# Nạp, kiểm tra và phân tích hai bộ dữ liệu thực.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

ENERGY_COLUMNS <- c(
  "relative_compactness", "surface_area", "wall_area", "roof_area",
  "overall_height", "orientation", "glazing_area",
  "glazing_area_distribution", "heating_load", "cooling_load"
)
CONCRETE_COLUMNS <- c(
  "cement", "blast_furnace_slag", "fly_ash", "water",
  "superplasticizer", "coarse_aggregate", "fine_aggregate", "age",
  "compressive_strength"
)

DEFAULT_DATA_PATHS <- list(
  energy = file.path("data", "ENB2012_data.xlsx"),
  concrete = file.path("data", "Concrete_Data.xls")
)

resolve_data_file <- function(project_root, data_path) {
  candidate <- path.expand(data_path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", candidate)) {
    candidate <- file.path(project_root, candidate)
  }
  if (!file.exists(candidate)) {
    stop(
      "Không tìm thấy tệp dữ liệu: ", candidate,
      ". Cần cập nhật hằng số đường dẫn trong cell Cấu hình."
    )
  }
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}



file_sha256 <- function(path) {
  commands <- list(
    list(command = "sha256sum", args = c(shQuote(path))),
    list(command = "shasum", args = c("-a", "256", shQuote(path)))
  )
  for (entry in commands) {
    executable <- Sys.which(entry$command)
    if (!nzchar(executable)) next
    output <- tryCatch(
      system2(executable, args = entry$args, stdout = TRUE, stderr = FALSE),
      error = function(error) character(0)
    )
    if (length(output)) {
      candidate <- strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]]
      if (grepl("^[0-9a-fA-F]{64}$", candidate)) return(tolower(candidate))
    }
  }
  NA_character_
}

load_real_dataset <- function(project_root, dataset, data_path = NULL) {
  if (!dataset %in% names(DEFAULT_DATA_PATHS)) {
    stop("Dataset must be energy or concrete.")
  }
  if (is.null(data_path)) data_path <- DEFAULT_DATA_PATHS[[dataset]]
  path <- resolve_data_file(project_root, data_path)
  if (identical(dataset, "energy")) {
    frame <- as.data.frame(readxl::read_excel(path))
    if (!identical(dim(frame), c(768L, 10L))) {
      stop("Unexpected Energy shape: ", paste(dim(frame), collapse = " x "))
    }
    names(frame) <- ENERGY_COLUMNS
    response <- "heating_load"
    x_name <- "relative_compactness"
    z_names <- setdiff(ENERGY_COLUMNS[1:8], c(x_name, "roof_area"))
    dropped <- "roof_area"
    relation <- "surface_area = wall_area + 2 * roof_area"
  } else if (identical(dataset, "concrete")) {
    frame <- as.data.frame(readxl::read_excel(path))
    if (!identical(dim(frame), c(1030L, 9L))) {
      stop("Unexpected Concrete shape: ", paste(dim(frame), collapse = " x "))
    }
    names(frame) <- CONCRETE_COLUMNS
    response <- "compressive_strength"
    x_name <- "water"
    z_names <- setdiff(CONCRETE_COLUMNS[1:8], x_name)
    dropped <- character(0)
    relation <- ""
  }
  frame[] <- lapply(frame, function(column) as.numeric(column))
  if (anyNA(frame) || any(!is.finite(as.matrix(frame)))) {
    stop(dataset, " contains missing or non-finite values.")
  }
  metadata <- list(
    dataset = dataset,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    rows = nrow(frame),
    columns = ncol(frame),
    response = response,
    x = x_name,
    z = z_names,
    dropped_collinear_z = dropped,
    deterministic_relation = relation,
    column_names = names(frame),
    response_variance = stats::var(frame[[response]]),
    file_size_bytes = file.info(path)$size,
    sha256 = file_sha256(path)
  )
  list(frame = frame, metadata = metadata)
}

real_data_object <- function(frame, metadata) {
  list(
    x = as.matrix(frame[metadata$x]),
    z = as.matrix(frame[metadata$z]),
    y = as.numeric(frame[[metadata$response]]),
    basis = "identity",
    targets = list()
  )
}

cross_fitted_predictability <- function(data, folds, learner, config, seed) {
  x <- as_2d(data$x)
  z <- as_2d(data$z)
  predictions <- matrix(NA_real_, nrow(x), ncol(x))
  backend <- ""
  for (fold in seq_len(config$folds)) {
    train <- folds != fold
    test <- folds == fold
    model <- fit_regressor(
      learner,
      z[train, , drop = FALSE],
      x[train, 1L],
      seed_from(seed, fold),
      config$forest_trees,
      config$n_threads
    )
    predictions[test, 1L] <- model$predict(z[test, , drop = FALSE])
    backend <- model$label
  }
  residual <- x[, 1L] - predictions[, 1L]
  total <- sum((x[, 1L] - mean(x[, 1L]))^2)
  data.frame(
    learner = learner,
    backend = backend,
    predictability_r2 = 1 - sum(residual^2) / total,
    residual_variance_ratio = stats::var(residual) / stats::var(x[, 1L]),
    max_abs_pearson_r = max(abs(stats::cor(x[, 1L], z))),
    stringsAsFactors = FALSE
  )
}


summarize_real_diagnostics <- function(results, predictability, metadata) {
  diagnostics <- merge(
    results,
    predictability[c(
      "learner", "predictability_r2", "residual_variance_ratio",
      "max_abs_pearson_r"
    )],
    by = "learner",
    all.x = TRUE,
    sort = FALSE
  )
  scale_value <- metadata$response_variance
  diagnostics$estimate_scaled_by_var_y <- diagnostics$estimate / scale_value
  diagnostics$ci_low_cross_scaled_by_var_y <- diagnostics$ci_low_cross / scale_value
  diagnostics$ci_high_cross_scaled_by_var_y <- diagnostics$ci_high_cross / scale_value
  diagnostics$ci_width_cross_scaled_by_var_y <- diagnostics$ci_width_cross / scale_value
  diagnostics$correction_width_ratio <- diagnostics$ci_width /
    diagnostics$ci_width_cross
  diagnostics$diagnostic <- vapply(seq_len(nrow(diagnostics)), function(index) {
    row <- diagnostics[index, , drop = FALSE]
    flags <- character(0)
    if (
      !identical(as.character(row$identification_status), "identified") ||
        !is.finite(row$estimate) || row$finite_fraction < 1
    ) {
      flags <- c(
        flags,
        "ước lượng không định danh hoặc không hữu hạn trên toàn bộ fold"
      )
    }
    if (is.finite(row$predictability_r2) && row$predictability_r2 >= 0.90) {
      flags <- c(flags, "X được dự đoán rất mạnh từ Z")
    }
    if (
      is.finite(row$residual_variance_ratio) &&
        row$residual_variance_ratio < 0.10
    ) {
      flags <- c(flags, "tỷ lệ phương sai phần dư của X dưới 0,10")
    }
    if (
      is.finite(row$density_ratio_ess_fraction_mean) &&
        row$density_ratio_ess_fraction_mean < 0.10
    ) {
      flags <- c(flags, "ESS của tỷ số mật độ dưới 10% cỡ mẫu")
    }
    if (is.finite(row$correction_width_ratio) && row$correction_width_ratio > 10) {
      flags <- c(
        flags,
        "khoảng có hiệu chỉnh rộng hơn khoảng t-Cross trên 10 lần"
      )
    }
    if (length(flags)) {
      paste(flags, collapse = "; ")
    } else {
      "không phát hiện cảnh báo theo các ngưỡng định trước"
    }
  }, character(1))
  diagnostics
}

summarize_learner_sensitivity <- function(diagnostics) {
  groups <- split(diagnostics, diagnostics$estimator)
  result <- do.call(rbind, lapply(groups, function(group) {
    values <- group$estimate_scaled_by_var_y
    finite <- values[is.finite(values)]
    median_value <- if (length(finite)) stats::median(finite) else NA_real_
    range_value <- if (length(finite)) diff(range(finite)) else NA_real_
    relative_range <- if (
      is.finite(range_value) && is.finite(median_value)
    ) {
      range_value / max(abs(median_value), 0.01)
    } else NA_real_
    data.frame(
      estimator = group$estimator[[1L]],
      learners_available = length(finite),
      minimum_scaled_estimate = if (length(finite)) min(finite) else NA_real_,
      maximum_scaled_estimate = if (length(finite)) max(finite) else NA_real_,
      range_across_learners = range_value,
      median_scaled_estimate = median_value,
      relative_range_across_learners = relative_range,
      learner_sensitivity_flag = if (!length(finite)) {
        "không có ước lượng hữu hạn"
      } else if (length(finite) < 3L) {
        "không đủ ba learner để đánh giá độ nhạy"
      } else if (is.finite(relative_range) && relative_range > 1) {
        "độ biến thiên giữa learner vượt độ lớn trung vị"
      } else {
        "độ biến thiên giữa learner không vượt độ lớn trung vị"
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}

build_real_report_table <- function(
  diagnostics,
  learner_sensitivity,
  metadata
) {
  sensitivity_columns <- learner_sensitivity[c(
    "estimator", "range_across_learners",
    "relative_range_across_learners", "learner_sensitivity_flag"
  )]
  report <- merge(
    diagnostics,
    sensitivity_columns,
    by = "estimator",
    all.x = TRUE,
    sort = FALSE
  )
  report <- data.frame(
    dataset = metadata$dataset,
    response = metadata$response,
    x_feature = metadata$x,
    learner = report$learner,
    estimator = report$estimator,
    estimate_scaled_by_var_y = report$estimate_scaled_by_var_y,
    ci_low_cross_scaled_by_var_y = report$ci_low_cross_scaled_by_var_y,
    ci_high_cross_scaled_by_var_y = report$ci_high_cross_scaled_by_var_y,
    predictability_r2 = report$predictability_r2,
    residual_variance_ratio = report$residual_variance_ratio,
    density_ratio_ess_fraction_mean = report$density_ratio_ess_fraction_mean,
    finite_fraction = report$finite_fraction,
    identification_status = report$identification_status,
    range_across_learners = report$range_across_learners,
    relative_range_across_learners = report$relative_range_across_learners,
    learner_sensitivity_flag = report$learner_sensitivity_flag,
    diagnostic = report$diagnostic,
    stringsAsFactors = FALSE
  )
  report[
    order(report$learner, report$estimator),
    ,
    drop = FALSE
  ]
}

run_real_dataset <- function(
  project_root,
  dataset,
  config = get_profile(),
  learners = c("linear", "additive", "forest"),
  language = "R",
  data_path = NULL
) {
  learners <- as.character(learners)
  loaded <- load_real_dataset(project_root, dataset, data_path = data_path)
  frame <- loaded$frame
  metadata <- loaded$metadata
  data <- real_data_object(frame, metadata)
  result_dir <- file.path(project_root, "results", language)
  checkpoint_dir <- file.path(project_root, "checkpoints", language, "real_data")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  fold_seed <- seed_from(
    nrow(frame), if (dataset == "energy") 811L else 977L
  )
  folds_path <- file.path(checkpoint_dir, paste0(dataset, "_folds.rds"))
  folds <- if (file.exists(folds_path)) {
    readRDS(folds_path)
  } else {
    value <- make_folds(nrow(frame), config$folds, fold_seed)
    saveRDS(value, folds_path)
    value
  }
  cache_base <- file.path(
    checkpoint_dir, paste0(dataset, "_psi0_mc", config$psi0_mc_samples)
  )
  caches <- load_or_prepare_psi0_caches(
    cache_base, data, folds, config, fold_seed
  )
  output <- file.path(result_dir, paste0(dataset, "_", config$name, ".csv"))
  fold_output <- file.path(
    result_dir, paste0(dataset, "_", config$name, "_folds.csv")
  )
  existing <- if (file.exists(output)) read.csv(output, check.names = FALSE) else data.frame()
  existing_folds <- if (file.exists(fold_output)) {
    read.csv(fold_output, check.names = FALSE)
  } else data.frame()
  required_schema <- c(
    "identification_status", "density_ratio_ess_mean", "density_ratio_max",
    "shared_runtime_seconds", "covers_psi0_cross"
  )
  if (
    nrow(existing) &&
    (!all(required_schema %in% names(existing)) || !nrow(existing_folds))
  ) {
    existing <- data.frame()
    existing_folds <- data.frame()
  }
  completed <- if (nrow(existing)) unique(existing$learner) else character(0)
  completed <- intersect(completed, learners)
  pieces <- if (nrow(existing)) list(existing) else list()
  fold_pieces <- if (nrow(existing_folds)) list(existing_folds) else list()
  predictability <- list()
  progress_started <- proc.time()[["elapsed"]]
  progress_count <- length(completed)
  initial_completed <- progress_count
  show_progress(
    paste0(tools::toTitleCase(dataset), " dataset"),
    0L,
    length(learners),
    progress_started,
    initial_completed = initial_completed
  )
  for (learner in learners) {
    predictability[[length(predictability) + 1L]] <- cross_fitted_predictability(
      data, folds, learner, config, fold_seed
    )
    if (learner %in% completed) next
    result <- t_cross(
      data,
      learner,
      fold_seed,
      config,
      fold_assignments = folds,
      psi0_caches = caches,
      include_psi0 = TRUE,
      include_baselines = TRUE
    )
    summary <- data.frame(
      language = language,
      profile = config$name,
      task = "real_data",
      dataset = dataset,
      x_feature = metadata$x,
      response = metadata$response,
      seed = fold_seed,
      n = nrow(frame),
      result$summary,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    fold_frame <- data.frame(
      language = language,
      profile = config$name,
      task = "real_data",
      dataset = dataset,
      x_feature = metadata$x,
      response = metadata$response,
      seed = fold_seed,
      n = nrow(frame),
      result$folds,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    pieces[[length(pieces) + 1L]] <- summary
    fold_pieces[[length(fold_pieces) + 1L]] <- fold_frame
    write.csv(do.call(rbind, pieces), output, row.names = FALSE, na = "")
    write.csv(
      do.call(rbind, fold_pieces),
      fold_output,
      row.names = FALSE,
      na = ""
    )
    progress_count <- progress_count + 1L
    show_progress(
      paste0(tools::toTitleCase(dataset), " dataset"),
      progress_count,
      length(learners),
      progress_started,
      paste0("learner=", learner),
      initial_completed
    )
  }
  predictability_frame <- do.call(rbind, predictability)
  predictability_frame <- data.frame(
    dataset = dataset,
    x_feature = metadata$x,
    predictability_frame,
    stringsAsFactors = FALSE
  )
  write.csv(
    predictability_frame,
    file.path(result_dir, paste0(dataset, "_", config$name, "_predictability.csv")),
    row.names = FALSE
  )
  jsonlite::write_json(
    metadata,
    file.path(result_dir, paste0(dataset, "_metadata.json")),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  results_frame <- read.csv(output, check.names = FALSE)
  diagnostics <- summarize_real_diagnostics(
    results_frame, predictability_frame, metadata
  )
  sensitivity <- summarize_learner_sensitivity(diagnostics)
  report_table <- build_real_report_table(
    diagnostics, sensitivity, metadata
  )
  write.csv(
    diagnostics,
    file.path(result_dir, paste0(dataset, "_", config$name, "_diagnostics.csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    sensitivity,
    file.path(result_dir, paste0(dataset, "_", config$name, "_learner_sensitivity.csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    report_table,
    file.path(result_dir, paste0(dataset, "_", config$name, "_report_table.csv")),
    row.names = FALSE,
    na = ""
  )
  list(
    results = results_frame,
    predictability = predictability_frame,
    diagnostics = diagnostics,
    learner_sensitivity = sensitivity,
    report_table = report_table,
    metadata = metadata
  )
}
