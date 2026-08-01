# DGP của Kịch bản mô phỏng 1-5, checkpoint và các vòng chạy mô phỏng.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

analytical_targets <- function(example, rho = 0.7) {
  if (example == 1L) {
    psi_l <- 4 * (1 - rho^2)
    return(list(
      psi_L = psi_l,
      psi_0 = 4,
      psi_1 = if (rho <= 0.5) psi_l else 4,
      psi_2 = 4,
      psi_3 = 4
    ))
  }
  if (example == 2L) return(list(psi_0 = 60))
  if (example == 3L) return(list(psi_0 = 100))
  if (example == 4L) {
    return(list(psi_0 = 1 / 7 + 49 / 125 - 49 / 225))
  }
  if (example == 5L) return(list(psi_0 = 9.16))
  stop("Example must be an integer from one to five.")
}

generate_example <- function(example, n, seed, rho = 0.7, h = 5L) {
  set.seed(seed)
  noise <- stats::rnorm(n)
  if (example == 1L) {
    if (rho < 0 || rho >= 1) stop("Example 1 requires 0 <= rho < 1.")
    x <- stats::rnorm(n)
    delta <- rho / sqrt(1 - rho^2)
    z <- cbind(
      delta * x + stats::rnorm(n),
      matrix(stats::rnorm(n * (h - 1L)), nrow = n)
    )
    y <- 2 * x + noise
    basis <- "identity"
  } else if (example == 2L) {
    x <- stats::rnorm(n)
    z <- cbind(
      x + stats::rnorm(n, sd = 0.4),
      matrix(stats::rnorm(n * (h - 1L)), nrow = n)
    )
    y <- 2 * x^3 + noise
    basis <- "orthogonal_polynomial"
  } else if (example == 3L) {
    z <- matrix(stats::rnorm(n * h), nrow = n)
    x <- cbind(
      2 * z[, 1L] + stats::rnorm(n),
      2 * z[, 2L] + stats::rnorm(n)
    )
    y <- 2 * x[, 1L] * x[, 2L] + noise
    basis <- "identity"
  } else if (example == 4L) {
    x <- stats::runif(n, -1, 1)
    z <- matrix(stats::runif(n * h, -1, 1), nrow = n)
    y <- x^2 * (x + 7 / 5) + (25 / 9) * z[, 1L]^2 + noise
    basis <- "orthogonal_polynomial"
  } else if (example == 5L) {
    x <- stats::rnorm(n)
    z <- cbind(
      x + stats::rnorm(n, sd = 0.4),
      matrix(stats::rnorm(n * (h - 1L)), nrow = n)
    )
    y <- 2 * x^2 + x * z[, 1L] + noise
    basis <- "orthogonal_polynomial"
  } else {
    stop("Example must be an integer from one to five.")
  }
  list(
    x = as_2d(x),
    z = as_2d(z),
    y = as.numeric(y),
    basis = basis,
    targets = analytical_targets(example, rho),
    example = example,
    rho = rho
  )
}

data_cache_path <- function(
  root, language, example, repetition, seed, n, rho
) {
  rho_tag <- gsub("\\.", "p", sprintf("%.2f", rho))
  file.path(
    root, "checkpoints", language, "data",
    paste0(
      "example", example, "_rep", repetition, "_seed", seed,
      "_n", n, "_rho", rho_tag, ".rds"
    )
  )
}

load_or_generate_example <- function(
  root, language, example, repetition, seed, n, rho
) {
  path <- data_cache_path(root, language, example, repetition, seed, n, rho)
  if (file.exists(path)) return(readRDS(path))
  data <- generate_example(example, n, seed, rho)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data, path)
  data
}

load_or_prepare_psi0_caches <- function(
  base, data, folds, config, seed
) {
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  x <- as_2d(data$x)
  z <- as_2d(data$z)
  caches <- vector("list", config$folds)
  for (fold in seq_len(config$folds)) {
    path <- file.path(base, paste0("psi0_fold", fold, ".rds"))
    if (file.exists(path)) {
      caches[[fold]] <- readRDS(path)
    } else {
      train <- folds != fold
      test <- folds == fold
      caches[[fold]] <- tryCatch(
        prepare_psi0_cache(
          x[train, , drop = FALSE],
          x[test, , drop = FALSE],
          z[train, , drop = FALSE],
          z[test, , drop = FALSE],
          config$psi0_mc_samples,
          seed_from(seed, fold, 77L),
          config$density_ratio_cap
        ),
        error = function(error) list(error = conditionMessage(error))
      )
      saveRDS(caches[[fold]], path)
    }
  }
  caches
}

simulation_diagnostic_flag <- function(
  finite_fraction,
  relative_absolute_bias,
  coverage,
  coverage_trials,
  mean_ci_width,
  target,
  density_ratio_ess_fraction = NA_real_
) {
  flags <- character(0)
  if (is.finite(finite_fraction) && finite_fraction < 0.95) {
    flags <- c(flags, "tỷ lệ ước lượng hữu hạn dưới 0,95")
  }
  if (is.finite(relative_absolute_bias) && relative_absolute_bias > 0.20) {
    flags <- c(flags, "độ lệch tương đối tuyệt đối vượt 20%")
  }
  if (
    is.finite(density_ratio_ess_fraction) &&
      density_ratio_ess_fraction < 0.10
  ) {
    flags <- c(flags, "ESS của tỷ số mật độ dưới 10% cỡ mẫu")
  }
  if (is.finite(coverage_trials) && coverage_trials < 30L) {
    flags <- c(flags, "số lần lặp dưới 30; coverage chỉ mang tính mô tả")
  } else if (is.finite(coverage) && coverage < 0.90) {
    flags <- c(flags, "coverage t-Cross dưới 0,90")
  }
  if (
    is.finite(mean_ci_width) && is.finite(target) && abs(target) > 1e-12 &&
      mean_ci_width / abs(target) > 1
  ) {
    flags <- c(
      flags,
      "độ rộng khoảng t-Cross trung bình vượt độ lớn của target"
    )
  }
  if (length(flags)) {
    paste(flags, collapse = "; ")
  } else {
    "không phát hiện cảnh báo theo các ngưỡng định trước"
  }
}

aggregate_simulation_results <- function(results) {
  if (!nrow(results)) return(results)
  keys <- c("language", "profile", "example", "rho", "learner", "estimator")
  split_rows <- split(results, interaction(results[keys], drop = TRUE))
  aggregate <- do.call(rbind, lapply(split_rows, function(group) {
    estimates <- as.numeric(group$estimate)
    target_values <- as.numeric(group$target_psi0)
    target_values <- target_values[is.finite(target_values)]
    target <- if (length(target_values)) target_values[[1L]] else NA_real_

    cross_values <- group$covers_psi0_cross[
      !is.na(group$covers_psi0_cross)
    ]
    cross_trials <- length(cross_values)
    cross_successes <- if (cross_trials) sum(cross_values) else 0L
    cross_interval <- exact_binomial_interval(cross_successes, cross_trials)
    cross_coverage <- if (cross_trials) {
      cross_successes / cross_trials
    } else NA_real_

    literal_values <- group$covers_psi0[!is.na(group$covers_psi0)]
    literal_trials <- length(literal_values)
    literal_successes <- if (literal_trials) sum(literal_values) else 0L
    literal_interval <- exact_binomial_interval(
      literal_successes, literal_trials
    )
    literal_coverage <- if (literal_trials) {
      literal_successes / literal_trials
    } else NA_real_

    mean_estimate <- safe_mean(estimates)
    bias <- if (is.finite(mean_estimate) && is.finite(target)) {
      mean_estimate - target
    } else NA_real_
    relative_absolute_bias <- if (is.finite(bias) && abs(target) > 1e-12) {
      abs(bias) / abs(target)
    } else NA_real_
    rmse <- if (is.finite(target)) {
      sqrt(safe_mean((estimates - target)^2))
    } else NA_real_
    relative_rmse <- if (is.finite(rmse) && abs(target) > 1e-12) {
      rmse / abs(target)
    } else NA_real_

    mean_width_cross <- safe_mean(group$ci_width_cross)
    mean_width_literal <- safe_mean(group$ci_width)
    width_inflation <- if (
      is.finite(mean_width_cross) && mean_width_cross > 0 &&
        is.finite(mean_width_literal)
    ) {
      mean_width_literal / mean_width_cross
    } else NA_real_
    ess_fraction <- safe_mean(group$density_ratio_ess_fraction_mean)

    data.frame(
      language = group$language[[1L]],
      profile = group$profile[[1L]],
      example = group$example[[1L]],
      rho = group$rho[[1L]],
      learner = group$learner[[1L]],
      estimator = group$estimator[[1L]],
      repetitions = length(unique(group$repetition)),
      target_psi0 = target,
      mean_estimate = mean_estimate,
      bias = bias,
      absolute_bias = if (is.finite(bias)) abs(bias) else NA_real_,
      relative_absolute_bias = relative_absolute_bias,
      rmse = rmse,
      relative_rmse = relative_rmse,
      standard_deviation = safe_sd(estimates),

      # Khoảng t-Cross không hiệu chỉnh là quy ước trình bày chính.
      coverage = cross_coverage,
      coverage_ci_low = cross_interval[[1L]],
      coverage_ci_high = cross_interval[[2L]],
      coverage_mcse = if (cross_trials) {
        sqrt(cross_coverage * (1 - cross_coverage) / cross_trials)
      } else NA_real_,
      coverage_trials = cross_trials,
      mean_ci_width = mean_width_cross,
      median_ci_width = safe_median(group$ci_width_cross),

      # Các cột tương thích được giữ lại để hỗ trợ mã hậu xử lý.
      coverage_cross = cross_coverage,
      coverage_cross_ci_low = cross_interval[[1L]],
      coverage_cross_ci_high = cross_interval[[2L]],
      mean_ci_width_cross = mean_width_cross,
      median_ci_width_cross = safe_median(group$ci_width_cross),

      # Khoảng theo cách diễn giải trực tiếp c = Var(Y)^2 là phân tích độ nhạy.
      coverage_paper_literal = literal_coverage,
      coverage_paper_literal_ci_low = literal_interval[[1L]],
      coverage_paper_literal_ci_high = literal_interval[[2L]],
      mean_ci_width_paper_literal = mean_width_literal,
      median_ci_width_paper_literal = safe_median(group$ci_width),
      interval_width_inflation_paper_literal = width_inflation,

      finite_fraction = safe_mean(group$finite_fraction),
      density_ratio_ess_fraction_mean = ess_fraction,
      shared_runtime_seconds = safe_mean(group$shared_runtime_seconds) *
        length(unique(group$repetition)),
      diagnostic = simulation_diagnostic_flag(
        safe_mean(group$finite_fraction),
        relative_absolute_bias,
        cross_coverage,
        cross_trials,
        mean_width_cross,
        target,
        ess_fraction
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(aggregate) <- NULL
  aggregate <- merge(
    aggregate,
    TABLE2_REFERENCE,
    by = c("example", "learner", "estimator"),
    all.x = TRUE,
    sort = FALSE
  )
  aggregate$coverage_minus_paper <- aggregate$coverage -
    aggregate$paper_coverage
  aggregate
}

build_simulation_diagnostic_table <- function(aggregate) {
  columns <- c(
    "example", "rho", "learner", "estimator", "repetitions",
    "target_psi0", "mean_estimate", "relative_absolute_bias",
    "relative_rmse", "standard_deviation",
    "coverage", "coverage_ci_low", "coverage_ci_high",
    "paper_coverage", "coverage_minus_paper",
    "mean_ci_width", "coverage_paper_literal",
    "coverage_paper_literal_ci_low", "coverage_paper_literal_ci_high",
    "interval_width_inflation_paper_literal", "finite_fraction", "density_ratio_ess_fraction_mean", "diagnostic"
  )
  aggregate[
    aggregate$example >= 2L & aggregate$example <= 5L,
    columns,
    drop = FALSE
  ]
}

run_simulations <- function(
  project_root,
  config = get_profile(),
  learners = c("linear", "additive", "forest"),
  language = "R"
) {
  learners <- as.character(learners)
  result_dir <- file.path(project_root, "results", language)
  checkpoint_dir <- file.path(project_root, "checkpoints", language)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(result_dir, paste0("simulation_", config$name, ".csv"))
  fold_output <- file.path(
    result_dir, paste0("simulation_", config$name, "_folds.csv")
  )
  existing <- if (file.exists(output)) {
    read.csv(output, check.names = FALSE)
  } else data.frame()
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
  completed <- character(0)
  if (nrow(existing)) {
    completed <- unique(paste(
      existing$example, existing$rho, existing$repetition, existing$learner,
      sep = "|"
    ))
  }
  pieces <- if (nrow(existing)) list(existing) else list()
  fold_pieces <- if (nrow(existing_folds)) list(existing_folds) else list()
  expected_keys <- unlist(lapply(config$examples, function(example) {
    correlations <- if (example == 1L) config$correlations else 0.7
    unlist(lapply(correlations, function(rho) {
      unlist(lapply(seq_len(config$repetitions), function(repetition) {
        paste(example, rho, repetition, learners, sep = "|")
      }))
    }))
  }))
  completed <- intersect(completed, expected_keys)
  total_tasks <- length(expected_keys)
  progress_started <- proc.time()[["elapsed"]]
  progress_count <- length(completed)
  initial_completed <- progress_count
  show_progress(
    "Paper simulations",
    progress_count,
    total_tasks,
    progress_started,
    initial_completed = initial_completed
  )

  for (example in config$examples) {
    correlations <- if (example == 1L) config$correlations else 0.7
    for (rho in correlations) {
      for (repetition in seq_len(config$repetitions)) {
        data_seed <- seed_from(
          example, repetition, config$n, as.integer(round(rho * 1000))
        )
        data <- load_or_generate_example(
          project_root, language, example, repetition, data_seed, config$n, rho
        )
        fold_seed <- seed_from(data_seed, 101L)
        folds_path <- file.path(
          checkpoint_dir, "folds",
          paste0("example", example, "_rep", repetition, "_rho", sprintf("%.2f", rho), ".rds")
        )
        dir.create(dirname(folds_path), recursive = TRUE, showWarnings = FALSE)
        folds <- if (file.exists(folds_path)) {
          readRDS(folds_path)
        } else {
          value <- make_folds(config$n, config$folds, fold_seed)
          saveRDS(value, folds_path)
          value
        }
        psi0_base <- file.path(
          checkpoint_dir, "psi0",
          paste0(
            "example", example, "_rep", repetition, "_rho",
            sprintf("%.2f", rho), "_mc", config$psi0_mc_samples
          )
        )
        caches <- load_or_prepare_psi0_caches(
          psi0_base, data, folds, config, fold_seed
        )
        for (learner in learners) {
          key <- paste(example, rho, repetition, learner, sep = "|")
          detail <- paste0(
            "example=", example,
            ", rho=", sprintf("%.2f", rho),
            ", repetition=", repetition,
            ", learner=", learner
          )
          if (key %in% completed) next
          result <- t_cross(
            data,
            learner,
            fold_seed,
            config,
            fold_assignments = folds,
            psi0_caches = caches,
            include_psi0 = TRUE
          )
          summary <- data.frame(
            language = language,
            profile = config$name,
            task = "paper_simulation",
            example = example,
            rho = rho,
            repetition = repetition,
            seed = data_seed,
            n = config$n,
            result$summary,
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
          fold_frame <- data.frame(
            language = language,
            profile = config$name,
            task = "paper_simulation",
            example = example,
            rho = rho,
            repetition = repetition,
            seed = data_seed,
            n = config$n,
            result$folds,
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
          pieces[[length(pieces) + 1L]] <- summary
          fold_pieces[[length(fold_pieces) + 1L]] <- fold_frame
          current <- do.call(rbind, pieces)
          write.csv(current, output, row.names = FALSE, na = "")
          write.csv(
            do.call(rbind, fold_pieces),
            fold_output,
            row.names = FALSE,
            na = ""
          )
          progress_count <- progress_count + 1L
          show_progress(
            "Paper simulations",
            progress_count,
            total_tasks,
            progress_started,
            detail,
            initial_completed
          )
        }
      }
    }
  }
  results <- if (file.exists(output)) read.csv(output, check.names = FALSE) else data.frame()
  aggregate <- aggregate_simulation_results(results)
  write.csv(
    aggregate,
    file.path(result_dir, paste0("simulation_", config$name, "_aggregate.csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    TABLE2_REFERENCE,
    file.path(result_dir, "table2_reference.csv"),
    row.names = FALSE
  )
  list(results = results, aggregate = aggregate)
}

baseline_method_role <- function(estimator) {
  roles <- c(
    psi_L = "LOCO; estimand chịu ảnh hưởng của phân phối chung X,Z",
    cpi_residual = "Xấp xỉ CPI bằng hoán vị phần dư của mô hình X|Z",
    pfi = "Marginal PFI; phép nhiễu không cùng estimand với tham số decorrelated"
  )
  unname(roles[estimator])
}

run_baseline_benchmark <- function(
  project_root,
  config = get_profile(),
  learners = c("linear", "additive", "forest"),
  language = "R"
) {
  result_dir <- file.path(project_root, "results", language)
  checkpoint_dir <- file.path(project_root, "checkpoints", language)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(
    result_dir, paste0("baseline_benchmark_", config$name, ".csv")
  )
  existing <- if (file.exists(output)) {
    read.csv(output, check.names = FALSE)
  } else data.frame()
  required_schema <- c(
    "example", "rho", "repetition", "learner", "estimator",
    "covers_psi0_cross", "shared_runtime_seconds"
  )
  if (nrow(existing) && !all(required_schema %in% names(existing))) {
    existing <- data.frame()
  }
  expected_estimators <- BASELINE_ESTIMATORS
  completed <- character(0)
  if (nrow(existing)) {
    key <- paste(
      existing$example, existing$rho, existing$repetition, existing$learner,
      sep = "|"
    )
    counts <- tapply(existing$estimator, key, function(x) {
      all(expected_estimators %in% unique(x))
    })
    completed <- names(counts)[counts]
    existing <- existing[key %in% completed, , drop = FALSE]
  }
  pieces <- if (nrow(existing)) list(existing) else list()
  task_grid <- do.call(rbind, lapply(config$baseline_examples, function(example) {
    rho <- if (example == 1L) config$baseline_rho else 0.7
    expand.grid(
      example = example,
      rho = rho,
      repetition = seq_len(config$baseline_repetitions),
      learner = learners,
      stringsAsFactors = FALSE
    )
  }))
  task_grid$key <- paste(
    task_grid$example, task_grid$rho, task_grid$repetition,
    task_grid$learner, sep = "|"
  )
  progress_started <- proc.time()[["elapsed"]]
  progress_count <- sum(task_grid$key %in% completed)
  initial_completed <- progress_count
  show_progress(
    "Baseline benchmark",
    progress_count,
    nrow(task_grid),
    progress_started,
    initial_completed = initial_completed
  )
  for (row_index in seq_len(nrow(task_grid))) {
    task <- task_grid[row_index, , drop = FALSE]
    if (task$key %in% completed) next
    example <- as.integer(task$example)
    rho <- as.numeric(task$rho)
    repetition <- as.integer(task$repetition)
    learner <- as.character(task$learner)
    data_seed <- seed_from(
      example, repetition, config$n, as.integer(round(rho * 1000))
    )
    data <- load_or_generate_example(
      project_root, language, example, repetition, data_seed, config$n, rho
    )
    fold_seed <- seed_from(data_seed, 101L)
    folds_path <- file.path(
      checkpoint_dir, "folds",
      paste0(
        "example", example, "_rep", repetition,
        "_rho", sprintf("%.2f", rho), ".rds"
      )
    )
    dir.create(dirname(folds_path), recursive = TRUE, showWarnings = FALSE)
    folds <- if (file.exists(folds_path)) {
      readRDS(folds_path)
    } else {
      value <- make_folds(config$n, config$folds, fold_seed)
      saveRDS(value, folds_path)
      value
    }
    result <- t_cross(
      data,
      learner,
      fold_seed,
      config,
      fold_assignments = folds,
      include_psi0 = FALSE,
      include_baselines = TRUE,
      correction_mode = "none"
    )
    summary <- result$summary[
      result$summary$estimator %in% BASELINE_ESTIMATORS,
      , drop = FALSE
    ]
    summary <- data.frame(
      language = language,
      profile = config$name,
      task = "baseline_benchmark",
      example = example,
      rho = rho,
      repetition = repetition,
      seed = data_seed,
      n = config$n,
      summary,
      method_role = baseline_method_role(summary$estimator),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    pieces[[length(pieces) + 1L]] <- summary
    write.csv(do.call(rbind, pieces), output, row.names = FALSE, na = "")
    progress_count <- progress_count + 1L
    show_progress(
      "Baseline benchmark",
      progress_count,
      nrow(task_grid),
      progress_started,
      paste0(
        "example=", example,
        ", repetition=", repetition,
        ", learner=", learner
      ),
      initial_completed
    )
  }
  results <- read.csv(output, check.names = FALSE)
  aggregate <- aggregate_simulation_results(results)
  aggregate$method_role <- baseline_method_role(aggregate$estimator)
  aggregate_path <- file.path(
    result_dir, paste0("baseline_benchmark_", config$name, "_aggregate.csv")
  )
  write.csv(aggregate, aggregate_path, row.names = FALSE, na = "")
  list(results = results, aggregate = aggregate)
}

assemble_method_comparison <- function(simulation_aggregate, baseline_aggregate) {
  dvi <- simulation_aggregate[
    simulation_aggregate$estimator == "psi_0" &
      simulation_aggregate$example %in% unique(baseline_aggregate$example),
    , drop = FALSE
  ]
  if (nrow(dvi)) {
    baseline_rho_lookup <- unique(baseline_aggregate[c("example", "rho")])
    dvi <- merge(
      dvi,
      baseline_rho_lookup,
      by = c("example", "rho"),
      all = FALSE,
      sort = FALSE
    )
  }
  comparison <- rbind(
    dvi,
    baseline_aggregate[names(dvi)]
  )
  comparison$method_label <- unname(c(
    psi_0 = "DVI psi_0",
    psi_L = "LOCO",
    cpi_residual = "Residual CPI",
    pfi = "Marginal PFI"
  )[comparison$estimator])
  comparison
}

run_correction_ablation <- function(
  project_root,
  config = get_profile(),
  language = "R",
  example = 2L,
  learner = "linear"
) {
  result_dir <- file.path(project_root, "results", language)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  raw_path <- file.path(
    result_dir, paste0("correction_ablation_", config$name, ".csv")
  )
  existing <- if (file.exists(raw_path)) {
    read.csv(raw_path, check.names = FALSE)
  } else data.frame()
  modes <- c("none", "variance_scale", "paper_literal")
  expected_count <- length(modes) * 4L
  completed <- integer(0)
  if (nrow(existing) && all(c("repetition", "correction_mode") %in% names(existing))) {
    counts <- table(existing$repetition)
    completed <- as.integer(names(counts)[counts >= expected_count])
    existing <- existing[existing$repetition %in% completed, , drop = FALSE]
  } else {
    existing <- data.frame()
  }
  pieces <- if (nrow(existing)) list(existing) else list()
  ablation_config <- config
  ablation_config$name <- paste0(config$name, "_ablation")
  ablation_config$n <- as.integer(config$ablation_n)
  ablation_config$forest_trees <- min(
    if (is.null(config$forest_trees)) 500L else config$forest_trees,
    100L
  )
  ablation_config$psi0_mc_samples <- min(config$psi0_mc_samples, 20L)
  progress_started <- proc.time()[["elapsed"]]
  progress_count <- length(completed)
  initial_completed <- progress_count
  show_progress(
    "Correction ablation",
    progress_count,
    config$ablation_repetitions,
    progress_started,
    initial_completed = initial_completed
  )
  for (repetition in seq_len(config$ablation_repetitions)) {
    if (repetition %in% completed) next
    data_seed <- seed_from(BASE_SEED, 900L, example, repetition, ablation_config$n)
    data <- generate_example(example, ablation_config$n, data_seed, rho = 0.7)
    fold_seed <- seed_from(data_seed, 901L)
    folds <- make_folds(ablation_config$n, ablation_config$folds, fold_seed)
    fit <- t_cross(
      data,
      learner,
      fold_seed,
      ablation_config,
      fold_assignments = folds,
      include_psi0 = FALSE,
      include_baselines = FALSE,
      correction_mode = "none"
    )
    for (mode in modes) {
      summary <- summarize_fold_estimates(
        fit$folds,
        data$y,
        ablation_config$folds,
        data$targets$psi_0,
        correction_mode = mode
      )
      summary <- summary[summary$estimator %in% c("psi_L", "psi_1", "psi_2", "psi_3"), ]
      summary <- data.frame(
        language = language,
        profile = config$name,
        task = "correction_ablation",
        example = example,
        learner = learner,
        repetition = repetition,
        seed = data_seed,
        n = ablation_config$n,
        summary,
        stringsAsFactors = FALSE
      )
      pieces[[length(pieces) + 1L]] <- summary
    }
    write.csv(do.call(rbind, pieces), raw_path, row.names = FALSE, na = "")
    progress_count <- progress_count + 1L
    show_progress(
      "Correction ablation",
      progress_count,
      config$ablation_repetitions,
      progress_started,
      paste0("repetition=", repetition),
      initial_completed
    )
  }
  raw <- read.csv(raw_path, check.names = FALSE)
  groups <- split(
    raw,
    interaction(raw$correction_mode, raw$estimator, drop = TRUE)
  )
  summary <- do.call(rbind, lapply(groups, function(group) {
    coverage_values <- group$covers_psi0[!is.na(group$covers_psi0)]
    trials <- length(coverage_values)
    successes <- if (trials) sum(coverage_values) else 0L
    interval <- exact_binomial_interval(successes, trials)
    target_values <- group$target_psi0[is.finite(group$target_psi0)]
    target <- if (length(target_values)) target_values[[1L]] else NA_real_
    mean_estimate <- safe_mean(group$estimate)
    mean_width <- safe_mean(group$ci_width)
    data.frame(
      correction_mode = group$correction_mode[[1L]],
      estimator = group$estimator[[1L]],
      repetitions = length(unique(group$repetition)),
      target_psi0 = target,
      coverage = if (trials) successes / trials else NA_real_,
      coverage_ci_low = interval[[1L]],
      coverage_ci_high = interval[[2L]],
      mean_estimate = mean_estimate,
      absolute_bias = if (
        is.finite(mean_estimate) && is.finite(target)
      ) abs(mean_estimate - target) else NA_real_,
      empirical_sd = safe_sd(group$estimate),
      mean_width = mean_width,
      median_width = safe_median(group$ci_width),
      width_p90 = safe_quantile(group$ci_width, 0.90),
      mean_width_to_target = if (
        is.finite(mean_width) && is.finite(target) && abs(target) > 1e-12
      ) mean_width / abs(target) else NA_real_,
      finite_fraction = safe_mean(group$finite_fraction),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL

  reference <- summary[
    summary$correction_mode == "none",
    c("estimator", "coverage", "mean_width"),
    drop = FALSE
  ]
  names(reference) <- c(
    "estimator", "coverage_none", "mean_width_none"
  )
  summary <- merge(summary, reference, by = "estimator", all.x = TRUE, sort = FALSE)
  summary$coverage_change_vs_none <- summary$coverage - summary$coverage_none
  summary$width_inflation_vs_none <- summary$mean_width / summary$mean_width_none
  summary$coverage_none <- NULL
  summary$mean_width_none <- NULL

  summary_path <- file.path(
    result_dir, paste0("correction_ablation_", config$name, "_summary.csv")
  )
  write.csv(summary, summary_path, row.names = FALSE, na = "")
  list(raw = raw, summary = summary)
}
