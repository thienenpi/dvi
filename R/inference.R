# Tổng hợp theo fold và thủ tục suy luận t-Cross (mục 3.6 của bài báo).
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

summarize_fold_estimates <- function(
  fold_frame,
  y,
  folds,
  target_psi0 = NA_real_,
  correction_mode = "paper"
) {
  critical <- stats::qt(0.975, df = folds - 1L)
  variance_y <- stats::var(y)
  normalized_mode <- switch(
    correction_mode,
    paper = "paper_literal",
    variance_scale_sensitivity = "variance_scale",
    correction_mode
  )
  c_value <- switch(
    normalized_mode,
    paper_literal = variance_y^2,
    variance_scale = variance_y,
    none = 0,
    stop("Unknown correction mode: ", correction_mode)
  )
  rows <- lapply(unique(fold_frame$estimator), function(estimator) {
    group <- fold_frame[fold_frame$estimator == estimator, , drop = FALSE]
    values <- group$estimate
    finite_fraction <- mean(is.finite(values))
    if (sum(is.finite(values)) == folds) {
      center <- mean(values)
      fold_variance <- stats::var(values)
      se_cross <- sqrt(fold_variance / folds)
      se <- sqrt(se_cross^2 + c_value^2 / length(y))
      low_cross <- center - critical * se_cross
      high_cross <- center + critical * se_cross
      low <- center - critical * se
      high <- center + critical * se
    } else {
      center <- se_cross <- se <- low_cross <- high_cross <- low <- high <- NA_real_
    }
    statuses <- if ("identification_status" %in% names(group)) {
      unique(as.character(group$identification_status))
    } else character(0)
    statuses <- sort(statuses[
      nzchar(statuses) & !statuses %in% c("identified", "NA")
    ])
    finite_metric <- function(name, reducer) {
      if (!name %in% names(group)) return(NA_real_)
      metric_values <- as.numeric(group[[name]])
      metric_values <- metric_values[is.finite(metric_values)]
      if (length(metric_values)) reducer(metric_values) else NA_real_
    }
    data.frame(
      estimator = estimator,
      estimate = center,
      standard_error = se,
      standard_error_cross = se_cross,
      ci_low = low,
      ci_high = high,
      ci_width = if (is.finite(low)) high - low else NA_real_,
      ci_low_cross = low_cross,
      ci_high_cross = high_cross,
      ci_width_cross = if (is.finite(low_cross)) high_cross - low_cross else NA_real_,
      target_psi0 = target_psi0,
      covers_psi0 = if (is.finite(low) && is.finite(target_psi0)) {
        low <= target_psi0 && target_psi0 <= high
      } else NA,
      covers_psi0_cross = if (
        is.finite(low_cross) && is.finite(target_psi0)
      ) {
        low_cross <= target_psi0 && target_psi0 <= high_cross
      } else NA,
      finite_fraction = finite_fraction,
      identification_status = if (length(statuses)) {
        paste(statuses, collapse = "; ")
      } else "identified",
      residual_variance_ratio_min = finite_metric("residual_variance_ratio", min),
      density_ratio_ess_mean = finite_metric("density_ratio_ess", mean),
      density_ratio_ess_fraction_mean = finite_metric(
        "density_ratio_ess_fraction", mean
      ),
      density_ratio_p99_max = finite_metric("density_ratio_p99", max),
      density_ratio_max = finite_metric("density_ratio_max", max),
      correction_mode = normalized_mode,
      c_value = c_value,
      fold_estimates = paste(format(values, digits = 12), collapse = ";"),
      # Thời gian này được dùng chung bởi mọi estimator trong cùng một lần t-Cross.
      shared_runtime_seconds = sum(group$total_fold_seconds),
      shared_nuisance_fit_seconds = sum(group$nuisance_fit_seconds),
      shared_psi0_seconds = sum(group$psi0_seconds),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

t_cross <- function(
  data,
  learner,
  seed,
  config = get_profile(),
  fold_assignments = NULL,
  psi0_caches = NULL,
  include_psi0 = TRUE,
  include_baselines = FALSE,
  correction_mode = "paper"
) {
  y <- as.numeric(data$y)
  folds <- if (is.null(fold_assignments)) {
    make_folds(length(y), config$folds, seed)
  } else as.integer(fold_assignments)
  if (length(folds) != length(y)) stop("Fold assignment length mismatch.")
  rows <- list()
  row_index <- 1L
  for (fold in seq_len(config$folds)) {
    result <- estimate_fold(
      data,
      folds,
      fold,
      learner,
      seed_from(seed, fold - 1L),
      config,
      if (is.null(psi0_caches)) NULL else psi0_caches[[fold]],
      include_psi0,
      include_baselines
    )
    for (estimator in names(result$estimates)) {
      status <- "identified"
      if (estimator == "psi_0" && !is.null(result$diagnostics$psi0_error)) {
        status <- paste0(
          "weak_overlap_density_failure: ", result$diagnostics$psi0_error
        )
      } else if (estimator == "psi_2" && !isTRUE(result$diagnostics$psi2_valid)) {
        status <- "non_identified_or_singular_psi2_moment"
      } else if (estimator == "psi_3" && !isTRUE(result$diagnostics$psi3_valid)) {
        status <- "non_identified_or_singular_psi3_moment"
      }
      diagnostic_value <- function(name, only_psi0 = FALSE) {
        if (only_psi0 && estimator != "psi_0") return(NA_real_)
        value <- result$diagnostics[[name]]
        if (is.null(value)) NA_real_ else as.numeric(value)
      }
      rows[[row_index]] <- data.frame(
        fold = fold,
        learner = learner,
        estimator = estimator,
        estimate = unname(result$estimates[[estimator]]),
        backend = result$diagnostics$mu_xz_backend,
        identification_status = status,
        residual_variance_ratio = diagnostic_value("residual_variance_ratio"),
        density_ratio_ess = diagnostic_value("density_ratio_ess", TRUE),
        density_ratio_ess_fraction = diagnostic_value(
          "density_ratio_ess_fraction", TRUE
        ),
        density_ratio_p99 = diagnostic_value("density_ratio_p99", TRUE),
        density_ratio_max = diagnostic_value("density_ratio_max", TRUE),
        nuisance_fit_seconds = result$timings$nuisance_fit_seconds,
        psi0_seconds = result$timings$psi0_seconds,
        estimator_postprocess_seconds = result$timings$estimator_postprocess_seconds,
        total_fold_seconds = result$timings$total_fold_seconds,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  fold_frame <- do.call(rbind, rows)
  target <- if (!is.null(data$targets$psi_0)) data$targets$psi_0 else NA_real_
  summary <- summarize_fold_estimates(
    fold_frame, y, config$folds, target, correction_mode
  )
  summary <- data.frame(
    learner = learner,
    backend = fold_frame$backend[[1L]],
    summary,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(summary = summary, folds = fold_frame)
}
