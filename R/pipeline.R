# Các bước (stage) của quy trình thực nghiệm.
# Mỗi stage nhận và trả về `state`, nhờ đó main.R chạy được từng phần độc lập
# mà vẫn giữ đúng thứ tự phụ thuộc như khi chạy notebook tuần tự.

# Trong notebook thì hiển thị hình trực tiếp; khi chạy bằng Rscript thì chỉ in
# đường dẫn vì không có kênh hiển thị.
display_figure <- function(path) {
  if (interactive() && requireNamespace("IRdisplay", quietly = TRUE)) {
    IRdisplay::display_png(file = path)
  } else {
    cat("Hình:", path, "\n")
  }
  invisible(path)
}

new_pipeline_state <- function(project_root, data_paths, language = "R") {
  list(
    project_root = project_root,
    data_paths = data_paths,
    language = language,
    figures = list()
  )
}

results_path <- function(state, filename) {
  path <- file.path(state$project_root, "results", state$language, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

# --- Stage: kiểm tra môi trường, cấu hình và tính toàn vẹn dữ liệu ------------

stage_check_environment <- function(state, profile, n_threads) {
  required_packages <- c("jsonlite", "readxl", "mgcv", "grf")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop(
      "Thiếu R package: ", paste(missing_packages, collapse = ", "),
      ". Cần cài đặt package trước khi chạy."
    )
  }
  if (!requireNamespace("IRdisplay", quietly = TRUE)) {
    message("IRdisplay không có sẵn; hình sẽ chỉ được ghi ra tệp.")
  }

  state$config <- get_profile(profile, n_threads)
  state$backend <- check_reproduction_backends()
  cat("Project root:", state$project_root, "\n")
  print(state$config)
  print(state$backend)
  for (dataset_name in c("energy", "concrete")) {
    loaded <- load_real_dataset(
      state$project_root,
      dataset_name,
      data_path = state$data_paths[[dataset_name]]
    )
    metadata <- loaded$metadata
    cat(
      dataset_name,
      ": path=", metadata$path,
      ", shape=", metadata$rows, " x ", metadata$columns,
      ", response=", metadata$response,
      ", X=", metadata$x,
      ", Z=", paste(metadata$z, collapse = ", "),
      ", SHA-256=", metadata$sha256,
      "\n",
      sep = ""
    )
  }
  state
}

# --- Stage: giá trị mục tiêu giải tích ----------------------------------------

stage_targets <- function(state) {
  target_table <- do.call(rbind, lapply(1:5, function(example) {
    values <- analytical_targets(example, rho = 0.7)
    data.frame(
      simulation_scenario = example,
      psi_L = if (!is.null(values$psi_L)) values$psi_L else NA_real_,
      psi_0 = values$psi_0,
      psi_1 = if (!is.null(values$psi_1)) values$psi_1 else NA_real_,
      psi_2 = if (!is.null(values$psi_2)) values$psi_2 else NA_real_,
      psi_3 = if (!is.null(values$psi_3)) values$psi_3 else NA_real_
    )
  }))
  print(target_table)
  write.csv(
    target_table,
    results_path(state, "analytical_targets.csv"),
    row.names = FALSE,
    na = ""
  )
  state$target_table <- target_table
  state$figures$target <- plot_analytic_targets(
    state$project_root, state$language
  )
  display_figure(state$figures$target)
  state
}

# --- Stage: mô phỏng chính và phân tích điều kiện hoạt động -------------------

stage_simulations <- function(state) {
  config <- state$config
  simulation <- run_simulations(
    state$project_root, config, language = state$language
  )
  main_coverage <- build_simulation_diagnostic_table(simulation$aggregate)
  main_coverage <- main_coverage[
    order(
      main_coverage$example,
      main_coverage$learner,
      main_coverage$estimator
    ),
    ,
    drop = FALSE
  ]
  print(main_coverage)
  write.csv(
    main_coverage,
    results_path(
      state, paste0("simulation_", config$name, "_report_table.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  assessment_groups <- split(
    main_coverage,
    interaction(
      main_coverage$learner,
      main_coverage$estimator,
      drop = TRUE
    )
  )
  method_condition_summary <- do.call(rbind, lapply(
    assessment_groups,
    function(group) {
      relative_rmse <- group$relative_rmse[is.finite(group$relative_rmse)]
      relative_bias <- group$relative_absolute_bias[
        is.finite(group$relative_absolute_bias)
      ]
      ess <- group$density_ratio_ess_fraction_mean[
        is.finite(group$density_ratio_ess_fraction_mean)
      ]
      finite_fraction <- group$finite_fraction[is.finite(group$finite_fraction)]
      data.frame(
        learner = group$learner[[1L]],
        estimator = group$estimator[[1L]],
        evaluated_scenarios = nrow(group),
        median_relative_rmse = if (length(relative_rmse)) {
          stats::median(relative_rmse)
        } else NA_real_,
        maximum_relative_absolute_bias = if (length(relative_bias)) {
          max(relative_bias)
        } else NA_real_,
        scenarios_with_relative_bias_at_most_20_percent = sum(
          group$relative_absolute_bias <= 0.20,
          na.rm = TRUE
        ),
        scenarios_with_finite_fraction_at_least_0_95 = sum(
          group$finite_fraction >= 0.95,
          na.rm = TRUE
        ),
        scenarios_with_density_ratio_ess_below_0_10 = if (length(ess)) {
          sum(ess < 0.10)
        } else NA_integer_,
        condition_assessment = if (
          length(finite_fraction) &&
            any(finite_fraction < 0.95)
        ) {
          "có kịch bản không hữu hạn hoặc không định danh"
        } else if (
          length(relative_bias) && max(relative_bias) > 1
        ) {
          "có kịch bản có độ lệch vượt độ lớn của target"
        } else if (
          length(relative_rmse) && stats::median(relative_rmse) <= 0.20
        ) {
          "RMSE tương đối trung vị không vượt 0,20"
        } else {
          "kết quả phụ thuộc đáng kể vào kịch bản hoặc learner"
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(method_condition_summary) <- NULL
  method_condition_summary <- method_condition_summary[
    order(
      method_condition_summary$learner,
      method_condition_summary$estimator
    ),
    ,
    drop = FALSE
  ]
  print(method_condition_summary)
  write.csv(
    method_condition_summary,
    results_path(
      state, paste0("method_condition_summary_", config$name, ".csv")
    ),
    row.names = FALSE,
    na = ""
  )

  state$simulation <- simulation
  state$main_coverage <- main_coverage
  state$method_condition_summary <- method_condition_summary
  state$figures$simulation <- plot_simulation_summary(
    state$project_root, simulation$aggregate, config, language = state$language
  )
  state$figures$example1 <- plot_example1_empirical(
    state$project_root, simulation$aggregate, config, language = state$language
  )
  state$figures$psi0_diagnostic <- plot_psi0_overlap_diagnostic(
    state$project_root, simulation$aggregate, config, language = state$language
  )
  state$figures$coverage <- plot_coverage_heatmap(
    state$project_root, simulation$aggregate, config, language = state$language
  )
  state$figures$interval_endpoints <- plot_interval_endpoints(
    state$project_root, simulation$aggregate, config, language = state$language
  )
  for (figure in c(
    "simulation", "example1", "psi0_diagnostic", "coverage",
    "interval_endpoints"
  )) {
    display_figure(state$figures[[figure]])
  }
  state
}

# --- Stage: benchmark baseline có oracle --------------------------------------

stage_baseline <- function(state) {
  if (is.null(state$simulation)) {
    stop("Stage `baseline` cần kết quả của stage `simulations`.")
  }
  config <- state$config
  baseline_benchmark <- run_baseline_benchmark(
    state$project_root, config, language = state$language
  )
  baseline_table <- baseline_benchmark$aggregate[
    order(
      baseline_benchmark$aggregate$example,
      baseline_benchmark$aggregate$learner,
      baseline_benchmark$aggregate$estimator
    ),
    c(
      "example", "rho", "learner", "estimator", "repetitions",
      "target_psi0", "mean_estimate", "relative_absolute_bias",
      "relative_rmse", "standard_deviation",
      "coverage", "coverage_ci_low", "coverage_ci_high",
      "mean_ci_width", "finite_fraction", "method_role"
    ),
    drop = FALSE
  ]
  names(baseline_table)[names(baseline_table) == "coverage"] <-
    "oracle_containment_rate"
  names(baseline_table)[names(baseline_table) == "coverage_ci_low"] <-
    "oracle_containment_ci_low"
  names(baseline_table)[names(baseline_table) == "coverage_ci_high"] <-
    "oracle_containment_ci_high"
  print(baseline_table)
  write.csv(
    baseline_table,
    results_path(state, paste0("baseline_", config$name, "_report_table.csv")),
    row.names = FALSE,
    na = ""
  )

  method_comparison <- assemble_method_comparison(
    state$simulation$aggregate, baseline_benchmark$aggregate
  )
  write.csv(
    method_comparison,
    results_path(state, paste0("method_comparison_", config$name, ".csv")),
    row.names = FALSE,
    na = ""
  )

  comparison_groups <- split(
    method_comparison,
    method_comparison$method_label
  )
  baseline_assessment <- do.call(rbind, lapply(
    comparison_groups,
    function(group) {
      relative_rmse <- group$relative_rmse[is.finite(group$relative_rmse)]
      relative_bias <- group$relative_absolute_bias[
        is.finite(group$relative_absolute_bias)
      ]
      data.frame(
        method = group$method_label[[1L]],
        evaluated_combinations = nrow(group),
        median_relative_rmse = if (length(relative_rmse)) {
          stats::median(relative_rmse)
        } else NA_real_,
        maximum_relative_rmse = if (length(relative_rmse)) {
          max(relative_rmse)
        } else NA_real_,
        median_relative_absolute_bias = if (length(relative_bias)) {
          stats::median(relative_bias)
        } else NA_real_,
        finite_fraction = safe_mean(group$finite_fraction),
        interpretation_scope = if (
          identical(group$method_label[[1L]], "Marginal PFI")
        ) {
          "phép nhiễu biên; không cùng estimand với tham số decorrelated"
        } else if (
          identical(group$method_label[[1L]], "LOCO")
        ) {
          "estimand phụ thuộc phân phối chung của X và Z"
        } else if (
          identical(group$method_label[[1L]], "Residual CPI")
        ) {
          "xấp xỉ CPI phụ thuộc mô hình X|Z"
        } else {
          "ước lượng trực tiếp tham số decorrelated"
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(baseline_assessment) <- NULL
  baseline_assessment <- baseline_assessment[
    order(baseline_assessment$median_relative_rmse),
    ,
    drop = FALSE
  ]
  print(baseline_assessment)
  write.csv(
    baseline_assessment,
    results_path(state, paste0("baseline_assessment_", config$name, ".csv")),
    row.names = FALSE,
    na = ""
  )

  state$baseline_benchmark <- baseline_benchmark
  state$method_comparison <- method_comparison
  state$baseline_assessment <- baseline_assessment
  state$figures$baseline <- plot_baseline_comparison(
    state$project_root, method_comparison, config, language = state$language
  )
  display_figure(state$figures$baseline)
  state
}

# --- Stage: ablation của số hạng hiệu chỉnh t-Cross ---------------------------

stage_ablation <- function(state) {
  config <- state$config
  ablation <- run_correction_ablation(
    state$project_root, config, language = state$language,
    example = 2L, learner = "linear"
  )
  ablation_table <- ablation$summary[
    order(ablation$summary$estimator, ablation$summary$correction_mode),
    c(
      "estimator", "correction_mode", "repetitions", "target_psi0",
      "coverage", "coverage_ci_low", "coverage_ci_high",
      "coverage_change_vs_none", "mean_estimate", "absolute_bias",
      "empirical_sd", "mean_width", "median_width", "width_p90",
      "mean_width_to_target", "width_inflation_vs_none",
      "finite_fraction"
    ),
    drop = FALSE
  ]
  print(ablation_table)
  state$ablation <- ablation
  state$ablation_table <- ablation_table
  state$figures$ablation <- plot_ablation_summary(
    state$project_root, ablation$summary, config, language = state$language
  )
  display_figure(state$figures$ablation)
  state
}

# --- Stage: phân tích hai bộ dữ liệu thực -------------------------------------

stage_real_data <- function(state) {
  real_outputs <- list()
  for (dataset_name in c("energy", "concrete")) {
    analysis <- run_real_dataset(
      state$project_root,
      dataset_name,
      state$config,
      language = state$language,
      data_path = state$data_paths[[dataset_name]]
    )
    real_outputs[[dataset_name]] <- analysis
    metadata <- analysis$metadata
    cat(
      "\n", dataset_name,
      ": n=", metadata$rows,
      ", response=", metadata$response,
      ", X=", metadata$x,
      ", SHA-256=", metadata$sha256,
      "\n",
      sep = ""
    )
    print(analysis$predictability)
    print(analysis$report_table)
    print(analysis$learner_sensitivity)

    importance_figure <- plot_real_results(
      state$project_root, dataset_name, analysis$diagnostics,
      language = state$language
    )
    predictability_figure <- plot_real_predictability(
      state$project_root,
      dataset_name,
      analysis$predictability,
      diagnostics = analysis$diagnostics,
      language = state$language
    )
    state$figures[[paste0("real_importance_", dataset_name)]] <-
      importance_figure
    state$figures[[paste0("real_predictability_", dataset_name)]] <-
      predictability_figure
    display_figure(importance_figure)
    display_figure(predictability_figure)
  }
  state$real_outputs <- real_outputs
  state
}

# --- Stage: metadata môi trường, runtime và manifest --------------------------

stage_reporting <- function(state) {
  config <- state$config
  environment_path <- write_environment_metadata(
    state$project_root, config, language = state$language
  )
  readme_path <- write_experiment_readme(
    state$project_root, config, state$data_paths
  )
  cat("Environment metadata:", environment_path, "\n")
  cat("Execution guide:", readme_path, "\n")

  if (!is.null(state$simulation)) {
    # Mỗi lần t-Cross xuất một dòng cho nhiều estimator. Dòng LOCO được dùng
    # làm đại diện để không cộng lặp cùng một thời gian tính toán.
    runtime_rows <- state$simulation$results[
      state$simulation$results$estimator == "psi_L",
      ,
      drop = FALSE
    ]
    runtime_summary <- aggregate(
      shared_runtime_seconds ~ example + learner,
      data = runtime_rows,
      FUN = sum
    )
    names(runtime_summary)[
      names(runtime_summary) == "shared_runtime_seconds"
    ] <- "wall_clock_proxy_seconds"
    print(runtime_summary)
    write.csv(
      runtime_summary,
      results_path(state, paste0("runtime_summary_", config$name, ".csv")),
      row.names = FALSE
    )
    state$runtime_summary <- runtime_summary
  }

  manifest_path <- write_output_manifest(
    state$project_root, language = state$language
  )
  cat("Output manifest:", manifest_path, "\n")
  state$environment_path <- environment_path
  state$readme_path <- readme_path
  state$manifest_path <- manifest_path
  state
}

# --- Stage: đối chiếu với yêu cầu thực nghiệm của đồ án -----------------------

check_status <- function(condition) {
  if (isTRUE(condition)) "Đạt" else "Chưa đạt"
}

stage_compliance <- function(state) {
  dataset_count <- length(state$real_outputs)
  baseline_count <- if (!is.null(state$baseline_benchmark)) {
    length(unique(state$baseline_benchmark$aggregate$estimator))
  } else 0L
  has_metrics <- !is.null(state$main_coverage) && all(c(
    "relative_rmse", "coverage", "coverage_ci_low", "coverage_ci_high",
    "mean_ci_width", "finite_fraction"
  ) %in% names(state$main_coverage))
  required_figures <- c(
    "simulation", "example1", "psi0_diagnostic", "coverage",
    "interval_endpoints", "baseline", "ablation"
  )
  has_figures <- all(required_figures %in% names(state$figures)) &&
    all(file.exists(unlist(state$figures[required_figures])))
  has_method_analysis <- !is.null(state$method_condition_summary) &&
    !is.null(state$baseline_assessment) &&
    dataset_count >= 2L &&
    all(vapply(
      state$real_outputs,
      function(item) {
        nrow(item$report_table) > 0L &&
          nrow(item$learner_sensitivity) > 0L
      },
      logical(1)
    ))
  has_ablation <- !is.null(state$ablation) &&
    min(state$ablation$summary$repetitions, na.rm = TRUE) >= 30L
  has_reproducibility <- !is.null(state$manifest_path) &&
    all(file.exists(c(
      state$environment_path, state$readme_path, state$manifest_path
    ))) && dataset_count >= 1L && all(vapply(
      state$real_outputs,
      function(item) {
        is.character(item$metadata$sha256) &&
          nzchar(item$metadata$sha256) &&
          !is.na(item$metadata$sha256)
      },
      logical(1)
    ))

  coursework_checklist <- data.frame(
    requirement = c(
      "Ít nhất 2 bộ dữ liệu chuẩn",
      "Ít nhất 2 baseline",
      "Độ đo phù hợp và khoảng bất định",
      "Bảng và biểu đồ phục vụ báo cáo",
      "Phân tích điều kiện phương pháp hoạt động ổn định hoặc bất ổn",
      "Ablation study",
      "Tính tái thực nghiệm"
    ),
    status = c(
      check_status(dataset_count >= 2L),
      check_status(baseline_count >= 2L),
      check_status(has_metrics),
      check_status(has_figures),
      check_status(has_method_analysis),
      check_status(has_ablation),
      check_status(has_reproducibility)
    ),
    evidence = c(
      paste0(dataset_count, " bộ dữ liệu thực đã được phân tích"),
      paste0(baseline_count, " baseline đã được đánh giá"),
      "Độ lệch, RMSE, coverage, khoảng nhị thức, độ rộng và tỷ lệ hữu hạn",
      "PNG 300 DPI, PDF vector và các bảng CSV",
      "Bảng độ nhạy theo learner, ESS, định danh và cờ chẩn đoán",
      paste0(
        if (!is.null(state$ablation)) {
          min(state$ablation$summary$repetitions)
        } else 0L,
        " lần lặp với cùng dữ liệu và fold giữa các chế độ"
      ),
      "Seed, checkpoint, package metadata, README, manifest và SHA-256 dữ liệu"
    ),
    stringsAsFactors = FALSE
  )
  print(coursework_checklist)
  write.csv(
    coursework_checklist,
    results_path(state, "coursework_checklist.csv"),
    row.names = FALSE
  )
  manifest_path <- write_output_manifest(
    state$project_root, language = state$language
  )
  cat("Updated output manifest:", manifest_path, "\n")
  state$coursework_checklist <- coursework_checklist
  state$manifest_path <- manifest_path
  state
}

# --- Điều phối ----------------------------------------------------------------

DVI_STAGES <- c(
  "targets", "simulations", "baseline", "ablation", "real", "report", "check"
)

run_pipeline <- function(state, stages, profile, n_threads) {
  unknown <- setdiff(stages, DVI_STAGES)
  if (length(unknown)) {
    stop(
      "Stage không hợp lệ: ", paste(unknown, collapse = ", "),
      ". Các stage hợp lệ: ", paste(DVI_STAGES, collapse = ", ")
    )
  }
  if ("baseline" %in% stages && !"simulations" %in% stages) {
    stop("Stage `baseline` phải được chạy cùng stage `simulations`.")
  }
  state <- stage_check_environment(state, profile, n_threads)
  runners <- list(
    targets = stage_targets,
    simulations = stage_simulations,
    baseline = stage_baseline,
    ablation = stage_ablation,
    real = stage_real_data,
    report = stage_reporting,
    check = stage_compliance
  )
  for (stage in DVI_STAGES[DVI_STAGES %in% stages]) {
    started <- proc.time()[[3L]]
    cat("\n==== Stage:", stage, "====\n")
    state <- runners[[stage]](state)
    cat(
      "==== Stage", stage, "hoàn tất sau",
      format_progress_time(proc.time()[[3L]] - started), "====\n"
    )
  }
  invisible(state)
}
