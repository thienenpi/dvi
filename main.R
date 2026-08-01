#!/usr/bin/env Rscript
# Điểm vào của quy trình tái thực nghiệm Decorrelated Variable Importance.
#
# Ví dụ:
#   Rscript main.R --profile reduced --threads 4
#   Rscript main.R --profile paper --threads 16 --stages simulations,baseline
#   Rscript main.R --profile smoke --stages real,report,check
#
# Toàn bộ mã phương pháp nằm trong R/; tệp này chỉ đọc tham số, nạp module theo
# đúng thứ tự phụ thuộc và gọi các stage.

parse_arguments <- function(argv) {
  options <- list(
    profile = Sys.getenv("DVI_PROFILE", "reduced"),
    threads = Sys.getenv("DVI_N_THREADS", "4"),
    dpi = Sys.getenv("DVI_REPORT_DPI", "300"),
    interval = Sys.getenv("DVI_PRIMARY_INTERVAL", "paper"),
    project_root = Sys.getenv("DVI_PROJECT_ROOT", ""),
    energy_data = "",
    concrete_data = "",
    stages = paste(
      c("targets", "simulations", "baseline", "ablation", "real",
        "report", "check"),
      collapse = ","
    ),
    language = "R"
  )
  known <- names(options)
  index <- 1L
  while (index <= length(argv)) {
    token <- argv[[index]]
    if (identical(token, "--help") || identical(token, "-h")) {
      cat(usage_text())
      quit(save = "no", status = 0L)
    }
    if (!startsWith(token, "--")) {
      stop("Tham số không hợp lệ: ", token, "\n", usage_text())
    }
    body <- substring(token, 3L)
    if (grepl("=", body, fixed = TRUE)) {
      key <- sub("=.*$", "", body)
      value <- sub("^[^=]*=", "", body)
      index <- index + 1L
    } else {
      key <- body
      if (index + 1L > length(argv)) stop("Thiếu giá trị cho --", key)
      value <- argv[[index + 1L]]
      index <- index + 2L
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    if (!key %in% known) {
      stop("Tham số không được hỗ trợ: --", key, "\n", usage_text())
    }
    options[[key]] <- value
  }
  options$stages <- trimws(strsplit(options$stages, ",", fixed = TRUE)[[1L]])
  options$stages <- options$stages[nzchar(options$stages)]
  options
}

usage_text <- function() {
  paste(
    "Cách dùng: Rscript main.R [tham số]",
    "",
    "  --profile        smoke | reduced | paper (mặc định: reduced)",
    "  --threads        số luồng BLAS/OpenMP (mặc định: 4)",
    "  --dpi            độ phân giải PNG cho hình báo cáo (mặc định: 300)",
    "  --interval       paper | cross: khoảng tin cậy dùng làm quy ước chính",
    "                   (mặc định: paper, tức se^2 = s^2/B + c^2/n với c = Var(Y)^2)",
    "  --project-root   thư mục gốc chứa data/, results/, figures/",
    "  --energy-data    đường dẫn ENB2012_data.xlsx",
    "  --concrete-data  đường dẫn Concrete_Data.xls",
    "  --stages         danh sách stage, phân tách bằng dấu phẩy:",
    "                   targets,simulations,baseline,ablation,real,report,check",
    "  --language       nhãn thư mục đầu ra (mặc định: R)",
    "",
    sep = "\n"
  )
}

script_directory <- function() {
  argv <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", argv, value = TRUE)
  if (length(file_argument)) {
    return(dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]))))
  }
  normalizePath(getwd())
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  options <- parse_arguments(argv)
  root_directory <- script_directory()

  # Biến môi trường phải được đặt trước khi nạp module vì constants.R và
  # models.R đọc chúng ngay tại thời điểm source.
  Sys.setenv(DVI_PROFILE = tolower(options$profile))
  Sys.setenv(DVI_REPORT_DPI = options$dpi)
  Sys.setenv(DVI_PRIMARY_INTERVAL = tolower(options$interval))
  if (nzchar(options$project_root)) {
    Sys.setenv(DVI_PROJECT_ROOT = options$project_root)
  }

  source(file.path(root_directory, "R", "config.R"), encoding = "UTF-8")
  n_threads <- apply_thread_limits(options$threads)
  project_root <- detect_project_root(root_directory)
  data_paths <- resolve_data_paths(options$energy_data, options$concrete_data)
  load_dvi_modules(file.path(root_directory, "R"))

  cat(
    "Profile:", tolower(options$profile), "| threads:", n_threads,
    "| interval:", tolower(options$interval), "\n"
  )
  cat("Stages:", paste(options$stages, collapse = ", "), "\n")
  started <- proc.time()[[3L]]
  state <- new_pipeline_state(project_root, data_paths, options$language)
  state <- run_pipeline(
    state, options$stages, tolower(options$profile), n_threads
  )
  cat(
    "\nTổng thời gian:",
    format_progress_time(proc.time()[[3L]] - started), "\n"
  )
  invisible(state)
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
