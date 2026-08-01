# Metadata môi trường, hướng dẫn thực thi và manifest đầu ra.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

write_environment_metadata <- function(
  project_root, config, language = "R"
) {
  metadata <- c(
    list(
      language = "R",
      profile = config$name,
      dvi_n_threads = config$n_threads,
      base_seed = BASE_SEED,
      platform = R.version$platform,
      operating_system = paste(
        Sys.info()[c("sysname", "release", "machine")],
        collapse = " | "
      ),
      report_dpi = REPORT_DPI,
      primary_interval = if (primary_interval_is_paper()) {
        "t_cross_paper_c_equals_variance_y_squared"
      } else "t_cross_uncorrected",
      sensitivity_intervals = c(
        "t_cross_uncorrected",
        "c_equals_variance_y",
        "c_equals_variance_y_squared"
      ),
      peak_memory = "not recorded by this notebook",
      configuration = config
    ),
    check_reproduction_backends()
  )
  path <- file.path(
    project_root, "results", language,
    paste0("environment_", config$name, ".json")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(metadata, path, pretty = TRUE, auto_unbox = TRUE)
  path
}

write_experiment_readme <- function(project_root, config, data_paths) {
  loaded_metadata <- lapply(names(data_paths), function(dataset) {
    loaded <- load_real_dataset(
      project_root, dataset, data_path = data_paths[[dataset]]
    )
    loaded$metadata
  })
  names(loaded_metadata) <- names(data_paths)
  lines <- c(
    "# Hướng dẫn thực thi thực nghiệm DVI",
    "",
    "Repository tái thực nghiệm một số thí nghiệm trong bài Decorrelated Variable Importance.",
    "Profile mặc định là cấu hình giảm tải và phải được mô tả đúng phạm vi này trong báo cáo.",
    "",
    "## Package R bắt buộc",
    "",
    "`jsonlite`, `readxl`, `mgcv` và `grf`. `IRdisplay` chỉ cần khi chạy notebook.",
    "",
    "## Quy trình thực thi",
    "",
    "1. Đặt hai tệp dữ liệu vào `data/`, hoặc khai báo `DVI_ENERGY_DATA_PATH` và `DVI_CONCRETE_DATA_PATH`.",
    "2. Chọn profile `smoke`, `reduced` hoặc `paper`.",
    "3. Chạy `Rscript main.R --profile <profile> --threads <n>`; trên cluster dùng script trong `jobs/`.",
    "   Hoặc chạy `src/dvi-r-reproduce.ipynb` tuần tự từ đầu đến cuối; hai đường chạy dùng chung module trong `R/`.",
    "4. Kiểm tra bảng CSV chẩn đoán trước khi sử dụng hình trong báo cáo.",
    "",
    "## Chạy từng phần",
    "",
    "`--stages` nhận danh sách phân tách bằng dấu phẩy trong",
    "`targets,simulations,baseline,ablation,real,report,check`.",
    "Stage `baseline` phải chạy cùng `simulations` vì bảng so sánh cần kết quả của mô phỏng chính.",
    "Mô phỏng chính ghi kết quả tăng dần nên lần chạy bị ngắt có thể tiếp tục bằng đúng lệnh cũ.",
    "",
    "## Cấu hình đang sử dụng",
    "",
    paste0("- Profile: `", config$name, "`"),
    paste0("- Số lần lặp mô phỏng chính: ", config$repetitions),
    paste0("- Số lần lặp benchmark baseline: ", config$baseline_repetitions),
    paste0("- Số lần lặp ablation: ", config$ablation_repetitions),
    paste0("- Số fold cross-fitting: ", config$folds),
    "",
    "## Quy ước suy luận",
    "",
    paste0("- Khoảng trình bày chính: ", primary_interval_label(), "."),
    "- Mặc định là khoảng của mục 3.6: `se^2 = s^2/B + c^2/n` với `c = Var(Y)^2`; đây là khoảng so trực tiếp với Table 2.",
    "- Khoảng không hiệu chỉnh và khoảng với `c = Var(Y)` được giữ lại cho phân tích độ nhạy.",
    "- Độ bao phủ phải được đọc cùng khoảng nhị thức chính xác và độ rộng khoảng.",
    "",
    "## Tính toàn vẹn dữ liệu",
    "",
    paste0(
      "- Energy Efficiency: ", loaded_metadata$energy$rows, " mẫu; SHA-256 `",
      loaded_metadata$energy$sha256, "`"
    ),
    paste0(
      "- Concrete Strength: ", loaded_metadata$concrete$rows, " mẫu; SHA-256 `",
      loaded_metadata$concrete$sha256, "`"
    ),
    "",
    "## Thư mục đầu ra",
    "",
    "- `results/R/`: kết quả theo fold, bảng tổng hợp, chẩn đoán và metadata.",
    "- `figures/R/`: hình PNG 300 DPI và PDF vector.",
    "- `checkpoints/R/`: dữ liệu mô phỏng, fold assignment và cache tỷ số mật độ.",
    "",
    "## Giới hạn diễn giải",
    "",
    "Coverage với dưới 30 lần lặp chỉ mang tính mô tả. Marginal PFI không cùng estimand với tham số decorrelated. Kết quả dữ liệu thực không có oracle và phải được đọc cùng chẩn đoán định danh, overlap và độ nhạy theo learner."
  )
  path <- file.path(project_root, "README_EXPERIMENTS.md")
  writeLines(lines, path, useBytes = TRUE)
  path
}

write_output_manifest <- function(project_root, language = "R") {
  roots <- c(
    file.path(project_root, "results", language),
    file.path(project_root, "figures", language)
  )
  files <- unlist(lapply(roots, function(root) {
    if (dir.exists(root)) list.files(root, recursive = TRUE, full.names = TRUE) else character(0)
  }))
  files <- files[file.info(files)$isdir %in% FALSE]
  manifest <- data.frame(
    relative_path = substring(files, nchar(project_root) + 2L),
    size_bytes = file.info(files)$size,
    modified_utc = format(file.info(files)$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
  path <- file.path(project_root, "results", language, "output_manifest.csv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(manifest, path, row.names = FALSE)
  path
}
