#!/usr/bin/env Rscript
# Cài các package R cần thiết cho thực nghiệm.
#
#   Rscript scripts/install_packages.R              # cài phần bắt buộc
#   Rscript scripts/install_packages.R --notebook   # cài thêm phần cho notebook
#
# `jsonlite`, `readxl`, `mgcv`, `grf` là bắt buộc. `IRdisplay` và `IRkernel` chỉ
# cần khi chạy src/dvi-r-reproduce.ipynb.

argv <- commandArgs(trailingOnly = TRUE)
include_notebook <- "--notebook" %in% argv

required <- c("jsonlite", "readxl", "mgcv", "grf")
if (include_notebook) required <- c(required, "IRdisplay", "IRkernel")

repository <- Sys.getenv("DVI_CRAN_REPOSITORY", "https://cloud.r-project.org")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  cat("Cài đặt:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = repository)
} else {
  cat("Tất cả package đã có sẵn.\n")
}

cat("\n", R.version.string, "\n", sep = "")
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
for (package in required) {
  cat(
    sprintf(
      "%-10s %s\n",
      package,
      if (status[[package]]) {
        as.character(utils::packageVersion(package))
      } else "CHƯA CÀI ĐƯỢC"
    )
  )
}
if (!all(status)) {
  quit(save = "no", status = 1L)
}
cat("\nMôi trường đã sẵn sàng. Chạy thử: Rscript main.R --profile smoke\n")
