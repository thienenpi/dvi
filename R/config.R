# Cấu hình môi trường chạy: thư mục gốc, đường dẫn dữ liệu và giới hạn luồng.
# Các hàm ở đây phải được gọi trước khi source những module còn lại vì
# constants.R và models.R đọc biến môi trường ngay tại thời điểm nạp.

# Thư mục gốc của đồ án (chứa `data/`, `docs/`, `src/`, `R/`). Ưu tiên
# DVI_PROJECT_ROOT; nếu không có thì dò từ `start` lên tối đa ba cấp để tìm thư
# mục có `data/`.
detect_project_root <- function(start = getwd()) {
  explicit <- Sys.getenv("DVI_PROJECT_ROOT", "")
  if (nzchar(explicit)) {
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }
  candidate <- normalizePath(start, winslash = "/", mustWork = TRUE)
  for (level in 0:3) {
    if (dir.exists(file.path(candidate, "data"))) return(candidate)
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  normalizePath(start, winslash = "/", mustWork = TRUE)
}

# Đường dẫn hai bộ dữ liệu thực. Mặc định là đường dẫn tương đối theo
# PROJECT_ROOT. Trên Kaggle hoặc cluster có layout khác, đặt
# DVI_ENERGY_DATA_PATH và DVI_CONCRETE_DATA_PATH bằng đường dẫn tuyệt đối.
resolve_data_paths <- function(energy = NULL, concrete = NULL) {
  if (is.null(energy) || !nzchar(energy)) {
    energy <- Sys.getenv(
      "DVI_ENERGY_DATA_PATH",
      file.path("data", "ENB2012_data.xlsx")
    )
  }
  if (is.null(concrete) || !nzchar(concrete)) {
    concrete <- Sys.getenv(
      "DVI_CONCRETE_DATA_PATH",
      file.path("data", "Concrete_Data.xls")
    )
  }
  list(energy = energy, concrete = concrete)
}

# Giới hạn luồng BLAS/OpenMP phải được đặt trước khi backend tính toán được nạp.
apply_thread_limits <- function(n_threads) {
  n_threads <- as.integer(n_threads)
  if (!is.finite(n_threads) || n_threads < 1L) {
    stop("Số luồng phải là số nguyên không nhỏ hơn một.")
  }
  Sys.setenv(
    DVI_N_THREADS = n_threads,
    OMP_NUM_THREADS = n_threads,
    MKL_NUM_THREADS = n_threads,
    OPENBLAS_NUM_THREADS = n_threads,
    NUMEXPR_NUM_THREADS = n_threads
  )
  n_threads
}

# Thứ tự nạp module là bắt buộc: constants.R và models.R phụ thuộc vào biến môi
# trường đã đặt, các module sau phụ thuộc vào hàm của module trước.
DVI_MODULES <- c(
  "constants.R",
  "utils.R",
  "models.R",
  "estimators.R",
  "inference.R",
  "simulation.R",
  "real_data.R",
  "plots.R",
  "reporting.R",
  "pipeline.R"
)

load_dvi_modules <- function(module_dir, modules = DVI_MODULES) {
  for (module in modules) {
    path <- file.path(module_dir, module)
    if (!file.exists(path)) stop("Không tìm thấy module: ", path)
    source(path, local = globalenv(), encoding = "UTF-8")
  }
  invisible(modules)
}
