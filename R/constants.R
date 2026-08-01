# Hằng số toàn cục và profile thực nghiệm.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

# Đọc lại từ biến môi trường đã được main.R (hoặc cell cấu hình của notebook)
# đặt trước khi nạp module. Giá trị này là mặc định cho `get_profile()`.
DVI_N_THREADS <- as.integer(Sys.getenv("DVI_N_THREADS", "4"))
if (!is.finite(DVI_N_THREADS) || DVI_N_THREADS < 1L) DVI_N_THREADS <- 4L
DVI_PROFILE <- tolower(Sys.getenv("DVI_PROFILE", "reduced"))

# "paper": se^2 = s^2/B + c^2/n với c = Var(Y)^2 (mục 3.6). "cross": se^2 = s^2/B.
PRIMARY_INTERVAL <- tolower(Sys.getenv("DVI_PRIMARY_INTERVAL", "paper"))
if (!PRIMARY_INTERVAL %in% c("paper", "cross")) {
  stop("DVI_PRIMARY_INTERVAL must be paper or cross.")
}

primary_interval_is_paper <- function() identical(PRIMARY_INTERVAL, "paper")

primary_interval_label <- function() {
  if (primary_interval_is_paper()) {
    "t-Cross theo bài báo"
  } else {
    "t-Cross không hiệu chỉnh"
  }
}

BASE_SEED <- 2026L
ESTIMATORS <- c("psi_L", "psi_0", "psi_1", "psi_2", "psi_3")
BASELINE_ESTIMATORS <- c("psi_L", "cpi_residual", "pfi")
ALL_ESTIMATORS <- c(ESTIMATORS, "cpi_residual", "pfi")
REPORT_DPI <- as.integer(Sys.getenv("DVI_REPORT_DPI", "300"))
if (!is.finite(REPORT_DPI) || REPORT_DPI < 120L) REPORT_DPI <- 300L
ESTIMATOR_COLORS <- c(
  psi_L = "#0072B2",
  psi_0 = "#D55E00",
  psi_1 = "#009E73",
  psi_2 = "#CC79A7",
  psi_3 = "#E69F00",
  cpi_residual = "#56B4E9",
  pfi = "#6A6A6A"
)
CORRECTION_COLORS <- c(
  none = "#0072B2",
  variance_scale = "#009E73",
  paper_literal = "#D55E00"
)

get_profile <- function(name = Sys.getenv("DVI_PROFILE", "smoke"), n_threads = DVI_N_THREADS) {
  name <- tolower(name)
  profiles <- list(
    smoke = list(
      n = 500L,
      repetitions = 2L,
      examples = 1:2,
      correlations = c(0, 0.6),
      forest_trees = 100L,
      psi0_mc_samples = 20L,
      permutation_repetitions = 5L,
      ablation_n = 500L,
      ablation_repetitions = 5L,
      baseline_repetitions = 2L,
      baseline_examples = c(1L, 2L),
      baseline_rho = 0.75
    ),
    reduced = list(
      n = 10000L,
      repetitions = 10L,
      examples = 1:5,
      correlations = c(0, 0.25, 0.49, 0.51, 0.75, 0.95),
      forest_trees = 300L,
      psi0_mc_samples = 30L,
      permutation_repetitions = 20L,
      ablation_n = 1000L,
      ablation_repetitions = 30L,
      baseline_repetitions = 10L,
      baseline_examples = c(1L, 2L, 5L),
      baseline_rho = 0.75
    ),
    paper = list(
      n = 10000L,
      repetitions = 100L,
      examples = 1:5,
      correlations = c(seq(0, 0.9, by = 0.1), 0.95),
      forest_trees = NULL,
      psi0_mc_samples = 100L,
      permutation_repetitions = 20L,
      ablation_n = 10000L,
      ablation_repetitions = 100L,
      baseline_repetitions = 100L,
      baseline_examples = c(1L, 2L, 5L),
      baseline_rho = 0.75
    )
  )
  if (!name %in% names(profiles)) {
    stop("DVI_PROFILE must be smoke, reduced, or paper.")
  }
  c(
    list(
      name = name,
      folds = 5L,
      n_threads = as.integer(n_threads),
      psi1_threshold = 0.5,
      density_ratio_cap = NULL,
      condition_limit = 1e12,
      residual_variance_ratio_min = 1e-8,
      psi0_chunk_size = 256L
    ),
    profiles[[name]]
  )
}

TABLE2_REFERENCE <- data.frame(
  example = rep(2:5, each = 15L),
  learner = rep(rep(c("linear", "additive", "forest"), each = 5L), 4L),
  estimator = rep(c("psi_L", "psi_0", "psi_1", "psi_2", "psi_3"), 12L),
  paper_coverage = c(
    1, .84, 1, 1, 1, 0, .79, 1, 1, .97, 0, .75, 1, .99, .30,
    0, 0, 0, 0, .99, 0, .88, 0, 0, .92, 0, 0, 0, 0, .91,
    1, .87, 1, 1, 1, .98, .20, .98, .98, .98, 0, .21, 0, .86, .85,
    0, .01, 0, 0, 0, 0, .83, 0, 0, .85, 0, .05, 0, 0, 1
  ),
  stringsAsFactors = FALSE
)
