# Tiện ích số học, chia fold, basis đa thức, KDE và thống kê an toàn.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

seed_from <- function(...) {
  values <- as.numeric(c(...))
  state <- as.numeric(BASE_SEED)
  for (value in values) {
    state <- (state * 48271 + value * 104729 + 1) %% 2147483647
  }
  as.integer(max(1, floor(state)))
}

as_2d <- function(value) {
  if (is.data.frame(value)) value <- data.matrix(value)
  if (is.null(dim(value))) value <- matrix(as.numeric(value), ncol = 1L)
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  if (length(value) && any(!is.finite(value))) stop("Input contains non-finite values.")
  value
}

population_covariance <- function(value) {
  value <- as_2d(value)
  centered <- sweep(value, 2L, colMeans(value), FUN = "-")
  crossprod(centered) / nrow(value)
}

ordered_sum <- function(values) {
  total <- 0
  for (value in as.numeric(values)) total <- total + value
  total
}

invert_upper_triangular <- function(matrix) {
  dimension <- nrow(matrix)
  inverse <- matrix(0, nrow = dimension, ncol = dimension)
  for (column in seq_len(dimension)) {
    solution <- numeric(dimension)
    for (row in dimension:1L) {
      rhs <- if (row == column) 1 else 0
      if (row < dimension) {
        for (index in seq.int(row + 1L, dimension)) {
          rhs <- rhs - matrix[row, index] * solution[index]
        }
      }
      solution[row] <- rhs / matrix[row, row]
    }
    inverse[, column] <- solution
  }
  inverse
}

make_folds <- function(n, folds, seed) {
  if (folds < 2L || n < folds) stop("Invalid fold count.")
  assignment <- rep(seq_len(folds), length.out = n)
  set.seed(seed)
  sample(assignment, size = n, replace = FALSE)
}

fit_basis <- function(raw, kind) {
  raw <- as_2d(raw)
  if (identical(kind, "identity")) {
    transform <- function(new_raw) as_2d(new_raw)
    return(list(
      kind = kind,
      train = raw,
      transform = transform,
      center = rep(0, ncol(raw)),
      scale = rep(1, ncol(raw))
    ))
  }
  if (!identical(kind, "orthogonal_polynomial") || ncol(raw) != 1L) {
    stop("Unsupported basis: ", kind)
  }
  center <- ordered_sum(raw[, 1L]) / nrow(raw)
  centered <- raw[, 1L] - center
  scale <- sqrt(ordered_sum(centered * centered) / (nrow(raw) - 1L))
  if (scale[1L] <= 1e-14) stop("Cannot fit polynomial basis to constant X.")
  q <- (raw[, 1L] - center[1L]) / scale[1L]
  vandermonde <- cbind(q, q^2, q^3)
  orthonormal <- matrix(0, nrow(vandermonde), ncol(vandermonde))
  r <- matrix(0, nrow = 3L, ncol = 3L)
  for (column in seq_len(3L)) {
    vector <- vandermonde[, column]
    if (column > 1L) {
      for (previous in seq_len(column - 1L)) {
        coefficient <- ordered_sum(orthonormal[, previous] * vector)
        r[previous, column] <- coefficient
        vector <- vector - coefficient * orthonormal[, previous]
      }
    }
    norm <- sqrt(ordered_sum(vector * vector))
    if (norm <= 1e-14) stop("Polynomial basis is rank deficient.")
    r[column, column] <- norm
    orthonormal[, column] <- vector / norm
  }
  transform_matrix <- invert_upper_triangular(r)
  transform <- function(new_raw) {
    new_raw <- as_2d(new_raw)
    q_new <- (new_raw[, 1L] - center[1L]) / scale[1L]
    cbind(q_new, q_new^2, q_new^3) %*% transform_matrix * sqrt(nrow(raw))
  }
  list(
    kind = kind,
    train = transform(raw),
    transform = transform,
    center = center,
    scale = scale,
    transform_matrix = transform_matrix
  )
}

fit_standardized_gaussian_kde <- function(training) {
  training <- as_2d(training)
  location <- colMeans(training)
  scale <- apply(training, 2L, stats::sd)
  scale[scale < 1e-12] <- 1
  standardized <- sweep(sweep(training, 2L, location, "-"), 2L, scale, "/")
  dimension <- ncol(standardized)
  bandwidth <- nrow(training)^(-1 / (dimension + 4))
  kernel_covariance <- bandwidth^2 * population_covariance(standardized) +
    diag(1e-8, dimension)
  cholesky <- chol(kernel_covariance)
  log_determinant <- 2 * sum(log(diag(cholesky)))
  log_jacobian <- sum(log(scale))

  logpdf <- function(query, query_chunk = 128L) {
    query <- as_2d(query)
    query <- sweep(sweep(query, 2L, location, "-"), 2L, scale, "/")
    constant <- -0.5 * (dimension * log(2 * pi) + log_determinant) -
      log(nrow(standardized)) - log_jacobian
    output <- numeric(nrow(query))
    starts <- seq.int(1L, nrow(query), by = query_chunk)
    for (start in starts) {
      stop_index <- min(nrow(query), start + query_chunk - 1L)
      for (i in start:stop_index) {
        differences <- sweep(standardized, 2L, query[i, ], FUN = "-")
        whitened <- forwardsolve(t(cholesky), t(differences))
        log_kernel <- -0.5 * colSums(whitened^2)
        maximum <- max(log_kernel)
        output[i] <- maximum + log(sum(exp(log_kernel - maximum))) + constant
      }
    }
    output
  }

  sample_kde <- function(count, seed) {
    set.seed(seed)
    centers <- standardized[sample.int(nrow(standardized), count, replace = TRUE), , drop = FALSE]
    noise <- matrix(stats::rnorm(count * dimension), nrow = count) %*% cholesky
    sweep(sweep(centers + noise, 2L, scale, "*"), 2L, location, "+")
  }

  list(
    logpdf = logpdf,
    sample = sample_kde,
    location = location,
    scale = scale,
    kernel_covariance = kernel_covariance
  )
}

safe_abs_correlation <- function(x, z) {
  if (stats::sd(x) < 1e-14 || stats::sd(z) < 1e-14) return(0)
  abs(stats::cor(x, z))
}

row_kronecker <- function(x, z) {
  x <- as_2d(x)
  z <- as_2d(z)
  output <- matrix(0, nrow(x), ncol(x) * ncol(z))
  for (j in seq_len(ncol(x))) {
    columns <- ((j - 1L) * ncol(z) + 1L):(j * ncol(z))
    output[, columns] <- z * x[, j]
  }
  output
}

solve_moment <- function(matrix, vector, condition_limit) {
  dimension <- nrow(matrix)
  rank <- qr(matrix)$rank
  condition <- tryCatch(kappa(matrix, exact = TRUE), error = function(e) Inf)
  valid <- rank == dimension && is.finite(condition) && condition <= condition_limit
  list(
    coefficient = if (valid) as.numeric(solve(matrix, vector)) else NULL,
    rank = rank,
    dimension = dimension,
    condition_number = condition,
    valid = valid
  )
}

safe_mean <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values)) mean(values) else NA_real_
}

safe_sd <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) >= 2L) stats::sd(values) else NA_real_
}

safe_median <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values)) stats::median(values) else NA_real_
}

safe_quantile <- function(values, probability) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values)) {
    as.numeric(stats::quantile(values, probability, names = FALSE))
  } else NA_real_
}

exact_binomial_interval <- function(successes, trials, level = 0.95) {
  if (!is.finite(trials) || trials <= 0L) return(c(NA_real_, NA_real_))
  interval <- stats::binom.test(
    as.integer(successes), as.integer(trials), conf.level = level
  )$conf.int
  as.numeric(interval)
}

format_progress_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("--:--:--")
  seconds <- as.integer(round(seconds))
  sprintf(
    "%02d:%02d:%02d",
    seconds %/% 3600L,
    (seconds %% 3600L) %/% 60L,
    seconds %% 60L
  )
}

show_progress <- function(
  label,
  completed,
  total,
  started,
  detail = "",
  initial_completed = 0L
) {
  elapsed <- proc.time()[["elapsed"]] - started
  fraction <- if (total > 0L) completed / total else 1
  measured <- completed - initial_completed
  eta <- if (measured > 0L) elapsed * (total - completed) / measured else Inf
  if (completed >= total) eta <- 0
  width <- 24L
  filled <- min(width, floor(width * fraction))
  bar <- paste0(
    strrep("=", filled),
    if (filled < width) ">" else "",
    strrep(" ", max(0L, width - filled - as.integer(filled < width)))
  )
  cat(sprintf(
    "\r%s [%s] %3d%% | elapsed %s | ETA %s | %s",
    label,
    bar,
    round(100 * fraction),
    format_progress_time(elapsed),
    format_progress_time(eta),
    detail
  ))
  flush.console()
  if (completed >= total) cat("\n")
}
