# Bộ ước lượng psi_L, psi_0, psi_1, psi_2, psi_3 và các baseline theo fold.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

prepare_psi0_cache <- function(
  x_train,
  x_test,
  z_train,
  z_test,
  mc_samples,
  seed,
  density_ratio_cap = NULL
) {
  started <- proc.time()[[3L]]
  kde_x <- fit_standardized_gaussian_kde(x_train)
  kde_z <- fit_standardized_gaussian_kde(z_train)
  kde_xz <- fit_standardized_gaussian_kde(cbind(x_train, z_train))
  x_star <- kde_x$sample(mc_samples, seed_from(seed, 1L))
  z_star <- kde_z$sample(mc_samples, seed_from(seed, 2L))
  log_ratio <- kde_x$logpdf(x_test) + kde_z$logpdf(z_test) -
    kde_xz$logpdf(cbind(x_test, z_test))
  if (max(log_ratio) > log(.Machine$double.xmax)) {
    stop("Density ratio overflow; overlap is too weak.")
  }
  raw_density_ratio <- exp(log_ratio)
  density_ratio <- raw_density_ratio
  clipped_fraction <- 0
  if (!is.null(density_ratio_cap)) {
    clipped_fraction <- mean(density_ratio > density_ratio_cap)
    density_ratio <- pmin(density_ratio, density_ratio_cap)
  }
  square_sum <- sum(density_ratio^2)
  ess <- if (square_sum > 0) sum(density_ratio)^2 / square_sum else 0
  list(
    x_star = x_star,
    z_star = z_star,
    raw_density_ratio = raw_density_ratio,
    density_ratio = density_ratio,
    density_ratio_mean = mean(density_ratio),
    density_ratio_median = stats::median(density_ratio),
    density_ratio_p95 = unname(stats::quantile(density_ratio, 0.95)),
    density_ratio_p99 = unname(stats::quantile(density_ratio, 0.99)),
    density_ratio_max = max(density_ratio),
    density_ratio_ess = ess,
    density_ratio_ess_fraction = ess / length(density_ratio),
    density_ratio_clipped_fraction = clipped_fraction,
    psi0_kde_seconds = proc.time()[[3L]] - started
  )
}

predict_cross_x <- function(model, basis, x_star, z_values, chunk_size) {
  x_design <- basis$transform(x_star)
  output <- matrix(NA_real_, nrow(z_values), nrow(x_star))
  starts <- seq.int(1L, nrow(z_values), by = chunk_size)
  for (start in starts) {
    stop_index <- min(nrow(z_values), start + chunk_size - 1L)
    z_block <- z_values[start:stop_index, , drop = FALSE]
    design <- cbind(
      x_design[rep(seq_len(nrow(x_design)), times = nrow(z_block)), , drop = FALSE],
      z_block[rep(seq_len(nrow(z_block)), each = nrow(x_design)), , drop = FALSE]
    )
    output[start:stop_index, ] <- matrix(
      model$predict(design),
      nrow = nrow(z_block),
      byrow = TRUE
    )
  }
  output
}

predict_cross_z <- function(model, x_design, z_star, chunk_size) {
  output <- matrix(NA_real_, nrow(x_design), nrow(z_star))
  starts <- seq.int(1L, nrow(x_design), by = chunk_size)
  for (start in starts) {
    stop_index <- min(nrow(x_design), start + chunk_size - 1L)
    x_block <- x_design[start:stop_index, , drop = FALSE]
    design <- cbind(
      x_block[rep(seq_len(nrow(x_block)), each = nrow(z_star)), , drop = FALSE],
      z_star[rep(seq_len(nrow(z_star)), times = nrow(x_block)), , drop = FALSE]
    )
    output[start:stop_index, ] <- matrix(
      model$predict(design),
      nrow = nrow(x_block),
      byrow = TRUE
    )
  }
  output
}

estimate_psi0_fold <- function(
  x_test_design,
  z_test,
  y_test,
  basis,
  mu_xz_model,
  mu_xz_observed,
  cache,
  chunk_size
) {
  started <- proc.time()[[3L]]
  mu_xstar_z <- predict_cross_x(
    mu_xz_model, basis, cache$x_star, z_test, chunk_size
  )
  mu0_z <- rowMeans(mu_xstar_z)
  first <- rowMeans((mu_xstar_z - mu0_z)^2)
  mu_x_zstar <- predict_cross_z(
    mu_xz_model, x_test_design, cache$z_star, chunk_size
  )
  mu0_zstar <- rowMeans(predict_cross_x(
    mu_xz_model, basis, cache$x_star, cache$z_star, chunk_size
  ))
  second <- rowMeans((mu_x_zstar - rep(mu0_zstar, each = nrow(mu_x_zstar)))^2)
  correction <- cache$density_ratio *
    (mu_xz_observed - mu0_z) *
    (y_test - mu_xz_observed)
  estimate <- 0.5 * mean(first + second) + mean(correction)
  diagnostics <- cache[!vapply(cache, is.matrix, logical(1))]
  diagnostics$raw_density_ratio <- NULL
  diagnostics$density_ratio <- NULL
  diagnostics$psi0_first_integral <- mean(first)
  diagnostics$psi0_second_integral <- mean(second)
  diagnostics$psi0_correction <- mean(correction)
  diagnostics$psi0_estimator_seconds <- proc.time()[[3L]] - started
  list(estimate = estimate, diagnostics = diagnostics)
}

estimate_fold <- function(
  data,
  fold_assignments,
  fold,
  learner,
  seed,
  config,
  psi0_cache = NULL,
  include_psi0 = TRUE,
  include_baselines = FALSE
) {
  total_started <- proc.time()[[3L]]
  test <- which(fold_assignments == fold)
  train <- which(fold_assignments != fold)
  x_raw <- as_2d(data$x)
  z_all <- as_2d(data$z)
  y_all <- as.numeric(data$y)
  x_raw_train <- x_raw[train, , drop = FALSE]
  x_raw_test <- x_raw[test, , drop = FALSE]
  z_train <- z_all[train, , drop = FALSE]
  z_test <- z_all[test, , drop = FALSE]
  y_train <- y_all[train]
  y_test <- y_all[test]
  basis <- fit_basis(x_raw_train, data$basis)
  x_train <- basis$train
  x_test <- basis$transform(x_raw_test)

  nuisance_started <- proc.time()[[3L]]
  mu_z_model <- fit_regressor(
    learner,
    z_train,
    y_train,
    seed_from(seed, 1L),
    config$forest_trees,
    config$n_threads
  )
  mu_xz_model <- fit_regressor(
    learner,
    cbind(x_train, z_train),
    y_train,
    seed_from(seed, 2L),
    config$forest_trees,
    config$n_threads
  )
  mu_z <- mu_z_model$predict(z_test)
  mu_xz <- mu_xz_model$predict(cbind(x_test, z_test))
  nu_models <- fit_multivariate_regression(
    learner,
    z_train,
    x_train,
    seed_from(seed, 5L),
    config$forest_trees,
    config$n_threads
  )
  nu <- predict_multivariate(nu_models, z_test)

  # Residual CPI và marginal PFI được xây dựng trên X gốc. Với basis đa thức,
  # việc hoán vị trực tiếp từng cột của basis có thể tạo tổ hợp không tương ứng
  # với bất kỳ giá trị X nào trong không gian đầu vào.
  if (isTRUE(include_baselines)) {
    if (identical(data$basis, "identity")) {
      nu_raw <- nu
    } else {
      nu_raw_models <- fit_multivariate_regression(
        learner,
        z_train,
        x_raw_train,
        seed_from(seed, 15L),
        config$forest_trees,
        config$n_threads
      )
      nu_raw <- predict_multivariate(nu_raw_models, z_test)
    }
  }
  nuisance_seconds <- proc.time()[[3L]] - nuisance_started

  estimator_started <- proc.time()[[3L]]
  psi_l <- mean((y_test - mu_z)^2 - (y_test - mu_xz)^2)
  dependence <- vapply(seq_len(ncol(z_train)), function(j) {
    sum(vapply(seq_len(ncol(x_train)), function(k) {
      safe_abs_correlation(x_train[, k], z_train[, j])
    }, numeric(1)))
  }, numeric(1))
  selected <- which(dependence <= config$psi1_threshold)
  if (length(selected) == ncol(z_train)) {
    mu_v <- mu_z
    mu_xv <- mu_xz
  } else {
    v_train <- z_train[, selected, drop = FALSE]
    v_test <- z_test[, selected, drop = FALSE]
    mu_v_model <- fit_regressor(
      learner,
      v_train,
      y_train,
      seed_from(seed, 3L),
      config$forest_trees,
      config$n_threads
    )
    mu_xv_model <- fit_regressor(
      learner,
      cbind(x_train, v_train),
      y_train,
      seed_from(seed, 4L),
      config$forest_trees,
      config$n_threads
    )
    mu_v <- mu_v_model$predict(v_test)
    mu_xv <- mu_xv_model$predict(cbind(x_test, v_test))
  }
  psi_1 <- mean((y_test - mu_v)^2 - (y_test - mu_xv)^2)

  residual_x <- x_test - nu
  residual_y <- y_test - mu_z
  m <- length(test)
  a2 <- crossprod(residual_x) / m
  b2 <- as.numeric(crossprod(residual_x, residual_y) / m)
  sigma_x <- population_covariance(x_test)
  largest_x_variance <- max(eigen(sigma_x, symmetric = TRUE, only.values = TRUE)$values)
  smallest_residual_variance <- min(eigen(a2, symmetric = TRUE, only.values = TRUE)$values)
  residual_variance_ratio <- if (largest_x_variance > 0) {
    smallest_residual_variance / largest_x_variance
  } else 0
  solution2 <- solve_moment(a2, b2, config$condition_limit)
  psi2_valid <- isTRUE(solution2$valid) &&
    residual_variance_ratio >= config$residual_variance_ratio_min
  if (psi2_valid) {
    beta <- solution2$coefficient
    residual2 <- residual_y - as.numeric(residual_x %*% beta)
    phi_beta <- (residual_x * residual2) %*% t(solve(a2))
    centered_x <- sweep(x_test, 2L, colMeans(x_test), "-")
    psi_2 <- mean(
      as.numeric(centered_x %*% beta)^2 +
        2 * as.numeric(phi_beta %*% (sigma_x %*% beta))
    )
  } else {
    psi_2 <- NA_real_
  }

  z_tilde <- cbind(1, z_test)
  residual_xz <- row_kronecker(residual_x, z_tilde)
  a3 <- crossprod(residual_xz) / m
  b3 <- as.numeric(crossprod(residual_xz, residual_y) / m)
  solution3 <- solve_moment(a3, b3, config$condition_limit)
  psi3_valid <- isTRUE(solution3$valid) &&
    residual_variance_ratio >= config$residual_variance_ratio_min
  if (psi3_valid) {
    theta <- solution3$coefficient
    omega <- kronecker(sigma_x, crossprod(z_tilde) / m)
    residual3 <- residual_y - as.numeric(residual_xz %*% theta)
    phi_theta <- (residual_xz * residual3) %*% t(solve(a3))
    psi_3 <- as.numeric(t(theta) %*% omega %*% theta) +
      mean(2 * as.numeric(phi_theta %*% (omega %*% theta)))
  } else {
    psi_3 <- NA_real_
  }

  estimates <- c(
    psi_L = psi_l,
    psi_1 = psi_1,
    psi_2 = psi_2,
    psi_3 = psi_3
  )
  diagnostics <- list(
    selected_z_columns = selected,
    dependence_scores = dependence,
    residual_variance_ratio = residual_variance_ratio,
    psi2_valid = psi2_valid,
    psi3_valid = psi3_valid,
    psi2_condition_number = solution2$condition_number,
    psi3_condition_number = solution3$condition_number,
    mu_z_backend = mu_z_model$label,
    mu_xz_backend = mu_xz_model$label
  )

  psi0_seconds <- 0
  if (isTRUE(include_psi0)) {
    psi0_result <- tryCatch({
      if (!is.null(psi0_cache$error)) stop(psi0_cache$error)
      if (is.null(psi0_cache)) {
        psi0_cache <- prepare_psi0_cache(
          x_raw_train,
          x_raw_test,
          z_train,
          z_test,
          config$psi0_mc_samples,
          seed_from(seed, 6L),
          config$density_ratio_cap
        )
      }
      estimate_psi0_fold(
        x_test,
        z_test,
        y_test,
        basis,
        mu_xz_model,
        mu_xz,
        psi0_cache,
        config$psi0_chunk_size
      )
    }, error = function(error) error)
    if (inherits(psi0_result, "error")) {
      estimates <- c(estimates, psi_0 = NA_real_)
      diagnostics$psi0_error <- conditionMessage(psi0_result)
    } else {
      estimates <- c(estimates, psi_0 = psi0_result$estimate)
      diagnostics <- c(diagnostics, psi0_result$diagnostics)
      psi0_seconds <- psi0_result$diagnostics$psi0_estimator_seconds
    }
  }

  if (isTRUE(include_baselines)) {
    set.seed(seed_from(seed, 19L))
    base_loss <- (y_test - mu_xz)^2
    residual_x_raw <- x_raw_test - nu_raw
    cpi <- pfi <- numeric(config$permutation_repetitions)
    for (permutation in seq_len(config$permutation_repetitions)) {
      order <- sample.int(m, m, replace = FALSE)
      x_conditional_raw <- nu_raw + residual_x_raw[order, , drop = FALSE]
      x_marginal_raw <- x_raw_test[order, , drop = FALSE]
      x_conditional_design <- basis$transform(x_conditional_raw)
      x_marginal_design <- basis$transform(x_marginal_raw)
      conditional_prediction <- mu_xz_model$predict(
        cbind(x_conditional_design, z_test)
      )
      marginal_prediction <- mu_xz_model$predict(
        cbind(x_marginal_design, z_test)
      )
      cpi[permutation] <- mean((y_test - conditional_prediction)^2 - base_loss)
      pfi[permutation] <- mean((y_test - marginal_prediction)^2 - base_loss)
    }
    estimates <- c(
      estimates,
      cpi_residual = mean(cpi),
      pfi = mean(pfi)
    )
  }

  requested_order <- c(
    "psi_L", "psi_0", "psi_1", "psi_2", "psi_3", "cpi_residual", "pfi"
  )
  estimates <- estimates[requested_order[requested_order %in% names(estimates)]]
  estimator_seconds <- proc.time()[[3L]] - estimator_started
  timings <- list(
    nuisance_fit_seconds = nuisance_seconds,
    psi0_seconds = psi0_seconds,
    estimator_postprocess_seconds = estimator_seconds,
    total_fold_seconds = proc.time()[[3L]] - total_started
  )
  list(estimates = estimates, diagnostics = diagnostics, timings = timings)
}
