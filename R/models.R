# Mô hình nuisance: linear, mgcv::gam, grf::regression_forest.
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

MODEL_N_THREADS <- as.integer(Sys.getenv("DVI_N_THREADS", "4"))
if (!is.finite(MODEL_N_THREADS) || MODEL_N_THREADS < 1L) {
  stop("DVI_N_THREADS must be at least one.")
}

fit_additive_gam <- function(x_train, y_train) {
  x_train <- as_2d(x_train)
  frame <- as.data.frame(x_train)
  names(frame) <- paste0("x", seq_len(ncol(x_train)))
  frame$y <- as.numeric(y_train)
  terms <- vapply(seq_len(ncol(x_train)), function(j) {
    name <- paste0("x", j)
    if (length(unique(x_train[, j])) >= 10L) paste0("s(", name, ")") else name
  }, character(1))
  formula <- stats::as.formula(paste("y ~", paste(terms, collapse = " + ")))
  model <- mgcv::gam(formula, data = frame)
  list(
    label = "mgcv_gam",
    predict = function(new_x) {
      new_frame <- as.data.frame(as_2d(new_x))
      names(new_frame) <- paste0("x", seq_len(ncol(new_frame)))
      as.numeric(stats::predict(model, newdata = new_frame, type = "response"))
    }
  )
}

fit_regressor <- function(
  learner,
  x_train,
  y_train,
  seed,
  forest_trees,
  n_threads = MODEL_N_THREADS
) {
  x_train <- as_2d(x_train)
  y_train <- as.numeric(y_train)
  if (nrow(x_train) != length(y_train)) {
    stop("Response and design matrix have different row counts.")
  }

  active <- if (ncol(x_train)) {
    apply(x_train, 2L, stats::sd) > 1e-14
  } else {
    logical(0)
  }
  if (!length(active) || !any(active)) {
    mean_value <- mean(y_train)
    return(list(
      label = "constant_mean",
      predict = function(new_x) rep(mean_value, nrow(as_2d(new_x)))
    ))
  }
  active_index <- which(active)
  x_active <- x_train[, active_index, drop = FALSE]

  if (identical(learner, "linear")) {
    fit <- stats::lm.fit(cbind(1, x_active), y_train)
    coefficient <- fit$coefficients
    coefficient[is.na(coefficient)] <- 0
    raw_predict <- function(new_x) {
      as.numeric(cbind(1, as_2d(new_x)) %*% coefficient)
    }
    label <- "linear_ols"
  } else if (identical(learner, "additive")) {
    fitted <- fit_additive_gam(x_active, y_train)
    raw_predict <- fitted$predict
    label <- fitted$label
  } else if (identical(learner, "forest")) {
    set.seed(as.integer(seed))
    arguments <- list(
      X = x_active,
      Y = y_train,
      num.threads = as.integer(n_threads),
      seed = as.integer(seed)
    )
    if (!is.null(forest_trees) && is.finite(forest_trees)) {
      arguments$num.trees <- as.integer(forest_trees)
    }
    model <- do.call(grf::regression_forest, arguments)
    raw_predict <- function(new_x) {
      as.numeric(predict(
        model,
        newdata = as_2d(new_x),
        num.threads = as.integer(n_threads)
      )$predictions)
    }
    label <- "grf_regression_forest"
  } else {
    stop("Unsupported learner: ", learner)
  }

  list(
    label = label,
    predict = function(new_x) {
      new_x <- as_2d(new_x)
      raw_predict(new_x[, active_index, drop = FALSE])
    }
  )
}

fit_multivariate_regression <- function(
  learner,
  x_train,
  y_train,
  seed,
  forest_trees,
  n_threads = MODEL_N_THREADS
) {
  y_train <- as_2d(y_train)
  lapply(seq_len(ncol(y_train)), function(column) {
    fit_regressor(
      learner,
      x_train,
      y_train[, column],
      seed + 1009L * column,
      forest_trees,
      n_threads
    )
  })
}

predict_multivariate <- function(models, new_x) {
  do.call(cbind, lapply(models, function(model) model$predict(new_x)))
}

check_reproduction_backends <- function() {
  required <- c("mgcv", "grf")
  available <- vapply(
    required,
    function(package) requireNamespace(package, quietly = TRUE),
    logical(1)
  )
  missing <- required[!available]
  if (length(missing)) {
    stop("Paper reproduction requires R packages: ", paste(missing, collapse = ", "))
  }
  list(
    r = R.version.string,
    mgcv = as.character(utils::packageVersion("mgcv")),
    grf = as.character(utils::packageVersion("grf"))
  )
}
