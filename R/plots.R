# Hình cho báo cáo (PNG 300 DPI + PDF vector).
# Tách từ notebook src/dvi-r-reproduce.ipynb; nội dung hàm giữ nguyên.

save_report_plot <- function(path, width, height, draw) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  png_arguments <- list(
    filename = path,
    width = width * REPORT_DPI,
    height = height * REPORT_DPI,
    res = REPORT_DPI,
    pointsize = 11,
    family = "sans"
  )
  if (isTRUE(capabilities("cairo"))) {
    png_arguments$type <- "cairo"
  }
  do.call(grDevices::png, png_arguments)
  draw()
  grDevices::dev.off()

  pdf_path <- sub("\\.png$", ".pdf", path)
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      pdf_path, width = width, height = height,
      family = "sans", pointsize = 11
    )
  } else {
    grDevices::pdf(
      pdf_path, width = width, height = height,
      family = "sans", pointsize = 11, useDingbats = FALSE
    )
  }
  draw()
  grDevices::dev.off()
  path
}

estimator_math_labels <- function(estimator_names, include_method = FALSE) {
  if (isTRUE(include_method)) {
    label_map <- list(
      psi_L = bquote("LOCO (" * psi[L] * ")"),
      psi_0 = bquote("DVI (" * psi[0] * ")"),
      psi_1 = bquote(psi[1]),
      psi_2 = bquote(psi[2]),
      psi_3 = bquote(psi[3]),
      cpi_residual = "Residual CPI",
      pfi = "Marginal PFI"
    )
  } else {
    label_map <- list(
      psi_L = bquote(psi[L]),
      psi_0 = bquote(psi[0]),
      psi_1 = bquote(psi[1]),
      psi_2 = bquote(psi[2]),
      psi_3 = bquote(psi[3])
    )
  }
  as.expression(unname(label_map[estimator_names]))
}

scenario_axis_labels <- function(indices, compact = FALSE) {
  if (isTRUE(compact)) {
    paste0("KB ", indices)
  } else {
    paste("Kịch bản mô phỏng", indices)
  }
}

plot_analytic_targets <- function(project_root, language = "R") {
  path <- file.path(project_root, "figures", language, "example1_analytic_targets.png")
  rho <- seq(0, 0.95, length.out = 200L)
  psi_l <- vapply(
    rho, function(value) analytical_targets(1L, value)$psi_L, numeric(1)
  )
  draw <- function() {
    par(
      mar = c(4.5, 4.5, 2.8, 1),
      family = "sans",
      bty = "l",
      las = 1
    )
    plot(
      rho, psi_l, type = "n",
      xlab = expression(paste("Tương quan ", rho)), ylab = "Giá trị mục tiêu",
      main = "Kịch bản mô phỏng 1: giá trị mục tiêu giải tích",
      xlim = c(-0.02, 0.97), ylim = c(-0.05, 4.2),
      xaxs = "i", yaxs = "i", font.main = 1
    )
    grid(col = "#E5E5E5", lwd = 0.7)
    lines(
      rho, rep(4, length(rho)), lwd = 2.2, lty = 4,
      col = ESTIMATOR_COLORS[["psi_0"]]
    )
    lines(rho, psi_l, lwd = 2.2, col = ESTIMATOR_COLORS[["psi_L"]])
    left <- rho <= 0.5
    right <- rho > 0.5
    lines(
      rho[left], psi_l[left], lwd = 2.2, lty = 2,
      col = ESTIMATOR_COLORS[["psi_1"]]
    )
    lines(
      rho[right], rep(4, sum(right)), lwd = 2.2, lty = 2,
      col = ESTIMATOR_COLORS[["psi_1"]]
    )
    abline(v = 0.5, col = "#777777", lty = 3, lwd = 1)
    text(0.53, 0.45, labels = expression(t == 0.5), col = "#555555", cex = 0.85)
    legend(
      "bottomleft",
      inset = c(0.02, 0.03),
      legend = c(
        expression(psi[L]),
        expression(psi[1]),
        expression(paste(psi[0], " = ", psi[2], " = ", psi[3]))
      ),
      col = c(
        ESTIMATOR_COLORS[["psi_L"]],
        ESTIMATOR_COLORS[["psi_1"]],
        ESTIMATOR_COLORS[["psi_0"]]
      ),
      lty = c(1, 2, 4),
      lwd = 2.2,
      seg.len = 2.8,
      bty = "n"
    )
  }
  save_report_plot(path, 7.2, 4.2, draw)
}

plot_simulation_summary <- function(
  project_root, aggregate, config, language = "R"
) {
  path <- file.path(
    project_root, "figures", language,
    paste0("simulation_", config$name, ".png")
  )
  subset <- aggregate[aggregate$example == 2L, , drop = FALSE]
  learners <- c("linear", "additive", "forest")
  offsets <- seq(-0.24, 0.24, length.out = length(ESTIMATORS))
  draw <- function() {
    par(
      mfrow = c(1, 2), mar = c(5, 4.5, 3, 1),
      oma = c(2.8, 0, 2.6, 7.5), family = "sans",
      bty = "l", las = 1
    )
    estimate_limits <- range(
      c(
        subset$mean_estimate - subset$standard_deviation,
        subset$mean_estimate + subset$standard_deviation,
        60
      ),
      na.rm = TRUE
    )
    plot(
      NA, xlim = c(0.6, 3.4), ylim = estimate_limits,
      xaxt = "n", xlab = "Mô hình nuisance",
      ylab = "Ước lượng trung bình và độ lệch chuẩn thực nghiệm",
      main = "Kịch bản mô phỏng 2: ước lượng điểm"
    )
    axis(1, at = 1:3, labels = tools::toTitleCase(learners))
    grid(col = "#E2E2E2", lwd = 0.7)
    abline(h = 60, lty = 2, col = "#555555")
    for (index in seq_along(ESTIMATORS)) {
      estimator <- ESTIMATORS[[index]]
      part <- subset[subset$estimator == estimator, , drop = FALSE]
      part <- part[match(learners, part$learner), , drop = FALSE]
      positions <- 1:3 + offsets[[index]]
      spread <- ifelse(
        is.finite(part$standard_deviation),
        part$standard_deviation,
        0
      )
      arrows(
        positions, part$mean_estimate - spread,
        positions, part$mean_estimate + spread,
        angle = 90, code = 3, length = 0.04,
        col = ESTIMATOR_COLORS[[estimator]]
      )
      points(
        positions, part$mean_estimate, pch = 16,
        col = ESTIMATOR_COLORS[[estimator]]
      )
    }

    plot(
      NA, xlim = c(0.6, 3.4), ylim = c(0, 1),
      xaxt = "n", xlab = "Mô hình nuisance",
      ylab = bquote("Độ bao phủ đối với" ~ psi[0]),
      main = "Độ bao phủ t-Cross không hiệu chỉnh"
    )
    axis(1, at = 1:3, labels = tools::toTitleCase(learners))
    grid(col = "#E2E2E2", lwd = 0.7)
    abline(h = 0.95, lty = 2, col = "#555555")
    for (index in seq_along(ESTIMATORS)) {
      estimator <- ESTIMATORS[[index]]
      part <- subset[subset$estimator == estimator, , drop = FALSE]
      part <- part[match(learners, part$learner), , drop = FALSE]
      positions <- 1:3 + offsets[[index]]
      arrows(
        positions, part$coverage_ci_low,
        positions, part$coverage_ci_high,
        angle = 90, code = 3, length = 0.035,
        col = ESTIMATOR_COLORS[[estimator]], lwd = 1.1
      )
      points(
        positions, part$coverage, pch = 16,
        col = ESTIMATOR_COLORS[[estimator]]
      )
    }

    legend(
      "right", inset = c(-0.20, 0), xpd = NA,
      legend = c(
        expression(psi[L]), expression(psi[0]), expression(psi[1]),
        expression(psi[2]), expression(psi[3])
      ),
      col = ESTIMATOR_COLORS[ESTIMATORS], pch = 16, bty = "n", cex = 0.82
    )
    mtext(
      "Kịch bản mô phỏng 2: tái thực nghiệm giảm tải",
      side = 3, outer = TRUE, line = 0.8, cex = 1.05
    )
    mtext(
      paste0(
        config$repetitions,
        " lần lặp; đoạn thẳng đứng trong đồ thị độ bao phủ là khoảng nhị thức chính xác 95%"
      ),
      side = 1, outer = TRUE, line = 0.8, cex = 0.70, col = "#555555"
    )
  }
  save_report_plot(path, 13.2, 4.8, draw)
}

plot_real_results <- function(
  project_root, dataset, diagnostics, language = "R"
) {
  path <- file.path(
    project_root, "figures", language, paste0(dataset, "_importance.png")
  )
  learners <- c("linear", "additive", "forest")
  learners <- learners[learners %in% unique(diagnostics$learner)]
  estimator_order <- c(
    "psi_L", "psi_0", "psi_1", "psi_2", "psi_3", "cpi_residual", "pfi"
  )
  estimator_labels <- estimator_math_labels(
    estimator_order, include_method = TRUE
  )
  finite_limits <- unlist(diagnostics[c(
    "ci_low_cross_scaled_by_var_y", "ci_high_cross_scaled_by_var_y"
  )])
  finite_limits <- finite_limits[is.finite(finite_limits)]
  x_limits <- if (length(finite_limits)) {
    limits <- range(c(0, finite_limits))
    padding <- 0.06 * max(diff(limits), 0.1)
    limits + c(-padding, padding)
  } else c(-0.1, 0.1)
  draw <- function() {
    par(
      mfrow = c(1, length(learners)), mar = c(4.5, 8.2, 3.2, 1),
      oma = c(3.8, 0, 3, 0), family = "sans", bty = "l", las = 1,
      cex.axis = 0.86, cex.lab = 0.92
    )
    for (learner in learners) {
      part <- diagnostics[diagnostics$learner == learner, , drop = FALSE]
      part <- part[match(estimator_order, part$estimator), , drop = FALSE]
      part <- part[!is.na(part$estimator), , drop = FALSE]
      y <- rev(seq_len(nrow(part)))
      colors <- ESTIMATOR_COLORS[part$estimator]
      plot(
        part$estimate_scaled_by_var_y, y,
        pch = 16, col = colors, xlim = x_limits, yaxt = "n",
        xlab = "", ylab = "", main = tools::toTitleCase(learner)
      )
      axis(
        2, at = y, labels = estimator_labels[match(part$estimator, estimator_order)],
        las = 1, tick = FALSE, cex.axis = 0.76
      )
      grid(col = "#E8E8E8", lwd = 0.7)
      abline(v = 0, col = "#777777", lty = 3, lwd = 0.9)
      segments(
        part$ci_low_cross_scaled_by_var_y, y,
        part$ci_high_cross_scaled_by_var_y, y,
        col = colors, lwd = 2.1
      )
      points(part$estimate_scaled_by_var_y, y, pch = 16, col = colors)
      missing <- which(!is.finite(part$estimate_scaled_by_var_y))
      if (length(missing)) {
        text(
          x_limits[[1L]] + 0.02 * diff(x_limits), y[missing],
          labels = "không định danh", col = colors[missing],
          adj = c(0, 0.5), cex = 0.70
        )
      }
    }
    mtext(
      paste0(tools::toTitleCase(dataset), ": độ quan trọng của biến"),
      side = 3, outer = TRUE, line = 1, cex = 1.08
    )
    mtext(
      bquote("Ước lượng và khoảng t-Cross 95% không hiệu chỉnh, chuẩn hóa theo" ~ Var(Y)),
      side = 1, outer = TRUE, line = 2.0, cex = 0.88
    )
    mtext(
      "Các mô hình nuisance sử dụng cùng một thang ngang.",
      side = 1, outer = TRUE, line = 0.35, cex = 0.70, col = "#555555"
    )
  }
  save_report_plot(path, 12.0, 5.0, draw)
}

plot_example1_empirical <- function(
  project_root, aggregate, config, language = "R"
) {
  path <- file.path(
    project_root, "figures", language,
    paste0("example1_empirical_", config$name, ".png")
  )
  part <- aggregate[aggregate$example == 1L, , drop = FALSE]
  learners <- c("linear", "additive", "forest")
  learners <- learners[learners %in% unique(part$learner)]
  stable_estimators <- c("psi_L", "psi_1", "psi_2", "psi_3")

  draw <- function() {
    layout(
      matrix(
        c(
          1, 2, 3, 7,
          4, 5, 6, 8
        ),
        nrow = 2,
        byrow = TRUE
      ),
      widths = c(1, 1, 1, 0.40),
      heights = c(1, 1)
    )
    
    par(
      mar = c(4.3, 4.5, 3.0, 2.8),
      oma = c(3.0, 0, 3.0, 0),
      family = "sans",
      bty = "l",
      las = 1,
      cex.axis = 0.80,
      cex.lab = 0.86
    )

    for (learner in learners) {
      panel <- part[part$learner == learner, , drop = FALSE]
      rho_grid <- sort(unique(panel$rho))
      stable <- panel[panel$estimator %in% stable_estimators, , drop = FALSE]
      finite_values <- stable$mean_estimate[is.finite(stable$mean_estimate)]
      y_limits <- range(c(0, 4.2, finite_values), na.rm = TRUE)
      padding <- 0.06 * max(diff(y_limits), 1)
      y_limits <- y_limits + c(-padding, padding)

      plot(
        NA, xlim = range(rho_grid), ylim = y_limits,
        xlab = expression(paste("Tương quan ", rho)),
        ylab = "Ước lượng trung bình",
        main = tools::toTitleCase(learner)
      )
      grid(col = "#E8E8E8", lwd = 0.7)
      lines(
        rho_grid, 4 * (1 - rho_grid^2),
        col = "#444444", lty = 3, lwd = 1.4
      )
      abline(h = 4, col = "#444444", lty = 2, lwd = 1.4)
      for (estimator in stable_estimators) {
        series <- panel[panel$estimator == estimator, , drop = FALSE]
        series <- series[order(series$rho), , drop = FALSE]
        if (!nrow(series)) next
        lines(
          series$rho, series$mean_estimate,
          col = ESTIMATOR_COLORS[[estimator]], lwd = 1.8
        )
        points(
          series$rho, series$mean_estimate,
          col = ESTIMATOR_COLORS[[estimator]], pch = 16, cex = 0.70
        )
      }
    }

    for (learner_index in seq_along(learners)) {
      learner <- learners[[learner_index]]
      panel <- part[part$learner == learner, , drop = FALSE]
      psi0 <- panel[panel$estimator == "psi_0", , drop = FALSE]
      psi0 <- psi0[order(psi0$rho), , drop = FALSE]
      finite_estimates <- psi0$mean_estimate[is.finite(psi0$mean_estimate)]
      y_limits <- range(c(4, finite_estimates), na.rm = TRUE)
      padding <- 0.08 * max(diff(y_limits), 1)
      y_limits <- y_limits + c(-padding, padding)

      plot(
        psi0$rho, psi0$mean_estimate,
        type = "b", pch = 16, lwd = 1.8,
        col = ESTIMATOR_COLORS[["psi_0"]],
        xlab = expression(paste("Tương quan ", rho)),
        ylab = bquote("Ước lượng trung bình của" ~ psi[0]),
        main = bquote(psi[0] ~ "và ESS của tỷ số mật độ"),
        ylim = y_limits
      )
      grid(col = "#E8E8E8", lwd = 0.7)
      abline(h = 4, col = "#444444", lty = 2, lwd = 1.4)

      par(new = TRUE)
      plot(
        psi0$rho, psi0$density_ratio_ess_fraction_mean,
        type = "b", pch = 17, lty = 3, lwd = 1.4,
        col = "#6A3D9A", axes = FALSE, xlab = "", ylab = "",
        ylim = c(0, 1)
      )
      if (learner_index == length(learners)) {
        axis(
          4,
          at = seq(0, 1, by = 0.2),
          las = 1,
          cex.axis = 0.74
        )
    
        mtext(
          "Tỷ lệ ESS",
          side = 4,
          line = 2.0,
          cex = 0.78
        )
      }
    }

    # Panel chú giải cho hàng trên
    plot.new()
    
    legend(
      "center",
      legend = as.expression(list(
        bquote(psi[L]),
        bquote(psi[1]),
        bquote(psi[2]),
        bquote(psi[3]),
        bquote("Giá trị giải tích của" ~ psi[L]),
        bquote("Giá trị giải tích của" ~ psi[0])
      )),
      col = c(
        ESTIMATOR_COLORS[stable_estimators],
        "#444444",
        "#444444"
      ),
      lty = c(
        rep(1, length(stable_estimators)),
        3,
        2
      ),
      lwd = c(
        rep(1.8, length(stable_estimators)),
        1.4,
        1.4
      ),
      pch = c(
        rep(16, length(stable_estimators)),
        NA,
        NA
      ),
      bty = "n",
      cex = 0.82,
      y.intersp = 1.15,
      x.intersp = 0.8
    )
    
    # Panel chú giải cho hàng dưới
    plot.new()
    
    legend(
      "center",
      legend = as.expression(list(
        bquote("Ước lượng" ~ psi[0]),
        "Tỷ lệ ESS của tỷ số mật độ",
        bquote("Giá trị chuẩn" ~ psi[0] == 4)
      )),
      col = c(
        ESTIMATOR_COLORS[["psi_0"]],
        "#6A3D9A",
        "#444444"
      ),
      lty = c(1, 3, 2),
      pch = c(16, 17, NA),
      lwd = c(1.8, 1.4, 1.4),
      bty = "n",
      cex = 0.82,
      y.intersp = 1.25,
      x.intersp = 0.8
    )
    mtext(
      bquote("Kịch bản mô phỏng 1: sai lệch do tương quan và độ ổn định của" ~ psi[0]),
      side = 3, outer = TRUE, line = 1, cex = 1.08
    )
    mtext(
      paste0(
        "Cấu hình ", config$name, "; ", config$repetitions,
        " lần lặp cho mỗi mức tương quan"
      ),
      side = 1, outer = TRUE, line = 0.9, cex = 0.70, col = "#555555"
    )
  }
  save_report_plot(path, 16.0, 8.5, draw)
}

plot_psi0_overlap_diagnostic <- function(
  project_root, aggregate, config, language = "R"
) {
  path <- file.path(
    project_root,
    "figures",
    language,
    paste0("psi0_overlap_diagnostic_", config$name, ".png")
  )

  part <- aggregate[
    aggregate$estimator == "psi_0" &
      is.finite(aggregate$relative_absolute_bias) &
      is.finite(aggregate$density_ratio_ess_fraction_mean),
    ,
    drop = FALSE
  ]

  if (!nrow(part)) {
    stop("Không có quan sát hợp lệ để vẽ chẩn đoán psi_0 và ESS.")
  }

  learner_colors <- c(
    linear = "#1B9E77",
    additive = "#D95F02",
    forest = "#7570B3"
  )

  learner_labels <- c(
    linear = "Tuyến tính",
    additive = "Cộng tính",
    forest = "Rừng"
  )

  scenario_symbols <- c(
    `1` = 16,
    `2` = 17,
    `3` = 15,
    `4` = 18,
    `5` = 8
  )

  x <- part$density_ratio_ess_fraction_mean
  y <- pmax(part$relative_absolute_bias, 1e-4)

  # Giới hạn trục x bám theo dữ liệu nhưng vẫn chứa ngưỡng ESS = 0,10.
  x_upper <- min(
    1,
    max(
      0.30,
      max(x, na.rm = TRUE) * 1.12
    )
  )

  # Giới hạn trục y trên thang logarit.
  y_limits <- c(
    max(1e-4, min(y, na.rm = TRUE) / 1.35),
    max(y, na.rm = TRUE) * 1.35
  )

  # Các mốc trục y dễ đọc trên thang logarit.
  possible_y_ticks <- c(
    0.01, 0.02, 0.05,
    0.10, 0.20, 0.50,
    1, 2, 5, 10, 20, 50
  )

  y_ticks <- possible_y_ticks[
    possible_y_ticks >= y_limits[[1L]] &
      possible_y_ticks <= y_limits[[2L]]
  ]

  draw <- function() {
    # Cột bên phải được dành riêng cho chú giải.
    layout(
      matrix(c(1, 2), nrow = 1L),
      widths = c(4.8, 1.55)
    )

    par(
      oma = c(2.2, 0, 0.4, 0),
      family = "sans",
      bty = "l",
      las = 1
    )

    # ------------------------------------------------------------------
    # Panel chính
    # ------------------------------------------------------------------
    par(
      mar = c(5.0, 5.5, 3.8, 1.0),
      cex.axis = 0.84,
      cex.lab = 0.90
    )

    plot(
      x,
      y,
      type = "n",
      log = "y",
      xlim = c(0, x_upper),
      ylim = y_limits,
      yaxt = "n",
      xlab = "Tỷ lệ ESS trung bình của tỷ số mật độ",
      ylab = expression(
        paste("Độ lệch tương đối tuyệt đối của ", psi[0])
      ),
      main = expression(
        paste("Chẩn đoán độ ổn định của ", psi[0])
      ),
      font.main = 1
    )

    axis(
      2,
      at = y_ticks,
      labels = format(
        y_ticks,
        trim = TRUE,
        scientific = FALSE
      ),
      las = 1
    )

    # Vẽ lưới trước điểm dữ liệu để không che các điểm.
    grid(
      col = "#E8E8E8",
      lwd = 0.7
    )

    # Các ngưỡng chẩn đoán tham chiếu.
    abline(
      v = 0.10,
      lty = 2,
      lwd = 1.0,
      col = "#666666"
    )

    abline(
      h = 0.20,
      lty = 3,
      lwd = 1.0,
      col = "#666666"
    )

    points(
      x,
      y,
      pch = scenario_symbols[as.character(part$example)],
      col = learner_colors[part$learner],
      cex = 1.15,
      lwd = 1.2
    )

    # Gắn nhãn trực tiếp cho hai đường tham chiếu.
    text(
      x = 0.10,
      y = y_limits[[2L]] / 1.20,
      labels = "ESS = 0,10",
      srt = 90,
      pos = 2,
      cex = 0.68,
      col = "#555555"
    )

    text(
      x = x_upper * 0.97,
      y = 0.20,
      labels = "Sai lệch = 0,20",
      pos = 3,
      cex = 0.68,
      col = "#555555"
    )

    mtext(
      expression(
        paste(
          "Vùng thuận lợi hơn nằm về phía ESS lớn và sai lệch của ",
          psi[0],
          " nhỏ."
        )
      ),
      side = 3,
      line = 0.35,
      cex = 0.70,
      col = "#555555"
    )

    # ------------------------------------------------------------------
    # Panel chú giải
    # ------------------------------------------------------------------
    par(
      mar = c(0.5, 0.5, 0.5, 0.5)
    )

    plot.new()

    legend(
      "topleft",
      inset = c(0.02, 0.08),
      title = "Mô hình nuisance",
      legend = learner_labels[names(learner_colors)],
      col = learner_colors,
      pch = 16,
      pt.cex = 1.0,
      bty = "n",
      cex = 0.82,
      y.intersp = 1.20,
      x.intersp = 0.75
    )

    legend(
      "bottomleft",
      inset = c(0.02, 0.08),
      title = "Kịch bản mô phỏng",
      legend = scenario_axis_labels(
        names(scenario_symbols),
        compact = TRUE
      ),
      pch = scenario_symbols,
      col = "#333333",
      pt.cex = 1.0,
      bty = "n",
      cex = 0.82,
      y.intersp = 1.20,
      x.intersp = 0.75
    )

    mtext(
      "Màu biểu diễn mô hình nuisance; hình dạng điểm biểu diễn kịch bản mô phỏng.",
      side = 1,
      outer = TRUE,
      line = 0.65,
      cex = 0.68,
      col = "#555555"
    )
  }

  save_report_plot(
    path,
    width = 10.8,
    height = 5.8,
    draw = draw
  )
}

plot_coverage_heatmap <- function(
  project_root, aggregate, config, language = "R"
) {
  path <- file.path(
    project_root, "figures", language,
    paste0("coverage_heatmap_", config$name, ".png")
  )
  part <- aggregate[
    aggregate$example %in% 2:5 & aggregate$estimator %in% ESTIMATORS,
    , drop = FALSE
  ]
  row_keys <- c(
    "KB 2 - Linear", "KB 2 - Additive", "KB 2 - Forest",
    "KB 3 - Linear", "KB 3 - Additive", "KB 3 - Forest",
    "KB 4 - Linear", "KB 4 - Additive", "KB 4 - Forest",
    "KB 5 - Linear", "KB 5 - Additive", "KB 5 - Forest"
  )
  learner_labels <- c(linear = "Linear", additive = "Additive", forest = "Forest")
  part$row_key <- paste0(
    "KB ", part$example, " - ", learner_labels[part$learner]
  )

  coverage_matrix <- matrix(
    NA_real_, nrow = length(row_keys), ncol = length(ESTIMATORS),
    dimnames = list(row_keys, ESTIMATORS)
  )
  difference_matrix <- coverage_matrix
  for (index in seq_len(nrow(part))) {
    key <- part$row_key[[index]]
    estimator <- part$estimator[[index]]
    coverage_matrix[key, estimator] <- part$coverage[[index]]
    difference_matrix[key, estimator] <- part$coverage_minus_paper[[index]]
  }

  draw_matrix <- function(values, palette, z_limits, title, digits = 2) {
    image(
      seq_len(ncol(values)), seq_len(nrow(values)),
      t(values[nrow(values):1, , drop = FALSE]),
      col = palette, zlim = z_limits, axes = FALSE,
      xlab = "Bộ ước lượng", ylab = "", main = title
    )
    axis(
      1, at = seq_len(ncol(values)),
      labels = estimator_math_labels(colnames(values))
    )
    axis(
      2, at = seq_len(nrow(values)),
      labels = rev(rownames(values)), las = 1, tick = FALSE
    )
    for (row in seq_len(nrow(values))) {
      for (column in seq_len(ncol(values))) {
        value <- values[row, column]
        if (!is.finite(value)) next
        y_position <- nrow(values) - row + 1L
        text(
          column, y_position,
          labels = formatC(value, format = "f", digits = digits),
          col = if (
            abs(value) > 0.55 && identical(z_limits, c(-1, 1))
          ) "white" else if (
            value >= 0.65 && identical(z_limits, c(0, 1))
          ) "white" else "#222222",
          cex = 0.70
        )
      }
    }
  }

  draw <- function() {
    par(
      mfrow = c(1, 2), mar = c(5.4, 8.6, 3.5, 1.2),
      oma = c(3.3, 0, 3.0, 0), family = "sans", bty = "n",
      las = 1, cex.axis = 0.78
    )
    coverage_palette <- grDevices::colorRampPalette(
      c("#F7FBFF", "#6BAED6", "#08306B")
    )(101)
    difference_palette <- grDevices::colorRampPalette(
      c("#B2182B", "#F7F7F7", "#2166AC")
    )(101)

    draw_matrix(
      coverage_matrix, coverage_palette, c(0, 1),
      expression(paste("Độ bao phủ t-Cross giảm tải đối với ", psi[0]))
    )
    draw_matrix(
      difference_matrix, difference_palette, c(-1, 1),
      "Chênh lệch so với độ bao phủ trong Bảng 2"
    )
    mtext(
      paste0(
        "Mỗi ô sử dụng ", config$repetitions,
        " lần lặp; khoảng nhị thức chính xác được lưu trong tệp CSV."
      ),
      side = 1, outer = TRUE, line = 1.0, cex = 0.70, col = "#555555"
    )
    mtext(
      "So sánh độ bao phủ cho các Kịch bản mô phỏng 2–5",
      side = 3, outer = TRUE, line = 1, cex = 1.08
    )
  }
  save_report_plot(path, 13.0, 6.8, draw)
}

plot_baseline_comparison <- function(
  project_root, comparison, config, language = "R"
) {
  path <- file.path(
    project_root, "figures", language,
    paste0("baseline_comparison_", config$name, ".png")
  )
  methods <- c("DVI psi_0", "LOCO", "Residual CPI", "Marginal PFI")
  learners <- c("linear", "additive", "forest")
  learners <- learners[learners %in% unique(comparison$learner)]
  examples <- sort(unique(comparison$example))
  offsets <- seq(-0.24, 0.24, length.out = length(methods))
  method_colors <- c(
    "DVI psi_0" = ESTIMATOR_COLORS[["psi_0"]],
    "LOCO" = ESTIMATOR_COLORS[["psi_L"]],
    "Residual CPI" = ESTIMATOR_COLORS[["cpi_residual"]],
    "Marginal PFI" = ESTIMATOR_COLORS[["pfi"]]
  )
  finite_error <- comparison$relative_rmse[
    is.finite(comparison$relative_rmse) &
      comparison$relative_rmse > 0
  ]
  y_limits <- if (length(finite_error)) {
    range(c(1e-3, finite_error), na.rm = TRUE)
  } else c(1e-3, 1)
  y_limits[[2L]] <- max(y_limits[[2L]] * 1.35, 0.3)

  draw <- function() {
    par(
      mfrow = c(1, length(learners)), mar = c(4.8, 4.8, 3.2, 1),
      oma = c(3.2, 0, 3, 7.4), family = "sans", bty = "l", las = 1,
      cex.axis = 0.82, cex.lab = 0.88
    )
    for (learner in learners) {
      panel <- comparison[comparison$learner == learner, , drop = FALSE]
      plot(
        NA, xlim = c(0.6, length(examples) + 0.4), ylim = y_limits,
        log = "y", xaxt = "n", xlab = "Kịch bản mô phỏng",
        ylab = bquote(RMSE / abs(psi[0])),
        main = tools::toTitleCase(learner)
      )
      axis(1, at = seq_along(examples), labels = scenario_axis_labels(examples, compact = TRUE))
      grid(col = "#E8E8E8", lwd = 0.7)
      abline(h = c(0.1, 0.2), lty = c(2, 3), col = "#777777")
      for (method_index in seq_along(methods)) {
        method <- methods[[method_index]]
        method_part <- panel[
          panel$method_label == method,
          ,
          drop = FALSE
        ]
        method_part <- method_part[
          match(examples, method_part$example),
          ,
          drop = FALSE
        ]
        values <- pmax(method_part$relative_rmse, 1e-3)
        points(
          seq_along(examples) + offsets[[method_index]], values,
          pch = 16, col = method_colors[[method]], cex = 0.9
        )
      }
    }
    legend(
      "right", inset = c(-0.21, 0), xpd = NA, bty = "n",
      legend = as.expression(list(
        bquote("DVI" ~ psi[0]), "LOCO", "Residual CPI", "Marginal PFI"
      )),
      col = method_colors[methods], pch = 16, cex = 0.78
    )
    mtext(
      "So sánh với oracle trên dữ liệu mô phỏng",
      side = 3, outer = TRUE, line = 1, cex = 1.08
    )
    mtext(
      "RMSE phản ánh đồng thời độ lệch hệ thống và mức biến thiên. Marginal PFI là chẩn đoán có estimand khác.",
      side = 1, outer = TRUE, line = 1, cex = 0.70, col = "#555555"
    )
  }
  save_report_plot(path, 13.6, 4.9, draw)
}

plot_ablation_summary <- function(
  project_root, summary, config, language = "R"
) {
  path <- file.path(
    project_root, "figures", language,
    paste0("correction_ablation_", config$name, ".png")
  )
  estimators <- c("psi_L", "psi_1", "psi_2", "psi_3")
  modes <- c("none", "variance_scale", "paper_literal")
  mode_labels <- list(
    none = "Không hiệu chỉnh",
    variance_scale = bquote(c == Var(Y)),
    paper_literal = bquote(c == Var(Y)^2)
  )
  offsets <- c(-0.22, 0, 0.22)
  draw <- function() {
    par(
      mfrow = c(1, 2), mar = c(5.0, 4.8, 3.2, 1),
      oma = c(2.9, 0, 3, 7.4), family = "sans", bty = "l", las = 1,
      cex.axis = 0.82, cex.lab = 0.90
    )
    plot(
      NA, xlim = c(0.6, length(estimators) + 0.4), ylim = c(0, 1),
      xaxt = "n", xlab = "Bộ ước lượng", ylab = "Độ bao phủ",
      main = "Độ bao phủ và khoảng nhị thức chính xác 95%"
    )
    axis(1, at = seq_along(estimators), labels = estimator_math_labels(estimators))
    grid(col = "#E8E8E8", lwd = 0.7)
    abline(h = 0.95, lty = 2, col = "#555555")
    for (mode_index in seq_along(modes)) {
      mode <- modes[[mode_index]]
      part <- summary[summary$correction_mode == mode, , drop = FALSE]
      part <- part[match(estimators, part$estimator), , drop = FALSE]
      x <- seq_along(estimators) + offsets[[mode_index]]
      arrows(
        x, part$coverage_ci_low, x, part$coverage_ci_high,
        angle = 90, code = 3, length = 0.035,
        col = CORRECTION_COLORS[[mode]], lwd = 1.4
      )
      points(x, part$coverage, pch = 16, col = CORRECTION_COLORS[[mode]])
    }

    finite_inflation <- summary$width_inflation_vs_none[
      is.finite(summary$width_inflation_vs_none) &
        summary$width_inflation_vs_none > 0
    ]
    inflation_limits <- if (length(finite_inflation)) {
      c(
        max(min(finite_inflation) * 0.8, 0.5),
        max(finite_inflation) * 1.25
      )
    } else c(0.5, 2)
    plot(
      NA, xlim = c(0.6, length(estimators) + 0.4),
      ylim = inflation_limits, log = "y", xaxt = "n",
      xlab = "Bộ ước lượng", ylab = "Độ rộng / độ rộng không hiệu chỉnh",
      main = "Hệ số phóng đại độ rộng khoảng"
    )
    axis(1, at = seq_along(estimators), labels = estimator_math_labels(estimators))
    grid(col = "#E8E8E8", lwd = 0.7)
    abline(h = 1, lty = 2, col = "#555555")
    for (mode_index in seq_along(modes)) {
      mode <- modes[[mode_index]]
      part <- summary[summary$correction_mode == mode, , drop = FALSE]
      part <- part[match(estimators, part$estimator), , drop = FALSE]
      points(
        seq_along(estimators) + offsets[[mode_index]],
        part$width_inflation_vs_none,
        pch = 16, col = CORRECTION_COLORS[[mode]], cex = 0.9
      )
    }
    legend(
      "right", inset = c(-0.21, 0), xpd = NA, bty = "n",
      legend = as.expression(unname(mode_labels[modes])),
      col = CORRECTION_COLORS[modes],
      pch = 16, cex = 0.78
    )
    mtext(
      "Ablation số hạng hiệu chỉnh bảo thủ của t-Cross",
      side = 3, outer = TRUE, line = 1, cex = 1.08
    )
    mtext(
      paste0(
        config$ablation_repetitions,
        " lần lặp; dữ liệu và fold giống nhau giữa các chế độ hiệu chỉnh"
      ),
      side = 1, outer = TRUE, line = 0.8, cex = 0.70, col = "#555555"
    )
  }
  save_report_plot(path, 13.2, 4.9, draw)
}

plot_real_predictability <- function(
  project_root,
  dataset,
  predictability,
  diagnostics = NULL,
  language = "R"
) {
  path <- file.path(
    project_root, "figures", language, paste0(dataset, "_predictability.png")
  )
  learners <- predictability$learner
  ess <- rep(NA_real_, length(learners))
  if (!is.null(diagnostics) && nrow(diagnostics)) {
    psi0 <- diagnostics[
      diagnostics$estimator == "psi_0",
      c("learner", "density_ratio_ess_fraction_mean"),
      drop = FALSE
    ]
    ess <- psi0$density_ratio_ess_fraction_mean[
      match(learners, psi0$learner)
    ]
  }

  draw <- function() {
    par(
      mfrow = c(1, 3), mar = c(5.0, 4.7, 3.2, 1),
      oma = c(2.6, 0, 3, 0), family = "sans", bty = "l", las = 1,
      cex.axis = 0.80, cex.lab = 0.86
    )
    colors <- ESTIMATOR_COLORS[c("psi_L", "psi_1", "psi_0")][
      seq_along(learners)
    ]
    barplot(
      predictability$predictability_r2,
      names.arg = tools::toTitleCase(learners), ylim = c(0, 1),
      ylab = bquote(R^2 ~ "cross-fitted của" ~ X ~ "|" ~ Z),
      main = bquote("Khả năng dự đoán" ~ X), col = colors, border = NA
    )
    grid(nx = NA, ny = NULL, col = "#E8E8E8", lwd = 0.7)
    abline(h = 0.90, lty = 2, col = "#555555")

    ratio <- pmax(predictability$residual_variance_ratio, 1e-16)
    barplot(
      ratio,
      names.arg = tools::toTitleCase(learners), log = "y",
      ylab = "Tỷ lệ phương sai phần dư",
      main = bquote("Biến thiên phần dư của" ~ X), col = colors, border = NA
    )
    grid(nx = NA, ny = NULL, col = "#E8E8E8", lwd = 0.7)
    abline(h = c(0.10, 1e-8), lty = c(2, 3), col = "#555555")

    finite_ess <- ess[is.finite(ess)]
    ess_limits <- if (length(finite_ess)) c(0, max(1, finite_ess)) else c(0, 1)
    barplot(
      ess,
      names.arg = tools::toTitleCase(learners), ylim = ess_limits,
      ylab = "Tỷ lệ ESS",
      main = "Độ ổn định của tỷ số mật độ", col = colors, border = NA
    )
    grid(nx = NA, ny = NULL, col = "#E8E8E8", lwd = 0.7)
    abline(h = 0.10, lty = 2, col = "#555555")

    mtext(
      paste(tools::toTitleCase(dataset), "- chẩn đoán điều kiện định danh"),
      side = 3, outer = TRUE, line = 1, cex = 1.08
    )
    mtext(
      "Khả năng dự đoán cao, tỷ lệ phần dư nhỏ hoặc ESS thấp cho thấy overlap hiệu dụng yếu.",
      side = 1, outer = TRUE, line = 0.8, cex = 0.70, col = "#555555"
    )
  }
  save_report_plot(path, 12.2, 4.7, draw)
}
