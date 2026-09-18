col_sr   <- "green4"        # the proposed, specification-robust interval
col_pt   <- "dodgerblue2"   # point estimates
col_ax   <- "grey35"        # axes and frames
col_kl   <- "saddlebrown"   # the KL panel, kept off the estimate palette
col_href <- grDevices::adjustcolor("red4", alpha.f = 0.55)   # hull rules
col_fill <- grDevices::adjustcolor("red4", alpha.f = 0.08)   # hull band
x_label     <- "moment order k"

display <- c(set1 = "Collection 1", set2 = "Collection 2")
disp <- function(lab) if (lab %in% names(display)) display[[lab]] else lab

kl_column <- function(L, type)
  if (type == "in-sample") L$kl_in_sample else L$kl

kl_axis <- function(L, type) {
  r   <- range(kl_column(L, type))
  pad <- if (diff(r) > 0) 0.08*diff(r) else 0.05*max(abs(r[1]), 1)
  yl  <- c(r[1] - pad, r[2] + pad)
  yt  <- pretty(r, 3)
  list(ylim = yl, ticks = yt[yt >= yl[1] & yt <= yl[2]])
}

estimate_panel <- function(L) {
  k  <- L$k
  hl <- L$hull_lo[1]; hh <- L$hull_hi[1]
  yr <- range(c(L$ci_lo, L$ci_hi, hl, hh))
  yr <- yr + c(-0.05, 0.035) * diff(yr)
  plot(NA, xlim = range(k) + c(-0.18, 0.75), ylim = yr, axes = FALSE,
       xlab = "", ylab = "ATE (in 1000 USD)")

  rect(par("usr")[1], hl, par("usr")[2], hh, col = col_fill, border = NA)
  abline(h = c(hl, hh), col = col_href, lty = 2, lwd = 2)
  text(par("usr")[1], hh, "convex hull", col = col_href,
       adj = c(-0.08, 2.0), cex = 1.0)

  cw <- 0.12
  segments(k, L$ci_lo, k, L$ci_hi, col = col_sr, lwd = 2.4, lend = 1)
  segments(k - cw, L$ci_lo, k + cw, L$ci_lo, col = col_sr, lwd = 2.4)
  segments(k - cw, L$ci_hi, k + cw, L$ci_hi, col = col_sr, lwd = 2.4)
  points(k, L$estimate, pch = 19, col = col_pt, cex = 1.2)

  text(k, L$ci_hi, sprintf("%.2f", L$width), col = col_sr,
       adj = c(-0.12, 2.1), cex = 0.88)

  box(col = col_ax)
  axis(1, at = k, col = col_ax, col.axis = col_ax, labels = FALSE)
  axis(2, col = col_ax, col.axis = col_ax, las = 1)
}

kl_panel <- function(L, type) {
  k <- L$k; ka <- kl_axis(L, type); y <- kl_column(L, type)
  plot(NA, xlim = range(k) + c(-0.18, 0.75), ylim = ka$ylim, axes = FALSE,
       xlab = x_label, ylab = paste0("KL div (", type, ")"))
  lines(k, y, col = col_kl, lwd = 2)
  points(k, y, pch = 19, col = col_kl, cex = 1.2)
  box(col = col_ax)
  axis(1, at = k, col = col_ax, col.axis = col_ax)
  axis(2, at = ka$ticks, col = col_ax, col.axis = col_ax, las = 1)
}

write_final_figures <- function(LAD, OUTDIR,
                                type = c("out-of-fold", "in-sample")) {
  type <- match.arg(type)
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
  for (lab in names(LAD)) {
    L <- LAD[[lab]]; if (is.null(L)) next
    fp <- file.path(OUTDIR, sprintf("ladder_%s%s.pdf", lab,
                                    if (type == "in-sample") "_in_sample_KL" else ""))
    pdf(fp, width = 10, height = 7.5, pointsize = 14)
    op <- par(mgp = c(2.6, 0.6, 0), oma = c(0, 0, 0.6, 0.6),
              cex.axis = 1.05, cex.lab = 1.12, col.lab = col_ax)
    layout(matrix(1:2, ncol = 1), heights = c(2, 1.35))
    par(mar = c(0, 4.6, 0.4, 0.6));   estimate_panel(L)
    par(mar = c(4.2, 4.6, 0, 0.6));   kl_panel(L, type)
    layout(1); par(op); invisible(dev.off())
    cat("Wrote ", fp, "\n", sep = "")
  }
}

if (sys.nframe() == 0) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1) {
    rds <- a[1]
    out <- if (length(a) >= 2) a[2] else dirname(rds)
    type <- if (length(a) >= 3) a[3] else "out-of-fold"
    write_final_figures(readRDS(rds), out, type)
  } else {
    cat("usage: Rscript ladder_final_plot.R <ladder .rds> [outdir] [out-of-fold|in-sample]\n")
  }
}
