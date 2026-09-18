COL_SR   <- "green4"        # the proposed, specification-robust interval
COL_PT   <- "dodgerblue2"   # point estimates
COL_AX   <- "grey35"        # axes and frames
COL_KL   <- "saddlebrown"   # the KL panel, kept off the estimate palette
COL_HREF <- grDevices::adjustcolor("red4", alpha.f = 0.55)   # hull rules
COL_FILL <- grDevices::adjustcolor("red4", alpha.f = 0.08)   # hull band
XLAB     <- "moment order k"

DISPLAY <- c(set1 = "Collection 1", set2 = "Collection 2")
disp <- function(lab) if (lab %in% names(DISPLAY)) DISPLAY[[lab]] else lab

kl_axis <- function(L) {
  r   <- range(L$kl)
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
  plot(NA, xlim = range(k) + c(-0.18, 0.18), ylim = yr, axes = FALSE,
       xlab = "", ylab = "ATE (in 1000 USD)")

  rect(par("usr")[1], hl, par("usr")[2], hh, col = COL_FILL, border = NA)
  abline(h = c(hl, hh), col = COL_HREF, lty = 2, lwd = 2)
  text(par("usr")[1], hh, "convex hull", col = COL_HREF,
       adj = c(-0.08, 2.0), cex = 1.0)

  cw <- 0.12
  segments(k, L$ci_lo, k, L$ci_hi, col = COL_SR, lwd = 2.4, lend = 1)
  segments(k - cw, L$ci_lo, k + cw, L$ci_lo, col = COL_SR, lwd = 2.4)
  segments(k - cw, L$ci_hi, k + cw, L$ci_hi, col = COL_SR, lwd = 2.4)
  points(k, L$estimate, pch = 19, col = COL_PT, cex = 1.2)

  text(k, L$ci_hi, sprintf("%.2f", L$width), col = COL_SR,
       adj = c(0.5, -0.9), cex = 0.88)

  box(col = COL_AX)
  axis(1, at = k, col = COL_AX, col.axis = COL_AX, labels = FALSE)
  axis(2, col = COL_AX, col.axis = COL_AX, las = 1)
}

kl_panel <- function(L) {
  k <- L$k; ka <- kl_axis(L)
  plot(NA, xlim = range(k) + c(-0.18, 0.18), ylim = ka$ylim, axes = FALSE,
       xlab = XLAB, ylab = "KL divergence")
  lines(k, L$kl, col = COL_KL, lwd = 2)
  points(k, L$kl, pch = 19, col = COL_KL, cex = 1.2)
  box(col = COL_AX)
  axis(1, at = k, col = COL_AX, col.axis = COL_AX)
  axis(2, at = ka$ticks, col = COL_AX, col.axis = COL_AX, las = 1)
}

write_final_figures <- function(LAD, OUTDIR) {
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
  for (lab in names(LAD)) {
    L <- LAD[[lab]]; if (is.null(L)) next
    fp <- file.path(OUTDIR, sprintf("ladder_%s.pdf", lab))
    pdf(fp, width = 10, height = 7.5, pointsize = 14)
    op <- par(mgp = c(2.6, 0.6, 0), oma = c(0, 0, 0.6, 0.6),
              cex.axis = 1.05, cex.lab = 1.12, col.lab = COL_AX)
    layout(matrix(1:2, ncol = 1), heights = c(2, 1))
    par(mar = c(0, 4.6, 0.4, 0.6));   estimate_panel(L)
    par(mar = c(4.2, 4.6, 0, 0.6));   kl_panel(L)
    layout(1); par(op); invisible(dev.off())
    cat("Wrote ", fp, "\n", sep = "")
  }
}

if (!exists(".LADDER_FINAL_SOURCED")) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1) {
    rds <- a[1]
    out <- if (length(a) >= 2) a[2] else dirname(rds)
    write_final_figures(readRDS(rds), out)
  } else {
    cat("usage: Rscript ladder_final_plot.R <ladder .rds> [outdir]\n")
  }
}
