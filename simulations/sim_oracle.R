# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# Reviewer material for Section 5.1: the oracle reweighted estimand, the
# covariate shift for both examples in the style of the paper, and the
# moments of the original and reweighted populations.  Run from inside the
# simulations folder:
#   Rscript sim_oracle.R
# tau_R goes to tauR.txt, which sim_eval.R reads; the supporting numbers to
# one_run.txt, and the figures alongside.
# ==============================================================================

rm(list=ls())
library(specrobust)

results_dir <- "."

n_large <- 5e6
n_small <- 1000
alpha = 0.05
adj_sets <- list(c("X1"), c("X1","X2"))

# the paper's overlay: 200 bins, both populations, and a smoothed outline
overlay <- function(data, weights, breaks = 200, hcol = "orange",
                    hcol2 = "skyblue2", xlim = NULL, ylim = NULL, xlab = NULL,
                    legend.text = c("original population",
                                    "new target population"),
                    legend.pos = "topleft") {
  h1 <- hist(data, breaks = breaks, plot = FALSE)
  bw <- diff(h1$breaks)
  d0 <- h1$counts / sum(h1$counts) / bw
  wc <- sapply(seq_along(h1$counts), function(i)
    sum(weights[data >= h1$breaks[i] & data < h1$breaks[i + 1]]))
  d1 <- wc / sum(weights) / bw
  if (is.null(xlim)) xlim <- range(data)
  if (is.null(ylim)) ylim <- c(0, 1.05 * max(c(d0, d1)))
  plot(h1$mids, d0, type = "n", xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = "Density", main = NULL)
  legend(legend.pos, legend = legend.text,
         fill = c(adjustcolor(hcol, alpha.f = 0.6),
                  adjustcolor(hcol2, alpha.f = 0.6)),
         col = c(adjustcolor(hcol, alpha.f = 0.6),
                 adjustcolor(hcol2, alpha.f = 0.6)), bty = "n", cex = 1)
  rect(h1$breaks[-length(h1$breaks)], 0, h1$breaks[-1], d0,
       col = adjustcolor(hcol, alpha.f = 0.6),
       border = adjustcolor(hcol, alpha.f = 0.2))
  rect(h1$breaks[-length(h1$breaks)], 0, h1$breaks[-1], d1,
       col = adjustcolor(hcol2, alpha.f = 0.6),
       border = adjustcolor(hcol2, alpha.f = 0.2))
  k <- 5
  lines(h1$mids, stats::filter(d0, rep(1/k, k), sides = 2) + 5e-4,
        col = hcol, lwd = 3)
  lines(h1$mids, stats::filter(d1, rep(1/k, k), sides = 2) + 5e-4,
        col = hcol2, lwd = 3)
  box()
}

draw_eg1 <- function(n) {
  set.seed(123)
  X1 = rnorm(n); X2 = rnorm(n)
  A = as.numeric(runif(n) <= 1/(exp(5*X1+5*X2)+1))
  Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
  list(Y = Y, A = A, X = data.frame(X1, X2), tau = 1)
}
draw_eg2 <- function(n) {
  set.seed(123)
  U1 = rnorm(n); U2 = rnorm(n)
  X1 = rnorm(n)
  A = as.numeric(U1 + X1 > 0)
  X2 = U1 + U2
  Y = ifelse(A==1, 1 + 0.5*X1 - U2 + rnorm(n), 5*U2 + rnorm(n))
  list(Y = Y, A = A, X = data.frame(X1, X2), tau = 1)
}

wmean <- function(x, w) sum(w*x)/sum(w)
wvar  <- function(x, w) sum(w*(x - wmean(x, w))^2)/sum(w)

out_lines <- character(0)
say <- function(...) out_lines <<- c(out_lines, sprintf(...))
tauR <- c()

for (eg in c("Example 1", "Example 2")) {
  gen <- if (eg == "Example 1") draw_eg1 else draw_eg2
  lim <- if (eg == "Example 1")
    list(X1 = list(xlim = c(-4.5, 4.8), ylim = c(0, 0.45), pos = "topleft"),
         X2 = list(xlim = c(-5, 5),     ylim = c(0, 0.45), pos = "topleft"))
  else
    list(X1 = list(xlim = c(-4.8, 4.8), ylim = c(0, 0.45), pos = "topright"),
         X2 = list(xlim = c(-7, 7),     ylim = c(0, 0.32), pos = "topleft"))

  say("================ %s ================", eg)

  d <- gen(n_small)
  s <- specrobust(d$Y, d$A, d$X, adj_sets, reg_mode = "lm", alpha = alpha,
                  n_boot = 1000, seed = 123)
  say("n = %d, one replication (seed 123)", n_small)
  say("  true ATE  tau_bar               : %.4f", d$tau)
  for (k in 1:2)
    say("  adj set %d  estimate of tau_bar   : %8.4f   CI [%8.4f, %8.4f]",
        k, s$candidate_estimates[k], s$candidate_ci[k,1], s$candidate_ci[k,2])
  say("  convex hull of the above        :            CI [%8.4f, %8.4f]",
      s$hull_ci[1], s$hull_ci[2])
  say("  proposed estimand  tau_R        : %8.4f   CI [%8.4f, %8.4f]",
      s$estimate, s$ci[1], s$ci[2])
  say("  lambda*                         : %8.4f", s$lambda)

  d <- gen(n_large)
  o <- specrobust(d$Y, d$A, d$X, adj_sets, reg_mode = "lm", alpha = alpha,
                  n_boot = 0, seed = 123)
  w <- o$weights; X1 <- d$X$X1; X2 <- d$X$X2
  say("n = %g, oracle (seed 123)", n_large)
  say("  oracle tau_R                    : %8.4f", o$estimate)
  tauR[eg] <- o$estimate
  say("  lambda*                         : %8.4f", o$lambda)
  say("  reweighted estimates per set    : %s",
      paste(sprintf("%.4f", o$reweighted_estimates), collapse = ", "))
  say("  moments of the shared covariate X1 and of X2:")
  say("    X1  original   mean %8.4f  var %8.4f", mean(X1), var(X1))
  say("    X1  reweighted mean %8.4f  var %8.4f", wmean(X1, w), wvar(X1, w))
  say("    X2  original   mean %8.4f  var %8.4f", mean(X2), var(X2))
  say("    X2  reweighted mean %8.4f  var %8.4f", wmean(X2, w), wvar(X2, w))
  say("  implied mean shift of X1        : %8.4f", wmean(X1, w) - mean(X1))
  say("")

  tag <- if (eg == "Example 1") "eg1" else "eg2"
  for (v in c("X1", "X2")) {
    f <- file.path(results_dir, sprintf("oracle_shift_%s_%s.pdf", tag, v))
    pdf(f, width = 5.5, height = 4.2)
    par(mar = c(4.0, 4.2, 1.0, 1.0), mgp = c(2.3, 0.8, 0))
    overlay(d$X[[v]], w, breaks = 200, xlab = v,
            xlim = lim[[v]]$xlim, ylim = lim[[v]]$ylim,
            legend.pos = lim[[v]]$pos)
    invisible(dev.off())
    cat("Wrote ", f, "\n", sep = "")
  }
}

f <- file.path(results_dir, "tauR.txt")
writeLines(c("example,tauR", sprintf("%s,%.10f", names(tauR), tauR)), f)
cat("Wrote ", f, "\n", sep = "")

f <- file.path(results_dir, "one_run.txt")
writeLines(out_lines, f)
cat("Wrote ", f, "\n", sep = "")
cat(paste(out_lines, collapse = "\n"), "\n")
