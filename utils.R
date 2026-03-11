# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# This contains all helper functions used in the main file/examples.
# ==============================================================================

suppressPackageStartupMessages({
  library(ranger)
  library(grf)
  library(dplyr)
  library(plotrix)
  library(parallel)
  library(pbapply)
  library(pbmcapply)
  library(foreign)
  library(MASS)
})

# ===========================================================================
# Plots: Function for comparing histograms side by side
# ===========================================================================

compare_hists <- function(data, weights, breaks = 20, 
                          hcol = "orange", hcol2 = "skyblue2", 
                          xlab = NULL, xlim = NULL,
                          mainOG = "Original population", 
                          mainAR = "New target population", 
                          align = "vertical") {
  if(is.null(xlim)) xlim <- range(data)
  if(align == "vertical"){
    par(mar = c(5, 5, 4, 2) + 0.1)
    par(oma = rep(0.1, 4)) #margins below, left, top, right 
  }  
  else{
    par(mar = c(5, 5, 4, 2) + 0.1)
    par(oma = rep(0.1, 4))
  } 
  hist(data, breaks = breaks, col = hcol, main = mainOG, xlim = xlim, 
       xlab = xlab, ylab = "Frequency")
  box()
  hist_data <- hist(data, breaks = breaks, plot = FALSE)
  counts <- hist_data$counts
  breaks <- hist_data$breaks
  midpoints <- hist_data$mids
  weighted_counts <- sapply(seq_along(counts), function(i) {
    idx <- which(data >= breaks[i] & data < breaks[i + 1])
    sum(weights[idx])
  })
  plot(midpoints, weighted_counts, type = "n", main = mainAR, 
       xlab = xlab, ylab = "Frequency", xlim = xlim)
  rect(breaks[-length(breaks)], 0, breaks[-1], weighted_counts, col = hcol2)
}

# ===========================================================================
# Plots: Function for comparing histograms overlay (as in the paper)
# ===========================================================================

compare_hists_overlay <- function(data, weights, breaks = 25, 
                                  hcol = "orange", 
                                  hcol2 = "skyblue2", 
                                  xlim = NULL, ylim = NULL,
                                  xlab = NULL, main = NULL,
                                  legend.text = c("original population", 
                                                  "new target population"),
                                  legend.pos = "topleft",
                                  draw.curve = F) {
  h1 <- hist(data, breaks = breaks, plot = FALSE)
  bin_widths <- diff(h1$breaks)
  h1_density <- h1$counts / sum(h1$counts) / bin_widths
  
  weighted_counts <- sapply(seq_along(h1$counts), function(i) {
    idx <- which(data >= h1$breaks[i] & data < h1$breaks[i + 1])
    sum(weights[idx])
  })
  
  weighted_density <- weighted_counts / sum(weights) / bin_widths
  if(is.null(xlim)) xlim <- range(data)
  if(is.null(ylim)) ylim <- range(c(h1_density, weighted_density))
  
  plot(h1$mids, h1_density, type = "n", 
       xlim = xlim, ylim = ylim, 
       xlab = xlab, ylab = "Density", 
       main = main)
  legend(legend.pos, legend = legend.text, 
         fill = c(adjustcolor(hcol, alpha.f = 0.6), 
                  adjustcolor(hcol2, alpha.f = 0.6)), 
         col = c(adjustcolor(hcol, alpha.f = 0.6), 
                 adjustcolor(hcol2, alpha.f = 0.6)), 
         bty = "n", cex = 1)
  rect(h1$breaks[-length(h1$breaks)], 0, 
       h1$breaks[-1], h1_density, 
       col = adjustcolor(hcol, alpha.f = 0.6), 
       border = adjustcolor(hcol, alpha.f = 0.2))
  rect(h1$breaks[-length(h1$breaks)], 0, 
       h1$breaks[-1], weighted_density, 
       col = adjustcolor(hcol2, alpha.f = 0.6), 
       border = adjustcolor(hcol2, alpha.f = 0.2))
  if(draw.curve == T){
    k <- 5
    h1_density_smooth <- stats::filter(h1_density, rep(1/k, k), sides = 2)
    weighted_density_smooth <- stats::filter(weighted_density, rep(1/k, k), sides = 2)
    lines(h1$mids, h1_density_smooth+5e-4, col = hcol, lwd = 3)
    lines(h1$mids, weighted_density_smooth+5e-4, col = hcol2, lwd = 3)
  }
  box()
}

# ===========================================================================
# Various helper functions
# ===========================================================================

get_ci_from_ests <- function(tau.hat, tau.se, alpha = 0.05) {
  z <- qnorm(alpha / 2, lower.tail = FALSE)
  c(tau.hat - z * tau.se, tau.hat + z * tau.se)
}

to_mm <- function(df) {
  mm <- stats::model.matrix(~ . - 1, data = as.data.frame(df))
  storage.mode(mm) <- "double"
  mm
}

prep_covariates <- function(covariates) {
  X <- as.data.frame(covariates)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  for (j in seq_along(X)) {
    xj <- X[[j]]
    if (is.logical(xj) || is.character(xj)) { X[[j]] <- as.factor(xj); next }
    if (is.numeric(xj) || is.integer(xj)) {
      if (length(unique(stats::na.omit(xj))) <= 2) X[[j]] <- factor(xj)
    }
  }
  X
}

make_folds_3 <- function(n, seed = 123) {
  set.seed(seed)
  sample(rep(1:3, length.out = n), size = n, replace = FALSE)
}

validate_inputs <- function(response, treatment, covariates, adj_sets) {
  stopifnot(length(response) == length(treatment),
            nrow(as.data.frame(covariates)) == length(response),
            all(treatment %in% c(0, 1)),
            is.list(adj_sets), length(adj_sets) >= 2)
  nm <- colnames(as.data.frame(covariates))
  bad <- setdiff(unique(unlist(adj_sets)), nm)
  if (length(bad) > 0) stop("Variables not in covariates: ", paste(bad, collapse = ", "))
  if (length(Reduce(intersect, adj_sets)) == 0)
    stop("Adjustment sets must have a non-empty intersection.")
}

# ===========================================================================
# Truncation helpers (used for numerical stability)
# ===========================================================================

clip01 <- function(p, eps = 1e-2) pmin(pmax(as.numeric(p), eps), 1 - eps)

# Truncate weights at quantile q, renormalise to mean 1.
truncate_weights <- function(w, q = 0.99) {
  cap <- quantile(w, q); w <- pmin(w, cap)
  w <- w / mean(w)
  ess <- sum(w)^2 / sum(w^2)
  list(w = w, ess = ess)
}

# Winsorise a vector at symmetric quantiles
winsorise <- function(x, q = 0.995) {
  lo <- quantile(x, 1 - q); hi <- quantile(x, q)
  pmin(pmax(x, lo), hi)
}

# ===========================================================================
# Main nuisance learners (mu_a, propensity) using grf
# ===========================================================================

make_learners <- function(num.trees = 400, seed = 123) {
  list(
    fit_outcome = function(x, y)
      regression_forest(to_mm(x), y, num.trees = num.trees, seed = seed),
    predict_outcome = function(fit, x)
      as.numeric(predict(fit, to_mm(x))$predictions),
    fit_propensity = function(x, a)
      probability_forest(to_mm(x), factor(a, levels = c(0, 1)),
                         num.trees = num.trees, seed = seed),
    predict_propensity = function(fit, x) {
      pred <- predict(fit, to_mm(x))$predictions
      clip01(if (is.matrix(pred)) as.numeric(pred[, ncol(pred)]) else as.numeric(pred))
    }
  )
}

# ===========================================================================
# g/h learner functions (h_k=E[tau(X_{S_k}) | X_common], g_k = h_1 - h_{k+1})
# ===========================================================================

make_g_learners_grf <- function(num.trees = 400, seed = 123) {
  list(
    fit = function(x, y) regression_forest(to_mm(x), y, num.trees = num.trees, seed = seed),
    predict = function(fit, x) as.numeric(predict(fit, to_mm(x))$predictions)
  )
}

make_g_learners_ranger <- function(num.trees = 400, seed = 123) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("ranger not installed")
  list(
    fit = function(x, y) {
      df <- data.frame(y = y, x)
      ranger::ranger(y ~ ., data = df, num.trees = num.trees,
                     respect.unordered.factors = "partition", seed = seed)
    },
    predict = function(fit, x) as.numeric(predict(fit, data = as.data.frame(x))$predictions)
  )
}

make_g_learners_lm <- function() {
  list(
    fit = function(x, y) lm(y ~ ., data = data.frame(y = y, x)),
    predict = function(fit, x) as.numeric(predict(fit, newdata = as.data.frame(x)))
  )
}

# ===========================================================================
# Helpers for the affine weights (nu)
# ===========================================================================

# Ridge-stabilized linear solve (fallback for safe_qsolve)
safe_qsolve <- function(M, b, ridge = 1e-8) {
  M <- as.matrix(M); b <- as.numeric(b)
  if (length(b) == 0) return(numeric(0))
  out <- tryCatch(qr.solve(M, b), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out)))
    out <- as.numeric(MASS::ginv(M + diag(ridge, ncol(M))) %*% b)
  out
}

# Regularized solve for nu_{2:K}, shrinking toward equal weights (1/K).
#   (M + ridge*I)^{-1} (rhs + ridge*target)
# Adaptive ridge = nu_regularize * max(trace(M)/(K-1), floor).
# When g has signal: ridge << diag(M), so nearly unbiased.
# When g ~ 0: ridge dominates, nu -> 1/K (equal weights).
solve_nu <- function(M, rhs, K, nu_regularize = 0.1) {
  M <- as.matrix(M); d <- ncol(M); rhs <- as.numeric(rhs)
  if (d == 0) return(numeric(0))
  target <- rep(1 / K, d)
  if (nu_regularize <= 0) return(safe_qsolve(M, rhs))
  ridge <- nu_regularize * max(sum(diag(M)) / d, 1e-6)
  M_reg <- M + diag(ridge, d)
  rhs_reg <- rhs + ridge * target
  out <- tryCatch(as.numeric(solve(M_reg, rhs_reg)), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) out <- target
  out
}