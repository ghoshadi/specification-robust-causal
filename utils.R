# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# Helper functions for main.R, protect_covariates.R and the examples.
#
# No package is attached; contributed functions are called with their
# namespace, so grf, MASS and withr need only be installed.
# ==============================================================================
plot_hists <- function(data, weights, vars, labels = NULL, breaks = 25,
                       legend.pos = "topleft", hcol = "orange",
                       hcol2 = "skyblue2",
                       legend.text = c("original population",
                                       "new target population")) {
  if (!is.list(weights)) weights <- list(weights)
  nr <- length(weights); nc <- length(vars)
  per <- function(a, j) if (length(a) == 1L) a[[1L]] else a[[j]]

  cell <- function(i, j) {
    w <- weights[[i]]; ok <- !is.na(w)
    x <- data[[vars[j]]][ok]; w <- w[ok]
    b <- per(breaks, j)
    if (identical(b, "unit"))
      b <- seq(floor(min(x)) - 0.5, ceiling(max(x)) + 0.5, by = 1)
    h  <- hist(x, breaks = b, plot = FALSE)
    bw <- diff(h$breaks)
    wc <- as.numeric(tapply(w, cut(x, h$breaks, include.lowest = TRUE), sum))
    wc[is.na(wc)] <- 0
    list(breaks = h$breaks, xlim = range(x),
         d = list(h$counts/sum(h$counts)/bw, wc/sum(w)/bw))
  }
  cells <- lapply(seq_len(nr), function(i)
    lapply(seq_len(nc), function(j) cell(i, j)))
  ylim <- lapply(seq_len(nc), function(j)
    range(unlist(lapply(cells, function(row) row[[j]]$d))))

  op <- par(mfrow = c(nr, nc), mar = c(4.0, 4.2, 2.4, 0.8),
            mgp = c(2.3, 0.8, 0))
  on.exit(par(op))
  col <- adjustcolor(c(hcol, hcol2), alpha.f = 0.6)
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    z <- cells[[i]][[j]]; b <- z$breaks
    plot(NA, xlim = z$xlim, ylim = ylim[[j]], xlab = vars[j], ylab = "Density")
    # the row is named once, over its leftmost panel
    if (!is.null(labels) && j == 1L)
      title(main = labels[i], adj = 0, font.main = 1)
    legend(per(legend.pos, j), legend = legend.text, fill = col, col = col,
           bty = "n")
    for (k in 1:2)
      rect(b[-length(b)], 0, b[-1], z$d[[k]], col = col[k],
           border = adjustcolor(c(hcol, hcol2)[k], alpha.f = 0.2))
    box()
  }
  invisible(NULL)
}

# ===========================================================================
# Various helper functions
# ===========================================================================

get_ci_from_ests <- function(tau.hat, tau.se, alpha = 0.05) {
  z <- qnorm(alpha / 2, lower.tail = FALSE)
  c(tau.hat - z * tau.se, tau.hat + z * tau.se)
}

model_mat_grf <- function(df) {
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

# ===========================================================================
# Main nuisance learners (mu_a, propensity) using grf
# ===========================================================================

# ===========================================================================
# g/m learner functions (m_k=E[tau(X_{S_k}) | X_common], g_k = m_1 - m_k)
# ===========================================================================

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

# ===========================================================================
# Exponential tilting for the transfer weights
# ===========================================================================
# Exponential tilting: argmin_lambda mean(exp(G %*% lambda))
solve_lam <- function(G, reltol = 1e-10, maxit = 5000) {
  G <- as.matrix(G)
  if (!is.numeric(G)) storage.mode(G) <- "double"
  p <- ncol(G)
  if (p == 0) return(list(lambda = numeric(0), weights = rep(1, nrow(G)), convergence = 0))
  sc <- apply(G, 2, stats::sd); sc[!is.finite(sc) | sc <= 0] <- 1
  Gs <- sweep(G, 2, sc, "/")
  obj <- function(par) { eta <- drop(Gs %*% par); m <- max(eta); mean(exp(eta - m)) * exp(m) }
  gr  <- function(par) { eta <- drop(Gs %*% par); m <- max(eta); as.numeric(colMeans(Gs * exp(eta - m)) * exp(m)) }
  opt <- optim(rep(0, p), obj, gr, method = "BFGS",
               control = list(reltol = reltol, maxit = maxit))
  lambda <- opt$par / sc
  w_raw <- exp(drop(G %*% lambda))
  list(lambda = as.numeric(lambda), weights = as.numeric(w_raw / mean(w_raw)),
       convergence = opt$convergence)
}

# ===========================================================================
# Machinery for specification_robust()
# ===========================================================================

layer1_learners <- function(seed = 123, num.trees = 400) list(
  fit_outcome        = function(x, y) grf::regression_forest(model_mat_grf(x), y,
                          num.trees = num.trees, seed = seed),
  predict_outcome    = function(fit, x) as.numeric(predict(fit, model_mat_grf(x))$predictions),
  predict_outcome_oob    = function(fit) as.numeric(predict(fit)$predictions),
  predict_propensity_oob = function(fit) {
    p <- predict(fit)$predictions
    if (is.matrix(p)) as.numeric(p[, ncol(p)]) else as.numeric(p) },
  fit_propensity     = function(x, a) grf::probability_forest(model_mat_grf(x),
                          factor(a, levels = c(0, 1)), num.trees = num.trees, seed = seed),
  predict_propensity = function(fit, x) {
    p <- predict(fit, model_mat_grf(x))$predictions
    if (is.matrix(p)) as.numeric(p[, ncol(p)]) else as.numeric(p) })

layer2_learners <- function(seed = 123, num.trees = 2000) list(
  fit     = function(x, y) grf::regression_forest(model_mat_grf(x), y,
                num.trees = num.trees, seed = seed),
  predict = function(fit, x) as.numeric(predict(fit, model_mat_grf(x))$predictions),
  predict_oob = function(fit) as.numeric(predict(fit)$predictions))

make_folds <- function(n, num_folds, seed = 123) {
  withr::with_seed(seed,
    sample(rep(seq_len(num_folds), length.out = n), size = n, replace = FALSE))
}

apply_tilt <- function(G, lam) {
  e <- drop(as.matrix(G) %*% lam)
  v <- if (max(e) > 300) exp(e - max(e)) else exp(e)   # shift only if needed
  list(weights = v/mean(v))
}
