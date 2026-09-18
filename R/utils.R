get_ci_from_ests = function(tau.hat, tau.se, alpha = 0.05) {
  z = stats::qnorm(alpha / 2, lower.tail = FALSE)
  c(tau.hat - z * tau.se, tau.hat + z * tau.se)
}

model_mat_grf = function(df) {
  mm = stats::model.matrix(~ . - 1, data = as.data.frame(df))
  storage.mode(mm) = "double"
  mm
}

prep_covariates = function(covariates) {
  X = as.data.frame(covariates)
  if (is.null(colnames(X))) colnames(X) = paste0("X", seq_len(ncol(X)))
  for (j in seq_along(X)) {
    xj = X[[j]]
    if (is.logical(xj) || is.character(xj)) { X[[j]] = as.factor(xj); next }
    if (is.numeric(xj) || is.integer(xj)) {
      if (length(unique(stats::na.omit(xj))) <= 2) X[[j]] = factor(xj)
    }
  }
  X
}

validate_inputs = function(response, treatment, covariates, adj_sets) {
  stopifnot(length(response) == length(treatment),
            nrow(as.data.frame(covariates)) == length(response),
            all(treatment %in% c(0, 1)),
            is.list(adj_sets), length(adj_sets) >= 2)
  nm = colnames(as.data.frame(covariates))
  bad = setdiff(unique(unlist(adj_sets)), nm)
  if (length(bad) > 0) stop("Variables not in covariates: ", paste(bad, collapse = ", "))
  if (length(Reduce(intersect, adj_sets)) == 0)
    stop("Adjustment sets must have a non-empty intersection.")
}

clip01 = function(p, eps = 1e-2) pmin(pmax(as.numeric(p), eps), 1 - eps)

safe_qsolve = function(M, b, ridge = 1e-8) {
  M = as.matrix(M); b = as.numeric(b)
  if (length(b) == 0) return(numeric(0))
  out = tryCatch(qr.solve(M, b), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out)))
    out = as.numeric(MASS::ginv(M + diag(ridge, ncol(M))) %*% b)
  out
}

solve_nu_protected = function(M, rhs, K, d, nu_regularize = 0.1) {
  M = as.matrix(M); p = ncol(M); rhs = as.numeric(rhs)
  if (p == 0) return(numeric(0))
  target = c(rep(1/K, K - 1), rep(0, d))
  if (nu_regularize <= 0) return(safe_qsolve(M, rhs))
  ridge = nu_regularize * max(sum(diag(M)) / p, 1e-6)
  out = tryCatch(as.numeric(solve(M + diag(ridge, p), rhs + ridge * target)),
                 error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) out = target
  out
}

solve_lam = function(G, reltol = 1e-10, maxit = 5000) {
  G = as.matrix(G)
  if (!is.numeric(G)) storage.mode(G) = "double"
  p = ncol(G)
  if (p == 0) return(list(lambda = numeric(0), weights = rep(1, nrow(G)), convergence = 0))
  sc = apply(G, 2, stats::sd); sc[!is.finite(sc) | sc <= 0] = 1
  Gs = sweep(G, 2, sc, "/")
  obj = function(par) { eta = drop(Gs %*% par); m = max(eta); mean(exp(eta - m)) * exp(m) }
  gr  = function(par) { eta = drop(Gs %*% par); m = max(eta); as.numeric(colMeans(Gs * exp(eta - m)) * exp(m)) }
  opt = stats::optim(rep(0, p), obj, gr, method = "BFGS",
                     control = list(reltol = reltol, maxit = maxit))
  lambda = opt$par / sc
  w_raw = exp(drop(G %*% lambda))
  list(lambda = as.numeric(lambda), weights = as.numeric(w_raw / mean(w_raw)),
       convergence = opt$convergence)
}

apply_tilt = function(G, lam) {
  e = drop(as.matrix(G) %*% lam)
  v = if (max(e) > 300) exp(e - max(e)) else exp(e)
  list(weights = v / mean(v))
}

layer1_learners = function(seed = 123, num.trees = 400) list(
  fit_outcome        = function(x, y) grf::regression_forest(model_mat_grf(x), y,
                          num.trees = num.trees, seed = seed),
  predict_outcome    = function(fit, x) as.numeric(stats::predict(fit, model_mat_grf(x))$predictions),
  predict_outcome_oob    = function(fit) as.numeric(stats::predict(fit)$predictions),
  predict_propensity_oob = function(fit) {
    p = stats::predict(fit)$predictions
    if (is.matrix(p)) as.numeric(p[, ncol(p)]) else as.numeric(p) },
  fit_propensity     = function(x, a) grf::probability_forest(model_mat_grf(x),
                          factor(a, levels = c(0, 1)), num.trees = num.trees, seed = seed),
  predict_propensity = function(fit, x) {
    p = stats::predict(fit, model_mat_grf(x))$predictions
    if (is.matrix(p)) as.numeric(p[, ncol(p)]) else as.numeric(p) })

layer2_learners = function(seed = 123, num.trees = 2000) list(
  fit     = function(x, y) grf::regression_forest(model_mat_grf(x), y,
                num.trees = num.trees, seed = seed),
  predict = function(fit, x) as.numeric(stats::predict(fit, model_mat_grf(x))$predictions),
  predict_oob = function(fit) as.numeric(stats::predict(fit)$predictions))

make_folds = function(n, num_folds, seed = 123) {
  withr::with_seed(seed,
    sample(rep(seq_len(num_folds), length.out = n), size = n, replace = FALSE))
}

build_protect_block = function(X, common, protect_vars = character(0),
                               protect_fun = NULL, orthonormalize = FALSE) {
  n = nrow(X)
  blocks = list()

  protect_vars = unique(protect_vars)
  if (length(protect_vars) > 0) {
    bad = setdiff(protect_vars, common)
    if (length(bad) > 0)
      stop("Protected variables must lie in the intersection of the adjustment sets: ",
           paste(bad, collapse = ", "))
    blocks$vars = data.matrix(X[, protect_vars, drop = FALSE])
  }
  if (!is.null(protect_fun)) {
    if (!is.function(protect_fun)) stop("protect_fun must be a function.")
    fx = as.matrix(protect_fun(X[, common, drop = FALSE]))
    storage.mode(fx) = "double"
    if (nrow(fx) != n)
      stop("protect_fun must return one row per observation (got ", nrow(fx),
           ", expected ", n, ").")
    if (is.null(colnames(fx)))
      colnames(fx) = if (ncol(fx) == 1) "f" else paste0("f", seq_len(ncol(fx)))
    blocks$fun = fx
  }

  if (length(blocks) == 0) {
    empty = matrix(NA_real_, n, 0)
    return(list(F = empty, raw = empty, center = numeric(0),
                names = character(0), d = 0))
  }

  f_raw = do.call(cbind, blocks)
  storage.mode(f_raw) = "double"
  if (anyNA(f_raw)) stop("Protected columns must not contain missing values.")
  center = colMeans(f_raw)
  Fmat = sweep(f_raw, 2, center, "-")

  if (orthonormalize && ncol(Fmat) > 1) {
    P = matrix(0, n, 0); keep = integer(0)
    nrm = function(u) sqrt(mean(u^2))
    for (j in seq_len(ncol(Fmat))) {
      v0 = Fmat[, j] - mean(Fmat[, j]); s0 = nrm(v0)
      if (!is.finite(s0) || s0 <= 0) next
      v = v0
      if (ncol(P) > 0) for (pass in 1:2) {
        v = as.numeric(v - P %*% (crossprod(P, v) / n)); v = v - mean(v)
      }
      s = nrm(v)
      if (is.finite(s) && s > 1e-8 * s0) { P = cbind(P, v / s); keep = c(keep, j) }
    }
    colnames(P) = colnames(Fmat)[keep]
    Fmat = P
  }

  list(F = Fmat, raw = f_raw, center = as.numeric(center),
       names = colnames(Fmat), d = ncol(Fmat))
}
