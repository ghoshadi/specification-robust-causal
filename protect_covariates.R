# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# Algorithm A.1: specification_robust() of main.R with additional covariate
# constraints, protecting a known R^d-valued function f of the common
# covariates X_{S_cap}.  The additions are
#
#   step 5    G -> [g, F] with F = f(X_common) - P_n[f(X_common)], so the tilt
#             balances the protected moments as well as the disagreements
#   step 9    the ridge target of the nu system is extended by d zeros
#   step 10   the bias correction and eta carry the protected term F nu_protect
#   step 11   the influence function carries -(w - 1) F nu_protect
#
# f is a known function of X_common, so the protected block needs no
# out-of-bag treatment.  With protect_vars = character(0) and
# protect_fun = NULL we have d = 0 and the procedure reduces exactly to
# specification_robust().
# ==============================================================================

# ===========================================================================
# Step 5: the protected block  f(X_{S_cap}) - P_n[f(X_{S_cap})]
# ===========================================================================

build_protect_block <- function(X, common, protect_vars = character(0),
                                protect_fun = NULL, orthonormalize = FALSE) {
  n <- nrow(X)
  blocks <- list()

  protect_vars <- unique(protect_vars)
  if (length(protect_vars) > 0) {
    bad <- setdiff(protect_vars, common)
    if (length(bad) > 0)
      stop("Protected variables must lie in the intersection of the adjustment sets: ",
           paste(bad, collapse = ", "))
    blocks$vars <- data.matrix(X[, protect_vars, drop = FALSE])
  }
  if (!is.null(protect_fun)) {
    if (!is.function(protect_fun)) stop("protect_fun must be a function.")
    fx <- as.matrix(protect_fun(X[, common, drop = FALSE]))
    storage.mode(fx) <- "double"
    if (nrow(fx) != n)
      stop("protect_fun must return one row per observation (got ", nrow(fx),
           ", expected ", n, ").")
    if (is.null(colnames(fx)))
      colnames(fx) <- if (ncol(fx) == 1) "f" else paste0("f", seq_len(ncol(fx)))
    blocks$fun <- fx
  }

  if (length(blocks) == 0) {
    empty <- matrix(NA_real_, n, 0)
    return(list(F = empty, raw = empty, center = numeric(0),
                names = character(0), d = 0L))
  }

  f_raw <- do.call(cbind, blocks)
  storage.mode(f_raw) <- "double"
  if (anyNA(f_raw)) stop("Protected columns must not contain missing values.")
  center <- colMeans(f_raw)
  Fmat <- sweep(f_raw, 2, center, "-")

  if (orthonormalize && ncol(Fmat) > 1) {
    P <- matrix(0, n, 0); keep <- integer(0)
    nrm <- function(u) sqrt(mean(u^2))
    for (j in seq_len(ncol(Fmat))) {
      v0 <- Fmat[, j] - mean(Fmat[, j]); s0 <- nrm(v0)
      if (!is.finite(s0) || s0 <= 0) next
      v <- v0
      if (ncol(P) > 0) for (pass in 1:2) {
        v <- as.numeric(v - P %*% (crossprod(P, v) / n)); v <- v - mean(v)
      }
      s <- nrm(v)
      if (is.finite(s) && s > 1e-8 * s0) { P <- cbind(P, v / s); keep <- c(keep, j) }
    }
    colnames(P) <- colnames(Fmat)[keep]
    Fmat <- P
  }

  list(F = Fmat, raw = f_raw, center = as.numeric(center),
       names = colnames(Fmat), d = ncol(Fmat))
}

# ===========================================================================
# Step 9: the nu system with the protected block
# ===========================================================================

solve_nu_protected <- function(M, rhs, K, d, nu_regularize = 0.1) {
  M <- as.matrix(M); p <- ncol(M); rhs <- as.numeric(rhs)
  if (p == 0) return(numeric(0))
  target <- c(rep(1/K, K - 1), rep(0, d))
  if (nu_regularize <= 0) return(safe_qsolve(M, rhs))
  ridge <- nu_regularize * max(sum(diag(M)) / p, 1e-6)
  out <- tryCatch(as.numeric(solve(M + diag(ridge, p), rhs + ridge * target)),
                  error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) out <- target
  out
}

# ===========================================================================

specification_robust_protect <- function(response,
                                 treatment,
                                 covariates,
                                 adj_sets,
                                 protect_vars     = character(0),
                                 protect_fun      = NULL,
                                 protect_orthonormalize = FALSE,
                                 ref_index        = 1L,
                                 verbose          = TRUE,
                                 alpha            = 0.05,
                                 seed             = 42,
                                 num_folds        = 2,
                                 propensity_clip  = 0,
                                 num_trees        = 400,
                                 nu_regularize    = 1e-6,
                                 bias_corr        = TRUE,
                                 aipw_trim        = 0.01,
                                 weight_trim      = 0.05,
                                 return_full      = FALSE) {

  stopifnot(length(bias_corr) == 1L, is.logical(bias_corr), !is.na(bias_corr))
  stopifnot(length(weight_trim) == 1L, is.finite(weight_trim),
            weight_trim >= 0, weight_trim < 1)
  stopifnot(length(aipw_trim) == 1L, is.finite(aipw_trim),
            aipw_trim >= 0, aipw_trim < 0.5)
  num_folds <- as.integer(num_folds)
  stopifnot(length(num_folds) == 1L, !is.na(num_folds), num_folds >= 2)

  validate_inputs(response, treatment, covariates, adj_sets)
  K <- length(adj_sets)
  if (is.null(ref_index)) ref_index <- 1L
  if (is.character(ref_index))
    ref_index <- switch(ref_index,
      first    = 1L,
      largest  = which.max(lengths(adj_sets)),
      smallest = which.min(lengths(adj_sets)),
      stop("ref_index must be a numeric index, or one of ",
           "\"first\", \"largest\", \"smallest\""))
  ref_index <- as.integer(ref_index)
  stopifnot(length(ref_index) == 1L, ref_index >= 1, ref_index <= K)
  if (ref_index != 1) {
    perm <- c(ref_index:K, if (ref_index > 1) 1:(ref_index - 1))
    adj_sets <- adj_sets[perm]
    if (verbose) cat("Reference: original S", ref_index, " moved to index 1.\n", sep = "")
  }

  y <- as.numeric(response); a <- as.integer(treatment)
  X <- prep_covariates(covariates)
  learners  <- layer1_learners(seed = seed, num.trees = num_trees)
  g_learner <- layer2_learners(seed = seed, num.trees = num_trees)

  n <- length(y); common <- Reduce(intersect, adj_sets)

  Fb <- build_protect_block(X, common, protect_vars, protect_fun,
                            orthonormalize = protect_orthonormalize)
  pd <- Fb$d                       # protected dimension
  fold_id <- make_folds(n, num_folds, seed)
  cycles <- cbind(eval_fold = seq_len(num_folds),
                  n_eval  = as.integer(table(factor(fold_id, seq_len(num_folds)))),
                  n_train = n - as.integer(table(factor(fold_id, seq_len(num_folds)))))

  if (verbose)
    cat(sprintf("K=%d  common={%s}  n=%d\n  %d folds  trees=%d  protected={%s} d=%d\n\n",
        K, paste(common, collapse = ", "), n, num_folds, num_trees,
        if (pd == 0) "none" else paste(Fb$names, collapse = ", "), pd))

  FD <- vector("list", num_folds); nu_sys <- vector("list", num_folds)
  aipw_oof <- tau_hat_oof <- matrix(NA_real_, n, K)
  g_hat_oof <- matrix(NA_real_, n, K - 1 + pd)
  m1_oof <- weights_oof <- w_raw_oof <- eta_oof <- psi_oof <- bc_oof <- rep(NA_real_, n)
  lambda_by_fold <- matrix(NA_real_, num_folds, K - 1 + pd)
  balance_by_fold <- balance_train_by_fold <-
    matrix(NA_real_, num_folds, K - 1 + pd)
  nu_protect_by_fold <- matrix(NA_real_, num_folds, pd)
  nu_by_fold <- reweighted_by_fold <- matrix(NA_real_, num_folds, K)
  fold_estimate <- fold_se <- fold_ess <- rep(NA_real_, num_folds)

  for (o in seq_len(num_folds)) {
    TR <- which(fold_id != o); EV <- which(fold_id == o)
    Xtr <- X[TR, , drop = FALSE]; Xev <- X[EV, , drop = FALSE]
    XcTR <- Xtr[, common, drop = FALSE]; XcEV <- Xev[, common, drop = FALSE]
    n_ev <- length(EV)

    P <- Tt <- R1 <- R0 <- vector("list", K); ED <- vector("list", K)
    for (k in seq_len(K)) {
      vk  <- adj_sets[[k]]
      xtr <- Xtr[, vk, drop = FALSE]; itr <- TR
      xev <- Xev[, vk, drop = FALSE]; iev <- EV

      e_fit   <- learners$fit_propensity(xtr, a[itr])
      mu0_fit <- learners$fit_outcome(xtr[a[itr] == 0, , drop = FALSE], y[itr][a[itr] == 0])
      mu1_fit <- learners$fit_outcome(xtr[a[itr] == 1, , drop = FALSE], y[itr][a[itr] == 1])

      i1 <- which(a[itr] == 1); i0 <- which(a[itr] == 0)
      u1_tr <- u0_tr <- numeric(length(itr))
      u1_tr[i1] <- learners$predict_outcome_oob(mu1_fit)
      u1_tr[i0] <- learners$predict_outcome(mu1_fit, xtr[i0, , drop = FALSE])
      u0_tr[i0] <- learners$predict_outcome_oob(mu0_fit)
      u0_tr[i1] <- learners$predict_outcome(mu0_fit, xtr[i1, , drop = FALSE])
      # the layer-2 response: the layer-1 arm predictions, not psi
      R1[[k]] <- u1_tr; R0[[k]] <- u0_tr

      # out-of-sample on EVAL
      eraw <- learners$predict_propensity(e_fit, xev)
      e_ev <- clip01(eraw, eps = propensity_clip)
      u0   <- learners$predict_outcome(mu0_fit, xev)
      u1   <- learners$predict_outcome(mu1_fit, xev)
      ED[[k]] <- c(min = min(eraw), max = max(eraw),
                   n_clip = sum(eraw < propensity_clip | eraw > 1 - propensity_clip))
      Tt[[k]] <- u1 - u0
      P[[k]]  <- (u1 - u0) + a[iev]/e_ev*(y[iev] - u1) -
                 (1 - a[iev])/(1 - e_ev)*(y[iev] - u0)
    }
    psi_ev <- do.call(cbind, P); tau_ev <- do.call(cbind, Tt)
    r1_tr  <- do.call(cbind, R1); r0_tr <- do.call(cbind, R0)

    if (aipw_trim > 0) {
      ak <- rep(TRUE, nrow(psi_ev))
      for (k in seq_len(K)) {
        qq <- stats::quantile(psi_ev[, k], c(aipw_trim, 1 - aipw_trim))
        ak <- ak & psi_ev[, k] >= qq[1] & psi_ev[, k] <= qq[2]
      }
      EV     <- EV[ak]
      Xev    <- Xev[ak, , drop = FALSE]; XcEV <- Xev[, common, drop = FALSE]
      psi_ev <- psi_ev[ak, , drop = FALSE]; tau_ev <- tau_ev[ak, , drop = FALSE]
      n_ev   <- length(EV)
    }

    Me <- Mt <- vector("list", K)
    for (k in seq_len(K)) {
      q1 <- g_learner$fit(XcTR, r1_tr[, k])
      q0 <- g_learner$fit(XcTR, r0_tr[, k])
      Me[[k]] <- g_learner$predict(q1, XcEV) - g_learner$predict(q0, XcEV)
      # the same projection on the TRAINING rows, out of bag
      Mt[[k]] <- g_learner$predict_oob(q1) - g_learner$predict_oob(q0)
    }
    m_ev <- do.call(cbind, Me); m_tr <- do.call(cbind, Mt)
    for (k in seq_len(K))
      if (setequal(adj_sets[[k]], common)) {
        m_ev[, k] <- tau_ev[, k]
        m_tr[, k] <- r1_tr[, k] - r0_tr[, k]
      }
    G_ev <- do.call(cbind, lapply(2:K, function(k) m_ev[, 1] - m_ev[, k]))
    G_tr <- do.call(cbind, lapply(2:K, function(k) m_tr[, 1] - m_tr[, k]))

    F_ev  <- Fb$F[EV, , drop = FALSE]
    F_tr  <- Fb$F[TR, , drop = FALSE]
    Gd_ev <- G_ev
    G_ev  <- cbind(G_ev, F_ev)
    G_tr  <- cbind(G_tr, F_tr)

    G_fit <- G_tr
    w_obj <- solve_lam(as.matrix(G_fit))
    lam   <- as.numeric(w_obj$lambda)
    w_raw_full <- as.numeric(apply_tilt(G_ev, lam)$weights)
    w_tr   <- as.numeric(apply_tilt(G_tr, lam)$weights)
    sd_tr  <- apply(G_tr, 2, stats::sd); sd_tr[sd_tr <= 0] <- 1
    bal_tr <- colMeans(w_tr * G_tr)/sd_tr
    sd_ev <- apply(G_ev, 2, stats::sd); sd_ev[sd_ev <= 0] <- 1
    bal_ev <- colMeans(w_raw_full * G_ev)/sd_ev

    ev_all <- EV
    aipw_oof[ev_all, ] <- psi_ev; tau_hat_oof[ev_all, ] <- tau_ev
    g_hat_oof[ev_all, ] <- G_ev; m1_oof[ev_all] <- m_ev[, 1]

    if (weight_trim > 0) {
      keep <- which(w_raw_full <= stats::quantile(w_raw_full, 1 - weight_trim))
      w <- w_raw_full[keep]/mean(w_raw_full[keep])
    } else {
      keep <- seq_along(w_raw_full)
      w <- w_raw_full                      # apply_tilt already gives mean 1
    }
    ess_used <- sum(w)^2/sum(w^2)
    psi_ev <- psi_ev[keep, , drop = FALSE]; tau_ev <- tau_ev[keep, , drop = FALSE]
    m_ev   <- m_ev[keep, , drop = FALSE];   G_ev   <- G_ev[keep, , drop = FALSE]
    Gd_ev  <- Gd_ev[keep, , drop = FALSE];  F_ev   <- F_ev[keep, , drop = FALSE]
    w_raw  <- w_raw_full[keep]; n_ev <- length(keep)

    tau_R1 <- mean(w * psi_ev[, 1])
    M   <- crossprod(G_ev, G_ev * w)/n_ev
    rhs <- colMeans(w * G_ev * (m_ev[, 1] - tau_R1))
    nu_all <- solve_nu_protected(M, rhs, K, pd, nu_regularize)
    nu2 <- nu_all[seq_len(K - 1)]
    nu_protect <- if (pd > 0) nu_all[K - 1 + seq_len(pd)] else numeric(0)
    nu  <- c(1 - sum(nu2), nu2)
    protect_ev <- if (pd > 0) drop(F_ev %*% nu_protect) else rep(0, n_ev)
    nu_sys[[o]] <- list(nu = nu, M = M, rhs = rhs, tau_R1 = tau_R1, n = n_ev)

    m_nu  <- drop(m_ev %*% nu)
    phi   <- m_nu - tau_R1 - protect_ev
    dpsi  <- do.call(cbind, lapply(2:K, function(k) psi_ev[, 1] - psi_ev[, k]))
    R     <- drop((dpsi - Gd_ev) %*% lam[seq_len(K - 1)])
    if (bias_corr) {
      bc    <- w * R * phi
      bc_if <- w * R
    } else {
      bc    <- rep(0, n_ev)
      bc_if <- rep(0, n_ev)
    }
    eta   <- w * (drop(psi_ev %*% nu) - protect_ev) + bc
    est   <- mean(eta)
    if_bc <- bc_if * (m_nu - protect_ev - est)
    psi_i <- w * (drop(psi_ev %*% nu) - est) - (w - 1) * protect_ev + if_bc

    ev <- EV[keep]
    weights_oof[ev] <- w
    w_raw_oof[EV] <- w_raw_full
    eta_oof[ev] <- eta; psi_oof[ev] <- psi_i; bc_oof[ev] <- bc
    lambda_by_fold[o, ] <- lam; nu_by_fold[o, ] <- nu
    if (pd > 0) nu_protect_by_fold[o, ] <- nu_protect
    balance_by_fold[o, ] <- bal_ev; balance_train_by_fold[o, ] <- bal_tr
    reweighted_by_fold[o, ] <- colMeans(w * psi_ev)
    fold_estimate[o] <- est
    fold_se[o] <- sqrt(stats::var(psi_i)/n_ev)
    fold_ess[o] <- ess_used

    FD[[o]] <- list(
      idx = list(ev = ev), n_ev = n_ev, n_ca = n_ev,
      psi_ev = psi_ev, tau_ev = tau_ev, m_ev = m_ev, m_ca = m_ev,
      G_ev = G_ev, G_ca = G_ev,
      w_ev = w, w_ev_raw = w_raw, ess_ev = ess_used, keep = keep,
      mom_post_trim = colMeans(w * G_ev)/apply(G_ev, 2, sd),
      lambda = lam, lambda_conv = w_obj$convergence,
      mom_train = bal_tr,          # exact: the fold lambda was solved on
      mom_ev    = bal_ev,          # the held-out residual, the honest number
      G_tr = G_tr, w_tr = w_tr,
      e_diag  = do.call(rbind, ED),
      nu = nu, nu_protect = nu_protect, protect_ev = protect_ev,
      F_ev = F_ev, Gd_ev = Gd_ev,
      phi = phi, R = R, bc = bc, eta = eta, psi_if = psi_i,
      tau_R1 = tau_R1)

    if (verbose) {
      ef <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
      cat(sprintf("Fold %d of %d  (EVAL = fold %d, n = %d; nuisances fitted on the other %d fold(s), n = %d)\n",
                  o, num_folds, o, length(EV), num_folds - 1L, length(TR)))
      cat("  lambda :", paste(sprintf("%.3e", lam), collapse = " "),
          sprintf("  (solved on the other %d fold(s))", num_folds - 1L), "\n")
      cat(sprintf("  balance: %.2e on the fold lambda came from, %.2e on this one\n",
                  max(abs(bal_tr)), max(abs(bal_ev))))
      cat("  nu     :", paste(sprintf("%8.4f", nu), collapse = " "), "\n")
      if (pd > 0)
        cat("  nu_prot:", paste(sprintf("%8.4f", nu_protect), collapse = " "), "\n")
      cat(sprintf("  cond(M) %.1f   fold est %.4f   fold s.e. %.4f   ESS %.0f/%d\n\n",
                  max(ef)/min(ef), est, fold_se[o], ess_used, n_ev))
    }
  }

  used <- !is.na(eta_oof); n_used <- sum(used)
  estimate <- mean(eta_oof[used]); se <- sqrt(stats::var(psi_oof[used])/n_used)
  ci <- get_ci_from_ests(estimate, se, alpha)
  a_ok <- stats::complete.cases(aipw_oof); n_a <- sum(a_ok)
  aipw_estimates <- colMeans(aipw_oof[a_ok, , drop = FALSE])
  aipw_se <- apply(aipw_oof[a_ok, , drop = FALSE], 2, stats::sd)/sqrt(n_a)
  reweighted_estimates <- colMeans(weights_oof[used] * aipw_oof[used, , drop = FALSE])
  aipw_cis <- lapply(seq_len(K), function(k)
    get_ci_from_ests(aipw_estimates[k], aipw_se[k], alpha))
  convex_hull <- c(min(sapply(aipw_cis, `[`, 1)), max(sapply(aipw_cis, `[`, 2)))
  nev <- sapply(FD, `[[`, "n_ev"); nev <- nev/sum(nev)
  nu_effective <- as.numeric(colSums(nu_by_fold * nev))

  colnames(tau_hat_oof) <- paste0("tau_hat_", seq_len(K))
  colnames(aipw_oof) <- paste0("aipw_", seq_len(K))
  colnames(g_hat_oof) <- c(paste0("g_", 2:K),
    if (pd > 0) paste0("protect_", Fb$names) else character(0))

  # What the protection achieved: P_n[w f_j] against P_n[f_j], in sd(f_j).
  # The weights here are the fitted weight function before the trim, which is
  # what the tilt balanced; the trim is a separate, later choice.
  protected_summary <- data.frame()
  if (pd > 0) {
    wok <- !is.na(w_raw_oof)
    fw  <- as.numeric(colSums(w_raw_oof[wok] * Fb$raw[wok, , drop = FALSE]) /
                      sum(w_raw_oof[wok]))
    fsd <- apply(Fb$raw, 2, stats::sd); fsd[!is.finite(fsd) | fsd <= 0] <- 1
    protected_summary <- data.frame(
      variable = colnames(Fb$raw), original_mean = Fb$center,
      weighted_mean = fw, std_difference = (fw - Fb$center)/fsd,
      row.names = NULL)
  }

  if (verbose) {
    cat("================ FINAL SUMMARY ================\n")
    cat("AIPW estimates      :", paste(round(aipw_estimates, 4), collapse = ", "), "\n")
    cat("AIPW s.e.           :", paste(round(aipw_se, 4), collapse = ", "), "\n")
    cat("Convex hull CI      : [", round(convex_hull[1], 4), ", ", round(convex_hull[2], 4),
        "]  (width = ", round(diff(convex_hull), 4), ")\n", sep = "")
    cat("nu (effective)      :", paste(sprintf("%8.4f", nu_effective), collapse = " "), "\n")
    cat("Reweighted estimates:", paste(round(reweighted_estimates, 4), collapse = ", "), "\n")
    cat("Spec-robust est.    :", round(estimate, 4), "\n")
    cat("Spec-robust s.e.    :", round(se, 4), "\n")
    cat(sprintf("Rows used           : %d of %d\n", n_used, n))
    cat(sprintf("Folds               : %d   (per-fold lambda and nu)\n", num_folds))
    cat("Spec-robust CI      : [", round(ci[1], 4), ", ", round(ci[2], 4),
        "]  (width = ", round(diff(ci), 4), ")\n", sep = "")
    cat(sprintf("Width reduction     : %.1f%%\n", 100*(1 - diff(ci)/diff(convex_hull))))
    if (pd > 0) {
      cat("---- protected moments (weighted mean vs original) ----\n")
      for (j in seq_len(nrow(protected_summary)))
        cat(sprintf("  %-18s %10.4f -> %10.4f   (%+.2e sd)\n",
            protected_summary$variable[j], protected_summary$original_mean[j],
            protected_summary$weighted_mean[j], protected_summary$std_difference[j]))
    }
    cat("================================================\n")
  }

  out <- list(
    estimate = estimate, se = se, ci = ci,
    aipw_estimates = aipw_estimates, aipw_se = aipw_se,
    reweighted_estimates = reweighted_estimates, convex_hull_ci = convex_hull,
    nu = nu_effective, nu_by_fold = nu_by_fold, lambda_by_fold = lambda_by_fold,
    balance_by_fold = balance_by_fold,
    balance_train_by_fold = balance_train_by_fold,
    nu_protect_by_fold = nu_protect_by_fold, protected_summary = protected_summary,
    protect_names = Fb$names, protect_d = pd,
    reweighted_by_fold = reweighted_by_fold,
    fold_estimate = fold_estimate, fold_se = fold_se, fold_ess = fold_ess,
    weights = weights_oof, weights_raw = w_raw_oof,
    fold_id = fold_id, n_used = n_used,
    common_covariates = common, adj_sets = adj_sets, ref_index = ref_index, alpha = alpha,
    settings = list(n_folds = num_folds,
      g_construction = "diff_of_proj", lambda_source = "train",
      layer1_oob = TRUE, exact_common = TRUE, nu_use_m1 = TRUE,
      nu_source = "eval_fold",
      seed = seed,
      propensity_clip = propensity_clip, num_trees = num_trees,
      nu_regularize = nu_regularize, bias_corr = bias_corr,
      aipw_trim = aipw_trim, weight_trim = weight_trim,
      protect_vars = protect_vars, protect_d = pd,
      protect_orthonormalize = protect_orthonormalize))

  if (return_full) out$full <- list(
    folds = FD, nu_sys = nu_sys, cycles = cycles, K = K,
    tau_hat_oof = tau_hat_oof, aipw_oof = aipw_oof, g_hat_oof = g_hat_oof,
    m1_oof = m1_oof, w_raw_oof = w_raw_oof, eta_oof = eta_oof, Fb = Fb,
    psi_oof = psi_oof, bc_oof = bc_oof)
  out
}
