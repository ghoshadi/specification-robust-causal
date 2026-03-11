# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# This uses 3-fold cyclic cross-fitting with grf nuisance learners.
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

source('./utils.R')

# ===========================================================================
# Function for finding transfer weights
# ===========================================================================

# Exponential tilting: argmin_lambda mean(exp(G %*% lambda))
get_weights <- function(G, reltol = 1e-10, maxit = 5000) {
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
# Assumption-robust inference on ATE using linear models
# ===========================================================================

specification_robust_lm <- function(response, treatment, covariates, adj_sets, verbose = FALSE) {
  k <- length(adj_sets)
  fits <- lapply(adj_sets, function(adj) {
    formula <- as.formula(paste("response ~ treatment * (", 
                                paste(adj, collapse = " + "), ")"))
    if(length(adj)==1){
      df <- data.frame(response, treatment, covariates[,adj]); names(df)[3] = adj
      lm(formula, data = df)
    } else{
      lm(formula, data = data.frame(cbind(response, treatment, covariates[, adj])))
    }
  })
  
  int_adj <- Reduce(intersect, adj_sets)
  data_for_pred_0 <- data.frame(treatment = rep(0, nrow(covariates)), covariates)
  data_for_pred_1 <- data.frame(treatment = rep(1, nrow(covariates)), covariates)
  common_covariates <- as.data.frame(covariates[, int_adj, drop = FALSE])
  
  g <- lapply(seq_len(k), function(i) {
    pred_treat_0 <- predict(fits[[i]], newdata = data_for_pred_0)
    pred_treat_1 <- predict(fits[[i]], newdata = data_for_pred_1)
    cate <- pred_treat_1 - pred_treat_0
    lm(cate ~ ., data = common_covariates)$fitted.values
  })
  
  Delta_g <- do.call(cbind, lapply(seq_len(k)[-1], function(i) g[[1]] - g[[i]]))
  out <- get_weights(Delta_g)
  if (verbose) print(paste0("lambda_hat: ",round(out$lambda, 3)))
  
  rewt_est <- sapply(g, function(gi) mean(out$weights * gi))
  theta_hats <- sapply(fits, function(fit) coef(fit)["treatment"])
  
  list(rewt.estimates = rewt_est, 
       ireg.estimates = theta_hats, 
       weights = out$weights, 
       lambda = out$lambda)
}

# ===========================================================================
# Assumption-robust inference on ATE using generalized random forests
# ===========================================================================

specification_robust <- function(response,
                                 treatment,
                                 covariates,
                                 adj_sets,
                                 ref_index       = NULL,
                                 verbose         = TRUE,
                                 alpha           = 0.05,
                                 seed            = 123,
                                 # --- Numerical stability knobs ---
                                 propensity_clip = 0.01,  # clip extreme propensity scores
                                 aipw_quantile   = 0.99,  # winsorise extreme AIPW summands
                                 weight_quantile = 0.99,  # truncate extreme transfer weights
                                 nu_regularize     = 1,     # ridge shrinkage toward 1/K
                                 nu_use_h1       = TRUE,  # use projected contrasts while finding nu
                                 bc_scale        = 1,     # use bias-correction? 1 = Yes, 0 = No 
                                 bc_g_from_cal   = FALSE, # independent g for bias-correction residual
                                 psi_quantile    = 0.99,  # winsorise extreme influence function summands
                                 g_method        = "grf", # method for double conditional, supports "ranger" and "grf"
                                 learners        = NULL,
                                 return_oof      = TRUE) {
  
  validate_inputs(response, treatment, covariates, adj_sets)
  if(is.null(ref_index)) ref_index = which.max(lengths(adj_sets))
  
  K <- length(adj_sets)
  ref_index <- as.integer(ref_index)
  stopifnot(ref_index >= 1, ref_index <= K)
  if (ref_index != 1) {
    perm <- c(ref_index:K, if (ref_index > 1) 1:(ref_index - 1))
    adj_sets <- adj_sets[perm]
    if (verbose) cat("Reference index: original S", ref_index,
                     " used as the reference.\n\n", sep = "")
  }
  
  y <- as.numeric(response)
  a <- as.integer(treatment)
  X <- prep_covariates(covariates)
  if (is.null(learners)) learners <- make_learners(seed = seed)
  if(g_method == "grf") g_learner = make_g_learners_grf(seed = seed)
  else g_learner = make_g_learners_ranger(seed = seed)
  
  n <- length(y)
  common <- Reduce(intersect, adj_sets)
  fold_id <- make_folds_3(n, seed)
  cycles <- rbind(c(1, 2, 3), c(2, 3, 1), c(3, 1, 2))
  
  if (verbose) {
    # cat(sprintf("g_method=%s  nu_use_h1=%s  bc_g_from_cal=%s  nu_regularize=%.2f  bc_scale=%.2f\n",
    #             g_method, nu_use_h1, bc_g_from_cal, nu_regularize, bc_scale))
    cat(sprintf("K=%d  common={%s}  n=%d\n\n", K, paste(common, collapse=", "), n))
  }
  
  # Out-of-fold storage
  tau_hat_oof <- matrix(NA_real_, n, K)
  aipw_oof    <- matrix(NA_real_, n, K)
  m_hat_oof   <- matrix(NA_real_, n, K)
  g_hat_oof   <- matrix(NA_real_, n, K - 1)
  weights_oof <- rep(NA_real_, n)
  eta_oof     <- rep(NA_real_, n)
  psi_oof     <- rep(NA_real_, n)
  
  lambda_by_fold     <- matrix(NA_real_, 3, K - 1)
  nu_by_fold         <- matrix(NA_real_, 3, K)
  reweighted_by_fold <- matrix(NA_real_, 3, K)
  fold_estimate      <- fold_se <- fold_ess <- rep(NA_real_, 3)
  
  for (cyc in 1:3) {
    tr_idx  <- which(fold_id == cycles[cyc, 1])
    cal_idx <- which(fold_id == cycles[cyc, 2])
    ev_idx  <- which(fold_id == cycles[cyc, 3])
    
    X_tr  <- X[tr_idx, , drop = FALSE]
    X_cal <- X[cal_idx, , drop = FALSE]
    X_ev  <- X[ev_idx, , drop = FALSE]
    y_tr <- y[tr_idx]; y_cal <- y[cal_idx]; y_ev <- y[ev_idx]
    a_tr <- a[tr_idx]; a_cal <- a[cal_idx]; a_ev <- a[ev_idx]
    n_ev <- length(ev_idx)
    
    # ==========================================================
    # Fit mu_a, propensity on TRAIN.  Predict on all folds.
    # ==========================================================
    tau_tr_list <- tau_cal_list <- tau_ev_list <- vector("list", K)
    aipw_ev_list <- m_ev_list <- vector("list", K)
    
    for (k in seq_len(K)) {
      vars_k <- adj_sets[[k]]
      xk_tr  <- X_tr[, vars_k, drop = FALSE]
      xk_cal <- X_cal[, vars_k, drop = FALSE]
      xk_ev  <- X_ev[, vars_k, drop = FALSE]
      
      e_fit <- learners$fit_propensity(xk_tr, a_tr)
      e_ev  <- clip01(learners$predict_propensity(e_fit, xk_ev), eps = propensity_clip)
      
      mu0_fit <- learners$fit_outcome(xk_tr[a_tr == 0, , drop = FALSE], y_tr[a_tr == 0])
      mu1_fit <- learners$fit_outcome(xk_tr[a_tr == 1, , drop = FALSE], y_tr[a_tr == 1])
      
      mu0_tr  <- learners$predict_outcome(mu0_fit, xk_tr)
      mu1_tr  <- learners$predict_outcome(mu1_fit, xk_tr)
      mu0_cal <- learners$predict_outcome(mu0_fit, xk_cal)
      mu1_cal <- learners$predict_outcome(mu1_fit, xk_cal)
      mu0_ev  <- learners$predict_outcome(mu0_fit, xk_ev)
      mu1_ev  <- learners$predict_outcome(mu1_fit, xk_ev)
      
      tau_tr_list[[k]]  <- mu1_tr - mu0_tr
      tau_cal_list[[k]] <- mu1_cal - mu0_cal
      tau_ev_list[[k]]  <- mu1_ev - mu0_ev
      
      aipw_raw <- tau_ev_list[[k]] +
        a_ev / e_ev * (y_ev - mu1_ev) -
        (1 - a_ev) / (1 - e_ev) * (y_ev - mu0_ev)
      
      # Winsorise AIPW summands to tame heavy tails
      aipw_ev_list[[k]] <- winsorise(aipw_raw, q = aipw_quantile)
      
      # h_k = E[tau_k | X_common], trained on TRAIN, predicted on EVAL
      h_fit <- g_learner$fit(X_tr[, common, drop = FALSE], tau_tr_list[[k]])
      m_ev_list[[k]] <- g_learner$predict(h_fit, X_ev[, common, drop = FALSE])
    }
    
    # ==========================================================
    # g = E[Delta_tau | X_common], fit ONCE on TRAIN, predict on CAL and EVAL
    # ==========================================================
    g_fits_train <- lapply(2:K, function(k)
      g_learner$fit(X_tr[, common, drop = FALSE],
                    tau_tr_list[[1]] - tau_tr_list[[k]])
    )
    G_cal <- do.call(cbind, lapply(g_fits_train, function(f)
      g_learner$predict(f, X_cal[, common, drop = FALSE])))
    G_ev_from_train <- do.call(cbind, lapply(g_fits_train, function(f)
      g_learner$predict(f, X_ev[, common, drop = FALSE])))
    
    # Optional: independent g for BC residual from CAL fold
    if (bc_g_from_cal) {
      g_fits_cal <- lapply(2:K, function(k)
        g_learner$fit(X_cal[, common, drop = FALSE],
                      tau_cal_list[[1]] - tau_cal_list[[k]]))
      G_ev_for_bc <- do.call(cbind, lapply(g_fits_cal, function(f)
        g_learner$predict(f, X_ev[, common, drop = FALSE])))
    } else {
      G_ev_for_bc <- G_ev_from_train
    }
    
    # ==========================================================
    # Find Lambda (exponential tilting on CAL fold)
    # ==========================================================
    w_obj      <- get_weights(as.matrix(G_cal))
    lambda_hat <- as.numeric(w_obj$lambda)
    
    # Transfer weights on EVAL with truncation
    w_ev_raw <- as.numeric(exp(drop(as.matrix(G_ev_from_train) %*% lambda_hat)))
    w_ev_raw <- w_ev_raw / mean(w_ev_raw)
    wt <- truncate_weights(w_ev_raw, q = weight_quantile)
    w_ev <- wt$w
    
    # ==========================================================
    # Affine weights nu (with regularization)
    # ==========================================================
    tau_ev_mat  <- do.call(cbind, tau_ev_list)
    aipw_ev_mat <- do.call(cbind, aipw_ev_list)
    m_ev_mat    <- do.call(cbind, m_ev_list)
    
    tau_R1 <- mean(w_ev * aipw_ev_mat[, 1])
    G_ev   <- G_ev_from_train
    M_nu   <- crossprod(G_ev, G_ev * w_ev) / n_ev
    
    if (nu_use_h1) {
      rhs_nu <- colMeans(w_ev * G_ev * (m_ev_mat[, 1] - tau_R1))
    } else {
      rhs_nu <- colMeans(w_ev * G_ev * (tau_ev_mat[, 1] - tau_R1))
    }
    nu_2K  <- solve_nu(M_nu, rhs_nu, K, nu_regularize)
    nu_hat <- c(1 - sum(nu_2K), nu_2K)
    
    # ==========================================================
    # Bias correction summands
    # ==========================================================
    m_nu_ev <- drop(m_ev_mat %*% nu_hat)
    phi_ev  <- m_nu_ev - tau_R1
    
    delta_aipw_ev <- do.call(cbind, lapply(2:K, function(k)
      aipw_ev_mat[, 1] - aipw_ev_mat[, k]))
    
    bc_residual <- drop((delta_aipw_ev - G_ev_for_bc) %*% lambda_hat)
    bc_ev <- bc_scale * w_ev * bc_residual * phi_ev 
    
    eta_ev <- w_ev * drop(aipw_ev_mat %*% nu_hat) + bc_ev
    est_ev <- mean(eta_ev)
    
    # ==========================================================
    # Influence function
    # ==========================================================
    psi_ev <- w_ev * (drop(aipw_ev_mat %*% nu_hat) - est_ev) +
      bc_scale * w_ev * bc_residual * (m_nu_ev - est_ev)
    
    # Winsorise psi to prevent a few extreme values dominating var
    psi_ev <- winsorise(psi_ev, q = psi_quantile)
    
    # ==========================================================
    # Store
    # ==========================================================
    tau_hat_oof[ev_idx, ] <- tau_ev_mat
    aipw_oof[ev_idx, ]    <- aipw_ev_mat
    m_hat_oof[ev_idx, ]   <- m_ev_mat
    g_hat_oof[ev_idx, ]   <- G_ev_from_train
    weights_oof[ev_idx]   <- w_ev
    eta_oof[ev_idx]       <- eta_ev
    psi_oof[ev_idx]       <- psi_ev
    
    lambda_by_fold[cyc, ]     <- lambda_hat
    nu_by_fold[cyc, ]         <- nu_hat
    reweighted_by_fold[cyc, ] <- colMeans(w_ev * aipw_ev_mat)
    fold_estimate[cyc]        <- est_ev
    fold_se[cyc]              <- sqrt(stats::var(psi_ev) / n_ev)
    fold_ess[cyc]             <- wt$ess
    
    if (verbose) {
      cat(sprintf("Cycle %d (train/cal/eval = %d/%d/%d)\n",
                  cyc, cycles[cyc, 1], cycles[cyc, 2], cycles[cyc, 3]))
      cat("  lambda:", paste(round(lambda_hat, 6), collapse = ", "), "\n")
      cat("  nu    :", paste(round(nu_hat, 4), collapse = ", "), "\n")
      cat(sprintf("  fold est: %.4f   fold s.e.: %.4f   ESS: %.0f / %d",
                  est_ev, fold_se[cyc], wt$ess, n_ev))
      if (wt$ess < n_ev / 3) cat("  ** LOW ESS **")
      if (cyc > 1 && fold_se[cyc] > 3 * min(fold_se[1:(cyc-1)], na.rm = TRUE))
        cat("  ** HIGH SE **")
      cat("\n\n")
    }
  }
  
  # ================================================================
  # Aggregate estimators
  # ================================================================
  estimate <- mean(eta_oof)
  se <- sqrt(stats::var(psi_oof) / n)
  ci <- get_ci_from_ests(estimate, se, alpha)
  
  aipw_estimates <- colMeans(aipw_oof)
  aipw_se <- apply(aipw_oof, 2, stats::sd) / sqrt(n)
  reweighted_estimates <- colMeans(weights_oof * aipw_oof)
  
  # Pooled lambda/nu (for reporting only)
  sol_pool <- get_weights(g_hat_oof)
  w_pool   <- as.numeric(exp(drop(g_hat_oof %*% sol_pool$lambda)))
  w_pool   <- truncate_weights(w_pool / mean(w_pool), q = weight_quantile)$w
  tau_R1_pool <- mean(w_pool * aipw_oof[, 1])
  M_pool <- crossprod(g_hat_oof, g_hat_oof * w_pool) / n
  if (nu_use_h1) {
    rhs_pool <- colMeans(w_pool * g_hat_oof * (m_hat_oof[, 1] - tau_R1_pool))
  } else {
    rhs_pool <- colMeans(w_pool * g_hat_oof * (tau_hat_oof[, 1] - tau_R1_pool))
  }
  nu_pool <- c(1 - sum(solve_nu(M_pool, rhs_pool, K, nu_regularize)),
               solve_nu(M_pool, rhs_pool, K, nu_regularize))
  
  # Convex hull CI
  aipw_cis <- lapply(seq_len(K), function(k)
    get_ci_from_ests(aipw_estimates[k], aipw_se[k], alpha))
  convex_hull <- c(min(sapply(aipw_cis, `[`, 1)), max(sapply(aipw_cis, `[`, 2)))
  
  colnames(tau_hat_oof) <- paste0("tau_hat_", seq_len(K))
  colnames(aipw_oof) <- paste0("aipw_", seq_len(K))
  colnames(m_hat_oof) <- paste0("m_hat_", seq_len(K))
  colnames(g_hat_oof) <- paste0("g_", 2:K)
  
  if (verbose) {
    cat("================ FINAL SUMMARY ================\n")
    cat("AIPW estimates      :", paste(round(aipw_estimates, 4), collapse = ", "), "\n")
    cat("AIPW s.e.           :", paste(round(aipw_se, 4), collapse = ", "), "\n")
    cat("Convex hull CI      : [", round(convex_hull[1], 4), ", ",
        round(convex_hull[2], 4), "]  (width = ", round(diff(convex_hull), 4), ")\n", sep = "")
    cat("lambda (exp tilt)   :", paste(round(sol_pool$lambda, 6), collapse = ", "), "\n")
    cat("nu (affine weights) :", paste(round(nu_pool, 4), collapse = ", "), "\n")
    cat("Reweighted estimates:", paste(round(reweighted_estimates, 4), collapse = ", "), "\n")
    cat("Spec-robust est.    :", round(estimate, 4), "\n")
    cat("Spec-robust s.e.    :", round(se, 4), "\n")
    cat("Spec-robust CI      : [", round(ci[1], 4), ", ", round(ci[2], 4),
        "]  (width = ", round(diff(ci), 4), ")\n", sep = "")
    cat(sprintf("Width reduction     : %.0f%%\n", 100 * (1 - diff(ci) / diff(convex_hull))))
    cat("================================================\n")
  }
  
  out <- list(
    estimate = estimate, se = se, ci = ci,
    aipw_estimates = aipw_estimates, aipw_se = aipw_se,
    reweighted_estimates = reweighted_estimates,
    convex_hull_ci = convex_hull,
    lambda = sol_pool$lambda, nu = nu_pool,
    lambda_by_fold = lambda_by_fold, nu_by_fold = nu_by_fold,
    reweighted_by_fold = reweighted_by_fold,
    fold_estimate = fold_estimate, fold_se = fold_se, fold_ess = fold_ess,
    weights = weights_oof, fold_id = fold_id,
    common_covariates = common, adj_sets = adj_sets,
    ref_index = ref_index, alpha = alpha,
    settings = list(g_method = g_method, nu_use_h1 = nu_use_h1,
                    bc_g_from_cal = bc_g_from_cal, nu_regularize = nu_regularize,
                    bc_scale = bc_scale, propensity_clip = propensity_clip,
                    aipw_quantile = aipw_quantile, weight_quantile = weight_quantile,
                    psi_quantile = psi_quantile)
  )
  if (return_oof) {
    out$tau_hat_oof <- tau_hat_oof
    out$aipw_summands <- aipw_oof
    out$m_hat_oof <- m_hat_oof
    out$g_hat_oof <- g_hat_oof
    out$estimator_summand <- eta_oof
    out$influence <- psi_oof
  }
  out
}
