library(ranger)
library(grf)
library(dplyr)
library(plotrix)
library(parallel)
library(pbapply)
library(pbmcapply)
library(foreign)
library(MASS)


## Function for finding weights

get_weights <- function(g, scale_factor = 1) {
  g <- as.matrix(g)
  if(is.null(scale_factor)) scale_factor = diff(range(g))
  g <- g / scale_factor
  objective <- function(lmbda) (mean(exp(g %*% lmbda)))
  result <- optim(par = rep(0, ncol(g)), fn = objective, 
                  method = "BFGS", control = list(reltol = 1e-12, maxit = 1e4))
  opt_lmbda <- result$par
  if (!anyNA(opt_lmbda)) {
    opt_w <- exp(g %*% opt_lmbda)
    opt_w <- opt_w / mean(opt_w)
  } else {
    opt_w <- rep(1, nrow(g))
    opt_lmbda <- NA
    print("no solution found, defaulting to equal weights")
  }
  list(weights = as.numeric(opt_w), lambda = opt_lmbda/scale_factor)
}

library(nleqslv)

get_weights_newton <- function(g, tol = 1e-14, maxit = 2000) {
  g <- as.matrix(g)
  n <- nrow(g); K <- ncol(g)
  # F(lambda):  1/n * t(g) %*% exp(g %*% lambda)
  F <- function(lambda) {
    w_vec <- as.numeric(exp(g %*% lambda))
    crossprod(g, w_vec) / n
  }
  # J(lambda):  1/n * t(g * w) %*% g
  J <- function(lambda) {
    w_vec <- as.numeric(exp(g %*% lambda))
    crossprod(g * w_vec, g) / n
  }
  sol <- nleqslv(x = rep(0, K), fn = F, jac = J,
                 method = "Newton",
                 control = list(ftol = tol,xtol = tol, maxit = maxit))
  
  lambda_hat <- sol$x
  w_raw      <- as.numeric(exp(g %*% lambda_hat))
  w_norm     <- w_raw / mean(w_raw)
  list(weights = w_norm, lambda = lambda_hat)
}


## Assumption-robust inference on ATE using linear models

assump_robust_lm <- function(response, treatment, covariates, adj_sets, verbose = FALSE) {
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

## Assumption-robust inference on ATE using flexible nonparametric regressors
assump_robust <- function(response, treatment, covariates, adj_sets, verbose = FALSE) {
  
  # Step 1: Estimate the contrasts E[Y | A = 1, X_{S_k}] - E[Y | A = 0, X_{S_k}] for each adjustment set S_k
  k <- length(adj_sets); n <- length(response)
  cfs <- lapply(adj_sets, function(adj){
    causal_forest(covariates[, adj, drop = FALSE], response, treatment)
  })
  cates <- lapply(cfs, function(cf) predict(cf)$predictions)
  
  # Step 2: Regress the contrasts on the common covariates 
  common_covariates <- covariates[, Reduce(intersect, adj_sets), drop = FALSE]
  # # g records cates regressed on common covariates
  g <- lapply(seq_len(k), function(i) {
    fit <- lm(cates[[i]] ~ ., data = common_covariates)
    predict(fit)
  })
  
  # Step 3: Obtain weights that empirically solve the KL optimization problem
  Delta_g <- do.call(cbind, lapply(seq_len(k)[-1], function(i) g[[1]] - g[[i]]))
  out <- get_weights_newton(Delta_g)
  if(verbose){
    cat("lambda_hat: "); print(round(out$lambda,4))
  }
  
  # Step 4: Compute AIPW summands with standard errors
  prop_scores = list()
  aipw_computation <- sapply(seq_len(k), function(i) {
    ps_model <- glm(treatment ~ ., data = covariates[, adj_sets[[i]], drop = FALSE], 
                    family = binomial)
    ps <- predict(ps_model, type = "response")
    prop_scores[[i]] = ps
    rf_mu0 <- ranger(response ~ ., 
                     data = na.omit(data.frame(response = response[treatment == 0],
                                               covariates[treatment == 0, adj_sets[[i]], drop = FALSE])
                     )
    )
    
    mu0 <- predict(rf_mu0, data = covariates[, adj_sets[[i]], drop = FALSE])$predictions
    mu1 <- mu0 + cates[[i]]
    
    aipw_cate <- ((treatment * (response - mu1)) / ps) - ((1 - treatment) * (response - mu0) / (1 - ps)) + mu1 - mu0
    aipw_ests <- mean(aipw_cate)
    rewt_aipw <- mean(out$weights * aipw_cate)
    aipw_se <- sd(aipw_cate)/sqrt(n)
    rewt_aipw_se <- sd(out$weights * aipw_cate)/sqrt(n)
    return(list(aipw_ests = aipw_ests, 
                rewt_aipw = rewt_aipw, 
                prop_scores = prop_scores,
                aipw_se = aipw_se,
                rewt_aipw_se = rewt_aipw_se,
                aipw_cate = aipw_cate))
  })
  
  # Step 5: Aggregate estimators and compute assumption-robust estimator without bias-correction
  est_grf <- sapply(g, function(g_i) mean(g_i))
  rewt_est_grf <- sapply(cates, function(tau_i) mean(out$weights * tau_i))
  aipw_ests <- unlist(aipw_computation["aipw_ests",])
  rewt_aipw <- unlist(aipw_computation["rewt_aipw",])
  prop_scores <- unlist(aipw_computation["prop_scores",])
  boxplot(cbind(est_grf, rewt_est_grf, aipw_ests, rewt_aipw))
  assump_robust_est <- mean(rewt_est_grf) # equal weights, can serve as tau_init
  
  # Step 6: Finding the nu's for the optimal convex combination
  nu_2toK = lm(cates[[2]] ~ Delta_g -1, weights = out$weights)$coef
  nu = c(1-sum(nu_2toK), nu_2toK)
  if(verbose){
    cat("nu_hat: "); print(round(nu,4))
  }
  
  # Step 7: Bias correction step
  aipw_cates = do.call(cbind, lapply(seq_len(k), function(i) unlist(aipw_computation['aipw_cate',][i])))
  Delta_cates <- do.call(cbind, lapply(seq_len(k)[-1], function(i) cates[[1]] - cates[[i]]))
  bc <- out$weights * (Reduce(`+`, Map(`*`, nu, g)) - assump_robust_est) * (Delta_cates - Delta_g) %*% out$lambda
  
  assump_robust_est_bc <- weighted.mean(rewt_aipw, nu) + mean(bc) # check again about this sign
  
  # Step 8: estimating the standard error of the bias-corrected estimator
  # write code for std err here and return all the std errs
  new_se = sd(out$weights * aipw_cates %*% nu + bc)/sqrt(n)
  
  return(list(aipw_estimates = unlist(aipw_ests),
              rewt.estimates = rewt_est_grf,
              rewt_aipw_estimates = unlist(rewt_aipw), 
              assump_robust_est = assump_robust_est,
              assump_robust_est_bc = assump_robust_est_bc,
              aipw_se = unlist(aipw_computation['aipw_se',]),
              rewt_aipw_se = new_se,
              prop_scores = prop_scores,
              nu = nu, cates = g,
              grf_estimates = est_grf,
              weights = out$weights, 
              lambda = out$lambda)
  )
}

## Auxiliary functions
get_ci_from_ests <- function(tau.hat, tau.se, alpha = 0.05){
  return(c(tau.hat - qnorm(alpha/2, lower.tail = F)*tau.se, 
           tau.hat + qnorm(alpha/2, lower.tail = F)*tau.se))
}


## Function for comparing histograms side by side

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

clip <- function(p) pmin(pmax(p, 0.01), 0.99)

cross_fit_assump_robust <- function(response,
                                    treatment,
                                    covariates,
                                    adj_sets,
                                    verbose = FALSE,
                                    seed = 42) {
  set.seed(seed)
  n <- length(response)
  # random two‐fold split
  fold_id <- sample(rep(1:2, length.out = n))
  fold_est = fold_se = rep(0, 2)
  weights <- rep(0, n)
  aipw_summands <- vector("list", length(adj_sets))
  for(j in seq_along(adj_sets)){
    aipw_summands[[j]] <- rep(0, n)
  }
  for (fold in 1:2) {
    # define train/test indices
    tr_idx <- which(fold_id != fold)
    te_idx <- which(fold_id == fold)
    
    # subset data
    y_tr <- response[tr_idx];  y_te <- response[te_idx]
    A_tr <- treatment[tr_idx]; A_te <- treatment[te_idx]
    X_tr <- covariates[tr_idx, , drop = FALSE]
    X_te <- covariates[te_idx, , drop = FALSE]
    
    ## STEP 1: compute AIPW estimate on test fold (using models trained on training fold)
    
    # 1) Propensity scores
    ps_models <- lapply(adj_sets, function(adj)
      glm(A_tr ~ ., data = data.frame(X_tr[, adj, drop = FALSE]), family = binomial)
    )
    ps_te <- lapply(seq_along(ps_models), function(i)
      predict(ps_models[[i]], data.frame(X_te[, adj_sets[[i]], drop = FALSE]), type = "response")
    )
    ps_tr <- lapply(seq_along(ps_models), function(i)
      predict(ps_models[[i]], X_tr[, adj_sets[[i]], drop = FALSE], type = "response")
    )
    ps_te = lapply(ps_te, clip)
    ps_tr = lapply(ps_tr, clip)
    
    # 2) Outcome models
    mu0_mods <- lapply(adj_sets, function(adj) {
      df0 <- data.frame(y = y_tr[A_tr==0], X_tr[A_tr==0, adj, drop = FALSE])
      ranger(y ~ ., data = df0)
    })
    mu0_tr  <- mapply(function(mod, adj)
      predict(mod, X_tr[, adj, drop = FALSE])$predictions,
      mu0_mods, adj_sets, SIMPLIFY = FALSE
    )
    mu0_te  <- mapply(function(mod, adj)
      predict(mod, X_te[, adj, drop = FALSE])$predictions,
      mu0_mods, adj_sets, SIMPLIFY = FALSE
    )
    
    mu1_mods <- lapply(adj_sets, function(adj) {
      df1 <- data.frame(y = y_tr[A_tr==1], X_tr[A_tr==1, adj, drop = FALSE])
      ranger(y ~ ., data = df1)
    })
    mu1_tr  <- mapply(function(mod, adj)
      predict(mod, X_tr[, adj, drop = FALSE])$predictions,
      mu1_mods, adj_sets, SIMPLIFY = FALSE
    )
    mu1_te  <- mapply(function(mod, adj)
      predict(mod, X_te[, adj, drop = FALSE])$predictions,
      mu1_mods, adj_sets, SIMPLIFY = FALSE
    )
    
    # 3) AIPW summands
    tau_tr <- mapply(function(ps, y, mu0, mu1) {
      ((A_tr * (y - mu1)) / ps) -
        (((1-A_tr) * (y - mu0)) / (1-ps)) + (mu1 - mu0)
    }, ps_tr, MoreArgs = list(y = y_tr), mu0 = mu0_tr, mu1 = mu1_tr,
    SIMPLIFY = FALSE)
    
    tau_te <- mapply(function(ps, y, mu0, mu1) {
      ((A_te * (y - mu1)) / ps) -
        (((1-A_te) * (y - mu0)) / (1-ps)) + (mu1 - mu0)
    }, ps_te, MoreArgs = list(y = y_te), mu0 = mu0_te, mu1 = mu1_te,
    SIMPLIFY = FALSE)
    for(j in seq_along(adj_sets)) {
      aipw_summands[[j]][te_idx] <- tau_te[[j]]
    }
    
    ## STEP 2: regress contrasts on common covariates (again on the training fold)
    common <- Reduce(intersect, adj_sets)
    
    G_tr <- lapply(tau_tr, function(tau)
      predict(lm(tau ~ ., data = X_tr[, common, drop = FALSE]))
    )
    
    G_te <- lapply(tau_tr, function(tau)
      predict(lm(tau ~ ., data = X_tr[, common, drop = FALSE]),
              newdata = X_te[, common, drop = FALSE])
    )
    
    ## STEP 3: solve for weights on training fold
    Delta_G_tr <- do.call(cbind,
                          lapply(2:length(G_tr), function(i) G_tr[[1]] - G_tr[[i]])
    )
    w_out <- get_weights(Delta_G_tr)
    w_tr  <- w_out$weights
    lambda <- w_out$lambda
    Delta_G_te <- do.call(cbind,
                          lapply(2:length(G_te), function(i) G_te[[1]] - G_te[[i]])
    )
    w_te <- exp(Delta_G_te %*% lambda)
    w_te <- w_te / mean(w_te)
    if (verbose) cat("Fold", fold, "lambda:", round(lambda,4), "\n")
    
    unlist(lapply(G_te, function(z) weighted.mean(z, w_te)))
    unlist(lapply(tau_te, function(z) weighted.mean(z, w_te)))
    unlist(lapply(G_tr, function(z) weighted.mean(z, w_tr)))
    unlist(lapply(tau_tr, function(z) weighted.mean(z, w_tr)))
    # bias‐correction on test
    nu2 <- lm(I(G_te[[1]]-mean(G_te[[1]] * w_te)) ~ Delta_G_te - 1, weights = w_te)$coef
    #nu2 <- lm(I(G_tr[[1]]-mean(G_tr[[1]] * w_tr)) ~ Delta_G_tr - 1, weights = w_tr)$coef
    nu  <- c(1 - sum(nu2), nu2)
    nu[is.na(nu)] = 0
    bkt_1 <- (sapply(2:length(adj_sets), function(i) tau_te[[1]]-tau_te[[i]] - Delta_G_te[,i-1]) * as.vector(w_te)) %*% lambda
    bkt_2 <- scale(apply(sapply(seq_along(G_te), function(i) nu[i] * G_te[[i]]), 1, sum), center = T, scale = F)
    bc_term <- bkt_1 * bkt_2
    
    # final fold estimate & SE
    est_bc <- sum(unlist(lapply(tau_te, function(z) weighted.mean(z, w_te))) * nu) + mean(bc_term)
    se_bc  <- sd((do.call(cbind, tau_te)%*% nu)* w_te  + bc_term) / sqrt(length(te_idx))
    
    fold_est[fold] <- est_bc
    fold_se[fold]  <- se_bc
    weights[te_idx] <- w_te
  }
  
  list(
    weights = weights,
    assump_robust_est = mean(fold_est),
    assump_robust_se = sqrt(mean(fold_se^2)),
    aipw_estimates = lapply(seq_along(adj_sets), function(j) { mean(aipw_summands[[j]]) }),
    aipw_se = lapply(seq_along(adj_sets), function(j) { sd(aipw_summands[[j]]) / sqrt(n)})
  )
}
