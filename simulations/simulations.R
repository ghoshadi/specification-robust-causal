# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# This contains all simulations used in the main paper.  Run from the folder
# above simulations:
#   Rscript simulations/simulations.R
# The per-replication results are written to simulations/.
#
# The contrasts here come from interacted linear regression, not from the
# cross-fitted forests of main.R; the only piece of the shared code this needs
# is solve_lam() in utils.R, so specification_robust_lm() is defined below.
# ==============================================================================

rm(list=ls())
source('./utils.R')

RESULTS <- "simulations"

n_samples <- 1000; n_boot <- 100; n_sims <- 100
n_cores = parallel::detectCores()
alpha = 0.05; z = qnorm(alpha/2, lower.tail = F)

# ===========================================================================
# Specification-robust inference on the ATE using linear models
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
  out <- solve_lam(Delta_g)
  if (verbose) print(paste0("lambda_hat: ",round(out$lambda, 3)))
  
  rewt_est <- sapply(g, function(gi) mean(out$weights * gi))
  theta_hats <- sapply(fits, function(fit) coef(fit)["treatment"])
  
  list(rewt.estimates = rewt_est, 
       ireg.estimates = theta_hats, 
       weights = out$weights, 
       lambda = out$lambda)
}

##-------------------------------------------------------
## Example: Two confounders
##-------------------------------------------------------

set.seed(123)
n = n_samples
tau = 1 # true ATE
X1 = rnorm(n); X2 = rnorm(n)
A = as.numeric(runif(n) <= 1/(exp(5*X1+5*X2)+1))
Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
print(summary(fit1 <- lm(Y ~ A*X1))$coef["A",], digits=3) # incorrect 
print(summary(fit2 <- lm(Y ~ A*(X1+X2)))$coef["A",], digits=3) # correct

# lm-based confidence intervals for each adjustment set
ci_1 = get_ci_from_ests(summary(fit1)$coef["A",1], summary(fit1)$coef["A",2], alpha)
ci_2 = get_ci_from_ests(summary(fit2)$coef["A",1], summary(fit2)$coef["A",2], alpha)

# Naive confidence interval
ci_naive = range(c(ci_1,ci_2))

# Applying our reweighting approach
out = specification_robust_lm(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")), verbose = T)

# bootstrap estimate of the variance of reweighted estimators
bootstrap_specification_robust_lm <- function(iter) {
  indices <- sample(1:n, replace = TRUE)
  Y_boot <- Y[indices]
  A_boot <- A[indices]
  X_boot <- cbind(X1[indices], X2[indices])
  colnames(X_boot) = c("X1", "X2")
  result <- specification_robust_lm(Y_boot, A_boot, X_boot, list(c("X1"), c("X1","X2")))
  return(list(rewt = result$rewt.estimates, ireg = result$ireg.estimates))
}
boot_results <- pbmcapply::pbmclapply(1:n_boot, bootstrap_specification_robust_lm, mc.cores = n_cores)
rewt.estimates <- do.call(rbind, lapply(boot_results, function(x) x$rewt))
rewt_se <- mean(apply(rewt.estimates, 2, sd)) 

# specification-robust confidence intervals
rewt_ci = get_ci_from_ests(mean(out$rewt.estimates), rewt_se)

# Printing the results
print(paste0("C.I. for adj set 1: [",round(ci_1, 3)[1],", ", round(ci_1, 3)[2],"]"))
print(paste0("C.I. for adj set 2: [",round(ci_2, 3)[1],", ", round(ci_2, 3)[2],"]"))
print(paste0("naive C.I. (range of prev C.I.'s): [",round(ci_naive, 3)[1],", ", 
             round(ci_naive, 3)[2],"]", " has length ",round(diff(ci_naive),3)))
print(paste0("specification-robust C.I.: [",round(rewt_ci, 3)[1],", ", 
             round(rewt_ci, 3)[2],"]", " has length ",round(diff(rewt_ci),3)))
print(out$lambda) # \lambda* 
mean(out$rewt.estimates) # tauR = \E[w\tau_1]= \E[w\tau_2]
mean(out$weights*X1) # \E_w[X_1]

#-------------------------------------------------------
# Preparing a plot for the paper
#-------------------------------------------------------

set.seed(123)
n <- 1e6
tau = 1 # true ATE
X1 = rnorm(n); X2 = rnorm(n)
A = as.numeric(runif(n) <= 1/(exp(5*X1+5*X2)+1))
Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
print(summary(fit1 <- lm(Y ~ A*X1))$coef["A",], digits=3) # incorrect 
print(summary(fit2 <- lm(Y ~ A*(X1+X2)))$coef["A",], digits=3) # correct

# lm-based confidence intervals for each adjustment set
ci_1 = get_ci_from_ests(summary(fit1)$coef["A",1], summary(fit1)$coef["A",2], alpha)
ci_2 = get_ci_from_ests(summary(fit2)$coef["A",1], summary(fit2)$coef["A",2], alpha)

# Naive confidence interval
ci_naive = range(c(ci_1,ci_2))

# Applying our reweighting approach
out = specification_robust_lm(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")), verbose = T)

# Plotting the histograms
plot_hists(data.frame(X1, X2), out$weights, vars = c("X1", "X2"), breaks = 200)

#-------------------------------------------------------
# Empirical coverage and average lenght of CIs
#-------------------------------------------------------

print(tauR <- mean(out$rewt.estimates))

run_experiment <- function(itr) {
  set.seed(itr)
  n = n_samples
  tau = 1 # true ATE 
  X1 = rnorm(n); X2 = rnorm(n)
  A = as.numeric(runif(n) <= 1/(exp(5*X1+5*X2)+1))
  Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
  out <- specification_robust_lm(Y, A, cbind(X1, X2), list(c("X1"), c("X1","X2")))
  bootstrap_specification_robust_lm <- function(iter) {
    indices <- sample(1:length(Y), replace = TRUE)
    Y_boot <- Y[indices]
    A_boot <- A[indices]
    X_boot <- cbind(X1[indices], X2[indices])
    colnames(X_boot) = c("X1", "X2")
    result <- specification_robust_lm(Y_boot, A_boot, X_boot, list(c("X1"), c("X1", "X2")))
    return(list(rewt = result$rewt.estimates, ireg = result$ireg.estimates))
  }
  boot_results <- parallel::mclapply(1:n_boot, bootstrap_specification_robust_lm, mc.cores = n_cores)
  rewt.estimates <- do.call(rbind, lapply(boot_results, function(x) x$rewt))
  ireg_estimates <- do.call(rbind, lapply(boot_results, function(x) x$ireg))
  
  rewt_variance <- apply(rewt.estimates, 2, var)
  ireg_variance <- apply(ireg_estimates, 2, var)
  
  rewt.ci.low = mean(out$rewt.estimates - z*sqrt(rewt_variance))
  rewt.ci.upp = mean(out$rewt.estimates + z*sqrt(rewt_variance))
  ireg.ci.low = out$ireg.estimates - z*sqrt(ireg_variance)
  ireg.ci.upp = out$ireg.estimates + z*sqrt(ireg_variance)
  
  naive_CI = range(c(ireg.ci.low, ireg.ci.upp))
  
  ireg.ecov = as.numeric(ireg.ci.low <= tau)*as.numeric(tau <= ireg.ci.upp)
  ireg.len = as.numeric(ireg.ci.upp-ireg.ci.low)
  naive.ecov = as.numeric(naive_CI[1] <= tau)*as.numeric(tau <= naive_CI[2])
  naive.len = naive_CI[2] - naive_CI[1]
  
  rewt.ecov = as.numeric(rewt.ci.low <= tauR)*as.numeric(tauR <= rewt.ci.upp)
  rewt.len = as.numeric(rewt.ci.upp-rewt.ci.low)
  
  res <- c(ireg.ecov, naive.ecov, rewt.ecov, ireg.len, naive.len, rewt.len)
  names(res) <- c(paste0("ireg", seq_along(ireg.ecov), ".ecov"), "naive.ecov", "rewt.ecov",
                  paste0("ireg", seq_along(ireg.len), ".len"), "naive.len", "rewt.len")
  return(res)
}

set.seed(123)
results <- pbapply::pblapply(1:n_sims, run_experiment)
write.csv(do.call(rbind, results),
          file.path(RESULTS, paste0("Eg1_",n_sims,"sims_",n_boot,"boot.csv")), row.names = F)

results_matrix = matrix(apply(read.csv(file.path(RESULTS,
                        paste0("Eg1_",n_sims,"sims_",n_boot,"boot.csv"))), 2, mean),
                        nrow = 4, byrow = F,
                        dimnames = list(c("C.I. using adj set 1", 
                                          "C.I. using adj set 2", 
                                          "naive C.I. (range of prev C.I.'s)", 
                                          "specification-robust C.I."),
                                        c("emp coverage", "avg width")
                        )
)
print(results_matrix, digits = 3)

##-------------------------------------------------------
## Example: M-bias
##-------------------------------------------------------

set.seed(123)
n = n_samples
tau = 1 # true ATE 
U1 = rnorm(n); U2 = rnorm(n) # unobserved
X1 = rnorm(n)
A = as.numeric(U1 + X1 > 0)
X2 = U1 + U2 # adjusting for X2 introduces M-bias
Y1 = tau + X1 - U2 + rnorm(n); Y0 = 5*U2 + rnorm(n)
Y = ifelse(A==1, Y1, Y0) 
print(summary(fit1 <- lm(Y ~ A*X1))$coef["A",], digits=3) # correct 
print(summary(fit2 <- lm(Y ~ A*(X1+X2)))$coef["A",], digits=3) # incorrect

# lm-based confidence intervals for each adjustment set
ci_1 = get_ci_from_ests(summary(fit1)$coef["A",1], summary(fit1)$coef["A",2], alpha)
ci_2 = get_ci_from_ests(summary(fit2)$coef["A",1], summary(fit2)$coef["A",2], alpha)

# Naive confidence interval
ci_naive = range(c(ci_1,ci_2))

# Applying our reweighting approach
out = specification_robust_lm(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")), verbose = T)

# bootstrap estimate of the variance of reweighted estimators
bootstrap_specification_robust_lm <- function(iter) {
  indices <- sample(1:n, replace = TRUE)
  Y_boot <- Y[indices]
  A_boot <- A[indices]
  X_boot <- cbind(X1[indices], X2[indices])
  colnames(X_boot) = c("X1", "X2")
  result <- specification_robust_lm(Y_boot, A_boot, X_boot, list(c("X1"), c("X1","X2")))
  return(list(rewt = result$rewt.estimates, ireg = result$ireg.estimates))
}
boot_results <- pbmcapply::pbmclapply(1:n_boot, bootstrap_specification_robust_lm, mc.cores = n_cores)
rewt.estimates <- do.call(rbind, lapply(boot_results, function(x) x$rewt))
rewt_se <- mean(apply(rewt.estimates, 2, sd)) 

# specification-robust confidence intervals
rewt_ci = get_ci_from_ests(mean(out$rewt.estimates), rewt_se)

# Printing the results
print(paste0("C.I. for adj set 1: [",round(ci_1, 3)[1],", ", round(ci_1, 3)[2],"]"))
print(paste0("C.I. for adj set 2: [",round(ci_2, 3)[1],", ", round(ci_2, 3)[2],"]"))
print(paste0("naive C.I. (range of prev C.I.'s): [",round(ci_naive, 3)[1],", ", 
             round(ci_naive, 3)[2],"]", " has length ",round(diff(ci_naive),3)))
print(paste0("specification-robust C.I.: [",round(rewt_ci, 3)[1],", ", 
             round(rewt_ci, 3)[2],"]", " has length ",round(diff(rewt_ci),3)))
print(out$lambda) # \lambda* 
mean(out$rewt.estimates) # tauR = \E[w\tau_1]= \E[w\tau_2]
mean(out$weights*X1) # \E_w[X_1]

#-------------------------------------------------------
# Preparing a plot for the paper
#-------------------------------------------------------

set.seed(123)
n <- 1e6
tau = 1 # true ATE 
U1 = rnorm(n); U2 = rnorm(n) # unobserved
X1 = rnorm(n)
A = as.numeric(U1 + X1 > 0)
X2 = U1 + U2 # adjusting for X2 introduces M-bias
Y1 = tau + X1 - U2 + rnorm(n); Y0 = 5*U2 + rnorm(n)
Y = ifelse(A==1, Y1, Y0) 
print(summary(fit1 <- lm(Y ~ A*X1))$coef["A",], digits=3) # incorrect 
print(summary(fit2 <- lm(Y ~ A*(X1+X2)))$coef["A",], digits=3) # correct

# lm-based confidence intervals for each adjustment set
ci_1 = get_ci_from_ests(summary(fit1)$coef["A",1], summary(fit1)$coef["A",2], alpha)
ci_2 = get_ci_from_ests(summary(fit2)$coef["A",1], summary(fit2)$coef["A",2], alpha)

# Naive confidence interval
ci_naive = range(c(ci_1,ci_2))

# Applying our reweighting approach
out = specification_robust_lm(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")), verbose = T)

# Plotting the histograms
plot_hists(data.frame(X1, X2), out$weights, vars = c("X1", "X2"), breaks = 200)

#-------------------------------------------------------
# Empirical coverage and average lenght of CIs
#-------------------------------------------------------

print(tauR <- mean(out$rewt.estimates))

run_experiment <- function(itr) {
  set.seed(itr)
  n = n_samples
  tau = 1 # true ATE 
  U1 = rnorm(n); U2 = rnorm(n) # unobserved
  X1 = rnorm(n)
  A = as.numeric(U1 + X1 > 0)
  X2 = U1 + U2 # adjusting for X2 introduces M-bias
  Y1 = tau + X1 - U2 + rnorm(n); Y0 = 5*U2 + rnorm(n)
  Y = ifelse(A==1, Y1, Y0) 
  out <- specification_robust_lm(Y, A, cbind(X1, X2), list(c("X1"), c("X1","X2")))
  bootstrap_specification_robust_lm <- function(iter) {
    indices <- sample(1:length(Y), replace = TRUE)
    Y_boot <- Y[indices]
    A_boot <- A[indices]
    X_boot <- cbind(X1[indices], X2[indices])
    colnames(X_boot) = c("X1", "X2")
    result <- specification_robust_lm(Y_boot, A_boot, X_boot, list(c("X1"), c("X1", "X2")))
    return(list(rewt = result$rewt.estimates, ireg = result$ireg.estimates))
  }
  boot_results <- parallel::mclapply(1:n_boot, bootstrap_specification_robust_lm, mc.cores = n_cores)
  rewt.estimates <- do.call(rbind, lapply(boot_results, function(x) x$rewt))
  ireg_estimates <- do.call(rbind, lapply(boot_results, function(x) x$ireg))
  
  rewt_variance <- apply(rewt.estimates, 2, var)
  ireg_variance <- apply(ireg_estimates, 2, var)
  
  rewt.ci.low = mean(out$rewt.estimates - z*sqrt(rewt_variance))
  rewt.ci.upp = mean(out$rewt.estimates + z*sqrt(rewt_variance))
  ireg.ci.low = out$ireg.estimates - z*sqrt(ireg_variance)
  ireg.ci.upp = out$ireg.estimates + z*sqrt(ireg_variance)
  
  naive_CI = range(c(ireg.ci.low, ireg.ci.upp))
  
  ireg.ecov = as.numeric(ireg.ci.low <= tau)*as.numeric(tau <= ireg.ci.upp)
  ireg.len = as.numeric(ireg.ci.upp-ireg.ci.low)
  naive.ecov = as.numeric(naive_CI[1] <= tau)*as.numeric(tau <= naive_CI[2])
  naive.len = naive_CI[2] - naive_CI[1]
  
  rewt.ecov = as.numeric(rewt.ci.low <= tauR)*as.numeric(tauR <= rewt.ci.upp)
  rewt.len = as.numeric(rewt.ci.upp-rewt.ci.low)
  
  res <- c(ireg.ecov, naive.ecov, rewt.ecov, ireg.len, naive.len, rewt.len)
  names(res) <- c(paste0("ireg", seq_along(ireg.ecov), ".ecov"), "naive.ecov", "rewt.ecov",
                  paste0("ireg", seq_along(ireg.len), ".len"), "naive.len", "rewt.len")
  return(res)
}

set.seed(123)
results <- pbapply::pblapply(1:n_sims, run_experiment)
write.csv(do.call(rbind, results),
          file.path(RESULTS, paste0("Eg2_",n_sims,"sims_",n_boot,"boot.csv")), row.names = F)

results_matrix = matrix(apply(read.csv(file.path(RESULTS,
                        paste0("Eg2_",n_sims,"sims_",n_boot,"boot.csv"))), 2, mean),
                        nrow = 4, byrow = F,
                        dimnames = list(c("C.I. using adj set 1", 
                                          "C.I. using adj set 2", 
                                          "naive C.I. (range of prev C.I.'s)", 
                                          "specification-robust C.I."),
                                        c("emp coverage", "avg width")
                        )
)
print(results_matrix, digits = 3)


#-------------------------------------------------------
# n_sims = 1000; n_boot = 1000
#-------------------------------------------------------

results_matrix = matrix(apply(read.csv(file.path(RESULTS,
                        "Eg1_1000sims_1000boot.csv")), 2, mean),
                        nrow = 4, byrow = F,
                        dimnames = list(c("C.I. using adj set 1", 
                                          "C.I. using adj set 2", 
                                          "naive C.I. (range of prev C.I.'s)", 
                                          "specification-robust C.I."),
                                        c("emp coverage", "avg width")
                        )
)
print(results_matrix, digits = 3)

results_matrix = matrix(apply(read.csv(file.path(RESULTS,
                        "Eg2_1000sims_1000boot.csv")), 2, mean),
                        nrow = 4, byrow = F,
                        dimnames = list(c("C.I. using adj set 1", 
                                          "C.I. using adj set 2", 
                                          "naive C.I. (range of prev C.I.'s)", 
                                          "specification-robust C.I."),
                                        c("emp coverage", "avg width")
                        )
)
print(results_matrix, digits = 3)
