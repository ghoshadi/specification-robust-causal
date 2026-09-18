# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# This contains all simulations used in the main paper.  Run from the folder
# above simulations:
#   Rscript simulations/simulations.R
# The per-replication results are written to simulations/.
#
# The contrasts here come from interacted linear regression rather than the
# cross-fitted forests, so every fit is specrobust(..., reg_mode = "lm").
# ==============================================================================

rm(list=ls())
library(specrobust)
if (!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")

results_dir <- "simulations"

n_samples <- 1000; n_boot <- 100; n_sims <- 100
n_cores = parallel::detectCores()
alpha = 0.05

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

# Applying our reweighting approach
out = specrobust(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")),
                 reg_mode = "lm", alpha = alpha, n_boot = n_boot,
                 n_cores = n_cores, verbose = T)

# lm-based confidence intervals for each adjustment set
ci_1 = out$candidate_ci[1,]
ci_2 = out$candidate_ci[2,]

# Naive confidence interval
ci_naive = out$hull_ci

# specification-robust confidence intervals
rewt_ci = out$ci

# Printing the results
print(paste0("C.I. for adj set 1: [",round(ci_1, 3)[1],", ", round(ci_1, 3)[2],"]"))
print(paste0("C.I. for adj set 2: [",round(ci_2, 3)[1],", ", round(ci_2, 3)[2],"]"))
print(paste0("naive C.I. (range of prev C.I.'s): [",round(ci_naive, 3)[1],", ", 
             round(ci_naive, 3)[2],"]", " has length ",round(diff(ci_naive),3)))
print(paste0("specification-robust C.I.: [",round(rewt_ci, 3)[1],", ", 
             round(rewt_ci, 3)[2],"]", " has length ",round(diff(rewt_ci),3)))
print(out$lambda) # \lambda* 
mean(out$reweighted_estimates) # tauR = \E[w\tau_1]= \E[w\tau_2]
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

# Applying our reweighting approach
out = specrobust(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")),
                 reg_mode = "lm", alpha = alpha, n_boot = 0, verbose = T)

# lm-based confidence intervals for each adjustment set
ci_1 = out$candidate_ci[1,]
ci_2 = out$candidate_ci[2,]

# Naive confidence interval
ci_naive = out$hull_ci

# Plotting the histograms
plot(out, covariates = c("X1", "X2"), breaks = 200)

#-------------------------------------------------------
# Empirical coverage and average lenght of CIs
#-------------------------------------------------------

print(tauR <- mean(out$reweighted_estimates))

run_experiment <- function(itr) {
  set.seed(itr)
  n = n_samples
  tau = 1 # true ATE 
  X1 = rnorm(n); X2 = rnorm(n)
  A = as.numeric(runif(n) <= 1/(exp(5*X1+5*X2)+1))
  Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
  out <- specrobust(Y, A, data.frame(X1, X2), list(c("X1"), c("X1","X2")),
                    reg_mode = "lm", alpha = alpha, n_boot = n_boot,
                    n_cores = n_cores, seed = itr)

  rewt.ci.low = out$ci[1]
  rewt.ci.upp = out$ci[2]
  ireg.ci.low = out$candidate_ci_boot[, 1]
  ireg.ci.upp = out$candidate_ci_boot[, 2]
  
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
          file.path(results_dir, paste0("Eg1_",n_sims,"sims_",n_boot,"boot.csv")), row.names = F)

results_matrix = matrix(apply(read.csv(file.path(results_dir,
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

# Applying our reweighting approach
out = specrobust(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")),
                 reg_mode = "lm", alpha = alpha, n_boot = n_boot,
                 n_cores = n_cores, verbose = T)

# lm-based confidence intervals for each adjustment set
ci_1 = out$candidate_ci[1,]
ci_2 = out$candidate_ci[2,]

# Naive confidence interval
ci_naive = out$hull_ci

# specification-robust confidence intervals
rewt_ci = out$ci

# Printing the results
print(paste0("C.I. for adj set 1: [",round(ci_1, 3)[1],", ", round(ci_1, 3)[2],"]"))
print(paste0("C.I. for adj set 2: [",round(ci_2, 3)[1],", ", round(ci_2, 3)[2],"]"))
print(paste0("naive C.I. (range of prev C.I.'s): [",round(ci_naive, 3)[1],", ", 
             round(ci_naive, 3)[2],"]", " has length ",round(diff(ci_naive),3)))
print(paste0("specification-robust C.I.: [",round(rewt_ci, 3)[1],", ", 
             round(rewt_ci, 3)[2],"]", " has length ",round(diff(rewt_ci),3)))
print(out$lambda) # \lambda* 
mean(out$reweighted_estimates) # tauR = \E[w\tau_1]= \E[w\tau_2]
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

# Applying our reweighting approach
out = specrobust(Y, A, data.frame(X1, X2), list(c("X1"),c("X1","X2")),
                 reg_mode = "lm", alpha = alpha, n_boot = 0, verbose = T)

# lm-based confidence intervals for each adjustment set
ci_1 = out$candidate_ci[1,]
ci_2 = out$candidate_ci[2,]

# Naive confidence interval
ci_naive = out$hull_ci

# Plotting the histograms
plot(out, covariates = c("X1", "X2"), breaks = 200)

#-------------------------------------------------------
# Empirical coverage and average lenght of CIs
#-------------------------------------------------------

print(tauR <- mean(out$reweighted_estimates))

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
  out <- specrobust(Y, A, data.frame(X1, X2), list(c("X1"), c("X1","X2")),
                    reg_mode = "lm", alpha = alpha, n_boot = n_boot,
                    n_cores = n_cores, seed = itr)

  rewt.ci.low = out$ci[1]
  rewt.ci.upp = out$ci[2]
  ireg.ci.low = out$candidate_ci_boot[, 1]
  ireg.ci.upp = out$candidate_ci_boot[, 2]
  
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
          file.path(results_dir, paste0("Eg2_",n_sims,"sims_",n_boot,"boot.csv")), row.names = F)

results_matrix = matrix(apply(read.csv(file.path(results_dir,
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

results_matrix = matrix(apply(read.csv(file.path(results_dir,
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

results_matrix = matrix(apply(read.csv(file.path(results_dir,
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
