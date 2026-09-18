lm_one_fit <- function(response, treatment, covariates, adj_sets) {
  K <- length(adj_sets)
  n <- nrow(covariates)

  fits <- lapply(adj_sets, function(adj) {
    form <- stats::as.formula(paste("response ~ treatment * (",
                                    paste(adj, collapse = " + "), ")"))
    stats::lm(form, data = data.frame(response = response, treatment = treatment,
                                      covariates[, adj, drop = FALSE]))
  })

  int_adj <- Reduce(intersect, adj_sets)
  pred0 <- data.frame(treatment = rep(0, n), covariates)
  pred1 <- data.frame(treatment = rep(1, n), covariates)
  common_covariates <- as.data.frame(covariates[, int_adj, drop = FALSE])

  g <- lapply(seq_len(K), function(k) {
    cate <- stats::predict(fits[[k]], newdata = pred1) -
            stats::predict(fits[[k]], newdata = pred0)
    stats::lm(cate ~ ., data = common_covariates)$fitted.values
  })

  Delta_g <- do.call(cbind, lapply(seq_len(K)[-1], function(k) g[[1]] - g[[k]]))
  tilt <- solve_lam(Delta_g)

  cf <- lapply(fits, function(f) summary(f)$coefficients)
  list(reweighted = vapply(g, function(gk) mean(tilt$weights * gk), numeric(1)),
       candidate = vapply(cf, function(m) m["treatment", 1], numeric(1)),
       candidate_se = vapply(cf, function(m) m["treatment", 2], numeric(1)),
       weights = tilt$weights, lambda = tilt$lambda)
}

fit_lm <- function(response, treatment, covariates, adj_sets, ref_index,
                   verbose, alpha, seed, n_boot, n_cores) {

  K <- length(adj_sets)
  if (is.null(ref_index)) ref_index <- 1
  if (is.character(ref_index))
    ref_index <- switch(ref_index,
      first    = 1,
      largest  = which.max(lengths(adj_sets)),
      smallest = which.min(lengths(adj_sets)),
      stop("ref_index must be a numeric index, or one of ",
           "\"first\", \"largest\", \"smallest\""))
  ref_index <- as.integer(ref_index)
  stopifnot(length(ref_index) == 1, ref_index >= 1, ref_index <= K)
  if (ref_index != 1) {
    perm <- c(ref_index:K, if (ref_index > 1) 1:(ref_index - 1))
    adj_sets <- adj_sets[perm]
    if (verbose) cat("Reference: original S", ref_index, " moved to index 1.\n", sep = "")
  }

  y <- as.numeric(response); a <- as.integer(treatment)
  X <- as.data.frame(covariates)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  n <- length(y); common <- Reduce(intersect, adj_sets)

  if (verbose)
    cat(sprintf("K=%d  common={%s}  n=%d\n  %d bootstrap replications\n\n",
                K, paste(common, collapse = ", "), n, n_boot))

  point <- lm_one_fit(y, a, X, adj_sets)

  estimate <- mean(point$reweighted)
  candidate_estimates <- point$candidate
  candidate_se <- point$candidate_se
  candidate_ci <- t(vapply(seq_len(K), function(k)
    get_ci_from_ests(candidate_estimates[k], candidate_se[k], alpha), numeric(2)))
  hull_ci <- c(min(candidate_ci[, 1]), max(candidate_ci[, 2]))

  if (n_boot == 0) {
    n_rep <- 0
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
    candidate_se_boot <- rep(NA_real_, K)
    candidate_ci_boot <- matrix(NA_real_, K, 2)
  } else {
    idx <- withr::with_seed(seed,
      lapply(seq_len(n_boot), function(b) sample.int(n, n, replace = TRUE)))
    one_boot <- function(i)
      tryCatch(lm_one_fit(y[i], a[i], X[i, , drop = FALSE], adj_sets),
               error = function(e) NULL)
    reps <- if (n_cores > 1 && .Platform$OS.type == "unix")
              parallel::mclapply(idx, one_boot, mc.cores = n_cores)
            else lapply(idx, one_boot)

    ok <- vapply(reps, function(r)
      is.list(r) && all(is.finite(r$reweighted)) && all(is.finite(r$candidate)),
      logical(1))
    if (!any(ok)) stop("Every bootstrap replication failed.")
    if (!all(ok))
      warning(sum(!ok), " of ", n_boot,
              " bootstrap replications were dropped as non-finite.", call. = FALSE)
    reps <- reps[ok]
    n_rep <- length(reps)

    boot_reweighted <- do.call(rbind, lapply(reps, `[[`, "reweighted"))
    boot_candidate  <- do.call(rbind, lapply(reps, `[[`, "candidate"))

    se <- mean(apply(boot_reweighted, 2, stats::sd))
    ci <- get_ci_from_ests(estimate, se, alpha)
    candidate_se_boot <- apply(boot_candidate, 2, stats::sd)
    candidate_ci_boot <- t(vapply(seq_len(K), function(k)
      get_ci_from_ests(candidate_estimates[k], candidate_se_boot[k], alpha), numeric(2)))
  }

  names(candidate_estimates) <- names(candidate_se) <-
    names(candidate_se_boot) <- paste0("S", seq_len(K))

  list(estimate = estimate, se = se, ci = ci,
       candidate_estimates = candidate_estimates, candidate_se = candidate_se,
       candidate_ci = candidate_ci, hull_ci = hull_ci,
       candidate_se_boot = candidate_se_boot, candidate_ci_boot = candidate_ci_boot,
       reweighted_estimates = point$reweighted, lambda = point$lambda,
       weights = point$weights, n_used = n, n_boot = n_rep,
       common_covariates = common, adj_sets = adj_sets, covariates = X,
       ref_index = ref_index, alpha = alpha,
       settings = list(seed = seed, n_boot = n_boot, n_cores = n_cores))
}
