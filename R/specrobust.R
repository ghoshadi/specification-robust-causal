#' specrobust: Specification-Robust Causal Inference
#'
#' Estimation and inference for the average treatment effect when several
#' candidate covariate adjustment sets are plausible and it is not known which
#' of them satisfies ignorability, as proposed by Ghosh and Rothenhaeusler
#' (2025+). The population is reweighted, as little as possible in
#' Kullback-Leibler divergence, so that all candidate estimands agree, and the
#' effect on that reweighted population is reported against the convex hull of
#' the per-candidate intervals.
#'
#' @param response The outcomes.
#' @param treatment The binary treatment indicator, coded 0 and 1.
#' @param covariates A data frame or matrix of covariates with named columns.
#' @param adj_sets A list of at least two candidate adjustment sets, each a
#'   character vector of column names of `covariates`. Their intersection must
#'   be non-empty; the transfer weights are a function of that intersection.
#' @param reg_mode How the contrasts are estimated. `"grf"` uses cross-fitted
#'   generalized random forests with augmented inverse probability weighting;
#'   `"lm"` uses linear regression with treatment-covariate interactions and
#'   bootstraps the standard errors.
#' @param protect_vars Character vector of variables in the intersection of
#'   the adjustment sets whose means are held fixed by the reweighting. Only
#'   available when `reg_mode = "grf"`.
#' @param protect_fun A function of the shared covariates returning a matrix
#'   with one row per observation, whose column means are held fixed by the
#'   reweighting. Only available when `reg_mode = "grf"`.
#' @param protect_orthonormalize Whether the protected block is orthonormalized
#'   before the tilt is solved.
#' @param alpha Significance level, so the intervals have nominal coverage
#'   `1 - alpha`.
#' @param ref_index Which candidate serves as the reference. A numeric index,
#'   or one of `"first"`, `"largest"`, `"smallest"`.
#' @param num_folds Number of cross-fitting folds. `reg_mode = "grf"` only.
#' @param num_trees Number of trees in every forest, both layers.
#'   `reg_mode = "grf"` only.
#' @param propensity_clip Propensity scores are clipped to
#'   `[propensity_clip, 1 - propensity_clip]`. `reg_mode = "grf"` only.
#' @param nu_regularize Ridge parameter of the system that combines the
#'   candidates, shrinking towards equal weights. `reg_mode = "grf"` only.
#' @param bias_corr Whether the second-order bias correction is applied.
#'   `reg_mode = "grf"` only.
#' @param aipw_trim Evaluation rows outside the `aipw_trim` and
#'   `1 - aipw_trim` quantiles of any candidate pseudo-outcome are dropped.
#'   `reg_mode = "grf"` only.
#' @param weight_trim Evaluation rows above the `1 - weight_trim` quantile of
#'   the transfer weights are dropped. `reg_mode = "grf"` only.
#' @param n_boot Number of bootstrap replications, or `0` to fit the point
#'   estimates and the transfer weights only, leaving `se` and `ci` as `NA`.
#'   `reg_mode = "lm"` only.
#' @param n_cores Number of cores used for the bootstrap, by forking, so more
#'   than one core has no effect on Windows. `reg_mode = "lm"` only.
#' @param seed Random seed for reproducibility.
#' @param verbose Whether progress information is printed while fitting.
#' @param return_full Whether the per-fold intermediates are attached to the
#'   result as `$full`. `reg_mode = "grf"` only.
#'
#' @return A `specrobust` object. The fields common to both modes are
#'   `estimate`, `se` and `ci` for the reweighted population;
#'   `candidate_estimates`, `candidate_se` and `candidate_ci` for the
#'   individual adjustment sets; `hull_ci` for their convex hull; `weights`
#'   for the transfer weights; and `adj_sets`, `common_covariates`,
#'   `covariates`, `alpha` and `n_used`. With `reg_mode = "grf"` it also
#'   carries `nu`, the per-fold diagnostics, and, when something is protected,
#'   `protected_summary`. With `reg_mode = "lm"` it also carries
#'   `candidate_se_boot` and `candidate_ci_boot`.
#'
#' @references Ghosh, A., & Rothenhaeusler, D. (2025+).
#' Which Covariates to Adjust for? Specification-Robust Causal Inference in
#' Observational Studies.
#'
#' @examples
#' \donttest{
#' # Example 1, two confounders: {X1} is invalid, {X1, X2} is valid
#' set.seed(123)
#' n = 500
#' X1 = rnorm(n); X2 = rnorm(n)
#' A = as.numeric(runif(n) <= 1/(exp(5*X1 + 5*X2) + 1))
#' Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
#' out = specrobust(Y, A, data.frame(X1, X2), list("X1", c("X1", "X2")))
#' print(out)
#' plot(out)
#'
#' # Example 2, M-bias: {X1} is valid, {X1, X2} is invalid
#' set.seed(123)
#' n = 500
#' U1 = rnorm(n); U2 = rnorm(n); X1 = rnorm(n)
#' A = as.numeric(U1 + X1 > 0)
#' X2 = U1 + U2
#' Y = ifelse(A == 1, 1 + X1 - U2 + rnorm(n), 5*U2 + rnorm(n))
#' out = specrobust(Y, A, data.frame(X1, X2), list("X1", c("X1", "X2")),
#'                  reg_mode = "lm", n_boot = 200)
#' print(out)
#' }
#'
#' @export
specrobust <- function(response,
                       treatment,
                       covariates,
                       adj_sets,
                       reg_mode = c("grf", "lm"),
                       protect_vars = character(0),
                       protect_fun = NULL,
                       protect_orthonormalize = FALSE,
                       alpha = 0.05,
                       ref_index = 1,
                       num_folds = 2,
                       num_trees = 400,
                       propensity_clip = 1e-3,
                       nu_regularize = 1e-6,
                       bias_corr = TRUE,
                       aipw_trim = 0,
                       weight_trim = 0,
                       n_boot = 1000,
                       n_cores = 1,
                       seed = 42,
                       verbose = FALSE,
                       return_full = FALSE) {

  reg_mode <- match.arg(reg_mode)
  validate_inputs(response, treatment, covariates, adj_sets)
  stopifnot(length(alpha) == 1, is.finite(alpha), alpha > 0, alpha < 1)

  if (reg_mode == "lm") {
    if (length(protect_vars) > 0 || !is.null(protect_fun))
      stop("Covariate protection is available only with reg_mode = \"grf\".")
    n_boot <- as.integer(n_boot)
    n_cores <- as.integer(n_cores)
    if (!(length(n_boot) == 1) || is.na(n_boot) || n_boot == 1 || n_boot < 0)
      stop("n_boot must be 0, for no bootstrap, or at least 2.")
    stopifnot(length(n_cores) == 1, !is.na(n_cores), n_cores >= 1)
    out <- fit_lm(response, treatment, covariates, adj_sets,
                  ref_index = ref_index, verbose = verbose, alpha = alpha,
                  seed = seed, n_boot = n_boot, n_cores = n_cores)
  } else {
    stopifnot(length(bias_corr) == 1, is.logical(bias_corr), !is.na(bias_corr))
    stopifnot(length(weight_trim) == 1, is.finite(weight_trim),
              weight_trim >= 0, weight_trim < 1)
    stopifnot(length(aipw_trim) == 1, is.finite(aipw_trim),
              aipw_trim >= 0, aipw_trim < 0.5)
    num_folds <- as.integer(num_folds)
    stopifnot(length(num_folds) == 1, !is.na(num_folds), num_folds >= 2)
    out <- fit_grf(response, treatment, covariates, adj_sets,
                   protect_vars = protect_vars, protect_fun = protect_fun,
                   protect_orthonormalize = protect_orthonormalize,
                   ref_index = ref_index, verbose = verbose, alpha = alpha,
                   seed = seed, num_folds = num_folds,
                   propensity_clip = propensity_clip, num_trees = num_trees,
                   nu_regularize = nu_regularize, bias_corr = bias_corr,
                   aipw_trim = aipw_trim, weight_trim = weight_trim,
                   return_full = return_full)
  }

  out$reg_mode <- reg_mode
  out$width_reduction <- 1 - diff(out$ci)/diff(out$hull_ci)
  class(out) <- "specrobust"
  out
}
