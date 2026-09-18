#' Summarize a specrobust object
#'
#' @param object specrobust object
#' @param ... Additional arguments (currently ignored).
#' @return A data frame with one row per interval: the individual adjustment
#'   sets, their convex hull, and the specification-robust interval. With
#'   `reg_mode = "lm"` each adjustment set appears twice, once with the linear
#'   regression standard error and once with the bootstrap standard error.
#' @export
summary.specrobust = function(object, ...) {
  K = length(object$adj_sets)
  lab = paste0("S", seq_len(K))

  blocks = list(data.frame(
    Estimate = object$candidate_estimates,
    SE = object$candidate_se,
    `CI Lower` = object$candidate_ci[, 1],
    `CI Upper` = object$candidate_ci[, 2],
    row.names = if (object$reg_mode == "lm") paste0(lab, " (lm)") else lab,
    check.names = FALSE))

  if (object$reg_mode == "lm")
    blocks[[length(blocks) + 1]] = data.frame(
      Estimate = object$candidate_estimates,
      SE = object$candidate_se_boot,
      `CI Lower` = object$candidate_ci_boot[, 1],
      `CI Upper` = object$candidate_ci_boot[, 2],
      row.names = paste0(lab, " (bootstrap)"),
      check.names = FALSE)

  blocks[[length(blocks) + 1]] = data.frame(
    Estimate = NA_real_, SE = NA_real_,
    `CI Lower` = object$hull_ci[1], `CI Upper` = object$hull_ci[2],
    row.names = "Convex hull", check.names = FALSE)

  blocks[[length(blocks) + 1]] = data.frame(
    Estimate = object$estimate, SE = object$se,
    `CI Lower` = object$ci[1], `CI Upper` = object$ci[2],
    row.names = "Specification-robust", check.names = FALSE)

  do.call(rbind, blocks)
}

#' Print a specrobust object
#'
#' @param x specrobust object
#' @param digits number of digits to print
#' @param ... Additional arguments passed to print methods.
#' @return `x`, invisibly.
#' @export
print.specrobust = function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Specification-robust causal inference (reg_mode = \"", x$reg_mode, "\")\n\n", sep = "")
  cat("Candidate adjustment sets:\n")
  for (k in seq_along(x$adj_sets))
    cat(paste0("  S", k, ": ", paste(x$adj_sets[[k]], collapse = ", "), "\n"))
  cat(paste0("Shared covariates: ", paste(x$common_covariates, collapse = ", "), "\n"))
  cat(paste0("Observations used: ", x$n_used, "\n"))
  cat(paste0("Confidence level: ", (1 - x$alpha) * 100, "%", "\n"))
  if (!is.null(x$protect_d) && x$protect_d > 0)
    cat(paste0("Protected: ", paste(x$protect_names, collapse = ", "), "\n"))
  if (is.finite(x$width_reduction))
    cat(paste0("Width reduction against the convex hull: ",
               signif(100 * x$width_reduction, 3), "%\n"))
  cat("\n")
  print(summary(x), digits = digits, ...)
  if (!is.null(x$protected_summary) && nrow(x$protected_summary) > 0) {
    cat("\nProtected moments, original against reweighted mean:\n")
    print(x$protected_summary, digits = digits, row.names = FALSE)
  }
  invisible(x)
}

#' Plot a specrobust object
#'
#' @param x specrobust object
#' @param type Which plot to draw. `"covariates"`, the default, overlays the
#'   covariate distributions under the original and the reweighted population,
#'   one panel per covariate, at most three panels per row. `"weights"` draws
#'   a histogram of the transfer weights.
#' @param covariates Which covariates to panel when `type = "covariates"`.
#'   `"common"`, the default, uses the intersection of the adjustment sets, the
#'   covariates the weights are a function of. `"all"` uses every supplied
#'   covariate. A character vector selects those covariates by name.
#' @param breaks Number of histogram breaks, or `"unit"` for unit-width bins.
#'   For `type = "covariates"` a vector or list supplies one value per panel.
#' @param hcol Fill for the original population.
#' @param hcol2 Fill for the reweighted population.
#' @param legend.pos Legend position, recycled or given once per panel.
#' @param legend.text Legend labels.
#' @param ... Additional arguments (currently ignored).
#' @return `NULL`, invisibly.
#' @export
plot.specrobust = function(x, type = c("covariates", "weights"),
                           covariates = "common", breaks = 25,
                           hcol = "orange", hcol2 = "skyblue2",
                           legend.pos = "topleft",
                           legend.text = c("original population",
                                           "new target population"), ...) {
  type = match.arg(type)
  w = x$weights
  ok = !is.na(w)
  w = w[ok]

  if (type == "weights") {
    op = graphics::par(mar = c(4.0, 4.2, 2.4, 0.8), mgp = c(2.3, 0.8, 0))
    on.exit(graphics::par(op))
    b = if (identical(breaks, "unit")) "Sturges" else breaks
    graphics::hist(w, breaks = b, xlab = "Transfer weight", ylab = "Frequency",
                   main = NULL, col = grDevices::adjustcolor(hcol2, alpha.f = 0.6),
                   border = grDevices::adjustcolor(hcol2, alpha.f = 0.2))
    graphics::box()
    return(invisible(NULL))
  }

  avail = colnames(x$covariates)
  vars = if (identical(covariates, "common")) x$common_covariates
         else if (identical(covariates, "all")) avail
         else as.character(unlist(covariates))
  bad = setdiff(vars, avail)
  if (length(bad) > 0)
    stop("Covariates not available in the fit: ", paste(bad, collapse = ", "))
  if (length(vars) == 0) stop("No covariates to plot.")

  data = x$covariates[ok, , drop = FALSE]
  per = function(a, j) if (length(a) == 1) a[[1]] else a[[j]]

  as_num = function(v, nm) {
    if (is.numeric(v)) return(as.numeric(v))
    z = suppressWarnings(as.numeric(as.character(v)))
    if (anyNA(z)) stop("Covariate ", nm, " is not numeric and cannot be plotted.")
    z
  }

  cells = lapply(seq_along(vars), function(j) {
    xv = as_num(data[[vars[j]]], vars[j])
    b = per(breaks, j)
    if (identical(b, "unit"))
      b = seq(floor(min(xv)) - 0.5, ceiling(max(xv)) + 0.5, by = 1)
    h = graphics::hist(xv, breaks = b, plot = FALSE)
    bw = diff(h$breaks)
    wc = as.numeric(tapply(w, cut(xv, h$breaks, include.lowest = TRUE), sum))
    wc[is.na(wc)] = 0
    list(breaks = h$breaks, xlim = range(xv),
         d = list(h$counts/sum(h$counts)/bw, wc/sum(w)/bw))
  })

  nc = min(3, length(vars))
  nr = ceiling(length(vars)/nc)
  op = graphics::par(mfrow = c(nr, nc), mar = c(4.0, 4.2, 2.4, 0.8),
                     mgp = c(2.3, 0.8, 0))
  on.exit(graphics::par(op))
  col = grDevices::adjustcolor(c(hcol, hcol2), alpha.f = 0.6)

  for (j in seq_along(vars)) {
    z = cells[[j]]; b = z$breaks
    graphics::plot(NA, xlim = z$xlim, ylim = range(unlist(z$d)),
                   xlab = vars[j], ylab = "Density")
    graphics::legend(per(legend.pos, j), legend = legend.text, fill = col,
                     col = col, bty = "n")
    for (i in 1:2)
      graphics::rect(b[-length(b)], 0, b[-1], z$d[[i]], col = col[i],
                     border = grDevices::adjustcolor(c(hcol, hcol2)[i], alpha.f = 0.2))
    graphics::box()
  }
  invisible(NULL)
}
