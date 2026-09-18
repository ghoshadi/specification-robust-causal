# ==============================================================================
# The 401(k) example under Algorithm A.1: unprotected, then with the mean of
# age protected, then with its mean & var protected.  Run from the
# folder above 401k_example:
#   Rscript 401k_example/scripts/401k_example_protect.R
# Fits are written to 401k_example/results, figures to 401k_example/plots.
# ==============================================================================
rm(list = ls())

library(specrobust)
if (!requireNamespace("hdm", quietly = TRUE)) install.packages("hdm")

results_dir <- "401k_example/results"
plots_dir   <- "401k_example/plots"

specs <- list(
  list(tag = "unprot",   label = "no protection",
       vars = character(0), fun = NULL),
  list(tag = "mean",     label = "mean of age protected",
       vars = "age",        fun = NULL),
  list(tag = "moments2", label = "mean & var of age protected",
       vars = character(0),
       fun = function(xc) {
         m <- cbind(xc$age, xc$age^2); colnames(m) <- c("age", "age^2"); m }))

utils::data("pension", package = "hdm", envir = environment())
df <- as.data.frame(pension)
names(df)[names(df) == "inc"]     <- "income"
names(df)[names(df) == "fsize"]   <- "family_size"
names(df)[names(df) == "marr"]    <- "marital_status"
names(df)[names(df) == "twoearn"] <- "two_earner"
names(df)[names(df) == "hown"]    <- "home_ownership"
names(df)[names(df) == "db"]      <- "defined_pension"
names(df)[names(df) == "pira"]    <- "ira_participation"

response  <- df$net_tfa / 1000          # thousands of dollars
treatment <- df$e401

base_covs <- c("age", "educ", "family_size", "marital_status", "two_earner",
               "home_ownership", "defined_pension")

adj_sets_set1 <- list(
  c("age", "educ"),
  c("age", "educ", "family_size", "income"),
  c("age", "educ", "family_size", "marital_status", "two_earner", "income"),
  c("age", "educ", "family_size", "marital_status", "two_earner", "income",
    "home_ownership", "ira_participation", "defined_pension"))

adj_sets_set2 <- list(
  base_covs,
  c(base_covs, "income"),
  c(base_covs, "ira_participation"),
  c(base_covs, "income", "ira_participation"))

stack_weight_grid <- function(data, weights, vars, labels = NULL, breaks = 25,
                       legend.pos = "topleft", hcol = "orange",
                       hcol2 = "skyblue2",
                       legend.text = c("original population",
                                       "new target population")) {
  if (!is.list(weights)) weights <- list(weights)
  nr <- length(weights); nc <- length(vars)
  per <- function(a, j) if (length(a) == 1) a[[1]] else a[[j]]

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
    if (!is.null(labels) && j == 1)
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

run_collection <- function(label, adj_sets) {
  dir.create(results_dir, showWarnings = FALSE)
  dir.create(plots_dir,   showWarnings = FALSE)

  outs <- lapply(specs, function(s) {
    cat("\n================ ", label, ", ", s$label,
        " ================\n", sep = "")
    o <- specrobust(
      response = response, treatment = treatment,
      covariates = df[, unique(unlist(adj_sets))],
      adj_sets = adj_sets, protect_vars = s$vars, protect_fun = s$fun,
      weight_trim = 0.05, aipw_trim = 0.01, verbose = TRUE)
    saveRDS(o, file.path(results_dir, sprintf("401k_%s_%s.rds", s$tag, label)))
    o
  })
  names(outs) <- vapply(specs, `[[`, character(1), "tag")

  f <- file.path(plots_dir, sprintf("401k_protect_%s_weights.pdf", label))
  pdf(f, width = 12, height = 3.9 * length(outs))
  stack_weight_grid(df, lapply(outs, `[[`, "weights"),
                    vars = c("age", "educ", "income"),
                    labels = vapply(specs, `[[`, character(1), "label"),
                    breaks = list(25, "unit", 25),
                    legend.pos = c("topright", "topleft", "topright"))
  invisible(dev.off())
  cat("Wrote ", f, "\n", sep = "")
  outs
}

fits <- list(set1 = run_collection("set1", adj_sets_set1),
             set2 = run_collection("set2", adj_sets_set2))

cat("\n================ SUMMARY ================\n")
cat(sprintf("%-6s %-34s %9s %8s %22s %11s\n", "", "protected", "estimate",
            "s.e.", "spec-robust CI", "reduction"))
for (nm in names(fits)) for (i in seq_along(specs)) {
  o <- fits[[nm]][[i]]
  cat(sprintf("%-6s %-34s %9.4f %8.4f  [%8.4f, %8.4f] %10.1f%%\n",
              nm, specs[[i]]$label, o$estimate, o$se, o$ci[1], o$ci[2],
              100 * (1 - diff(o$ci)/diff(o$hull_ci))))
}

cat("\n---- the protected moments ----\n")
cat(sprintf("%-6s %-34s %-10s %14s %14s %14s\n", "", "protected", "variable",
            "original", "reweighted", "difference/sd"))
for (nm in names(fits)) for (i in seq_along(specs)) {
  ps <- fits[[nm]][[i]]$protected_summary
  if (is.null(ps) || nrow(ps) == 0) next
  for (j in seq_len(nrow(ps)))
    cat(sprintf("%-6s %-34s %-10s %14.4f %14.4f %14.2e\n", nm,
                specs[[i]]$label, ps$variable[j], ps$original_mean[j],
                ps$weighted_mean[j], ps$std_difference[j]))
}
