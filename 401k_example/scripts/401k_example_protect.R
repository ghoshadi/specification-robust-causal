# ==============================================================================
# The 401(k) example under Algorithm A.1: unprotected, then with the mean of
# age protected, then with its mean & var protected.  Run from the
# folder above 401k_example:
#   Rscript 401k_example/scripts/401k_example_protect.R
# Fits are written to 401k_example/results, figures to 401k_example/plots.
# ==============================================================================
rm(list = ls())

source('./utils.R')
source('./protect_covariates.R')

RESULTS <- "401k_example/results"
PLOTS   <- "401k_example/plots"

SPECS <- list(
  list(tag = "unprot",   label = "no protection",
       vars = character(0), fun = NULL),
  list(tag = "mean",     label = "mean of age protected",
       vars = "age",        fun = NULL),
  list(tag = "moments2", label = "mean & var of age protected",
       vars = character(0),
       fun = function(xc) {
         m <- cbind(xc$age, xc$age^2); colnames(m) <- c("age", "age^2"); m }))

VARS <- c("age", "educ", "income")
BINS <- list(25, "unit", 25)
LEG  <- c("topright", "topleft", "topright")

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

run_collection <- function(label, adj_sets) {
  dir.create(RESULTS, showWarnings = FALSE)
  dir.create(PLOTS,   showWarnings = FALSE)

  outs <- lapply(SPECS, function(s) {
    cat("\n================ ", label, ", ", s$label,
        " ================\n", sep = "")
    o <- specification_robust_protect(
      response = response, treatment = treatment,
      covariates = df[, unique(unlist(adj_sets))],
      adj_sets = adj_sets, protect_vars = s$vars, protect_fun = s$fun,
      verbose = TRUE)
    saveRDS(o, file.path(RESULTS, sprintf("401k_%s_%s.rds", s$tag, label)))
    o
  })
  names(outs) <- vapply(SPECS, `[[`, character(1), "tag")

  f <- file.path(PLOTS, sprintf("401k_protect_%s_weights.pdf", label))
  pdf(f, width = 12, height = 3.9 * length(outs))
  plot_hists(df, lapply(outs, `[[`, "weights"), vars = VARS,
             labels = vapply(SPECS, `[[`, character(1), "label"),
             breaks = BINS, legend.pos = LEG)
  invisible(dev.off())
  cat("Wrote ", f, "\n", sep = "")
  outs
}

fits <- list(set1 = run_collection("set1", adj_sets_set1),
             set2 = run_collection("set2", adj_sets_set2))

cat("\n================ SUMMARY ================\n")
cat(sprintf("%-6s %-34s %9s %8s %22s %11s\n", "", "protected", "estimate",
            "s.e.", "spec-robust CI", "reduction"))
for (nm in names(fits)) for (i in seq_along(SPECS)) {
  o <- fits[[nm]][[i]]
  cat(sprintf("%-6s %-34s %9.4f %8.4f  [%8.4f, %8.4f] %10.1f%%\n",
              nm, SPECS[[i]]$label, o$estimate, o$se, o$ci[1], o$ci[2],
              100 * (1 - diff(o$ci)/diff(o$convex_hull_ci))))
}

cat("\n---- the protected moments ----\n")
cat(sprintf("%-6s %-34s %-10s %14s %14s %14s\n", "", "protected", "variable",
            "original", "reweighted", "difference/sd"))
for (nm in names(fits)) for (i in seq_along(SPECS)) {
  ps <- fits[[nm]][[i]]$protected_summary
  if (is.null(ps) || nrow(ps) == 0) next
  for (j in seq_len(nrow(ps)))
    cat(sprintf("%-6s %-34s %-10s %14.4f %14.4f %14.2e\n", nm,
                SPECS[[i]]$label, ps$variable[j], ps$original_mean[j],
                ps$weighted_mean[j], ps$std_difference[j]))
}
