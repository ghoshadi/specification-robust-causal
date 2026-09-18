# ==============================================================================
# The 401(k) example, on two candidate collections, at the implementation's
# defaults.  Run from the folder above 401k_example:
#   Rscript 401k_example/scripts/401k_example.R
# Fits are written to 401k_example/results, figures to 401k_example/plots.
# ==============================================================================
rm(list = ls())

library(specrobust)

results_dir <- "401k_example/results"
plots_dir   <- "401k_example/plots"

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
  cat("\n================ ", label, " ================\n", sep = "")
  out <- specrobust(
    response = response, treatment = treatment,
    covariates = df[, unique(unlist(adj_sets))],
    adj_sets = adj_sets, weight_trim = 0.05, aipw_trim = 0.01,
    verbose = TRUE)

  dir.create(results_dir, showWarnings = FALSE)
  dir.create(plots_dir,   showWarnings = FALSE)
  saveRDS(out, file.path(results_dir, sprintf("401k_default_%s.rds", label)))

  f <- file.path(plots_dir, sprintf("401k_default_%s_weights.pdf", label))
  pdf(f, width = 12, height = 3.9)
  plot(out, covariates = c("age", "educ", "income"),
       breaks = list(25, "unit", 25),
       legend.pos = c("topright", "topleft", "topright"))
  invisible(dev.off())
  cat("Wrote ", f, "\n", sep = "")
  out
}

out_set1 <- run_collection("set1", adj_sets_set1)
out_set2 <- run_collection("set2", adj_sets_set2)

cat("\n================ SUMMARY ================\n")
cat(sprintf("%-6s %9s %8s %22s %11s\n",
            "", "estimate", "s.e.", "spec-robust CI", "reduction"))
for (nm in c("set1", "set2")) {
  o <- get(paste0("out_", nm))
  cat(sprintf("%-6s %9.4f %8.4f  [%8.4f, %8.4f] %10.1f%%\n", nm, o$estimate, o$se,
              o$ci[1], o$ci[2], 100 * (1 - diff(o$ci)/diff(o$hull_ci))))
}
