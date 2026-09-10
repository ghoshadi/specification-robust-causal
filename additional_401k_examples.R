# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# The 401(k) example under alternative collections of candidate adjustment sets.
# ==============================================================================

rm(list=ls())
source('./utils.R')
source('./main.R')

SAVE_PLOTS <- TRUE
OUTDIR     <- "results_more_examples"
SEED       <- 123

# ============================================================
# Preparing the dataset
# ============================================================

library(hdm)
data("pension")
df <- as.data.frame(pension)
names(df)[names(df) == "inc"]     <- "income"
names(df)[names(df) == "fsize"]   <- "family_size"
names(df)[names(df) == "marr"]    <- "marital_status"
names(df)[names(df) == "twoearn"] <- "two_earner"
names(df)[names(df) == "hown"]    <- "home_ownership"
names(df)[names(df) == "db"]      <- "defined_pension"
names(df)[names(df) == "pira"]    <- "ira_participation"

# The nine covariates of Benjamin (2003, J. Public Econ. 87, 1259-1290), also
# used in the 401(k) application of Chernozhukov, Chetverikov, Demirer, Duflo,
# Hansen, Newey and Robins (2018, Econometrics J. 21, C1-C68).
dml_covs <- c("age", "educ", "income", "family_size", "marital_status",
              "two_earner", "defined_pension", "ira_participation",
              "home_ownership")

# The same, minus the two indicators of saving behavior.
base_covs <- setdiff(dml_covs, c("ira_participation", "home_ownership"))

# ============================================================
# Candidate collections
# ============================================================

adj_collections <- list(

  # Nested sequence, as in 401k_example.R.
  "C1_paper" = list(
    c("age", "educ"),
    c("age", "educ", "family_size", "income"),
    c("age", "educ", "family_size", "marital_status", "two_earner", "income"),
    c("age", "educ", "family_size", "marital_status", "two_earner", "income",
      "home_ownership", "ira_participation", "defined_pension")
  ),

  # As C1, with income added to the smallest set, following Poterba, Venti and
  # Wise (1996, J. Econ. Perspect. 10(4), 91-112).
  "C2_paper_income_in_S1" = list(
    c("age", "educ", "income"),
    c("age", "educ", "family_size", "income"),
    c("age", "educ", "family_size", "marital_status", "two_earner", "income"),
    c("age", "educ", "family_size", "marital_status", "two_earner", "income",
      "home_ownership", "ira_participation", "defined_pension")
  ),

  # One candidate per contested covariate, dropped from the nine.
  "C3_leave_one_out" = list(
    dml_covs,
    setdiff(dml_covs, "ira_participation"),
    setdiff(dml_covs, "home_ownership"),
    setdiff(dml_covs, "defined_pension")
  ),

  # Two-by-two factorial over the two indicators of saving behavior, which
  # Benjamin (2003) argues proxy an unobserved taste for saving.
  "C4_saving_proxy_factorial" = list(
    base_covs,
    c(base_covs, "ira_participation"),
    c(base_covs, "home_ownership"),
    c(base_covs, "ira_participation", "home_ownership")
  )
)

if (SAVE_PLOTS) dir.create(OUTDIR, showWarnings = FALSE)

# ============================================================
# Run the specification-robust procedure on each collection
# ============================================================

results <- list()

for (nm in names(adj_collections)) {

  adj_sets <- adj_collections[[nm]]
  xvars    <- unique(unlist(adj_sets))

  cat("\n\n----------------------------------------------------------\n")
  cat("---  ", nm, "\n")
  cat("----------------------------------------------------------\n\n")

  set.seed(SEED)
  out <- specification_robust(
    response   = df$net_tfa,
    treatment  = df$e401,
    covariates = df[, xvars],
    adj_sets   = adj_sets,
    verbose    = TRUE
  )

  results[[nm]] <- out

  if (SAVE_PLOTS) pdf(file.path(OUTDIR, paste0(nm, "_weights.pdf")),
                      width = 6, height = 5)
  compare_hists_overlay(df$age,    out$weights, xlab = "age",
                        draw.curve = FALSE, legend.pos = "topright")
  compare_hists_overlay(df$educ,   out$weights, xlab = "educ",
                        draw.curve = FALSE, legend.pos = "topleft")
  compare_hists_overlay(df$income, out$weights, xlab = "income",
                        draw.curve = FALSE, legend.pos = "topright")
  if (SAVE_PLOTS) dev.off()
}

# ============================================================
# Comparison across collections
# ============================================================

tab <- do.call(rbind, lapply(names(results), function(nm) {
  o <- results[[nm]]
  data.frame(
    collection = nm,
    K          = length(o$adj_sets),
    common     = length(o$common_covariates),
    estimate   = round(o$estimate),
    se         = round(o$se),
    ci_lo      = round(o$ci[1]),
    ci_hi      = round(o$ci[2]),
    ci_width   = round(diff(o$ci)),
    hull_width = round(diff(o$convex_hull_ci)),
    width_redn = sprintf("%.0f%%", 100 * (1 - diff(o$ci) / diff(o$convex_hull_ci))),
    row.names  = NULL
  )
}))

cat("\n\n================ COMPARISON ACROSS COLLECTIONS ================\n")
print(tab, row.names = FALSE)
cat("===============================================================\n")

if (SAVE_PLOTS) {
  saveRDS(results, file.path(OUTDIR, "results.rds"))
  write.csv(tab, file.path(OUTDIR, "comparison.csv"), row.names = FALSE)
}
