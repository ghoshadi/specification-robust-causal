# ==============================================================================
# Specification-robust Causal Inference (Ghosh & Rothenhaeusler 2026)
# This illustrates the proposed method using the 401(k) real-data example.
# ==============================================================================

rm(list=ls())
source('./utils.R')
source('./main.R')

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

adj_sets <- list(
  c("age", "educ"),
  c("age", "educ", "family_size", "income"),
  c("age", "educ", "family_size", "marital_status", "two_earner", "income"),
  c("age", "educ", "family_size", "marital_status", "two_earner", "income", 
    "home_ownership", "ira_participation", "defined_pension")
)

xvars <- unique(unlist(adj_sets))

# ============================================================
# Run the specification-robust procedure (Algorithms 1 & 2)
# ============================================================

out = specification_robust(
  response = df$net_tfa, 
  treatment = df$e401,
  covariates = df[, xvars], 
  adj_sets = adj_sets,
  verbose = TRUE
)
compare_hists_overlay(df$age, out$weights, xlab="age", draw.curve=F, legend.pos = "topright")
compare_hists_overlay(df$educ, out$weights, xlab="educ", draw.curve=F, legend.pos = "topleft")
compare_hists_overlay(df$inc, out$weights, xlab="inc", draw.curve=F, legend.pos = "topright")
