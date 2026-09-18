# ==============================================================================
# The moment ladder for the 401(k) example: mixed moments of the common
# covariates up to order k are protected, k = 0, 1, 2, ...  Run from the folder
# above 401k_example:
#   Rscript 401k_example/scripts/moment_ladder.R
# The fit is written to 401k_example/results, the figures to 401k_example/plots.
# ==============================================================================
rm(list = ls())
library(specrobust)

WQ_USED <- eval(formals(specrobust)$weight_trim)

.LADDER_FINAL_SOURCED <- TRUE
source('./401k_example/scripts/ladder_final_plot.R')

RESULTS <- Sys.getenv("SPECROBUST_RESULTS", unset = "401k_example/results")
PLOTS   <- Sys.getenv("SPECROBUST_PLOTS",   unset = "401k_example/plots")

NF        <- suppressWarnings(as.integer(Sys.getenv("LADDER_FOLDS", unset = NA)))
BCT       <- Sys.getenv("LADDER_BC",     unset = "")
NUREG     <- suppressWarnings(as.numeric(Sys.getenv("LADDER_NUREG", unset = NA)))
KMAX      <- as.integer(Sys.getenv("LADDER_KMAX",  unset = "6"))
D_MAX     <- as.integer(Sys.getenv("LADDER_DMAX",  unset = "40"))
SE_JUMP   <- 2.0     # stop when the s.e. doubles against the previous order
ESS_FLOOR <- 0.05    # stop when the pre-trimming effective sample size falls

utils::data("pension", package = "hdm", envir = environment())
df <- as.data.frame(pension)
names(df)[names(df) == "inc"]     <- "income"
names(df)[names(df) == "fsize"]   <- "family_size"
names(df)[names(df) == "marr"]    <- "marital_status"
names(df)[names(df) == "twoearn"] <- "two_earner"
names(df)[names(df) == "hown"]    <- "home_ownership"
names(df)[names(df) == "db"]      <- "defined_pension"
names(df)[names(df) == "pira"]    <- "ira_participation"

response  <- df$net_tfa / 1000
treatment <- df$e401

base_covs <- c("age", "educ", "family_size", "marital_status", "two_earner",
               "home_ownership", "defined_pension")

COLLECTIONS <- list(
  set1 = list(
    adj_sets = list(
      c("age", "educ"),
      c("age", "educ", "family_size", "income"),
      c("age", "educ", "family_size", "marital_status", "two_earner", "income"),
      c("age", "educ", "family_size", "marital_status", "two_earner", "income",
        "home_ownership", "ira_participation", "defined_pension")),
    protect = c("age", "educ")),
  set2 = list(
    adj_sets = list(
      base_covs,
      c(base_covs, "income"),
      c(base_covs, "ira_participation"),
      c(base_covs, "income", "ira_participation")),
    protect = c("age", "educ", "income")))

nrm <- function(u) sqrt(mean(u^2))

compositions <- function(deg, p) {
  if (p == 1L) return(matrix(deg, 1L, 1L))
  do.call(rbind, lapply(0:deg, function(v) cbind(v, compositions(deg - v, p - 1L))))
}
monomial_exponents <- function(p, k) {
  if (k < 1) return(matrix(integer(0), 0L, p))
  E <- do.call(rbind, lapply(seq_len(k), function(deg) compositions(deg, p)))
  dimnames(E) <- NULL
  E[order(rowSums(E), -E[, 1]), , drop = FALSE]
}
mono <- function(Z, e) {
  out <- rep(1, nrow(Z))
  for (j in seq_along(e)) if (e[j] > 0) out <- out * Z[, j]^e[j]
  out
}
mono_name <- function(e, nms) {
  kp <- e > 0
  paste(sprintf("%s^%d", nms[kp], e[kp]), collapse = ":")
}

onb_nested <- function(Xp, kmax, tol = 1e-8) {
  nms <- colnames(Xp); p <- ncol(Xp); n <- nrow(Xp)
  Z <- vapply(seq_len(p), function(j) {
    x <- as.numeric(Xp[[j]]); (x - mean(x))/stats::sd(x) }, numeric(n))
  E <- monomial_exponents(p, kmax); deg <- as.integer(rowSums(E))
  P <- matrix(0, n, 0); keep <- integer(0)
  for (i in seq_len(nrow(E))) {
    v0 <- mono(Z, E[i, ]); v0 <- v0 - mean(v0); s0 <- nrm(v0)
    if (!is.finite(s0) || s0 <= 0) next
    v <- v0
    if (ncol(P) > 0) for (pass in 1:2) {
      v <- as.numeric(v - P %*% (crossprod(P, v)/n)); v <- v - mean(v) }
    s <- nrm(v)
    if (is.finite(s) && s > tol*s0) { P <- cbind(P, v/s); keep <- c(keep, i) }
  }
  colnames(P) <- apply(E[keep, , drop = FALSE], 1, mono_name, nms = nms)
  list(Q = P,
       dim_k = vapply(seq_len(kmax), function(k) sum(deg[keep] <= k), integer(1)),
       req_k = vapply(seq_len(kmax), function(k) sum(deg      <= k), integer(1)))
}

bal_avg <- function(Xp, w, deg) {
  if (deg < 1) return(NA_real_)
  E  <- monomial_exponents(ncol(Xp), deg)
  Xm <- as.matrix(vapply(Xp, as.numeric, numeric(nrow(Xp))))
  ok <- !is.na(w)
  mean(abs(vapply(seq_len(nrow(E)), function(i) {
    xj <- mono(Xm, E[i, ])
    sdj <- stats::sd(xj); if (!is.finite(sdj) || sdj <= 0) sdj <- 1
    (sum(w[ok] * xj[ok])/sum(w[ok]) - mean(xj))/sdj }, numeric(1))))
}

kl_div  <- function(w) { w <- w[!is.na(w)]; w <- w/mean(w); mean(ifelse(w > 0, w*log(w), 0)) }
ess_of  <- function(w) { w <- w[!is.na(w)]; sum(w)^2/sum(w^2) }

run_ladder <- function(label, adj_sets, protect_req) {
  common  <- Reduce(intersect, adj_sets)
  protect <- intersect(protect_req, common)
  dropped <- setdiff(protect_req, common)

  cat("\n================================ ", disp(label),
      " ================================\n", sep = "")
  cat("weight_trim: ", WQ_USED, "   num_folds: ",
      if (is.na(NF)) "implementation default" else NF,
      "   bias_corr: ", if (nzchar(BCT)) BCT else "default",
      "   nu_regularize: ", if (is.na(NUREG)) "default" else NUREG, "\n", sep = "")
  cat("common     : ", paste(common, collapse = ", "), "\n", sep = "")
  cat("requested  : ", paste(protect_req, collapse = ", "), "\n", sep = "")
  cat("protected  : ", paste(protect, collapse = ", "), "\n", sep = "")
  if (length(dropped))
    cat("DROPPED    : ", paste(dropped, collapse = ", "),
        "  -- not in the intersection, so no function of it is admissible\n", sep = "")
  if (!length(protect)) { cat("nothing admissible to protect; ladder skipped\n"); return(NULL) }

  Xp <- df[, protect, drop = FALSE]
  ob <- onb_nested(Xp, KMAX)
  cat(sprintf("dim V_k    : %s   (of %s requested; the gap is linear dependence\n",
              paste(ob$dim_k, collapse = ", "), paste(ob$req_k, collapse = ", ")))
  cat("             on the sample support)\n")

  covs <- df[, unique(unlist(adj_sets))]
  rows <- list(); W <- list()
  se_prev <- NA_real_; stop_at <- NA_integer_; stop_why <- ""

  for (k in 0:KMAX) {
    d <- if (k == 0) 0L else ob$dim_k[k]
    if (d > D_MAX) { stop_at <- k; stop_why <- sprintf("dim V_%d = %d exceeds D_MAX = %d", k, d, D_MAX); break }

    pf <- if (d == 0) NULL else local({ QQ <- ob$Q[, seq_len(d), drop = FALSE]; function(xc) QQ })

    t0  <- proc.time()[["elapsed"]]
    fargs <- list(response, treatment, covs, adj_sets, protect_fun = pf,
                  verbose = FALSE)
    if (!is.na(NF)) fargs$num_folds <- NF
    if (nzchar(BCT)) fargs$bias_corr <- BCT != "0"
    if (!is.na(NUREG)) fargs$nu_regularize <- NUREG
    fit <- tryCatch(do.call(specrobust, fargs),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "err"))
    el <- proc.time()[["elapsed"]] - t0

    if (inherits(fit, "err")) {
      cat(sprintf("k=%d  d=%3d  FAILED: %s\n", k, d, fit$msg))
      stop_at <- k; stop_why <- paste("the fit failed:", fit$msg); break
    }

    red <- 100 * (1 - diff(fit$ci)/diff(fit$hull_ci))
    rows[[length(rows) + 1L]] <- data.frame(
      k = k, d = d, estimate = fit$estimate, se = fit$se,
      ci_lo = fit$ci[1], ci_hi = fit$ci[2], width = diff(fit$ci),
      reduction = red,
      hull_lo = fit$hull_ci[1], hull_hi = fit$hull_ci[2],
      hull_width = diff(fit$hull_ci),
      ess_frac = ess_of(fit$weights)/fit$n_used,
      kl = kl_div(fit$weights_raw),
      lam_norm = mean(sqrt(rowSums(fit$lambda_by_fold^2))),
      max_w = max(fit$weights, na.rm = TRUE),
      n_used = fit$n_used, secs = el)
    W[[length(rows)]] <- fit$weights_raw
    r <- rows[[length(rows)]]
    cat(sprintf("k=%d  d=%3d  est %8.4f  se %6.4f  width %6.4f  red %5.1f%%  ESS %4.1f%%  KL %.4f  %4.0fs\n",
                r$k, r$d, r$estimate, r$se, r$width, r$reduction,
                100*r$ess_frac, r$kl, r$secs))

    if (!is.na(se_prev) && r$se > SE_JUMP * se_prev) {
      stop_at <- k; stop_why <- sprintf("the s.e. rose by more than %.1fx against order %d", SE_JUMP, k - 1); break }
    if (r$ess_frac < ESS_FLOOR) {
      stop_at <- k; stop_why <- sprintf("the effective sample size fell below %.0f%%", 100*ESS_FLOOR); break }
    se_prev <- r$se
  }

  out <- do.call(rbind, rows)
  out$bal_1   <- vapply(seq_len(nrow(out)), function(i)
    bal_avg(Xp, W[[i]], 1L),   numeric(1))
  out$bal_all <- vapply(seq_len(nrow(out)), function(i)
    bal_avg(Xp, W[[i]], KMAX), numeric(1))
  out$bal_1_unprot   <- bal_avg(Xp, W[[1]], 1L)
  out$bal_all_unprot <- bal_avg(Xp, W[[1]], KMAX)
  if (!is.na(stop_at)) cat(sprintf("\nladder stopped before order %d: %s\n", stop_at, stop_why))
  attr(out, "stop_at") <- stop_at; attr(out, "stop_why") <- stop_why
  attr(out, "protect") <- protect;  attr(out, "dropped")  <- dropped
  attr(out, "dim_k")   <- ob$dim_k
  attr(out, "weights") <- W
  attr(out, "Xp")      <- Xp
  attr(out, "trim")    <- WQ_USED
  out
}

LAD <- lapply(names(COLLECTIONS), function(lab)
  run_ladder(lab, COLLECTIONS[[lab]]$adj_sets, COLLECTIONS[[lab]]$protect))
names(LAD) <- names(COLLECTIONS)

cat("\n\n================ LADDER SUMMARY ================\n")
for (lab in names(LAD)) {
  L <- LAD[[lab]]; if (is.null(L)) next
  cat("\n", disp(lab), "   protected: ", paste(attr(L, "protect"), collapse = ", "),
      if (length(attr(L, "dropped")))
        paste0("   (dropped: ", paste(attr(L, "dropped"), collapse = ", "), ")") else "",
      "\n", sep = "")
  cat(sprintf("%2s %4s %9s %7s %7s %9s %10s %10s %7s %8s\n",
              "k", "d", "estimate", "s.e.", "width", "reduction",
              "bal 1st", sprintf("bal 1-%d", KMAX), "ESS", "KL"))
  for (i in seq_len(nrow(L)))
    cat(sprintf("%2d %4d %9.4f %7.4f %7.4f %8.1f%% %10.2e %10.2e %6.1f%% %8.4f\n",
                L$k[i], L$d[i], L$estimate[i], L$se[i], L$width[i],
                L$reduction[i], L$bal_1[i], L$bal_all[i], 100*L$ess_frac[i], L$kl[i]))
  cat(sprintf("%2s %4s %9s %7s %7s %9s %10.2e %10.2e\n", "--", "", "no protection",
              "", "", "", L$bal_1_unprot[1], L$bal_all_unprot[1]))
  b <- which.min(L$width)
  cat(sprintf("narrowest interval at k = %d (d = %d): [%.4f, %.4f], %.1f%% narrower than the hull\n",
              L$k[b], L$d[b], L$ci_lo[b], L$ci_hi[b], L$reduction[b]))
}

dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS,   showWarnings = FALSE, recursive = TRUE)
f <- file.path(RESULTS, "moment_ladder.rds")
saveRDS(LAD, f); cat("\nWrote ", f, "\n", sep = "")
write_final_figures(LAD, PLOTS)
