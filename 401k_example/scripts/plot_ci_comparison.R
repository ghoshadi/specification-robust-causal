# ==============================================================================
# The confidence-interval figure of the 401(k) example, reproducing the TikZ
# picture of Figure 4 of the paper: the K per-candidate AIPW intervals in
# black, their convex hull in red, the specification-robust intervals in green,
# point estimates as blue rings, the adjustment sets described on the right.
# Lengths below are in the TikZ coordinate system, x = 0.7 cm and y = 3.5 cm.
#   Rscript 401k_example/scripts/plot_ci_comparison.R
# Reads 401k_example/results, writes 401k_example/plots.
# ==============================================================================
rm(list = ls())

TAG     <- Sys.getenv("SPECROBUST_TAG",     unset = "unprot")
RESULTS <- Sys.getenv("SPECROBUST_RESULTS", unset = "401k_example/results")
PLOTS   <- Sys.getenv("SPECROBUST_PLOTS",   unset = "401k_example/plots")

EXTRA <- Sys.getenv("SPECROBUST_EXTRA",
  unset = paste("mean=mean of age protected",
                "moments2=mean & var of age protected", sep = ";"))
EXTRA <- if (nzchar(EXTRA)) lapply(strsplit(EXTRA, ";", fixed = TRUE)[[1]],
  function(x) list(tag   = trimws(sub("=.*", "", x)),
                   label = trimws(sub("^[^=]*=", "", x)))) else list()

STEM <- Sys.getenv("SPECROBUST_STEM",
                   unset = if (length(EXTRA)) "protect" else TAG)
BASE <- Sys.getenv("SPECROBUST_BASE",
                   unset = if (length(EXTRA)) "no protection" else "Proposed CI")

LABELS <- c(age = "age", educ = "education", family_size = "family size",
            marital_status = "marital status", two_earner = "two-earner household",
            home_ownership = "home-ownership", defined_pension = "defined pension",
            ira_participation = "IRA participation", income = "income")
pretty_names <- function(v) paste(ifelse(v %in% names(LABELS), LABELS[v], v),
                                  collapse = ", ")

set_descriptions <- function(adj_sets) {
  K <- length(adj_sets)
  if (all(vapply(seq_len(K - 1), function(k)
          all(adj_sets[[k]] %in% adj_sets[[k + 1]]), logical(1))))
    return(vapply(seq_len(K), function(k)
      if (k == 1) pretty_names(adj_sets[[1]]) else
        paste0("+ ", pretty_names(setdiff(adj_sets[[k]], adj_sets[[k - 1]]))),
      character(1)))
  top <- which(vapply(seq_len(K), function(k)
    all(vapply(adj_sets, function(s) all(s %in% adj_sets[[k]]), logical(1))),
    logical(1)))
  base <- if (length(top)) adj_sets[[top[1]]] else adj_sets[[1]]
  vapply(seq_len(K), function(k) {
    if (length(top) && k == top[1]) return(pretty_names(base))
    drop <- setdiff(base, adj_sets[[k]])
    if (!length(top)) paste0("+ ", pretty_names(setdiff(adj_sets[[k]], base)))
    else sprintf("S%d \\ {%s}", top[1], pretty_names(drop))
  }, character(1))
}

wrap_descriptions <- function(desc) {
  ne   <- nzchar(desc)
  room <- vapply(seq_along(desc), function(i) {
    r <- 1L; j <- i + 1L
    while (j <= length(ne) && !ne[j]) { r <- r + 1L; j <- j + 1L }
    r
  }, 1L)
  ok <- function(w) all(vapply(seq_along(desc),
    function(i) length(strwrap(desc[i], w)) <= room[i], logical(1)))
  w <- Find(ok, seq(42, 90, by = 2))
  lapply(lapply(desc, strwrap, width = if (is.null(w)) 90 else w),
         function(v) lapply(v, mathify))
}

mathify <- function(s) {
  m <- regmatches(s, regexec("^S([0-9]+) \\\\ \\{(.*)\\}$", s))[[1]]
  if (length(m) == 3L)
    bquote(S[.(as.integer(m[2]))] ~ "\\" ~ .(paste0("{", m[3], "}"))) else s
}

COL_AIPW <- "black"
COL_HULL <- "#CC4C4C"   # grayred
COL_SR   <- "#007300"   # darkgreen
COL_PT   <- "#0050FF"   # myblue

UX  <- 0.7/2.54         # inches per x unit
UY  <- 3.5/2.54         # inches per y unit
PT  <- 1/72.27          # inches per TeX point
FS  <- 8                # \scriptsize
LW      <- 2.0/72.27*96 # line width=2pt, in R lwd units (1/96 inch)
LW_AX   <- 0.8/72.27*96 # thick
LW_TICK <- 0.4/72.27*96 # TikZ default
R_OUT   <- 3.00*PT      # blue disc
R_IN    <- 1.25*PT      # white disc inside it
SEP     <- 0.3333*FS*PT # TikZ default inner sep

ring <- function(x, y) {
  th <- seq(0, 2*pi, length.out = 72)
  for (r in list(c(R_OUT, COL_PT), c(R_IN, "white")))
    polygon(x + as.numeric(r[[1]])/UX * cos(th), y + as.numeric(r[[1]])/UY * sin(th),
            col = r[[2]], border = NA)
}

geom <- function(o, rng, extra = list()) {
  z   <- stats::qnorm(1 - o$alpha/2)
  K   <- length(o$aipw_estimates)
  ex  <- lapply(extra, function(e) e$fit)
  lo  <- c(o$aipw_estimates - z * o$aipw_se, o$convex_hull_ci[1], o$ci[1],
           vapply(ex, function(q) q$ci[1], 0))
  hi  <- c(o$aipw_estimates + z * o$aipw_se, o$convex_hull_ci[2], o$ci[2],
           vapply(ex, function(q) q$ci[2], 0))
  est <- c(o$aipw_estimates, NA, o$estimate,
           vapply(ex, function(q) q$estimate, 0))

  left  <- c(lapply(seq_len(K), function(k) bquote(AIPW ~ (S[.(k)]))),
             list("Convex hull", BASE),
             lapply(extra, function(e) e$label))
  right <- c(wrap_descriptions(c(set_descriptions(o$adj_sets), "", "")),
             rep(list(character(0)), length(extra)))

  nl <- vapply(right, function(v) max(1L, length(v)), 1L)
  ne <- vapply(right, function(v) length(v) > 0 && any(v != ""), TRUE)
  nr <- length(lo)
  y  <- numeric(nr)
  for (i in seq_len(nr - 1) + 1) y[i] <- y[i-1] - 0.1*if (ne[i]) nl[i-1] else 1
  y  <- y - y[nr] - 0.1

  ticks <- pretty(rng, 5)
  xlo   <- rng[1] - 0.08 * diff(rng)
  xhi   <- rng[2] + 0.05 * diff(rng)
  ticks <- ticks[ticks >= xlo & ticks <= xhi]

  grDevices::pdf(file = nullfile(), family = "Times", pointsize = FS)
  wl <- max(vapply(left, function(s) strwidth(s, units = "inches"), 0))
  wr <- max(0, unlist(lapply(right, function(v)
    vapply(seq_along(v), function(k)
      strwidth(v[[k]], units = "inches") + (k > 1)*0.5*UX, 0))))
  ht <- strheight("0", units = "inches")
  invisible(grDevices::dev.off())

  ybot <- -0.22 - (SEP + ht)/UY
  ytop <- y[1] + (R_OUT + 0.5*ht)/UY
  list(lo = lo, hi = hi, est = est, y = y, K = K, left = left, right = right,
       col = c(rep(COL_AIPW, K), COL_HULL, rep(COL_SR, 1 + length(extra))),
       xlo = xlo, xhi = xhi, ticks = ticks, ybot = ybot, ytop = ytop,
       mai_l = wl + 0.05*UX + SEP, mai_r = wr + 0.10*UX + SEP,
       pw = (xhi - xlo)*UX, ph = (ytop - ybot)*UY)
}

panel <- function(g, slack = 0) {
  par(mai = c(0, g$mai_l, 0, g$mai_r + slack), ps = FS, family = "Times",
      xpd = NA, lend = "butt")
  plot(NA, xlim = c(g$xlo, g$xhi), ylim = c(g$ybot, g$ytop),
       xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")

  arrows(g$xlo, -0.2, g$xhi, -0.2, code = 3, length = 0.055, angle = 22,
         lwd = LW_AX)
  segments(g$ticks, -0.18, g$ticks, -0.22, lwd = LW_TICK)
  text(g$ticks, -0.22 - SEP/UY, g$ticks, adj = c(0.5, 1))

  segments(g$lo, g$y, g$hi, g$y, lwd = LW, col = g$col)
  for (i in seq_along(g$y)) if (!is.na(g$est[i])) ring(g$est[i], g$y[i])

  for (i in seq_along(g$left))
    text(g$xlo - 0.05 - SEP/UX, g$y[i], g$left[[i]], adj = c(1, 0.5))
  for (i in seq_along(g$right)) {
    v <- g$right[[i]]
    for (k in seq_along(v))
      text(g$xhi + 0.10 + SEP/UX + (k > 1)*0.5, g$y[i] - 0.1*(k - 1), v[[k]],
           adj = c(0, 0.5))
  }
}

rd   <- function(tag, lab)
  readRDS(file.path(RESULTS, sprintf("401k_%s_%s.rds", tag, lab)))
SETS <- strsplit(Sys.getenv("SPECROBUST_SETS", unset = "set1,set2"), ",")[[1]]
fits <- lapply(stats::setNames(SETS, SETS), function(lab) rd(TAG, lab))
have <- function(tag, lab)
  file.exists(file.path(RESULTS, sprintf("401k_%s_%s.rds", tag, lab)))
exs  <- lapply(names(fits), function(lab)
  lapply(Filter(function(e) have(e$tag, lab), EXTRA),
         function(e) c(e, list(fit = rd(e$tag, lab)))))
names(exs) <- names(fits)

rng <- range(unlist(lapply(fits, function(o) {
  z <- stats::qnorm(1 - o$alpha/2)
  c(o$aipw_estimates - z*o$aipw_se, o$aipw_estimates + z*o$aipw_se,
    o$convex_hull_ci, o$ci)
})), unlist(lapply(exs, function(v)
       vapply(v, function(e) e$fit$ci, numeric(2)))))
XR <- Sys.getenv("SPECROBUST_XRANGE", unset = "")
if (nzchar(XR)) rng <- as.numeric(strsplit(XR, ",", fixed = TRUE)[[1]])
gs  <- lapply(names(fits), function(lab)
  geom(fits[[lab]], rng = rng, extra = exs[[lab]]))
names(gs) <- names(fits)

wid <- max(vapply(gs, function(g) g$mai_l + g$pw + g$mai_r, 0))
for (lab in names(gs)) {
  g <- gs[[lab]]
  dir.create(PLOTS, showWarnings = FALSE)
  f <- file.path(PLOTS, sprintf("401k_%s_%s_ci.pdf", STEM, lab))
  pdf(f, width = wid, height = g$ph, family = "Times", pointsize = FS)
  panel(g, slack = wid - (g$mai_l + g$pw + g$mai_r)); invisible(dev.off())
  cat("Wrote ", f, "\n", sep = "")
}
cat("\n")

for (lab in names(fits)) {
  o <- fits[[lab]]; K <- length(o$aipw_estimates)
  z <- stats::qnorm(1 - o$alpha/2)
  cat("================ ", lab, " ================\n", sep = "")
  cat(sprintf("%-36s %9s %8s %22s %8s %8s\n", "", "estimate", "s.e.",
              "95% CI", "width", "red."))
  for (k in seq_len(K))
    cat(sprintf("%-36s %9.4f %8.4f  [%8.4f, %8.4f] %8.4f\n",
                sprintf("AIPW (S%d)", k), o$aipw_estimates[k], o$aipw_se[k],
                o$aipw_estimates[k] - z*o$aipw_se[k],
                o$aipw_estimates[k] + z*o$aipw_se[k], 2*z*o$aipw_se[k]))
  cat(sprintf("%-36s %9s %8s  [%8.4f, %8.4f] %8.4f\n", "Convex hull", "", "",
              o$convex_hull_ci[1], o$convex_hull_ci[2], diff(o$convex_hull_ci)))
  for (e in c(list(list(label = BASE, fit = o)), exs[[lab]])) {
    q <- e$fit
    cat(sprintf("%-36s %9.4f %8.4f  [%8.4f, %8.4f] %8.4f %7.1f%%\n",
                if (is.character(e$label)) e$label else BASE,
                q$estimate, q$se, q$ci[1], q$ci[2], diff(q$ci),
                100 * (1 - diff(q$ci)/diff(o$convex_hull_ci))))
  }
  cat("\n")
}
