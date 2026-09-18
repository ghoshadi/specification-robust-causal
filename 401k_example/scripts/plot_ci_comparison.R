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

results_dir <- "401k_example/results"
plots_dir   <- "401k_example/plots"

tag        <- "unprot"
file_stem  <- "protect"
base_label <- "no protection"
extra <- list(list(tag = "mean",     label = "mean of age protected"),
              list(tag = "moments2", label = "mean & var of age protected"))

var_labels <- c(age = "age", educ = "education", family_size = "family size",
            marital_status = "marital status", two_earner = "two-earner household",
            home_ownership = "home-ownership", defined_pension = "defined pension",
            ira_participation = "IRA participation", income = "income")
pretty_names <- function(v) paste(ifelse(v %in% names(var_labels), var_labels[v], v),
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
    r <- 1; j <- i + 1
    while (j <= length(ne) && !ne[j]) { r <- r + 1; j <- j + 1 }
    r
  }, 1)
  ok <- function(w) all(vapply(seq_along(desc),
    function(i) length(strwrap(desc[i], w)) <= room[i], logical(1)))
  w <- Find(ok, seq(42, 90, by = 2))
  lapply(lapply(desc, strwrap, width = if (is.null(w)) 90 else w),
         function(v) lapply(v, mathify))
}

mathify <- function(s) {
  m <- regmatches(s, regexec("^S([0-9]+) \\\\ \\{(.*)\\}$", s))[[1]]
  if (length(m) == 3)
    bquote(S[.(as.integer(m[2]))] ~ "\\" ~ .(paste0("{", m[3], "}"))) else s
}

col_aipw <- "black"
col_hull <- "#CC4C4C"   # grayred
col_sr   <- "#007300"   # darkgreen
col_pt   <- "#0050FF"   # myblue

ux  <- 0.7/2.54         # inches per x unit
uy  <- 3.5/2.54         # inches per y unit
pt  <- 1/72.27          # inches per TeX point
fs  <- 8                # \scriptsize
lw      <- 2.0/72.27*96 # line width=2pt, in R lwd units (1/96 inch)
lw_ax   <- 0.8/72.27*96 # thick
lw_tick <- 0.4/72.27*96 # TikZ default
r_out   <- 3.00*pt      # blue disc
r_in    <- 1.25*pt      # white disc inside it
inner_sep     <- 0.3333*fs*pt # TikZ default inner sep

ring <- function(x, y) {
  th <- seq(0, 2*pi, length.out = 72)
  for (r in list(c(r_out, col_pt), c(r_in, "white")))
    polygon(x + as.numeric(r[[1]])/ux * cos(th), y + as.numeric(r[[1]])/uy * sin(th),
            col = r[[2]], border = NA)
}

geom <- function(o, rng, extra = list()) {
  z   <- stats::qnorm(1 - o$alpha/2)
  K   <- length(o$candidate_estimates)
  ex  <- lapply(extra, function(e) e$fit)
  lo  <- c(o$candidate_estimates - z * o$candidate_se, o$hull_ci[1], o$ci[1],
           vapply(ex, function(q) q$ci[1], 0))
  hi  <- c(o$candidate_estimates + z * o$candidate_se, o$hull_ci[2], o$ci[2],
           vapply(ex, function(q) q$ci[2], 0))
  est <- c(o$candidate_estimates, NA, o$estimate,
           vapply(ex, function(q) q$estimate, 0))

  left  <- c(lapply(seq_len(K), function(k) bquote(AIPW ~ (S[.(k)]))),
             list("Convex hull", base_label),
             lapply(extra, function(e) e$label))
  right <- c(wrap_descriptions(c(set_descriptions(o$adj_sets), "", "")),
             rep(list(character(0)), length(extra)))

  nl <- vapply(right, function(v) max(1, length(v)), 1)
  ne <- vapply(right, function(v) length(v) > 0 && any(v != ""), TRUE)
  nr <- length(lo)
  y  <- numeric(nr)
  for (i in seq_len(nr - 1) + 1) y[i] <- y[i-1] - 0.1*if (ne[i]) nl[i-1] else 1
  y  <- y - y[nr] - 0.1

  ticks <- pretty(rng, 5)
  xlo   <- rng[1] - 0.08 * diff(rng)
  xhi   <- rng[2] + 0.05 * diff(rng)
  ticks <- ticks[ticks >= xlo & ticks <= xhi]

  grDevices::pdf(file = nullfile(), family = "Times", pointsize = fs)
  wl <- max(vapply(left, function(s) strwidth(s, units = "inches"), 0))
  wr <- max(0, unlist(lapply(right, function(v)
    vapply(seq_along(v), function(k)
      strwidth(v[[k]], units = "inches") + (k > 1)*0.5*ux, 0))))
  ht <- strheight("0", units = "inches")
  invisible(grDevices::dev.off())

  ybot <- -0.22 - (inner_sep + ht)/uy
  ytop <- y[1] + (r_out + 0.5*ht)/uy
  list(lo = lo, hi = hi, est = est, y = y, K = K, left = left, right = right,
       col = c(rep(col_aipw, K), col_hull, rep(col_sr, 1 + length(extra))),
       xlo = xlo, xhi = xhi, ticks = ticks, ybot = ybot, ytop = ytop,
       mai_l = wl + 0.05*ux + inner_sep, mai_r = wr + 0.10*ux + inner_sep,
       pw = (xhi - xlo)*ux, ph = (ytop - ybot)*uy)
}

panel <- function(g, slack = 0) {
  par(mai = c(0, g$mai_l, 0, g$mai_r + slack), ps = fs, family = "Times",
      xpd = NA, lend = "butt")
  plot(NA, xlim = c(g$xlo, g$xhi), ylim = c(g$ybot, g$ytop),
       xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")

  arrows(g$xlo, -0.2, g$xhi, -0.2, code = 3, length = 0.055, angle = 22,
         lwd = lw_ax)
  segments(g$ticks, -0.18, g$ticks, -0.22, lwd = lw_tick)
  text(g$ticks, -0.22 - inner_sep/uy, g$ticks, adj = c(0.5, 1))

  segments(g$lo, g$y, g$hi, g$y, lwd = lw, col = g$col)
  for (i in seq_along(g$y)) if (!is.na(g$est[i])) ring(g$est[i], g$y[i])

  for (i in seq_along(g$left))
    text(g$xlo - 0.05 - inner_sep/ux, g$y[i], g$left[[i]], adj = c(1, 0.5))
  for (i in seq_along(g$right)) {
    v <- g$right[[i]]
    for (k in seq_along(v))
      text(g$xhi + 0.10 + inner_sep/ux + (k > 1)*0.5, g$y[i] - 0.1*(k - 1), v[[k]],
           adj = c(0, 0.5))
  }
}

rd   <- function(tag, lab)
  readRDS(file.path(results_dir, sprintf("401k_%s_%s.rds", tag, lab)))
sets <- c("set1", "set2")
fits <- lapply(stats::setNames(sets, sets), function(lab) rd(tag, lab))
have <- function(tag, lab)
  file.exists(file.path(results_dir, sprintf("401k_%s_%s.rds", tag, lab)))
exs  <- lapply(names(fits), function(lab)
  lapply(Filter(function(e) have(e$tag, lab), extra),
         function(e) c(e, list(fit = rd(e$tag, lab)))))
names(exs) <- names(fits)

rng <- range(unlist(lapply(fits, function(o) {
  z <- stats::qnorm(1 - o$alpha/2)
  c(o$candidate_estimates - z*o$candidate_se, o$candidate_estimates + z*o$candidate_se,
    o$hull_ci, o$ci)
})), unlist(lapply(exs, function(v)
       vapply(v, function(e) e$fit$ci, numeric(2)))))
gs  <- lapply(names(fits), function(lab)
  geom(fits[[lab]], rng = rng, extra = exs[[lab]]))
names(gs) <- names(fits)

wid <- max(vapply(gs, function(g) g$mai_l + g$pw + g$mai_r, 0))
for (lab in names(gs)) {
  g <- gs[[lab]]
  dir.create(plots_dir, showWarnings = FALSE)
  f <- file.path(plots_dir, sprintf("401k_%s_%s_ci.pdf", file_stem, lab))
  pdf(f, width = wid, height = g$ph, family = "Times", pointsize = fs)
  panel(g, slack = wid - (g$mai_l + g$pw + g$mai_r)); invisible(dev.off())
  cat("Wrote ", f, "\n", sep = "")
}
cat("\n")

for (lab in names(fits)) {
  o <- fits[[lab]]; K <- length(o$candidate_estimates)
  z <- stats::qnorm(1 - o$alpha/2)
  cat("================ ", lab, " ================\n", sep = "")
  cat(sprintf("%-36s %9s %8s %22s %8s %8s\n", "", "estimate", "s.e.",
              "95% CI", "width", "red."))
  for (k in seq_len(K))
    cat(sprintf("%-36s %9.4f %8.4f  [%8.4f, %8.4f] %8.4f\n",
                sprintf("AIPW (S%d)", k), o$candidate_estimates[k], o$candidate_se[k],
                o$candidate_estimates[k] - z*o$candidate_se[k],
                o$candidate_estimates[k] + z*o$candidate_se[k], 2*z*o$candidate_se[k]))
  cat(sprintf("%-36s %9s %8s  [%8.4f, %8.4f] %8.4f\n", "Convex hull", "", "",
              o$hull_ci[1], o$hull_ci[2], diff(o$hull_ci)))
  for (e in c(list(list(label = base_label, fit = o)), exs[[lab]])) {
    q <- e$fit
    cat(sprintf("%-36s %9.4f %8.4f  [%8.4f, %8.4f] %8.4f %7.1f%%\n",
                if (is.character(e$label)) e$label else base_label,
                q$estimate, q$se, q$ci[1], q$ci[2], diff(q$ci),
                100 * (1 - diff(q$ci)/diff(o$hull_ci))))
  }
  cat("\n")
}
