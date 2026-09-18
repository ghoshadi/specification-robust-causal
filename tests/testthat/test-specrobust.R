example_1 = function(n = 500, seed = 123) {
  set.seed(seed)
  X1 = rnorm(n); X2 = rnorm(n)
  A = as.numeric(runif(n) <= 1/(exp(5*X1 + 5*X2) + 1))
  Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
  list(Y = Y, A = A, X = data.frame(X1, X2), adj_sets = list("X1", c("X1", "X2")))
}

example_2 = function(n = 500, seed = 123) {
  set.seed(seed)
  U1 = rnorm(n); U2 = rnorm(n); X1 = rnorm(n)
  A = as.numeric(U1 + X1 > 0)
  X2 = U1 + U2
  Y = ifelse(A == 1, 1 + X1 - U2 + rnorm(n), 5*U2 + rnorm(n))
  list(Y = Y, A = A, X = data.frame(X1, X2), adj_sets = list("X1", c("X1", "X2")))
}

test_that("grf mode reproduces the reference implementation", {
  for (mk in list(example_1, example_2)) {
    d = mk()
    out = specrobust(d$Y, d$A, d$X, d$adj_sets)
    expect_s3_class(out, "specrobust")
    expect_equal(out$protect_d, 0L)
    expect_length(out$estimate, 1L)
    expect_true(is.finite(out$estimate), is.finite(out$se))
    expect_equal(out$ci, c(out$estimate - stats::qnorm(0.975) * out$se,
                           out$estimate + stats::qnorm(0.975) * out$se))
    expect_equal(sum(out$nu), 1)
    expect_equal(out$hull_ci[1], min(out$candidate_ci[, 1]))
    expect_equal(out$hull_ci[2], max(out$candidate_ci[, 2]))
  }
})

test_that("grf mode is deterministic given the seed", {
  d = example_1()
  a = specrobust(d$Y, d$A, d$X, d$adj_sets)
  b = specrobust(d$Y, d$A, d$X, d$adj_sets)
  expect_identical(a$estimate, b$estimate)
  expect_identical(a$se, b$se)
  expect_identical(a$weights, b$weights)
})

test_that("protection engages only when asked", {
  d = example_1()
  plain = specrobust(d$Y, d$A, d$X, d$adj_sets)
  prot = specrobust(d$Y, d$A, d$X, d$adj_sets, protect_vars = "X1")

  expect_equal(plain$protect_d, 0L)
  expect_equal(prot$protect_d, 1L)
  expect_equal(prot$protect_names, "X1")
  expect_false(isTRUE(all.equal(plain$estimate, prot$estimate)))
  expect_equal(nrow(prot$protected_summary), 1L)
  expect_equal(prot$protected_summary$original_mean, mean(d$X$X1))

  expect_error(specrobust(d$Y, d$A, d$X, d$adj_sets, protect_vars = "X2"),
               "intersection")
  expect_error(specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm",
                          protect_vars = "X1"),
               "reg_mode")
})

test_that("the tilt balances the moments it is solved on", {
  # a feasible constraint: mild confounding, so protecting X1 and forcing the
  # candidates to agree do not pull against each other
  set.seed(123)
  n = 800
  X1 = rnorm(n); X2 = rnorm(n)
  A = as.numeric(runif(n) <= 1/(1 + exp(-(0.5*X1 + 0.5*X2))))
  Y = A * (1 + X1 - X2) + X2 + rnorm(n)
  X = data.frame(X1, X2); adj_sets = list("X1", c("X1", "X2"))

  plain = specrobust(Y, A, X, adj_sets)
  expect_lt(max(abs(plain$balance_train_by_fold)), 1e-4)

  prot = specrobust(Y, A, X, adj_sets, protect_vars = "X1")
  expect_lt(max(abs(prot$balance_train_by_fold)), 1e-4)
  expect_lt(abs(prot$protected_summary$std_difference), 0.05)
})

test_that("lm mode matches the interacted regressions it is built on", {
  d = example_1()
  out = specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm", n_boot = 50)

  f1 = stats::lm(d$Y ~ d$A * d$X$X1)
  f2 = stats::lm(d$Y ~ d$A * (d$X$X1 + d$X$X2))
  expect_equal(unname(out$candidate_estimates[1]),
               unname(summary(f1)$coefficients[2, 1]))
  expect_equal(unname(out$candidate_estimates[2]),
               unname(summary(f2)$coefficients[2, 1]))
  expect_equal(unname(out$candidate_se[1]),
               unname(summary(f1)$coefficients[2, 2]))
  expect_equal(unname(out$candidate_se[2]),
               unname(summary(f2)$coefficients[2, 2]))

  expect_equal(mean(out$weights), 1)
  expect_true(all(is.finite(out$candidate_se_boot)))
})

test_that("the lm bootstrap is reproducible and core-count independent", {
  d = example_1(n = 300)
  a = specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm", n_boot = 40)
  b = specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm", n_boot = 40)
  expect_identical(a$se, b$se)
  expect_identical(a$candidate_se_boot, b$candidate_se_boot)

  if (.Platform$OS.type == "unix") {
    p = specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm", n_boot = 40,
                   n_cores = 2L)
    expect_identical(a$se, p$se)
    expect_identical(a$candidate_se_boot, p$candidate_se_boot)
  }
})

test_that("the bootstrap does not disturb the caller's random stream", {
  d = example_1(n = 200)
  set.seed(99); before = runif(1)
  set.seed(99)
  specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm", n_boot = 20)
  expect_identical(runif(1), before)
})

test_that("summary has the documented shape in both modes", {
  d = example_1(n = 300)
  g = summary(specrobust(d$Y, d$A, d$X, d$adj_sets))
  expect_equal(colnames(g), c("Estimate", "SE", "CI Lower", "CI Upper"))
  expect_equal(rownames(g), c("S1", "S2", "Convex hull", "Specification-robust"))

  l = summary(specrobust(d$Y, d$A, d$X, d$adj_sets, reg_mode = "lm", n_boot = 20))
  expect_equal(rownames(l), c("S1 (lm)", "S2 (lm)", "S1 (bootstrap)",
                              "S2 (bootstrap)", "Convex hull",
                              "Specification-robust"))
  expect_equal(l["S1 (lm)", "Estimate"], l["S1 (bootstrap)", "Estimate"])
})

test_that("plot selects covariates as documented and rejects bad input", {
  d = example_1(n = 200)
  out = specrobust(d$Y, d$A, d$X, d$adj_sets)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(out))
  expect_silent(plot(out, covariates = "all"))
  expect_silent(plot(out, covariates = c("X1", "X2")))
  expect_silent(plot(out, type = "weights"))
  expect_error(plot(out, covariates = "X3"), "not available")
  expect_error(plot(out, type = "nonsense"))
})

test_that("input validation rejects malformed problems", {
  d = example_1(n = 200)
  expect_error(specrobust(d$Y, d$A, d$X, list("X1")))
  expect_error(specrobust(d$Y, d$A, d$X, list("X1", "X9")), "not in covariates")
  expect_error(specrobust(d$Y, d$A, d$X, list("X1", "X2")), "intersection")
  expect_error(specrobust(d$Y, d$A + 1, d$X, d$adj_sets))
})
