# Specification-robust Causal Inference
Which Covariates to Adjust for? Specification-robust Causal Inference in Observational Studies, as proposed by Ghosh and Rothenhäusler (2025+).

In observational causal inference, domain knowledge often leaves multiple covariate adjustments plausible, yet which sets satisfy ignorability is untestable. Different adjustment sets can yield conflicting estimates of the average treatment effect, and standard remedies (adjusting for their union or intersection, or reporting the union or convex hull of confidence intervals) can fail or produce intervals whose width does not vanish with sample size. We propose a specification-robust procedure that returns a single point estimate and a confidence interval that is valid as long as at least one candidate adjustment set is valid and has width shrinking at the parametric $n^{-1/2}$ rate. Our approach mirrors how trimming and overlap weighting handle overlap violations: We shift the target to a reweighted population, closest in KL-divergence to the original population, for which credible, specification-robust inference is feasible. We also provide diagnostic plots to assess the population shift and an extension to protect any function of the covariates used for reweighting, similar to calipers in matching.

The development version of our package `specrobust` can be installed using devtools:

```R
devtools::install_github("ghoshadi/specrobust")
```

Example usage:

```R
library(specrobust)
set.seed(42)
n = 5000
X1 = rnorm(n); X2 = rnorm(n)
A = as.numeric(runif(n) <= 1/(exp(5*X1 + 5*X2) + 1))
Y = A * (1 + X1 - 5*X2) + 4*X2 + rnorm(n)
out = specrobust(Y, A, data.frame(X1, X2), list("X1", c("X1", "X2")))
print(out)
plot(out)
plot(out, type = "weights")
```

The contrasts are estimated by cross-fitted generalized random forests by
default. Passing `reg_mode = "lm"` estimates them by linear regression with
treatment-covariate interactions instead and bootstraps the standard errors:

```R
out_lm = specrobust(Y, A, data.frame(X1, X2), list("X1", c("X1", "X2")),
                    reg_mode = "lm", n_boot = 1000)
```

Any variable in the intersection of the adjustment sets, or any known function
of those variables, can be held fixed under the reweighting:

```R
out_prot = specrobust(Y, A, data.frame(X1, X2), list("X1", c("X1", "X2")),
                      protect_vars = "X1")
```

#### References
Aditya Ghosh and Dominik Rothenhäusler.
<b>Which Covariates to Adjust for? Specification-robust Causal Inference in Observational Studies.</b>, [arXiv preprint arXiv:2505.08729](https://arxiv.org/abs/2505.08729).
