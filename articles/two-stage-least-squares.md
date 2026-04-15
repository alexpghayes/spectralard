# Linear-in-Means Estimation with 2SLS and ARD

``` r
library(spectralard)
library(Matrix)
```

## Introduction

In this vignette, we demonstrate how to perform linear-in-means
estimation of peer effects when the network adjacency matrix $A$ is
unknown. Instead of $A$, we only observe Aggregated Relational Data
(ARD), which we use to estimate the network’s spectral structure. We
then use this estimate to construct instruments for a two-stage least
squares (2SLS) estimation.

## The Linear-in-Means Model

Consider the standard linear-in-means model for an outcome $y$:

$$y = \alpha\iota + \beta Gy + \gamma X + \delta GX + \epsilon$$

where $G$ is the row-normalized adjacency matrix, $X$ are node-level
exogenous covariates, and $Gy$ represents the endogenous peer effect.

## Simulation Setup

First, we simulate a network and ARD using the package’s built-in
simulation tools.

``` r
set.seed(42)
n <- 500
k <- 2

# Simulate ARD data
sim <- simulate_ard_data(
  n = n, k = k, corr = 0.8,
  num_good_traits = 5, num_bad_traits = 2
)

A <- sim$A
traits <- sim$traits
Y_ard <- sim$Y

# True parameters for linear-in-means model
alpha <- 1
beta <- 0.5
gamma <- 2
delta <- 1

# Generate exogenous covariate X
X <- rnorm(n)

# Row-normalize A to get G
row_sums <- rowSums(A)
row_sums[row_sums == 0] <- 1
G <- A / row_sums

# Solve for endogenous y: y = (I - beta * G)^{-1} * (alpha + gamma * X + delta * G * X + eps)
eps <- rnorm(n, sd = 0.5)
I <- Diagonal(n)
rhs <- alpha + gamma * X + delta * (G %*% X) + eps
y <- as.numeric(solve(I - beta * G, rhs))
```

## Spectral ARD Estimation

Next, we estimate the network spectral structure from the ARD.

``` r
estimate <- estimate_spectrum(Y_ard, traits)

# Extract estimated components
Utilde <- estimate$Utilde
Sigmatilde <- estimate$Sigmatilde

# Construct P_tilde = Utilde %*% Sigmatilde %*% t(Utilde)
# This serves as a plug-in for the unknown adjacency matrix A
P_tilde <- tcrossprod(Utilde %*% Sigmatilde, Utilde)

# Row-normalize P_tilde to get G_hat
p_row_sums <- rowSums(P_tilde)
p_row_sums[p_row_sums == 0] <- 1
G_hat <- P_tilde / p_row_sums
```

## 2SLS Estimation

We use ${\widehat{G}}^{2}X$ as an instrument for the endogenous peer
effect $Gy$, since we don’t observe $G$.

``` r
# Endogenous variable (using true G as is standard in peer effect literature, 
# but in practice we might only have an estimate)
Gy <- as.numeric(G %*% y)
Gx <- as.numeric(G %*% X)

# Instruments using G_hat
# We use G_hat^2 * X as an instrument for G * y
G2x <- as.numeric(G_hat %*% (G_hat %*% X))

# First stage: Regress Gy on G2x, X, and Gx
first_stage <- lm(Gy ~ G2x + X + Gx)
Gy_hat <- predict(first_stage)

# Second stage: Regress y on Gy_hat, X, and Gx
second_stage <- lm(y ~ Gy_hat + X + Gx)

summary(second_stage)
#> 
#> Call:
#> lm(formula = y ~ Gy_hat + X + Gx)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.39294 -0.34450  0.05311  0.32990  1.19379 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)  
#> (Intercept) -30.0378    68.4703  -0.439   0.6611  
#> Gy_hat       14.5099    30.9363   0.469   0.6393  
#> X             1.6892     0.7262   2.326   0.0204 *
#> Gx          -29.3164    67.5012  -0.434   0.6643  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.4858 on 496 degrees of freedom
#> Multiple R-squared:  0.9465, Adjusted R-squared:  0.9462 
#> F-statistic:  2924 on 3 and 496 DF,  p-value: < 2.2e-16
```

## Conclusion

By substituting the spectral estimate $\widehat{G}$ for the unknown
row-normalized adjacency matrix $G$, we can successfully recover peer
effect estimates even when the underlying network structure is not
directly observed.
