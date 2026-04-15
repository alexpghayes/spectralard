# Linear-in-Means Estimation with 2SLS and ARD

``` r
library(spectralard)
library(Matrix)
library(ivreg)
```

## Introduction

In this vignette, we demonstrate how to perform linear-in-means
estimation of peer effects when the network adjacency matrix $G$ is
unknown. Instead of $G$, we only observe Aggregated Relational Data
(ARD), which we use to estimate the network’s spectral structure. We
then use this estimate to construct instruments for a two-stage least
squares (2SLS) estimation using the `ivreg` package.

## The Linear-in-Means Model

Consider the standard linear-in-means model for an outcome $y$:

$$y = \alpha\iota + \beta Gy + \gamma W + \delta GW + \eta X + \epsilon$$

where: - $G$ is the row-normalized adjacency matrix. - $W$ is a vector
of node-level exogenous covariates. - $X$ are the latent positions from
the underlying network model. - $Gy$ represents the endogenous peer
effect. - $GW$ represents the contextual effects.

## Simulation Setup

First, we simulate a network and ARD using the package’s built-in
simulation tools. We use $n = 1000$ nodes and $k = 2$ latent dimensions.

``` r
set.seed(42)
n <- 1000
k <- 2

# Simulate ARD data
sim <- simulate_ard_data(
  n = n, k = k, corr = 0.8,
  num_good_traits = 5, num_bad_traits = 0
)

A <- sim$A
traits <- sim$traits
ard <- sim$ard
X <- sim$s_pop$u %*% diag(sqrt(sim$s_pop$d)) # True latent positions

# True parameters for linear-in-means model
alpha <- 1
beta <- 0.5
gamma <- 2
delta <- 1
eta <- c(0.5, -0.5) # Coefficients for latent positions

# Generate exogenous covariate W
W <- rnorm(n)

# Row-normalize A to get G
row_sums <- rowSums(A)
row_sums[row_sums == 0] <- 1
G <- A / row_sums

# Solve for endogenous y: 
# y = (I - beta * G)^{-1} * (alpha + gamma * W + delta * G * W + X * eta + eps)
eps <- rnorm(n, sd = 0.5)
I <- Diagonal(n)
rhs <- alpha + gamma * W + delta * (G %*% W) + (as.matrix(X) %*% eta) + eps
y <- as.numeric(solve(I - beta * G, rhs))
```

## Spectral ARD Estimation

Next, we estimate the network spectral structure from the ARD. We assume
that the node degrees (row sums of $A$) are known and use them to
normalize our spectral estimate.

``` r
estimate <- estimate_spectrum(ard, traits)

# Extract estimated components
Utilde <- estimate$Utilde
Sigmatilde <- estimate$Sigmatilde

# Construct P_tilde = Utilde %*% Sigmatilde %*% t(Utilde)
# This serves as a plug-in for the unknown adjacency matrix A
P_tilde <- tcrossprod(Utilde %*% Sigmatilde, Utilde)

# Normalize P_tilde using KNOWN degrees to get G_hat
G_hat <- P_tilde / row_sums
```

## 2SLS Estimation

We use ${\widehat{G}}^{2}W$ and $\widehat{G}X$ as instruments for the
endogenous peer effect $Gy$. We also include $W$, $GW$, and the latent
positions $X$ as covariates.

``` r
# Prepare peer effects and instruments
Gy <- as.numeric(G %*% y)
Gw <- as.numeric(G %*% W)
G2w_hat <- as.numeric(G_hat %*% (G_hat %*% W))

# GX latent positions as instruments
Gx_hat <- as.matrix(G_hat %*% X)
colnames(Gx_hat) <- c("Gx1", "Gx2")

# Latent positions as covariates
X1 <- X[, 1]
X2 <- X[, 2]

# Perform 2SLS using ivreg
# Model: y ~ Gy + W + Gw + X1 + X2
# Instruments for Gy: G2w_hat + Gx1 + Gx2
model_2sls <- ivreg(
  y ~ Gy + W + Gw + X1 + X2 | 
    G2w_hat + Gx_hat + W + Gw + X1 + X2
)

summary(model_2sls)
#> 
#> Call:
#> ivreg(formula = y ~ Gy + W + Gw + X1 + X2 | G2w_hat + Gx_hat + 
#>     W + Gw + X1 + X2)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.87373 -0.33386  0.01529  0.33654  1.62678 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  0.60807    0.31066   1.957 0.050586 .  
#> Gy           0.79920    0.21050   3.797 0.000156 ***
#> W            1.97456    0.01660 118.961  < 2e-16 ***
#> Gw           0.36282    0.45950   0.790 0.429954    
#> X1           0.63789    0.08537   7.472 1.73e-13 ***
#> X2          -0.39107    0.09033  -4.329 1.65e-05 ***
#> 
#> Diagnostic tests:
#>                  df1 df2 statistic p-value    
#> Weak instruments   3 992    725.40  <2e-16 ***
#> Wu-Hausman         1 993      0.04   0.842    
#> Sargan             2  NA      0.35   0.839    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.502 on 994 degrees of freedom
#> Multiple R-Squared: 0.9398,  Adjusted R-squared: 0.9395 
#> Wald test:  3102 on 5 and 994 DF,  p-value: < 2.2e-16
```

## Conclusion

By substituting the spectral estimate $\widehat{G}$ for the unknown
row-normalized adjacency matrix $G$, and leveraging known degrees, we
can successfully recover peer effect estimates. Including latent
positions and their estimated neighborhood averages as instruments
significantly aids identification.
