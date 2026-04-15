# spectralard

The goal of `spectralard` is to provide efficient spectral estimation of
low-rank network models from Aggregated Relational Data (ARD).

## Installation

You can install the development version of `spectralard` from
[GitHub](https://github.com/alexpghayes/spectralard) with:

``` r
# install.packages("devtools")
devtools::install_github("alexpghayes/spectralard")
```

## Features

- **Efficient Spectral Estimation**: Estimate latent positions and
  community structures directly from ARD using
  [`estimate_spectrum()`](https://alexpghayes.github.io/spectralard/reference/estimate_spectrum.md).
- **Sparse Graph Sampling**: Generate representative graphs from
  estimated models using
  [`sample_graph()`](https://alexpghayes.github.io/spectralard/reference/sample_graph.md),
  with support for large networks through memory-efficient chunked
  sampling.
- **Model-Based Simulation**: Generate synthetic ARD datasets using
  [`simulate_ard_data()`](https://alexpghayes.github.io/spectralard/reference/simulate_ard_data.md).
  This function generates a network from a Degree-Corrected Stochastic
  Block Model (DCSBM) and simulates traits that are correlated with the
  underlying latent positions. ARD is then formed by aggregating these
  traits over the neighborhood of each node.

## Quick Start

This example demonstrates how to simulate ARD data, estimate the
underlying network model, and analyze the results.

``` r
library(spectralard)

# 1. Simulate ARD data
# Uses a DCSBM and generates traits correlated with latent positions
sim <- simulate_ard_data(
  n = 200,              # Number of nodes
  k = 3,                 # Latent dimensions
  corr = 0.8,            # Correlation between traits and latent positions
  num_good_traits = 5,   # Informative traits
  num_bad_traits = 0     # Uninformative traits
)
```

The estimation returns a spectral decomposition of the underlying
network model.

``` r
# 2. Estimate the network model
estimate <- estimate_spectrum(sim$ard, sim$traits)

# The estimate contains the spectral components
estimate$Utilde[1:5, ]
#>            [,1]       [,2]        [,3]        [,4]         [,5]
#> [1,] -0.1415037 0.11690606 -0.05471458 -0.08199877 -0.058363992
#> [2,] -0.1380731 0.11706110 -0.06413892 -0.02666455  0.018530022
#> [3,] -0.1370887 0.11595520 -0.05625102 -0.08994069  0.026554500
#> [4,] -0.1148585 0.11914494 -0.04518377  0.13510443 -0.021374103
#> [5,] -0.1147252 0.09921956 -0.06781349 -0.05768563  0.004314305
estimate$Sigmatilde
#>             [,1]        [,2]       [,3]       [,4]         [,5]
#> [1,] 40.73364425 -0.06716709 -0.7911385   1.133723  -0.06352423
#> [2,] -0.06716709 29.23102358  0.2159438   1.743138  -0.52952347
#> [3,] -0.79113848  0.21594380 26.5554528  -2.805831   2.30528812
#> [4,]  1.13372256  1.74313757 -2.8058308  44.448773 -14.76779592
#> [5,] -0.06352423 -0.52952347  2.3052881 -14.767796 -11.99922151
```

We can use the estimated latent positions `Utilde` to recover community
structure (e.g., via k-means clustering) or sample new graphs.

``` r
# 3. Find blocks using k-means on Utilde
clusters <- kmeans(estimate$Utilde, centers = 3)

# 4. Evaluate clustering accuracy (Adjusted Rand Index)
# A helper to calculate ARI manually to avoid extra dependencies
calc_ari <- function(x, y) {
  tab <- table(x, y)
  sum_comb_n_ij <- sum(choose(tab, 2))
  sum_comb_a_i <- sum(choose(rowSums(tab), 2))
  sum_comb_b_j <- sum(choose(colSums(tab), 2))
  n <- sum(tab)
  comb_n_2 <- choose(n, 2)
  prod_comb <- (sum_comb_a_i * sum_comb_b_j) / comb_n_2
  (sum_comb_n_ij - prod_comb) / (0.5 * (sum_comb_a_i + sum_comb_b_j) - prod_comb)
}

calc_ari(clusters$cluster, sim$model$z)
#> [1] 1

# 5. Sample a new graph from the estimate
g <- sample_graph(
  estimate, 
  allow_self_loops = FALSE, 
  poisson_edges = FALSE
)
g
#> IGRAPH aa0efaa D--- 200 8883 -- 
#> + edges from aa0efaa:
#>  [1]   2->1   3->1   4->1   5->1   7->1   8->1   9->1  10->1  11->1  13->1
#> [11]  14->1  15->1  16->1  17->1  18->1  19->1  21->1  22->1  24->1  25->1
#> [21]  26->1  28->1  29->1  32->1  34->1  35->1  36->1  37->1  38->1  39->1
#> [31]  41->1  43->1  44->1  46->1  47->1  50->1  54->1  56->1  58->1  59->1
#> [41]  61->1  62->1  64->1  65->1  68->1  72->1  73->1  74->1  79->1  80->1
#> [51]  85->1  87->1  94->1 108->1 109->1 115->1 119->1 129->1 132->1 134->1
#> [61] 138->1 144->1 154->1 158->1 160->1 163->1 164->1 167->1 170->1 176->1
#> [71] 184->1 190->1 191->1 192->1 194->1 199->1   1->2   3->2   4->2   5->2
#> [81]   6->2   7->2   8->2   9->2  10->2  11->2  13->2  15->2  16->2  17->2
#> + ... omitted several edges
```

The Adjusted Rand Index (ARI) of 1.0 indicates that the clustering
perfectly recovers the true community structure from the noisy ARD. Even
with only $0.8$ correlation between traits and latent positions, and the
presence of uninformative traits, the spectral estimation accurately
identifies the underlying blocks.
