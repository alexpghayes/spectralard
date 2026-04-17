# Simulate Aggregated Relational Data (ARD)

Simulate Aggregated Relational Data (ARD)

## Usage

``` r
simulate_ard_data(n, k, corr, num_good_traits, num_bad_traits)
```

## Arguments

- n:

  Number of nodes.

- k:

  Number of communities/latent dimensions.

- corr:

  Correlation between traits and latent positions.

- num_good_traits:

  Number of traits correlated with latent positions.

- num_bad_traits:

  Number of traits uncorrelated with latent positions.

## Value

A list with the following elements:

- `A`: The adjacency matrix.

- `ard`: The ARD matrix.

- `traits`: The node-level traits.

- `s_pop`: The population spectral decomposition.

## Examples

``` r
set.seed(42)

# 1. Simulate ARD data
sim <- simulate_ard_data(
  n = 100, k = 2, corr = 0.7,
  num_good_traits = 2, num_bad_traits = 2
)

# 2. Estimate the network model
estimate <- estimate_spectrum(sim$ard, sim$traits)
estimate
#> 
#> ── Spectral ARD Estimate ───────────────────────────────────────────────────────
#> Nodes: 100
#> Latent dimensions: 4

# 3. Sample a new graph from the estimate
g <- sample_graph(estimate)
g
#> IGRAPH 98a23fd D--- 100 2454 -- 
#> + edges from 98a23fd:
#>   [1]  2->1  3->1  4->1  5->1  6->1  7->1  8->1  9->1 10->1 11->1 12->1 13->1
#>  [13] 14->1 15->1 16->1 18->1 21->1 24->1 25->1 26->1 27->1 28->1 29->1 30->1
#>  [25] 31->1 32->1 33->1 37->1 38->1 40->1 42->1 43->1 44->1 45->1 47->1 48->1
#>  [37] 56->1 59->1 66->1 68->1 71->1 73->1 93->1  1->2  3->2  4->2  5->2  6->2
#>  [49]  7->2  8->2  9->2 10->2 11->2 12->2 13->2 14->2 15->2 16->2 17->2 18->2
#>  [61] 19->2 20->2 21->2 23->2 24->2 25->2 26->2 27->2 29->2 30->2 33->2 37->2
#>  [73] 38->2 40->2 42->2 43->2 45->2 48->2 51->2 54->2 56->2 65->2 76->2 79->2
#>  [85] 80->2 82->2 84->2 87->2  1->3  2->3  4->3  5->3  7->3  8->3  9->3 10->3
#>  [97] 11->3 12->3 13->3 15->3 16->3 21->3 22->3 23->3 24->3 25->3 26->3 27->3
#> + ... omitted several edges
```
