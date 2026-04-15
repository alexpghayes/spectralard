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

- `Y`: The ARD matrix.

- `traits`: The node-level traits.

- `s_pop`: The population spectral decomposition.

## Examples

``` r
sim <- simulate_ard_data(
  n = 100, k = 2, corr = 0.7,
  num_good_traits = 2, num_bad_traits = 2
)
```
