library(testthat)
library(spectralard)

test_that("estimate_spectrum and sample_graph work", {
  n <- 100
  k <- 2
  sim <- simulate_ard_data(
    n = n, k = k, corr = 0.7,
    num_good_traits = 2, num_bad_traits = 2
  )

  expect_s4_class(sim$Y, "dgeMatrix")
  expect_equal(nrow(sim$Y), n)

  estimate <- estimate_spectrum(sim$Y, sim$traits)
  expect_s3_class(estimate, "SpectralArdEstimate")
  expect_equal(nrow(estimate$Utilde), n)
  expect_equal(ncol(estimate$Utilde), ncol(sim$Y))

  # Test Bernoulli sampling without self-loops
  g <- sample_graph(estimate, allow_self_loops = FALSE, poisson_edges = FALSE)
  expect_s3_class(g, "igraph")
  expect_equal(igraph::vcount(g), n)
  expect_false(igraph::any_loop(g))

  # Test Poisson sampling with self-loops
  g_pois <- sample_graph(estimate, allow_self_loops = TRUE, poisson_edges = TRUE)
  expect_s3_class(g_pois, "igraph")
  expect_equal(igraph::vcount(g_pois), n)
})

test_that("estimate_spectrum handles invalid inputs", {
  expect_error(estimate_spectrum(matrix(NA), matrix(1)))
  expect_error(estimate_spectrum(matrix(1), matrix(NA)))
  expect_error(estimate_spectrum(matrix(1, 2, 2), matrix(1, 2, 3)))
})
