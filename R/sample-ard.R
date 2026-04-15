#' Left-pad an integer sequence with zeros
#'
#' @param x A vector of integers.
#' @return A character vector of the same length as `x`, left-padded with zeros.
#' @keywords internal
left_padded_sequence <- function(x) {
  original <- withr::with_options(
    c(scipen = 999),
    as.character(x)
  )

  max_digits <- max(vapply(original, nchar, integer(1)))
  formatC(x, width = max_digits, format = "d", flag = "0")
}

#' Generate data with specified correlation to a vector
#'
#' @param y Vector to correlate with.
#' @param corr Target correlation value (between -1 and 1).
#' @return A vector with the specified correlation to `y`.
#' @export
#' @examples
#' y <- rnorm(100)
#' x <- complement(y, corr = 0.8)
#' cor(x, y)
complement <- function(y, corr = 0.5) {
  x <- stats::rnorm(length(y))
  y.perp <- stats::residuals(stats::lm(x ~ y))
  corr * stats::sd(y.perp) * y + y.perp * stats::sd(y) * sqrt(1 - corr^2)
}

#' Generate traits correlated with latent positions
#'
#' @param U A matrix of latent positions.
#' @param corr A numeric vector (or scalar) of correlations.
#' @param num_traits The number of traits to generate.
#' @return A matrix of traits.
#' @keywords internal
generate_traits <- function(U, corr, num_traits) {
  k <- ncol(U)

  if (length(corr) == 1) {
    corr <- rep(corr, num_traits)
  } else {
    # must specify corr either once as a scalar, or as
    # a vector for all traits all at once
    stopifnot(length(corr) == num_traits)
  }

  traits <- matrix(0, nrow = nrow(U), ncol = num_traits)

  for (trait in 1:num_traits) {
    # repeat traits if num_traits > k
    trait_index <- trait %% k
    if (trait_index == 0) {
      trait_index <- k
    }
    u <- U[, trait_index, drop = TRUE]
    traits[, trait] <- complement(u, corr = corr[trait])
  }

  traits
}

#' Simulate Aggregated Relational Data (ARD)
#'
#' @param n Number of nodes.
#' @param k Number of communities/latent dimensions.
#' @param corr Correlation between traits and latent positions.
#' @param num_good_traits Number of traits correlated with latent positions.
#' @param num_bad_traits Number of traits uncorrelated with latent positions.
#' @return A list with the following elements:
#'   - `A`: The adjacency matrix.
#'   - `Y`: The ARD matrix.
#'   - `traits`: The node-level traits.
#'   - `s_pop`: The population spectral decomposition.
#' @export
#' @examples
#' sim <- simulate_ard_data(
#'   n = 100, k = 2, corr = 0.7,
#'   num_good_traits = 2, num_bad_traits = 2
#' )
simulate_ard_data <- function(n, k, corr, num_good_traits, num_bad_traits) {
  stopifnot(num_good_traits + num_bad_traits > 0)

  diagonal <- 0.5
  off_diagonal <- 0.05
  B <- matrix(off_diagonal, nrow = k, ncol = k)
  diag(B) <- diagonal
  pi <- rep(1, k) / k
  theta <- stats::rexp(n) + 1

  network_model <- fastRG::dcsbm(
    theta = theta, # heterogeneous degrees
    B = B,
    pi = pi,
    allow_self_loops = FALSE,
    poisson_edges = FALSE,
    expected_degree = 2 * n^0.6
  )

  # RSpectra::svds might be better but for now let's assume it's available
  # or use fastRG equivalent if exists. fastRG models have a custom svds method.
  s_pop <- RSpectra::svds(network_model, k = k)

  X <- s_pop$u %*% diag(sqrt(s_pop$d))

  traits <- matrix(nrow = n, ncol = 0)

  if (num_good_traits > 0) {
    good_traits <- generate_traits(X, corr = corr, num_traits = num_good_traits)

    colnames(good_traits) <- paste0(
      "good",
      left_padded_sequence(1:num_good_traits)
    )

    traits <- cbind(traits, good_traits)
  }

  if (num_bad_traits > 0) {
    bad_traits <- matrix(
      stats::rnorm(n * num_bad_traits),
      nrow = n,
      ncol = num_bad_traits
    )

    colnames(bad_traits) <- paste0(
      "bad",
      left_padded_sequence(1:num_bad_traits)
    )

    traits <- cbind(traits, bad_traits)
  }

  A <- fastRG::sample_sparse(network_model)
  Y <- A %*% traits

  list(
    A = A,
    Y = Y,
    traits = traits,
    s_pop = s_pop
  )
}
