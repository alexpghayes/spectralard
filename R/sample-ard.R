left_padded_sequence <- function(x) {
  original <- withr::with_options(
    c(scipen = 999),
    as.character(x)
  )

  max_digits <- max(vapply(original, nchar, integer(1)))
  formatC(x, width = max_digits, format = "d", flag = "0")
}

#' Generate data with specified correlation to X
#' @param y Vector to correlate with
#' @param corr Target correlation value
#' @return Vector with specified correlation to y
complement <- function(y, corr = 0.5) {
  x <- rnorm(length(y))
  y.perp <- residuals(lm(x ~ y))
  corr * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - corr^2)
}

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

simulate_ard_data <- function(sim_params) {
  n <- sim_params$n
  k <- sim_params$k
  corr <- sim_params$corr
  num_good_traits <- sim_params$num_good_traits
  num_bad_traits <- sim_params$num_bad_traits

  stopifnot(num_good_traits + num_bad_traits > 0)

  diagonal <- 0.5
  off_diagonal <- 0.05
  B <- matrix(off_diagonal, nrow = k, ncol = k)
  diag(B) <- diagonal
  pi <- rep(1, k) / k
  theta <- rexp(n) + 1

  network_model <- dcsbm(
    theta = theta, # heterogeneous degrees
    B = B,
    pi = pi,
    allow_self_loops = FALSE,
    poisson_edges = FALSE,
    expected_degree = 2 * n^0.6
  )

  s_pop <- svds(network_model)

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
      rnorm(n * num_bad_traits),
      nrow = n,
      ncol = num_bad_traits
    )

    colnames(bad_traits) <- paste0(
      "bad",
      left_padded_sequence(1:num_bad_traits)
    )

    traits <- cbind(traits, bad_traits)
  }

  A <- sample_sparse(network_model)
  Y <- A %*% traits

  list(
    A = A,
    Y = Y,
    traits = traits,
    s_pop = s_pop
  )
}
