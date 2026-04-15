SpectralArdEstimate <- function(Utilde, Sigmatilde) {
  structure(
    list(Utilde = Utilde, Sigmatilde = Sigmatilde),
    class = "SpectralArdEstimate"
  )
}

#' Estimate a low-rank network model from Aggregated Relational Data (ARD)
#'
#' @param ard A `matrix` containing Aggregated Relational Data. Each row of `ard`
#'   corresponds to one node in the network, and each column represents an
#'   aggregated trait. The order of the nodes (rows) and aggregated traits (columns)
#'   should match `traits`. Missing data is not allowed.
#'
#' @param traits A `matrix` containing node-level trait information. Each row of `traits`
#'   corresponds to one node in the network, and each column represents a particular trait.
#'   The order of the nodes (rows) and aggregated traits (columns) should match `ard`.
#'   Missing data is not allowed.
#'
#' @returns A `SpectralArdEstimate` object.
#'
#' @export
#' @examples
#'
#' sim <- simulate_ard_data(n = 100, k = 2, corr = 0.7,
#'   num_good_traits = 2, num_bad_traits = 2)
#'
#' estimate <- estimate_spectrum(sim$Y, sim$traits)
#' estimate
estimate_spectrum <- function(ard, traits) {
  if (anyNA(ard)) {
    cli::cli_abort("{.arg ard} must be matrix without any `NA` entries.")
  }

  if (anyNA(traits)) {
    cli::cli_abort("{.arg traits} must be matrix without any `NA` entries.")
  }

  if (ncol(ard) != ncol(traits) || nrow(ard) != nrow(traits)) {
    cli::cli_abort(
      "{.arg ard} and {.arg traits} must have the same number of rows and columns"
    )
  }
  s_ard <- svd(ard)

  Utilde <- s_ard$u

  UtY <- crossprod(Utilde, ard)
  UtW <- crossprod(Utilde, traits)

  Ntilde <- tcrossprod(UtY, UtW)
  Dtilde <- tcrossprod(UtW, UtW)
  Sigmatilde <- as.matrix(Ntilde %*% solve(Dtilde))

  SpectralArdEstimate(Utilde, Sigmatilde)
}
