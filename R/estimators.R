SpectralArdEstimate <- function(Utilde, Sigmatilde) {
  structure(
    list(Utilde = Utilde, Sigmatilde = Sigmatilde),
    class = "SpectralArdEstimate"
  )
}

#' @export
#' @method print SpectralArdEstimate
print.SpectralArdEstimate <- function(x, ...) {
  cli::cli_h1("Spectral ARD Estimate")
  cli::cli_bullets(c(
    "i" = "Nodes: {.val {nrow(x$Utilde)}}",
    "i" = "Latent dimensions: {.val {ncol(x$Utilde)}}",
    "v" = "Estimated components (access with {.code $}):",
    "$" = "{.field Utilde}: {.val {nrow(x$Utilde)}} x {.val {ncol(x$Utilde)}} matrix",
    "$" = "{.field Sigmatilde}: {.val {ncol(x$Utilde)}} x {.val {ncol(x$Utilde)}} matrix"
  ))
  invisible(x)
}

#' Sample a graph from an estimated spectral ARD model
#'
#' @param object A `SpectralArdEstimate` object.
#' @param allow_self_loops Logical. Should self-loops be allowed? Defaults to `FALSE`.
#' @param poisson_edges Logical. Should the number of edges between nodes be
#'   sampled from a Poisson distribution? If `FALSE` (the default), edges are
#'   sampled from a Bernoulli distribution.
#' @param ... Currently ignored.
#'
#' @return An `igraph` object.
#' @export
#' @inherit estimate_spectrum examples
sample_graph <- function(
  object,
  allow_self_loops = FALSE,
  poisson_edges = FALSE,
  ...
) {
  UseMethod("sample_graph")
}

#' @export
sample_graph.SpectralArdEstimate <- function(
  object,
  allow_self_loops = FALSE,
  poisson_edges = FALSE,
  ...
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort(
      "The {.pkg igraph} package is required to use {.fn sample_graph}."
    )
  }

  n <- nrow(object$Utilde)

  # To avoid an n x n dense matrix, we sample in chunks.
  # We target a chunk size that keeps memory usage reasonable (~500MB).
  # A 65M element double matrix is about 500MB.
  chunk_size <- max(1, floor(65e6 / n))

  # Precompute USigma to avoid redundant work in the loop
  USigma <- object$Utilde %*% object$Sigmatilde
  U <- object$Utilde

  edge_list <- list()

  for (start in seq(1, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n)
    n_chunk <- end - start + 1

    # Compute chunk of intensities: P[i, j] = (USigma[i, ] %*% t(U))[i, j]
    # P_chunk is n_chunk x n
    P_chunk <- tcrossprod(USigma[start:end, , drop = FALSE], U)

    # For both Bernoulli and Poisson, we need non-negative rates/probabilities.
    # Estimated values may be negative; we clamp them to 0.
    P_chunk[P_chunk < 0] <- 0

    if (poisson_edges) {
      # In the Poisson case, intensities can be > 1.
      # Sample edge counts from Poisson distribution.
      adj_vec <- stats::rpois(n_chunk * n, lambda = P_chunk)
    } else {
      # In the Bernoulli case, clamp probabilities to [0, 1].
      P_chunk[P_chunk > 1] <- 1
      adj_vec <- stats::rbinom(n_chunk * n, size = 1, prob = P_chunk)
    }
    rm(P_chunk)

    # Find indices of edges
    hits <- which(adj_vec > 0)

    if (length(hits) > 0) {
      counts <- adj_vec[hits]
      rm(adj_vec)

      # Convert linear indices to (row, col)
      # Col-major: index = (col-1) * n_chunk + row
      rows_in_chunk <- ((hits - 1) %% n_chunk) + 1
      cols <- ((hits - 1) %/% n_chunk) + 1
      rows <- rows_in_chunk + start - 1

      # Handle self-loops
      if (!allow_self_loops) {
        valid <- rows != cols
        rows <- rows[valid]
        cols <- cols[valid]
        counts <- counts[valid]
      }

      if (length(rows) > 0) {
        # For multi-graphs (Poisson), we repeat the edges based on their counts.
        # igraph::graph_from_edgelist handles multi-edges.
        edge_list[[length(edge_list) + 1]] <- cbind(
          rep(rows, times = counts),
          rep(cols, times = counts)
        )
      }
    } else {
      rm(adj_vec)
    }
  }

  if (length(edge_list) == 0) {
    return(igraph::make_empty_graph(n, directed = TRUE))
  }

  all_edges <- do.call(rbind, edge_list)
  igraph::graph_from_edgelist(all_edges, directed = TRUE)
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
#' set.seed(42)
#'
#' # 1. Simulate ARD data
#' sim <- simulate_ard_data(
#'   n = 100, k = 2, corr = 0.7,
#'   num_good_traits = 2, num_bad_traits = 2
#' )
#'
#' # 2. Estimate the network model
#' estimate <- estimate_spectrum(sim$ard, sim$traits)
#' estimate
#'
#' estimate_spectrum(sim$ard, sim$traits, khat = 2)
#'
#' # 3. Sample a new graph from the estimate
#' g <- sample_graph(estimate)
#' g
#'
estimate_spectrum <- function(ard, traits, khat = NULL) {
  if ((!is.matrix(ard) && !inherits(ard, "Matrix")) || anyNA(ard)) {
    cli::cli_abort("{.arg ard} must be matrix without any `NA` entries.")
  }

  if ((!is.matrix(traits) && !inherits(traits, "Matrix")) || anyNA(traits)) {
    cli::cli_abort("{.arg traits} must be matrix without any `NA` entries.")
  }

  if (ncol(ard) != ncol(traits) || nrow(ard) != nrow(traits)) {
    cli::cli_abort(
      "{.arg ard} and {.arg traits} must have the same number of rows and columns"
    )
  }

  if (is.null(khat)) {
    khat <- ncol(traits)
  }

  s_ard <- svd(ard, nu = khat, nv = khat)

  Utilde <- s_ard$u

  UtY <- crossprod(Utilde, ard)
  UtW <- crossprod(Utilde, traits)

  Ntilde <- tcrossprod(UtY, UtW)
  Dtilde <- tcrossprod(UtW, UtW)
  Sigmatilde <- as.matrix(Ntilde %*% solve(Dtilde))

  SpectralArdEstimate(Utilde, Sigmatilde)
}
