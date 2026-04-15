generate_centering_mx <- function(n) {
  # Construct and return the centering matrix for n data points
  diag(nrow = n) - matrix(rep(1, n * n), c(n, n)) / n
}

construct_chunk <- function(X) {
  # Construct X (I-n^{-1}J)/sqrt(n),
  # the matrix that shows up a bunch in our proofs.
  nxx <- nrow(X)
  Cxx <- generate_centering_mx(nxx)
  Cxx %*% X / sqrt(nxx)
}

construct_cca_part <- function(X, regu = 0) {
  # Construct the product of subspaces for X that would be used if we were
  # to CCA X against some other traits.
  # Since CCA( X,traits ) = SigmaX^{-1/2} \frac{1}{n} X^T M traits Sigmatraits^{-1/2},
  # this means constructing the left and right subspaces of
  # ( n^{-1} X^T X + regu )^{-1/2} X^T M / sqrt(n).
  # fs encodes whether or not to apply the Fisher-Sun covar regularization.

  Xchunk <- construct_chunk(X)
  Xsvd <- svd(Xchunk, ncol(X))

  UX <- Xsvd$u
  VX <- Xsvd$v
  RX <- diag(Xsvd$d, ncol = length(Xsvd$d), nrow = length(Xsvd$d))

  # Otherwise, just use the given regularization.
  dd <- length(Xsvd$d)
  RXreguinv <- diag(1 / sqrt(Xsvd$d^2 + regu), ncol = dd, nrow = dd)

  # RXregu used to NOT be an inverse as it is above. I think that was a bug,
  # but maybe I'm not thinking right and we have to switch it back?
  # Just a note in case things are behaving weird.
  UX %*% RXreguinv %*% RX %*% t(VX)
}

v_leverage <- function(C) {
  rowSums(svd(C)$v^2)
}

# this maybe has a tendency to underselect? should investigate
cca_ase <- function(A, traits, dhat = NULL, num_perms = 100) {
  # A : adjacency matrix n-by-n
  # traits : node features, n-by-p
  # B : number of MC iterates in permutation testing.
  # dhat : estimate of the network dimension.
  #		Defaults to using Sussman's estimate.
  #	max(c(1,which(svdA$d>=3^(1/4)*n^(3/4)*log(n)^(1/4))))
  #	see AK's my_ASE function.
  # fs : boolean for using Fisher-Sun or not.

  if (is.null(dhat)) {
    stop("Must specify `dhat`.", call. = FALSE)
  }

  s_A <- irlba(A, dhat)
  Uhat <- s_A$u
  Xhat <- s_A$u %*% diag(sqrt(s_A$d))

  n <- nrow(traits)
  if (nrow(Xhat) != n) {
    stop('Number of rows in data matrices should agree.')
  }

  Xpart <- construct_cca_part(Xhat, regu = 0)
  traitspart <- construct_cca_part(traits, regu = 0)

  cca_observed <- crossprod(Xpart, traitspart)
  stat_observed <- v_leverage(cca_observed)

  stat_permuted <- matrix(0, nrow = num_perms, ncol = length(stat_observed))

  for (i in 1:num_perms) {
    permutation <- sample(n, n, replace = FALSE)
    cca_permuted <- crossprod(Xpart, traitspart[permutation, ])

    # note comparison here
    stat_permuted[i, ] <- stat_observed < v_leverage(cca_permuted)
  }

  pvalues <- colMeans(stat_permuted)
  names(pvalues) <- colnames(traits)
  pvalues
}


cca_full <- function(
  A,
  traits,
  num_perms = 100,
  regu = NULL,
  alpha = 0.05,
  dhat = NULL
) {
  n <- nrow(traits)
  num_traits <- ncol(traits)

  if (nrow(A) != n) {
    stop('Number of rows in data matrices should agree.')
  }

  if (is.null(regu)) {
    regu <- sqrt(n)
  }

  Xpart <- construct_cca_part(A, regu = regu)
  traitspart <- construct_cca_part(traits, regu = 0)

  cca_observed <- crossprod(Xpart, traitspart)
  cca_svd <- svd(cca_observed)

  singular_observed <- cca_svd$d
  singular_permuted <- matrix(
    0,
    nrow = num_perms,
    ncol = length(singular_observed)
  )

  svd_permuted <- vector(mode = "list", length = length(num_perms))

  for (i in 1:num_perms) {
    permutation <- sample(n, n, replace = FALSE)
    cca_permuted <- crossprod(Xpart, traitspart[permutation, ])
    svd_permuted[[i]] <- svd(cca_permuted)

    # note comparison here
    singular_permuted[i, ] <- singular_observed < svd_permuted[[i]]$d
  }

  singular_pvalues <- colMeans(singular_permuted)

  if (all(singular_pvalues <= alpha)) {
    rank_estimate <- num_traits
  } else {
    rank_estimate <- min(which(singular_pvalues > alpha)) - 1
  }

  if (is.null(dhat)) {
    dhat <- rank_estimate
  }

  v_leverage_observed <- rowSums(
    cca_svd$v[, 1:min(dhat, num_traits), drop = FALSE]^2
  )
  v_leverage_permuted <- matrix(
    0,
    nrow = num_perms,
    ncol = length(v_leverage_observed)
  )

  for (i in 1:num_perms) {
    permuted_leverage <- rowSums(
      svd_permuted[[i]]$v[, 1:min(dhat, num_traits), drop = FALSE]^2
    )
    v_leverage_permuted[i, ] <- v_leverage_observed < permuted_leverage
  }

  trait_pvalues <- colMeans(v_leverage_permuted)
  names(trait_pvalues) <- colnames(traits)

  list(
    cca_spectrum_pvalues = singular_pvalues,
    rank_estimate = rank_estimate,
    trait_pvalues = trait_pvalues,
    regu = regu
  )
}


get_pilot_traits <- function(ard_data) {
  A <- ard_data$A
  traits <- ard_data$traits

  cf <- cca_full(A, traits)
  # ca <- cca_ase(A, traits, dhat = cf$rank_estimate)

  dhat <- max(2, cf$rank_estimate)

  gl <- grplasso_select(A, traits, dhat = dhat)

  marginal_cfs <- map_dfr(seq_len(ncol(traits)), \(j) {
    cca_full(A, traits[, j, drop = FALSE])
  })

  list(
    cf = cf,
    # ca = ca,
    marginal_cf = marginal_cfs,
    gl = gl
  )
}
