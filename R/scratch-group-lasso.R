grplasso_select <- function(A, traits, dhat = NULL) {
    if (is.null(dhat)) {
        stop("Must specify `dhat`.", call. = FALSE)
    }

    if (dhat <= 1) {
        stop("`dhat` must be at least two, and an integer.", call. = FALSE)
    }

    s_A <- svds(A, dhat)
    Xhat <- s_A$u %*% diag(sqrt(s_A$d))

    gl_path <- cv.glmnet(
        x = traits,
        y = Xhat,
        family = "mgaussian",
        alpha = 1,
        intercept = FALSE # ?
    )

    # variables selected

    beta_1se <- coef(gl_path, s = "lambda.1se")
    beta_min <- coef(gl_path, s = "lambda.min")

    traits_1se <- names(which(beta_1se$y1[, 1] != 0))
    traits_min <- names(which(beta_min$y1[, 1] != 0))

    fit <- gl_path$glmnet.fit

    # assumes dhat >= 2, otherwise indexing into beta bad. i <3 lack of type stability
    betas <- as.matrix(fit$beta[[1]])
    varlist <- row.names(betas)
    which_step <- rep(NA, length(varlist))
    for (i in seq_along(varlist)) {
        which_step[i] <- which.max(betas[i, ] != 0)
    }
    entry_order <- varlist[order(which_step)]

    list(
        traits_1se = traits_1se,
        traits_min = traits_min,
        entry_order = entry_order
    )
}
