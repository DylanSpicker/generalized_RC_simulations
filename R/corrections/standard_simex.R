library(MASS)

simex <- function(proxies, 
                  Sigma, 
                  est.fun=mean, 
                  B=200, 
                  M=20, 
                  max.lambda=2, 
                  seed=NULL, 
                  models=list(Y~a+b*X), 
                  starts=list(list(a=1, b=1)), 
                  nls.control=list(warnOnly=TRUE),
                  which.parms = NULL) {
    if(!is.null(seed)) set.seed(seed)
    
    proxies <- lapply(proxies, as.matrix)
    eigen.d <- eigen(Sigma)
    neg_ev <- which(eigen.d$values < 0)
    q <- ncol(proxies[[1]])
    N <- nrow(as.matrix(proxies[[1]]))

    if(length(neg_ev) > 0 && length(neg_ev) < q) {
        # Compute the closest positive definite matrix to Sigma
        # if this is not possible, as all eigen values are zero and so it
        # would return the zero matrix, we simply use the naive estimator
        # the naive estimator, in that case, should be approximately equal 
        # to the truth.
        k_tr <- sum(eigen.d$values)
        eigen.d$values[neg_ev] <- 0
        eigen.d$values <- eigen.d$values*k_tr/sum(eigen.d$values)
        Sigma <- eigen.d$vectors %*% diag(eigen.d$values) %*% solve(eigen.d$vectors)
        warn("The provided covariance matrix is not positive definite. An approximate positive definite matrix is computed to approximate the provided covariance matrix.")
    } else if (length(neg_ev) > 0) {
        stop("The covariance matrix provided is not positive definite, and all eigenvalues are negative, meaning a positive definite matrix cannot be approximated.")
    }

    idx_mat <- do.call(cbind, lapply(proxies, function(Wj){ apply(! is.na(Wj), FUN=all, MARGIN=1) }))
    unique_idxs <- unique(idx_mat, MARGIN=1)
    lambdas <- seq(0, max.lambda, length.out=M)
    
    sim.dfs <- list()

    for (i_row in 1:nrow(unique_idxs)) {
        r_filter <- unique_idxs[i_row, ]
        ki <- sum(r_filter)
        y.filter <- which(apply(idx_mat, FUN=function(r){ all(r == r_filter) }, MARGIN=1))

        X <- as.matrix(Reduce("+", lapply(1:ki, function(j){
            (1/ki)*proxies[r_filter][[j]][y.filter, ]
        })))

        sim.df <- do.call(rbind, lapply(lambdas, function(l){
            if(l == 0){
                est.fun(X, y.filter=y.filter)
            } else {
                theta <- 0

                for(ii in 1:B){
                    rv <- mvrnorm(nrow(X), mu=rep(0, q), Sigma=l*(1/ki)*Sigma)
                    Xb <- rv + X
                    theta <- theta + est.fun(Xb, y.filter=y.filter)
                }
            
                (1/B)*theta
            }
        }))

        colnames(sim.df) <- 1:ncol(sim.df)
        sim.df <- data.frame(Y=sim.df, X=lambdas)
        sim.dfs[[length(sim.dfs)+1]] <- list(df = sim.df, 
                                             comb = r_filter, 
                                             n = length(y.filter))
    }
    ests <- lapply(sim.dfs, function(sim.df){
        (sim.df$n/N)*simex.ext(sim.df$df, models=models, starts=starts, nls.control=nls.control, which.parms = which.parms)
    })

    Reduce("+", ests)
}

empirical_simex <- function(proxies, Z=NULL,
                            est.fun = mean, 
                            B=200, 
                            M=20, 
                            max.lambda=2, 
                            seed=NULL, 
                            models=list(Y~a+b*X), 
                            starts=list(list(a=1, b=1)), 
                            nls.control=list(warnOnly=TRUE),
                            which.parms = NULL) {

    if(!is.null(seed)) set.seed(seed)

    proxies <- lapply(proxies, as.matrix)
    if(!is.null(Z)) parms <- get_correction_params(proxies, as.matrix(Z))
    else parms <- get_correction_params(proxies)
    
    idx_mat <- do.call(cbind, lapply(proxies, function(Wj){ apply(! is.na(Wj), FUN=all, MARGIN=1) }))

    ki <- parms$k.i
    Wbar.i <- parms$Wbar.i
    y.filter <- which(ki > 1)

    idx_mat <- idx_mat[y.filter, ]
    Wbar.i <- as.matrix(Wbar.i[y.filter, ])
    ki <- ki[y.filter]

    lambdas <- seq(0, max.lambda, length.out=M)
    total_obs <- sum(ki)
    sel_ends <- cumsum(ki)
    sel_starts <- c(1, sel_ends[1:(length(sel_ends)-1)]+1)

    contrast_means <- lapply(1:B, function(ii){
        pseudo_errors <- rnorm(total_obs)
        W.comb <- array(0, dim(Wbar.i))

        for (jj in 1:nrow(Wbar.i)) {
            Zi <- pseudo_errors[sel_starts[jj]:sel_ends[jj]]
            Zbar <- mean(Zi)

            contrast <- (Zi - Zbar)/sqrt(sum((Zi-Zbar)^2))
            W.comb[jj, ] <- Reduce("+", lapply(1:ki[jj], function(num){
                contrast[num]*proxies[idx_mat[jj,]][[num]][y.filter[jj], ]
            }))
        }

        W.comb
    })


    sim.df <- do.call(rbind, lapply(lambdas, function(l){
        if(l == 0){
            est.fun(Wbar.i, y.filter=y.filter)
        } else {
            (1/B)*Reduce("+", lapply(1:B, function(b){
                Xb <- Wbar.i + sqrt(l/ki)*contrast_means[[b]]
                est.fun(Xb, y.filter=y.filter)
            }))
        }
    }))

    sim.df <- data.frame(Y=sim.df, X=lambdas)
    simex.ext(sim.df, models=models, starts=starts, nls.control=nls.control, which.parms = which.parms)

}