library(MASS)

# Start by assuming that X is scalar
generalized_simex.gen.combine_first <- function(proxies, J0=NULL, J1=NULL, use.proxies=NULL, Z=NULL, est.fun=mean, B=200, M=5, max.lambda=2, seed=NULL, ...) {
    
    if(length(proxies) < 2) stop("At least two proxies are required to use the generalized method.")
    proxies <- lapply(proxies, as.matrix)

    if(is.null(J0)) J0 <- 1:length(proxies)
    if(is.null(J1)) J1 <- 1:length(proxies)
    if(is.null(use.proxies)) use.proxies <- 1:length(proxies)

    if(!is.null(seed)) set.seed(seed)
    
    if(!is.null(Z)) {
        parms <- compute_correction_params(list(proxies = proxies, error_free=as.matrix(Z), J0=J0, J1=J1))
    } else {
        parms <- compute_correction_params(list(proxies = proxies, J0=J0, J1=J1))
    }


    lambdas <- seq(0, max.lambda, length.out=M)
    q <- ncol(proxies[[1]])
    
    idx_mat <- do.call(cbind, lapply(proxies, function(Wj){ apply(! is.na(Wj), FUN=all, MARGIN=1) }))
    unique_idxs <- unique(idx_mat, MARGIN=1)
    
    data_frames <- list()

    for (i_row in 1:nrow(unique_idxs)) {
        r_filter <- unique_idxs[i_row, ]
        ki <- sum(r_filter)
        y.filter <- which(apply(idx_mat, FUN=function(r){ all(r == r_filter) }, MARGIN=1))

        #TODO: Update this to take weights, once we decide how they can be computed.
        X <- as.matrix(Reduce("+", lapply(1:ki, function(j){
            (1/ki)*proxies[r_filter][[j]][y.filter, ]
        })))

        eta0 <- Reduce("+", lapply(1:ki, function(j){
            (1/ki)*parms$Eta0j[r_filter][[j]]
        }))
        eta1 <- Reduce("+", lapply(1:ki, function(j){
            (1/ki)*parms$Eta1j[r_filter][[j]]
        }))
        Mj <- Reduce("+", lapply(1:ki, function(j){ 
            (1/ki)^2*diag(parms$Eta1j[r_filter][[j]], nrow=length(parms$Eta1j[r_filter][[j]]))%*%parms$Mj[r_filter][[j]]%*%diag(parms$Eta1j[r_filter][[j]], nrow=length(parms$Eta1j[r_filter][[j]]))
        }))

        thetas <- list()

        for (l in lambdas){

            theta <- 0
            Sigma <- l*Mj
            continue <- FALSE 

            if(l != 0) {
                eigen.d <- eigen(Sigma)
                neg_ev <- which(eigen.d$values < 0)

                if(length(neg_ev) > 0 && length(neg_ev) < q) {
                    # Compute the closest positive definite matrix to Sigma
                    # if this is not possible, as all eigen values are zero and so it
                    # would return the zero matrix, we simply use the naive estimator
                    # the naive estimator, in that case, should be approximately equal 
                    # to the truth.
                    continue <- TRUE 
                    k_tr <- sum(eigen.d$values)
                    eigen.d$values[neg_ev] <- 0
                    eigen.d$values <- eigen.d$values*k_tr/sum(eigen.d$values)
                    Sigma <- eigen.d$vectors %*% diag(eigen.d$values) %*% solve(eigen.d$vectors)
                } else if (length(neg_ev) == 0) {
                    continue <- TRUE 
                }
            }
            
            if (continue) {
                for(ii in 1:B){
                    rv <- mvrnorm(nrow(X), mu=rep(0, q), Sigma=Sigma)
                    Xb <- 1/as.numeric(eta1)*(rv + X - as.numeric(eta0))
                    theta <- theta + est.fun(Xb, y.filter=y.filter)
                }
                
                theta <- (1/B)*theta

            } else {

                Xb <- 1/as.numeric(eta1)*(X - as.numeric(eta0))

                theta <- est.fun(Xb, y.filter=y.filter)
            }

            thetas[[length(thetas)+1]] <- theta
        }

        thetas <- do.call(rbind, thetas)
        colnames(thetas) <- 1:ncol(thetas)
        data_frames[[length(data_frames) + 1]] <- list(df = data.frame(Y = thetas, X = lambdas), 
                                                       comb = r_filter, 
                                                       n = length(y.filter))

    }
    data_frames
}

# Start by assuming that X is scalar
generalized_simex.gen <- function(proxies, J0=NULL, J1=NULL, use.proxies=NULL, Z=NULL, est.fun=mean, B=200, M=5, max.lambda=2, seed=NULL, ...) {
    
    if(length(proxies) < 2) stop("At least two proxies are required to use the generalized method.")
    proxies <- lapply(proxies, as.matrix)

    if(is.null(J0)) J0 <- 1:length(proxies)
    if(is.null(J1)) J1 <- 1:length(proxies)
    if(is.null(use.proxies)) use.proxies <- 1:length(proxies)

    if(!is.null(seed)) set.seed(seed)
    
    if(!is.null(Z)) {
        parms <- compute_correction_params(list(proxies = proxies, error_free=as.matrix(Z), J0=J0, J1=J1))
    } else {
        parms <- compute_correction_params(list(proxies = proxies, J0=J0, J1=J1))
    }


    lambdas <- seq(0, max.lambda, length.out=M)
    n <- nrow(proxies[[1]])
    q <- ncol(proxies[[1]])

    data_frames <- list()

    for(idx in use.proxies){
        X <- proxies[[idx]]
        idxs <- complete.cases(X)

        if(sum(idxs) == nrow(X)) {
            y.filter <- NULL
        } else {
            X <- as.matrix(X[idxs, ])
            y.filter <- which(idxs)
        }

        Mj <- parms$Mj[[idx]]
        eta0 <- parms$Eta0j[[idx]]
        eta1 <- parms$Eta1j[[idx]]
        
        thetas <- list()

        for (l in lambdas){

            theta <- 0
            Sigma <- l*Mj
            continue <- FALSE 

            if(l != 0) {
                eigen.d <- eigen(Sigma)
                neg_ev <- which(eigen.d$values < 0)

                if(length(neg_ev) > 0 && length(neg_ev) < q) {
                    # Compute the closest positive definite matrix to Sigma
                    # if this is not possible, as all eigen values are zero and so it
                    # would return the zero matrix, we simply use the naive estimator
                    # the naive estimator, in that case, should be approximately equal 
                    # to the truth.
                    continue <- TRUE 
                    k_tr <- sum(eigen.d$values)
                    eigen.d$values[neg_ev] <- 0
                    eigen.d$values <- eigen.d$values*k_tr/sum(eigen.d$values)
                    Sigma <- eigen.d$vectors %*% diag(eigen.d$values) %*% solve(eigen.d$vectors)
                } else if (length(neg_ev) == 0) {
                    continue <- TRUE 
                }
            }
            
            if (continue) {
                ii <- 1
                while(ii <= B){
                    rv <- mvrnorm(nrow(X), mu=rep(0, q), Sigma=Sigma)
                    Xb <- 1/as.numeric(eta1)*(rv + X - as.numeric(eta0))
                    new_cont <- est.fun(Xb, y.filter=y.filter)
                    ii <- ii + 1
                    theta <- theta + new_cont
                }
                
                theta <- (1/B)*theta

            } else {

                Xb <- 1/as.numeric(eta1)*(X - as.numeric(eta0))

                theta <- est.fun(Xb, y.filter=y.filter)
            }

            thetas[[length(thetas)+1]] <- theta
        }

        thetas <- do.call(rbind, thetas)
        colnames(thetas) <- 1:ncol(thetas)
        data_frames[[length(data_frames) + 1]] <- data.frame(Y = thetas, X = lambdas)
    }

    data_frames
}

simex.ext <- function(data, 
                      models=list(Y~a + b*X), 
                      starts = list(list(a=1,b=1)),
                      which.parms = NULL,
                      nls.control=list(warnOnly=TRUE),
                      ...) {
    
    

    unlist(lapply(1:length(models), function(ii) {
        p <- c()

        if(length(models) == 1){
            for(param in 1:(ncol(data)-1)) {
                dsub <- data[,c(param,ncol(data))]
                names(dsub) <- c("Y","X")
                
                if (class(starts[[ii]]) == "list") {
                    start.vals <- starts[[ii]]
                } else if (class(starts[[ii]]) == "function") {
                    start.vals <- starts[[ii]](dsub)
                }

                p <- c(p, predict(suppressWarnings(nls(models[[ii]], 
                                            data=dsub, 
                                            start=start.vals,
                                            control=nls.control,
                                            ...)), newdata=data.frame(X=c(-1))))
            }
        } else if (class(which.parms) == 'list' && length(which.parms) == length(models)) {
            for(param in which.parms[[ii]]) {
                dsub <- data[,c(param,ncol(data))]
                names(dsub) <- c("Y","X")
                
                if (class(starts[[ii]]) == "list") {
                    start.vals <- starts[[ii]]
                } else if (class(starts[[ii]]) == "function") {
                    start.vals <- starts[[ii]](dsub)
                }

                p <- c(p, predict(suppressWarnings(nls(models[[ii]], 
                                            data=dsub, 
                                            start=start.vals,
                                            control=nls.control,
                                            ...)), newdata=data.frame(X=c(-1))))
            }        
        } else {
            dsub <- data[,c(ii,ncol(data))]
            names(dsub) <- c("Y","X")
            
            if (class(starts[[ii]]) == "list") {
                start.vals <- starts[[ii]]
            } else if (class(starts[[ii]]) == "function") {
                start.vals <- starts[[ii]](dsub)
            }

            p <- predict(suppressWarnings(nls(models[[ii]], 
                                              data=dsub, 
                                              start=start.vals,
                                              control=nls.control,
                                              ...)), newdata=data.frame(X=c(-1)))
        } 
        
        p
    }))
}

generalized_simex <- function(proxies, 
                              J0=NULL, 
                              J1=NULL, 
                              use.proxies=NULL, 
                              Z=NULL, 
                              est.fun=mean, 
                              B=200, 
                              M=5, 
                              max.lambda=2,
                              seed=NULL, 
                              save.object = FALSE,
                              models=list(Y~a + b/(c+X)), 
                              starts = list(c(a=1,b=1,c=1)), 
                              nls.control=list(warnOnly=TRUE), 
                              combine.first = TRUE,
                              which.parms = NULL,
                              ...) {

    
    if(combine.first) {
        sim.df.list <- generalized_simex.gen.combine_first(proxies, J0, J1, use.proxies, Z, est.fun, B=B, M=M, max.lambda=max.lambda,seed=seed, ...)
        N <- nrow(as.matrix(proxies[[1]]))
        ests <- lapply(sim.df.list, function(sim.df){
            (sim.df$n/N)*simex.ext(sim.df$df, models=models, starts=starts, nls.control=nls.control, which.parms = which.parms, ...)
        })

    } else {
        sim.df.list <- generalized_simex.gen(proxies, J0, J1, use.proxies, Z, est.fun, B=B, M=M, max.lambda=max.lambda,seed=seed, ...)        
        ests <- lapply(sim.df.list, function(sim.df){
            (1/length(sim.df.list))*simex.ext(sim.df, models=models, starts=starts, nls.control=nls.control, which.parms = which.parms,  ...)
        })
    }

    
    ret_ests <- Reduce("+", ests)

    if(save.object) {
        list(estimates = ret_ests, sim.df = sim.df.list)
    } else {
        list(estimates = ret_ests)
    }
}
