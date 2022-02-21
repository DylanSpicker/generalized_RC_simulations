library(rootSolve)

cov_total <- function(X, Y=NULL, meanX, meanY = NULL, scale = NULL) {
    X <- as.matrix(X)

    if(is.null(Y)) Y <- X
    else Y <- as.matrix(Y)

    if(is.null(meanY)) meanY = meanX
    if(is.null(scale)) scale <- rep(1, nrow(X))

    ret <- matrix(0, nrow = ncol(X), ncol = ncol(Y))

    for (ii in 1:nrow(X)) {
        if(!any(is.na(X[ii,])) && !any(is.na(Y[ii, ]))) {
            ret <- ret + scale[ii]*(X[ii, ] - meanX)%*%t(Y[ii, ] - meanY)
        }
    }

    ret
}


compute_correction_params <- function(obs) {
    if(length(obs$proxies) <= 1) stop("At least two proxy measurements are needed.") 
    if(length(obs$J1) <= 1) stop("At least 2 observations need to be in J1.")
    if (length(obs$J0) < 1) stop("At least 1 observation needs to be in J0.")

    use_Z <- !is.null(obs$error_free)

    obs$proxies <- lapply(obs$proxies, as.matrix)

    obs$J01 <- intersect(obs$J0, obs$J1)
    obs$J0_ <- setdiff(obs$J0, obs$J01)
    obs$J1_ <- setdiff(obs$J1, obs$J01)
    obs$K <- length(obs$proxies)
    obs$dimX <- ncol(obs$proxies[[1]])

    if(use_Z) {
        obs$error_free <- as.matrix(obs$error_free)
        obs$dimZ <- ncol(obs$error_free)
    }


    # Build Return Object
    return_obj <- list()


    # Collect the Building Blocks
    obs$covariance <- list()
    obs$mean <- list()
    obs$cross_terms.Xj.Z <- list() 

    for (ii in 1:obs$K) {
        obs$covariance[[ii]] <- list()
        obs$mean[[ii]] <- colMeans(obs$proxies[[ii]], na.rm=T)

        for (jj in 1:obs$K) {
            obs$covariance[[ii]][[jj]] <- cov(obs$proxies[[ii]], obs$proxies[[jj]], use="complete.obs")
        }
    }

    # Compute Parameters 

    if(use_Z) {
        return_obj$SigmaZXj <- list()

        return_obj$SigmaZ <- cov(obs$error_free)
        return_obj$MuZ <- colMeans(obs$error_free)

        # Return obj 
        return_obj$SigmaZXj <- lapply(obs$proxies, function(Wj){
            # cov(obs$error_free, Wj, use='complete.obs')
            1/(sum(complete.cases(cbind(obs$error_free, Wj)))-1)*cov_total(X = obs$error_free, Y = Wj, meanX = return_obj$MuZ, meanY = colMeans(Wj, na.rm=T))
        })

        return_obj$SigmaZX <- (1/obs$K)*Reduce("+", return_obj$SigmaZXj)

        return_obj$Eta1j <- list()
        for (jj in 1:obs$K) {
            if (jj %in% obs$J1) {
                return_obj$Eta1j[[jj]] <- rep(1, obs$dimX)
            } else {
                v1 <- cov(obs$error_free, obs$proxies[[jj]], use='complete.obs')
                return_obj$Eta1j[[jj]] <- as.numeric(sqrt(diag(solve(t(return_obj$SigmaZX)%*%return_obj$SigmaZX)%*%t(v1)%*%v1)))
            }
        }
        
        # SigmaXXj
        return_obj$SigmaXXj <- list()
        for(jj in 1:obs$K) {
            return_obj$SigmaXXj[[jj]] <- (1/(obs$K - 1))*Reduce("+", lapply(
                setdiff(1:obs$K, jj),
                function(idx){
                    diag(1/return_obj$Eta1j[[idx]], nrow=length(return_obj$Eta1j[[idx]]))%*%obs$covariance[[idx]][[jj]]
                }
            ))
        }

    } else {
        return_obj$MuZ <- NULL
        return_obj$SigmaZX <- NULL
 
        # SigmaXXj
        return_obj$SigmaXXj <- list()
        for (jj in 1:obs$K) {
            res <- 0
            idxs <- setdiff(obs$J1, jj)
            
            for(ii in idxs) {
                res <- res + obs$covariance[[ii]][[jj]]
            }

            return_obj$SigmaXXj[[jj]] <- res*(1/length(idxs))
        }

        # Eta1J
        return_obj$Eta1j <- list()
        for (jj in 1:obs$K){
            if (jj %in% obs$J1) {
                return_obj$Eta1j[[jj]] <- rep(1, obs$dimX)
            } else {
                res <- 0
                idxs <- setdiff(1:obs$K, jj)
                for (ii in idxs) {
                    res <- res + (1/length(idxs))*obs$covariance[[jj]][[ii]]%*%solve(return_obj$SigmaXXj[[ii]])
                }

                return_obj$Eta1j[[jj]] <- as.numeric(res)
            }
        }

    }

    # MuX
    total.X <- 0
    n.X <- 0
    total_len <- length(obs$J01) + length(obs$J0_)

    for (ii in 1:obs$K) {
        if (ii %in% obs$J01) {
            total.X <- total.X + colSums(obs$proxies[[ii]], na.rm=T)
            n.X <- n.X + sum(rowSums(is.na(obs$proxies[[ii]])) == 0)
        } else if (ii %in% obs$J0_) {
            total.X <- total.X + solve(diag(return_obj$Eta1j[[ii]], nrow=length(return_obj$Eta1j[[ii]])))%*%colSums(obs$proxies[[ii]], na.rm=T)
            n.X <- n.X + sum(rowSums(is.na(obs$proxies[[ii]])) == 0)
        }
    }
    
    return_obj$MuX <- total.X/n.X

    # Eta0J
    return_obj$Eta0j <- list()
    for (jj in 1:obs$K) {
        if (jj %in% obs$J0){
            return_obj$Eta0j[[jj]] <- rep(0, obs$dimX)
        } else {
            return_obj$Eta0j[[jj]] <- obs$mean[[jj]] - diag(return_obj$Eta1j[[jj]], nrow=length(return_obj$Eta1j[[jj]]))%*%return_obj$MuX
        }
    }

    # SigmaX 
    return_obj$SigmaX <- 0
    for (ii in 1:obs$K) {
        return_obj$SigmaX <- return_obj$SigmaX + (1/obs$K)*solve(diag(return_obj$Eta1j[[ii]], nrow=length(return_obj$Eta1j[[ii]])))%*%return_obj$SigmaXXj[[ii]]
    }

    # Mj 
    return_obj$Mj <- list()
    for (jj in 1:obs$K){
        return_obj$Mj[[jj]] <- obs$covariance[[jj]][[jj]] - diag(return_obj$Eta1j[[jj]], nrow=length(return_obj$Eta1j[[jj]]))%*%return_obj$SigmaX%*%diag(return_obj$Eta1j[[jj]], nrow=length(return_obj$Eta1j[[jj]]))
    }

    # MuJ
    return_obj$Muj <- obs$mean

    return(return_obj)
}

generalized_RC.weights <- function(proxies, Z=NULL, J0=NULL, J1=NULL,return_parms=FALSE) {
    if(length(proxies) <= 1) stop("At least two proxy measurements are needed.")
    if(is.null(J0)) J0 <- 1:length(proxies)
    if(is.null(J1)) J1 <- 1:length(proxies)

    if(length(J1) <= 1) stop("At least 2 observations need to be in J1.")
    if (length(J0) < 1) stop("At least 1 observation needs to be in J0.")

    proxies <- lapply(proxies, as.matrix)
    

    if(!is.null(Z)){ 
        Z <- as.matrix(Z)
        corr <- compute_correction_params(list(proxies=proxies, error_free=Z, J0=J0, J1=J1))
    } else {
        corr <- compute_correction_params(list(proxies=proxies, error_free=NULL, J0=J0, J1=J1))
    }

    mat <- array(0, dim(proxies[[1]]))

    idx_mat <- do.call(cbind, lapply(proxies, function(Wj){ apply(! is.na(Wj), FUN=all, MARGIN=1) }))
    unique_idxs <- unique(idx_mat, MARGIN=1)

    if(return_parms) {
        parm_list <- list()
    }
    
    for (i_row in 1:nrow(unique_idxs)) {
        r_filter <- unique_idxs[i_row, ]
        in_set <- which(apply(idx_mat, FUN=function(r){ all(r == r_filter) }, MARGIN=1))

        if(is.null(Z)) {
            start.SigmaX <- (1/sum(r_filter))*Reduce("+", corr$SigmaXXj[r_filter])
            start.M <- (1/sum(r_filter))^2*Reduce("+", corr$Mj[r_filter])
            start.beta <- start.SigmaX%*%solve(start.SigmaX + start.M)
            start.mu.star <- (1/sum(r_filter))*Reduce("+", corr$Muj[r_filter])

            starts <- list(
                mu = as.numeric(corr$MuX - start.beta%*%start.mu.star),
                beta = as.numeric(start.beta),
                alpha = rep(1/sum(r_filter), sum(r_filter)),
                lambda = 0
            )
            
            names <- list(
                mu = 1:length(starts$mu),
                beta = (length(starts$mu)+1):(length(starts$mu)+length(starts$beta)),
                alpha = (length(starts$mu) + length(starts$beta)+1):(length(starts$mu) + length(starts$beta)+length(starts$alpha)),
                lambda = length(starts$mu) + length(starts$beta) + length(starts$alpha) + 1
            )

            g <- function(theta, names) {
                
                p.mu <- matrix(theta[names$mu])
                p.beta <- matrix(theta[names$beta], nrow=sqrt(length(theta[names$beta])))
                p.lambda <- theta[names$lambda]

                mu_X.star <- Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]*corr$Muj[r_filter][[ii]]
                }))

                Sigma.XX.star <- Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]*corr$SigmaXXj[r_filter][[ii]]
                }))

                Sigma.X.star <- corr$SigmaX + Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]^2*corr$Mj[r_filter][[ii]]
                }))


                component_1 <- corr$MuX - p.mu - p.beta%*%mu_X.star
                component_2 <- Sigma.XX.star - p.beta%*%Sigma.X.star

                components_3 <- lapply(1:sum(r_filter), function(ii){
                    Sigma.XsXj <- corr$SigmaX + theta[names$alpha[ii]]*corr$Mj[r_filter][[ii]]
                    sum(diag((corr$SigmaXXj[r_filter][[ii]] - p.beta%*%Sigma.XsXj)%*%t(p.beta))) - p.lambda
                })
                component_4 <- 1 - sum(theta[names$alpha])

                c(component_1, component_2, unlist(components_3), component_4)
            }

        } else {
            start.SigmaX <- (1/sum(r_filter))*Reduce("+", corr$SigmaXXj[r_filter])
            start.M <- (1/sum(r_filter))^2*Reduce("+", corr$Mj[r_filter])
            start.SigmaZXs <- (1/sum(r_filter))*Reduce("+", corr$SigmaZXj[r_filter])
            start.inv.SigmaZ <- solve(corr$SigmaZ)

            start.beta <- (start.SigmaX - t(corr$SigmaZX)%*%start.inv.SigmaZ%*%start.SigmaZXs)%*%solve(start.SigmaX + start.M - t(start.SigmaZXs)%*%start.inv.SigmaZ%*%start.SigmaZXs)
            start.gamma <- (t(corr$SigmaZX) - start.beta%*%t(start.SigmaZXs))%*%start.inv.SigmaZ
            start.mu.star <- (1/sum(r_filter))*Reduce("+", corr$Muj[r_filter])

            starts <- list(
                mu = as.numeric(corr$MuX - start.beta%*%start.mu.star - start.gamma%*%corr$MuZ),
                beta = as.numeric(start.beta),
                gamma = as.numeric(start.gamma),
                alpha = rep(1/sum(r_filter), sum(r_filter)),
                lambda = 0
            )
            
            names <- list(
                mu = 1:length(starts$mu),
                beta = (length(starts$mu)+1):(length(starts$mu)+length(starts$beta)),
                gamma = (length(starts$mu) + length(starts$beta)+1):(length(starts$mu) + length(starts$beta)+length(starts$gamma)),
                alpha = (length(starts$mu) + length(starts$beta)+length(starts$gamma)+1):(length(starts$mu) + length(starts$beta)+length(starts$alpha)+length(starts$gamma)),
                lambda = length(starts$mu) + length(starts$beta)+length(starts$alpha)+length(starts$gamma) + 1
            )


            g <- function(theta, names) {

                p.mu <- matrix(theta[names$mu])
                p.beta <- matrix(theta[names$beta], nrow=sqrt(length(theta[names$beta])))
                p.gamma <- matrix(theta[names$gamma], nrow=sqrt(length(theta[names$gamma])))
                p.lambda <- theta[names$lambda]

                mu_X.star <- Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]*corr$Muj[r_filter][[ii]]
                }))

                Sigma.XX.star <- Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]*corr$SigmaXXj[r_filter][[ii]]
                }))

                Sigma.X.star <- corr$SigmaX + Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]^2*corr$Mj[r_filter][[ii]]
                }))

                Sigma.ZX.star <-  Reduce("+", lapply(1:sum(r_filter), function(ii){
                    theta[names$alpha[ii]]*corr$SigmaZXj[r_filter][[ii]]
                }))

                component_1 <- corr$MuX - p.mu - p.beta%*%mu_X.star - p.gamma%*%corr$MuZ
                component_2 <- Sigma.XX.star - p.mu%*%mu_X.star - p.beta%*%Sigma.X.star - p.gamma%*%Sigma.ZX.star
                component_3 <- t(corr$SigmaZX) - p.mu%*%corr$MuZ - p.beta%*%t(Sigma.ZX.star) - p.gamma%*%corr$SigmaZ
                components_4 <- lapply(1:sum(r_filter), function(ii){
                    Sigma.XsXj <- corr$SigmaX + theta[names$alpha[ii]]*corr$Mj[r_filter][[ii]]
                    sum(diag((corr$SigmaXXj[r_filter][[ii]] - p.mu%*%t(corr$Muj[r_filter][[ii]]) - p.beta%*%Sigma.XsXj - p.gamma%*%corr$SigmaZXj[r_filter][[ii]])%*%p.beta)) - p.lambda
                })
                component_5 <- 1 - sum(theta[names$alpha])

                c(component_1, component_2, component_3, unlist(components_4), component_5)
            }
        }

        solve_parms <- tryCatch({
            solve_parms <- multiroot(f=g, start=unlist(starts), parms=names)
            solve_parms$root
        }, warning = function(w){
            unlist(starts)
        })

        mu <- matrix(solve_parms[names$mu])
        beta <- matrix(solve_parms[names$beta], nrow=sqrt(length(names$beta)))
        if(!is.null(Z)) gamma <- matrix(solve_parms[names$gamma], nrow=sqrt(length(names$gamma)))
        alpha <- solve_parms[names$alpha]

        if(return_parms) {
            parm_list = list(mu=mu, beta=beta, alpha=alpha)
        }

        Xs <- Reduce("+", lapply(1:length(alpha), function(j){
            alpha[j]*proxies[r_filter][[j]]
        }))

        for(ii in in_set) {
            
            if(!is.null(Z)) mat[ii, ] <- mu + beta%*%Xs[ii,] + gamma%*%Z[ii,]
            else mat[ii, ] <- mu + beta%*%Xs[ii,]
        }

    }

    if(return_parms) {
        list(mat = mat, parms = parm_list)
    } else {
        mat
    }

}

generalized_RC <- function(proxies, Z=NULL, J0=NULL, J1=NULL) {
    if(length(proxies) <= 1) stop("At least two proxy measurements are needed.")
    if(is.null(J0)) J0 <- 1:length(proxies)
    if(is.null(J1)) J1 <- 1:length(proxies)

    if(length(J1) <= 1) stop("At least 2 observations need to be in J1.")
    if (length(J0) < 1) stop("At least 1 observation needs to be in J0.")

    proxies <- lapply(proxies, as.matrix)
    

    if(!is.null(Z)){ 
        Z <- as.matrix(Z)
        corr <- compute_correction_params(list(proxies=proxies, error_free=Z, J0=J0, J1=J1))
    } else {
        corr <- compute_correction_params(list(proxies=proxies, error_free=NULL, J0=J0, J1=J1))
    }

    mat <- array(0, dim(proxies[[1]]))

    for(ii in 1:nrow(proxies[[1]])) {
        
        idx <- which(unlist(lapply(proxies, function(Wj) { any(is.na(Wj[ii, ])) })) == 0)
        
        ki <- length(idx)
        Sigma.XXs <- (1/ki)*Reduce("+", corr$SigmaXXj[idx])
        Mu.Xs <- (1/ki)*Reduce("+", corr$Muj[idx])

        Sigma.X <- (1/ki)^2*Reduce("+", lapply(idx, function(ii){
            Reduce("+", lapply(idx, function(jj){
                diag(corr$Eta1j[[ii]], nrow=length(corr$Eta1j[[ii]]))%*%corr$SigmaX%*%diag(corr$Eta1j[[jj]], nrow=length(corr$Eta1j[[jj]]))
            }))
        }))

        Sigma.U <- (1/ki)^2*Reduce("+", corr$Mj[idx])
        Xs <- (1/ki)*Reduce("+", lapply(proxies[idx], function(Wj){ Wj[ii,] }))

        if(is.null(Z)){ 
            mat[ii, ] <- t(corr$MuX + Sigma.XXs%*%solve(Sigma.X + Sigma.U)%*%(Xs - Mu.Xs))
        } else{ 
            SigmaZXs <- (1/ki)*Reduce("+", corr$SigmaZXj[idx])

            mat[ii, ] <- t(corr$MuX + cbind(Sigma.XXs, t(corr$SigmaZX)) %*% 
                                     solve(rbind(
                                         cbind(Sigma.X + Sigma.U, t(SigmaZXs)),
                                         cbind(SigmaZXs, corr$SigmaZ))) %*%
                                    as.matrix(c(Xs - Mu.Xs, Z[ii,] - corr$MuZ)))
        }
    }

    mat

}

