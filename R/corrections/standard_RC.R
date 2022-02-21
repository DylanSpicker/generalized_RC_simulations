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

get_correction_params <- function(W, Z=NULL) {
    W <- lapply(W, as.matrix)
    if(!is.null(Z)) Z <- as.matrix(Z)

    get_k <- function(W) {
        length(W) - Reduce("+", lapply(W,
            function(Wj) {
                rowSums(is.na(Wj)) > 0
            }
        ))
    }
    get_inner_sum <- function(W){
        Reduce("+", lapply(W, function(Wj){
            Wj[!complete.cases(Wj), ] <- 0
            Wj
        }))
    }
    get_MuX <- function(k.i, T.i) {
        colSums(T.i)/sum(k.i)
    }
    get_v <- function(k.i){
        k.i <- get_k(W)
        sum(k.i) - sum(k.i^2)/sum(k.i)
    }
    get_sigma_uu <- function(Wbar.i, W, k.i) {
        Reduce("+", lapply(W, function(Wj){
            mat <- 0
            for (ii in 1:nrow(Wj)) {
                if(!any(is.na(Wj[ii,]))) mat <- mat + (Wj[ii,] - Wbar.i[ii,])%*%t(Wj[ii,] - Wbar.i[ii,])
            }
            mat 
        }))/sum(k.i-1)
    }
    
    k.i <- get_k(W)
    T.i <- get_inner_sum(W)
    v <- get_v(W) 

    Wbar.i <- T.i/k.i 
    mu.X <- get_MuX(k.i, T.i)
    mu.Z <- NULL
    Sigma.UU <- get_sigma_uu(Wbar.i, W, k.i)
    Sigma.ZZ <- NULL
    Sigma.XZ <- NULL

    if(!is.null(Z)) {
        Sigma.ZZ <- cov(Z)
        mu.Z <- colMeans(Z)
        Sigma.XZ <- cov_total(Wbar.i, Z, mu.X, mu.Z, scale = k.i)/v
    }

    Sigma.XX <- (cov_total(Wbar.i, meanX=mu.X, scale = k.i) - (nrow(Wbar.i) - 1)*Sigma.UU)/v

    list(
        n=nrow(W[[1]]),
        k.i = k.i,
        Wbar.i = Wbar.i,
        mu.X = mu.X,
        mu.W = mu.X,
        mu.Z = mu.Z,
        Sigma.UU = Sigma.UU,
        Sigma.XX = Sigma.XX,
        Sigma.ZZ = Sigma.ZZ,
        Sigma.XZ = Sigma.XZ
    )
}

RC <- function(W, Z=NULL) {
    params <- get_correction_params(W, Z)
    corrected <- array(0, dim(as.matrix(W[[1]])))

    if (is.null(Z)){
        for(ii in 1:params$n){
            corrected[ii, ] <- t(params$mu.W + 
            params$Sigma.XX%*%solve(params$Sigma.XX + params$Sigma.UU/params$k.i[ii])%*%(params$Wbar.i[ii,] - params$mu.W))
        }
    } else {
        Z <- as.matrix(Z)
        for(ii in 1:params$n){
            corrected[ii, ] <- params$mu.W + 
            cbind(params$Sigma.XX,
                  params$Sigma.XZ)%*%solve(
                      rbind(cbind(params$Sigma.XX + params$Sigma.UU/params$k.i[ii],
                            params$Sigma.XZ),
                      cbind(t(params$Sigma.XZ), params$Sigma.ZZ)))%*%as.matrix(c(params$Wbar.i[ii,] - params$mu.W, Z[ii,] - params$mu.Z))
        }
    }

    corrected
}
