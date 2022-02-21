experiment_2_truth <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y
    est_coefs.true <- coef(glm(Y~Z+X, family=Gamma(link='log')))
    names(est_coefs.true) <- c("Intercept (Truth)", "Z (Truth)", "X (Truth)")
    data.frame(t(est_coefs.true), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_RC <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y

    Wstar.standard <- RC(W, Z)
    est_coefs.s <- coef(glm(Y~Z+Wstar.standard, family=Gamma(link='log')))
    names(est_coefs.s) <- c("Intercept (RC)", "Z (RC)", "X (RC)")
    data.frame(t(est_coefs.s), scenario = sim2.data$scenario[1], check.names=FALSE)

}

experiment_2_generalized_RC <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y
    Wstar.general <- generalized_RC(W, Z=as.matrix(Z), J0=c(1,2,3), J1=c(1,2,3))
    est_coefs.g <- coef(glm(Y~Z+Wstar.general, family=Gamma(link='log')))
    names(est_coefs.g) <- c("Intercept (Generalized RC)", "Z (Generalized RC)", "X (Generalized RC)")
    data.frame(t(est_coefs.g), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_generalized_RC_noZ <- function(sim2.data){
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y
    
    Wstar.general.noZ <- generalized_RC(W, Z=NULL, J0=c(1,2,3), J1=c(1,2,3))
    est_coefs.g.noZ <- coef(glm(Y~Z+Wstar.general.noZ, family=Gamma(link='log')))
    names(est_coefs.g.noZ) <- c("Intercept (Generalized RC (no Z))", "Z (Generalized RC (no Z))", "X (Generalized RC (no Z))")
    data.frame(t(est_coefs.g.noZ), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_generalized_RC_wtd <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y

    Wstar.general.wt <- generalized_RC.weights(W, Z=as.matrix(Z), J0=c(1,2,3), J1=c(1,2,3))
    est_coefs.g.wt <- coef(glm(Y~Z+Wstar.general.wt, family=Gamma(link='log')))
    names(est_coefs.g.wt) <- c("Intercept (Generalized RC Weighted)", "Z (Generalized RC Weighted)", "X (Generalized RC Weighted)")
    data.frame(t(est_coefs.g.wt), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_generalized_RC_wtd_noZ <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y

    Wstar.general.noZ.wt <- generalized_RC.weights(W, Z=NULL, J0=c(1,2,3), J1=c(1,2,3))
    est_coefs.g.noZ.wt <- coef(glm(Y~Z+Wstar.general.noZ.wt, family=Gamma(link='log')))
    names(est_coefs.g.noZ.wt) <- c("Intercept (Generalized RC Weighted (no Z))", "Z (Generalized RC Weighted (no Z))", "X (Generalized RC Weighted (no Z))")
    data.frame(t(est_coefs.g.noZ.wt), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_generalized_SIMEX <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y


    # SIMEX
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X + Z[y.filter], family=Gamma(link='log')))
        } else {
            coef(glm(Y~X + Z, family=Gamma(link='log')))
        }
    }

    sim.general <- generalized_simex(W,
                                    Z = as.matrix(Z),
                                    J0=c(1,2,3),
                                    J1=c(1,2,3),
                                    use.proxies=c(1,2,3), 
                                    est.fun=get_coef,
                                    M=10,
                                    B=50,
                                    models=list(Y~a+b/(c+X)), 
                                    starts=list(list(a=1,b=1,c=1)),
                                    save.object=FALSE,
                                    combine.first = FALSE)
    sim.ests <- unlist(sim.general$estimates)
    names(sim.ests) <- c("Intercept (Generalized SIMEX (A))", "X (Generalized SIMEX (A))", "Z (Generalized SIMEX (A))")
    data.frame(t(sim.ests), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_generalized_SIMEX_combine_first <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y


    # SIMEX
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X + Z[y.filter], family=Gamma(link='log')))
        } else {
            coef(glm(Y~X + Z, family=Gamma(link='log')))
        }
    }


    sim.general.c <- generalized_simex(W,
                                    Z = as.matrix(Z),
                                    J0=c(1,2,3),
                                    J1=c(1,2,3),
                                    use.proxies=c(1,2,3), 
                                    est.fun=get_coef,
                                    M=10,
                                    B=50,
                                    models=list(Y~a+b/(c+X)), 
                                    starts=list(list(a=1,b=1,c=1)),
                                    save.object=FALSE,
                                    combine.first = TRUE)

    sim.ests.c <- unlist(sim.general.c$estimates)
    names(sim.ests.c) <- c("Intercept (Generalized SIMEX (B))", "X (Generalized SIMEX (B))", "Z (Generalized SIMEX (B))")
    data.frame(t(sim.ests.c), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_standard_simex <- function(sim2.data) {
    Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y

    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X + Z[y.filter], family=Gamma(link='log')))
        } else {
            coef(glm(Y~X + Z, family=Gamma(link='log')))
        }
    }

    parms <- get_correction_params(W, Z)
    standard.simex <- simex(W,
                            parms$Sigma.UU,
                            est.fun=get_coef,
                            M=15,
                            models=list(Y~a+b/(c+X)), 
                            starts=list(list(a=1, b=1, c=1)))

    names(standard.simex) <- c("Intercept (SIMEX)", "X (SIMEX)", "Z (SIMEX)")
    data.frame(t(standard.simex), scenario = sim2.data$scenario[1], check.names=FALSE)
}

experiment_2_empirical_simex <- function(sim2.data) {
        Z <- sim2.data$Z
    W <- list(sim2.data$W1,
              sim2.data$W2,
              sim2.data$W3)
    X <- sim2.data$X
    Y <- sim2.data$Y

    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X + Z[y.filter], family=Gamma(link='log')))
        } else {
            coef(glm(Y~X + Z, family=Gamma(link='log')))
        }
    }
    
    empirical.simex.est <- empirical_simex(W, Z=Z,
                                            est.fun=get_coef,
                                            M=15,
                                            models=list(Y~a+b/(c+X)), 
                                            starts=list(list(a=1, b=1, c=1)))

    names(empirical.simex.est) <- c("Intercept (E-SIMEX)", "X (E-SIMEX)", "Z (E-SIMEX)")
    data.frame(t(empirical.simex.est), scenario = sim2.data$scenario[1], check.names=FALSE)
}