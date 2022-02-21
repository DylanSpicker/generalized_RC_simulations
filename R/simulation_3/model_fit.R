experiment_3_all <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.all <- list(sim3.data$W1, sim3.data$W2, sim3.data$W3)

    Wstar.all <- generalized_RC(W.all, Z=NULL, J0=c(1, 2), J1=c(1,2))
    fit.mod.all <- glm(Y~Wstar.all, family=binomial(link='logit'))
    est_coefs.all <- coef(fit.mod.all)
    names(est_coefs.all) <- c("Intercept (Generalized RC)", "X (Generalized RC)")
    
    data.frame(t(est_coefs.all), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_reps <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.reps <- list(sim3.data$W1, sim3.data$W2)
    
    Wstar.reps <- generalized_RC(W.reps, Z=NULL, J0=c(1, 2), J1=c(1,2))
    fit.mod.reps <- glm(Y~Wstar.reps, family=binomial(link='logit'))
    est_coefs.reps <- coef(fit.mod.reps)
    names(est_coefs.reps) <- c("Intercept (Generalized RC Reps. only)", "X (Generalized RC Reps. only)")

    data.frame(t(est_coefs.reps), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_all_wt <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.all <- list(sim3.data$W1, sim3.data$W2, sim3.data$W3)

    Wstar.all.wt <- generalized_RC.weights(W.all, Z=NULL, J0=c(1, 2), J1=c(1,2))
    fit.mod.all.wt <- glm(Y~Wstar.all.wt, family=binomial(link='logit'))
    est_coefs.all.wt <- coef(fit.mod.all.wt)
    names(est_coefs.all.wt) <- c("Intercept (Generalized RC Weighted)", "X (Generalized RC Weighted)")

    data.frame(t(est_coefs.all.wt), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_reps_wt <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.reps <- list(sim3.data$W1, sim3.data$W2)
    Wstar.reps.wt <- generalized_RC.weights(W.reps, Z=NULL, J0=c(1, 2), J1=c(1,2))
    fit.mod.reps.wt <- glm(Y~Wstar.reps.wt, family=binomial(link='logit'))
    est_coefs.reps.wt <- coef(fit.mod.reps.wt)
    names(est_coefs.reps.wt) <- c("Intercept (Generalized RC Weighted Reps. only)", "X (Generalized RC Weighted Reps. only)")


    data.frame(t(est_coefs.reps.wt), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_generalized_sim <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.all <- list(sim3.data$W1, sim3.data$W2, sim3.data$W3)
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X, family=binomial(link='logit')))
            # coef(lm(Y[y.filter] ~ X))
        } else {
            coef(glm(Y~X, family=binomial(link='logit')))
            # coef(lm(Y~X))
        }
    }
    sim.general <- generalized_simex(W.all, 
                                     Z = NULL, 
                                     J0=c(1, 2), 
                                     J1=c(1, 2), 
                                     use.proxies=c(1,2,3), 
                                     est.fun=get_coef,
                                     M=10, 
                                     B=50, 
                                     models=list(Y~a+b/(c+X)), 
                                     starts=list(list(a=1, b=1, c=1)),  
                                     save.object=FALSE,
                                     combine.first = FALSE)
    sim.ests <- unlist(sim.general$estimates)
    names(sim.ests) <- c("Intercept (Generalized SIMEX (A))", "X (Generalized SIMEX (A))")

    data.frame(t(sim.ests), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_generalized_sim_reps <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.reps <- list(sim3.data$W1, sim3.data$W2)

    # SIMEX
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X, family=binomial(link='logit')))
            # coef(lm(Y[y.filter] ~ X))
        } else {
            coef(glm(Y~X, family=binomial(link='logit')))
            # coef(lm(Y~X))
        }
    }

    sim.reps <- generalized_simex(W.reps, 
                                  Z = NULL, 
                                  J0=c(1,2), 
                                  J1=c(1,2), 
                                  use.proxies=c(1,2), 
                                  est.fun=get_coef, 
                                  M=10, 
                                  B=50, 
                                  models=list(Y~a+b/(c+X)), 
                                  starts=list(list(a=1, b=1, c=1)),  
                                  save.object=FALSE, 
                                  combine.first = FALSE)
    sim.reps <- unlist(sim.reps$estimates)
    names(sim.reps) <- c("Intercept (Generalized SIMEX (A) Reps. only)", "X (Generalized SIMEX (A) Reps. only)")

    data.frame(t(sim.reps), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_generalized_sim_combine_first <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.all <- list(sim3.data$W1, sim3.data$W2, sim3.data$W3)

    # SIMEX
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X, family=binomial(link='logit')))
            # coef(lm(Y[y.filter] ~ X))
        } else {
            coef(glm(Y~X, family=binomial(link='logit')))
            # coef(lm(Y~X))
        }
    }

    sim.general.c <- generalized_simex(W.all, 
                                     Z = NULL, 
                                     J0=c(1, 2), 
                                     J1=c(1, 2), 
                                     use.proxies=c(1,2,3), 
                                     est.fun=get_coef,
                                     M=10, 
                                     B=50, 
                                     models=list(Y~a+b/(c+X)), 
                                     starts=list(list(a=1, b=1, c=1)),  
                                     save.object=FALSE,
                                     combine.first = TRUE)

    sim.ests.c <- unlist(sim.general.c$estimates)
    names(sim.ests.c) <- c("Intercept (Generalized SIMEX (B))", "X (Generalized SIMEX (B))")

    data.frame(t(sim.ests.c), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_generalized_sim_combine_first_reps <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.reps <- list(sim3.data$W1, sim3.data$W2)
    
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X, family=binomial(link='logit')))
            # coef(lm(Y[y.filter] ~ X))
        } else {
            coef(glm(Y~X, family=binomial(link='logit')))
            # coef(lm(Y~X))
        }
    }
    sim.reps.c <- generalized_simex(W.reps, 
                                  Z = NULL, 
                                  J0=c(1,2), 
                                  J1=c(1,2), 
                                  use.proxies=c(1,2), 
                                  est.fun=get_coef, 
                                  M=10, 
                                  B=50, 
                                  models=list(Y~a+b/(c+X)), 
                                  starts=list(list(a=1, b=1, c=1)),  
                                  save.object=FALSE, 
                                  combine.first = TRUE)
    sim.reps.c <- unlist(sim.reps.c$estimates)
    names(sim.reps.c) <- c("Intercept (Generalized SIMEX (B) Reps. only)", "X (Generalized SIMEX (B) Reps. only)")

    data.frame(t(sim.reps.c), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_simex <- function(sim3.data) {
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.all <- list(sim3.data$W1, sim3.data$W2, sim3.data$W3)

    # SIMEX
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X, family=binomial(link='logit')))
        } else {
            coef(glm(Y~X, family=binomial(link='logit')))
        }
    }

    parms <- get_correction_params(W.all)
    sim.standard <- simex(W.all, parms$Sigma.UU, est.fun=get_coef, M=10, B=50, models=list(Y~a+b/(c+X)), starts=list(list(a=1, b=1, c=1)))
    names(sim.standard) <- c("Intercept (SIMEX)", "X (SIMEX)")
    data.frame(t(sim.standard), scenario = sim3.data$scenario[1], check.names=FALSE)
}

experiment_3_empirical_simex <- function(sim3.data){
    X <- sim3.data$X
    Y <- sim3.data$Y
    W.all <- list(sim3.data$W1, sim3.data$W2, sim3.data$W3)

    # SIMEX
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(glm(Y[y.filter] ~ X, family=binomial(link='logit')))
        } else {
            coef(glm(Y~X, family=binomial(link='logit')))
        }
    }
    sim.empirical <- empirical_simex(W.all,Z = NULL, est.fun=get_coef, M=10, B=50, models=list(Y~a+b/(c+X)), starts=list(list(a=1, b=1, c=1)))
    names(sim.empirical) <- c("Intercept (E-SIMEX)", "X (E-SIMEX)")
    data.frame(t(sim.empirical), scenario = sim3.data$scenario[1], check.names=FALSE)
}