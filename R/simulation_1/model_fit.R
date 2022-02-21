get_w_matrix <- function(sim1.df) {
    list(sim1.df[,c("W1.1", "W1.2", "W1.3")], 
         sim1.df[,c("W2.1", "W2.2", "W2.3")], 
         sim1.df[,c("W3.1", "W3.2", "W3.3")])
}

# Perform the analysis using standard regression calibration
experiment_1_RC <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)
    Wstar.standard <- RC(W)
    standard <- coef(lm(sim1.df$Y~Wstar.standard))
    names(standard) <- c("Intercept (RC)", "X1 (RC)", "X2 (RC)", "X3 (RC)")
    data.frame(t(standard), scenario = sim1.df$scenario[1], check.names=FALSE)
}

# Perform the analysis using generalized regression calibration
experiment_1_generalized_RC <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)
    Wstar.general <- generalized_RC(W, Z=NULL, J0=c(1,2,3), J1=c(1,2,3))
    general <- coef(lm(sim1.df$Y~Wstar.general))
    names(general) <- c("Intercept (Generalized RC)", "X1 (Generalized RC)", "X2 (Generalized RC)", "X3 (Generalized RC)")
    data.frame(t(general), scenario = sim1.df$scenario[1], check.names=FALSE)
}

# Perform the analysis using generalized regression calibration,
# with optimal weighting
experiment_1_weighted_generalized_RC <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)
    Wstar.general.wt <- generalized_RC.weights(W, Z=NULL, J0=c(1,2,3), J1=c(1,2,3), return_parms=FALSE)
    general.wt <- coef(lm(sim1.df$Y~Wstar.general.wt))

    names(general.wt) <- c("Intercept (Generalized RC Weighted)", "X1 (Generalized RC Weighted)", "X2 (Generalized RC Weighted)", "X3 (Generalized RC Weighted)")
    data.frame(t(general.wt), scenario = sim1.df$scenario[1], check.names=FALSE)
}

# Perform the analysis using standard simulation extrapolation.
experiment_1_SIMEX <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)

    # Defined within the function to be based on the correct 'Y'
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(lm(sim1.df$Y[y.filter]~X)) 
        } else {
            coef(lm(sim1.df$Y~X)) 
        }
    }

    parms <- get_correction_params(W)
    standard.simex <- simex(W,
                            parms$Sigma.UU,
                            est.fun=get_coef,
                            M=15,
                            models=list(Y~a+b/(c+X)), 
                            starts=list(list(a=1, b=1, c=1)))
    names(standard.simex) <- c("Intercept (SIMEX)", "X1 (SIMEX)", "X2 (SIMEX)", "X3 (SIMEX)")
    data.frame(t(standard.simex), scenario = sim1.df$scenario[1], check.names=FALSE)
}

# Perform the analysis using standard 'emprical' simulation extrapolation
experiment_1_ESIMEX <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(lm(sim1.df$Y[y.filter]~X)) 
        } else {
            coef(lm(sim1.df$Y~X)) 
        }
    }

    empirical.simex.est <- empirical_simex(W, 
                                           est.fun=get_coef,
                                           M=15,
                                           models=list(Y~a+b/(c+X)), 
                                           starts=list(list(a=1, b=1, c=1)))
    names(empirical.simex.est) <- c("Intercept (E-SIMEX)", "X1 (E-SIMEX)", "X2 (E-SIMEX)", "X3 (E-SIMEX)")
    data.frame(t(empirical.simex.est), scenario = sim1.df$scenario[1], check.names=FALSE)
}

# Perform the analysis using generalized SIMEX where averaging 
# takes place *after* 
experiment_1_generalized_SIMEX <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)    
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(lm(sim1.df$Y[y.filter]~X)) 
        } else {
            coef(lm(sim1.df$Y~X)) 
        }
    }

    sim.general <- generalized_simex(W,
                                J0=c(1,2,3),
                                J1=c(1,2,3),
                                use.proxies=c(1,2,3), 
                                est.fun=get_coef,
                                M=15,
                                models=list(Y~a+b/(c+X)), 
                                starts=list(list(a=1, b=1, c=1)), 
                                save.object=TRUE,
                                combine.first=FALSE)

    sim.ests <- unlist(sim.general$estimates)
    names(sim.ests) <- c("Intercept (Generalized SIMEX (A))", "X1 (Generalized SIMEX (A))", "X2 (Generalized SIMEX (A))", "X3 (Generalized SIMEX (A))")
    data.frame(t(sim.ests), scenario = sim1.df$scenario[1], check.names=FALSE)
}

# Perform the analysis using generalized SIMEX where averaging 
# takes place *before*
experiment_1_combine_first_generalized_SIMEX <- function(sim1.df) {
    W <- get_w_matrix(sim1.df)    
    get_coef <- function(X, y.filter=NULL){ 
        if(!is.null(y.filter)) {
            coef(lm(sim1.df$Y[y.filter]~X)) 
        } else {
            coef(lm(sim1.df$Y~X)) 
        }
    }

    sim.general.combine <- generalized_simex(W,
                                            J0=c(1,2,3),
                                            J1=c(1,2,3),
                                            use.proxies=c(1,2,3), 
                                            est.fun=get_coef,
                                            M=15,
                                            models=list(Y~a+b/(c+X)), 
                                            starts=list(list(a=1, b=1, c=1)), 
                                            save.object=TRUE)
    sim.ests.com <- unlist(sim.general.combine$estimates)
    names(sim.ests.com) <- c("Intercept (Generalized SIMEX (B))", "X1 (Generalized SIMEX (B))", "X2 (Generalized SIMEX (B))", "X3 (Generalized SIMEX (B))")
    data.frame(t(sim.ests.com), scenario = sim1.df$scenario[1], check.names=FALSE)
}