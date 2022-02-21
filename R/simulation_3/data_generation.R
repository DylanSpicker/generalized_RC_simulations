generate_data_simulation_3 <- function(n, eta0, eta1) {
    X <- rnorm(n, mean=3, sd=1)
    Y <- rbinom(n, 1, expit(0.5 - 0.5*X))

    W1 <- X + rnorm(n)
    W2 <- X + rnorm(n)
    W3 <- eta0 + eta1*X + runif(n, -0.5, 0.5)
    
    W2[sample(1:n, floor(0.8*n))] <- NA

    data.frame(Y=Y, X=X, W1=W1, W2=W2, W3=W3)
}