generate_data_simulation_2 <- function(n, ev1, ev2) {
    Z <- rbinom(n, size=1, prob=0.3)
    X <- rnorm(n, Z/50, sqrt(.5))
    Y <- rgamma(n, shape=1, scale=exp(2 - 3*Z + 2*X))

    a <- 0.7
    b <- sqrt(12*ev1) + a
    W1 <- X*runif(n, a, b)
    W2 <- X + rnorm(n, 0, sqrt(ev2))
    W3 <- X + rnorm(n, 0, sqrt(ev2))

    # Make W3 only available for 50% of individuals
    W3[sample(1:n, floor(0.5*n))] <- NA

    data.frame(Y,Z,X,W1,W2,W3)
}