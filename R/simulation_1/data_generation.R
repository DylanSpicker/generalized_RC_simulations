generate_data_simulation_1 <- function(n) {
    X <- cbind(rnorm(n), rnorm(n, 3, sqrt(2)), rnorm(n, 1, sqrt(3)))
    Y <- 2 + X%*%c(-1,2,0.5) + rnorm(n)

    W1 <- X + cbind(rnorm(n), rnorm(n), rnorm(n))
    W2 <- X + cbind(rnorm(n), rnorm(n, 0, sqrt(4)), rnorm(n, 0, sqrt(3)))
    W3 <- X + cbind(rnorm(n, 0, sqrt(2)), rnorm(n, 0, sqrt(2)), rnorm(n, 0, sqrt(5)))

    # Missingness
    W2[1:(n%/%2), ] <- NA
    W3[1:floor(0.2*n), ] <- NA 
    
    data.frame(Y=Y,X=X,W1=W1,W2=W2,W3=W3)
   
}
