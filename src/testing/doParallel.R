library(rstan)
library(doParallel)

coin_tosses = c(1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
coin_toss_dat = list(N = length(coin_tosses), 
                    y = coin_tosses) 
stan_model = "
data {
    // # of obs; constrained to be greater than 0
    int<lower=0> N;  

    // define Y as an array of integers length N;
    //  each element either 0 and 1
    int<lower=0,upper=1> y[N]; 
    
}
parameters {
    real<lower=0,upper=1> theta;
}
model {
    // our prior for theta
    theta ~ beta(2,6); 

    // our likelihood for y
    y ~ bernoulli(theta);
}"

library(rstan)
library(doParallel)
cl = makeCluster(4)
registerDoParallel(cl)
output = foreach(i = 1:100,
        .packages = c("rstan", "rethinking")) %dopar% {

    f <- alist(
    y ~ dnorm( mu , sigma ),
    mu ~ dnorm( 0 , 10 ),
    sigma ~ dexp( 1 )
    )
    fit_stan = ulam(f , data=list(y=c(-1,1)) )
    return(fit_stan)

}
stopCluster(cl)

