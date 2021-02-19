
library(data.table)
library(rethinking)
library(texreg)
library(brms)

source("src/utils.R")

library(rethinking)
data(WaffleDivorce)
d <- data.table(WaffleDivorce)

d[, D_obs := standardize(d$Divorce)]
d[, D_sd := d$Divorce.SE/sd(d$Divorce) * 1.7]
d[, M := standardize(Marriage)]
d[, A := standardize(MedianAgeMarriage)]

prior = c(prior(normal(0.0, 0.2), class = Intercept),
    prior(normal(0, 0.5), class = b),
    prior(exponential(1), class = sigma))

m1.2 = brm(D_obs | mi(D_sd) ~ 1 + A, family = gaussian,
    save_pars = save_pars(latent = TRUE),
    data = d, prior = prior, chains = 4, cores =4)
m1.1 = brm(D_obs ~ 1 + A, family = gaussian,
    data = d, prior = prior, chains = 4, cores =4)

screenreg(list(m1.1, m1.2))


## R code 15.3
dlist <- list(
    D_obs = standardize(d$Divorce),
    D_sd = d$Divorce.SE/sd(d$Divorce),
    M = standardize( d$Marriage ),
    A = standardize( d$MedianAgeMarriage ),
    N = nrow(d)
)

m2.1 <- ulam(
        alist(
            D_obs ~ dnorm(mu, sigma),
            mu <- a + bA*A + bM*M,
            a ~ dnorm(0,0.2),
            bA ~ dnorm(0,0.5),
            bM ~ dnorm(0,0.5),
            sigma ~ dexp(1)
        ) , data=dlist , chains=4 , cores=4 )

m2.2 <- ulam(
        alist(
            D_obs ~ dnorm( D_true , D_sd ),
            vector[N]:D_true ~ dnorm( mu , sigma ),
            mu <- a + bA*A + bM*M,
            a ~ dnorm(0,0.2),
            bA ~ dnorm(0,0.5),
            bM ~ dnorm(0,0.5),
            sigma ~ dexp(1)
        ) , data=dlist , chains=4 , cores=4 )

precis(m2.1)
precis(m2.2)

summary(m1.1)

prior = c(prior(normal(9.0, 0.2), class = Intercept),
    prior(normal(0, 0.5), class = b),
    prior(exponential(1), class = sigma))


m1.1 = brm(Divorce | mi(Divorce.SE) ~ 1 + A + M,
    data = dat, prior = prior, chains = 4, cores =4)
summary(m1.1)


summary(d)
