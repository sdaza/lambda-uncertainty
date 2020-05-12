# project: lambda
# author: sebastian
# run models in iteration

# run at linstat
# R < /home/s/sdaza/00projects/lambda/src/le_boostrap_1900_uniform_prob.R > /home/s/sdaza/00projects/lambda/src/le_boostrap_1900_uniform_prob.log  --no-save  &


library(doParallel)
library(data.table)

cl = makeCluster(15)
registerDoParallel(cl)
seed = 103231

df = fread('/home/s/sdaza/00projects/lambda/data/bs_samples_uniform.csv')
setorder(df, ctry, year)

country_labels = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
               "Costa_Rica", "Cuba", "Dominican_Republic", "Ecuador",
               "El_Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua",
               "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela")

sample_size = max(df$sample_index)


results = foreach(i=1:sample_size) %dopar% {

    library(data.table)
    library(brms)
    library(loo)
    library(stringr)
    library(texreg)
    source('/home/s/sdaza/00projects/lambda/src/functions.R')

    prior = set_prior("normal(0, 5)", class = "b")

    test = copy(df[sample_index==i])

    # adjustment is by country!
    test[, y := le/max(le+1.05), by = ctry]
    test[, wy := log(-log(1-y))]

    # to recover values later
    test[, max_le := max(le+1.05), by = ctry]

    m1.1 = brm(formula = wy ~ 1 + igdp_log  + (igdp_log|ctry_year),
           data = test,
           iter = 2000,
           chains = 2,
           seed = seed,
           prior=prior)

    m1.2 = brm(formula = wy ~ 1 + igdp_log  + iurban_log +  (igdp_log|ctry_year),
          data = test,
          iter = 2000,
          chains = 2,
          seed = seed,
          prior = prior)

    m1.3 = brm(formula = wy ~ 1 + igdp_log  + ilit_log +  (igdp_log|ctry_year),
          data = test,
          iter = 2000,
          chains = 2,
          seed = seed,
          prior = prior)

    m1.4 = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + (igdp_log|ctry_year),
          data = test,
          iter = 2000,
          chains = 2,
          seed = seed,
          prior = prior,
          control=  list(adapt_delta=0.90))

    loo1.1 = loo(m1.1, reloo=TRUE)
    loo1.2 = loo(m1.2, reloo=TRUE)
    loo1.3 = loo(m1.3, reloo=TRUE)
    loo1.4 = loo(m1.4, reloo=TRUE)

    models  = list(m1.1, m1.2, m1.3, m1.4)
    loo_list = list(loo1.1, loo1.2, loo1.3, loo1.4)
    model_weights = as.vector(loo_model_weights(loo_list))

    compute_shifts(models = models,
                   weights = model_weights,
                   posterior_nsample = 100,
                   data = test,
                   obs_var = 'le',
                   transform = TRUE,
                   countries = country_labels,
                   years = c(1930, 1950, 1970, 1990, 2010))

}

r = rbindlist(results, idcol='sample_index')
saveRDS(r, "/home/s/sdaza/00projects/lambda/output/shifts_1900_uniform.rds")
