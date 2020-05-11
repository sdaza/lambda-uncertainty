#####################
# le models
# create shifts
#####################

# R < /home/s/sdaza/00projects/lambda/src/le_models.R > /home/s/sdaza/00projects/lambda/output/le_models.log  --no-save  &


# libraries, functions and options
library(brms)
library(data.table)
library(texreg)
library(stringr)

source('/home/s/sdaza/00projects/lambda/src/functions.R')
options(mc.cores = parallel::detectCores()-5)

seed = 103231

# read data
df = fread('/home/s/sdaza/00projects/lambda/data/featured_LE_data.csv')
setorder(df, ctry, year)

country_labels = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
               "Costa_Rica", "Cuba", "Dominican_Republic", "Ecuador",
               "El_Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua",
               "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela")

################################
# 1900 no time adjustment
################################

# set prior of betas to normal(0,5)
prior = set_prior("normal(0, 5)", class = "b")

m1.1 = brm(formula = wy ~ 1 + igdp_log  + (igdp_log|ctry_year),
           data = df,
           iter = 10000,
           chains = 5,
           seed = seed,
           prior=prior)

m1.2 = brm(formula = wy ~ 1 + igdp_log  + iurban_log +  (igdp_log|ctry_year),
          data = df,
          iter = 10000,
          chains = 5,
          seed = seed,
          prior = prior)

m1.3 = brm(formula = wy ~ 1 + igdp_log  + ilit_log +  (igdp_log|ctry_year),
          data = df,
          iter = 10000,
          chains = 5,
          seed = seed,
          prior = prior)

m1.4 = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + (igdp_log|ctry_year),
          data = df,
          iter = 10000,
          chains = 5,
          seed = seed,
          prior = prior,
          control=  list(adapt_delta=0.90))

loo1.1 = loo(m1.1, reloo=TRUE)
loo1.2 = loo(m1.2, reloo=TRUE)
loo1.3 = loo(m1.3, reloo=TRUE)
loo1.4 = loo(m1.4, reloo=TRUE)

# stacking
loo_list = list(loo1.1, loo1.2, loo1.3, loo1.4)
weights_models = as.vector(loo_model_weights(loo_list))
print(weights_models)

models = list(m1.1, m1.2, m1.3, m1.4)

# stacking

est_shifts = compute_shifts(models = models,
                        weights = weights_models,
                        data = df,
                        obs_var = 'Ex',
                        transform = TRUE,
                        posterior_nsample = 10000,
                        countries = country_labels,
                        years = c(1930, 1950, 1970, 1990, 2010))

fwrite(est_shifts, '/home/s/sdaza/00projects/lambda/output/shift_1900_stacking.csv')
saveRDS(weights_models, '/home/s/sdaza/00projects/lambda/output/model_weights_shift_1900_stacking.rds')

# # gdp only
# est_shifts = compute_shifts(models = list(m1.1),
#                         data = df,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1930, 1950, 1970, 1990, 2010))

# fwrite(est_shifts, '/home/s/sdaza/00projects/lambda/output/shift_1900_gdponly.csv')

# ################################
# # 1900 time adjustment
# ################################

# m1.1.year = brm(formula = wy ~ 1 + igdp_log + zyear + (igdp_log|ctry_year),
#            data = df,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# m1.2.year = brm(formula = wy ~ 1 + igdp_log  + iurban_log + zyear + (igdp_log|ctry_year),
#           data = df,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior)

# m1.3.year = brm(formula = wy ~ 1 + igdp_log  + ilit_log + zyear +  (igdp_log|ctry_year),
#           data = df,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior)

# m1.4.year = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + zyear + (igdp_log|ctry_year),
#           data = df,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior)


# loo1.1.year = loo(m1.1.year, reloo=TRUE)
# loo1.2.year = loo(m1.2.year, reloo=TRUE)
# loo1.3.year = loo(m1.3.year, reloo=TRUE)
# loo1.4.year = loo(m1.4.year, reloo=TRUE)

# # stacking
# models_year = list(m1.1.year, m1.2.year, m1.3.year, m1.4.year)

# loo_list_year = list(loo1.1.year, loo1.2.year, loo1.3.year, loo1.4.year)
# weights_models_year = as.vector(loo_model_weights(loo_list_year))

# # stacking
# est_shifts_year = compute_shifts(models = models_year,
#                         weights = weights_models_year,
#                         data = df,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1930, 1950, 1970, 1990, 2010))

# fwrite(est_shifts_year, '/home/s/sdaza/00projects/lambda/output/shift_1900_stacking_year.csv')

# # gdp only
# est_shifts_year = compute_shifts(models = list(m1.1.year),
#                         data = df,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1930, 1950, 1970, 1990, 2010))

# fwrite(est_shifts_year, '/home/s/sdaza/00projects/lambda/output/shift_1900_gdponly_year.csv')

# ################################
# # 1900 autocorrelation
# ################################

# m1.1.corr = brm(formula = wy ~ 1 + igdp_log + (igdp_log|ctry_year),
#           autocor = cor_ar(~wy|ctry, p = 1),
#           data = df,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior,
#           control=  list(adapt_delta=0.90))

# m1.4.corr = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + (igdp_log|ctry_year),
#           autocor = cor_ar(~wy|ctry, p = 1),
#           data = df,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior,
#           control=  list(adapt_delta=0.90))

# est_shifts_corr = compute_shifts(models = list(m1.1.corr),
#                         data = df,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1930, 1950, 1970, 1990, 2010))

# fwrite(est_shifts_corr, '/home/s/sdaza/00projects/lambda/output/shift_1900_gdponly_autocorrelation.csv')

# est_shifts_corr = compute_shifts(models = list(m1.4.corr),
#                         data = df,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1930, 1950, 1970, 1990, 2010))

# fwrite(est_shifts_corr, '/home/s/sdaza/00projects/lambda/output/shift_1900_full_autocorrelation.csv')


# # table

# ncountries_1900 = length(unique(df$ctry))
# avg_years_1900 = round(mean(df[, .N, ctry]$N), 0)

# name.map  <- list('Intercept' = "Intercept",
#                  igdp_log = "Log GDP",
#                  iurban_log = "Log Urbanity Rate",
#                  ilit_log = "Log Literacy Rate",
#                  zyear = 'Standardized Year',
#                  'sd(Intercept)' = 'Intercept SD',
#                  'sd(igdp_log)' = 'Log GDP Coef SD',
#                  'cor(Intercept,igdp_log)' = 'Correlation Intercept - GDP Coef'
#                  )

# models <- list(m1.1, m1.2, m1.3, m1.4, m1.4.year)

# texreg(models, include.r2 = TRUE, include.loo = TRUE,
#     stars = 0,
# #     custom.gof.names =c("Deaths", "Person-years"),
# #     custom.model.names = c("M1", "M1 MSM", "M2", "M2 MSM"),
#     groups = list("Fixed Effects" = 1:5, "Random Effects" = 6:8),
#     custom.coef.map = name.map,
#     custom.note = "Year blocks: $<$1950, 1950-1969, 1970-1989, 1990$+$. Ony the GPD coefficient is random.",
#     booktabs = TRUE,
#     dcolumn = TRUE,
#     use.packages = FALSE,
#     label = "le_1900",
#     caption = paste("Multilevel Models Weibull Transformation of Life Expectancy 1900,\\newline", ncountries_1900, "countries,",
#                     avg_years_1900, "average time points"),
#     caption.above = TRUE,
#     fontsize = "scriptsize",
#     float.pos = "htp",
#     file = "/home/s/sdaza/00projects/lambda/output/tables/le_1900.tex"
#     )


# ################################
# # 1950 no time adjustment
# ################################

# # select data since 1950
# df50 = copy(df[year>=1950])

# # set prior of betas to normal(0,5)
# prior = set_prior("normal(0, 5)", class = "b")

# m2.1 = brm(formula = wy ~ 1 + igdp_log + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# m2.2 = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior,
#            control=list(max_treedepth=15))

# m2.3 = brm(formula = wy ~ 1 + igdp_log + iwater_log + isewage_log + ielec_log + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# m2.4 = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log +
#             isewage_log + ielec_log + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# loo2.1 = loo(m2.1, reloo=TRUE)
# loo2.2 = loo(m2.2, reloo=TRUE)
# loo2.3 = loo(m2.3, reloo=TRUE)
# loo2.4 = loo(m2.4, reloo=TRUE)

# loo_2_list = list(loo2.1, loo2.2, loo2.3, loo2.4)
# weights_models_2 = as.vector(loo_model_weights(loo_2_list))
# models_2 = list(m2.1, m2.2, m2.3, m2.4)


# # stacking
# est_shifts = compute_shifts(models = models_2,
#                         weights = weights_models_2,
#                         data = df50,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1950, 1970, 1990, 2010))

# fwrite(est_shifts, '/home/s/sdaza/00projects/lambda/output/shift_1950_stacking.csv')

# # gdp only
# est_shifts = compute_shifts(models = list(m2.1),
#                         data = df50,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1950, 1970, 1990, 2010))

# fwrite(est_shifts, '/home/s/sdaza/00projects/lambda/output/shift_1950_gdponly.csv')


# ################################
# # 1950 time adjustment
# ################################

# m2.1.year = brm(formula = wy ~ 1 + igdp_log + zyear + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# m2.2.year = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + zyear + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# m2.3.year = brm(formula = wy ~ 1 + igdp_log + iwater_log + isewage_log + ielec_log + zyear + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# m2.4.year = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log +
#             isewage_log + ielec_log + zyear + (igdp_log|ctry_year),
#            data = df50,
#            iter = 2000,
#            chains = 2,
#            seed = seed,
#            prior=prior)

# loo2.1.year = loo(m2.1.year, reloo=TRUE)
# loo2.2.year = loo(m2.2.year, reloo=TRUE)
# loo2.3.year = loo(m2.3.year, reloo=TRUE)
# loo2.4.year = loo(m2.4.year, reloo=TRUE)

# loo_2_list_year = list(loo2.1.year, loo2.2.year, loo2.3.year, loo2.4.year)
# weights_models_2_year = as.vector(loo_model_weights(loo_2_list_year))
# models_2_year = list(m2.1.year, m2.2.year, m2.3.year, m2.4.year)


# # stacking
# est_shifts_year = compute_shifts(models = models_2_year,
#                         weights = weights_models_2_year,
#                         data = df50,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1950, 1970, 1990, 2010))

# fwrite(est_shifts_year, '/home/s/sdaza/00projects/lambda/output/shift_1950_stacking_year.csv')

# # gdp only
# est_shifts_year = compute_shifts(models = list(m2.1.year),
#                         data = df50,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1950, 1970, 1990, 2010))


# fwrite(est_shifts_year, '/home/s/sdaza/00projects/lambda/output/shift_1950_gdponly_year.csv')

# ################################
# # 1950 autocorrelation
# ################################

# m2.1.corr = brm(formula = wy ~ 1 + igdp_log + (igdp_log|ctry_year),
#           autocor = cor_ar(~wy|ctry, p = 1),
#           data = df50,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior,
#           control=  list(adapt_delta=0.90))

# m2.4.corr = brm(formula = wy ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log +
#             isewage_log + ielec_log + (igdp_log|ctry_year),
#           autocor = cor_ar(~wy|ctry, p = 1),
#           data = df50,
#           iter = 2000,
#           chains = 2,
#           seed = seed,
#           prior = prior,
#           control=  list(adapt_delta=0.90))

# # gdp only
# est_shifts_corr = compute_shifts(models = list(m2.1.corr),
#                         data = df50,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1950, 1970, 1990, 2010))

# fwrite(est_shifts_corr, '/home/s/sdaza/00projects/lambda/output/shift_1950_gdponly_autocorrelation.csv')

# # full
# est_shifts_corr = compute_shifts(models = list(m2.4.corr),
#                         data = df50,
#                         transform = TRUE,
#                         obs_var = 'Ex',
#                         countries = country_labels,
#                         years = c(1950, 1970, 1990, 2010))

# fwrite(est_shifts_corr, '/home/s/sdaza/00projects/lambda/output/shift_1950_full_autocorrelation.csv')


# # table

# ncountries_1950 = length(unique(df50$ctry))
# avg_years_1950 = round(mean(df50[, .N, ctry]$N), 0)

# name.map  <- list('Intercept' = "Intercept",
#                  igdp_log = "Log GDP",
#                  iurban_log = "Log Urbanity Rate",
#                  ilit_log = "Log Literacy Rate",
#                  iwater_log = "Log Water Rate",
#                  isewage_log = "Log Sewage Rate",
#                  ielec_log = "Log Electricity Rate",
#                  zyear = 'Standardized Year',
#                  'sd(Intercept)' = 'Intercept SD',
#                  'sd(igdp_log)' = 'Log GDP Coef SD',
#                  'cor(Intercept,igdp_log)' = 'Correlation Intercept - GDP Coef'
#                  )

# # length(name.map)

# models <- list(m2.1, m2.2, m2.3, m2.4, m2.4.year)

# texreg(models, include.r2 = TRUE, include.loo = TRUE,
#     stars = 0,
# #     custom.gof.names =c("Deaths", "Person-years"),
# #     custom.model.names = c("M1", "M1 MSM", "M2", "M2 MSM"),
#     groups = list("Fixed Effects" = 1:8, "Random Effects" = 9:11),
#     custom.coef.map = name.map,
#     custom.note = "Year blocks: 1950-1969, 1970-1989, 1990$+$. Ony the GPD coefficient is random.",
#     booktabs = TRUE,
#     dcolumn = TRUE,
#     use.packages = FALSE,
#     label = "le_1950",
#     caption = paste("Multilevel Models Weibull Transformation of Life Expectancy 1950,\\newline", ncountries_1950, "countries,",
#                     avg_years_1950, "average time points"),
#     caption.above = TRUE,
#     fontsize = "scriptsize",
#     float.pos = "htp",
#     file = "/home/s/sdaza/00projects/lambda/output/tables/le_1950.tex"
#     )
