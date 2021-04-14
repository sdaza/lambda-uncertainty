############################
# life expectancy models
# author: sebastian daza
############################

# R < /home/s/sdaza/00projects/lambda/src/le_models.R > /home/s/sdaza/00projects/lambda/output/le_models.log  --no-save  &

# libraries, functions and options
library(haven)
library(data.table)
library(stringr)
library(brms)
library(loo)
library(texreg)
library(mice)
library(miceadds)

library(ggplot2)
library(patchwork)

options(scipen=999)

source('src/utils.R')

options(mc.cores = 10)
seed = 103231

# decades
adat = data.table(read_stata("data/Ex_LA1840-2020_UncertaintyFile_bydecades.dta"))
adat
setnames(adat, names(adat), tolower(names(adat)))
setnames(adat, "myear", "year")

adat = adat[sex == 1 & age == 0 & year >= 1900]

adat[ctry == 2020 & year < 1950]
adat[year == 1905 & ctry == 2020, .(ctry, year, sex, age, ex)]

adat[, y := ex / max(ex + 1.05), by = ctry] # adjustment is by country!
adat[, wy := log(-log(1-y))]
adat[, max_le := max( ex + 1.05), by = ctry] # to recover values later

sadat = adat[, .(
    ex_mean = mean(ex),
    ex_random = getSample(ex),
    ex_sd = sd(ex),
    wy_mean = mean(wy),
    wy_sd = sd(wy), N = .N,
    gdp = getMin(gdp_pc)),
    .(ctry, year)]

dim(sadat)

# some descriptives
summary(sadat$ex_mean)
summary(sadat$wy_mean)

# log gdp
sadat[, log_gdp := scale(log(gdp), scale = FALSE)]
sadat[, zyear := scale(year)]

# recode year
sadat[year < 1950, gyear := "1950"]
sadat[year >= 1950 & year < 1970, gyear := "1950-1969"]
sadat[year >= 1970 & year < 1990, gyear := "1970-1989"]
sadat[year >= 1990, gyear := "1990"]
sadat[, ctryear := paste0(ctry,".", gyear)]
sadat[, year1950 := ifelse(year < 1950, "1950", "1950+")]

# flag missing records
sadat[, gdp_missing := ifelse(is.na(gdp), "missing", "observed")]
sadat[, sd_missing := ifelse(is.na(ex_sd), "missing", "observed")]

# correlations
cor(sadat[, .(wy_mean, wy_sd, ex_mean, ex_sd, log_gdp)])

# missing data
countmis(sadat)

# imputation using decade data
setorder(sadat, ctry, year)

imp = mice(sadat, maxit = 0)
meth = imp$method
meth[] = ""
pred = imp$pred
pred[,] = 0

countmis(sadat)
imp$loggedEvents

meth["log_gdp"] = "2l.pmm"
meth["ex_sd"] = "2l.pmm"

pred["ex_sd", c("zyear", "log_gdp", "ex_mean")] = 1
pred["ex_sd", c("zyear")] = 1
pred["ex_sd", c("ctry")] = -2

pred["log_gdp", c("zyear", "ex_mean", "ex_sd")] = 1
pred["log_gdp", c("zyear")] = 1
pred["log_gdp", c("ctry")] = -2

pred["log_gdp", ]

# five datasets and 10 iterations
imps = mice(sadat,
    m = 5,
    method = meth,
    maxit = 10,
    predictorMatrix = pred)

# create imputation plots
savepdf("manuscript/figures/mice_decade")
    plot(imps)
    densityplot(imps, ~ ex_sd)
    densityplot(imps, ~ log_gdp)
dev.off()

iadat = data.table(complete(imps, action = 1))
countries = unique(sadat[gdp_missing == "missing", ctry])

savepdf("manuscript/figures/imputation_check_gdp_decade")
for (i in countries) {
    print(ggplot(data = iadat[ctry == i, .(log_gdp, year, gdp_missing)],
        aes(year, log_gdp, color = gdp_missing)) +
        geom_point() + labs(title = i,  x = "Year", y = "Log GDP", color = NULL)
    )

}
dev.off()

countries = unique(sadat[sd_missing == "missing", ctry])
savepdf("manuscript/figures/imputation_check_sd_decade")
for (i in countries) {
    print(ggplot(data = iadat[ctry == i, .(ex_sd, year, sd_missing)],
        aes(year, ex_sd, color = sd_missing)) +
        geom_point() + labs(title = i, x = "Year", y = "SD life expectancty", color = NULL color = NULL)
    )

}
dev.off()

table(iadat$sd_missing)
cor(iadat[, .(ex_mean, log_gdp, ex_sd, N)])
cor(iadat[sd_missing != "missing", .(ex_mean, log_gdp, ex_sd, N)])

savepdf("manuscript/figures/gdp_sd_decade")    
    ggplot(iadat, aes(log_gdp, ex_sd, color = sd_missing)) + geom_point() + theme_minimal() +
    labs(x = "Log GDP", y = "SD life expectancy", color = NULL)
dev.off()

# models

# using year data
m1 = brm(ex_mean ~ log_gdp + (1 | ctryear),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_all_pars = TRUE,
    control = list(adapt_delta = 0.80)
)

summary(m1)

m2 = brm(ex_mean ~ log_gdp + (log_gdp | ctryear),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_all_pars = TRUE,
    control = list(adapt_delta = 0.99)
)

summary(m2)

m3 = brm(ex_mean ~ log_gdp + zyear + (zyear|ctry),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_all_pars = TRUE,
    control = list(adapt_delta = 0.99)
)

summary(m3)

m4 = brm(ex_mean ~ log_gdp + zyear + (zyear + log_gdp | ctry),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_all_pars = TRUE,
    control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(m4)

baseline_models = list(m1, m2, m3, m4)
loo_list = list()
for (i in seq_along(baseline_models)) {
    loo_list[[i]] = loo(baseline_models[[i]], reloo=TRUE)
}
model_weights = as.vector(loo_model_weights(loo_list))
round(model_weights, 3)


# predictive checks
ctrys = unique(sadat$ctry)
vars = c("ctry", "year", "log_gdp", "ex_mean")
plots_checks = prediction_checks_ex(m4, sadat, ctrys, vars, y = "ex_mean", x = "year")

savepdf("manuscript/figures/fit_no_error_m4_decades")
    print(plots_checks)
dev.off()


# create table
cnames = paste0("Model ", 1:4,  " (", round(model_weights, 4), ")")
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", "zyear" = "Year (standardized)")

texreg(baseline_models,
    caption = "Models for LE and log GPD (no measurement error), decades", 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:ex_no_error",
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = "manuscript/tables/models_no_error_decades.tex"
)



# error models 
e1 = brm(ex_mean | mi(ex_sd) ~ log_gdp + (1|ctryear),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.80)
)

summary(e1)

e2 = brm(ex_mean | mi(ex_sd) ~ log_gdp + (log_gdp|ctryear),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.99)
)

summary(e2)

e3 = brm(ex_mean | mi(ex_sd) ~ log_gdp + zyear + (zyear|ctry),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.99)
)

summary(e3)

e4 = brm(ex_mean | mi(ex_sd)  ~ log_gdp + zyear + (zyear + log_gdp|ctry),
    family = gaussian, data = iadat,
    iter = 15000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(e4)

baseline_models_e = list(e1, e2, e3, e4)
loo_list_e = list()
for (i in seq_along(baseline_models_e)) {
    loo_list_e[[i]] = loo(baseline_models_e[[i]], reloo = TRUE)
}
model_weights = as.vector(loo_model_weights(loo_list_e))
round(model_weights, 3)



# checking fit
vars = c("ctry", "year", "log_gdp", "ex_mean")

plots1 = prediction_checks_ex(m1, sdat, ctrys, vars, y = "ex_mean", x = "year")
plots2 = prediction_checks_ex(m2, sdat, ctrys, vars, y = "ex_mean", x = "year")

plots3 = prediction_checks_ex(m3, sdat, ctrys, vars, y = "ex_mean", x = "year")
plots4 = prediction_checks_ex(m4, sdat, ctrys, vars, y = "ex_mean", x = "year")


savepdf("fit_m1_m2")
print(plots1[[1]])
print(plots2[[1]])
print(plots3[[1]])
print(plots4[[1]])
# wrap_plots(plots, ncol = 3)
dev.off()

sdat[, .(sum(is.na(ex_sd)), mean(ex_sd, na.rm = TRUE)), ctry]


est_shifts = compute_shifts(models = list(m1),
                        data = sdat,
                        obs_var = 'ex_mean',
                        transform = FALSE,
                        posterior_nsample = 100,
                        years = c(1950, 1970, 1990))

p1 = posterior_samples(m1, pars = "gdp")
p2 = posterior_samples(m2, pars = "gdp")
p3 = posterior_samples(m3, pars = "gdp")
p4 = posterior_samples(m4, pars = "gdp")
test1 =  p2$b_log_gdp - p1$b_log_gdp
test2 =  p4$b_log_gdp - p3$b_log_gdp

quantile(test1, prob = c(0.05, 0.95))
quantile(test2, prob = c(0.05, 0.95))

savepdf("test_m1_m2")
    hist(test1)
dev.off()

savepdf("test_m3_m4")
    hist(test2)
dev.off()
