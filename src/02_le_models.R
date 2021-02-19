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
setnames(adat, names(adat), tolower(names(adat)))
setnames(adat, "myear", "year")
adat = adat[sex == 1 & age == 0 & year >= 1900]
adat[, y := ex / max(ex + 1.05), by = ctry] # adjustment is by country!
adat[, wy := log(-log(1-y))]
adat[, max_le := max( ex + 1.05), by = ctry] # to recover values later

# year dataset
dat = data.table(read_stata("data/Ex_LA1850-2020_SES_FULL_Jan25-2021.dta"))
dat = dat[sex == 1 & age == 0 & year >= 1900 & year < 2020]
setnames(dat, names(dat), tolower(names(dat)))
dat[, y := ex / max(ex + 1.05), by = ctry] # adjustment is by country!
dat[, wy := log(-log(1-y))]
dat[, max_le := max(ex + 1.05), by = ctry] # to recover values later

# estimate standard deviation based on estimate
sdat = dat[, .(
    ex_mean = mean(ex),
    ex_random = getSample(ex),
    ex_sd = sd(ex),
    wy_mean = mean(wy),
    wy_sd = sd(wy), N = .N,
    gdp = getMin(gdp_pc),
    urban = getMin(urban),
    pop = getMin(pop)),
    .(ctry, year)]

sadat = adat[, .(
    ex_mean = mean(ex),
    ex_random = getSample(ex),
    ex_sd = sd(ex),
    wy_mean = mean(wy),
    wy_sd = sd(wy), N = .N,
    gdp = getMin(gdp_pc)),
    .(ctry, year)]

dim(sdat)
dim(sadat)

# some descriptives
summary(sdat$ex_mean)
summary(sdat$wy_mean)
summary(sadat$ex_mean)
summary(sadat$wy_mean)

# log gdp
sdat[, log_gdp := scale(log(gdp), scale = FALSE)]
sdat[, zpop := scale(pop)]
sdat[, zyear := scale(year)]

sadat[, log_gdp := scale(log(gdp), scale = FALSE)]
sadat[, zyear := scale(year)]

# recode year
sadat[year < 1950, gyear := "1950"]
sadat[year >= 1950 & year < 1970, gyear := "1950-1969"]
sadat[year >= 1970 & year < 1990, gyear := "1970-1989"]
sadat[year >= 1990, gyear := "1990"]
sadat[, ctryear := paste0(ctry,".", gyear)]
sadat[, year1950 := ifelse(year < 1950, "1950", "1950+")]

sdat[year < 1950, gyear := "1950"]
sdat[year >= 1950 & year < 1970, gyear := "1950-1969"]
sdat[year >= 1970 & year < 1990, gyear := "1970-1989"]
sdat[year >= 1990, gyear := "1990"]
sdat[, ctryear := paste0(ctry,".", gyear)]
sdat[, ctryearg := .GRP, ctryear]
# no gdp records before 1950 for country 2170
sdat[ctry == 2170 & ctryearg == 29, ctryearg := 30]
sdat[, year1950 := ifelse(year < 1950, "1950", "1950+")]

# flag missing records
sdat[, gdp_missing := ifelse(is.na(gdp), "missing", "observed")]
sdat[, sd_missing := ifelse(is.na(ex_sd), "missing", "observed")]
sadat[, gdp_missing := ifelse(is.na(gdp), "missing", "observed")]
sadat[, sd_missing := ifelse(is.na(ex_sd), "missing", "observed")]

# correlations
cor(sdat[, .(wy_mean, wy_sd, ex_mean, ex_sd, log_gdp)])
cor(sadat[, .(wy_mean, wy_sd, ex_mean, ex_sd, log_gdp)])

# missing data
countmis(sdat)
countmis(sadat)

# impute some missing data

# year data
setorder(sdat, ctry, year)

imp = mice(sdat, maxit = 0)
meth = imp$method
meth[] = ""
pred = imp$pred
pred[,] = 0

imp$loggedEvents

meth["zpop"] = "2l.pmm"
meth["log_gdp"] = "2l.pmm"
meth["ex_sd"] = "2l.pmm"

pred["zpop", c("zyear", "log_gdp", "ex_mean")] = 1
pred["zpop", c("ctryearg")] = -2

pred["ex_sd", c("zyear", "log_gdp", "ex_mean", "zpop")] = 1
pred["ex_sd", c("zyear")] = 2
pred["ex_sd", c("ctryearg")] = -2

pred["log_gdp", c("zyear", "ex_mean", "ex_sd", "zpop")] = 1
pred["log_gdp", c("zyear")] = 2
pred["log_gdp", c("ctryearg")] = -2

pred["log_gdp", ]

# five datasets and 10 iterations
imps = mice(sdat,
    m = 5,
    method = meth,
    maxit = 10,
    predictorMatrix = pred)

# create imputation plots
savepdf("manuscript/figures/mice")
    plot(imps)
    densityplot(imps, ~ ex_sd)
    densityplot(imps, ~ log_gdp)
    densityplot(imps, ~ zpop)
dev.off()

idat = data.table(complete(imps, action = 5))
countries = unique(sdat[gdp_missing == "missing", ctry])
savepdf("manuscript/figures/imputation_check_gdp")
for (i in countries) {
    print(ggplot(data = idat[ctry == i, .(log_gdp, year, gdp_missing)],
        aes(year, log_gdp, color = gdp_missing)) +
        geom_point() + labs(title = i))

}
dev.off()

countries = unique(sdat[sd_missing == "missing", ctry])
savepdf("output/plots/imputation_check_sd")
for (i in countries) {
    print(ggplot(data = idat[ctry == i, .(ex_sd, year, sd_missing)],
        aes(year, ex_sd, color = sd_missing)) +
        geom_point() + labs(title = i))

}
dev.off()

cor(idat[, .(ex_mean, wy_mean, log_gdp, ex_sd, wy_sd, N)])

savepdf("output/plots/gdp_sd")
    ggplot(idat, aes(log_gdp, wy_sd)) + geom_point() + theme_minimal()
    ggplot(idat, aes(log_gdp, ex_sd)) + geom_point() + theme_minimal()
dev.off()

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
savepdf("output/plots/mice_decade")
    plot(imps)
    densityplot(imps, ~ ex_sd)
    densityplot(imps, ~ log_gdp)
dev.off()

iadat = data.table(complete(imps, action = 5))
countries = unique(sadat[gdp_missing == "missing", ctry])

savepdf("output/plots/imputation_check_gdp_decade")
for (i in countries) {
    print(ggplot(data = iadat[ctry == i, .(log_gdp, year, gdp_missing)],
        aes(year, log_gdp, color = gdp_missing)) +
        geom_point() + labs(title = i))

}
dev.off()

countries = unique(sadat[sd_missing == "missing", ctry])
savepdf("output/plots/imputation_check_sd_decade")
for (i in countries) {
    print(ggplot(data = iadat[ctry == i, .(ex_sd, year, sd_missing)],
        aes(year, ex_sd, color = sd_missing)) +
        geom_point() + labs(title = i))

}
dev.off()

cor(iadat[, .(ex_mean, wy_mean, log_gdp, ex_sd, wy_sd, N)])

savepdf("output/plots/gdp_sd_decade")
    ggplot(iadat, aes(log_gdp, wy_sd)) + geom_point() + theme_minimal()
    ggplot(iadat, aes(log_gdp, ex_sd)) + geom_point() + theme_minimal()
dev.off()

# models using year data
m1 = brm(ex_mean ~ log_gdp + (1|ctryear),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(m1)

m2 = brm(ex_mean ~ log_gdp + (log_gdp|ctryear),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(m2)

m3 = brm(ex_mean ~ log_gdp + zyear + (zyear|ctry),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(m3)

m4 = brm(ex_mean ~ log_gdp + zyear + (log_gdp|ctry),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(m4)

baseline_models = list(m1, m2, m3, m4)
screenreg(baseline_models)

# error models 
e1 = brm(ex_mean | mi(ex_sd) ~ log_gdp + (1|ctryear),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(e1)

e2 = brm(ex_mean | mi(ex_sd) ~ log_gdp + (log_gdp|ctryear),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(e2)

e3 = brm(ex_mean | mi(ex_sd) ~  log_gdp + zyear + (1|ctry),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(e3)

e4 = brm(ex_mean | mi(ex_sd)  ~ log_gdp + zyear + (log_gdp|ctry),
    family = gaussian, data = idat,
    iter = 11000, 
    warmup = 1000, 
    chains = 10, 
    seed = 483892929,
    refresh = 11000, 
    cores = 10, 
    control = list(adapt_delta = 0.80)
)

summary(e4)

# decade







m2 = brm(ex_mean | mi(ex_sd) ~ log_gdp + (1|ctryear),
    family = gaussian, data = idat,
    iter = 10000,
    cores = 10,
    control = list(adapt_delta = 0.99),
    save_pars = save_pars(latent = TRUE)
    )

m3 = brm(ex_mean ~ log_gdp + (log_gdp|ctryear),
    family = gaussian, data = idat,
    iter = 10000,
    cores = 10,
    control = list(adapt_delta = 0.99)
    )

m4 = brm(ex_mean | mi(ex_sd) ~ log_gdp + (log_gdp|ctryear),
    family = gaussian, data = idat,
    iter = 10000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(latent = TRUE)
    )

screenreg(list(m1, m2, m3, m4), include.r2 = FALSE)

ctrys = unique(sdat$ctry)

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
                        co
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
