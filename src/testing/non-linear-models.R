############################
# non-linear models
############################


# libraries
library(data.table)
library(haven)
library(ggplot2)
library(brms)
source("src/utils.R")

# read data
dat = data.table(read_stata("data/Ex_LA1850-2020_SES_FULL_Jan25-2021.dta"))
setnames(dat, names(dat), tolower(names(dat)))

ctrylabs = attr(dat$ctry, "labels")
lab_list = as.list(setNames(attr(ctrylabs, "names"), as.numeric(ctrylabs)))
levels = as.numeric(ctrylabs)
labs = attr(ctrylabs, "names")
dat[, ctryf := factor(ctry, labels = labs, level = levels)]

# male
dat = dat[sex == 1 & age == 0 & year >= 1900 & year < 2010]
# remove LE duplicates
dat[, N := 1:.N, .(ctry, year, ex)]
dat = dat[N == 1]
anyDuplicated(dat[, .(ctry, year, ex)])
nrow(dat)

# select estimates
dat[, N := 1:.N, .(ctry, year)]
dat[, lambda := 0][tseries2 == 1, lambda := 1]
dat = dat[lambda == 1]
dat[, N := .N, .(ctry, year)]
dat = dat[N == 1]
dat[, max_ex := max(ex)]
dat[, tle := ex / (1.05 * max(ex))]
dat[, wle := transWeibull(ex, max_ex)]

dat[, gdp_pc2 := gdp_pc^2]
dat[, gdp_pc3 := gdp_pc^3]
dat[, log_gdp := scale(log(gdp_pc), scale = FALSE)]
dat[, log_gdp2 := log_gdp^2]
dat[, ctry50 := ctry *10 + ifelse(year < 1950, 0, 1)]

m1 = brm(wle ~ 1 + t2(log_gdp) + (1 + log_gdp|ctry50), data = dat, 
     control = list(adapt_delta = 0.95), 
     iter = 6000)
summary(m1)

savepdf("output/plots/preston_conditional")
plot(conditional_effects(m1), points = TRUE)da
dev.off()

conditions  = data.frame(ctry = unique(dat$ctry))
rownames(conditions) = unique(dat$ctry)

me_ctry = conditional_effects(m1, conditions = conditions,
    re_formula = NULL, method = "predict")

savepdf("output/plots/m1_cond_ctry")
plot(me_ctry, points = TRUE, cols = 2)
dev.off()

# model specification 
f = bf(tle ~ alpha / (1 + exp(beta + exp(log(C) * gdp_pc))), 
         alpha ~ 1, 
         beta ~ 1 + (1|ID1|year), 
         C ~ 1 + (1|ID1|year), nl = TRUE)

priors = c(prior(normal(5.4, 0.1), nlpar='beta'),
            prior(uniform(0.5,1), lb = 0.5, ub = 1, nlpar = 'C'), 
            prior(normal(201.8, 0.1), nlpar = 'alpha'))

preston = brm(f, data = dat, prior = priors, 
    iter = 15000,
    chains = 1, 
    cores = 1,
    control = list(adapt_delta = 0.90))

summary(preston)

savepdf("output/plots/preston_conditional")
plot(conditional_effects(preston), points = TRUE)
dev.off()

selected_years = c(1940, 1950, 1960, 1970, 1980, 1990)
conditions  = data.frame(year = unique(dat$ctry))
rownames(conditions) = unique(dat$ctry)

me_ctry = conditional_effects(preston, conditions = conditions,
    re_formula = NULL, method = "predict")

savepdf("output/plots/preston_conditional_year")
plot(me_year, points = TRUE)
dev.off()



selected_years = c(1940, 1950, 1960, 1970, 1980, 1990)
conditions  = data.frame(year = selected_years)
rownames(conditions) = selected_years

me_year = conditional_effects(preston, conditions = conditions,
    re_formula = NULL, method = "predict")

savepdf("output/plots/preston_conditional_year")
plot(me_year, points = TRUE)
dev.off()
