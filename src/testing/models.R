##############################
# testing models
##############################

# libraries, functions and options
library(data.table)
library(stringr)
library(brms)

seed = 103231

library(texreg)
library(ggplot2)
source("src/utils.R")

seed = 19380302
set.seed(seed)

# f (false) or t (true)
select_estimates = "t"

# paths
plots_path = "output/plots/"
tables_path = "output/tables/"
data_path = "output/data/"
manus_plots = "manuscript/plots"
manus_tables  = "manuscript/tables"

# read data
data_list = readRDS(paste0(data_path, select_estimates, "datalist.rds"))
timps = data_list[["imputations"]]
idat = data_list[["single-imputation"]]
country_labs = data_list[["ctrylabels"]]

table(idat$ctry)

iterations = list(
    "f" = list(6000, 6000, 9000, 25000, 25000, 25500),
    "t" = list(4000, 4000, 5000, 12000, 12000, 12000)
)

# check multiple wy values
length(timps)
dat = timps[[2]]

bprior <- c(prior(normal(0,3), class = b),
            prior(cauchy(0,0.2), class = sd, 
                   group = year1950, coef = Intercept), 
            prior(cauchy(0,0.2), class = sd, 
                   group = year1950, coef = log_gdp))



m2 = brm(wy_mean ~ 1 + log_gdp + (1 | ctry50), data = dat, 
        iter = 6000, 
        chains = 1,
        warmup = 1000, 
        prior = bprior)

m2 = brm(wy_mean ~ 1 + log_gdp + (1 + log_gdp | ctry50), data = dat, 
        iter = 6000, 
        chains = 1,
        warmup = 1000, 
        prior = bprior)

m2 = brm(wy_mean ~ 1 + log_gdp + (log_gdp | ctryear), data = dat, 
        iter = 6000, 
        chains = 1,
        warmup = 1000, 
        prior = bprior)

m3 = brm(wy_mean ~ 1 + log_gdp + zyear + (log_gdp | ctryear), data = dat, 
        iter = 6000, 
        chains = 1,
        warmup = 1000, 
        prior = bprior)

summary(m2)
prior_summary(m2)

countries = unique(dat$ctry)
vars = c("ctry", "ctry50", "ctryear", "zyear", "year1950")
pred = predict(m2, newdata = dat, summary = TRUE)

plots_checks = prediction_check_plots(pred, dat, countries, 
    vars, y = "ex_mean", x = "year", 
    transform = TRUE, 
    country_labels = country_labs)

savepdf(paste0(plots_path, select_estimates, "fit_no_error_testing"))
    print(plots_checks)
dev.off()
