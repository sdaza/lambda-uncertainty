############################
# life expectancy models
# author: sebastian daza
############################

# R < /home/s/sdaza/00projects/lambda/src/le_models.R > /home/s/sdaza/00projects/lambda/output/le_models.log  --no-save  &
# R < src/03_le_models_err.R > output/log/03_le_models_err.log  --no-save  &


# libraries, functions and options
library(haven)
library(data.table)
library(stringr)
library(brms)
library(loo)
library(future)
plan(multiprocess, workers = 15)
options(future.globals.maxSize = 8000 * 1024^2)
options(future.rng.onMisuse="ignore")

library(texreg)
library(ggplot2)
library(patchwork)
source("src/utils.R")
# options(mc.cores = 15)

# slack setup
slackr::slackr_setup(config_file = ".slackr")
seed = 19380302
set.seed(seed)

# read data
# sdat = readRDS("output/data/aggregate-data.rds")
idat = readRDS("output/data/single-imputation.rds")
timps = readRDS("output/data/imputations.rds")
nsamples = length(timps)

# check values
# values = NULL
# for (i in 1:100) {
    # values = c(values, timps[[i]][ctry == 2140 & year == 1950, wy])
# }

# models 
m1 = brm_multiple(wy ~ log_gdp + (1 | ctry50), data = timps, 
    iter = 4000, 
    warmup = 1000, 
    chains = 1, 
    combine = FALSE)

m2 = brm_multiple(wy ~ log_gdp + (log_gdp | ctry50), data = timps, 
    iter = 4000, 
    warmup = 1000, 
    chains = 1,
    combine = FALSE)

m3 = brm_multiple(wy ~ log_gdp + zyear + (log_gdp | ctry50), data = timps, 
    iter = 4000, 
    warmup = 1000, 
    chains = 1,  
    control = list(max_treedepth = 15), 
    combine = FALSE)

m4 = brm_multiple(wy ~ log_gdp + (1 | ctryear), data = timps, 
    iter = 12000, 
    warmup = 1000, 
    chains = 1, 
    combine = FALSE)

m5 = brm_multiple(wy ~ log_gdp + (log_gdp | ctryear), data = timps, 
    iter = 12000, 
    warmup = 1000, 
    chains = 1,  
    control = list(max_treedepth = 15), 
    combine = FALSE)

m6 = brm_multiple(wy ~ log_gdp + zyear + (log_gdp | ctryear), data = timps, 
    iter = 12000, 
    warmup = 1000, 
    chains = 1,  
    control = list(max_treedepth = 15),
    combine = FALSE)

# create datasets to contrast
countries = unique(idat$ctry)
years = c(1950)
cyears = list(c("1950", "1950-1969"))
dyears = list(c("1950", "1950+"))

datasets = list()
for (i in 1:nsamples) {
    datasets[[i]] = createComparisonData(timps[[i]], countries, 
    years, cyears, dyears)
}

# countries with values before and after 1950
countries = unique(datasets[[1]]$ctry)
model_replicates = list(m1, m2, m3, m4, m5, m6)
rm(m1, m2, m3, m4, m5, m6)
# saveRDS(model_replicates, "output/models/model_replicates_error.rds")
slackr::slackr_msg(txt = paste0("LE error replicates: ", Sys.time()))

# compute shifts
shifts = createShifts(model_replicates, datasets, countries = countries, 
    nsamples = 3000) 
saveRDS(shifts, "output/models/shifts_error.rds")
slackr::slackr_msg(txt = paste0("LE error saved shifts: ", Sys.time()))

cm1 = combine_models(mlist = model_replicates[[1]], check_data = FALSE)
cm2 = combine_models(mlist = model_replicates[[2]], check_data = FALSE)
cm3 = combine_models(mlist = model_replicates[[3]], check_data = FALSE)
cm4 = combine_models(mlist = model_replicates[[4]], check_data = FALSE)
cm5 = combine_models(mlist = model_replicates[[5]], check_data = FALSE)
cm6 = combine_models(mlist = model_replicates[[6]], check_data = FALSE)
model_list = list(cm1, cm2, cm3, cm4, cm5, cm6)
rm(model_replicates, cm1, cm2, cm3, cm4, cm5, cm6)
# saveRDS(model_list, "output/models/model_combined_error.rds")

# table
cnames = paste0("Model ", 1:length(model_list),  " (", round(shifts[["avg_weights"]], 2), ")")
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", "zyear" = "Year (standardized)")
caption = paste0("Models for LE and log GPD measurement error, stacking weight in parenthesis, ", 
    nsamples, " replicates")

texreg(model_list, 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:ex_error",
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    include.loo.ic = FALSE,
    inclue.rsquared = TRUE,
    include.waic = FALSE,
    include.random = FALSE,
    reloo = FALSE,
    file = "output/tables/models_error.tex"
)    
file.copy("output/tables/models_error.tex", "manuscript/tables/", recursive = TRUE)    


# send message to slack
slackr::slackr_msg(txt = paste0("LE models no error finished at: ", Sys.time()))
