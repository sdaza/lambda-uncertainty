########################################
# life expectancy models without error
# author: sebastian daza
#######################################

# R < /home/s/sdaza/00projects/lambda/src/le_models.R > /home/s/sdaza/00projects/lambda/output/le_models.log  --no-save  &
# R < src/02_le_models.R > output/log/02_le_models.log  --no-save  &


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
# options(mc.cores = 30, scipen=999)
# options(mc.cores = 10)
slackr::slackr_setup(config_file = ".slackr")
seed = 19380302
set.seed(seed)

# f (false) or t (true)
select_estimates = "f"

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
nsamples = length(timps)
nimputations = 10
timps = timps[1:nimputations]

# models 
m1 = brm_multiple(wy_mean ~ log_gdp + (1 | ctry50), data = timps, 
    iter = 6000, 
    warmup = 1000, 
    chains = 1, 
    combine = FALSE)

m2 = brm_multiple(wy_mean ~ log_gdp + (log_gdp | ctry50), data = timps, 
    iter = 6000, 
    warmup = 1000, 
    chains = 1, 
    control = list(max_treedepth = 15),
    combine = FALSE)

m3 = brm_multiple(wy_mean ~ log_gdp + zyear + (log_gdp | ctry50), data = timps, 
    iter = 9000, 
    warmup = 1000, 
    chains = 1, 
    control = list(max_treedepth = 15), 
    combine = FALSE)

m4 = brm_multiple(wy_mean ~ log_gdp + (1 | ctryear), data = timps, 
    iter = 25000, 
    warmup = 1000, 
    chains = 1, 
    combine = FALSE)

m5 = brm_multiple(wy_mean ~ log_gdp + (log_gdp | ctryear), data = timps, 
    iter = 25000, 
    warmup = 1000, 
    chains = 1, 
    control = list(max_treedepth = 15), 
    combine = FALSE)

m6 = brm_multiple(wy_mean ~ log_gdp + zyear + (log_gdp | ctryear), data = timps, 
    iter = 25000, 
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
for (i in 1:nimputations) {
    datasets[[i]] = createComparisonData(timps[[i]], countries, 
    years, cyears, dyears)
}

# countries with values before and after 1950
countries = unique(datasets[[1]]$ctry)
model_replicates = list(m1, m2, m3, m4, m5, m6)
rm(m1, m2, m3, m4, m5, m6)

# saveRDS(model_replicates, "output/models/model_replicates_error.rds")
slackr::slackr_msg(txt = paste0("LE no error replicates done!: ", Sys.time()))

# compute shifts
shifts = createShifts(model_replicates, datasets, countries = countries, 
    nsamples = 5000) 
saveRDS(shifts, paste0(data_path, select_estimates, "shifts_no_error.rds"))

# combine models
cm1 = combine_models(mlist = model_replicates[[1]], check_data = FALSE)
cm2 = combine_models(mlist = model_replicates[[2]], check_data = FALSE)
cm3 = combine_models(mlist = model_replicates[[3]], check_data = FALSE)
cm4 = combine_models(mlist = model_replicates[[4]], check_data = FALSE)
cm5 = combine_models(mlist = model_replicates[[5]], check_data = FALSE)
cm6 = combine_models(mlist = model_replicates[[6]], check_data = FALSE)
model_list = list(cm1, cm2, cm3, cm4, cm5, cm6)
rm(model_replicates, cm1, cm2, cm3, cm4, cm5, cm6)
# saveRDS(model_list, "output/models/model_combined_no_error.rds")

# create predicitive checks 
countries = unique(idat$ctry)
preferred_model = which.max(shifts[["avg_weights"]])

vars = c("ctry")

# stakcing
pred = rlang::invoke(pp_average, model_list, weights = shifts[["avg_weights"]], 
    newdata = idat, summary = TRUE, nsamples = 10000)
plots_checks = prediction_check_plots(pred, idat, countries, vars, 
    y = "ex_mean", x = "year", transform = TRUE, 
    country_labels = country_labs)
savepdf(paste0(plots_path, select_estimates, "fit_no_error_stacking")
    print(plots_checks)
dev.off()
file.copy(paste0(plots_path, select_estimates, "fit_no_error_stacking.pdf", 
    manus_plots, recursive = TRUE)    

# preferred model 
pred = predict(model_list[[preferred_model]], newdata = idat, summary = TRUE)
plots_checks = prediction_checks_plots(pred, idat, ctrys, 
    vars, y = "ex_mean", x = "year", 
    transform = TRUE, 
    country_labels = country_labs)

savepdf(paste0(plots_path, select_estimates, "fit_no_error_preferred_model")
    print(plots_checks)
dev.off()
file.copy(paste0(plots_path, select_estimates, "fit_no_error_preferred_model.pdf", 
    manus_plots, recursive = TRUE)    

# table
cnames = paste0("Model ", 1:length(model_list),  " (", round(shifts[["avg_weights"]], 2), ")")
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", 
    "zyear" = "Year (standardized)")
caption = paste0("Models for LE and log GPD no measurement error, stacking weight in parenthesis, ", 
    nimputations, " imputations")

texreg(model_list,
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:ex_no_error",
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
    file = paste0(tables_path, select_estimates, "models_no_error.tex")
)  

file.copy(paste0(tables_path, select_estimates, "models_no_error.tex"), manus_path, 
        recursive = TRUE)

# send message to slack
slackr::slackr_msg(txt = paste0("LE models no error done!: ", Sys.time()))