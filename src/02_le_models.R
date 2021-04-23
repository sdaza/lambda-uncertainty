########################################
# life expectancy models without error
# author: sebastian daza
#######################################

# R < src/02_le_models.R > output/log/02_le_models.log  --no-save  &


# libraries, functions and options
library(data.table)
library(stringr)
library(brms)

library(doParallel)
cl = makeCluster(20)
registerDoParallel(cl)
seed = 103231

library(texreg)
library(ggplot2)
source("src/utils.R")
slackr::slackr_setup(config_file = ".slackr")
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
nsamples = length(timps)

# number of imputation to consider
nimputations = 10
timps = timps[1:nimputations]

iterations = list(
    "f" = list(6000, 6000, 9000, 25000, 25000, 25500),
    "t" = list(4000, 4000, 5000, 12000, 12000, 12000)
)

multiResultClass = function(models = NULL, shifts = NULL) {
  me = list(models = models, shifts = shifts)
  class(me) = append(class(me),"multiResultClass")
  return(me)
}

output = foreach(i = 1:nimputations) %dopar% {

    results = multiResultClass()
    library(data.table)
    library(brms)
    library(loo)
    source("src/utils.R")
    
    models = list()
    dat = timps[[i]]

    models[[1]] = brm(wy_mean ~ log_gdp + (1 | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[1]], 
        warmup = 1000, 
        chains = 1)

    models[[2]] = brm(wy_mean ~ log_gdp + (log_gdp | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[2]], 
        warmup = 1000, 
        chains = 1, 
        control = list(max_treedepth = 15))

    models[[3]] = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[3]], 
        warmup = 1000, 
        chains = 1, 
        control = list(max_treedepth = 15))

    models[[4]] = brm(wy_mean ~ log_gdp + (1 | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[4]], 
        warmup = 1000, 
        chains = 1)

    models[[5]] = brm(wy_mean ~ log_gdp + (log_gdp | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[5]], 
        warmup = 1000, 
        chains = 1, 
        control = list(max_treedepth = 15))

    models[[6]] = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[6]], 
        warmup = 1000, 
        chains = 1,
        control = list(max_treedepth = 15))

    results$models = models

    countries = unique(idat$ctry)
    years = c(1950)
    cyears = list(c("1950", "1950-1969"))
    dyears = list(c("1950", "1950+"))

    newdata = createComparisonData(dat, countries, 
        years, cyears, dyears)
    countries = unique(newdata$ctry)

    results$shifts = createShifts(models, newdata, countries = countries, 
        nsamples = 5000) 

    return(results)
}

# extract results
model_replicates = list()
shifts = list()
weights = list()
nmodels = length(output[[1]]$models)
for (i in seq_along(output)) {
   shifts[[i]] = output[[i]][["shifts"]][["values"]]
   weights[[i]] = output[[i]][["shifts"]][["weights"]]
}

shifts = rbindlist(shifts, idcol = "replicate")
avg_weights = apply(do.call(rbind, weights), 2, mean)

model_list = list()
lmodels = list()
for (i in 1:nmodels) {  
    for (h in 1:nimputations) {
        lmodels[[h]] = output[[h]][["models"]][[i]]
    }
    model_list[[i]] = combine_models(mlist = lmodels, check_data = FALSE)
    lmodels = list()
}
rm(lmodels, output)

saveRDS(list("shifts" = shifts, "avg_weights" = avg_weights),
    paste0(data_path, select_estimates, "shifts_no_error.rds"))
rm(shifts)
slackr::slackr_msg(txt = paste0("LE no error saved shifts: ", Sys.time()))

# create predicitive checks 
countries = unique(idat$ctry)
preferred_model = which.max(avg_weights)
vars = c("ctry", "ctry50", "ctryear", "zyear")

# stakcing
pred = rlang::invoke(pp_average, model_list, weights = avg_weights, 
    newdata = idat, summary = TRUE, nsamples = 10000)

plots_checks = prediction_check_plots(pred, idat, countries, vars, 
    y = "ex_mean", x = "year", 
    transform = TRUE, 
    country_labels = country_labs)
savepdf(paste0(plots_path, select_estimates, "fit_no_error_stacking"))
    print(plots_checks)
dev.off()
file.copy(paste0(plots_path, select_estimates, "fit_no_error_stacking.pdf"), 
    manus_plots, recursive = TRUE)    
rm(pred, plots_checks)

# preferred model 
pred = predict(model_list[[preferred_model]], newdata = idat, summary = TRUE)
plots_checks = prediction_check_plots(pred, idat, countries, 
    vars, y = "ex_mean", x = "year", 
    transform = TRUE, 
    country_labels = country_labs)
savepdf(paste0(plots_path, select_estimates, "fit_no_error_preferred_model"))
    print(plots_checks)
dev.off()
file.copy(paste0(plots_path, select_estimates, "fit_no_error_preferred_model.pdf"), 
    manus_plots, recursive = TRUE)    
rm(pred, plots_checks)

# table
cnames = paste0("Model ", 1:length(model_list),  " (", round(avg_weights, 2), ")")
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
file.copy(paste0(tables_path, select_estimates, "models_no_error.tex"), manus_tables, 
        recursive = TRUE)

# send message to slack
slackr::slackr_msg(txt = paste0("LE models error finished at: ", Sys.time()))