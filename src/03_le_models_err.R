######################################
# life expectancy models with error
# author: sebastian daza
#####################################3

# R < src/03_le_models_err.R > output/log/03_le_models_err_t.log  --no-save  &
# R < src/03_le_models_err.R > output/log/03_le_models_err_f.log  --no-save  &


# libraries, functions and options
library(data.table)
library(stringr)
library(brms)

library(doParallel)
cl = makeCluster(20, outfile="")
registerDoParallel(cl)
seed = 103231

library(texreg)
library(ggplot2)
source("src/utils.R")
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

table(idat$ctry)

# check multiple wy values
length(timps)
timps[[1]][ctry == 2020 & year == 1950, .(ctry, year, ex_mean, wy_mean, wy)]
timps[[3]][ctry == 2020 & year == 1950, .(ctry, year, ex_mean, wy_mean, wy)]

# values = NULL
# for (i in 1:100) {
#     values = c(values, timps[[1]][ctry == 2020 & year == 1950, wy])
# }
# summary(values)# replicates

nsamples = length(timps)
K = 10
print(paste0("Number of replicates: ", nsamples))
timps = timps[1:nsamples]

iterations = list(
    "f" = list(6000, 6000, 9000, 25000, 25000, 25500),
    "t" = list(4000, 4000, 5000, 12000, 12000, 12000)
)

multiResultClass = function(models = NULL, shifts = NULL) {
  me = list(models = models, shifts = shifts)
  class(me) = append(class(me), "multiResultClass")
  return(me)
}

#Progress combine function

# models 
output = foreach(i = 1:nsamples) %dopar% {

    results = multiResultClass()
    library(data.table)
    library(brms)
    library(loo)
    
    source("src/utils.R")
    
    print(paste0(":::::::: Running iteration ", i, " ::::::::"))
    models = list()
    dat = timps[[i]]

    models[[1]] = brm(wy ~ log_gdp + (1 | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[1]], 
        warmup = 1000, 
        chains = 1)

    models[[2]] = brm(wy ~ log_gdp + (log_gdp | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[2]], 
        warmup = 1000, 
        chains = 1, 
        control = list(max_treedepth = 15))

    models[[3]] = brm(wy ~ log_gdp + zyear + (log_gdp | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[3]], 
        warmup = 1000, 
        chains = 1, 
        control = list(max_treedepth = 15))

    models[[4]] = brm(wy ~ log_gdp + (1 | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[4]], 
        warmup = 1000, 
        chains = 1)

    models[[5]] = brm(wy ~ log_gdp + (log_gdp | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[5]], 
        warmup = 1000, 
        chains = 1, 
        control = list(max_treedepth = 15))

    models[[6]] = brm(wy ~ log_gdp + zyear + (log_gdp | ctryear), data = dat, 
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
        nsamples = 5000, K = K) 

    print(paste0(":::::::: Finished iteration ", i, " ::::::::"))
    return(results)
}

# extract results
model_replicates = list()
shifts = list()
weights = list()
nmodels = length(output[[1]]$models)
print(paste0("Number of models: ", nmodels))

for (i in seq_along(output)) {
   shifts[[i]] = output[[i]][["shifts"]][["values"]]
   weights[[i]] = output[[i]][["shifts"]][["weights"]]
}
shifts = rbindlist(shifts, idcol = "replicate")
avg_weights = apply(do.call(rbind, weights), 2, mean)

saveRDS(list("shifts" = shifts, "avg_weights" = avg_weights),
    paste0(data_path, select_estimates, "shifts_error.rds"))
rm(shifts)
slackr::slackr_msg(txt = paste0("LE error saved shifts: ", Sys.time()))

model_list = list()
lmodels = list()
for (i in 1:nmodels) {  
    for (h in 1:nsamples) {
        lmodels[[h]] = output[[h]][["models"]][[i]]
    }
    model_list[[i]] = combine_models(mlist = lmodels, check_data = FALSE)
    lmodels = list()
}
rm(lmodels, output)

tabs = list()
for (i in seq_along(model_list)) {
    print(paste0("Extracting model ", i))
    tabs[[i]] = extractBRMS(model_list[[1]], r2 = FALSE)
    model_list[[1]] = NULL
}
saveRDS(tabs, paste0(data_path, select_estimates, "tab_error.rds"))

# table
cnames = paste0("Model ", 1:length(model_list),  " (", round(avg_weights, 2), ")")
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", 
    "zyear" = "Year (standardized)")
caption = paste0("Models for LE and log GPD measurement error, stacking weight in parenthesis, ", 
    nsamples, " replicates")

texreg::texreg(tabs, 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = paste0("tab:", select_estimates, "ex_error"),
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = paste0(tables_path, select_estimates, "models_error.tex")
)    
file.copy(paste0(tables_path, select_estimates, "models_error.tex"), manus_tables, 
    recursive = TRUE)    

# send message to slack
slackr::slackr_msg(txt = paste0("LE models error finished at: ", Sys.time()))