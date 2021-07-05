########################################
# life expectancy models without error
# author: sebastian daza
#######################################

# R < src/02_le_models.R > output/log/02_le_models_t.log  --no-save  &
# R < src/02_le_models.R > output/log/02_le_models_f.log  --no-save  &


# libraries, functions and options
library(data.table)
library(stringr)
library(brms)

library(doParallel)
cl = makeCluster(10, outfile="")
registerDoParallel(cl)
seed = 103231

library(texreg)
library(ggplot2)
source("src/utils.R")
slackr::slackr_setup(config_file = ".slackr")
seed = 19380302
set.seed(seed)

# cross-validation
K = 10 

# stacking weights
stacking = TRUE

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
covs = data_list[["covs"]]
country_labs = data_list[["ctrylabels"]]
nsamples = length(timps)

# number of imputation to consider
nimputations = 10
timps = timps[1:nimputations]

iterations = list(
    "f" = list(6000, 6000, 6000, 14000, 12000, 10000),
    "t" = list(6000, 6000, 6000, 14000, 12000, 10000)
)

# # data exploration
# dat = timps[[1]]
# dat[, year1950 := factor(year1950)]
# dat[, gyear := factor(gyear)]

# savepdf("output/plots/scatter")
# print(
#     ggplot(dat, aes(x = log_gdp, y = wy_mean,color = year1950)) +
#     geom_point() + geom_smooth(method=lm, se = TRUE) + 
#     geom_vline(xintercept=mean(dat[year == 1950, log_gdp]), linetype="dashed", 
#                 color = "red", size=1))

# print(ggplot(dat, aes(x = log_gdp, y = wy_mean,color = gyear)) +
#     geom_point() + geom_smooth(method=lm, se = TRUE) + 
#     geom_vline(xintercept=mean(dat[year == 1950, log_gdp]), linetype="dashed", 
#                 color = "red", size=0.5) + 
#     geom_vline(xintercept=mean(dat[year == 1970, log_gdp]), linetype="dashed", 
#                 color = "blue", size=0.5) + 
#     geom_vline(xintercept=mean(dat[year == 1990, log_gdp]), linetype="dashed", 
#                 color = "green", size=0.5)
# )
# dev.off()

multiResultClass = function(models = NULL, shifts_baseline = NULL, 
    shifts_stacking = NULL) {
    me = list(models = models, shifts_baseline = shifts_baseline, 
    shifts_stacking = shifts_stacking)
    class(me) = append(class(me), "multiResultClass")
    return(me)
}

# loop over imputations
output = foreach(i = 1:nimputations) %dopar% {

    results = multiResultClass()
    library(data.table)
    library(brms)
    library(loo)
    source("src/utils.R")
    
    models = list()
    dat = timps[[i]]

    bprior = c(prior(normal(0,5), class = b), 
        prior(normal(0,5), class = Intercept))
    models[[1]] = brm(wy_mean ~ 1 + log_gdp * year1950, data = dat, 
        iter = iterations[[select_estimates]][[1]], 
        chains = 1,
        refresh = 5000,
        warmup = 1000, 
        prior = bprior)

    bprior = c(prior(normal(0,5), class = b),
        prior(normal(0,5), class = Intercept),
        prior(cauchy(0,2.5), class = sd, group = qyear))
   
    models[[1]] = brm(wy_mean ~ 1 + log_gdp + (log_gdp | qyear), data = dat, 
        iter = iterations[[select_estimates]][[1]], 
        chains = 1,
        refresh = 5000,
        warmup = 1000, 
        prior = bprior)

    # shifts
    pred = predict(models[[1]], newdata = dat, 
        summary = TRUE)

    countries = unique(dat$ctry)
    vars = c("ctry", "ctry50", "ctryear", "zyear", "year1950", "qyear")
    plots_checks = prediction_check_plots(pred, dat, countries, 
        vars,  y = "ex_mean", x = "year", 
        xlab = "\nYear", ylab = "Average life expectancy at age 0\n", 
        transform = TRUE, 
        country_labels = country_labs)

    savepdf(paste0(plots_path, select_estimates, "fit_no_error_yearly"))
        print(plots_checks)
    dev.off()

    countries = unique(covs$ctry)
    years = c(1950)
    cyears = list(c("1950", "1950-1969"))
    dyears = list(c("1950", "1950+"))

    newdata = createComparisonData(covs, countries, 
        years, cyears, dyears)
    countries = unique(newdata$ctry)

    shifts = createShifts(models[[1]], newdata, 
        countries = countries, nsamples = 5000, K = K) 
    
    summaryFun = function(x) {
        list(
            "Mean" = mean(x), 
            "Q2.5" = quantile(x, 0.025), 
            "Q97.5" = quantile(x, 0.975)
        )
    }
    tab1 = t(sapply(shifts$values, summaryFun))
    tab1

    bprior = c(prior(normal(0,5), class = b),
        prior(normal(0,5), class = Intercept),
        prior(cauchy(0,2.5), class = sd, group = ctry50))

    models[[2]] = brm(wy_mean ~ log_gdp + (1 | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[2]], 
        warmup = 1000, 
        chains = 1,
        refresh = 5000,
        prior = bprior)

    models[[3]] = brm(wy_mean ~ log_gdp + (log_gdp | ctry50), data = dat, 
        iter = iterations[[select_estimates]][[3]], 
        warmup = 1000, 
        chains = 1, 
        refresh = 5000,
        prior = bprior)

    # models[[4]] = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctry50), data = dat, 
    #     iter = iterations[[select_estimates]][[4]], 
    #     warmup = 1000, 
    #     chains = 1,
    #     prior = bprior)

    bprior = c(prior(normal(0,5), class = b),
        prior(normal(0,5), class = Intercept),
        prior(cauchy(0,2.5), class = sd, group = ctryear))

    models[[4]] = brm(wy_mean ~ log_gdp + (1 | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[4]], 
        warmup = 1000, 
        chains = 1, 
        refresh = 5000,
        prior = bprior)

    models[[5]] = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[5]], 
        warmup = 1000, 
        chains = 1, 
        refresh = 5000,
        prior = bprior)


    models[[5]] = brm(wy_mean ~ log_gdp + lit + (log_gdp | ctryear), data = dat, 
        iter = iterations[[select_estimates]][[5]], 
        warmup = 1000, 
        chains = 1, 
        refresh = 5000,
        prior = bprior)

    summary(models[[5]])
    pred = predict(models[[5]], newdata = dat, 
        summary = TRUE)

    countries = unique(dat$ctry)
    vars = c("ctry", "ctry50", "ctryear", "zyear", "year1950", "qyear")
    plots_checks = prediction_check_plots(pred, dat, countries, 
        vars,  y = "ex_mean", x = "year", 
        xlab = "\nYear", ylab = "Average life expectancy at age 0\n", 
        transform = TRUE, 
        country_labels = country_labs)

    savepdf(paste0(plots_path, select_estimates, "fit_no_error_preferred"))
        print(plots_checks)
    dev.off()

    countries = unique(covs$ctry)
    years = c(1950)
    cyears = list(c("1950", "1950-1969"))
    dyears = list(c("1950", "1950+"))

    newdata = createComparisonData(covs, countries, 
        years, cyears, dyears)
    countries = unique(newdata$ctry)

    shifts = createShifts(models[[5]], newdata, 
        countries = countries, nsamples = 5000, K = K) 
    
    tab3 = t(sapply(shifts$values, summaryFun))
    tab3
    # models[[6]] = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctryear), data = dat, 
    #     iter = iterations[[select_estimates]][[6]], 
    #     warmup = 1000, 
    #     chains = 1, 
    #     refresh = 5000,
    #     prior = bprior)

    results$models = models

    countries = unique(covs$ctry)
    years = c(1950)
    cyears = list(c("1950", "1950-1969"))
    dyears = list(c("1950", "1950+"))

    newdata = createComparisonData(covs, countries, 
        years, cyears, dyears)
    countries = unique(newdata$ctry)

    results$shifts_stacking = createShifts(models, newdata, 
        countries = countries, nsamples = 5000, K = K, stacking = stacking) 

    results$shifts_baseline = createShifts(models[[1]], newdata, 
        countries = countries, nsamples = 5000) 

    return(results)
}

# stop clusters
stopCluster(cl)

# extract results
shifts = list()
weights = list()
nmodels = length(output[[1]]$models)
for (i in seq_along(output)) {
   shifts[[i]] = output[[i]][["shifts_stacking"]][["values"]]
   weights[[i]] = output[[i]][["shifts_stacking"]][["weights"]]
}

shifts = rbindlist(shifts, idcol = "replicate")
avg_weights = apply(do.call(rbind, weights), 2, mean)
saveRDS(list("shifts" = shifts, "avg_weights" = avg_weights),
    paste0(data_path, select_estimates, "shifts_stacking_no_error.rds"))

# shifts = list()
# for (i in seq_along(output)) {
#    shifts[[i]] = output[[i]][["shifts_baseline"]][["values"]]
# }
# shifts = rbindlist(shifts, idcol = "replicate")
# saveRDS(shifts,
#     paste0(data_path, select_estimates, "shifts_baseline_no_error.rds"))
rm(shifts)
slackr::slackr_msg(txt = paste0("LE no error saved shifts: ", Sys.time()))

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

# create predicitive checks 
countries = unique(idat$ctry)
preferred_model = which.max(avg_weights)
vars = c("ctry", "ctry50", "ctryear", "zyear", "year1950")

# stacking
pred = rlang::invoke(pp_average, model_list, weights = avg_weights, 
    newdata = idat, summary = TRUE, nsamples = 5000)

plots_checks = prediction_check_plots(pred, idat, countries, vars, 
    y = "ex_mean", x = "year", 
    xlab = "\nYear", ylab = "Average life expectancy at age 0\n",
    transform = TRUE, 
    country_labels = country_labs)
savepdf(paste0(plots_path, select_estimates, "fit_no_error_stacking"))
    print(plots_checks)
dev.off()
file.copy(paste0(plots_path, select_estimates, "fit_no_error_stacking.pdf"), 
    manus_plots, recursive = TRUE)    
rm(pred, plots_checks)

# # baseline model
# pred = predict(model_list[[1]], newdata = idat, summary = TRUE)
# plots_checks = prediction_check_plots(pred, idat, countries, 
#     vars,  y = "ex_mean", x = "year", 
#     xlab = "\nYear", ylab = "Average life expectancy at age 0\n",
#     transform = TRUE, 
#     country_labels = country_labs)
# savepdf(paste0(plots_path, select_estimates, "fit_no_error_baseline_model"))
#     print(plots_checks)
# dev.off()
# file.copy(paste0(plots_path, select_estimates, "fit_no_error_baseline_model.pdf"), 
#     manus_plots, recursive = TRUE)    
# rm(pred, plots_checks)

# preferred model 
pred = predict(model_list[[preferred_model]], newdata = idat, 
    summary = TRUE, nsamples = 10000)
plots_checks = prediction_check_plots(pred, idat, countries, 
    vars,  y = "ex_mean", x = "year", 
    xlab = "\nYear", ylab = "Average life expectancy at age 0\n", 
    transform = TRUE, 
    country_labels = country_labs)
savepdf(paste0(plots_path, select_estimates, "fit_no_error_preferred_model"))
    print(plots_checks)
dev.off()
file.copy(paste0(plots_path, select_estimates, "fit_no_error_preferred_model.pdf"), 
    manus_plots, recursive = TRUE)    
rm(pred, plots_checks)

# table
tabs = list()
for (i in seq_along(model_list)) {
    print(paste0("Extracting model ", i))
    tabs[[i]] = extractBRMS(model_list[[1]], r2 = TRUE, weight = avg_weights[i])
    model_list[[1]] = NULL
}
saveRDS(tabs, paste0(data_path, select_estimates, "tab_no_error.rds"))

cnames = paste0("Model ", 1:nmodels)
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", year19501950P = "1950+", 
    "log_gdp:year19501950P" = "Log GDP x 1950+")
# "zyear" = "Year (standardized)"

caption = paste0("Models for LE and log GPD no measurement error, ", 
    nimputations, " imputations")

texreg(tabs,
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = paste0("tab:", select_estimates,"ex_no_error"),
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = paste0(tables_path, select_estimates, "models_no_error.tex")
)
file.copy(paste0(tables_path, select_estimates, "models_no_error.tex"), manus_tables, 
        recursive = TRUE)

# send message to slack
slackr::slackr_msg(txt = paste0("LE models no error finished at: ", Sys.time()))
