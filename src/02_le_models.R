############################
# life expectancy models
# author: sebastian daza
############################

# R < /home/s/sdaza/00projects/lambda/src/le_models.R > /home/s/sdaza/00projects/lambda/output/le_models.log  --no-save  &
# R < src/02_le_models.R > output/log/02_le_models.log  --no-save  &


# libraries, functions and options
library(haven)
library(data.table)
library(stringr)
library(brms)
library(loo)
library(future)
plan("multicore")

library(texreg)
library(ggplot2)
library(patchwork)
source("src/utils.R")

# options(mc.cores = 30, scipen=999)
# options(mc.cores = 10)
seed = 19380302
set.seed(seed)

estimate_models = TRUE

# read data
sdat = readRDS("output/data/aggregate-data.rds")
idat = readRDS("output/data/single-imputation.rds")
timps = readRDS("output/data/imputations.rds")
nsamples = length(timps)

idat[, ctry50 := paste0(ctry, ".", year1950)]
idat[, xy_sd := transWeibull(ex_sd, max_ex)]

head(idat)

# models 
m1 = brm(wy_mean ~ log_gdp + (1 | ctry50), data = idat, 
    iter = 15000, 
    refresh = 11000, 
    warmup = 1000, 
    chains = 6, 
    cores = 6, 
    save_pars = save_pars(all = TRUE))

m2 = brm(wy_mean ~ log_gdp + (log_gdp | ctry50), data = idat, 
    iter = 15000, 
    refresh = 11000, 
    warmup = 1000, 
    chains = 6, 
    cores = 6, 
    save_pars = save_pars(all = TRUE))

m3 = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctry50), data = idat, 
    iter = 15000, 
    refresh = 11000, 
    warmup = 1000, 
    chains = 6, 
    cores = 6, 
    save_pars = save_pars(all = TRUE), 
    control = list(max_treedepth = 15))

m4 = brm(wy_mean ~ log_gdp + (1 | ctryear), data = idat, 
    iter = 20000, 
    refresh = 11000, 
    warmup = 1000, 
    chains = 6, 
    cores = 6, 
    save_pars = save_pars(all = TRUE))

m5 = brm(wy_mean ~ log_gdp + (log_gdp | ctryear), data = idat, 
    iter = 20000, 
    refresh = 11000, 
    warmup = 1000, 
    chains = 6, 
    cores = 6, 
    save_pars = save_pars(all = TRUE), 
    control = list(max_treedepth = 15))

m6 = brm(wy_mean ~ log_gdp + zyear + (log_gdp | ctryear), data = idat, 
    iter = 15000, 
    refresh = 11000, 
    warmup = 1000, 
    chains = 6, 
    cores = 6, 
    save_pars = save_pars(all = TRUE), 
    control = list(max_treedepth = 15))

model_list = list(m1, m2, m3, m4, m5, m6)
loo_list =list()
kfold_list = list()
for (i in 1:6) {
    # loo_list[[i]] = loo(model_list[[i]], cores = 6, moment_match = TRUE)
    kfold_list[[i]] = kfold(model_list[[i]], K = 10, chains  = 1)
}

# lpd_point = NULL
lpd_kfold = NULL 
for (i in c(1, 2, 3, 4, 5, 6)) {
    # lpd_point = cbind(lpd_point, loo_list[[i]]$pointwise[, "elpd_loo"])
    lpd_kfold = cbind(lpd_kfold, kfold_list[[i]]$pointwise[, "elpd_kfold"])
}

# pbma_wts = pseudobma_weights(lpd_point, BB=FALSE)
# pbma_BB_wts = pseudobma_weights(lpd_point)
# stacking_wts = stacking_weights(lpd_point)
# round(cbind(pbma_wts, pbma_BB_wts, stacking_wts), 2)

pbma_wts = pseudobma_weights(lpd_kfold, BB=FALSE)
pbma_BB_wts = pseudobma_weights(lpd_kfold)
stacking_wts = stacking_weights(lpd_kfold)
round(cbind(pbma_wts, pbma_BB_wts, stacking_wts), 2)

ppred = pp_average(m1, m2, m3, m4, m5, m6, weights = stacking_wts)

savepdf("output/plots/pred_check_m3")
    print(pp_check(m6))
dev.off()

ctrys = unique(idat$ctry)
vars = c("ctryear", "ctry50", "ctry", "year", "zyear", "log_gdp", "ex_mean", "wy_mean")
plots_checks = prediction_checks_pp_ex(ppred, idat, ctrys, vars, y = "ex_mean", x = "year", 
    transform = TRUE)
savepdf("output/plots/fit_no_error_pp")
    print(plots_checks)
dev.off()
file.copy("output/plots/fit_no_error_pp.pdf", "manuscript/plots/", recursive = TRUE)    

# preferred model 
preferred_model = which.max(stacking_wts)
plots_checks = prediction_checks_ex(model_list[[preferred_model]], idat, ctrys, vars, y = "ex_mean", x = "year", transform = TRUE)

savepdf(paste0("output/plots/fit_no_error_m", preferred_model))
    print(plots_checks)
dev.off()
file.copy(paste0("output/plots/fit_no_error_m", preferred_model, ".pdf"), "manuscript/plots/", recursive = TRUE)   

# gpd
savepdf("output/plots/imputation_check_gdp")
for (i in countries) {
    print(ggplot(data = idat[ctry == i, .(log_gdp, year, gdp_missing)],
    aes(year, log_gdp, color = gdp_missing)) +
    geom_point() + labs(title = i, x = "Year", y = "Log GDP", color = NULL))        
}
dev.off()
file.copy("output/plots/imputation_check_gdp.pdf", "manuscript/plots", recursive = TRUE)

# table
cnames = paste0("Model ", 1:length(model_list),  " (", round(stacking_wts, 2), ")")
cnames
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", "zyear" = "Year (standardized)")

texreg(baseline_models,
    caption = "Models for LE and log GPD no measurement error, stacking weight in parenthesis", 
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
    file = "output/tables/models_no_error.tex"
)    
file.copy("output/tables/models_no_error.tex", "manuscript/tables/", recursive = TRUE)    

# shift 
countries = unique(idat$ctry)
years = c(1950)
cyears = list(c("1950", "1950-1969"))
dyears = list(c("1950", "1950+"))
datatest = createComparisonData(idat, countries, years, cyears, dyears)

shift = pp_average(m1, m2, m3, m4, m5, m6, weights = stacking_wts, newdata = datatest, summary = FALSE)
shift = data.table(shift)

savepdf("")
vars = names(shift)
shift[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]
datatest

ctry = unique(datatest$ctry)

savepdf("output/plots/shifts_no_error_pp")
    index = c(1, 2)    
    plots ()
    for (i in seq_along(ctry)) {
            out = data.table(est = shift[[paste0("V", index[1])]]  - shift[[paste0("V", index[2])]])
            print(ggplot(out, aes(x = est)) +
                geom_histogram(color = "black", fill = "white", binwidth = 0.5) + 
                theme_minimal() + labs(x = "Estimated shift", y = "Frequency", title = ctry[i])
            )
            index = index + 2
    }
dev.off()
file.copy("output/plots/shifts_no_error_pp.pdf", "manuscript/plots", recursive = TRUE)

# send message to slack
slackr::slackr_setup(config_file = ".slackr")
slackr::slackr_msg(txt = paste0("LE models no error: ", Sys.time()))