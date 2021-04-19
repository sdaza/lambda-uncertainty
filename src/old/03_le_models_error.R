
# life expectancy error models
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
library(slackr)

plan(multicore)

library(texreg)
library(mice)
library(miceadds)

library(ggplot2)
library(patchwork)
source("src/utils.R")

options(mc.cores = 30, scipen=999)
seed = 19380302
set.seed(seed)

estimate_models = FALSE 

idat = readRDS("output/data/single-imputation.rds")
sdat = readRDS("output/data/aggregate-data.rds")
imps = readRDS("output/data/imputations.rds")

if (estimate_models) {

    e1 = brm_multiple(ex ~ log_gdp + (1|ctryear),
        family = gaussian, data = imps,
        iter = 10000, 
        warmup = 1000, 
        chains = 5, 
        seed = seed,
        refresh = 11000, 
        cores = 10
    )
    
    max(as.vector(round(e1$rhats, 2)))
    pred = predict(e1)
    test = loo(e1)
    
    e2 = brm_multiple(ex ~ log_gdp + (log_gdp | ctryear),
        family = gaussian, data = timps,
        iter = 10000, 
        warmup = 1000, 
        chains = 5, 
        seed = seed,
        refresh = 11000, 
        cores = 10
    )
    
    max(as.vector(round(e2$rhats, 2)))
    summary(e2)

    e3 = brm_multiple(ex ~ log_gdp + zyear + (1|ctry),
        family = gaussian, data = timps,
        iter = 11000, 
        warmup = 1000, 
        chains = 10, 
        seed = seed,
        refresh = 11000, 
        cores = 10, 
    )

    max(as.vector(round(e3$rhats, 2)))
    summary(e3)

    e4 = brm_multiple(ex ~ log_gdp + zyear + (zyear + log_gdp|ctry),
        family = gaussian, data = timps,
        iter = 15000, 
        warmup = 1000, 
        chains = 10, 
        seed = seed,
        refresh = 11000, 
        cores = 10
    )

    max(as.vector(round(e4$rhats, 2)))
    baseline_models_e = list(e1, e2, e3, e4)
    saveRDS(baseline_models_e, file = "output/models/baseline_models_e.rds")
} else { baseline_models_e = readRDS("output/models/baseline_models_e.rds") }

# if (estimate_mdoels) {
    # loo_list_e = list()
    # for (i in seq_along(baseline_models_e)) {
    #        loo_list_e[[i]] = loo(baseline_models_e[[i]], reloo=TRUE)
    # }
    # saveRDS(loo_list_e, file = "output/models/loo_baseline_models_e.rds")
# } else { loo_list_e = readRDS("output/models/loo_baseline_models_e.rds") }


# model_weights_e = as.vector(loo_model_weights(loo_list_e))
# round(model_weights_e, 3)

# cnames = paste0("Model ", 1:4,  " (", round(model_weights_e, 4), ")")
cnames = paste0("Model ", 1:4)
custom_coeff_map = list(Intercept = "Constant", "log_gdp" = "Log GDP", "zyear" = "Year (standardized)")

texreg(baseline_models_e,
    caption = "Models for LE and log GPD (measurement error)", 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:ex_error",
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    include.loo.ic = FALSE,
    inclue.rsquared = FALSE,
    include.waic = FALSE,
    include.random = FALSE,
    reloo = FALSE,
    file = "output/tables/models_error.tex"
)
file.copy("output/tables/models_error.tex", "manuscript/tables/", recursive = TRUE)

# checking fit
ctrys = unique(sdat$ctry)
vars = c("ctry", "year", "log_gdp", "ex_mean")

test = predict(baseline_models_e[[2]])


plots_checks = prediction_checks_ex(e2, sdat, ctrys, vars, y = "ex_mean", x = "year")

plots_checks[[1]]

savepdf("output/plots/fit_error_m2")
    print(plots_checks)
dev.off()
file.copy("output/plots/fit_error_m2.pdf", "manuscript/plots/", recursive = TRUE)

# send message to slack
slackr::slackr_setup(config_file = ".slackr")
slackr::slackr_msg(paste0("LE models no error: ", Sys.time()))

