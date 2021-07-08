# life expectancy models with error
# author: sebastian daza
#####################################
#####################################


# R < src/03_le_models_err.R > output/log/03_le_models_err_t.log  --no-save  &
# R < src/03_le_models_err.R > output/log/03_le_models_err_f.log  --no-save  &
# install.packages("StanHeaders", repos="https://cloud.r-project.org")
# install.packages("rstan", repos="https://cloud.r-project.org")


# libraries, functions and options
library(data.table)
library(doParallel)
library(foreach)
library(texreg)
# library(rethinking)
# library(bayesplot)
library(ggplot2)
# library(ggridges)
# library(ggcorrplot)
# library(patchwork)
# rstan_options(auto_write = FALSE)

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
samples = data_list[["imputations"]]
idat = data_list[["single-imputation"]]
country_labs = data_list[["ctrylabels"]]
covs = data_list[["covs"]]
ex_max = data_list[["ex_max"]]
lctryear = levels(covs$ctryear)
formulas = data_list[["formulas"]]
newdata = data_list[["newdata"]]

select_estimates = ""

# replicates
nsamples = length(samples)

# create data object for testing
# idat = samples[[sample(1:10, 1)]]
# idat[, cy := .GRP, ctryear]
# mdata = list(
#     wy = idat$wy, 
#     zyear = idat$zyear,
#     zinfrastructure = idat$zinfrastructure, 
#     zpop= idat$zpop, 
#     zilit = idat$zilit, 
#     zius_aid_pc  = idat$zius_aid_pc, 
#     zigdp_pc = idat$zigdp_pc, 
#     ctryearg = idat$ctryearg)

# model = rethinking::ulam(formulas[[1]],
#   data = mdata, chains = 4, cores = 4, iter = 4000
# )

# rethinking::precis(model, depth = 1)


# table list
# tabs = list()
# tshifts = list()
iterations = c(4000, 2000, 3000, 3000, 3000)

# get model output
model_number = 4
output = runModel(formulas[[model_number]], samples, newdata = newdata, ex_max = ex_max, 
    iterations = iterations[model_number], clusters = 3)

print(paste0("Number of parameters with Rhat4 > 1.01: ", output$rhat))
print(paste0("Number of parameters with neff < 100: ", output$rhat))

# prediction plot
plots = predictionPlots(output$predictions, country_labels = country_labs)
plots = patchwork::wrap_plots(plots, ncol = 3)
savepdf(paste0(plots_path, select_estimates, 
    paste0("fit_check_m", model_number)), width = 30, height = 35)
    print(plots)
dev.off()
file.copy(paste0(plots_path, select_estimates, 
    paste0("fit_check_m", model_number, ".pdf")), manus_plots, recursive = TRUE)

# shift plot
shift_plot = plotShifts(output$shifts, country_labs)
savepdf(paste0(plots_path, select_estimates, 
    paste0("shifts_by_period_m", model_number)), height = 20)
    print(shift_plot)
dev.off()
file.copy(paste0(plots_path, select_estimates, 
    paste0("shifts_by_period_m", model_number, ".pdf")), manus_plots, recursive = TRUE)

# shifts
tabshift = output$shifts
tabshift = na.omit(tabshift, "shift")
v = unlist(country_labs)
tabshift = tabshift[, .(estimate = paste0(specify_decimal(mean(shift), 2), 
    " [", specify_decimal(quantile(shift, probs = 0.025), 2), ", ", 
    specify_decimal(quantile(shift, probs = 0.975), 2), "]")), 
    .(ctry, year)]
tabshift[, lctry := as.factor(v[as.character(ctry)])]
tabshift[, model := model_number]
tshifts = readRDS(paste0(tables_path, "tab_shifts.rds"))
tshifts[[model_number]] = tabshift
saveRDS(tshifts, paste0(tables_path, "tab_shifts.rds"))

createShiftTable(tshifts[[model_number]], paste0(tables_path, "shifts_m", model_number, ".tex"))
file.copy(paste0(tables_path, select_estimates, "shifts_m", model_number, ".tex"), manus_tables, 
    recursive = TRUE)    
    
# model table
tabs = readRDS(paste0(tables_path, "tabs.rds"))
tabs[[model_number]] = extractStan(output$fit, 
    n = nrow(idat), r2 = output$r2)
rm(output)
saveRDS(tabs, paste0(tables_path, "tabs.rds"))
slackr::slackr_msg(txt = paste0("Stan loop finish at: ", Sys.time()))


# create regression table
select_estimates = ""
tabs = readRDS(paste0(tables_path, "tabs.rds"))

# rename some coefficients
cnames = tabs[[1]]@coef.names
cnames[grepl("sigma\\_cy", cnames)] = "sigma_cy[1]"
tabs[[1]]@coef.names = cnames

cnames = paste0("Model ", 1:5)
custom_coeff_map = list(
    "a" = "Constant", 
    "b_gdp" = "GDP per capita", 
    "b_year" = "Year", 
    "b_lit" = "Literacy", 
    "b_infra" = "Infrastructure",
    "b_pop" = "Population",
    "b_us" = "US aid per capita",
    "sigma_cy[1]" = "Intercept", 
    "sigma_cy[2]"  = "GDP", 
    "Rho[1,2]" = "cor(Intercept, GDP)"
)

caption = paste0("Models for LE and GPD ", 
    100, " replicates")

groups = list("Fixed effects" = 1:7, "Random effects" = 8:10)

texreg::texreg(tabs, 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = paste0("tab:", select_estimates, "models"),
    groups = groups,
    custom.note = "\\item  $^*$ Null hypothesis value outside the confidence interval. 
        \\item All covariates are standardized. Life expectancy estimates were transformed using $ln\\left(-ln( 1-\\frac{e_0}{ (78.6 + 1.05)}\\right)$.",
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = paste0(tables_path, select_estimates, "models.tex")
)    
file.copy(paste0(tables_path, select_estimates, "models.tex"), manus_tables, 
    recursive = TRUE)    

screenreg(tabs, 
  omit.coef = "^a\\_cy\\[.+|^b_gdp\\_cy\\[.+|^Rho\\[1,1\\]|^Rho\\[2,2\\]" )

# testing parallel work windows
library(rstan)
options(mc.cores = 1)
library(doParallel)
cl = makeCluster(4)
registerDoParallel(cl)
output = foreach(i = 1:20,
        .packages = c("rstan", "rethinking")) %dopar% {

    options(mc.cores = 1)
    dat = samples[[i]]
    mdata = list(
                wy = dat$wy, 
                zyear = dat$zyear,
                zinfrastructure = dat$zinfrastructure, 
                zpop= dat$zpop, 
                zilit = dat$zilit, 
                zius_aid_pc = dat$zius_aid_pc, 
                zigdp_pc = dat$zigdp_pc, 
                ctryearg = dat$ctryearg
    )
    
        model = ulam(
            formulas[[1]], 
            data = mdata, chains = 1, cores = 1, 
            iter = 4000
        )
    return(model)
    Sys.sleep(samples(20:40, 1))
}

stopCluster(cl)
