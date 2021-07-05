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
library(stringr)
library(ggplot2)
library(ggridges)
library(ggcorrplot)
library(patchwork)
library(doParallel)
library(texreg)
library(bayesplot)
library(Rcpp)
library(rethinking)
rstan_options(auto_write = TRUE)


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

# replicates
nsamples = length(samples)
# nsamples = 20˜
print(paste0("Number of replicates: ", nsamples))
samples = samples[1:nsamples]

# correlation between indicators
vars = c("zpop", "zyear", "zigdp_pc", "ziurban", "zielec", "zilit", "ziwater",
     "zisewage", "zius_aid_pc", "zinfrastructure")
r = cor(idat[, ..vars])
savepdf(paste0(plots_path, select_estimates, "correlation_imputed_covs"))
ggcorrplot(r, 
    hc.order = TRUE, 
    type = "lower",
    lab = TRUE)
dev.off()
file.copy(paste0(plots_path, select_estimates, "correlation_imputed_covs.pdf"), 
    manus_plots, recursive = TRUE)    

# data for shifts
vars = c(vars, c("ctry", "ctry50", "ctryear", "zyear", "year1950"))
countries = unique(covs$ctry)
years = c(1950, 1970, 1990)
cyears = list(c("1950", "1950-1969"), c("1950-1969", "1970-1989"), 
    c("1970-1989", "1990"))
dyears = list(c("1950", "1950+"), c("1950", "1950+"), 
    c("1950", "1950+"))
newdata = createComparisonData(covs, countries, 
    years, cyears, dyears)
newdata[, ctryear := factor(ctryear, levels = lctryear)]
newdata[, ctryearg := as.numeric(ctryear)]

# # create data object for testing
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

# m1 = ulam(
#   alist(
#     wy ~ normal(mu, sigma),
#     mu <- a_cy[ctryearg] + b_gdp * zigdp_pc,
#     a_cy[ctryearg] ~ normal(a, sigma_cy),
#     a ~ normal(0, 1),
#     # c(sigma_cy, sigma) ~ half_normal(0,1),
#     c(sigma_cy, sigma) ~ exponential(1),
#     b_gdp ~ normal(0, 1), 
#     save> pred <- a_cy[ctryearg] + b_gdp * zigdp_pc
#   ),
#   data = mdata, chains = 1, cores = 1, iter = 2000
# )

select_estimates = ""

# model 1
model_number = 1
flist = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp * zigdp_pc,
    a_cy[ctryearg] ~ normal(a, sigma_cy),
    a ~ normal(0, 1),
    # c(sigma_cy, sigma) ~ half_normal(0,1),
    c(sigma_cy, sigma) ~ exponential(1),
    b_gdp ~ normal(0, 0.5), 
    save> pred <- a_cy[ctryearg] + b_gdp * zigdp_pc
)

tabs = list()
models = list()
predictions = list()
shifts = list()
r2 = NULL
for (i in seq_along(samples)) {
    print(paste0(":::::::: Running iteration ", i, " ::::::::"))
      dat = samples[[i]]
      dat[, cy := as.numeric(ctryear)]
      mdata = list(
          wy = dat$wy,
          zigdp_pc = dat$zigdp_pc, 
          ctryearg = dat$ctryearg
      )
      print(":::::::: running stan models")
      iterations = 4000
      model = ulam(
          flist, 
          data = mdata, chains = 1, cores = 1, 
          iter = iterations, chain_id = i
      )
      check = as.matrix(precis(model, depth = 3))
      print(paste0(":::::::: There are ", 
          sum(check[, "Rhat4"] > 1.01), 
          " parameters with Rhat4 > 1.01"))
      print(paste0(":::::::: There are ", 
          sum(check[, "n_eff"] < 100), 
          " parameters with n_eff < 100"))
      rm(check)
      print(":::::::: saving output")
      models[[i]] = model@stanfit
      predictions[[i]] = predictData(model, idat, ex_max, n = 100)
      shifts[[i]] = computeShift(model, newdata, ex_max, n = 100)
      fit_ss = rstan::extract(models[[i]], pars = c("pred", "sigma"))
      r2 =c(r2, bayes_R2(fit_ss$pred, fit_ss$sigma))
      rm(fit_ss, model)
}

shifts = rbindlist(shifts)
pred = rbindlist(predictions)
fit = sflist2stanfit(models)
tabs[[model_number]] = extractStan(fit, n = nrow(idat), r2 = median(r2))

plots = predictionPlots(pred, country_labels = country_labs)
plots = wrap_plots(plots, ncol = 3)
savepdf(paste0(plots_path, select_estimates, 
    paste0("fit_check_m", model_number)), width = 30, height = 35)
    print(plots)
dev.off()

shift_plot = plotShifts(shifts, country_labs)
savepdf(paste0(plots_path, select_estimates, 
    paste0("shifts_by_period_m", model_number)), height = 20)
    print(shift_plot)
dev.off()



cl = makeCluster(5)
registerDoParallel(cl)
test1 = runModel(flist, samples, iterations = 2500, 
    model_number = 1, ex_max = ex_max, country_labs = country_labs, 
    newdata = newdata, plots_path = plots_path, 
    select_estimates = "")


stopCluster(cl)
# model 2
flist = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 2),
    b_gdp ~ normal(0, 1),
    c(sigma, sigma_cy) ~ exponential(1),
    Rho ~ lkj_corr(2),
    save> pred <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc
)

test2 = runModel(flist, samples[1:3], clusters = 2, 
    iterations = 2000, chains = 1,
    model_number = 2, ex_max = ex_max, country_labs = country_labs, 
    newdata = newdata, plots_path = plots_path, 
    select_estimates = "")



# model 3
flist = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_year * zyear,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 2),
    c(b_gdp, b_year) ~ normal(0, 1),
    c(sigma, sigma_cy) ~ exponential(1),
    Rho ~ lkj_corr(2),
    save> pred <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_year * zyear
)

test3 = runModel(flist, samples[1:3], clusters = 2, iterations = 3000, 
    model_number = 3, ex_max = ex_max, country_labs = country_labs, 
    newdata = newdata, plots_path = plots_path, 
    select_estimates = "")

screenreg(list(test[[1]], test2[[1]]), 
  omit.coef = "a\\_cy\\[.+|gdp\\_cy\\[.+")

# baseline model
cl = makeCluster(2)
# cl = makeCluster(10, outfile="")
# redefine samples 
nsamples = 10
registerDoParallel(cl)
output = foreach(i = 1:nsamples) %dopar% {
        library(data.table)
        library(rethinking)
        results = multiResultClass()

        print(paste0(":::::::: Running iteration ", i, " ::::::::"))
        dat = samples[[i]]
        dat[, cy := as.numeric(ctryear)]
        mdata = list(
            wy = dat$wy, 
            zyear = dat$zyear,
            zinfrastructure = dat$zinfrastructure, 
            zpop= dat$zpop, 
            zilit = dat$zilit, 
            zius_aid_pc  = dat$zius_aid_pc, 
            zigdp_pc = dat$zigdp_pc, 
            ctryearg = dat$ctryearg
        )
        model = ulam(
            alist(
                wy ~ normal(mu, sigma),
                mu <- a_cy[ctryearg] + b_gdp * zigdp_pc,
                a_cy[ctryearg] ~ normal(a, sigma_cy),
                a ~ normal(0, 1),
                # c(sigma_cy, sigma) ~ half_normal(0,1),
                c(sigma_cy, sigma) ~ exponential(1),
                b_gdp ~ normal(0, 1), 
                save> pred <- a_cy[ctryearg] + b_gdp * zigdp_pc
            ), 
            data = mdata, chains = 1, cores = 1, 
            iter = 3000, chain_id = i
        )
        
        results$models = model@stanfit
        results$shifts = computeShift(model, newdata, ex_max)
        return(results)
}
stopCluster(cl)

# extract restuls
shifts = list()
models = list()
for (i in seq_along(output)) {
   shifts[[i]] = output[[i]][["shifts"]]
   models[[i]] = output[[i]][["models"]]
}

shifts = rbindlist(shifts)
fit = sflist2stanfit(models)
t = extractStan(fit, n = nrow(idat))
screenreg(t, omit.coef = "a\\_cy\\[.+")

pred = predictionData(fit, idat, ex_max = ex_max)
plots = predictionPlots(pred, country_labels = country_labs)
plots = wrap_plots(plots, ncol = 3)
model_number = 1
savepdf(paste0(plots_path, select_estimates, 
    paste0("fit_check_m", model_number)), width = 30, height = 35)
    print(plots)
dev.off()

xmin = min(shifts$shift, na.rm = TRUE) - 1.5
xmax = max(shifts$shift, na.rm = TRUE) + 1.5
v = unlist(country_labs)
shifts[, lctry := as.factor(v[as.character(ctry)])]
shifts[, year := factor(year)]
shifts[, ctry := factor(ctry)]
savepdf(paste0(plots_path, select_estimates, 
    paste0("shifts_by_period_m", model_number)), height = 20)
    print(ggplot(shifts, aes(y = lctry)) +
    geom_density_ridges(aes(x = shift, fill = year), 
        alpha = .45, color = "white", from = xmin, to = xmax, scale = 1) +
    labs(x = "Shift", y = "Country", subtitle = "", caption = "") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_cyclical(
        values = c("#EC7063", "#F7DC6F", "#229954"), guide = "legend") +
    coord_cartesian(clip = "off") +
    theme_ridges(grid = TRUE) +
    theme(legend.position="top", 
    legend.title = element_blank())
    )
dev.off()



m2 = ulam(
  alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 2),
    b_gdp ~ normal(0, 1),
    c(sigma, sigma_cy) ~ exponential(1),
    Rho ~ lkj_corr(2)
  ),
  data = mdata, chains = 1, cores = 1, iter = 2000
)

precis(m2, depth = 2)
test = m2@stanfit
etest = extract(test)
str(etest)
m3 = ulam(
  alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_cy[cy] * zigdp_pc + byear * zyear,
    c(a_cy, b_cy)[cy] ~ multi_normal(c(a, b), Rho, sigma_cy),
    a ~ normal(0, 2),
    c(b, byear) ~ normal(0, 1),
    c(sigma_cy, sigma) ~ exponential(1),
    Rho ~ lkj_corr(2), 
    save> pred <- a_cy[cy] + b_gdp_cy[cy] * zigdp_pc
  ),
  data = mdata, chains = 1, cores = 1, iter = 2000
)

m3 = ulam(
  alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[cy] + b_cy[cy] * zigdp_pc + byear[cy] * zyear,
    c(a_cy, b_cy, byear)[cy] ~ multi_normal(c(a, b, c), Rho, sigma_cy),
    a ~ normal(0, 2),
    c(b, c, byear) ~ normal(0, 1),
    c(sigma_cy, sigma) ~ exponential(1),
    Rho ~ lkj_corr(2)
  ),
  data = mdata, chains = 1, cores = 1, iter = 4000
)
precis(m3, depth = 2)

m4 = ulam(
  alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[cy] + b_cy[cy] * zigdp_pc + b_infra * zinfrastructure + 
        b_pop * zpop + b_lit * zilit + b_us * zius_aid_pc,
    c(a_cy, b_cy)[cy] ~ multi_normal(c(a, b), Rho, sigma_cy),
    a ~ normal(0, 2),
    b ~ normal(0, 0.5),
    c(b_infra, b_lit, b_lit, b_pop, b_us) ~ normal(0, 1),
    c(sigma_cy, sigma) ~ exponential(1),
    Rho ~ lkj_corr(2)
  ),
  data = mdata, chains = 1, cores = 1, iter = 2000
)

m5 = ulam(
  alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[cy] + b_cy[cy] * zigdp_pc + b_infra * zinfrastructure + b_year * zyear +
        b_pop * zpop + b_lit * zilit + b_us * zius_aid_pc,
    c(a_cy, b_cy)[cy] ~ multi_normal(c(a, b), Rho, sigma_cy),
    a ~ normal(0, 2),
    b ~ normal(0, 0.5),
    b_year ~ normal(0, 0.5),
    b_infra ~ normal(0, 0.5),
    b_lit ~ normal(0, 0.5),
    b_pop ~ normal(0, 0.5),
    b_us ~ normal(0, 0.5),
    sigma_cy ~ exponential(1),
    sigma ~ exponential(1),
    Rho ~ lkj_corr(2)
  ),
  data = mdata, chains = 1, cores = 1, iter = 2000
)


cnames = paste0("Model ", 1:length(models))
custom_coeff_map = list(
    Intercept = "Constant", 
    "zigdp_pc" = "GDP per capita", 
    "zinfrastructure" = "Infrastructure", 
    "zilit" = "Literacy", 
    "zius_aid_pc" = "US aid per capita",
    "zpop" = "Population",
    "zyear" = "Year",
    "ctryear::sd(Intercept)" = "Intercept", 
    "ctryear::sd(zigdp_pc)" = "GDP", 
    "ctryear::cor(Intercept,zigdp_pc)" = "cor(Intercept, GDP)"
)

caption = paste0("Models for LE and GPD", 
    nsamples, " replicates")

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

# send message to slack
slackr::slackr_msg(txt = paste0("LE models error finished at: ", Sys.time()))

library(rstan)
example(stan_model,run.dontrun = TRUE)