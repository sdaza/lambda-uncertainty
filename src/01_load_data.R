##########################
# load and impute data
# author: sebastian daza
##########################


# libraries, functions and options
library(haven)
library(data.table)
library(stringr)

library(texreg)
library(ggplot2)
# library(naniar)
library(imputeTS)

nsamples = 100
seed = 19380302
set.seed(seed)

source("src/utils.R")
le_estimate_limit = 10/2

# f (false) or t (true)
select_estimates = "t"

# paths
plots_path = "output/plots/"
tables_path = "output/tables/"
data_path = "output/data/"
manus_plots = "manuscript/plots"
manus_tables  = "manuscript/tables"

# year dataset
dat = data.table(read_stata("data/Ex_LA1850-2020_SES_FULL_Jan25-2021.dta"))
setnames(dat, names(dat), tolower(names(dat)))

ctrylabs = attr(dat$ctry, "labels")
lab_list = as.list(setNames(attr(ctrylabs, "names"), as.numeric(ctrylabs)))
levels = as.numeric(ctrylabs)
labs = attr(ctrylabs, "names")
dat[, ctryf := factor(ctry, labels = labs, level = levels)]

dat[, qyear := ifelse(as.numeric(substr(as.character(year), 4, 4)) < 5, 
    year - as.numeric(substr(as.character(year), 4, 4)),
    year - as.numeric(substr(as.character(year), 4, 4)) + 5)]

# average males and females
dat = dat[age == 0 & year >= 1900 & year < 2010]
ex = dat[, .(ex = mean(ex)), .(ctry, ctryf, tseries2, year, name)]
ex = ex[, N := .N, .(ctry, year, ex)][N == 1]

# recode year
covs  = copy(dat[, N := 1:.N, .(ctry, year)][N == 1])
covs[year < 1950, gyear := "1950"]
covs[year >= 1950 & year < 1970, gyear := "1950-1969"]
covs[year >= 1970 & year < 1990, gyear := "1970-1989"]
covs[year >= 1990, gyear := "1990"]
covs[, ctryear := factor(paste0(ctry,".", gyear))]
covs[, ctryearg := .GRP, ctryear]
covs[, year1950 := factor(ifelse(year < 1950, "1950", "1950+"))]
covs[, ctry50 := paste0(ctry, ".", year1950)]

# impute variables (interpolation)
setorder(covs, ctry, year)
ovars = c("gdp_pc", "urban", "elec", "lit", "water", "sewage", "us_aid_pc")
covs[, paste0("i", ovars) := lapply(.SD, na_interpolation, option = "stine"), 
    by = ctry, .SDcols = ovars]
covs[, log_gdp_pc := log(igdp_pc)]

covs[, paste0("zi", ovars) := lapply(.SD, scale), 
    .SDcols = paste0("i", ovars)]
covs[, zpop := scale(pop)]
covs[, zyear := scale(year)]
covs[, zinfrastructure := apply(.SD, 1, mean), .SDcols = c("ziwater", "zisewage", 
    "zielec")]

vars = c("ctry", "year", "qyear", "gyear", "ctryear", "ctry50", "ctryearg", 
    "year1950", "zpop", "zyear", paste0("zi", ovars), "zinfrastructure", 
    "log_gdp_pc", "igdp_pc")

covs = covs[, ..vars]
covs[, ctry := as.numeric(ctry)]

# savepdf(paste0(plots_path, "missing_pattern"))
# gg_miss_fct(covs[, c("year", ovars), with = FALSE], year)
# dev.off()
# file.copy(paste0(plots_path, "missing_pattern.pdf"), manus_plots, 
#     recursive = TRUE)
countmis(covs)

# country labels
ex[, ctry := as.numeric(ctry)]
ex[, ctryf := droplevels(ctryf)]


labels = names(attr(ex$name, "labels"))
levels = as.numeric((attr(ex$name, "labels")))
ex[, lnames := factor(name, levels = levels, labels = labels)]
ex[, lambda := 0][tseries2 == 1, lambda := 1]
ex[lambda == 1, lambda_ex := ex][, lambda_ex := mean(lambda_ex, na.rm = TRUE), .(ctry, year)]
ex[, nlambda := 1]
ex[, nlambda := ifelse(ex >= (lambda_ex + le_estimate_limit) | ex <= (lambda_ex - le_estimate_limit), 0, nlambda)]
ex[, selection := 0]
ex[lambda == 1 | nlambda == 1, selection := 1]
ex[, N := .N, .(ctry, year)]
ex[N == 1 & lambda == 1, selection := 0]

# plot with included values
countries = unique(dat$ctry)
ex[, fselection := factor(selection,  labels = c("Removed", "Included"))]
savepdf(paste0(plots_path, "le_estimate_selection"))
for (i in countries) {
    print(
        ggplot(ex[ctry == i], aes(year, ex, color = fselection)) +
        geom_jitter(size = 0.5) +
        theme_minimal() +
        labs(title = lab_list[[as.character(i)]], 
            x = "\nYear", y = "Life expectancy at age 0\n") +
        theme(legend.position = "top", legend.title=element_blank())
    )
}
dev.off()
file.copy(paste0(plots_path, "le_estimate_selection.pdf"), manus_plots, 
    recursive = TRUE)

# selection of estimates?
if (select_estimates == "t") {
    ex = data.table::copy(ex[selection == 1])
} 

ex[, N := .N, .(ctry, year)]
table(ex$N)
anyDuplicated(ex[, .(ctry, year, ex)])

# max ex = 78.6
ex_max = max(ex$ex)
ex[, max_ex := max(ex)]
ex[, wy := transWeibull(ex, max_ex)]
ex[, ex_mean := mean(ex), by = .(ctry, year)]
ex[, wy_mean := transWeibull(ex_mean, max_ex)]

ex = merge(ex, covs, by = c("ctry", "year"), x.all = TRUE)

samples = list()
for (i in 1:nsamples)  {
    samples[[i]] = ex[, .SD[sample(.N, min(1, .N))], .(ctry, year)]
}

# missing data
countmis(samples[[2]])

# random imputation
idat = samples[[sample(1:nsamples, 1)]]
names(idat)

# data for shifts
lctryear = levels(covs$ctryear)
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

# define models (formulas, rethinking)
formulas = list()

formulas[[1]] = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp * zigdp_pc,
    a_cy[ctryearg] ~ normal(a, sigma_cy),
    a ~ normal(0, 1),
    # c(sigma_cy, sigma) ~ half_normal(0,1),
    c(sigma_cy, sigma) ~ exponential(1),
    b_gdp ~ normal(0, 0.5), 
    save> pred <- a_cy[ctryearg] + b_gdp * zigdp_pc
)

formulas[[2]] = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 2),
    b_gdp ~ normal(0, 1),
    c(sigma, sigma_cy) ~ exponential(1),
    Rho ~ lkj_corr(2),
    save> pred <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc
)

formulas[[3]] = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_year * zyear,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 1),
    c(b_gdp, b_year) ~ normal(0, 0.5),
    c(sigma, sigma_cy) ~ exponential(1),
    Rho ~ lkj_corr(2),
    save> pred <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_year * zyear
)


formulas[[4]] = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_infra * zinfrastructure + 
        b_pop * zpop + b_lit * zilit + b_us * zius_aid_pc,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 1),
    c(b_gdp, b_infra, b_lit, b_pop, b_us) ~ normal(0, 0.5),
    c(sigma_cy, sigma) ~ exponential(1),
    Rho ~ lkj_corr(2),
    save> pred <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_infra * zinfrastructure + 
        b_pop * zpop + b_lit * zilit + b_us * zius_aid_pc
)

formulas[[5]] = alist(
    wy ~ normal(mu, sigma),
    mu <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_infra * zinfrastructure + 
        b_pop * zpop + b_lit * zilit + b_us * zius_aid_pc + b_year * zyear,
    c(a_cy, b_gdp_cy)[ctryearg] ~ multi_normal(c(a, b_gdp), Rho, sigma_cy),
    a ~ normal(0, 1),
    c(b_gdp, b_infra, b_lit, b_pop, b_us, b_year) ~ normal(0, 0.5),
    c(sigma_cy, sigma) ~ exponential(1),
    Rho ~ lkj_corr(2),
    save> pred <- a_cy[ctryearg] + b_gdp_cy[ctryearg] * zigdp_pc + b_infra * zinfrastructure + 
        b_pop * zpop + b_lit * zilit + b_us * zius_aid_pc + b_year * zyear
)

# correlation between indicators
vars = c("zpop", "zyear", "zigdp_pc", "ziurban", "zielec", "zilit", "ziwater",
     "zisewage", "zius_aid_pc", "zinfrastructure")
r = cor(idat[, ..vars])
savepdf(paste0(plots_path, select_estimates, "correlation_imputed_covs"))
ggcorrplot::ggcorrplot(r, 
    hc.order = TRUE, 
    type = "lower",
    lab = TRUE)
dev.off()
file.copy(paste0(plots_path, select_estimates, "correlation_imputed_covs.pdf"), 
    manus_plots, recursive = TRUE)    

# all data points
savepdf(paste0(plots_path, select_estimates, "le_gdp_data"), 24, 16)
ggplot(idat, aes(zigdp_pc, wy_mean, color= ctryf, shape = gyear))  +
 geom_point(size = 1)  +  
 theme_minimal() +
 theme(legend.position="right", 
        legend.title = element_blank()) + 
labs(x = "\nGDP (z-score)", y = "Transformed e0\n")
dev.off()
file.copy(paste0(plots_path, select_estimates, "le_gdp_data.pdf"), 
    manus_plots, recursive = TRUE)

# create data list
data_list = list("single-imputation" = idat, 
    "raw-data" = dat,  "imputations" = samples, "ctrylabels" = lab_list,
    "covs" = covs, ex_max = ex_max, newdata = newdata, formulas = formulas)
    
saveRDS(data_list, paste0(data_path, select_estimates, "datalist.rds"))
