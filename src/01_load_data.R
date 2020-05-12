##########################
# load and impute data
# author: sebastian daza
##########################


# load libraries
library(haven)
library(data.table)
library(imputeTS)
library(stringr)
library(texreg)

seed = 850091718
set.seed(seed)
source("src/utils.R")

# laod data
df = data.table(read_stata('data/le-1850-2013-abbr-2018-03-03.dta'))

# add country labels
country_labels = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
    "Costa_Rica", "Cuba", "Dominican_Republic", "Ecuador",
    "El_Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua",
    "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela")
df[, ctry := factor(ctry, labels = country_labels)]

# covariates
covariates = c("gdp_pc", "urban", "elec", "lit" ,"water" ,"sewage" ,"gini" ,
    "tfr" ,"bf" ,"extFund" ,"healthGdp" ,"lsi" ,"polio" ,"bcg" ,"dpt1" ,
    "dpt3" ,"mcv1", "us_aid" ,"lunion")

c = df[tseries2 == 1 & age == 0 & year >= 1900, lapply(.SD, getMax),
    .SDcols = covariates, .(ctry, year)]

# average men and women
le = df[tseries2==1 & age == 0 & year>=1900,
    .(Ex = mean(Ex, na.rm = TRUE)), .(ctry, year)]

df = merge(c, le, by = c('ctry', 'year'))

if (nrow(le) != nrow(df)) {
    stop("Number of rows is not the same after merging datasets!")
}

# select columns
df = df[year>=1900, .(ctry, year, gdp_pc, urban, lit, Ex, water, sewage, elec, us_aid, tfr)]

# create year groups
df[year < 1950, gyear := '1950']
df[year >= 1950 & year < 1970, gyear := '1950-1969']
df[year >= 1970 & year < 1990, gyear := '1970-1989']
df[year >= 1990, gyear := '1990']

# df[, gyear := factor(gyear, levels=1:4, labels=c('1950', '1950-1969', '1970-1989', '1990'))]

# transform variable: weibull
# adjustment is by country!
df[, y := Ex / max(Ex + 1.05), ctry]
df[, wy := log(-log(1 - y))]
# to recover values later
df[, max_le := max(Ex + 1.05), by = ctry]
df[, ctry_year := paste0(ctry,'.', gyear)]

# interpolation
setorder(df, year)
interpolation_vars = c("gdp_pc", "urban", "lit", "tfr", "water", "sewage",
    "elec")
df[, (paste0("i", interpolation_vars)) := lapply(.SD, na_interpolation,
    option = "stine"), ctry, .SDcols = interpolation_vars]

# create plots interpolation
savepdf("output/plots/interpolation_by_country")
for (i in country_labels) {
    for (ii in interpolation_vars) {
        variable = ii
        imputed_variable = paste0("i", variable)
        plotNA.imputations(df[ctry == i][[variable]],df[ctry == i][[imputed_variable]],
            legend = TRUE, main = paste0("Interpolation: ", i, " ", variable))
    }
}
dev.off()

# create log version of variables
vars = c("gdp_pc", "urban", "lit", "water", "sewage", "elec")
ivars = paste0("i", vars)
nvars = c("igdp_log", "iurban_log", "ilit_log", "iwater_log", "isewage_log", "ielec_log")
df[, (nvars) := lapply(.SD, function(x) scale(log(x), scale = FALSE)),
    .SDcols = ivars]

# z-score
df[, zyear := scale(year, center = TRUE, scale = TRUE)]
fwrite(df, 'data/featured_LE_data.csv', row.names = FALSE)