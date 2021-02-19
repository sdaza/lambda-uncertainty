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
df = data.table(read_stata("data/Ex_LA1840-2020_UncertaintyFile_bydecades.dta"))

# add country labels
# country_labels = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
#     "Costa_Rica", "Cuba", "Dominican_Republic", "Ecuador",
#     "El_Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua",
#     "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela")
# df[, ctry := factor(ctry, labels = country_labels)]
table(df$ctry)
setnames(df, "myear", "year")

# covariates
covs = df[age == 0 & sex == 1][, .(gdp = getMax(gdp_pc)), .(ctry, year)]
nrow(covs)
summary(covs$year)

# life expectancy
le = df[age == 0 & sex == 1]
le[, le_std := scale(Ex)]
le[, le_sd := sd(Ex), .(ctry, year)]
le[, le_sd := le_sd / sd(Ex)]

hist(le$le_sd)
hist(le$le_std)

# transform life expectancy variable: weibull
# adjustment is by country!
# le[, y := Ex / max(Ex + 1.05), ctry]
# le[, wy := log(-log(1 - y))]

# to recover values later
# le[, max_le := max(Ex + 1.05), by = ctry]

# hist(le$wy)

# interpolation covariates
setorder(covs, year)
interpolation_vars = c("gdp")
covs[, (paste0("i", interpolation_vars)) := lapply(.SD, na_interpolation,
    option = "stine"), ctry, .SDcols = interpolation_vars]
summary(covs)

# create plots interpolation
countries_with_missing = unique(covs[is.na(gdp), ctry])
savepdf("output/plots/interpolation_by_country")
for (i in countries_with_missing) {
    for (ii in interpolation_vars) {
        variable = ii
        imputed_variable = paste0("i", variable)
        print(
            ggplot_na_imputations(covs[ctry == i][[variable]], covs[ctry == i][[imputed_variable]],
            legend = TRUE, title = paste0("Interpolation: ", i, " ", variable))
        )
    }
}
dev.off()

# create log version of variables
vars = c("igdp")
covs[, (paste0(vars, "_log")) := lapply(.SD, function(x) log(x)),
    .SDcols = vars]
summary(covs)
hist(covs$igdp_log)
sd_igdp_log = sd(covs$igdp_log)
covs[, igdp_std := scale(igdp_log)]
hist(covs$igdp_std)

# merge datasets
df = merge(le, covs[, .(ctry, year, igdp, igdp_log, igdp_std)], by = c("ctry", "year"))

if (nrow(le) != nrow(df)) {
    stop("Number of rows is not the same after merging datasets!")
}

# compute standard deviation per estimates
adf = df[, .(le_std = mean(le_std),
        le_sd = mean(le_sd),
        le_avg = mean(Ex),
        igdp_log = mean(igdp_log),
        igdp_std = mean(igdp_std),
        .N
    ),
    .(ctry, year)]

# # impute standard deviation for only one record per year
# setorder(adf, ctry, year)
# interpolation_vars = c("sd_wy", "sd_ex")
# adf[, (paste0("i", interpolation_vars)) := lapply(.SD, na_interpolation,
#     option = "stine"), ctry, .SDcols = interpolation_vars]

# country_only_one_record = unique(adf[is.na(sd_wy), ctry])
# savepdf("output/plots/interpolation_standard_deviation")
# for (i in country_only_one_record ) {
#     for (ii in interpolation_vars) {
#         variable = ii
#         imputed_variable = paste0("i", variable)
#         print(
#             ggplot_na_imputations(adf[ctry == i][[variable]], adf[ctry == i][[imputed_variable]],
#             legend = TRUE, title = paste0("Interpolation: ", i, " ", variable))
#         )
#     }
# }
# dev.off()

# create group of years
adf[year < 1950, gyear := "1950"]
adf[year >= 1950 & year < 1970, gyear := "1950-1969"]
adf[year >= 1970 & year < 1990, gyear := "1970-1989"]
adf[year >= 1990, gyear := "1990"]
adf[, ctry_year := paste0(ctry,".", gyear)]
adf[, year1950 := ifelse(year < 1950, "1950", "1950+")]

fwrite(adf, "data/featured_LE_data.csv", row.names = FALSE)
