##########################
# load and impute data
# author: sebastian daza
##########################

# libraries, functions and options
library(haven)
library(data.table)
library(stringr)

library(texreg)
library(mice)
library(miceadds)

library(ggplot2)
library(patchwork)

nsamples = 100
seed = 19380302
set.seed(seed)

source("src/utils.R")

# f (false) or t (true)
select_estimates = "f"

# paths
plots_path = "output/plots/"
tables_path = "output/tables/"
data_path = "output/data/"
manus_plots = "manuscript/plots"
manus_tables  = "manuscript/tables"

# year dataset
dat = data.table(read_stata("data/Ex_LA1850-2020_SES_FULL_Mar2-2021.dta"))
setnames(dat, names(dat), tolower(names(dat)))

ctrylabs = attr(dat$ctry, "labels")
lab_list = as.list(setNames(attr(ctrylabs, "names"), as.numeric(ctrylabs)))
levels = as.numeric(ctrylabs)
labs = attr(ctrylabs, "names")
dat[, ctryf := factor(ctry, labels = labs, level = levels)]

# male
dat = dat[sex == 1 & age == 0 & year >= 1900 & year < 2020]
# remove LE duplicates
dat[, N := 1:.N, .(ctry, year, ex)]
dat = dat[N == 1]
anyDuplicated(dat[, .(ctry, year, ex)])

# country labels
dat[, ctry := as.numeric(ctry)]
dat[, ctryf := droplevels(ctryf)]
table(dat$ctry)
table(dat$ctryf)

labels = names(attr(dat$name, "labels"))
levels = as.numeric((attr(dat$name, "labels")))
dat[, lnames := factor(name, levels = levels, labels = labels)]

dat[, selection := 1]
dat[year < 1950, selection := 0]
dat[year < 1950 & piv == 2, selection := 1]
dat[year < 1950 & name == 1 & (piv == 0 | piv == 1), selection := 1]

countries = unique(dat$ctry)
dat[, fselection := factor(selection,  labels = c("Removed", "Included"))]

savepdf(paste0(plots_path, "selection"))
for (i in countries) {
    print(
        ggplot(dat[ctry == i], aes(year, ex, color = fselection)) +
        geom_jitter(size = 0.5) +
        theme_minimal() +
        labs(title = lab_list[[as.character(i)]], 
            x = "\nYear", y = "Life expectanty at age 0\n") +
        theme(legend.position = "top", legend.title=element_blank())
    )
}
dev.off()
file.copy(paste0(plots_path, "selection.pdf"), manus_plots, 
    recursive = TRUE)

# problematics rows 
# dat[, selection := 1]
# dat[ctry == 2340 & year == 1936, selection := 0]

# selection cases
if (select_estimates == "t") {
    dat = data.table::copy(dat[selection == 1])
} 

dat[, N := .N, .(ctry, year)]
table(dat$N)
anyDuplicated(dat[, .(ctry, year, ex)])

# # testing cases
# test = dat[, .(.N, sd = sd(ex)) ,.(ctry, year)]
# test[year < 1950 & ctry == 2020, .(ctry, year, sd, N)]
# test[year < 1950 & ctry == 2460, .(ctry, year, sd, N)]
# test[year < 1950 & ctry == 2060, .(ctry, year, sd, N)]
# test[year < 1950 & ctry == 2130, .(ctry, year, sd, N)]
# test[year < 1950 & ctry == 2280, .(ctry, year, sd, N)]

# # testing values
# test = dat[ctry == 2020]
# test[year == 1904, .(ctry, year, ex, tseries1)]
# mean(test[year == 1904, .(ctry, year, ex, tseries1)]$ex)

# create variables
# max ex = 78.6
dat[, max_ex := max(ex)]

# dat[, ex_mean := mean(ex), by = .(ctry, year)] # to recover values later
# dat[, wy_mean := mean(wy), by = .(ctry, year)] # to recover values later

# estimate standard deviation based on estimate
sdat = dat[, .(
    max_ex = max(max_ex),
    ex_mean = mean(ex),
    ex_sd = sd(ex),
    N = .N,
    gdp = getMin(gdp_pc),
    urban = getMin(urban),
    pop = getMin(pop)),
    .(ctry, ctryf, year)]

sdat[, wy_mean := transWeibull(ex_mean, max_ex)]

sampex = list()
for (i in 1:nsamples)  {
    sampex[[i]] = dat[, .(ex = getSample(ex)), .(ctry, year, max_ex)][, 
        wy := transWeibull(ex, max_ex)][, .imp := i]
}

# log gdp
sdat[, log_gdp := scale(log(gdp), scale = FALSE)]
sdat[, zpop := scale(pop)]
sdat[, zyear := scale(year)]

# recode year
sdat[year < 1950, gyear := "1950"]
sdat[year >= 1950 & year < 1970, gyear := "1950-1969"]
sdat[year >= 1970 & year < 1990, gyear := "1970-1989"]
sdat[year >= 1990, gyear := "1990"]
sdat[, ctryear := paste0(ctry,".", gyear)]
sdat[, ctryearg := .GRP, ctryear]

# no gdp records before 1950 for country 2170
sdat[ctry == 2170 & ctryearg == 29, ctryearg := 30]
sdat[, year1950 := ifelse(year < 1950, "1950", "1950+")]

# flag missing records
sdat[, gdp_missing := ifelse(is.na(gdp), "missing", "observed")]
sdat[, sd_missing := ifelse(is.na(ex_sd), "missing", "observed")]

sdat[, max_ex := NULL]

# correlations
cor(sdat[, .(ex_mean, ex_sd, log_gdp)])

# missing data```
countmis(sdat)

# impute some missing data

# year data
setorder(sdat, ctry, year)

imp = mice(sdat, maxit = 0)
meth = imp$method
meth[] = ""
pred = imp$pred
pred[,] = 0

imp$loggedEvents
meth["zpop"] = "2l.pmm"
meth["log_gdp"] = "2l.norm"
meth["ex_sd"] = "2l.norm"

pred["zpop", c("zyear", "log_gdp", "ex_mean")] = 1
pred["zpop", c("ctryearg")] = -2

pred["ex_sd", c("zyear", "log_gdp", "ex_mean", "zpop")] = 1
pred["ex_sd", c("zyear")] = 2
pred["ex_sd", c("ctryearg")] = -2

pred["log_gdp", c("zyear", "ex_mean", "ex_sd", "zpop")] = 1
pred["log_gdp", c("zyear")] = 2
pred["log_gdp", c("ctryearg")] = -2

pred["log_gdp", ]

# negative sd values
summary(sdat$ex_sd)
post = imp$post
post["ex_sd"] = "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0.001, 17.0))"

# 100 datasets and 10 iterations
# using average LE values
imps = parmice(sdat,
    m = nsamples,
    method = meth,
    post = post,
    maxit = 10,
    predictorMatrix = pred, 
    n.core = 10, 
    n.imp.core = 10,
    cluster.seed = seed)

# create imputation plots
savepdf(paste0(plots_path, select_estimates, "mice"))
    print(plot(imps))
    print(densityplot(imps, ~ ex_sd))
    print(densityplot(imps, ~ log_gdp))
    print(densityplot(imps, ~ zpop))
dev.off()
file.copy(paste0(plots_path, select_estimates, "mice.pdf"), manus_plots, 
    recursive = TRUE)

# combine samples of ex with imputed values
cimps = data.table(complete(imps, action = "long", include = TRUE))
fimps = merge(cimps, rbindlist(sampex), all.x = TRUE, by = c("ctry", "year", ".imp"))
test_list = list()
for (i in 1:nsamples) {
    test_list[[i]] = fimps[.imp == i][, ctry50 := paste0(ctry, ".", year1950)]
}
timps = test_list

# random imputation
idat = timps[[sample(1:nsamples, 1)]]
names(idat)

# GDP checks
countries = unique(sdat[, ctry])
savepdf(paste0(plots_path, select_estimates, "imputation_check_gdp"))
for (i in countries) {
    print(ggplot(data = idat[ctry == i, .(log_gdp, year, gdp_missing)],
        aes(year, log_gdp, color = gdp_missing)) +
        geom_point() +
        labs(title = lab_list[[as.character(i)]], 
        x = "\nYear", y = "Log GDP\n", color = NULL) +
        theme_minimal() + 
        theme(legend.position = "top", legend.title=element_blank())
    ) 
}
dev.off()
file.copy(paste0(plots_path, select_estimates, "imputation_check_gdp.pdf"), 
    manus_plots,, recursive = TRUE)

# create data list
data_list = list("single-imputation" = idat, "aggregate-data"  = sdat, 
    "raw-data" = dat,  "imputations" = timps, "ctrylabels" = lab_list)

saveRDS(data_list, paste0(data_path, select_estimates, "datalist.rds"))