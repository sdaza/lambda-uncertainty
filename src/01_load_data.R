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

# year dataset
dat = data.table(read_stata("data/Ex_LA1850-2020_SES_FULL_Mar2-2021.dta"))
setnames(dat, names(dat), tolower(names(dat)))

ctrylabs = attr(dat$ctry, "labels")
levels = as.numeric(ctrylabs)
labs = attr(ctrylabs, "names")
dat[, ctryf := factor(ctry, labels = labs, level = levels)]

# male
dat = dat[sex == 1 & age == 0 & year >= 1900 & year < 2020]
dat[, N := 1:.N, .(ctry, year, ex)]
dat = dat[N == 1]
anyDuplicated(dat[, .(ctry, year, ex)])
dat[, ctry := as.numeric(ctry)]
dat[, ctryf := droplevels(ctryf)]
table(dat$ctry)
table(dat$ctryf)

labels = names(attr(dat$name, "labels"))
levels = as.numeric((attr(dat$name, "labels")))
dat[, lnames := factor(name, levels = levels, labels = labels)]
# str(dat)

dat[, selection := 1]
dat[year < 1950, selection := 0]
dat[year < 1950 & piv == 2, selection := 1]
dat[year < 1950 & name == 1 & (piv == 0 | piv == 1), selection := 1]

# problematics rows 
# dat[, selection := 1]
# dat[ctry == 2340 & year == 1936, selection := 0]

# selection cases
dat = data.table::copy(dat[selection == 1])
dat[, N := .N, .(ctry, year)]
table(dat$N)
anyDuplicated(dat[, .(ctry, year, ex)])

# testing values
test = dat[ctry == 2020]
test[year == 1904, .(ctry, year, ex, tseries1)]
mean(test[year == 1904, .(ctry, year, ex, tseries1)]$ex)

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

# five datasets and 10 iterations
imps = mice(sdat,
    m = nsamples,
    method = meth,
    post = post,
    maxit = 10,
    predictorMatrix = pred
)

# create imputation plots
savepdf("output/plots/mice")
    plot(imps)
    densityplot(imps, ~ ex_sd)
    densityplot(imps, ~ log_gdp)
    densityplot(imps, ~ zpop)
dev.off()
file.copy("output/plots/mice.pdf", "manuscript/plots/", recursive = TRUE)

# combine samples of ex with imputed values
cimps = data.table(complete(imps, action = "long", include = TRUE))
fimps = merge(cimps, rbindlist(sampex), all.x = TRUE, by = c("ctry", "year", ".imp"))
test_list = list()
for (i in 1:nsamples) {
    test_list[[i]] = fimps[.imp == i][, ctry50 := paste0(ctry, ".", year1950)]
}
timps = test_list

idat = timps[[sample(1:nsamples, 1)]]
names(idat)

countries = unique(sdat[, ctry])
savepdf("output/plots/imputation_check_gdp")
for (i in countries) {
    print(ggplot(data = idat[ctry == i, .(log_gdp, year, gdp_missing)],
        aes(year, log_gdp, color = gdp_missing)) +
        geom_point() + labs(title = i, x = "Year", y = "Log GDP", color = NULL)) 

}
dev.off()
file.copy("output/plots/imputation_check_gdp.pdf", "manuscript/plots", recursive = TRUE)

# save files
saveRDS(idat, "output/data/single-imputation.rds")
saveRDS(sdat, "output/data/aggregate-data.rds")
saveRDS(dat, "output/data/raw-data.rds")
saveRDS(timps, "output/data/imputations.rds")