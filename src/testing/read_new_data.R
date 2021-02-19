
library(haven)
library(data.table)

old = data.table(read_dta("data/Ex-LA1850-2020-estimates.dta"))[age == 0][sex == 1]
dat = data.table(read_dta("data/Ex_LA1840-2020_UncertaintyFile_bydecades.dta"))[
    age == 0][sex == 1]


setnames(old, "myear", "year")

names(old)
test = dat[, .N, .(ctry, year)]
testo = old[, .N, .(ctry, year)]

length(unique(dat$ctry))
length(table(test[N == 1, ctry]))
table(test[N == 1, year])

sort(table(test[N == 1, ctry]))

table(test[ctry == 2370, year])
table(test[ctry == 2130, year])
