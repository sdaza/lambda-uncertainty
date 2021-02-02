
library(haven)
library(data.table)

dat = data.table(read_dta("./data/Ex_LA1850-2020_SES_FULL_Jan25-2021.dta"))[
    age == 0]

dat = dat[sex == 1]
test = dat[, .N, .(ctry, year)]
table(test$N)


length(unique(dat$ctry))
length(table(test[N == 1, ctry]))
table(test[N == 1, year])

sort(table(test[N == 1, ctry]))

table(test[ctry == 2370, year])
table(test[ctry == 2130, year])