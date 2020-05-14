##########################
# create bootstrap samples
# author: sebastian daza
##########################


library(data.table)
library(haven)

df = data.table(read_dta('data/probability-assignment-2018-11-19.dta'))
country_labels = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
                   "Costa_Rica", "Cuba", "Dominican_Republic", "Ecuador",
                   "El_Salvador", "Guatemala", "Honduras", "Mexico", "Nicaragua",
                   "Panama", "Paraguay", "Peru", "Uruguay", "Venezuela")
df[, ctry := factor(ctry, labels=country_labels)]
df = melt(df, id_vars = c('ctry', 'year'), measure = patterns("^Ex", "^pr"),
    value.name = c("le", "pr"))
df = df[!is.na(le)]
setorder(df, ctry, year)


df[ctry == "Argentina" & year == 1869][, mean(le)]
df[ctry == "Argentina" & year == 1900]
df[, le := round(le, 2)]


covs = fread('data/featured_LE_data.csv')[, .(ctry, year, Ex)]
setnames(covs, "Ex", "le")
covs[, le := round(le, 2)]
df[, max_pr := max(pr), .(ctry, year)]

test = merge(covs, df, all.x = TRUE, by = c("ctry", "year", "le"))
summary(test)
test[is.na(pr)]
test[pr < max_pr]

test[is.na(pr)]
covs[ctry == "Argentina" & year == 1900]
covs[, c("Ex", "y", "wy", "max_le") := NULL]
df = merge(df, covs, by = c('ctry','year'))

# test distributions
test = df[, .(min_le = min(le),
    max_le = max(le), min_pr = min(pr),
    max_pr = max(pr)),
    .(ctry, year)]

test[, diff_years := max_le - min_le]
test[, diff_prob := max_pr - min_pr]

# sampling
sample_size = 1000
set.seed(seed)

df = df[year >= 1900]
# uniform probability
df[, upr := 1 / .N, .(ctry, year)]
# expert probability```
df[, spr := pr / sum(pr), by=.(ctry, year)]
df[, N := .N, .(ctry, year)]

uniform_samples = df[,.SD[ sample(.N, sample_size, replace = TRUE, prob=upr)], by = .(ctry, year)]
uniform_samples[, sample_index := 1:.N, by=.(ctry, year)]
dim(uniform_samples)

samples = df[,.SD[ sample(.N, sample_size, replace = TRUE, prob = spr)], .(ctry,year)]
samples[, sample_index := 1:.N, by=.(ctry, year)]

# save files
fwrite(samples, '../data/bs_samples_no_uniform.csv', row.names=FALSE)
fwrite(uniform_samples, '../data/bs_samples_uniform.csv', row.names=FALSE)