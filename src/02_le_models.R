#####################################
# life expectancy models with error
# author: sebastian daza
#####################################


# R < src/03_le_models.R > output/log/03_le_models.log  --no-save  &

# libraries, functions and options
library(data.table)
library(doParallel)
library(foreach)
library(texreg)
library(rethinking)
library(ggplot2)
library(factoextra)
library(RColorBrewer)
library(rlist)
rstan_options(auto_write = FALSE)

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
#   data = mdata, chains = 4, cores = 4, iter = 4000, cmdstan = TRUE
# )

# rethinking::precis(model, depth = 1)

# loop by model
iterations = c(6000, 4000, 4000, 4000, 4000)

for (m in seq_along(iterations)) {

    # get model output
    model_number = m
    output = runModel(formulas[[model_number]], samples, newdata = newdata, ex_max = ex_max, 
        iterations = iterations[model_number], clusters = 10)

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

    # save shifts
    # shifts = list()
    shifts = readRDS("output/data/shifts.rds")
    shifts[[model_number]] = output$shifts
    saveRDS(shifts, "output/data/shifts.rds")
    rm(shifts)
    
    # shift plot
    shift_plot = plotShifts(output$shifts, country_labs)
    savepdf(paste0(plots_path, select_estimates, 
        paste0("shifts_by_period_m", model_number)), height = 20)
        print(shift_plot)
    dev.off()
    file.copy(paste0(plots_path, select_estimates, 
        paste0("shifts_by_period_m", model_number, ".pdf")), manus_plots, recursive = TRUE)

    # shifts
    v = unlist(country_labs)
    tabshift = output$shifts
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
    slackr::slackr_msg(txt = paste0("Stan loop ", m, " finished at: ", Sys.time()))

}

# create cross model shift plots
shifts = readRDS("output/data/shifts.rds")
shifts = rbindlist(shifts, idcol = "model")

for (i in c(1950, 1970, 1990)) {
    print(paste0(":::::: plotting year ", i))
    name = paste0(plots_path, select_estimates, "shifts_by_year_", i)
    savepdf(name, height = 20)
        print(plotShiftModels(shifts[year == i], country_labs))
    dev.off()
    file.copy(paste0(name, ".pdf"), manus_plots, 
        recursive = TRUE)
}

# cluster analysis
shifts = readRDS("output/data/shifts.rds")
shifts = rbindlist(shifts, idcol = "model")
shifts = shifts[, .(shift = median(shift)), .(model, year, ctry)]
shifts = shifts[model %in% c(3, 5)]
shifts = shifts[year %in% c(1950, 1970, 1990)]
shifts = dcast(shifts, ctry ~ year + model, value.var = "shift")
v = unlist(country_labs)
shifts[, lctry := as.factor(v[as.character(ctry)])]

dt = data.frame(scale(shifts[,.SD, .SDcols = !c('ctry', 'lctry')]))
rownames(dt) = shifts$lctry
dt = dt[complete.cases(dt),]

savepdf("output/plots/nclusters")
fviz_nbclust(x = dt, FUNcluster = kmeans, method = "silhouette", k.max = 15) +
  labs(title = "Number of clusters")
dev.off()

hc = hcut(dt, 6)
savepdf("output/plots/dendogram", width = 18, height = 18)
fviz_dend(x = hc, cex = 0.5, k = 6)
dev.off()

savepdf("output/plots/silhouette")
fviz_silhouette(hc)
dev.off()


# hk = hkmeans(x = dt, hc.metric = "euclidean", hc.method = "complete", k = 7)
hk = kmeans(x = dt, 7)
savepdf("output/plots/clusters")
fviz_cluster(object = hc, data = dt, palette = "Set2", repel = TRUE) +
  theme_bw() + labs(title = "Clustering")
dev.off()

# create plot grouping by cluster 
colors = brewer.pal(6, "Set2")
shifts = readRDS("output/data/shifts.rds")
cluster = data.table(lctry = names(hc$cluster), cluster = hc$cluster)
setorder(cluster, cluster)
lab_order = as.character(cluster$lctry)
cluster[, color :=  colors[cluster]]

vlist = list()
for (i in seq_along(lab_order)) {
    vlist = append(vlist, list.search(country_labs, . == lab_order[i])) 
}

shift_plot = plotShifts(shifts[[5]], vlist, 
    color = cluster$color)
savepdf(paste0(plots_path, select_estimates, 
    paste0("shifts_by_period_ordered_m", 5)), height = 20)
    print(shift_plot)
dev.off()
file.copy(paste0(plots_path, select_estimates, 
    paste0("shifts_by_period_ordered_m", 5, ".pdf")), manus_plots, recursive = TRUE)


# create regression table
select_estimates = ""
tabs = readRDS(paste0(tables_path, "tabs.rds"))

# testing 
# screenreg(tabs, 
#   omit.coef = "^a\\_cy\\[.+|^b_gdp\\_cy\\[.+|^Rho\\[1,1\\]|^Rho\\[2,2\\]" )

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
        \\item All covariates are standardized. Life expectancy estimates were transformed using $ln\\left(-ln( 1-\\frac{e_0}{ (79.45 + 1.05)}\\right)$.",
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

# final message
slackr::slackr_msg(txt = paste0("Overall code finished at: ", Sys.time()))