################################################
# create plots with shift under different scenarios
# author: sebastian daza
################################################


library(data.table)
library(ggplot2)
library(xtable)
library(stringr)
source("src/utils.R")

# paths
plots_path = "output/plots/"
tables_path = "output/tables/"
data_path = "output/data/"
manus_plots = "manuscript/plots"
manus_tables  = "manuscript/tables"

# read data
select_estimates = "t"
data_list = readRDS(paste0(data_path, select_estimates, "datalist.rds"))
idat = data_list[["single-imputation"]]
country_labs = data_list[["ctrylabels"]]

shifte = readRDS(paste0(data_path, select_estimates, "shifts_error.rds"))[[1]][,
    error := "Bootstrap"]
shiftne = readRDS(paste0(data_path, select_estimates, "shifts_no_error.rds"))[[1]][,
    error := "Mean"]

# create table
shifte = shifte[sample(.N, nrow(shiftne))]

b = data.table::copy(shifte)[, c("replicate", "error") := NULL]
a = data.table::copy(shiftne)[, c("replicate", "error") := NULL]
d = a - b 

summaryFun = function(x) {
    list("Mean" = mean(x), 
    "Q2.5" = quantile(x, 0.025), 
    "Q97.5" = quantile(x, 0.975))
}
sa = t(sapply(a, summaryFun))
sb = t(sapply(b, summaryFun))
sd = t(sapply(d, summaryFun))
tab = cbind(sa, sb, sd)
rownames(tab) = as.vector(unlist(country_labs[rownames(tab)]))
tab = data.table(tab, keep.rownames = TRUE)
ctab = print(xtable(tab), include.rownames = FALSE)

tab_list = list()
tab_list[[1]] = "\\renewcommand{\\arraystretch}{1.2}
\\setlength{\\tabcolsep}{10pt}
\\begin{table}[htp]
\\centering
\\caption{1950's shifts by country (selected estimates)} 
\\label{tab:tshift_1950}
\\scriptsize
\\begin{tabular}{lrrrrrrrrr}
  \\hline
  \\addlinespace
    & \\multicolumn{3}{c}{Average (A)}  & \\multicolumn{3}{c}{Error (B)} & \\multicolumn{3}{c}{Difference (A-B)}  \\\\
    \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}
 & Mean & Q2.5 & Q97.5 & Mean & Q2.5 & Q97.5 & Mean & Q2.5 & Q97.5 \\\\ 
 \\addlinespace
  \\hline
\\addlinespace"
tab_list[[2]] = str_split(ctab, "\\\\hline")[[1]][[3]]
tab_list[[3]] = "\\addlinespace
\\hline
\\end{tabular}
\\end{table}
"

cat(paste0(unlist(tab_list)), file = paste0(tables_path, select_estimates, "shift_1950.tex"))
file.copy(paste0(tables_path, select_estimates, "shift_1950.tex"), manus_tables, recursive = TRUE)

# plot with differences
shifts = rbind(shifte, shiftne)
countries = na.omit(as.numeric(names(shifts)))

savepdf(paste0(plots_path, select_estimates, "shifts_1950_comparison"))
    for (i in seq_along(countries)) {
        out = shifts[, c(as.character(countries[i]), "error"), with = FALSE]
        setnames(out, as.character(countries[i]), "est")
        print(ggplot(out, aes(x = est, color= error, fill = error)) +
            geom_density(alpha = 0.2) + 
            theme_minimal() + 
            labs(x = "Estimated shift", y = "Density", 
                title = country_labs[[as.character(countries[i])]]) + 
            theme(legend.position = "top") + 
            theme(legend.title=element_blank()) + 
            scale_color_manual(values=c("#2b8cbe", "#de2d26")) +
            scale_fill_manual(values=c("#2b8cbe", "#de2d26"))
        )
    }
dev.off()
file.copy(paste0(plots_path, select_estimates, "shifts_1950_comparison.pdf"), 
    manus_plots, recursive = TRUE)


# read data all estimates
select_estimates = "f"
data_list = readRDS(paste0(data_path, select_estimates, "datalist.rds"))
idat = data_list[["single-imputation"]]
country_labs = data_list[["ctrylabels"]]

shifte = readRDS(paste0(data_path, select_estimates, "shifts_error.rds"))[[1]][,
    error := "Bootstrap"]
shiftne = readRDS(paste0(data_path, select_estimates, "shifts_no_error.rds"))[[1]][,
    error := "Mean"]

# create table
shifte = shifte[sample(.N, nrow(shiftne))]

b = data.table::copy(shifte)[, c("replicate", "error") := NULL]
a = data.table::copy(shiftne)[, c("replicate", "error") := NULL]
d = a - b 

summaryFun = function(x) {
    list("Mean" = mean(x), 
    "Q2.5" = quantile(x, 0.025), 
    "Q97.5" = quantile(x, 0.975))
}
sa = t(sapply(a, summaryFun))
sb = t(sapply(b, summaryFun))
sd = t(sapply(d, summaryFun))
tab = cbind(sa, sb, sd)
rownames(tab) = as.vector(unlist(country_labs[rownames(tab)]))
tab = data.table(tab, keep.rownames = TRUE)
ctab = print(xtable(tab), include.rownames = FALSE)

tab_list = list()
tab_list[[1]] = "\\renewcommand{\\arraystretch}{1.2}
\\setlength{\\tabcolsep}{10pt}
\\begin{table}[htp]
\\centering
\\caption{1950's shifts by country (all estimates)} 
\\label{tab:fshift_1950}
\\scriptsize
\\begin{tabular}{lrrrrrrrrr}
  \\hline
  \\addlinespace
    & \\multicolumn{3}{c}{Average (A)}  & \\multicolumn{3}{c}{Error (B)} & \\multicolumn{3}{c}{Difference (A-B)}  \\\\
    \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}
 & Mean & Q2.5 & Q97.5 & Mean & Q2.5 & Q97.5 & Mean & Q2.5 & Q97.5 \\\\ 
 \\addlinespace
  \\hline
\\addlinespace"
tab_list[[2]] = str_split(ctab, "\\\\hline")[[1]][[3]]
tab_list[[3]] = "\\addlinespace
\\hline
\\end{tabular}
\\end{table}
"
cat(paste0(unlist(tab_list)), file = paste0(tables_path, select_estimates, "shift_1950.tex"))
file.copy(paste0(tables_path, select_estimates, "shift_1950.tex"), manus_tables, recursive = TRUE)

# plot with differences
shifts = rbind(shifte, shiftne)
countries = na.omit(as.numeric(names(shifts)))

savepdf(paste0(plots_path, select_estimates, "shifts_1950_comparison"))
    for (i in seq_along(countries)) {
        out = shifts[, c(as.character(countries[i]), "error"), with = FALSE]
        setnames(out, as.character(countries[i]), "est")
        print(ggplot(out, aes(x = est, color= error, fill = error)) +
            geom_density(alpha = 0.2) + 
            theme_minimal() + 
            labs(x = "Estimated shift", y = "Density", 
                title = country_labs[[as.character(countries[i])]]) + 
            theme(legend.position = "top") + 
            theme(legend.title=element_blank()) + 
            scale_color_manual(values=c("#2b8cbe", "#de2d26")) +
            scale_fill_manual(values=c("#2b8cbe", "#de2d26"))
        )
    }
dev.off()
file.copy(paste0(plots_path, select_estimates, "shifts_1950_comparison.pdf"), 
    manus_plots, recursive = TRUE)

