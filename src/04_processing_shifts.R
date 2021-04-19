# create plots with shift under different scenarios
# author: sebastian daza


library(data.table)
library(ggplot2)
source("src/utils.R")

shifte = readRDS("output/models/shifts_error.rds")[[1]][, error := "Bootstrap"]
shiftne = readRDS("output/models/shifts_no_error.rds")[[1]][, error := "Mean"]
shifts = rbind(shifte, shiftne)
countries = na.omit(as.numeric(names(shifts)))

savepdf("output/plots/shifts_comparison")
    for (i in seq_along(countries)) {
            out = shifts[, c(as.character(countries[i]), "error"), with = FALSE]
            setnames(out, as.character(countries[i]), "est")
            print(ggplot(out, aes(x = est, color= error, fill = error)) +
                geom_density(alpha = 0.2) + 
                theme_minimal() + 
                labs(x = "Estimated shift", y = "Density", 
                    title = countries[i]) + 
                theme(legend.position = "top") + 
                theme(legend.title=element_blank()) + 
                scale_color_manual(values=c("#2b8cbe", "#de2d26")) +
                scale_fill_manual(values=c("#2b8cbe", "#de2d26"))
            )
    }
dev.off()
file.copy("output/plots/shifts_comparison.pdf", "manuscript/plots", recursive = TRUE)


# fit checks plots 
# todo: move them to le files
# plots_checks = prediction_checks_pp_ex(ppred, idat, ctrys, vars, y = "ex_mean", x = "year", 
#     transform = TRUE)
# savepdf("output/plots/fit_no_error_pp")
#     print(plots_checks)
# dev.off()
# file.copy("output/plots/fit_no_error_pp.pdf", "manuscript/plots/", recursive = TRUE)    

# # preferred model 
# preferred_model = which.max(stacking_wts)
# plots_checks = prediction_checks_ex(model_list[[preferred_model]], idat, ctrys, vars, y = "ex_mean", x = "year", transform = TRUE)

# savepdf(paste0("output/plots/fit_no_error_m", preferred_model))
#     print(plots_checks)
# dev.off()
# file.copy(paste0("output/plots/fit_no_error_m", preferred_model, ".pdf"), "manuscript/plots/", recursive = TRUE)   