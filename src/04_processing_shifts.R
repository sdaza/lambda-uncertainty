ctrys = unique(idat$ctry)
vars = c("ctryear", "ctry50", "ctry", "year", "zyear", "log_gdp", "ex_mean", "wy_mean")


savepdf("output/plots/pred_check_m3")
    print(pp_check(m6))
dev.off()

plots_checks = prediction_checks_pp_ex(ppred, idat, ctrys, vars, y = "ex_mean", x = "year", 
    transform = TRUE)
savepdf("output/plots/fit_no_error_pp")
    print(plots_checks)
dev.off()
file.copy("output/plots/fit_no_error_pp.pdf", "manuscript/plots/", recursive = TRUE)    

# preferred model 
preferred_model = which.max(stacking_wts)
plots_checks = prediction_checks_ex(model_list[[preferred_model]], idat, ctrys, vars, y = "ex_mean", x = "year", transform = TRUE)

savepdf(paste0("output/plots/fit_no_error_m", preferred_model))
    print(plots_checks)
dev.off()
file.copy(paste0("output/plots/fit_no_error_m", preferred_model, ".pdf"), "manuscript/plots/", recursive = TRUE)   

# gpd
savepdf("output/plots/imputation_check_gdp")
for (i in countries) {
    print(ggplot(data = idat[ctry == i, .(log_gdp, year, gdp_missing)],
    aes(year, log_gdp, color = gdp_missing)) +
    geom_point() + labs(title = i, x = "Year", y = "Log GDP", color = NULL))        
}
dev.off()
file.copy("output/plots/imputation_check_gdp.pdf", "manuscript/plots", recursive = TRUE)


savepdf("")
vars = names(shift)
shift[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]
datatest

ctry = unique(datatest$ctry)

savepdf("output/plots/shifts_no_error_pp")
    index = c(1, 2)    
    plots ()
    for (i in seq_along(ctry)) {
            out = data.table(est = shift[[paste0("V", index[1])]]  - shift[[paste0("V", index[2])]])
            print(ggplot(out, aes(x = est)) +
                geom_histogram(color = "black", fill = "white", binwidth = 0.5) + 
                theme_minimal() + labs(x = "Estimated shift", y = "Frequency", title = ctry[i])
            )
            index = index + 2
    }
dev.off()
file.copy("output/plots/shifts_no_error_pp.pdf", "manuscript/plots", recursive = TRUE)
