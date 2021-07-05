runModel = function(flist, samples, clusters = 2, iterations = 2000, 
    model_number = 1, ex_max, country_labs, newdata, plots_path, 
    selected_estimates = "") {

    tabs = list()
    idat = samples[[1]]
    cl = makeCluster(clusters, outfile="")
    nsamples = length(samples)
    registerDoParallel(cl)
    output = foreach(i = 1:nsamples) %dopar% {
        library(data.table)
        library(rethinking)
        source("src/utils.R")
        results = multiResultClass()

        print(paste0(":::::::: Running iteration ", i, " ::::::::"))
        dat = samples[[i]]
        dat[, cy := as.numeric(ctryear)]
        mdata = list(
            wy = dat$wy, 
            zyear = dat$zyear,
            zinfrastructure = dat$zinfrastructure, 
            zpop= dat$zpop, 
            zilit = dat$zilit, 
            zius_aid_pc  = dat$zius_aid_pc, 
            zigdp_pc = dat$zigdp_pc, 
            ctryearg = dat$ctryearg
        )
        model = ulam(
            flist, 
            data = mdata, chains = 1, cores = 1, 
            iter = iterations, chain_id = i
        )
        
        results$models = model@stanfit
        results$shifts = computeShift(model, newdata, ex_max)
        return(results)
    }
    stopCluster(cl)

    # extract restuls
    shifts = list()
    models = list()
    for (i in seq_along(output)) {
        shifts[[i]] = output[[i]][["shifts"]]
        models[[i]] = output[[i]][["models"]]
    }

    shifts = rbindlist(shifts)
    fit = sflist2stanfit(models)
    tab = extractStan(fit, n = nrow(idat))

    # plots
    pred = predictionData(fit, idat, ex_max = ex_max)
    plots = predictionPlots(pred, country_labels = country_labs)
    plots = wrap_plots(plots, ncol = 3)
    savepdf(paste0(plots_path, select_estimates, 
        paste0("fit_check_m", model_number)), width = 30, height = 35)
        print(plots)
    dev.off()

    xmin = min(shifts$shift, na.rm = TRUE) - 1.5
    xmax = max(shifts$shift, na.rm = TRUE) + 1.5
    v = unlist(country_labs)
    shifts[, lctry := as.factor(v[as.character(ctry)])]
    shifts[, year := factor(year)]
    shifts[, ctry := factor(ctry)]
    savepdf(paste0(plots_path, select_estimates, 
        paste0("shifts_by_period_m", model_number)), height = 20)
        print(ggplot(shifts, aes(y = lctry)) +
        geom_density_ridges(aes(x = shift, fill = year), 
            alpha = .45, color = "white", from = xmin, to = xmax, scale = 1) +
        labs(x = "Shift", y = "Country", subtitle = "", caption = "") +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_cyclical(
            values = c("#EC7063", "#F7DC6F", "#229954"), guide = "legend") +
        coord_cartesian(clip = "off") +
        theme_ridges(grid = TRUE) +
        theme(legend.position="top", 
        legend.title = element_blank())
        )
    dev.off()

    tabs[[model_number]] = tab
    return(tabs)
}