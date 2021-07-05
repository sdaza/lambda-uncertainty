runModel = function(flist, samples, iterations = 2000, chains = 1, 
    model_number = 1, ex_max, country_labs, newdata, plots_path, 
    select_estimates = "") {

    tabs = list()
    idat = samples[[1]]
    nsamples = length(samples)

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
        print(":::::::::::: running stan models")
        model = ulam(
            flist, 
            data = mdata, chains = chains, cores = 1, 
            iter = iterations, chain_id = i
        )
        
        check = as.matrix(precis(model, depth = 3))
        print(paste0("::::::::: There are ", 
            sum(check[, "Rhat4"] > 1.01), 
            " parameters with Rhat4 > 1.01"))
            
        results$models = model@stanfit
        results$shifts = computeShift(model, newdata, ex_max)
        return(results)
    }

    # extract restuls
    shifts = list()
    models = list()
    for (i in seq_along(output)) {
        shifts[[i]] = output[[i]][["shifts"]]
        models[[i]] = output[[i]][["models"]]
    }

    print(":::::::::::: creating tabs")
    shifts = rbindlist(shifts)
    fit = sflist2stanfit(models)
    tab = extractStan(fit, n = nrow(idat))

    # plots

    print(":::::::::::: creating plots")
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
    return(tab)
}




# model function 
createModel = function(f, samples, iterations = 20000, data_for_shifts, 
    countries, country_labs, number_clusters = 20, bprior = NULL, 
    model_number = NULL, vars, plots_path, selected_estimates = "t", 
    manus_plots, vmean = NULL, vsd = NULL) {

    nsamples = length(samples)
    # cl = makeCluster(20, outfile="")
    cl = makeCluster(20)
    registerDoParallel(cl)

    print("::::::::::: running models")
    output = foreach(i = 1:nsamples) %dopar% {
        library(data.table)
        library(brms)
        print(paste0(":::::::: Running iteration ", i, " ::::::::"))
        dat = samples[[i]]
        dat[, zwy := (wy - vmean) / vsd]
        model = brms::brm(f, data = dat, 
            iter = iterations, 
            # warmup = 1000, 
            chains = 4, 
            refresh = 5000,
            prior = bprior
        )
        return(model)
    }
    stopCluster(cl)

    # combine models
    model = combine_models(mlist = output, check_data = FALSE)
    idat = samples[[1]]
    idat[, zwy := (wy - vmean) / vsd]
    countries = unique(data_for_shifts$ctry)

    # fit plots
    print("::::::::::: predicting")
    pred = predict(model, newdata = idat, 
        summary = TRUE, nsamples = 10000)
    plots_checks = prediction_check_plots(pred, idat, countries, 
        vars,  y = "ex_mean", x = "year", 
        xlab = "\nYear", ylab = "e0\n", 
        transform = TRUE, country_labels = country_labs, 
        vmean = vmean, vsd = vsd)
    plots = wrap_plots(plots_checks, ncol = 3)
    savepdf(paste0(plots_path, select_estimates, 
        paste0("fit_check_m", model_number)), width = 30, height = 35)
        print(plots)
    dev.off()
    file.copy(paste0(plots_path, select_estimates, 
        paste0("fit_check_m", model_number, ".pdf")), manus_plots, 
        recursive = TRUE)  

    # shift plots 
    print("::::::::::: computing shifts")
    mshifts = list()
    for (i in unique(data_for_shifts$year)) {
        temp = createShifts(model, data_for_shifts[year == i], 
            countries = unique(data_for_shifts[year == i]$ctry), 
            nsamples = 5000, K = K, vmean = vmean, vsd = vsd)
        temp = data.table(temp$values)
        if(length(temp) == 0) { print("No shift data!") }
        mshifts[[as.character(i)]] = temp[, year := i]
    }
    shifts = rbindlist(mshifts, fill = TRUE)
    shifts = melt(shifts, id.vars = "year", 
        variable.name = "ctry", 
        value.name = "shift")
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
    file.copy(paste0(plots_path, select_estimates, 
        paste0("shifts_by_period_m", model_number, ".pdf")), manus_plots, 
        recursive = TRUE)  
    print("::::::::::: creating tab")
    tab = extractBRMS(model, r2 = TRUE, weight = NULL)
    slackr::slackr_msg(txt = paste0("Script done at: ", Sys.time()))
    return(list("model" = model, "shifts" = shifts, "tab" = tab))
}

    # bayes_R2 = function(mu, sigma) {
    #     mu = transpose(data.table(mu))
    #     sigma = data.tab
    #     var_mu = as.vector(t(mu[, lapply(.SD, var)]))
    #     sigma2 <- sigma^2
    #     r2 = var_mu / (var_mu + sigma2)
    #     median(r2)
    # }
prediction_check_plots = function(posterior, data, countries, variables, y, x, 
    maxy = 2, transform = FALSE, country_labels, xlab = NULL, 
    ylab = NULL, vmean = NULL, vsd = NULL) {

    posterior = data.table(posterior)
    if (transform) {
        vars = c('Estimate', 'Q2.5', 'Q97.5')
        if (!is.null(vmean)) {
            posterior[, (vars) := lapply(.SD, recoverZ, vmean, vsd), .SDcols = vars] 
        }
        posterior[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]

    }
    variables = unique(c(variables, y, x))
    pred = cbind(data[, ..variables], posterior)
    setnames(pred, c('Estimate', 'Q2.5', 'Q97.5'), c('m', 'lo', 'hi'))

    # setup limitis of plots
    wsd = sd(data[[y]])
    max_ex = max(data[, y, with = FALSE]) + maxy * wsd
    min_ex = min(data[, y, with = FALSE]) - maxy * wsd
    max_year = max(data[, x, with = FALSE])
    min_year = min(data[, x, with = FALSE])

    plots = list()
    for (c in seq_along(countries)) {
        plots[[c]] = ggplot(pred[ctry == countries[c]], aes_string(x=x, y=y))+
            geom_line(aes(y = m), color='#2b8cbe', size = 0.4)  +
            geom_ribbon(aes(ymin = lo, ymax = hi), fill = '#a6bddb', alpha=0.2) +
            geom_point(size=0.3, color='#e34a33', alpha=0.4) +
            labs(title = country_labels[[as.character(countries[c])]], 
                x = ifelse(is.null(xlab), x, xlab),
                y = ifelse(is.null(ylab), y, ylab)) +
            ylim(min_ex, max_ex) +
            xlim(min_year, max_year) +
            theme_minimal() +
            geom_vline(xintercept = 1950, size=0.5, color='red', 
                alpha=0.8, linetype = 'dotted') +
            geom_vline(xintercept = 1970, size=0.5, color='red', 
                alpha=0.8, linetype = 'dotted') +
            geom_vline(xintercept = 1990, size=0.5, color='red', 
                alpha=0.8, linetype = 'dotted')
    }
    return(plots)
}



predictionData = function(model, data, pars = "pred", ex_max, prob = 0.95) {
    
    pred = data.table(recoverWeibull(
        as.matrix(rstan::extract(model, permuted = TRUE, pars = "pred"))[[1]], 
            ex_max))

    m = as.vector(t(pred[, lapply(.SD, mean)]))
    lo = as.vector(t(pred[, lapply(.SD, quantile, prob = (1 - prob) / 2)]))
    hi = as.vector(t(pred[, lapply(.SD, quantile, prob = (prob + (1-0.95)/2))]))
    data[, m := m]
    data[, lo := lo]
    data[, hi := hi]
    return(data)
}
