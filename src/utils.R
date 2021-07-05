########################
# utilities
# author: sebastian daza
########################


# decimals
specify_decimal = function(x, k) trimws(format(round(x, k), nsmall=k))

# get original values from weibull transformation
transWeibull = function(x, maxvalue = 100) {
    y = x / (maxvalue + 1.05)
    y = log( -log(1-y))
    return(y)
}


recoverZ = function(x, vmean, vsd) {
    return((x * vsd) + vmean)
}


recoverWeibull = function(x, maxvalue = 100) {
    y = (1 - exp(-exp(x))) * (maxvalue + 1.05)
    return(y)
}


getSample = function(x) {
    x = na.omit(as.numeric(x))
    if (length(x) > 0) {
        s = sample(x, 1)
    } else {
        s = as.numeric(NA)
    }
    return(s)
}


table = function (...) base::table(..., useNA = 'ifany')


cor = function (...) stats::cor(..., use = "complete.obs")


nlocf = function(x) {
    which.na = c(which(!is.na(x)), length(x) + 1)
    values = na.omit(x)
    if (which.na[1] != 1) {
        which.na = c(1, which.na)
        values = c(values[1], values)
    }
    diffs = diff(which.na)
    return(rep(values, times = diffs))
}


lookvar = function(dat, varnames) {
    n = names(dat)
    nn = list()
    for (i in 1:length(varnames)) {
        nn[[i]] = grep(varnames[i],n)
    }

    nn = unlist(nn)

    if ( length(nn) >0 )
    {
        r = n[nn]
        return(r)
    }
    else { return("No variables found") }
}


countmis = function(dat, vars = NULL, pct = TRUE, exclude.complete = TRUE) {

    if (is.null(vars)) {
        vars = names(dat)
    }

    mis = sort(sapply(dat[, vars, with = FALSE],
        function(x) sum(is.na(x))), decreasing = TRUE)

    if (exclude.complete == TRUE) {
         mis = mis[mis > 0]
    }

    if (pct == FALSE)
        { return(mis) }
    else if ( pct == TRUE ) {
        return( round(mis / nrow(dat), 3))
    }
    return(mis)
}


getMax = function(x) {
    x = na.omit(x)
    if (length(x) == 0) {
        return(NA_real_)
    } else {
        return(max(x))
    }
}


getMin = function(x) {
    x = na.omit(x)
    if (length(x) == 0) {
        return(NA_real_)
    } else {
        return(min(x))
  }
}


savepdf = function(file, width = 16, height = 10, mgp = c(2.2,0.45,0),
    tcl = -0.4, mar = c(3.3,3.6,1.1,1.1)) {
    fname = paste0(file, ".pdf")
    pdf(fname, width=width / 2.54, height = height / 2.54,
        pointsize = 10)
    par(mgp = mgp, tcl = tcl, mar = mar)
}


predictionPlots = function(pred, x = "year", y = "ex_mean", country_labels, 
    xlab = "Year", ylab = "e0") {
    
    countries = unique(pred$ctry)
    wsd = sd(pred[[y]])/ max(pred[, y, with = FALSE]) 
    max_y = max(pred[, y, with = FALSE]) + max(pred[, y, with = FALSE]) * wsd
    min_y = min(pred[, y, with = FALSE]) - min(pred[, y, with = FALSE]) * wsd
    max_x = max(pred[, x, with = FALSE])
    min_x = min(pred[, x, with = FALSE])

    plots = list()
    
    for (c in seq_along(countries)) {
        plots[[c]] = ggplot(pred[ctry == countries[c]], aes_string(x=x, y=y))+
            geom_line(aes(y = m), color='#2b8cbe', size = 0.4)  +
            geom_ribbon(aes(ymin = lo, ymax = hi), fill = '#a6bddb', alpha=0.2) +
            geom_point(size=0.3, color='#e34a33', alpha=0.4) +
            labs(title = country_labels[[as.character(countries[c])]], 
                x = ifelse(is.null(xlab), x, xlab),
                y = ifelse(is.null(ylab), y, ylab)) +
            ylim(min_y, max_y) +
            xlim(min_x, max_x) +
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


predictData = function(model, ex_max = 100, n = 1000) {
    dt = recoverWeibull(rethinking::sim(model), ex_max)
    return(data.table(dt))
}


createComparisonData = function(data, countries, years, cyears, dyears) {
    datalist = list()

    for (i in seq_along(countries)) {
        for (h in seq_along(cyears)) {
            values = paste0(countries[i], ".", cyears[[h]])
            if (sum(values %in% unique(data$ctryear)) < 2 ) { next }
            else { 
                a = data[ctry == countries[i] & year == years[h]]
                b = copy(a)
                b[, qyear := years[[h]] - 5]
                b[, ctryear := paste0(countries[i], ".", cyears[[h]][1])]
                b[, ctry50 := paste0(countries[i], ".", dyears[[h]][1])]
                b[, year1950 := dyears[[h]][1]]
                datalist[[paste0(i, ".", h)]] = rbind(a, b)
                rm(a, b)
            }
        }
    }
    df = rbindlist(datalist, idcol = "comp")
    df[, ctryearg := as.numeric(ctryear)]
    return(df)
}


iShift = function(predictions, countries) {
    output = list()
    index = c(1, 2)    
    for (i in seq_along(countries)) {
        output[[as.character(countries[i])]] = 
        predictions[[paste0("V", index[1])]] - predictions[[paste0("V", index[2])]]
        index = index + 2
    }
    return(setDT(output))
}   


computeShift = function(model, newdata, ex_max, n = 1000) {
    years = unique(newdata$year)
    mshifts = list()    
    for (i in years) {
        temp = data.table::copy(newdata[year == i])
        pred = data.table(recoverWeibull(
            rethinking::link(model, data = temp, n = n), ex_max))
        countries = unique(temp$ctry)
        shifts = iShift(pred, countries)
        mshifts[[as.character(i)]] = shifts[, year := i]
    }
    shifts = rbindlist(mshifts, fill = TRUE)
    shifts = melt(shifts, id.vars = "year", 
        variable.name = "ctry", 
        value.name = "shift")
    shifts = na.omit(shifts, cols="shift")
    return(shifts)
}


plotShifts = function(shifts, country_labs, title = "", xlab = "Shift", ylab = "Country") {
    xmin = min(shifts$shift, na.rm = TRUE) - 1.5
    xmax = max(shifts$shift, na.rm = TRUE) + 1.5
    v = unlist(country_labs)
    shifts[, lctry := as.factor(v[as.character(ctry)])]
    shifts[, year := factor(year)]
    shifts[, ctry := factor(ctry)]
    plot = ggplot2::ggplot(shifts, aes(y = lctry)) +
        ggridges::geom_density_ridges(aes(x = shift, fill = year), 
            alpha = .45, color = "white", from = xmin, to = xmax, scale = 1) +
        labs(x = xlab, y = ylab, title = title) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        ggridges::scale_fill_cyclical(
            values = c("#EC7063", "#F7DC6F", "#229954"), guide = "legend") +
        coord_cartesian(clip = "off") +
        ggridges::theme_ridges(grid = TRUE) +
        theme(legend.position="top", 
        legend.title = element_blank())
    return(plot)
}


multiResultClass = function(models = NULL, shifts = NULL, predictions = NULL, 
    r2 = NULL, rhat = NULL, neff = NULL) {
  me = list(models = models, shifts = shifts, predictions = predictions, 
    r2 = r2, rhat = rhat, neff = neff)
  class(me) = append(class(me), "multiResultClass")
  return(me)
}


runModel = function(flist, samples, chains = 1, iterations = 2000, 
    n = 100, ex_max = 100, newdata, clusters = 2) {

    nsamples = length(samples)

    cl = makeCluster(clusters)
    registerDoParallel(cl)
    
    output = foreach(i = 1:nsamples) %dopar% {
    
        library(data.table)
        library(rethinking)
        source("src/utils.R")
        results = multiResultClass()

        dat = samples[[i]]
        mdata = list(
                wy = dat$wy, 
                zyear = dat$zyear,
                zinfrastructure = dat$zinfrastructure, 
                zpop= dat$zpop, 
                zilit = dat$zilit, 
                zius_aid_pc = dat$zius_aid_pc, 
                zigdp_pc = dat$zigdp_pc, 
                ctryearg = dat$ctryearg
        )

        model = ulam(
            flist, 
            data = mdata, chains = chains, cores = 1, 
            iter = iterations, chain_id = i
        )

        check = as.matrix(rethinking::precis(model, depth = 3))
        results$rhat = sum(na.omit(check[, "Rhat4"]) > 1.01)
        results$neff = sum(na.omit(check[, "n_eff"]) < 100)
        rm(check)

        results$models = model@stanfit
        results$predictions = predictData(model, ex_max = ex_max, n = 100)
        results$shifts = computeShift(model, newdata, ex_max, n = 100)
        fit_ss = rstan::extract(model@stanfit, pars = c("pred", "sigma"))
        results$r2 = bayes_R2(fit_ss$pred, fit_ss$sigma)
        rm(fit_ss, model)
        return(results)
    }

    stopCluster(cl)
    
    # extract results
    models = list()
    predictions = list()
    shifts = list()
    r2 = NULL
    rhat = NULL
    neff = NULL
    for (i in seq_along(output)) {
        shifts[[i]] = output[[i]][["shifts"]]
        models[[i]] = output[[i]][["models"]]
        predictions[[i]] = output[[i]][["predictions"]]
        r2 = c(r2, output[[i]][["r2"]])
        rhat = c(rhat, output[[i]][["rhat"]])
        neff = c(neff, output[[i]][["neff"]])
    }
    
    predictions = rbindlist(predictions)
    pred = samples[[1]]
    pred[, m := apply(predictions, 2, median)]
    pred[, lo := apply(predictions, 2, quantile, probs = 0.025)]
    pred[, hi := apply(predictions, 2, quantile, probs = 0.975)]
    rm(predictions)
    
    return(
        list(
            "fit" = sflist2stanfit(models), 
            "shifts" = rbindlist(shifts), 
            "predictions" = pred,
            "r2" = median(r2),
            "rhat" = sum(rhat), 
            "neff" = sum(neff)
        )
    )
}


bayes_R2 = function(mu, sigma) {
    mu = transpose(data.table(mu))
    var_mu = as.vector(t(mu[, lapply(.SD, var)]))
    sigma2 = sigma^2
    r2 = var_mu / (var_mu + sigma2)
    mean(r2)
}


extractStan = function(model, r2 = TRUE, n = NULL) {

    t = as.matrix(rethinking::precis(model, depth = 3, prob = 0.95))
    t  = t[-grep("pred|lp__", rownames(t)),] 
    coefnames = rownames(t)
    coefs = t[, 1]
    se = t[, 2]
    ci.low = t[, 3]
    ci.up = t[, 4]

    gof = numeric()
    gof.names = character()
    gof.decimal = logical()

    if (!is.null(n)) {
        gof = c(gof, n)
        gof.names = c(gof.names, "Num. obs.")
        gof.decimal = c(gof.decimal, FALSE)
    }

    if (!is.null(r2)) {
        # print(":::::::: extracting pred and sigma")
        # fit_ss = rstan::extract(model, permuted = TRUE, pars = c("pred", "sigma"))
        # print(":::::::: computing R2")
        # rs = bayes_R2(fit_ss$pred, fit_ss$sigma)
        gof = c(gof, r2)
        gof.names = c(gof.names, "R$^2$")
        gof.decimal = c(gof.decimal, TRUE)
    }

    tr = texreg::createTexreg(
        coef.names = coefnames,
        coef = coefs,
        se = se,
        ci.low = ci.low,
        ci.up = ci.up,
        gof.names = gof.names,
        gof = gof,
        gof.decimal = gof.decimal)
    return(tr)
}