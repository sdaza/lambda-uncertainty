########################
# utilities
# author: sebastian daza
########################


library(texreg)


# get original values from weibull transformation
transWeibull = function(x, maxvalue = 100) {
    y = x / (maxvalue + 1.05)
    y = log( -log(1-y))
    return(y)
}


recoverWeibull = function(x, maxvalue = 100) {
    y = (1 - exp(-exp(x))) * (maxvalue + 1.05)
    return(y)
}


getSample = function(x) {
    if (length(x) > 1) {
        s = sample(x, 1)
    } else {
        s = x
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


prediction_check_plots = function(posterior, data, countries, variables, y, x, 
    maxy = 2, transform = FALSE, country_labels) {

    posterior = data.table(posterior)
    if (transform) {
        vars = c('Estimate', 'Q2.5', 'Q97.5')
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
            labs(title = country_labels[[as.character(countries[c])]]) +
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


# shift functions
# create data for comparison (predictive values)
createComparisonData = function(data, countries, years, cyears, dyears) {
    datalist = list()
    for (i in seq_along(countries)) {
        for (h in seq_along(cyears)) {
            values = paste0(countries[i], ".", cyears[[h]])
            if (sum(values %in% unique(idat$ctryear)) < 2 ) next 
            a = data[ctry == countries[i] & year == years[h] & 
                gyear == cyears[[h]][2]]
            b = copy(a)
            b[, ctryear := paste0(countries[i], ".", cyears[[h]][1])]
            b[, ctry50 := paste0(countries[i], ".", dyears[[h]][1])]
            datalist[[paste0(i, ".", h)]] = rbind(a, b)
        }
    }
    return(rbindlist(datalist, idcol = "comp"))
}


computeShift = function(predictions, countries) {
    output = list()
    index = c(1, 2)    
    for (i in seq_along(countries)) {
        output[[as.character(countries[i])]] = predictions[[paste0("V", index[1])]] - predictions[[paste0("V", index[2])]]
        index = index + 2
    }
    return(output)
}


createShifts = function(models, newdata, nsamples = 1000, countries = NULL, 
    transform = TRUE, K = 10) {
    
    lpd = NULL
    for (h in seq_along(models)) {
        print(paste0("::::::::: model ", h, " :::::::::"))
        cv = suppressWarnings(suppressMessages(suppressPackageStartupMessages(
                kfold(models[[h]], K = K, chains = 1, cores = 2))))
            lpd = cbind(lpd, cv$pointwise[, "elpd_kfold"])
    }
    
    weights= as.vector(stacking_weights(lpd))
    pred = rlang::invoke(pp_average, models, weights = weights, 
            newdata = newdata, 
            summary = FALSE, nsamples = nsamples)
    pred = data.table(pred)
    if (transform) {
        vars = names(pred)
        pred[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]
    }
    shifts = setDT(computeShift(pred, countries))
    return(list("values" = shifts, "weights" = weights)) 
}



# createShifts = function(model_replicates, newdata, nsamples = 1000, countries = NULL, 
#     transform = TRUE) {
    
#     models = length(model_replicates)
#     replicates = length(model_replicates[[1]])

#     # list to save output
#     weight_list = list()
#     shift_estimates = list()
    
#     for (i in 1:replicates) {

#         print(paste0("::::::::: replicate ", i, " :::::::::"))

#         cv_list = list()
#         lpd = NULL
        
#         # temporary model list
#         rmodels = list()

#         for (h in 1:models) {
#             print(paste0("::::::::: model ", h, " :::::::::"))
#             cv_list[[h]] = suppressWarnings(suppressMessages(suppressPackageStartupMessages(
#                 kfold(model_replicates[[h]][[i]], K = 10, chains = 1))))
#             lpd = cbind(lpd, cv_list[[h]]$pointwise[, "elpd_kfold"])

#             # cv_list[[h]] = brms::loo(rmodels[[h]], moment_match = TRUE, reloo = FALSE, 
#                 # cores = 20, reloo_args = list(chains = 1))
#             # lpd = cbind(lpd, cv_list[[h]]$pointwise[, "elpd_loo"])
#         }


#         weight_list[[i]]= as.vector(stacking_weights(lpd))
#         ndata = newdata[[i]]
#         pred = rlang::invoke(pp_average, 
#             getModels(model_replicates, models, i), weights = weight_list[[i]], newdata = ndata, 
#             summary = FALSE, nsamples = nsamples)
#         pred = data.table(pred)

#         if (transform) {
#             vars = names(pred)
#             pred[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]
#         }

#         shift_estimates[[i]] = setDT(computeShift(pred, countries))

#     }
#     shifts = rbindlist(shift_estimates, idcol = "replicate")
#     avg_weights = apply(do.call(rbind, weight_list), 2, mean)
    
#     return(list("shifts" = shifts, "weight_list" = weight_list, "avg_weights" = avg_weights))
       
# }


# imputation function
parmice = function(data, n.core = detectCores() - 1, n.imp.core = 2,
    seed = NULL, m = NULL, ...) {
        
    suppressMessages(require(parallel))
    cl = makeCluster(n.core, ...)
    clusterExport(cl, varlist = "data", envir = environment())
    clusterEvalQ(cl, library(miceadds)) # to use miceadds!
    if (!is.null(seed)) {
        clusterSetRNGStream(cl, seed)
    }
    if (!is.null(m)) {
        n.imp.core = ceiling(m / n.core)
    }
    imps = parLapply(cl = cl, X = 1:n.core, fun = function(i) {
        mice(data, print = FALSE, m = n.imp.core, ...)
        }
    )
    stopCluster(cl)
    imp = imps[[1]]
    if (length(imps) > 1) {
        for (i in 2:length(imps)) {
            imp = ibind(imp, imps[[i]])
        }
    } 
    return(imp)
}


extractBRMS = function(model, r2 = TRUE) {
    sf = fixef(model)
    coefnames = rownames(sf)
    coefs = sf[, 1]
    se = sf[, 2]
    ci.low = sf[, 3]
    ci.up = sf[, 4]

    gof = numeric()
    gof.names = character()
    gof.decimal = logical()

    n = stats::nobs(model)
    gof = c(gof, n)
    gof.names = c(gof.names, "Num. obs.")
    gof.decimal = c(gof.decimal, FALSE)

    if (r2) {
        rs = brms::bayes_R2(model)[1]
        gof = c(gof, rs)
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