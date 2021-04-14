########################
# utilities
# author: sebastian daza
########################


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


# create data for comparison (predictive values)
createComparisonData = function(data, countries, years, cyears, dyears) {
    datalist = list()
    for (i in seq_along(countries)) {
        for (h in seq_along(cyears)) {
            values = paste0(countries[i], ".", cyears[[h]])
            if (sum(values %in% unique(idat$ctryear)) < 2 ) next 
            a = data[ctry == countries[i] & year == years[h] & gyear == cyears[[h]][2]]
            b = copy(a)
            b[, ctryear := paste0(countries[i], ".", cyears[[h]][1])]
            b[, ctry50 := paste0(countries[i], ".", dyears[[h]][1])]
            datalist[[paste0(i, ".", h)]] = rbind(a, b)
        }
    }
    return(rbindlist(datalist, idcol = "comp"))
}


# shifts with stacking for only one country (auxiliary function)
estimate_shift = function(models=NULL, # list
                          ps=NULL, # list, posterior samples
                          posterior_nsample=500,
                          data=NULL, # data.table
                          country=NULL, # string
                          obs_var = NULL,
                          weights= NULL, # vector, length = number of models
                          cfyear=NULL, # numeric
                          transform = TRUE, # weibull transform
                          segment=NULL, # string representing period, valid values 1950, 1950-1969, 1970-1989, 1990
                          model_pred = list('1930' = '1950', '1950' = '1950',
                                '1970' = '1950-1969', '1990' = '1970-1989', '2010' = '1990')
                         )  {

    output = list()

    # equal weights (average) if they are not specified
    if (is.null(weights)) { weights = rep(1/length(models), length(models)) }

    # loop through models
    for (i in seq_along(models)) {

        if (is.null(ps)) {  s = data.table(
                            posterior_samples(models[[i]], subset = sample(
                              1:nsamples(models[[i]]), posterior_nsample, replace=FALSE))) }
        else { s = data.table(ps[[i]]) }
        # counterfactual using previous coefficients and intercepts (random effects)
        if (is.null(segment)) { igyear = model_pred[[as.character(cfyear)]] }
        else { igyear = segment }
        # print(paste('segment: ', segment))
        # print(paste('igyear: ', igyear))
        # print(unlist(model_pred))

        if (!igyear %in% as.vector(unlist(model_pred))) { stop('Segment (period) is not valid!') }

        colnames = names(s)
        betas = grep('^b_', colnames, value=TRUE)
        random = str_subset(colnames, paste0('^r_.+\\[', country, '.', igyear, ','))
        coef = c(betas, random)
        # print(coef)

        variables = c('ctry', 'year', sub('b_', '', betas))
        variables = variables[variables != 'Intercept']
        covariates = sub('b_', '', betas)

        data[, Intercept := 1]
        dt = data[ctry==country & year==cfyear, ..covariates]
#         print(head(dt))

        if (transform) {
            max_le = unique(data[ctry==country & year==cfyear, max_le])
        }

        ex_obs = unique(data[ctry==country & year==cfyear, ..obs_var])[[1]]


        st = s[, ..coef] # select coefficients
#         print(coef)
#         print('counterfactual')
#         print(head(st))
        for (h in seq_along(covariates)) {
            st[, covariates[h] := rowSums(.SD), .SDcols = grep(covariates[h], names(st), value=TRUE)]
        }

        mt = as.matrix(st[, ..covariates])

        # using observed values!
        cf = mt %*% as.vector(as.matrix(dt)) # counterfactual

        if (transform) {
            cf = unlist(sapply(cf, function(x) get_orig_values_weibull(x, max_value=max_le)))
                               }
        else { cf = cf[,1] }

        # using predicted values!
        pigyear =  data[ctry==country & year == cfyear, gyear]
#             print(igyear)
#             print('prediction')
        prandom = str_subset(colnames, paste0('^r_.+\\[', country, '.', pigyear, ','))
        pcoef = c(betas, prandom)
#             print(coef)
        pt = s[, ..pcoef]
#             print(head(pt))
#             print(covariates)

        for (h in seq_along(covariates)) {
            pt[, covariates[h] := rowSums(.SD), .SDcols = grep(covariates[h], names(pt), value=TRUE)]
        }

#             print(dt)
        mpt = as.matrix(pt[, ..covariates])
#             print(head(mpt))
        ex_pred = mpt %*% as.vector(as.matrix(dt))
        if (transform) {
            ex_pred = unlist(sapply(ex_pred, function(x) get_orig_values_weibull(x, max_value=max_le)))
        } else { ex_pred = ex_pred[,1]}

        # using obs and pred for generality
        output[[i]] = data.table(obs=ex_obs, pred=ex_pred, counterfactual=cf)
        output[[i]][, shift_obs := obs - counterfactual]
        output[[i]][, shift_pred := pred - counterfactual]

#                     'prediction' = as.vector(as.matrix(setDT(predictions)) %*% weights),
#                     'counterfactual' = as.vector(as.matrix(setDT(cfs)) %*% weights)))
    }
    if (length(output)==1) {
        return(output[[1]])
    }
    else {
        for (i in seq_along(models)) {
            output[[i]] = output[[i]] * weights[i]
            output[[i]][, model := i]
            output[[i]][, id := 1:.N]
        }
        results = rbindlist(output, fill=FALSE)[, lapply(.SD, sum, na.rm = FALSE), by = id]
        results[, model := NULL]
        return(results)
    }
}


# estimate shift for each country
compute_shifts = function(models = NULL, # list of models
                        ps = NULL, # list of posterior samples
                        posterior_nsample = 500,
                        weights = NULL, # vector with model weigths
                        data = NULL,
                        obs_var = NULL,
                        transform = TRUE,
                        countries = NULL,
                        years = NULL,
                        model_pred = list('1930' = '1950', '1950' = '1950', '1970' = '1950-1969',
                          '1990' = '1970-1989', '2010' = '1990')) {


    # print(model_pred)
    # list to save results
    results = list()

    for (c in countries) {

        iyears = as.numeric(unique(data[ctry==c & year %in% years, year]))
        segments = as.character(unique(data[ctry==c, gyear]))

    for (ys in iyears) {

        for (seg in segments) {

        est = estimate_shift(models = models,
            ps = ps,
            posterior_nsample = posterior_nsample,
            weights = weights,
            data= data,
            obs_var = obs_var,
            transform = transform,
            country = c,
            cfyear = ys,
            segment = seg,
            model_pred = model_pred)

       name = paste0(c(c,ys,seg), collapse='.')

       results[[paste0(c(c,ys,seg), collapse='.')]] = est[, name := name]

       }
     }
    }

    results = rbindlist(results)
    results[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)][,
                                        num_models := length(models)]

    return(results[, .(ctry, year, segment, num_models, obs, pred,
        counterfactual, shift_obs, shift_pred)])

}

prediction_checks_ex = function(model, data, countries, variables, y, x, maxy = 2, transform = FALSE) {

    variables = unique(c(variables, y, x))
    posterior = data.table(predict(model))
    if (transform) {
        vars = c('Estimate', 'Q2.5', 'Q97.5')
        posterior[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]

    }
    pred = data.table(cbind(data[, ..variables], posterior))

    setnames(pred, c('Estimate', 'Q2.5', 'Q97.5'), c('m', 'lo', 'hi'))

    wsd = sd(data[[y]])
    max_ex = max(data[, y, with = FALSE]) + maxy * wsd
    min_ex = min(data[, y, with = FALSE]) - maxy * wsd
    max_year = max(data[, x, with = FALSE])
    min_year = min(data[, x, with = FALSE])

    plots = list()
    for (c in seq_along(countries)) {
        plots[[c]] = ggplot(pred[ctry == countries[c]], aes_string(x=x, y=y))+
            geom_line(aes(y = m), color='#2b8cbe', size = 0.4)  +
            geom_ribbon(aes(ymin = lo, ymax = hi), fill = '#a6bddb', alpha=0.2)  +
            geom_point(size=0.3, color='#e34a33', alpha=0.4) +
            labs(title = countries[c]) +
            ylim(min_ex, max_ex) +
            xlim(min_year, max_year) +
            theme_minimal() +
            geom_vline(xintercept = 1950, size=0.5, color='red', alpha=0.8, linetype = 'dotted') +
            geom_vline(xintercept = 1970, size=0.5, color='red', alpha=0.8, linetype = 'dotted') +
            geom_vline(xintercept = 1990, size=0.5, color='red', alpha=0.8, linetype = 'dotted')
    }
    return(plots)
}


prediction_checks_pp_ex = function(posterior, data, countries, variables, y, x, maxy = 2, transform = FALSE) {

    posterior = data.table(posterior)
    if (transform) {
        vars = c('Estimate', 'Q2.5', 'Q97.5')
        posterior[, (vars) := lapply(.SD, recoverWeibull, maxvalue = 78.6), .SDcols = vars]

    }
    variables = unique(c(variables, y, x))
    pred = cbind(data[, ..variables], posterior)
    setnames(pred, c('Estimate', 'Q2.5', 'Q97.5'), c('m', 'lo', 'hi'))

    wsd = sd(data[[y]])
    max_ex = max(data[, y, with = FALSE]) + maxy * wsd
    min_ex = min(data[, y, with = FALSE]) - maxy * wsd
    max_year = max(data[, x, with = FALSE])
    min_year = min(data[, x, with = FALSE])

    plots = list()
    for (c in seq_along(countries)) {
        plots[[c]] = ggplot(pred[ctry == countries[c]], aes_string(x=x, y=y))+
            geom_line(aes(y = m), color='#2b8cbe', size = 0.4)  +
            geom_ribbon(aes(ymin = lo, ymax = hi), fill = '#a6bddb', alpha=0.2)  +
            geom_point(size=0.3, color='#e34a33', alpha=0.4) +
            labs(title = countries[c]) +
            ylim(min_ex, max_ex) +
            xlim(min_year, max_year) +
            theme_minimal() +
            geom_vline(xintercept = 1950, size=0.5, color='red', alpha=0.8, linetype = 'dotted') +
            geom_vline(xintercept = 1970, size=0.5, color='red', alpha=0.8, linetype = 'dotted') +
            geom_vline(xintercept = 1990, size=0.5, color='red', alpha=0.8, linetype = 'dotted')
    }
    return(plots)
}



extract.brmsfit <- function (model,
                             use.HDI = TRUE,
                             level = 0.9,
                             include.random = TRUE,
                             include.rsquared = TRUE,
                             include.nobs = TRUE,
                             include.loo.ic = TRUE,
                             reloo = FALSE,
                             include.waic = TRUE,
                             ...) {
  sf <- summary(model, ...)$fixed
  coefnames <- rownames(sf)
  coefs <- sf[, 1]
  se <- sf[, 2]
  if (isTRUE(use.HDI)) {
    hdis <- coda::HPDinterval(brms::as.mcmc(model, combine_chains = TRUE),
                              prob = level)
    hdis <- hdis[seq(1:length(coefnames)), ]
    ci.low = hdis[, "lower"]
    ci.up = hdis[, "upper"]
  } else { # default using 95% posterior quantiles from summary.brmsfit
    ci.low = sf[, 3]
    ci.up = sf[, 4]
  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (isTRUE(include.random) & isFALSE(!nrow(model$ranef))) {
    sr <- summary(model, ...)$random
    sd.names <- character()
    sd.values <- numeric()
    for (i in 1:length(sr)) {
      sd <- sr[[i]][, 1]
      sd.names <- c(sd.names, paste0("SD: ", names(sr)[[i]], names(sd)))
      sd.values <- c(sd.values, sd)
    }
    gof <- c(gof, sd.values)
    gof.names <- c(gof.names, sd.names)
    gof.decimal <- c(gof.decimal, rep(TRUE, length(sd.values)))
  }
  if (isTRUE(include.rsquared)) {
    rs <- brms::bayes_R2(model)[1]
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.nobs)) {
    n <- stats::nobs(model)
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (isTRUE(include.loo.ic)) {
    looic <- brms::loo(model, reloo = reloo)$estimates["looic", "Estimate"]
    gof <- c(gof, looic)
    gof.names <- c(gof.names, "loo IC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.waic)) {
    waic <- brms::waic(model)$estimates["waic", "Estimate"]
    gof <- c(gof, waic)
    gof.names <- c(gof.names, "WAIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(coef.names = coefnames,
                     coef = coefs,
                     se = se,
                     ci.low = ci.low,
                     ci.up = ci.up,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}

setMethod("extract",
          signature = className("brmsfit_multiple", "brms"),
          definition = extract.brmsfit)