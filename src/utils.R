# lambda project
# author: sebastian daza


# functions

# get original values from weibull transformation
get_orig_values_weibull = function(x, max_value) {
    return ( (1 - exp(-exp(x))) * max_value )
}

# shifts with stacking for only one country (auxiliary function)
estimate_shift = function(models=NULL, # list
                          ps=NULL, # list, posterior samples
                          posterior_nsample=500,
                          data=NULL, # data.table
                          country=NULL, # string
                          obs_var = NULL,
                          weights= NULL, # vector, lenght = number of models
                          cfyear=NULL, # numeric
                          transform = TRUE, # weibull transform
                          segment=NULL, # string representing period, valid values 1950, 1950-1969, 1970-1989, 1990
                          model_pred = list('1930' = '1950', '1950' = '1950', '1970' = '1950-1969', '1990' = '1970-1989', '2010' = '1990')
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

    for (c in countries ) {

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

    return(results[, .(ctry, year, segment, num_models, obs, pred, counterfactual, shift_obs, shift_pred)])

}




# lags with stacking for only one country (auxiliary function)
estimate_lag = function(models=NULL, # list
                          ps=NULL, # list, posterior samples
                          data=NULL, # data.table
                          country=NULL, # string
                          weights= NULL, # vector, lenght = number of models
                          cfyear=NULL, # numeric
                          segment=NULL, # string representing period, valid values 1950, 1950-1969, 1970-1989, 1990
                          predicted_values=FALSE # boolean
                         )  {


    setorder(data, year)

    year_values = data[ctry==country, year]
    ex_values = data[ctry==country, Ex]

    if(!(length(year_values) == length(ex_values))) { stop('LE values should have same lenght as years')}

    differences = list()

    model_pred = list('1950' = '1950', '1970' = '1950-1969', '1990' = '1970-1989', '2010' = '1990')

    # equal weights (average) if they are not specified
    if (is.null(weights)) { weights = rep(1/length(models), length(models)) }

    # loop through models
    for (i in seq_along(models)) {

        if (is.null(ps)) {  s = data.table(posterior_samples(models[[i]])) }
        else {  s = data.table(ps[[i]]) }

        # counterfactual using previous coefficients and intercepts (random effects)
        if (is.null(segment)) { igyear = model_pred[as.character(cfyear)] }
        else { igyear = segment }


        if (!igyear %in% as.vector(unlist(model_pred))) { stop('Segment (period) is not valid!') }

        colnames = names(s)
        betas = grep('^b_', colnames, value=TRUE)
        random = str_subset(colnames, paste0('^r_.+\\[', country, '.', igyear, ','))
        coef = c(betas, random)
#         print(coef)

        variables = c('ctry', 'year', sub('b_', '', betas))
        variables = variables[variables != 'Intercept']
        covariates = sub('b_', '', betas)

        data[, Intercept := 1]
        dt  = data[ctry==country & year==cfyear, ..covariates]
#         print(head(dt))

        max_le = unique(data[ctry==country & year==cfyear, max_le])
        ex_obs = unique(data[ctry==country & year==cfyear, Ex])

#     print(ex_obs)

        st = s[, ..coef] # select coefficients
#         print(coef)
#         print('counterfactual')
#         print(head(st))
        for (h in seq_along(covariates)) {
            st[, covariates[h] := rowSums(.SD), .SDcols = grep(covariates[h], names(st), value=TRUE)]
        }

        mt = as.matrix(st[, ..covariates])

#         print(head(mt))

        cf = mt %*% as.vector(as.matrix(dt)) # counterfactual
        cf = unlist(sapply(cf, function(x) get_orig_values_weibull(x, max_value=max_le)))

#         print(head(cf))

        if (predicted_values) { # using predicted values (all random effects) instead of observed Ex

           ex_values = predict(models[[i]], data[ctry==country], summary=FALSE)
           ex_values = apply(ex_values, 2, mean)
           ex_values  = unlist(lapply(ex_values,  function(x)  get_orig_values_weibull(x, max_le)))
#            print(ex_values)

#            print(length(year_values))
#            print(length(cf))

           ind = NULL
           for (j in 1:length(cf)) {
             ind[j] = which.min(abs(ex_values - cf[j]))
           }

           differences[[i]]  = year_values[ind] - cfyear

        } else {

           ind = NULL
           for (j in 1:length(cf)) {
             ind[j] = which.min(abs(ex_values - cf[j]))
           }

           differences[[i]] = year_values[ind] - cfyear
        } # using observed Ex values
#             print(head(differences[[i]]))
      }

    # combine values (differences) using weights
    return( as.vector(as.matrix(setDT(differences)) %*% weights) ) # return a vector
}

# compute lags for each country

compute_lags = function(models = NULL, # list of models
                        ps = NULL, # list of posterior samples
                        weights = NULL, # vector with model weigths
                        data = NULL,
                        countries = NULL,
                        years = NULL,
                        predicted_values=FALSE,
                        model_pred = list('1950' = '1950', '1970' = '1950-1969', '1990' = '1970-1989', '2010' = '1990')) {



    # list to save results
    lags = list()

    for (c in countries ) {

        iyears = as.numeric(unique(data[ctry==c & year %in% years, year]))
        segments = as.character(unique(data[ctry==c, gyear]))

    for (ys in iyears) {

        for (seg in segments) {

        est = estimate_lag(models = models,
            ps = ps,
            weights = weights,
            data= data,
            country = c,
            cfyear = ys,
            segment = seg,
            predicted_values=predicted_values)

        name = paste0(c(c,ys,seg), collapse='.')
        lags[[paste0(c(c,ys,seg), collapse='.')]] = data.table(name, pred_lag = est)

        }
    }
    }

    lags = rbindlist(lags)
    lags[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)][,
                                        num_models := length(models)]
    lags = lags[, .(ctry, year, segment, num_models, pred_lag)]

    return(lags)

}

# texreg function for brms

extract.brms = function(model, include.r2 = TRUE, include.loo = FALSE, ...) {

  s = summary(model)

  # fixed
  coefficient.names = names(s$fixed[,1])
  coefficients = s$fixed[,1]
  standard.errors = s$fixed[, 2]

  # random
  if ('random' %in% names(s)) {
   r = s$random[[1]]
   random.names = stringr::str_replace_all(rownames(r), stringr::fixed('\\_'), "\\_")
   random.estimates = r[,1]
   random.se = r[,2]

   coefficient.names = c(coefficient.names, random.names)
   coefficients = c(coefficients, random.estimates)
   standard.errors = c(standard.errors, random.se)

  }

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  gof = c(gof, s$nobs)
  gof.names <- c(gof.names, "Num.\ obs.")
  gof.decimal <- c(gof.decimal, FALSE)

  if ('ngrps' %in% names(s)) {
    for (i in seq_along(s$ngrps)) {
          gof = c(gof, s$ngrps[[i]])
          gof.names <- c(gof.names, paste('Num.\ obs. ',
            stringr::str_replace_all(names(s$ngrps[i]), stringr::fixed('_'), '\\_')))
          gof.decimal <- c(gof.decimal, FALSE)
    }
  }

  if (include.loo == TRUE) {
    loo = loo::loo(model, reloo=TRUE)
    gof <- c(gof, loo$estimates[3,1])
    gof.names <- c(gof.names, "LOO Information Criterion")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  if (include.r2 == TRUE) {
    r2 = bayes_R2(model)
    gof <- c(gof, round(r2[1, 1], 2))
    gof.names <- c(gof.names, "Bayes $R^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }


  tr <- createTexreg(
      coef.names = coefficient.names,
      coef = coefficients,
      se = standard.errors,
      gof.names = gof.names,
      gof = gof,
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("brmsfit", "brms"),
    definition = extract.brms)
