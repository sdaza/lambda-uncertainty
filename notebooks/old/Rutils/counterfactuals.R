# recover transformed values
get_orig_values_weibull = function(x, max_value) { 
    return ( (1 - exp(-exp(x))) * max_value )
}

# estimate shifts 

estimate_shift = function(samples,
                          gdp_value=NULL,
                          le_value=NULL,
                          country=NULL, 
                          model_year=NULL,
                          year=NULL,
                          max_le = NULL,
                          coefficients = c('Intercept', 'gdp_log')) {
    
    colnames = names(samples)
                          
    betas = paste0('b_', coefficients)
    
    random = str_subset(colnames, paste0('^r_.+\\[', country, '.', model_year, ','))
    s = samples[, c(betas, random)]  
    
    pred = (s[,1] + s[,3]) + (s[,2] + s[,4]) * gdp_value
    pred = unlist(sapply(pred, function(x) get_orig_values_weibull(x, max_value=max_le)))
    
    return(le_value-pred)           
}

# compute shifts over years and countries 
                         
compute_shifts = function(model, 
                          data, 
                          country_labels,
                          iyears, 
                          max_le,
                          coefficients,
                          predicted_le=FALSE) {
        
    samples = posterior_samples(model)
    colnames = names(samples)

    model_pred = list('1950' = '1950-1969', '1970' = '1970-1989', '1990' = '1990', '2010'= '1990')
    
    shifts = list()

    for (c in country_labels) {
        
        years = as.numeric(data[ctry==c & year %in% iyears, year])
        segments = as.character(unique(data[ctry==c, gyear])) 
        max_le_value = max_le[ctry==c, max_le]
        
        for (ys in years) {
            
             fsegments = segments[-which(segments == model_pred[[as.character(ys)]])]
            
             if (predicted_le) {
                 
                 betas = paste0('b_', coefficients)
                 my = model_pred[[as.character(ys)]]
     
                 random = str_subset(colnames, paste0('^r_.+\\[', c, '.', my, ','))
                 s = samples[, c(betas, random)]
    
                 gdp_log = data[ctry==c & year==ys, gdp_log]
                 
                 le_value = (s[,1] + s[,3]) + (s[,2] + s[,4]) * gdp_log[1]
                 
                 le_value = unlist(sapply(le_value, function(x) 
                     get_orig_values_weibull(x, max_value=max_le_value)))
            } 
            else {
                le_value = data[ctry==c & year==ys, le][1]
            }
                                          
            for (ysm in fsegments) {
               
               est = estimate_shift(samples=samples, 
                   le_value=le_value, 
                   gdp_value = data[ctry==c & year==ys, gdp_log],
                   country=c, 
                   year=ys, 
                   model_year=ysm,
                   max_le=max_le_value, 
                   coefficients=coefficients)
            
            name = paste0(c(c,ys,ysm), collapse='.')
            shifts[[paste0(c(c,ys,ysm), collapse='.')]] = data.table(name, pred_shift = est)
                
            }
        }   
        
     }
       shifts = rbindlist(shifts)
       shifts[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)]

       return(shifts)
}
                                          
# estimate lags

estimate_lags = function(samples, 
                         country=NULL, 
                         model_year=NULL,
                         le_values=NULL,
                         year_values=NULL,
                         year=NULL,
                         gdp_value = NULL, 
                         max_le_value = NULL,
                         coefficients = c('Intercept', 'gdp_log')) {
    
    colnames = names(samples)
    betas = paste0('b_', coefficients)
    random = str_subset(colnames, paste0('^r_.+\\[', country, '.', model_year, ','))

    s = samples[, c(betas, random)]

    pred = (s[,1] + s[,3]) + (s[,2] + s[,4]) * gdp_value
    pred = unlist(sapply(pred, function(x) 
            get_orig_values_weibull(x, max_value=max_le_value)))
     
    ind = NULL
    for (i in 1:length(pred)) {
        ind[i] = which.min(abs(le_values - pred[i]))
        }
           
    return(year_values[ind]-year)
}
                         
compute_lags = function(model,
                        data, 
                        country_labels, 
                        iyears,
                        max_le, 
                        predicted_le=FALSE, 
                        coefficients) {

    model_years = list('1950' = '1950-1969', '1970' = '1970-1989', '1990' = '1990', '2010' = '1990')
    samples = posterior_samples(model)
    
    lags = list()

    for (c in country_labels) {
    
    max_le_value = max_le[ctry==c, max_le]
    iyears = as.numeric(data[ctry==c & year %in% iyears, year])
    segments = as.character(unique(data[ctry==c, gyear]))
        
      if (predicted_le) {
          le_values = predict(model, data[ctry==c], summary=FALSE)
          le_values = apply(le_values, 2, mean)
          
          le_values = unlist(lapply(le_values,  function(x) 
              get_orig_values_weibull(x, max_le_value)))
          year_values = data[ctry==c, year]
        } 
        else { 
            le_values = data[ctry==c, le]
            year_values = data[ctry==c, year]       
        }
    
    for (ys in iyears) {
        
            fsegments = segments[-which(segments == model_years[[as.character(ys)]])]
            gdp_value= data[ctry==c & year==ys, gdp_log]
        
        for (ysm in fsegments) { 
         
            est = estimate_lags(samples=samples, 
                         country=c, 
                         model_year=ysm,
                         year=ys,
                         le_values=le_values,
                         year_values=year_values,
                         gdp_value=gdp_value, 
                         max_le_value=max_le_value,
                     )

                name = paste0(c(c,ys,ysm), collapse='.')
                lags[[paste0(c(c,ys,ysm), collapse='.')]] = data.table(name, pred_lag = est)

                  }
            }
    }
       lags = rbindlist(lags)
       lags[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)]

       return(lags)
}