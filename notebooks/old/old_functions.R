compute_shifts = function(model, 
                          data, 
                          country_labels,
                          iyears, 
                          max_le, 
                          predicted_le=FALSE) {

    model_years = list('1950' = '1950', '1970' = '1950-1969', '1990' = '1970-1989', '2010'= '1990')
    
    shifts = list()

    for (c in country_labels) {
        
#         print(c)
        max_le_value = max_le[ctry==c, max_le]
        years = as.numeric(data[ctry==c & year %in% iyears, year])
        segments = as.character(unique(data[ctry==c, gyear]))

        for (ys in years) {
            
            
#             print(ys)
            tv = model_years[[as.character(ys)]]
#             print(tv)
            if (!(tv %in% segments)) { next }
            
            est = estimate_shift(model=model, 
                                 data=data,
                                 country=c, 
                                 model_year=tv, 
                                 max_le_value=max_le_value, 
                                 predicted_le=predicted_le)
            
            name = paste0(c(c,ys,tv), collapse='.')
            shifts[[paste0(c(c,ys,tv), collapse='.')]] = data.table(name, pred_shift = est)

            }

       }

       shifts = rbindlist(shifts)
       shifts[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)]

       levels = c('1990', '1970-1989', '1950-1969', '1950')
       labels = c('>=1990', '1970-1989', '1950-1969', '<1950')
       shifts[, segment := factor(segment, levels=levels, 
                             labels=labels ) ] 

        return(shifts)
}



compute_shifts = function(model, 
                          data, 
                          country_labels,
                          iyears, 
                          max_le, 
                          predicted_le=FALSE) {

    model_years = list('1950' = '1950', '1970' = '1950-1969', '1990' = '1970-1989', '2010'= '1990')
    
    shifts = list()

    for (c in country_labels) {
        
#         print(c)
        max_le_value = max_le[ctry==c, max_le]
        years = as.numeric(data[ctry==c & year %in% iyears, year])
        segments = as.character(unique(data[ctry==c, gyear]))

        for (ys in years) {
            
            
#             print(ys)
            tv = model_years[[as.character(ys)]]
#             print(tv)
            if (!(tv %in% segments)) { next }
            
            est = estimate_shift(model=model, 
                                 data=data,
                                 country=c, 
                                 model_year=tv, 
                                 max_le_value=max_le_value, 
                                 predicted_le=predicted_le)
            
            name = paste0(c(c,ys,tv), collapse='.')
            shifts[[paste0(c(c,ys,tv), collapse='.')]] = data.table(name, pred_shift = est)

            }

       }

       shifts = rbindlist(shifts)
       shifts[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)]

       levels = c('1990', '1970-1989', '1950-1969', '1950')
       labels = c('>=1990', '1970-1989', '1950-1969', '<1950')
       shifts[, segment := factor(segment, levels=levels, 
                             labels=labels ) ] 

        return(shifts)
}

iyears = c(1950, 1970, 1990, 2010)
model_year = c('1950', '1950-1969', '1970-1989', '1990')

shift_total_obs = compute_shifts(samples=samples_total, 
                             data=total, 
                             model=t1,
                             country_labels=country_labels, 
                             iyears=iyears, 
                             model_year=model_year,
                             max_le=max_le_total, 
                             predicted_le=FALSE)

shift_total_pred = compute_shifts(samples=samples_total, 
                             data=total, 
                             model=t1,
                             country_labels=country_labels, 
                             iyears=iyears, 
                             model_year=model_year,
                             max_le=max_le_total, 
                             predicted_le=TRUE)

shift_total_obs[, sex := 'total']
shift_total_pred[, sex := 'total']


# estimate_shift = function(samples,
#                           le_value=NULL,
#                           country=NULL, 
#                           model_year=NULL,
#                           year=NULL,
#                           max_le = NULL,
#                           coefficients = c('Intercept', 'gdp_log')) {
    
#     colnames = names(samples)
                          
#     betas = paste0('b_', coefficients)
    
#     random = str_subset(colnames, paste0('^r_.+\\[', country, '.', model_year, ','))
#     s = samples[, c(betas, random)]
#     gdp_log = data[ctry==country & year==year, gdp_log]   
                                 
#     pred = (s[,1] + s[,3]) + (s[,2] + s[,4]) * gdp_log)
#     pred = unlist(sapply(pred, function(x) 
#             get_orig_values_weibull(x, max_value=max_le)))
    
#      return(le_value-pred)
            
# }

# compute_shifts = function(samples, model, data, 
#                           country_labels,
#                           predicted_le=FALSE,
#                           iyears,
#                           model_year, max_le) {

#     shifts = list()

#     for (c in country_labels) {
        
#         max_le_value = max_le[ctry==c, max_le]
#         years = as.numeric(data[ctry == c & year %in% iyears, year])
#         model_years = as.character(unique(data[ctry == c, gyear]))

#         for (ys in years) {

#             gdp_value = data[ctry==c & year==ys, gdp_log]
            
#             if (predicted_le) {
                
#                 le_value = predict(model, data[ctry==c & year==ys], summary=FALSE)
#                 le_value = unlist(sapply(le_value, function(x) 
#                                 get_orig_values_weibull(x,max_le_value)))    
#             } else {
#                 le_value =  data[ctry==c & year==ys, le]
#                 }
            
#             for (ysm in model_years) { 
    
#                 est = estimate_shift(samples, 
#                                      country = c, 
#                                      model_year = ysm, 
#                                      gdp_value = gdp_value,
#                                      le_value = le_value,
#                                      max_le_value=max_le_value)
#                 name = paste0(c(c,ys,ysm), collapse='.')
#                 shifts[[paste0(c(c,ys,ysm), collapse='.')]] = data.table(name, pred_shift = est)

#                   }

#             }

#         }

#        shifts = rbindlist(shifts)
#        shifts[, c('ctry', 'year', 'segment') := tstrsplit(name, ".", fixed=TRUE)]

#        levels = c('1990', '1970-1989', '1950-1969', '1950')
#        labels = c('>=1990', '1970-1989', '1950-1969', '<1950')
#        shifts[, segment := factor(segment, levels=levels, 
#                              labels=labels ) ] 

#         return(shifts)
# }


# estimate_shift = function(samples, 
#                           country=NULL, 
#                           model_year=NULL,
#                           le_value=NULL, # value or vector
#                           gdp_value = NULL, 
#                           max_le_value = NULL,
#                           coefficients = c('Intercept', 'gdp_log')) {
    
#     colnames = names(samples)
#     betas = paste0('b_', coefficients)
    
#     random = str_subset(colnames, paste0('^r_.+\\[', country, '.', model_year, ','))

#     s = samples[, c(betas, random)]

#     pp = (s[,1] + s[,3]) + (s[,2] + s[,4]) * gdp_value
#     pg = unlist(sapply(pp, function(x) 
#             get_orig_values_weibull(x, max_value=max_le_value)))
    
#      return(le_value-pg)
# }

estimate_shift = function(model=NULL, 
                          data=NULL,
                          country=NULL, 
                          model_year=NULL,
                          year=NULL,
                          max_le_value=NULL, 
                          predicted_le=FALSE) {
    
    if (predicted_le) {
        le_value = predict(model, data[ctry==country & year==year], summary=FALSE)
        le_value = unlist(sapply(le_value, function(x)
        get_orig_values_weibull(x,max_le_value)))    
    } 
    else {
        le_value =  data[ctry==country & year==year, le]
    }
                                         
    newdata = data[ctry==country & year==year]
    newdata[, ctry_year := paste0(c(country, model_year), collapse='.')]

    pred = predict(model, newdata, summary=FALSE, allow_new_levels=TRUE) # predict values
    pred = unlist(sapply(pred, function(x) 
            get_orig_values_weibull(x, max_value=max_le_value)))
    
     return(le_value-pred)
}
                         
                         