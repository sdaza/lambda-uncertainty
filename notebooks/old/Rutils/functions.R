# plot for frequentist approach
plot_predictions = function (x, y, pred) {

    #rmse = sqrt(mean(residuals(model)^2))
    ggplot(data.frame(x, y, pred), aes(x = x, y = y)) +
    geom_point(size=0.7, color='#e34a33', alpha=0.4) + 
    geom_line(aes(x = x, y = pred), color='#2b8cbe', size = 0.4) +
   # labs(title = paste0("RMSE: ", round(rmse, 2))) + 
    theme_classic()

}


# plot function for posterior
plot_posterior = function(model, x, y, pred_name='y_pred', ylabel='y ', xlabel='x',
                          prob=0.90, type='stan') {
    
    if (type == 'stan' ) {
        ypred = extract(model, pred_name)[[1]]
    }
    if (type == 'matrix') {
        ypred = model 
    }
    
    m =  as.vector(apply(ypred, 2, median, na.rm=TRUE))
    lo = as.vector(apply(ypred, 2, function(x) quantile(x, prob=(1-prob), na.rm = TRUE)))
    hi = as.vector(apply(ypred, 2, function(x) quantile(x, prob=prob, na.rm = TRUE)))
        
    dat = data.frame(y, x, m, lo, hi)
    
    ggplot(dat, aes(x=x, y=y)) + 
        geom_line(aes(y=m), color='#2b8cbe', size = 0.4)  +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = '#a6bddb', alpha=0.2)  + 
        geom_point(size=0.7, color='#e34a33', alpha=0.4) + 
        labs(x=xlabel, y=ylabel) + 
        theme_classic()

}
                         
                         
# function for box-cox transformation ofr prediction
pred_bc = function(x, lambda) {
    if (lambda!=0) {
        return(sapply(x, function(x) (x * lambda + 1)^(1/lambda)))
        }
    else {
        return(exp(x))
    }
}  
        
plot_posterior_w = function(stan_model_list, w, x, y, pred_name='y_pred', ylabel='y ', xlabel='x',
                          prob=0.90) {
    
    nmodels = length(stan_model_list)
    
    preds = list()
    for (i in 1:nmodels) {
        preds[[i]] =  extract(stan_model_list[[i]], pred_name)[[1]]    
    }
    
    ypred = matrix(0, nrow(preds[[1]]), ncol(preds[[1]]))

    for (i in 1:nmodels) { 
        temp =  w[i] * preds[[i]]
        ypred = ypred + temp
        
    }
    
    
    m =  as.vector(apply(ypred, 2, median))
    lo = as.vector(apply(ypred, 2, function(x) quantile(x, prob=(1-prob))))
    hi = as.vector(apply(ypred, 2, function(x) quantile(x, prob=prob)))
        
    dat = data.frame(y, x, m, lo, hi)
    
    ggplot(dat, aes(x=x, y=y)) + 
        geom_line(aes(y=m), color='#2b8cbe', size = 0.4)  +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = '#a6bddb', alpha=0.2)  + 
        geom_point(size=0.7, color='#e34a33', alpha=0.4) + 
        labs(x=xlabel, y=ylabel) + 
        theme_classic()

}
  