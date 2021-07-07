data{
    vector[1261] zius_aid_pc;
    vector[1261] zilit;
    vector[1261] zpop;
    vector[1261] zinfrastructure;
    vector[1261] zyear;
    vector[1261] wy;
    vector[1261] zigdp_pc;
    int ctryearg[1261];
}
parameters{
    vector[74] a_cy;
    real a;
    real<lower=0> sigma;
    real<lower=0> sigma_cy;
    real b_gdp;
}
model{
    vector[1261] mu;
    vector[1261] pred;
    for ( i in 1:1261 ) {
        pred[i] = a_cy[ctryearg[i]] + b_gdp * zigdp_pc[i];
    }
    b_gdp ~ normal( 0 , 0.5 );
    sigma_cy ~ exponential( 1 );
    sigma ~ exponential( 1 );
    a ~ normal( 0 , 1 );
    a_cy ~ normal( a , sigma_cy );
    for ( i in 1:1261 ) {
        mu[i] = a_cy[ctryearg[i]] + b_gdp * zigdp_pc[i];
    }
    wy ~ normal( mu , sigma );
}
generated quantities{
    vector[1261] pred;
    for ( i in 1:1261 ) {
        pred[i] = a_cy[ctryearg[i]] + b_gdp * zigdp_pc[i];
    }
}