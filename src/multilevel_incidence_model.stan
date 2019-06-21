functions {
  vector gelman_standardize(vector x){
    return((x - mean(x)) / (2 * sd(x)));
  }
}

data {

  int<lower=0> n_cities;
  int<lower=0> n_epidemics;
  vector<lower=0>[n_epidemics] incidences;
  int<lower=0, upper=n_cities> city[n_epidemics];

  // predictors
  vector<lower=0, upper=1>[n_epidemics] antigenic_change;
  vector<lower=0>[n_epidemics] abs_humidity;
  vector<lower=0>[n_epidemics] cumulative_prior_incidence;
  vector<lower=0>[n_epidemics] prior_activity;
}

transformed data {

  // center and scale predictors 
  vector[n_epidemics] abs_humidity_std;
  vector[n_epidemics] cumulative_prior_incidence_std;
  vector[n_epidemics] prior_activity_std;

  abs_humidity_std = gelman_standardize(abs_humidity);

  cumulative_prior_incidence_std =
    gelman_standardize(cumulative_prior_incidence);

  prior_activity_std =
    gelman_standardize(prior_activity);

}

  

parameters{
  real<lower=0> sd_incidences;

  real effect_antigenic_change;
  real effect_abs_hum;
  real effect_cumulative_prior_inc;
  real effect_prior_activity;

  real intercept;
  
}

transformed parameters{

  vector[n_epidemics] expected_incidences;

  // elementwise
  expected_incidences = intercept +
    effect_antigenic_change * antigenic_change +
    effect_abs_hum * abs_humidity_std +
    effect_cumulative_prior_inc * cumulative_prior_incidence_std +
    effect_prior_activity * prior_activity_std;
}

model {

  incidences ~ normal(expected_incidences, sd_incidences);

  effect_abs_hum ~ normal(0, 1);
  effect_cumulative_prior_inc ~ normal(0, 1);
  effect_prior_activity ~ normal(0, 1);

  intercept ~ normal(0, 3);
  
  sd_incidences ~ normal(0, 10);
}
