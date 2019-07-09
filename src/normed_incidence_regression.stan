functions {
  vector gelman_standardize(vector x){
    return((x - mean(x)) / (2 * sd(x)));
  }
}

data {

  int<lower=0> n_cities;
  int<lower=0> n_epidemics;
  vector[n_epidemics] normed_metric;
  int<lower=0, upper=n_cities> city[n_epidemics];

  // predictors
  vector<lower=0, upper=1>[n_epidemics] antigenic_change;
  vector<lower=0>[n_epidemics] abs_humidity;
  vector<lower=0>[n_epidemics] temperature;
  vector<lower=0>[n_epidemics] cumulative_prior_incidence;
  vector<lower=0>[n_epidemics] other_subtype_activity;
  vector<lower=0>[n_epidemics] start_date;

  // hyperparameters set at runtime
  real<lower=0> sd_sd_incidences;
}

transformed data {

  // center and scale predictors 
  vector[n_epidemics] abs_humidity_std;
  vector[n_epidemics] temperature_std;
  vector[n_epidemics] other_subtype_activity_std;
  vector[n_epidemics] cumulative_prior_incidence_std;
  
  abs_humidity_std = gelman_standardize(abs_humidity);

  temperature_std = gelman_standardize(temperature);
  
  other_subtype_activity_std =
    gelman_standardize(other_subtype_activity);

  cumulative_prior_incidence_std =
    gelman_standardize(cumulative_prior_incidence);
}

  

parameters{
  real<lower=0> sd_incidences;

  real intercept;
  real effect_antigenic_change;
  real effect_abs_humidity;
  real effect_temperature;
  real effect_cumulative_prior_inc;
  real effect_other_subtype_activity;
  real effect_start_date;
  
}

transformed parameters{

  vector[n_epidemics] expected_incidences;
  real log_intercept;

  log_intercept = log(intercept);  

  // sum log(reporting rate) and other stuff because log scale!
  expected_incidences =
    log_intercept +
    effect_antigenic_change * antigenic_change +
    effect_abs_humidity * abs_humidity_std +
    effect_start_date * start_date +
    effect_temperature * temperature_std +
    effect_cumulative_prior_inc * (1 - antigenic_change) .*
    cumulative_prior_incidence_std +
    effect_other_subtype_activity * other_subtype_activity_std;
}

model {
  
  normed_metric ~ normal(expected_incidences, sd_incidences);

  intercept ~ normal(0, 1);
  
  // weakly informative priors on regression coeffecients
  effect_antigenic_change ~ normal(0, 1);
  effect_abs_humidity ~ normal(0, 1);
  effect_cumulative_prior_inc ~ normal(0, 1);
  effect_other_subtype_activity ~ normal(0, 1);
  effect_start_date ~ normal(0, 1);
  effect_temperature ~ normal(0, 1);

  sd_incidences ~ normal(0, sd_sd_incidences);
}

generated quantities {
  vector[n_epidemics] possible_sizes_given_params;

  for(epi_id in 1:n_epidemics){
    possible_sizes_given_params[epi_id] = normal_rng(expected_incidences[epi_id], sd_incidences);
  }
}
