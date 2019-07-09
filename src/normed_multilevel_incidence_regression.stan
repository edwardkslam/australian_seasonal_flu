functions {
  vector gelman_standardize(vector x){
    return((x - mean(x)) / (2 * sd(x)));
  }
}

data {

  int<lower=0> n_cities;
  int<lower=0> n_subtypes;
  int<lower=0> n_epidemics;
  vector[n_epidemics] normed_metric;
  int<lower=0, upper=n_cities> city[n_epidemics];
  int<lower=0, upper=n_subtypes> subtype[n_epidemics];

  // predictors
  vector<lower=0, upper=1>[n_epidemics] antigenic_change;
  vector<lower=0>[n_epidemics] abs_humidity;
  vector<lower=0>[n_epidemics] temperature;
  vector<lower=0>[n_epidemics] cumulative_prior_incidence;
  vector<lower=0>[n_epidemics] other_subtype_activity;
  vector<lower=0>[n_epidemics] start_date;

  // hyperparameters set at runtime
  real<lower=0> sd_sd_incidences;
  real<lower=0> sd_mean_effect_sizes;
  real<lower=0> sd_sd_effect_sizes;
  real<lower=0> sd_mean_intercept;
  real<lower=0> sd_sd_intercept;
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

  vector[n_subtypes] intercepts;
  vector[n_subtypes] effect_antigenic_change;
  vector[n_subtypes] effect_abs_humidity;
  vector[n_subtypes] effect_temperature;
  vector[n_subtypes] effect_cumulative_prior_inc;
  vector[n_subtypes] effect_other_subtype_activity;
  vector[n_subtypes] effect_start_date;


  real mean_intercept;
  real<lower=0> sd_intercept;
  

  real mean_effect_antigenic_change;
  real<lower=0> sd_effect_antigenic_change;

  real mean_effect_abs_humidity;
  real<lower=0> sd_effect_abs_humidity;
  
  real mean_effect_temperature;
  real<lower=0> sd_effect_temperature;

  real mean_effect_cumulative_prior_inc;
  real<lower=0> sd_effect_cumulative_prior_inc;

  real mean_effect_other_subtype_activity;
  real<lower=0> sd_effect_other_subtype_activity;

  real mean_effect_start_date;
  real<lower=0> sd_effect_start_date;
}

transformed parameters{

  vector[n_epidemics] expected_incidences;


  for(epi_id in 1:n_epidemics){
    expected_incidences[epi_id] =
      intercepts[subtype[epi_id]] +
      effect_antigenic_change[subtype[epi_id]] *
      antigenic_change[epi_id] +
      effect_abs_humidity[subtype[epi_id]] *
      abs_humidity_std[epi_id] +
      effect_start_date[subtype[epi_id]] * start_date[epi_id] +
      effect_temperature[subtype[epi_id]] * temperature_std[epi_id] +
      effect_cumulative_prior_inc[subtype[epi_id]] *
      (1 - antigenic_change[epi_id]) *
      cumulative_prior_incidence_std[epi_id] +
      effect_other_subtype_activity[subtype[epi_id]] *
      other_subtype_activity_std[epi_id];
  }
}

model {
  
  normed_metric ~ normal(expected_incidences, sd_incidences);

  intercepts ~ normal(mean_intercept, sd_intercept);
  
  // regression coeffecients hierarchical by subtype
  effect_antigenic_change ~ normal(mean_effect_antigenic_change,
                                   sd_effect_antigenic_change);

  effect_abs_humidity ~ normal(mean_effect_abs_humidity,
                               sd_effect_abs_humidity);

  effect_cumulative_prior_inc ~ normal(mean_effect_cumulative_prior_inc,
                                       sd_effect_cumulative_prior_inc);

  effect_other_subtype_activity ~ normal(mean_effect_other_subtype_activity,
                                         sd_effect_other_subtype_activity);

  effect_start_date ~ normal(mean_effect_start_date,
                             sd_effect_start_date);

  effect_temperature ~ normal(mean_effect_temperature,
                              sd_effect_temperature);

  sd_incidences ~ normal(0, sd_sd_incidences);


  mean_intercept ~ normal(0, sd_mean_intercept);
  sd_intercept ~ normal(0, sd_sd_intercept);

  // weakly informative priors on mean regression effects
  
  mean_effect_antigenic_change ~ normal(0, sd_mean_effect_sizes);
  sd_effect_antigenic_change ~ normal(0, sd_sd_effect_sizes);

  mean_effect_abs_humidity ~ normal(0, sd_mean_effect_sizes);
  sd_effect_abs_humidity ~ normal(0, sd_sd_effect_sizes);
  
  mean_effect_temperature ~ normal(0, sd_mean_effect_sizes);
  sd_effect_temperature ~ normal(0, sd_sd_effect_sizes);

  mean_effect_cumulative_prior_inc ~ normal(0, sd_mean_effect_sizes);
  sd_effect_cumulative_prior_inc ~ normal(0, sd_sd_effect_sizes);

  mean_effect_other_subtype_activity ~ normal(0, sd_mean_effect_sizes);
  sd_effect_other_subtype_activity ~ normal(0, sd_sd_effect_sizes);

  mean_effect_start_date ~ normal(0, sd_mean_effect_sizes);
  sd_effect_start_date ~ normal(0, sd_sd_effect_sizes);

  
}

generated quantities {
  vector[n_epidemics] possible_sizes_given_params;

  for(epi_id in 1:n_epidemics){
    possible_sizes_given_params[epi_id] = normal_rng(expected_incidences[epi_id], sd_incidences);
  }
}
