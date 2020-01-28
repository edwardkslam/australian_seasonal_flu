functions {
  vector gelman_standardize(vector x){
    return((x - mean(x)) / (2 * sd(x)));
  }
}

data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower=0> n_cities;
  int<lower=0> n_subtypes;
  int<lower=0> n_epidemics;
  vector[n_epidemics] normed_metric;
  int<lower=0, upper=n_cities> city[n_epidemics];
  int<lower=0, upper=n_subtypes> subtype[n_epidemics];

  /////////////////////////////////
  // predictors
  /////////////////////////////////

  // antigenic change predictors
  vector<lower=0, upper=1>[n_epidemics] antigenic_change;
  vector[n_epidemics] cumulative_prior_incidence_std;

  // climate predictors
  vector<lower=0>[n_epidemics] abs_humidity;
  vector<lower=0>[n_epidemics] rainfall;

  // competition among epidemics predictors
  vector<lower=0, upper=1>[n_epidemics] is_first_of_season;
  vector[n_epidemics] prior_season_activity_std;

  // general seasonality predictors
  vector<lower=0>[n_epidemics] start_date;

  /////////////////////////////////
  // hyperparameters set at runtime
  /////////////////////////////////
  real<lower=0> sd_sd_incidences;
  real<lower=0> sd_mean_intercept;
  real<lower=0> sd_sd_intercept;
  real<lower=0> nu;
}

transformed data {

  // center and scale predictors 
  vector[n_epidemics] abs_humidity_std;
  vector[n_epidemics] rainfall_std;
  vector[n_epidemics] start_date_std;
  
  abs_humidity_std = gelman_standardize(abs_humidity);
  
  rainfall_std = gelman_standardize(rainfall);
  
  start_date_std = gelman_standardize(start_date);
  
}

  

parameters{
  real<lower=0> sd_incidences;


  // non centered hierarchical effects
  vector[n_subtypes] intercept_errors;

  vector[n_subtypes] effect_antigenic_change_errors;
  vector[n_subtypes] effect_cumulative_prior_inc_errors;

  vector[n_subtypes] effect_abs_humidity_errors;

  vector[n_subtypes] effect_is_first_of_season_errors;
  vector[n_subtypes] effect_prior_season_activity_errors;

  vector[n_subtypes] effect_start_date_errors;

  vector[n_subtypes] effect_rainfall_errors;


  real mean_intercept;
  real<lower=0> sd_intercept;

  
  real mean_effect_antigenic_change;
  real<lower=0> sd_effect_antigenic_change;
  real mean_effect_cumulative_prior_inc;
  real<lower=0> sd_effect_cumulative_prior_inc;

  real mean_effect_abs_humidity;
  real<lower=0> sd_effect_abs_humidity;

  real mean_effect_is_first_of_season;
  real<lower=0> sd_effect_is_first_of_season;
  real mean_effect_prior_season_activity;
  real<lower=0> sd_effect_prior_season_activity;

  real mean_effect_start_date;
  real<lower=0> sd_effect_start_date;

  real mean_effect_rainfall;
  real<lower=0> sd_effect_rainfall;

  real<lower=0> sd_sd_effect_sizes;
}

transformed parameters{
  vector[n_epidemics] expected_incidences;

  vector[n_subtypes] intercept;

  vector[n_subtypes] effect_antigenic_change;
  vector[n_subtypes] effect_cumulative_prior_inc;
  
  vector[n_subtypes] effect_abs_humidity;

  vector[n_subtypes] effect_is_first_of_season;
  vector[n_subtypes] effect_prior_season_activity;

  vector[n_subtypes] effect_start_date;

  vector[n_subtypes] effect_rainfall;

  intercept = mean_intercept + intercept_errors * sd_intercept;

  
  effect_antigenic_change = mean_effect_antigenic_change +
    effect_antigenic_change_errors * sd_effect_antigenic_change;

  effect_cumulative_prior_inc = mean_effect_cumulative_prior_inc +
    effect_cumulative_prior_inc_errors * sd_effect_cumulative_prior_inc;


  effect_abs_humidity = mean_effect_abs_humidity +
    effect_abs_humidity_errors * sd_effect_abs_humidity;  
  
  effect_is_first_of_season = mean_effect_is_first_of_season +
    effect_is_first_of_season_errors * sd_effect_is_first_of_season;
  
  effect_prior_season_activity = mean_effect_prior_season_activity +
    effect_prior_season_activity_errors * sd_effect_prior_season_activity;

  effect_start_date = mean_effect_start_date +
    effect_start_date_errors * sd_effect_start_date;

  effect_rainfall = mean_effect_rainfall +
    effect_rainfall_errors * sd_effect_rainfall;
  
  // calculate expected values of normed incidence
  for(epi_id in 1:n_epidemics){
    expected_incidences[epi_id] =
      intercept[subtype[epi_id]] +
      
      effect_antigenic_change[subtype[epi_id]] *
      antigenic_change[epi_id] +

      effect_cumulative_prior_inc[subtype[epi_id]] *
      (1 - antigenic_change[epi_id]) *
      cumulative_prior_incidence_std[epi_id] +

      effect_abs_humidity[subtype[epi_id]] *
      abs_humidity_std[epi_id] +
      
      effect_is_first_of_season[subtype[epi_id]] *
      is_first_of_season[epi_id] +

      effect_prior_season_activity[subtype[epi_id]] *
      (1 - is_first_of_season[epi_id]) * 
      prior_season_activity_std[epi_id] +

      effect_start_date[subtype[epi_id]] *
      start_date_std[epi_id] +

      effect_rainfall[subtype[epi_id]] *
      rainfall_std[epi_id];

  }
}

model {
  
  normed_metric ~ normal(expected_incidences, sd_incidences);
  
  // estimated intercept with non-centered intercept errors
  mean_intercept ~ normal(0, sd_mean_intercept);
  sd_intercept ~ normal(0, sd_sd_intercept);
  intercept_errors ~ normal(0, 1);

  // regression coeffecients hierarchical by subtype
  effect_antigenic_change_errors ~ normal(0, 1);
  effect_cumulative_prior_inc_errors ~ normal(0, 1);

  effect_abs_humidity_errors ~ normal(0, 1);

  effect_is_first_of_season_errors ~ normal(0, 1);
  effect_prior_season_activity_errors ~ normal(0, 1);

  effect_start_date_errors ~ normal(0, 1);
  
  effect_rainfall_errors ~ normal(0, 1);

  // Gaussian global errors
  sd_incidences ~ normal(0, sd_sd_incidences);

  // weakly informative priors on mean regression effects
  // pooling limited by an estimated sd with a gaussian prior
  
  mean_effect_antigenic_change ~ student_t(nu, 0, 2.5);
  sd_effect_antigenic_change ~ normal(0, sd_sd_effect_sizes);

  mean_effect_abs_humidity ~ student_t(nu, 0, 2.5);
  sd_effect_abs_humidity ~ normal(0, sd_sd_effect_sizes);
    
  mean_effect_cumulative_prior_inc ~ student_t(nu, 0, 2.5);
  sd_effect_cumulative_prior_inc ~ normal(0, sd_sd_effect_sizes);

  mean_effect_is_first_of_season ~ student_t(nu, 0, 2.5);
  sd_effect_is_first_of_season ~ normal(0, sd_sd_effect_sizes);
  
  mean_effect_prior_season_activity ~ student_t(nu, 0, 2.5);
  sd_effect_prior_season_activity ~ normal(0, sd_sd_effect_sizes);

  mean_effect_start_date ~ student_t(nu, 0, 2.5);
  sd_effect_start_date ~ normal(0, sd_sd_effect_sizes);

  mean_effect_rainfall ~ student_t(nu, 0, 2.5);
  sd_effect_rainfall ~ normal(0, sd_sd_effect_sizes);
  
  sd_sd_effect_sizes ~ normal(0, 1);
  
}

generated quantities {
  vector[n_epidemics] possible_sizes_given_params;

  for(epi_id in 1:n_epidemics){
    possible_sizes_given_params[epi_id] = normal_rng(expected_incidences[epi_id], sd_incidences);
  }
}
