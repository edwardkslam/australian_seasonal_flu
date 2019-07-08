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

  // hyperparameters set at runtime
  real<lower=0> sd_city_reporting_rates_per_hundred;
  real<lower=0,upper=100> mean_city_reporting_rates_per_hundred;
  real<lower=0> alpha_average_epi_attack_rate;
  real<lower=0> beta_average_epi_attack_rate;
  real<lower=0> sd_sd_incidences;
}

transformed data {

  // center and scale predictors 
  vector[n_epidemics] abs_humidity_std;
  vector[n_epidemics] cumulative_prior_incidence_std;
  vector[n_epidemics] prior_activity_std;

  abs_humidity_std = gelman_standardize(abs_humidity);

  cumulative_prior_incidence_std = cumulative_prior_incidence;

  prior_activity_std =
    gelman_standardize(prior_activity);

}

  

parameters{
  real<lower=0> sd_incidences;

  real effect_antigenic_change;
  real effect_abs_hum;
  real effect_cumulative_prior_inc;
  real effect_prior_activity;

  real<lower=0, upper=1> average_epi_attack_rate;

  vector<lower=0, upper=100>[n_cities] city_reporting_rates_per_hundred;
  
}

transformed parameters{

  vector[n_epidemics] epi_reporting_rates;
  vector[n_epidemics] expected_incidences;
  real log_intercept;

  log_intercept = log(average_epi_attack_rate * 1000000);

  // for now, reporting rate is fixed in time for a given city
  for(epi_id in 1:n_epidemics){
    epi_reporting_rates[epi_id] =
      city_reporting_rates_per_hundred[city[epi_id]] / 100;
  }
  


  // sum log(reporting rate) and other stuff because log scale!
  expected_incidences = log(epi_reporting_rates) +
    log_intercept +
    effect_antigenic_change * antigenic_change +
    effect_abs_hum * abs_humidity_std +
    effect_cumulative_prior_inc * cumulative_prior_incidence_std * (1 - antigenic_change) +
    effect_prior_activity * prior_activity_std;
}

model {

  city_reporting_rates_per_hundred ~ normal(mean_city_reporting_rates_per_hundred,
                                            sd_city_reporting_rates_per_hundred);
  
  incidences ~ normal(expected_incidences, sd_incidences);

  effect_antigenic_change ~ normal(0, 1);
  effect_abs_hum ~ normal(0, 1);
  effect_cumulative_prior_inc ~ normal(0, 1);
  effect_prior_activity ~ normal(0, 1);

  average_epi_attack_rate ~ beta(alpha_average_epi_attack_rate,
                                 beta_average_epi_attack_rate);
  
  sd_incidences ~ normal(0, sd_sd_incidences);
}

generated quantities {
  vector[n_epidemics] true_attack_rate;
  vector[n_epidemics] true_report;

  for(epi_id in 1:n_epidemics){
    true_report[epi_id] = normal_rng(expected_incidences[epi_id], sd_incidences);
  }

  true_attack_rate = exp(true_report) ./ (1e6 * epi_reporting_rates);
}
