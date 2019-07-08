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
  vector<lower=0>[n_epidemics] temperature;
  vector<lower=0>[n_epidemics] cumulative_prior_incidence;
  vector<lower=0>[n_epidemics] other_subtype_activity;
  vector<lower=0>[n_epidemics] start_date;

  // hyperparameters set at runtime
  real<lower=0> sd_city_reporting_rates_per_mil;
  real<lower=0,upper=1e6> mean_city_reporting_rates_per_mil;
  real<lower=0> alpha_average_epi_attack_rate;
  real<lower=0> beta_average_epi_attack_rate;
  real<lower=0> sd_sd_incidences;
}

transformed data {

  // center and scale predictors 
  vector[n_epidemics] abs_humidity_std;
  vector[n_epidemics] temperature_std;
  vector[n_epidemics] cumulative_prior_incidence_std;
  vector[n_epidemics] other_subtype_activity_std;

  abs_humidity_std = gelman_standardize(abs_humidity);

  temperature_std = gelman_standardize(temperature);
  
  cumulative_prior_incidence_std =
    gelman_standardize(cumulative_prior_incidence);

  other_subtype_activity_std =
    gelman_standardize(other_subtype_activity);

}

  

parameters{
  real<lower=0> sd_incidences;

  real effect_antigenic_change;
  real effect_abs_humidity;
  real effect_cumulative_prior_inc;
  real effect_other_subtype_activity;
  real effect_start_date;
  real effect_temperature;

  real<lower=0, upper=1> average_epi_attack_rate;

  vector<lower=0, upper=1e6>[n_cities] city_reporting_rates_per_mil;
  
}

transformed parameters{

  vector[n_epidemics] epi_reporting_rates;
  vector[n_epidemics] expected_incidences;
  real log_intercept;

  log_intercept = log(average_epi_attack_rate * 1e6);  

  // for now, reporting rate is fixed in time for a given city
  for(epi_id in 1:n_epidemics){
    epi_reporting_rates[epi_id] =
      city_reporting_rates_per_mil[city[epi_id]] / 1e6;
  }

  
  // sum log(reporting rate) and other stuff because log scale!
  expected_incidences =
    log(epi_reporting_rates) + 
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

  city_reporting_rates_per_mil ~ normal(mean_city_reporting_rates_per_mil,
                                        sd_city_reporting_rates_per_mil);
  
  incidences ~ normal(expected_incidences, sd_incidences);

  // weakly informative priors on regression coeffecients
  effect_antigenic_change ~ normal(0, 1);
  effect_abs_humidity ~ normal(0, 1);
  effect_cumulative_prior_inc ~ normal(0, 1);
  effect_other_subtype_activity ~ normal(0, 1);
  effect_start_date ~ normal(0, 1);
  effect_temperature ~ normal(0, 1);
  
  average_epi_attack_rate ~ beta(alpha_average_epi_attack_rate,
                                 beta_average_epi_attack_rate);
  
  sd_incidences ~ normal(0, sd_sd_incidences);
}

generated quantities {
  vector[n_epidemics] true_attack_rate;
  vector[n_epidemics] epi_with_change;
  vector[n_epidemics] epi_no_change;
  vector[n_epidemics] epi_with_change_attack;
  vector[n_epidemics] epi_no_change_attack;
  vector[n_epidemics] true_report;

  for(epi_id in 1:n_epidemics){
    true_report[epi_id] =
      normal_rng(expected_incidences[epi_id], sd_incidences);
    
    epi_with_change[epi_id] =
      normal_rng((expected_incidences[epi_id] -
                  effect_antigenic_change *
                  antigenic_change[epi_id] +
                  effect_antigenic_change),
                 sd_incidences);
    
    epi_no_change[epi_id] =
      normal_rng((expected_incidences[epi_id] -
                  effect_antigenic_change *
                  antigenic_change[epi_id]),
                 sd_incidences);
  }

  true_attack_rate = exp(true_report) ./ (1e6 * epi_reporting_rates);
  epi_with_change_attack = exp(epi_with_change) ./
    (1e6 * epi_reporting_rates);
  epi_no_change_attack = exp(epi_no_change) ./
    (1e6 * epi_reporting_rates);
}
