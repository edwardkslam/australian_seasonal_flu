functions {
  vector gelman_standardize(vector x){
    return((x - mean(x)) / (2 * sd(x)));
  }
}

data {

  int<lower=0> n_cities;
  int<lower=0> n_subtypes;
  int<lower=0> n_possible_epidemics;
  int<lower=0, upper=1> epi_occurred[n_possible_epidemics];
  int<lower=0, upper=n_cities> city[n_possible_epidemics];
  int<lower=0, upper=n_subtypes> subtype[n_possible_epidemics];

  // predictors
  vector<lower=0, upper=1>[n_possible_epidemics] antigenic_change;
  vector<lower=0>[n_possible_epidemics] cumulative_prior_incidence;

  // hyperparameters set at runtime
  real<lower=0> sd_mean_effect_sizes;
  real<lower=0> sd_sd_effect_sizes;
  real<lower=0> sd_mean_intercept;
  real<lower=0> sd_sd_intercept;
}

transformed data {

  // center and scale predictors 
  vector[n_possible_epidemics] cumulative_prior_incidence_std;
  
  cumulative_prior_incidence_std =
    gelman_standardize(cumulative_prior_incidence);
}

  

parameters{

  // non centered hierarchical effects
  vector[n_subtypes] intercept_errors;
  vector[n_subtypes] effect_antigenic_change_errors;
  vector[n_subtypes] effect_cumulative_prior_inc_errors;


  real mean_intercept;
  real<lower=0> sd_intercept;
  

  //real<lower=0> sd_effects;
  
  real mean_effect_antigenic_change;
  real<lower=0> sd_effect_antigenic_change;

  real mean_effect_cumulative_prior_inc;
  real<lower=0> sd_effect_cumulative_prior_inc;

}

transformed parameters{

  vector[n_possible_epidemics] logit_epi_prob;
  vector[n_possible_epidemics] epi_prob;

  for(epi_id in 1:n_possible_epidemics){
    logit_epi_prob[epi_id] =

      mean_intercept +
      intercept_errors[subtype[epi_id]] *
      sd_intercept +
      
      (mean_effect_antigenic_change +
       effect_antigenic_change_errors[subtype[epi_id]] *
       sd_effect_antigenic_change) *
      antigenic_change[epi_id] +

      (mean_effect_cumulative_prior_inc +
       effect_cumulative_prior_inc_errors[subtype[epi_id]] *
       sd_effect_cumulative_prior_inc) *
      (1 - antigenic_change[epi_id]) *
      cumulative_prior_incidence_std[epi_id];
  }

  epi_prob = inv_logit(logit_epi_prob);
}

model {
  
  epi_occurred ~ bernoulli(epi_prob);

  // non-centered intercept errors
  intercept_errors ~ normal(0, 1);

  
  // regression coeffecients hierarchical by subtype
  effect_antigenic_change_errors ~ normal(0, 1);

  effect_cumulative_prior_inc_errors ~ normal(0, 1);


  mean_intercept ~ normal(0, sd_mean_intercept);
  sd_intercept ~ normal(0, sd_sd_intercept);

  // weakly informative priors on mean regression effects
  
  mean_effect_antigenic_change ~ normal(0, sd_mean_effect_sizes);
  sd_effect_antigenic_change ~ normal(0, sd_sd_effect_sizes);

  mean_effect_cumulative_prior_inc ~ normal(0, sd_mean_effect_sizes);
  sd_effect_cumulative_prior_inc ~ normal(0, sd_sd_effect_sizes);
  
}
