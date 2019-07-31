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
  vector[n_possible_epidemics] cumulative_prior_incidence_std;

  // hyperparameters set at runtime
  real<lower=0> sd_sd_effect_sizes;
  real<lower=0> sd_mean_intercept;
  real<lower=0> sd_sd_intercept;
  real<lower=1> nu;
  real sd_sd_errors;

}
  

parameters{

  // non centered hierarchical effects
  vector[n_subtypes] intercept_errors;
  //vector[n_subtypes] effect_antigenic_change_errors;
  vector[n_subtypes] effect_cumulative_prior_inc_errors;
  vector[n_possible_epidemics] errors;

  real mean_intercept;
  real<lower=0> sd_intercept;
  real<lower=0> sd_error;
  
  //real mean_effect_antigenic_change;
  // real<lower=0> sd_effect_antigenic_change;

  real mean_effect_cumulative_prior_inc;
  real<lower=0> sd_effect_cumulative_prior_inc;

}

transformed parameters{

  vector[n_possible_epidemics] logit_epi_prob;
  vector<lower=0, upper=1>[n_possible_epidemics] epi_prob;

  vector[n_subtypes] intercept;
  //vector[n_subtypes] effect_antigenic_change;
  vector[n_subtypes] effect_cumulative_prior_inc;

  intercept = mean_intercept + intercept_errors * sd_intercept;

  //effect_antigenic_change = mean_effect_antigenic_change +
  // effect_antigenic_change_errors * sd_effect_antigenic_change;

  effect_cumulative_prior_inc = mean_effect_cumulative_prior_inc +
    effect_cumulative_prior_inc_errors * sd_effect_cumulative_prior_inc;

  // calculate expected values of normed incidence
  for(epi_id in 1:n_possible_epidemics){
    logit_epi_prob[epi_id] =
      intercept[subtype[epi_id]] +
      
      //effect_antigenic_change[subtype[epi_id]] *
      //antigenic_change[epi_id] +
      
      effect_cumulative_prior_inc[subtype[epi_id]] *
      //(1 - antigenic_change[epi_id]) *
      cumulative_prior_incidence_std[epi_id];
  }

  epi_prob = inv_logit(logit_epi_prob + errors * sd_error);
}

model {
  
  epi_occurred ~ bernoulli(epi_prob);

  // non-centered intercept errors
  intercept_errors ~ normal(0, 1);

  
  // regression coeffecients hierarchical by subtype
  //effect_antigenic_change_errors ~ normal(0, 1);

  effect_cumulative_prior_inc_errors ~ normal(0, 1);


  mean_intercept ~ normal(0, sd_mean_intercept);
  sd_intercept ~ normal(0, sd_sd_intercept);

  // weakly informative priors on mean regression effects
  
  //mean_effect_antigenic_change ~ student_t(nu, 0, 2.5);
  //sd_effect_antigenic_change ~ normal(0, sd_sd_effect_sizes);

  mean_effect_cumulative_prior_inc ~ student_t(nu, 0, 2.5);
  sd_effect_cumulative_prior_inc ~ normal(0, sd_sd_effect_sizes);


  errors ~ normal(0, 1);

  sd_error ~ normal(0, sd_sd_errors);
}

