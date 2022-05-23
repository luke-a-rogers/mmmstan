functions {
  #include functions.stan
}

data {
  // Model structure
  int<lower=0, upper=1> model_time;
  int<lower=0, upper=1> model_term;
  int<lower=0, upper=1> model_size;
  // Index limits
  int<lower=2> N; // Number of model steps (months/quarters/years) in the study
  int<lower=2> D; // Maximum duration at large in model steps
  int<lower=1> L; // Number of size or sex classes
  int<lower=2> X; // Number of geographic regions
  int<lower=1> T; // Number of years in the study
  int<lower=1> K; // Number of terms (quarters/trimesters/halfs) per year
  // Constants
  int<lower=1> H; // Number of steps per term
  int<lower=1> J; // Number of steps per year
  // Movement index arrays
  array[N] int<lower=1, upper=T> n_to_t; // Model step to year index
  array[N] int<lower=1, upper=K> n_to_k; // Model step to term index
  // Tag data
  array[N - 1, D, L, X, X] int<lower=0> tags;
  // Movement index (will be paired with a matrix version)
  array[X, X] int<lower=0, upper=1> movement_index;
  // Movement step mean priors
  vector<lower=0, upper=1>[X] mu_movement_diag;
  vector<lower=0>[X] sd_movement_diag;
  // Fishing rate priors
  array[T] vector<lower=0>[X] mu_fishing_rate;
  real<lower=0> cv_fishing_rate;
  // Selectivity priors
//  vector<lower=0, upper=1>[L - 1] mu_selectivity;
//  real<lower=0> cv_selectivity;
  // Fishing weight priors
//  array[K] vector<lower=0, upper=1>[X] mu_fishing_weight;
//  real<lower=0> cv_fishing_weight;
  // Natural mortality rate priors
  vector<lower=0>[X] mu_natural_mortality_rate;
  vector<lower=0>[X] sd_natural_mortality_rate;
  // Fractional (per tag) reporting rate priors
  vector<lower=0, upper=1>[X] mu_reporting_rate;
  vector<lower=0>[X] sd_reporting_rate;
  // Fractional (per tag) initial loss rate priors
  real<lower=0, upper=1> mu_initial_loss_rate;
  real<lower=0> sd_initial_loss_rate;
  // Instantaneous ongoing loss rate priors
  real<lower=0> mu_ongoing_loss_rate;
  real<lower=0> sd_ongoing_loss_rate;
  // Autoregression priors
//  real<lower=0, upper=1> mu_autoregress;
//  real<lower=0, upper=1> sd_autoregress;
//  real<lower=0> mu_sigma;
//  real<lower=0> sd_sigma;
  // Dispersion priors
  real<lower=0> mu_dispersion;
  real<lower=0> sd_dispersion;
  // Tolerance values
  real<lower=0> tolerance_expected;
  real<lower=0> tolerance_fishing;
}

transformed data {
  // Maximum number of observations (upper bound actually)
  int<lower=0> C = N * D * L * X * X;
  // Matrix version of movement index
  matrix[X, X] movement_matrix = to_matrix(movement_index);
  // Simplex dimensions
  array[6] int simplex_dims_mean = simplex_dims(movement_index);
  array[6] int simplex_dims_time = simplex_dims(movement_index, model_time, T);
  array[6] int simplex_dims_term = simplex_dims(movement_index, model_term, K);
  array[6] int simplex_dims_size = simplex_dims(movement_index, model_size, L);
  // Declare tags released
  array[N - 1, L] vector[X] tags_released = assemble_tags_released(tags);
  // Declare tags transpose (permute x and y)
  array[N - 1, D, L, X, X] int tags_transpose = assemble_tags_transpose(tags);
  // Declare movement possible values
  array[D] matrix[X, X] movement_possible = assemble_movement_possible(
    movement_matrix,
    D
  );
  // Declare partial sum index
  array[N - 1] int partial_sum_index = rep_array(1, N - 1);
  // Declare partial sum grainsize
  int grainsize = 1;
}

parameters {
  // Stepwise movement mean simplexes
  array[simplex_dims_mean[1]] simplex[1] m1;
  array[simplex_dims_mean[2]] simplex[2] m2;
  array[simplex_dims_mean[3]] simplex[3] m3;
  array[simplex_dims_mean[4]] simplex[4] m4;
  array[simplex_dims_mean[5]] simplex[5] m5;
  array[simplex_dims_mean[6]] simplex[6] m6;
  // Stepwise movement time simplexes
  array[simplex_dims_time[1]] simplex[1] t1;
  array[simplex_dims_time[2]] simplex[2] t2;
  array[simplex_dims_time[3]] simplex[3] t3;
  array[simplex_dims_time[4]] simplex[4] t4;
  array[simplex_dims_time[5]] simplex[5] t5;
  array[simplex_dims_time[6]] simplex[6] t6;
  // Stepwise movement term simplexes
  array[simplex_dims_term[1]] simplex[1] k1;
  array[simplex_dims_term[2]] simplex[2] k2;
  array[simplex_dims_term[3]] simplex[3] k3;
  array[simplex_dims_term[4]] simplex[4] k4;
  array[simplex_dims_term[5]] simplex[5] k5;
  array[simplex_dims_term[6]] simplex[6] k6;
  // Stepwise movement size simplexes
  array[simplex_dims_size[1]] simplex[1] l1;
  array[simplex_dims_size[2]] simplex[2] l2;
  array[simplex_dims_size[3]] simplex[3] l3;
  array[simplex_dims_size[4]] simplex[4] l4;
  array[simplex_dims_size[5]] simplex[5] l5;
  array[simplex_dims_size[6]] simplex[6] l6;
  // Instantaneous stepwise rates
  array[T] vector<lower=0>[X] fishing_step;
  vector<lower=0>[X] natural_mortality_step;
  real<lower=0> ongoing_loss_step;
  // Fractional (per tag) stepwise rates
  vector<lower=0, upper=1>[X] reporting_step;
  real<lower=0, upper=1> initial_loss_step;
//  // Selectivity (per fish)
//  vector<lower=0, upper=1>[L - 1] selectivity_short;
//  // Seasonal fishing weights
//  array[X] simplex[K] fishing_weight_transpose;
//  // Autoregression priors
//  vector<lower=0, upper=1>[form > 0 ? X : 0] autoregress;
//  vector<lower=0>[form > 0 ? X : 0] sigma;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Stepwise movement rate components
  matrix<lower=0, upper=1>[X, X] movement_step_mean;
  array[T, X] matrix<lower=0, upper=1>[X, X] movement_tran_time;
  array[K, X] matrix<lower=0, upper=1>[X, X] movement_tran_term;
  array[L, X] matrix<lower=0, upper=1>[X, X] movement_tran_size;
  // Stepwise movement rate
  array[T, K, L] matrix<lower=0, upper=1>[X, X] movement_step;
  // Stepwise survival rate
  array[T, K, L] vector<lower=0, upper=1>[X] survival_step;
  // Stepwise transition rate
  array[N, L] matrix<lower=0, upper=1>[X, X] transition_step;
  // Stepwise observation rate
  array[N, L] vector<lower=0, upper=1>[X] observation_step;

  // Declare instantaneous rates
  array[T] vector<lower=0>[X] fishing_rate;
  vector<lower=0>[X] natural_mortality_rate = natural_mortality_step * J;
  real<lower=0> ongoing_loss_rate = ongoing_loss_step * J;
  // Declare fractional (per tag) rates
  vector<lower=0, upper=1>[X] reporting_rate = reporting_step;
  real<lower=0, upper=1> initial_loss_rate = initial_loss_step;
//  // Selectivity
//  vector<lower=0, upper=1>[S] selectivity = append_row(selectivity_short, 1.0);
//  // Declare fishing weight
//  array[K] vector<lower=0, upper=1>[X] fishing_weight;

  // Assemble stepwise movement rate components
  movement_step_mean = assemble_movement_step_mean(
    m1, m2, m3, m4, m5, m6,
    movement_index
  );
  movement_tran_time = assemble_movement_transformation(
    t1, t2, t3, t4, t5, t6,
    movement_index,
    model_time,
    T
  );
  movement_tran_term = assemble_movement_transformation(
    k1, k2, k3, k4, k5, k6,
    movement_index,
    model_term,
    K
  );
  movement_tran_size = assemble_movement_transformation(
    l1, l2, l3, l4, l5, l6,
    movement_index,
    model_size,
    L
  );
  // Assemble stepwise movement rate [T, K, L][X, X]
  movement_step = assemble_movement_step(
    movement_step_mean,
    movement_tran_time,
    movement_tran_term,
    movement_tran_size,
    model_time,
    model_term,
    model_size
  );
//  // Assemble fishing weight [K][X]
//  fishing_weight = assemble_fishing_weight(fishing_weight_transpose);
  // Assemble stepwise survival rate [T, K, L][X]
  survival_step = assemble_survival_step(
    fishing_step,
//    fishing_weight,
//    selectivity,
    natural_mortality_step,
    ongoing_loss_step,
    K, L
  );
  // Assemble stepwise transition rate [N, L][X, X]
  transition_step = assemble_transition_step(
    movement_step,
    survival_step
  );
  // Assemble stepwize observation rate [N, L][X]
  observation_step = assemble_observation_step(
    fishing_step,
//    fishing_weight,
//    selectivity,
    reporting_step,
    K, L
  );
  // Assemble fishing rate [T][X]
  fishing_rate = assemble_fishing_rate(fishing_step, J);
}

model {
  // Stepwise movement mean priors
  diagonal(movement_step_mean) ~ normal(mu_movement_diag, sd_movement_diag);
  // Stepwise movement transformation priors
  if (model_time) {}
  if (model_term) {}
  if (model_size) {}
  // Fishing rate prior
  for (t in 1:T) {
    fishing_rate[t] ~ normal(
      mu_fishing_rate[t],
      mu_fishing_rate[t] * cv_fishing_rate + tolerance_fishing
    );
  }
  // Natural mortality rate prior
  natural_mortality_rate ~ normal(
    mu_natural_mortality_rate,
    sd_natural_mortality_rate
  );
  // Reporting rate prior
  reporting_rate ~ normal(mu_reporting_rate, sd_reporting_rate);
  // Ongoing loss rate prior
  ongoing_loss_rate ~ normal(mu_ongoing_loss_rate, sd_ongoing_loss_rate);
  // Initial loss rate prior
  initial_loss_rate ~ normal(mu_initial_loss_rate, sd_initial_loss_rate);
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, sd_dispersion);
  // Use reduce sum
  target += reduce_sum(
    partial_sum_lpmf, // Partial sum function
    partial_sum_index, // Each element corresponds to a summand
    grainsize, // Set to one to leave partitioning up to scheduler
    N, D, L, X, // Arguments shared by every term...
    tags_transpose,
    tags_released,
    transition_step,
    observation_step,
    movement_possible,
    initial_loss_step,
    tolerance_expected,
    dispersion
  );
}

generated quantities {
  // Movement mean
  matrix[X, X] movement_mean = matrix_power(movement_step_mean, J);
  // Movement facets
  array[T] matrix[X, X] movement_time;
  array[K] matrix[X, X] movement_term;
  array[L] matrix[X, X] movement_size;
  // Populate movement time
  movement_time = assemble_movement_rate(
    movement_step_mean,
    movement_tran_time,
    model_time,
    J
  );
  // Populate movement term
  movement_term = assemble_movement_rate(
    movement_step_mean,
    movement_tran_term,
    model_term,
    H
  );
  // Populate movement size
  movement_size = assemble_movement_rate(
    movement_step_mean,
    movement_tran_size,
    model_size,
    J
  );
}
