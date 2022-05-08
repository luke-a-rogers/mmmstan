functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=2> N; // Number of model steps (months) in the study
  int<lower=1> S; // Number of size or sex classes
  int<lower=2> L; // Number of maximum model steps at liberty
  int<lower=2> X; // Number of geographic regions
//  int<lower=1> T; // Number of fishing steps (years)
//  int<lower=1> W; // Number of fishing weight steps (seasons)
  int<lower=1> I; // Number of movement steps
  int<lower=1> D; // Number of movement sizes
  int<lower=1> P; // Number of movement step parameters in one [N, S] slice
  int<lower=1> K; // Matrix power to convert movement steps to rates
  // Index arrays
//  array[N] int<lower=1> n_to_t; // Model step to time (year) index
  array[N] int<lower=1> n_to_i; // Model step to movement step index
//  array[N] int<lower=1> n_to_w; // Model step to fishing weight step index
  array[S] int<lower=1> s_to_d; // Model size to movement size index
  // Tag data
  array[N - 1, S, L, X, X] int<lower=0> tags;
  // Indexes
  array[X, X] int<lower=0, upper=1> movement_index;
  // Prior means
  matrix<lower=0, upper=1>[X, X] mu_movement_mean; // Not used
  matrix<lower=0, upper=1>[X, X] mu_movement_total; // Not used
  real<lower=0> mu_dispersion;
  // Prior standard deviations
  matrix<lower=0, upper=1>[X, X] sd_movement_mean; // Not used
  matrix<lower=0, upper=1>[X, X] sd_movement_total; // Not used
  // Prior coefficients of variation
  real<lower=0> cv_dispersion;
  // Tolerance values
  real<lower=0> tolerance_expected;
}

transformed data {
  // Upper bound on number of observations
  int<lower=0> C = N * S * L * X * X;
  // Declare simplex dimensions
  array[6] int simplex_dimensions = assemble_simplex_dimensions(
    movement_index,
    I, D // Each set to one
  );
  // Declare tags released
  array[N - 1, S] vector[X] tags_released = assemble_tags_released(tags);
  // Declare tags transpose (permute x and y)
  array[N - 1, S, L, X, X] int tags_transpose = assemble_tags_transpose(tags);
  // Declare movement possible values
  array[L] matrix[X, X] movement_possible = assemble_movement_possible(
    movement_index,
    L
  );
  // Declare partial sum index
  array[N - 1] int partial_sum_index = rep_array(1, N - 1);
  // Declare partial sum grainsize
  int grainsize = 1;
}

parameters {
  // Movement simplexes
  array[simplex_dimensions[1]] simplex[1] a1;
  array[simplex_dimensions[2]] simplex[2] a2;
  array[simplex_dimensions[3]] simplex[3] a3;
  array[simplex_dimensions[4]] simplex[4] a4;
  array[simplex_dimensions[5]] simplex[5] a5;
  array[simplex_dimensions[6]] simplex[6] a6;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Fix some parameters
  real<lower=0, upper=1> initial_loss_step = 0.1;
  // Declare stepwise rates
  array[I, D] matrix<lower=0, upper=1>[X, X] movement_step; // [1, 1][X, X]
  array[N, S] vector<lower=0, upper=1>[X] survival_step; // [N, S][X]
  array[N, S] vector<lower=0, upper=1>[X] observed_step; // [N, S][X]

  // Assemble movement step [X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6,
    movement_index,
    I, D // Each set to one
  );
  // Assemble survival step
  survival_step = rep_array(rep_vector(0.85, X), N, S);
  // Assemble observed step
  observed_step = rep_array(rep_vector(0.05, X), N, S);
}

model {
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, cv_dispersion * mu_dispersion);
  // Use reduce sum
  target += reduce_sum(
    partial_sum_lpmf, // Partial sum function
    partial_sum_index, // Each element corresponds to a summand
    grainsize, // Set to one to leave partitioning up to scheduler
    N, S, L, X, // Arguments shared by every term...
    n_to_i,
    s_to_d,
    tags_transpose,
    tags_released,
    movement_step,
    survival_step,
    observed_step,
    movement_possible,
    initial_loss_step,
    tolerance_expected,
    dispersion
  );
}

generated quantities {
  // Assemble movement mean
  matrix[X, X] movement_mean = matrix_power(movement_step[1, 1], K); // [X, X]
  // Assemble movement total
  matrix[X, X] movement_total = matrix_power(movement_step[1, 1], N); // [X, X]
}
