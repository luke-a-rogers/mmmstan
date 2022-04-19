functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=0> X; // Number of geographic regions
  int<lower=0> T; // Number of times (usually years; release only)
  int<lower=0> I; // Number of terms per unit of time
  int<lower=0> S; // Number of size classes
  int<lower=0> N; // Number of study steps (T * I)
  int<lower=0> L; // Number of maximum steps at liberty
  int<lower=0> P; // Number of movement rate mean parameters
  // Tag data
  array[N - 1, S, L, X, X] int<lower=0> tags;
  // Movement index array
  array[X, X] int<lower=0, upper=1> movement_index;
  // Prior means
  array[T] vector<lower=0>[X] mu_fishing_rate;
  vector<lower=0>[X] mu_mortality_rate;
  vector<lower=0>[X] mu_reporting_rate;
  vector<lower=0, upper=1>[S - 1] mu_selectivity;
  real<lower=0, upper=1> mu_initial_loss_rate; // Fraction
  real<lower=0> mu_ongoing_loss_rate; // Instantaneous
  real<lower=0> mu_dispersion;
  // Prior coefficients of variation
  real<lower=0> cv_fishing_rate;
  real<lower=0> cv_mortality_rate;
  real<lower=0> cv_reporting_rate;
  real<lower=0> cv_selectivity;
  real<lower=0> cv_initial_loss_rate;
  real<lower=0> cv_ongoing_loss_rate;
  real<lower=0> cv_dispersion;
  // Movement prior coefficients of variation
  real<lower=0> cv_movement_time_deviation;
  real<lower=0> cv_movement_term_deviation;
  real<lower=0> cv_movement_size_deviation;
  // Tolerance values
  real<lower=0> tolerance_expected;
  real<lower=0> tolerance_movement;
  real<lower=0> tolerance_fishing;
}

transformed data {
  // Upper bound on number of observations
  int<lower=0> C = N * S * L * X * X;
  // Declare simplex dimensions
  array[6] int simplex_dimensions = assemble_simplex_dimensions(
    movement_index,
    T, I, S
  );
  // Declare partial sum index
  array[N - 1] int partial_sum_index = rep_array(1, N - 1);
  // Declare partial sum grainsize
  int grainsize = 1;
  // Declare movement possible transpose values
  array[L] matrix[X, X] movement_possible_transpose;
  // Populate movement possible transpose
  movement_possible_transpose = assemble_movement_possible_transpose(
    movement_index,
    L
  );
}

parameters {
  // Movement simplexes
  array[simplex_dimensions[1]] simplex[1] a1;
  array[simplex_dimensions[2]] simplex[2] a2;
  array[simplex_dimensions[3]] simplex[3] a3;
  array[simplex_dimensions[4]] simplex[4] a4;
  array[simplex_dimensions[5]] simplex[5] a5;
  array[simplex_dimensions[6]] simplex[6] a6;
  // Stepwise instantaneous rate parameters
  array[T] vector<lower=0>[X] fishing_step; // Instantaneous
  vector<lower=0>[X] mortality_step; // Instantaneous
  real<lower=0> ongoing_loss_step; // Instantaneous
  // Stepwise fractional (per tag) rate parameters
  vector<lower=0, upper=1>[X] reporting_step; // Fraction (per tag)
  real<lower=0, upper=1> initial_loss_step; // Fraction (per tag)
  // Selectivity
  vector<lower=0, upper=1>[S - 1] selectivity_short;
  // Fishing term weights
  // array[X] simplex[I] fishing_simplex;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Declare step values (for computation)
  array[N, S] matrix<lower=0, upper=1>[X, X] movement_step;
  array[N, S] matrix<lower=0, upper=1>[X, X] survival_step = rep_array(rep_matrix(0.0,X,X),N,S);
  array[N, S] matrix<lower=0, upper=1>[X, X] observed_step = rep_array(rep_matrix(0.0,X,X),N,S);
  // Declare movement means (for priors)
  matrix<lower=0, upper=1>[X, X] movement_mean_step = rep_matrix(0.0, X, X);
  matrix<lower=0, upper=1>[X, X] movement_mean_rate = rep_matrix(0.0, X, X);
  // Declare movement rates (for priors)
  array[T] matrix<lower=0, upper=1>[X, X] movement_time_rate = rep_array(rep_matrix(0.0,X,X),T);
  array[I] matrix<lower=0, upper=1>[X, X] movement_term_rate = rep_array(rep_matrix(0.0,X,X),I);
  array[S] matrix<lower=0, upper=1>[X, X] movement_size_rate = rep_array(rep_matrix(0.0,X,X),S);
  // Declare movement deviations (for priors)
  array[T] matrix<lower=-1, upper=1>[X, X] movement_time_deviation = rep_array(rep_matrix(0.0,X,X),T);
  array[I] matrix<lower=-1, upper=1>[X, X] movement_term_deviation = rep_array(rep_matrix(0.0,X,X),I);
  array[S] matrix<lower=-1, upper=1>[X, X] movement_size_deviation = rep_array(rep_matrix(0.0,X,X),S);
  // Declare fishing rates (for priors)
  array[T] vector<lower=0>[X] fishing_rate = rep_array(rep_vector(0.0, X), T);
  // Declare instantaneous rates (for priors)
  vector<lower=0>[X] mortality_rate = mortality_step * I; // Instantaneous
  real<lower=0> ongoing_loss_rate = ongoing_loss_step * I; // Instantaneous
  // Declare fractional (per tag) rates
  vector<lower=0, upper=1>[X] reporting_rate = reporting_step; // Fraction
  real<lower=0, upper=1> initial_loss_rate = initial_loss_step; // Fraction
  // Selectivity
  vector<lower=0, upper=1>[S] selectivity = append_row(selectivity_short, 1.0);
  // Fishing term weight
  // array[I] matrix[X, X] fishing_term = assemble_fishing_term(fishing_simplex);

  // Fix fishing_term
  array[I] matrix[X, X] fishing_term = rep_array(
    diag_matrix(rep_vector(1.0, X)),
    I
  );

  // Assemble movement step [N, S] [X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6,
    movement_index,
    T, I, S,
    tolerance_movement
  );
  // Assemble survival step [N, S] [X, X]
  survival_step = assemble_survival_step(
    fishing_step,
    fishing_term,
    selectivity,
    mortality_step,
    ongoing_loss_step
  );
  // Assemble observed step [N, S] [X, X]
  observed_step = assemble_observed_step(
    fishing_step,
    fishing_term,
    selectivity,
    reporting_step
  );
  // Assemble movement mean step and rate [X, X]
  movement_mean_step = assemble_movement_mean(movement_step, 1);
  movement_mean_rate = assemble_movement_mean(movement_step, I);
  // Assemble movement time [T] term [I] and size [S] rates [X, X]
  movement_time_rate = assemble_movement_time_rate(movement_step, T, I);
  movement_term_rate = assemble_movement_term_rate(movement_step, T, I);
  movement_size_rate = assemble_movement_size_rate(movement_step, T, I);
  // Assemble movement time deviation [T] [X, X]
  movement_time_deviation = assemble_movement_time_deviation(
    movement_time_rate,
    movement_mean_rate
  );
  // Assemble movement term deviation [I] [X, X]
  movement_term_deviation = assemble_movement_term_deviation(
    movement_term_rate,
    movement_mean_step
  );
  // Assemble movement size deviation [S] [X, X]
  movement_size_deviation = assemble_movement_size_deviation(
    movement_size_rate,
    movement_mean_rate
  );
  // Assemble fishing rate [T] [X]
  fishing_rate = assemble_fishing_rate(fishing_step, I);
}

model {

  // Priors --------------------------------------------------------------------

  // Movement time deviation prior
  for (t in 2:T) {
    for (y in 1:X) {
      for (x in 1:X) {
        if (movement_mean_rate[x, y] > 0) {
          movement_time_deviation[t, x, y] ~ normal(
            movement_time_deviation[t - 1, x, y],
            cv_movement_time_deviation * movement_mean_rate[x, y]
          );
        } else {
          movement_time_deviation[t, x, y] ~ normal(0.0, tolerance_movement);
        }
      }
    }
  }
  // Movement term deviation prior
  for (i in 2:I) {
    for (y in 1:X) {
      for (x in 1:X) {
        if (movement_mean_step[x, y] > 0) {
          movement_term_deviation[i, x, y] ~ normal(
            movement_term_deviation[i - 1, x, y],
            cv_movement_term_deviation * movement_mean_step[x, y]
          );
        } else {
          movement_term_deviation[i, x, y] ~ normal(0.0, tolerance_movement);
        }
      }
    }
  }
  for (y in 1:X) {
    for (x in 1:X) {
      if (movement_mean_rate[x, y] > 0) {
        movement_term_deviation[1, x, y] ~ normal(
          movement_term_deviation[I, x, y],
          cv_movement_term_deviation * movement_mean_rate[x, y]
        );
      } else {
          movement_term_deviation[1, x, y] ~ normal(0.0, tolerance_movement);
      }
    }
  }
  // Movement size deviation prior
  for (s in 2:S) {
    for (y in 1:X) {
      for (x in 1:X) {
        if (movement_mean_rate[x, y] > 0) {
          movement_size_deviation[s, x, y] ~ normal(
            movement_size_deviation[s - 1, x, y],
            cv_movement_size_deviation * movement_mean_rate[x, y]
          );
        } else {
          movement_size_deviation[s, x, y] ~ normal(0.0, tolerance_movement);
        }
      }
    }
  }
  // Fishing rate prior
  for (t in 1:T) {
    fishing_rate[t] ~ normal(
      mu_fishing_rate[t],
      cv_fishing_rate
      * mu_fishing_rate[t]
      + tolerance_fishing
    );
  }
  // Natural mortality rate prior
  mortality_rate ~ normal(
    mu_mortality_rate,
    cv_mortality_rate * mu_mortality_rate
  );
  // Reporting rate prior
  reporting_rate ~ normal(
    mu_reporting_rate,
    cv_reporting_rate * mu_reporting_rate
  );
  // Selectivity prior
  selectivity_short ~ normal(
    mu_selectivity,
    cv_selectivity * mu_selectivity
  );
  // Initial tag loss rate prior
  initial_loss_rate ~ normal(
    mu_initial_loss_rate,
    cv_initial_loss_rate * mu_initial_loss_rate
  );
  // Ongoing tag loss rate prior
  ongoing_loss_rate ~ normal(
    mu_ongoing_loss_rate,
    cv_ongoing_loss_rate * mu_ongoing_loss_rate
  );
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, cv_dispersion * mu_dispersion);

  // Sampling statement --------------------------------------------------------

  // Sampling statement (var = mu + mu^2 / dispersion)
  // observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);

  // Use reduce sum
  target += reduce_sum(
    partial_sum_lpmf, // Partial sum function
    partial_sum_index, // Each element corresponds to a summand
    grainsize, // Set to one to leave partitioning up to scheduler
    X, T, I, S, N, L, C,  // Arguments shared in every term...
    tags,
    movement_step,
    survival_step,
    observed_step,
    movement_possible_transpose,
    initial_loss_step,
    tolerance_expected,
    dispersion
  );
}
