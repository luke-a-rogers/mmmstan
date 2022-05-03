functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=2> N; // Number of model steps (months) in the study
  int<lower=1> T; // Number of times (years) in the study
  int<lower=1> I; // Number of unique terms (seasons) in the study
  int<lower=1> S; // Number of size classes
  int<lower=2> L; // Number of maximum model steps at liberty
  int<lower=2> X; // Number of geographic regions
  int<lower=1> P; // Number of movement step parameters in one [N, S] slice
  int<lower=1> K; // Matrix power to convert movement step to movement rate
  // Tag data
  array[N - 1, S, L, X, X] int<lower=0> tags;
  // Indexes
  array[X, X] int<lower=0, upper=1> movement_index;
  array[N] int<lower=1, upper=T> step_to_time;
  array[N] int<lower=1, upper=I> step_to_term;
  // Structure
  int<lower=0, upper=1> nest_term_within_time;
  // Prior means
  array[T] vector<lower=0>[X] mu_fishing_rate;
  vector<lower=0>[X] mu_mortality_rate;
  vector<lower=0>[X] mu_reporting_rate;
  vector<lower=0, upper=1>[S - 1] mu_selectivity;
  real<lower=0, upper=1> mu_initial_loss_rate; // Fraction
  real<lower=0> mu_ongoing_loss_rate; // Instantaneous
//  real<lower=0> mu_random_walk_sigma;
  real<lower=0> mu_dispersion;
  // Prior coefficients of variation
  real<lower=0> cv_fishing_rate;
  real<lower=0> cv_mortality_rate;
  real<lower=0> cv_reporting_rate;
  real<lower=0> cv_selectivity;
  real<lower=0> cv_initial_loss_rate;
  real<lower=0> cv_ongoing_loss_rate;
//  real<lower=0> sd_random_walk_sigma;
  real<lower=0> cv_dispersion;
  // Tolerance values
  real<lower=0> tolerance_expected;
  real<lower=0> tolerance_movement;
  real<lower=0> tolerance_fishing;
}

transformed data {
  // Upper bound on number of observations
  int<lower=0> C = N * S * L * X * X;
  // Declare simplex dimensions
  array[8] int simplex_dimensions = assemble_simplex_dimensions(
    movement_index,
    I, S, 8 // Currently ndims <= 8
  );
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
  array[simplex_dimensions[7]] simplex[7] a7;
  array[simplex_dimensions[8]] simplex[8] a8;

  /**

  // Stepwise instantaneous rate parameters
  array[T] vector<lower=0>[X] fishing_step; // Instantaneous
  vector<lower=0>[X] mortality_step; // Instantaneous
  real<lower=0> ongoing_loss_step; // Instantaneous
  // Stepwise fractional (per tag) rate parameters
  vector<lower=0, upper=1>[X] reporting_step; // Fraction (per tag)
  real<lower=0, upper=1> initial_loss_step; // Fraction (per tag)
  // Selectivity (per fish)
  vector<lower=0, upper=1>[S - 1] selectivity_short;

  */
  // Random walk sigma
  //real<lower=0> random_walk_sigma;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Declare step values (for computation)
  array[I, S] matrix<lower=0, upper=1>[X, X] movement_step;
  array[I, S] matrix<lower=0, upper=1>[X, X] movement_rate;

  // Assemble movement step [I, S] [X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6, a7, a8,
    movement_index,
    I, S, 8 // Currently ndims <= 8
  );
  // Assemble movement rate [I, S] [X, X]
  movement_rate = assemble_movement_rate(
    movement_step,
    K // Matrix power to convert movement step to movement rate
  );

  // array[N, S] matrix<lower=0, upper=1>[X, X] survival_step = rep_array(rep_matrix(0.0,X,X),N,S);
  // array[N, S] matrix<lower=0, upper=1>[X, X] observed_step = rep_array(rep_matrix(0.0,X,X),N,S);
  // Declare movement means (for priors)
  // matrix<lower=0, upper=1>[X, X] movement_mean_step = rep_matrix(0.0, X, X);
  // matrix<lower=0, upper=1>[X, X] movement_mean_rate = rep_matrix(0.0, X, X);
  // Declare movement rates (for priors)
  // array[T] matrix<lower=0, upper=1>[X, X] movement_time_rate = rep_array(rep_matrix(0.0,X,X),T);
  // array[I] matrix<lower=0, upper=1>[X, X] movement_term_rate = rep_array(rep_matrix(0.0,X,X),I);
  // array[S] matrix<lower=0, upper=1>[X, X] movement_size_rate = rep_array(rep_matrix(0.0,X,X),S);
  // Declare movement deviations (for priors)
  // array[T] matrix<lower=-1, upper=1>[X, X] movement_time_deviation = rep_array(rep_matrix(0.0,X,X),T);
  // array[I] matrix<lower=-1, upper=1>[X, X] movement_term_deviation = rep_array(rep_matrix(0.0,X,X),I);
  // array[S] matrix<lower=-1, upper=1>[X, X] movement_size_deviation = rep_array(rep_matrix(0.0,X,X),S);
  // Declare fishing rates (for priors)
  // array[T] vector<lower=0>[X] fishing_rate = rep_array(rep_vector(0.0, X), T);
  // Declare instantaneous rates (for priors)
  // vector<lower=0>[X] mortality_rate = mortality_step * I; // Instantaneous
  // real<lower=0> ongoing_loss_rate = ongoing_loss_step * I; // Instantaneous
  // Declare fractional (per tag) rates
  // vector<lower=0, upper=1>[X] reporting_rate = reporting_step; // Fraction
  // real<lower=0, upper=1> initial_loss_rate = initial_loss_step; // Fraction
  // Selectivity
  //vector<lower=0, upper=1>[S] selectivity = append_row(selectivity_short, 1.0);
  // Fishing term weight
  // array[I] matrix[X, X] fishing_term = assemble_fishing_term(fishing_simplex);

  // Fix fishing_term
  // array[I] matrix[X, X] fishing_term = rep_array(
  //   diag_matrix(rep_vector(1.0, X)),
  //   I
  // );


  // Assemble survival step [N, S] [X, X]
  // survival_step = assemble_survival_step(
  //   fishing_step,
  //   fishing_term,
  //   selectivity,
  //   mortality_step,
  //   ongoing_loss_step
  // );
  // Assemble observed step [N, S] [X, X]
  // observed_step = assemble_observed_step(
  //   fishing_step,
  //   fishing_term,
  //   selectivity,
  //   reporting_step
  // );
  // Assemble movement mean step and rate [X, X]
  // movement_mean_step = assemble_movement_mean(movement_step, 1);
  // movement_mean_rate = assemble_movement_mean(movement_step, I);
  // Assemble movement time [T] term [I] and size [S] rates [X, X]
  // movement_time_rate = assemble_movement_time_rate(movement_step, T, I);
  // movement_term_rate = assemble_movement_term_rate(movement_step, T, I);
  // movement_size_rate = assemble_movement_size_rate(movement_step, T, I);
  // Assemble movement time deviation [T] [X, X]
  // movement_time_deviation = assemble_movement_time_deviation(
  //   movement_time_rate,
  //   movement_mean_rate
  // );
  // Assemble movement term deviation [I] [X, X]
  // movement_term_deviation = assemble_movement_term_deviation(
  //   movement_term_rate,
  //   movement_mean_step
  // );
  // Assemble movement size deviation [S] [X, X]
  // movement_size_deviation = assemble_movement_size_deviation(
  //   movement_size_rate,
  //   movement_mean_rate
  // );
  // Assemble fishing rate [T] [X]
  // fishing_rate = assemble_fishing_rate(fishing_step, I);
}

model {
  // Declare enumeration values
  array[N-1,S,L] matrix[X,X] abundance = rep_array(rep_matrix(0.0,X,X),N-1,S,L);
  array[N-1,S,L] matrix[X,X] predicted = rep_array(rep_matrix(0.0,X,X),N-1,S,L);
  array[C] int observed = rep_array(0, C);
  array[C] real expected = rep_array(0.0, C);
  // Declare index values
  int count;
  // Populate released abundance
  for (n in 1:(N - 1)) { // Released step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        abundance[n, s, 1, x, x] = tags[n, s, 1, x, x]; // Integer scalar
//        * (1 - initial_loss_step); // Real scalar
      }
    }
  }
  // Initialize count
  count = 0;
  // Compute expected recoveries
  for (n in 1:(N - 1)) { // Released step
    for (s in 1:S) { // Released size
      for (l in 2:min(N - n + 1, L)) { // Liberty step
        // Propagate abundance
        abundance[n, s, l] = abundance[n, s, l - 1]
//        * survival_step[n + l - 2, s] // Previous step; diagonal matrix [X, X]
        * 0.85
        * movement_step[step_to_term[n + l - 2], s]; // Prevous step; square matrix [X, X]
        // Compute predicted
        predicted[n, s, l] = abundance[n, s, l] // Square matrix [X, X]
        * 0.05;
//        * observed_step[n + l - 1, s]; // Current step; diagonal matrix [X, X]
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags[n, s, 1, x, x] > 0) { // Were any tags released?
              if (movement_possible_transpose[l][y, x] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags[n, s, l, x, y]; // Integer
                expected[count] = predicted[n, s, l, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End s
  } // End n

/**
  random_walk_sigma ~ normal(mu_random_walk_sigma, sd_random_walk_sigma);

  if (I > 1) {
    // Cyclic movement prior deviation
    if (nest_term_within_time == 1) {
      for (s in 1:S) {
        for (y in 1:X) {
          for (x in 1:X) {
            if (movement_index[x, y] == 1) {
              movement_step[1, s][x, y] ~ normal(
                movement_step[I, s][x, y],
                cv_movement_deviation * movement_step[I, s][x, y]
              );
            }
          }
        }
      }
    }
    // Movement deviation prior
    for (i in 2:I) {
      for (s in 1:S) {
        for (y in 1:X) {
          for (x in 1:X) {
            if (movement_index[x, y] == 1) {
              movement_step[i, s][x, y] ~ normal(
                movement_step[i - 1, s][x, y],
                cv_movement_deviation * movement_step[i - 1, s][x, y]
              );
            }
          }
        }
      }
    }
  }
  */

  /**
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
  */
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, cv_dispersion * mu_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);
}

generated quantities {
  // Parse which versions of movement rates etc to generate
  // movement_mean_rate
  // movement_time_rate
  // movement_term_rate
  // movement_size_rate
}
