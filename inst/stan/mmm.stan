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
  vector<lower=0, upper=1>[X] mu_retention_rate;
  array[T] vector<lower=0>[X] mu_fishing_rate;
  vector<lower=0>[X] mu_mortality_rate;
  vector<lower=0>[X] mu_reporting_rate;
  vector<lower=0, upper=1>[S - 1] mu_selectivity;
  real<lower=0, upper=1> mu_initial_loss_rate; // Fraction
  real<lower=0> mu_ongoing_loss_rate; // Instantaneous
  real<lower=0> mu_dispersion;
  // Prior coefficients of variation
  real<lower=0> cv_retention_rate;
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
  // Declare movement possible transpose values
  array[L] matrix[X, X] movement_possible_transpose;
  // Populate movement possible transpose
  movement_possible_transpose = assemble_movement_possible_transpose(
    movement_index,
    L
  );
}

parameters {
  // Movement parameters
  array[P] real movement_parameter_mean_step; // Mean across dimensions
  array[T, P] real movement_parameter_time_step; // Deviation by 'year'
//  array[I, P] real movement_parameter_term_step; // Deviation by 'season'
//  array[S, P] real movement_parameter_size_step; // Deviation by 'size'
  // Stepwise instantaneous rate parameters
//  array[T] vector<lower=0>[X] fishing_step; // Instantaneous
//  vector<lower=0>[X] mortality_step; // Instantaneous
//  real<lower=0> ongoing_loss_step; // Instantaneous
  // Stepwise fractional (per tag) rate parameters
//  vector<lower=0, upper=1>[X] reporting_step; // Fraction (per tag)
//  real<lower=0, upper=1> initial_loss_step; // Fraction (per tag)
  // Selectivity
//  vector<lower=0, upper=1>[S - 1] selectivity_short;
  // Autoregressive process parameters
  real<lower=0, upper=1> ar_movement_time_deviation;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Fix non-movement parameters

/**
  // Stepwise instantaneous rate parameters
  array[T] vector<lower=0>[X] fishing_step = rep_array(rep_vector(0.015,X),T);
  vector<lower=0>[X] mortality_step = rep_vector(0.025, X);
  real<lower=0> ongoing_loss_step = 0.005;
  // Stepwise fractional (per tag) rate parameters
  vector<lower=0, upper=1>[X] reporting_step = rep_vector(0.5, X);
  real<lower=0, upper=1> initial_loss_step = 0.1;
  // Selectivity
  vector<lower=0, upper=1>[S - 1] selectivity_short = rep_vector(1.0, S - 1);
  // Fishing term weight
  array[I] vector[X] fishing_term = rep_array(rep_vector(1.0, X), I);
*/


  // Declare step values (for computation)
  array[N, S] matrix<lower=0, upper=1>[X,X] movement_step = rep_array(rep_matrix(0.0,X,X),N,S);
//  array[N, S] matrix<lower=0, upper=1>[X,X] survival_step = rep_array(rep_matrix(0.0,X,X),N,S);
//  array[N, S] matrix<lower=0, upper=1>[X,X] observed_step = rep_array(rep_matrix(0.0,X,X),N,S);
  // Declare movement values (for priors)
  matrix<lower=0, upper=1>[X, X] movement_mean_step = rep_matrix(0.0, X, X);
  matrix<lower=0, upper=1>[X, X] movement_mean_rate = rep_matrix(0.0, X, X);
  array[T] matrix<lower=0, upper=1>[X,X] movement_time_step = rep_array(rep_matrix(0.0, X, X), T);
  array[T] matrix<lower=0, upper=1>[X,X] movement_time_rate = rep_array(rep_matrix(0.0, X, X), T);
  matrix<lower=0, upper=1>[X, X] movement_time_product = rep_matrix(0.0, X, X);
//  array[T] matrix<lower=-1, upper=1>[X,X] movement_time_deviation=rep_array(rep_matrix(0.0,X,X),T);
//  array[I] matrix<lower=-1, upper=1>[X,X] movement_term_deviation=rep_array(rep_matrix(0.0,X,X),I);
//  array[S] matrix<lower=-1, upper=1>[X,X] movement_size_deviation=rep_array(rep_matrix(0.0,X,X),S);
  // Declare instantaneous rates (for priors)
//  array[T] vector[X] fishing_rate = rep_array(rep_vector(0.0, X), T);
//  vector<lower=0>[X] mortality_rate = mortality_step * I; // Instantaneous
//  real<lower=0> ongoing_loss_rate = ongoing_loss_step * I; // Instantaneous
  // Declare fractional (per tag) rates
//  vector<lower=0, upper=1>[X] reporting_rate = reporting_step; // Fraction
//  real<lower=0, upper=1> initial_loss_rate = initial_loss_step; // Fraction
  // Selectivity
//  vector<lower=0, upper=1>[S] selectivity = append_row(selectivity_short, 1.0);
  // Assemble movement step [N, S] [X, X]
  movement_step = assemble_movement_step(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    rep_array(0.0, 1, P), // movement_parameter_term_step,
    rep_array(0.0, 1, P), // movement_parameter_size_step,
    movement_index
  );
  // Assemble movement mean step [X, X]
  movement_mean_step = assemble_movement_mean_rate(
    movement_parameter_mean_step,
    movement_index,
    1 // nterm = 1 for movement mean step
  );
  // Assemble movement mean rate [X, X]
  movement_mean_rate = assemble_movement_mean_rate(
    movement_parameter_mean_step,
    movement_index,
    I // nterm = I for movement mean rate
  );
  // Assemble movement time step [T] [X, X]
//  movement_time_step = assemble_movement_time_rate(
//    movement_parameter_mean_step,
//    movement_parameter_time_step,
//    movement_index,
//    1 // nterm = 1 for movement time step
//  );
  // Assemble movement time rate [T] [X, X]
  movement_time_rate = assemble_movement_time_rate(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    movement_index,
    I
  );
  // Assemble movement time product [X, X]
  // movement_time_product = assemble_movement_time_product(movement_time_step);

  /**

  // Assemble movement time deviation [T] [X, X]
  movement_time_deviation = assemble_movement_time_deviation(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    movement_index,
    I
  );
  // Assemble movement term deviation [I] [X, X]
  movement_term_deviation = assemble_movement_term_deviation(
    movement_parameter_mean_step,
    movement_parameter_term_step,
    movement_index,
    1 // nterm = 1 for movement term deviations
  );
  // Assemble movement size deviation [S] [X, X]
  movement_size_deviation = assemble_movement_size_deviation(
    movement_parameter_mean_step,
    movement_parameter_size_step,
    movement_index,
    I
  );
  // Assemble survival step [N, S] [X, X]
  survival_step = assemble_survival_step(
    fishing_step,
    fishing_term,
    selectivity,
    mortality_step,
    ongoing_loss_step
  );
  // Assemble observal step [N, S] [X, X]
  observed_step = assemble_observed_step(
    fishing_step,
    fishing_term,
    selectivity,
    reporting_step
  );
  */

  print("movement_mean_rate: ", movement_mean_rate);

  /**
  // Assemble fishing rate [T] [X]
  fishing_rate = assemble_fishing_rate(
    fishing_step,
    T,
    I
  );
  // Assemble fishing term deviation [T, I] [X]
  fishing_term_deviation = assemble_fishing_term_deviation(
    fishing_step,
    fishing_rate,
    I
  );
  */
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
  for (n in 1:(N - 1)) { // Study step
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
        * movement_step[n + l - 2, s]; // Prevous step; square matrix [X, X]
        // Compute predicted
        predicted[n, s, l] = abundance[n, s, l] // Square matrix [X, X]
//        * observed_step[n + l - 1, s]; // Current step; diagonal matrix [X, X]
        * 0.05;
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

//  // Movement time prior: matrix product constraint
//  to_vector(movement_time_product) ~ normal(
//    to_vector(matrix_power(movement_mean_step, T)),
//    0.01
//  );
  // Movement time prior: AR1 process
  to_vector(movement_time_rate[1]) ~ normal(
    to_vector(movement_mean_rate),
    cv_movement_time_deviation
    * to_vector(movement_mean_rate)
    + tolerance_movement
  );
  for (t in 2:T) {
    to_vector(movement_time_rate[t]) ~ normal(
      to_vector(movement_mean_rate)
      + ar_movement_time_deviation
      * to_vector(movement_time_rate[t - 1]),
      cv_movement_time_deviation
      * to_vector(movement_mean_rate)
      + tolerance_movement
    );
  }


  // Movement retention prior
  /**
  diagonal(movement_mean_rate) ~ normal(
    mu_retention_rate,
    cv_retention_rate * mu_retention_rate
  );
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
*/

  /**
  // Fishing rate prior
  for (t in 1:T) {
    fishing_rate[t] ~ normal(
      mu_fishing_rate[t],
      cv_fishing_rate * mu_fishing_rate[t] + tolerance_fishing
    );
  }
  // Fishing term deviation prior across times
  for (t in 2:T) {
    for (i in 1:I) {
      fishing_term_deviation[t, i] ~ normal(
        fishing_term_deviation[t - 1, i],
        tolerance_fishing
      );
    }
  }
  // Fishing term deviation prior across terms
  for (t in 1:T) {
    for (i in 2:I) {
      fishing_term_deviation[t, i] ~ normal(
        fishing_term_deviation[t, i - 1],
        sd_fishing_term_deviation
      );
    }
    fishing_term_deviation[t, 1] ~ normal(
      fishing_term_deviation[t, I],
      sd_fishing_term_deviation
    );
  }
  */

  /**
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
  /**

  // Declare movement values
  array[T] matrix[X,X] movement_time_rate = rep_array(rep_matrix(0.0,X,X),T);
  array[I] matrix[X,X] movement_term_rate = rep_array(rep_matrix(0.0,X,X),I);
  array[S] matrix[X,X] movement_size_rate = rep_array(rep_matrix(0.0,X,X),S);
  // Declare fishing values
  array[N, S, X] real fishing_step = rep_array(0.0, N, S, X);
  array[T, X] real fishing_time_rate = rep_array(0.0, T, X);
  array[I, X] real fishing_term_rate = rep_array(0.0, I, X);
  array[S, X] real fishing_size_rate = rep_array(0.0, S, X);

  // Assemble movement term rate [I] [X, X]
  movement_term_rate = assemble_movement_term_rate(
    movement_parameter_mean_step,
    movement_parameter_term_step,
    movement_index,
    1 // nterm = 1 for movement term rates
  );
  // Assemble movement size rate [S] [X, X]
  movement_size_rate = assemble_movement_size_rate(
    movement_parameter_mean_step,
    movement_parameter_size_step,
    movement_index,
    I
  );
  // Assemble fishing mortality step [N, S, X]
  fishing_step = assemble_fishing_step(
    fishing_parameter_mean_step,
    fishing_parameter_time_step,
    fishing_parameter_term_step,
    fishing_parameter_size_step
  );
  // Assemble fishing time rate [T, X]
  fishing_time_rate = assemble_fishing_time_rate(
    fishing_parameter_mean_step,
    fishing_parameter_time_step,
    I
  );
  // Assemble fishing term rate [I, X]
  fishing_term_rate = assemble_fishing_term_rate(
    fishing_parameter_mean_step,
    fishing_parameter_term_step,
    1 // nterm = 1 for fishing term rates
  );
  // Assemble fishing size rate [S, X]
  fishing_size_rate = assemble_fishing_size_rate(
    fishing_parameter_mean_step,
    fishing_parameter_size_step,
    I
  );

  */
}
