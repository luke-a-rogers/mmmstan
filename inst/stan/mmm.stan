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
  array[N - 1, S, X] int<lower=0> tags_released;
  array[N - 1, S, X, L, X] int<lower=0> tags_recovered;
  // Movement index array
  array[X, X] int<lower=0, upper=1> movement_index;
  // Prior means
  array[X] real<lower=0> mu_mortality_rate;
  array[X] real<lower=0> mu_reporting_rate;
  array[X] real<lower=0> mu_fishing_mean_rate;
  real<lower=0, upper=1> mu_initial_loss_rate; // Fraction
  real<lower=0> mu_ongoing_loss_rate; // Instantaneous
  real<lower=0> mu_dispersion;
  // Fishing prior mean deviations by time
  array[T, X] real mu_fishing_time_deviation;
  // Prior standard deviations
  array[X] real<lower=0> sd_mortality_rate;
  array[X] real<lower=0> sd_reporting_rate;
  array[X] real<lower=0> sd_fishing_mean_rate;
  real<lower=0> sd_initial_loss_rate;
  real<lower=0> sd_ongoing_loss_rate;
  real<lower=0> sd_dispersion;
  // Movement prior coefficients of variation
  real<lower=0> cv_movement_time_deviation;
  real<lower=0> cv_movement_term_deviation;
  real<lower=0> cv_movement_size_deviation;
  // Fishing prior coefficients of variation
  real<lower=0> cv_fishing_time_deviation;
  real<lower=0> cv_fishing_term_deviation;
  real<lower=0> cv_fishing_size_deviation;
  // Tolerance values
  real<lower=0> tolerance_expected;
  real<lower=0> tolerance_movement;
}

transformed data {
  int<lower=0> C = 0; // Number of observations
  // Count observations
  for (n in 1:(N - 1)) { // Released step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        if (tags_released[n, s, x] > 0) {
          C += (min(N - n + 1, L) - 1) * X; // Realized steps liberty
        }
      }
    }
  }
}

parameters {
  // Movement parameters
  array[P] real movement_parameter_mean_step; // Mean across dimensions
  array[T, P] real movement_parameter_time_step; // Deviation by 'year'
  array[I, P] real movement_parameter_term_step; // Deviation by 'season'
  array[S, P] real movement_parameter_size_step; // Deviation by 'size'
  // Fishing mortality parameters
  array[X] real fishing_parameter_mean_step; // Mean across dimensions
  array[T, X] real fishing_parameter_time_step; // Deviation by 'year'
  array[I, X] real fishing_parameter_term_step; // Deviation by 'season'
  array[S, X] real fishing_parameter_size_step; // Deviation by 'size'
  // Natural mortality parameters
  array[X] real<lower=0> mortality_step; // Instantaneous
  // Tag reporting parameters
  array[X] real<lower=0, upper=1> reporting_step; // Fraction
  // Tag loss parameters
  real<lower=0, upper=1> initial_loss_rate; // Fraction
  real<lower=0> ongoing_loss_step; // Instantaneous
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Ongoing tag loss rate
  real<lower=0> ongoing_loss_rate = ongoing_loss_step * I;
  // Natural mortality rate
  array[X] real<lower=0> mortality_rate;
  // Tag reporting rate
  array[X] real<lower=0, upper=1> reporting_rate; // Fraction
  // Populate natural mortality rate
  for (x in 1:X) {
    mortality_rate[x] = mortality_step[x] * I;
  }
  // Populate tag reporting rate
  for (x in 1:X) {
    reporting_rate[x] = reporting_step[x]; // Fraction
  }
}

model {
  // Declare step values (for computation)
  array[N, S] matrix[X, X] movement_step;
  array[N, S, X] real fishing_step;
  array[N, S, X] real survival_step;
  // Declare movement rate values (for priors)
  matrix[X, X] movement_mean_step;
  matrix[X, X] movement_mean_rate;
  array[T] matrix[X, X] movement_time_deviation;
  array[I] matrix[X, X] movement_term_deviation;
  array[S] matrix[X, X] movement_size_deviation;
  // Declare fishing rate values (for priors)
  array[X] real fishing_mean_rate;
  array[T, X] real fishing_time_deviation;
  array[I, X] real fishing_term_deviation;
  array[S, X] real fishing_size_deviation;
  // Declare remaining values
  array[N - 1, S, X, L, X] real abundance;
  array[C] int observed;
  array[C] real expected;
  int count = 1;
  // Assemble movement step [N, S] [X, X]
  movement_step = assemble_movement_step(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    movement_parameter_term_step,
    movement_parameter_size_step,
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
    I
  );
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
  // Assemble fishing mortality step [N, S, X]
  fishing_step = assemble_fishing_step(
    fishing_parameter_mean_step,
    fishing_parameter_time_step,
    fishing_parameter_term_step,
    fishing_parameter_size_step
  );
  // Assemble fishing mean rate [X]
  fishing_mean_rate = assemble_fishing_mean_rate(
    fishing_parameter_mean_step,
    I
  );
  // Assemble fishing time deviation [T, X]
  fishing_time_deviation = assemble_fishing_time_deviation(
    fishing_parameter_mean_step,
    fishing_parameter_time_step,
    I
  );
  // Assemble fishing term deviation [I, X]
  fishing_term_deviation = assemble_fishing_term_deviation(
    fishing_parameter_mean_step,
    fishing_parameter_term_step,
    1 // nterm = 1 for fishing term deviations
  );
  // Assemble fishing size deviation [S, X]
  fishing_size_deviation = assemble_fishing_size_deviation(
    fishing_parameter_mean_step,
    fishing_parameter_size_step,
    I
  );
  // Compute survival step [N, S, X]
  for (n in 1:N) { // Study step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Region
        survival_step[n, s, x] = exp(
          -fishing_step[n, s, x]
          - mortality_step[x]
          - ongoing_loss_step
        );
      }
    }
  }
  // Populate abundance released [N - 1, S, X]
  for (n in 1:(N - 1)) { // Released step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Region
        abundance[n, s, x, 1, x] = (1 - initial_loss_rate)
        * tags_released[n, s, x];
      }
    }
  }
  // Compute expected recoveries
  for (n in 1:(N - 1)) { // Released step
    for (s in 1:S) { // Released size
      for (w in 1:X) { // Released region
        if (tags_released[n, s, w] > 0) {
          for (l in 2:min(N - n + 1, L)) { // Populate abundance
            for (y in 1:X) { // Current region
              for (x in 1:X) { // Previous region
                abundance[n, s, w, l, y] += abundance[n, s, w, l - 1, x]
                * survival_step[n + l - 2, s, x] // Previous step
                * movement_step[n + l - 2, s, x, y]; // Previous step
              } // End x
            } // End y
          } // End l
          for (l in 2:min(N - n + 1, L)) { // Populate 1D arrays
            for (y in 1:X) { // Current region
              observed[count] = tags_recovered[n, s, w, l, y];
              expected[count] = abundance[n, s, w, l, y]
              * (1 - exp(fishing_step[n + l - 1, s, y]))
              * reporting_step[y]
              + tolerance_expected;
              count += 1;
            } // End y
          } // End l
        } // End if
      } // End w
    } // End s
  } // End n
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
  // Fishing mean rate prior
  fishing_mean_rate ~ normal(mu_fishing_mean_rate, sd_fishing_mean_rate);
  // Fishing time deviation prior
  for (t in 1:T) {
    for (x in 1:X) {
      fishing_time_deviation[t, x] ~ normal(
        mu_fishing_time_deviation[t, x],
        cv_fishing_time_deviation * fishing_mean_rate[x]
      );
    }
  }
  // Fishing term deviation prior
  for (i in 2:I) {
    for (x in 1:X) {
      fishing_term_deviation[i, x] ~ normal(
        fishing_term_deviation[i - 1, x],
        cv_fishing_term_deviation * fishing_mean_rate[x]
      );
    }
  }
  for (x in 1:X) {
    fishing_term_deviation[1, x] ~ normal(
      fishing_term_deviation[I, x],
      cv_fishing_term_deviation * fishing_mean_rate[x]
    );
  }
  // Fishing size deviation prior
  for (s in 2:S) {
    for (x in 1:X) {
      fishing_size_deviation[s, x] ~ normal(
        fishing_size_deviation[s - 1, x],
        cv_fishing_size_deviation * fishing_mean_rate[x]
      );
    }
  }
  // Natural mortality rate prior
  mortality_rate ~ normal(mu_mortality_rate, sd_mortality_rate);
  // Reporting rate prior
  reporting_rate ~ normal(mu_reporting_rate, sd_reporting_rate);
  // Tag loss rate priors
  initial_loss_rate ~ normal(mu_initial_loss_rate, sd_initial_loss_rate);
  ongoing_loss_rate ~ normal(mu_ongoing_loss_rate, sd_ongoing_loss_rate);
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, sd_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed ~ neg_binomial_2(expected, dispersion);
}

generated quantities {
  // Declare movement values
  array[T] matrix[X,X] movement_time_rate = rep_array(rep_matrix(0.0,X,X),T);
  array[I] matrix[X,X] movement_term_rate = rep_array(rep_matrix(0.0,X,X),I);
  array[S] matrix[X,X] movement_size_rate = rep_array(rep_matrix(0.0,X,X),S);
  // Declare fishing values
  array[T, X] real fishing_time_rate = rep_array(0.0, T, X);
  array[I, X] real fishing_term_rate = rep_array(0.0, I, X);
  array[S, X] real fishing_size_rate = rep_array(0.0, S, X);
  // Assemble movement time rate [T] [X, X]
  movement_time_rate = assemble_movement_time_rate(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    movement_index,
    I
  );
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
}
