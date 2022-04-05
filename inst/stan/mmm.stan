functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=0> X; // Number of geographic regions
  int<lower=0> T; // Number of times (usually years; release only)
  int<lower=0> I; // Number of terms per unit of time
  int<lower=0> S; // Number of size classes
  int<lower=0> N; // Number of release steps ((T * I) - 1)
  int<lower=0> L; // Number of maximum steps at liberty
  int<lower=0> P; // Number of movement rate mean parameters
  // Tag data
  array[N, S, X] int<lower=0> tags_released;
  array[N, S, X, L, X] int<lower=0> tags_recovered;
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
  array[T, X] real mu_fishing_deviation_time_rate;
  // Prior standard deviations
  array[X] real<lower=0> sd_mortality_rate;
  array[X] real<lower=0> sd_reporting_rate;
  array[X] real<lower=0> sd_fishing_mean_rate;
  real<lower=0> sd_initial_loss_rate;
  real<lower=0> sd_ongoing_loss_rate;
  real<lower=0> sd_dispersion;
  // Movement prior coefficients of variation
  real<lower=0> cv_movement_deviation_time_rate;
  real<lower=0> cv_movement_deviation_term_rate;
  real<lower=0> cv_movement_deviation_size_rate;
  // Fishing prior coefficients of variation
  real<lower=0> cv_fishing_deviation_time_rate;
  real<lower=0> cv_fishing_deviation_term_rate;
  real<lower=0> cv_fishing_deviation_size_rate;
  // Fudge values
  real<lower=0> expected_fudge;
}

transformed data {
  int<lower=0> C = 0; // Number of observations
  // Count observations
  for (n in 1:N) { // Released step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        if (tags_released[n, s, x] > 0) {
          C += (min(N - n + 2, L) - 1) * X; // Realized steps liberty
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
  // Define step values (for computation)
  array[N, S] matrix[X, X] movement_step;
  array[N, S, X] real fishing_step;
  array[N, S, X] real survival_step;
  // Define movement rate values (for priors)
  matrix[X, X] movement_mean_rate;
  array[T] matrix[X, X] movement_deviation_time_rate;
  array[I] matrix[X, X] movement_deviation_term_rate;
  array[S] matrix[X, X] movement_deviation_size_rate;
  // Define fishing rate values (for priors)
  array[X] real fishing_mean_rate;
  array[T, X] real fishing_deviation_time_rate;
  array[I, X] real fishing_deviation_term_rate;
  array[S, X] real fishing_deviation_size_rate;
  // Define remaining values
  array[N, S, X, L, X] real abundance;
  array[C] int  observed;
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
  // Assemble movement mean rate
  movement_mean_rate = assemble_movement_rate(
    movement_parameter_mean_step,
    movement_index,
    I
  );
  // Assemble movement deviation time rate [T] [X, X]
  movement_deviation_time_rate = assemble_movement_deviation(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    movement_index,
    I
  );
  // Assemble movement deviation term rate [I] [X, X]
  movement_deviation_term_rate = assemble_movement_deviation(
    movement_parameter_mean_step,
    rep_array(0.0, 1, P),
    movement_parameter_term_step,
    rep_array(0.0, 1, P),
    movement_index,
    I
  );
  // Assemble movement deviation size rate [S] [X, X]
  movement_deviation_size_rate = assemble_movement_deviation(
    movement_parameter_mean_step,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
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
  // Assemble fishing mean rate
  fishing_mean_rate = assemble_fishing_rate(
    fishing_parameter_mean_step,
    I
  );
  // Assemble fishing deviation time rate [T, X]
  fishing_deviation_time_rate = assemble_fishing_deviation(
    fishing_parameter_mean_step,
    fishing_parameter_time_step,
    rep_array(0.0, 1, X),
    rep_array(0.0, 1, X),
    I
  );
  // Assemble fishing deviation term rate [I, X]
  fishing_deviation_term_rate = assemble_fishing_deviation(
    fishing_parameter_mean_step,
    rep_array(0.0, 1, X),
    fishing_parameter_term_step,
    rep_array(0.0, 1, X),
    I
  );
  // Assemble fishing deviation size rate [S, X]
  fishing_deviation_size_rate = assemble_fishing_deviation(
    fishing_parameter_mean_step,
    rep_array(0.0, 1, X),
    rep_array(0.0, 1, X),
    fishing_parameter_size_step,
    I
  );
  // Compute survival step [N, S, X]
  for (n in 1:N) { // Released step
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
  // Populate abundance released [N, S, X]
  for (n in 1:N) { // Released step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Region
        abundance[n, s, x, 1, x] = (1 - initial_loss_rate)
        * tags_released[n, s, x];
      }
    }
  }
  // Compute expected recoveries
  for (n in 1:N) { // Released step
    for (s in 1:S) { // Released size
      for (w in 1:X) { // Released region
        if (tags_released[n, s, w] > 0) {
          for (l in 2:min(N - n + 2, L)) { // Populate abundance
            for (y in 1:X) { // Current region
              for (x in 1:X) { // Previous region
                abundance[n, s, w, l, y] += abundance[n, s, w, l - 1, x]
                * survival_step[n + l - 2, s, x] // Previous step
                * movement_step[n + l - 2, s, x, y]; // Previous step
              } // End x
            } // End y
          } // End l
          for (l in 2:min(N - n + 2, L)) { // Populate 1D arrays
            for (y in 1:X) { // Current region
              observed[count] = tags_recovered[n, s, w, l, y];
              expected[count] = abundance[n, s, w, l, y]
              * (1 - exp(fishing_step[n + l - 1, s, y]))
              * reporting_step[y]
              + expected_fudge;
              count += 1;
            } // End y
          } // End l
        } // End if
      } // End w
    } // End s
  } // End n
  // Movement deviation time rate prior
  for (t in 2:T) {
    for (y in 1:X) {
      for (x in 1:X) {
        movement_deviation_time_rate[t, x, y] ~ normal(
          movement_deviation_time_rate[t - 1, x, y],
          cv_movement_deviation_time_rate * movement_mean_rate[x, y]
        );
      }
    }
  }
  // Movement deviation term rate prior
  for (i in 2:I) {
    for (y in 1:X) {
      for (x in 1:X) {
        movement_deviation_term_rate[i, x, y] ~ normal(
          movement_deviation_term_rate[i - 1, x, y],
          cv_movement_deviation_term_rate * movement_mean_rate[x, y]
        );
      }
    }
  }
  for (y in 1:X) {
    for (x in 1:X) {
      movement_deviation_term_rate[1, x, y] ~ normal(
        movement_deviation_term_rate[I, x, y],
        cv_movement_deviation_term_rate * movement_mean_rate[x, y]
      );
    }
  }
  // Movement deviation size rate prior
  for (s in 2:S) {
    for (y in 1:X) {
      for (x in 1:X) {
        movement_deviation_size_rate[s, x, y] ~ normal(
          movement_deviation_size_rate[s - 1, x, y],
          cv_movement_deviation_size_rate * movement_mean_rate[x, y]
        );
      }
    }
  }
  // Fishing mean rate prior
  fishing_mean_rate ~ normal(mu_fishing_mean_rate, sd_fishing_mean_rate);
  // Fishing deviation time rate prior
  for (t in 1:T) {
    for (x in 1:X) {
      fishing_deviation_time_rate[t, x] ~ normal(
        mu_fishing_deviation_time_rate[t, x],
        cv_fishing_deviation_time_rate * fishing_mean_rate[x]
      );
    }
  }
  // Fishing deviation term rate prior
  for (i in 2:I) {
    for (x in 1:X) {
      fishing_deviation_term_rate[i, x] ~ normal(
        fishing_deviation_term_rate[i - 1, x],
        cv_fishing_deviation_term_rate * fishing_mean_rate[x]
      );
    }
  }
  for (x in 1:X) {
    fishing_deviation_term_rate[1, x] ~ normal(
      fishing_deviation_term_rate[I, x],
      cv_fishing_deviation_term_rate * fishing_mean_rate[x]
    );
  }
  // Fishing deviation size rate prior
  for (s in 2:S) {
    for (x in 1:X) {
      fishing_deviation_size_rate[s, x] ~ normal(
        fishing_deviation_size_rate[s - 1, x],
        cv_fishing_deviation_size_rate * fishing_mean_rate[x]
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
  // Define movement values
  array[T] matrix[X,X] movement_time_rate = rep_array(rep_matrix(0.0,X,X),T);
  array[I] matrix[X,X] movement_term_rate = rep_array(rep_matrix(0.0,X,X),I);
  array[S] matrix[X,X] movement_size_rate = rep_array(rep_matrix(0.0,X,X),S);
  // Define fishing values
  array[T, X] real fishing_time_rate = rep_array(0.0, T, X);
  array[I, X] real fishing_term_rate = rep_array(0.0, I, X);
  array[S, X] real fishing_size_rate = rep_array(0.0, S, X);
  // Assemble movement time rate
  movement_time_rate = assemble_movement_rate(
    movement_parameter_mean_step,
    movement_parameter_time_step,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    movement_index,
    I
  );
  // Assemble movement term rate
  movement_term_rate = assemble_movement_rate(
    movement_parameter_mean_step,
    rep_array(0.0, 1, P),
    movement_parameter_term_step,
    rep_array(0.0, 1, P),
    movement_index,
    I
  );
  // Assemble movement size rate
  movement_size_rate = assemble_movement_rate(
    movement_parameter_mean_step,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    movement_parameter_size_step,
    movement_index,
    I
  );
  // Assemble fishing time rate
  fishing_time_rate = assemble_fishing_rate(
    fishing_parameter_mean_step,
    fishing_parameter_time_step,
    rep_array(0.0, 1, X),
    rep_array(0.0, 1, X),
    I
  );
  // Assemble fishing term rate
  fishing_term_rate = assemble_fishing_rate(
    fishing_parameter_mean_step,
    rep_array(0.0, 1, X),
    fishing_parameter_term_step,
    rep_array(0.0, 1, X),
    I
  );
  // Assemble fishing size rate
  fishing_size_rate = assemble_fishing_rate(
    fishing_parameter_mean_step,
    rep_array(0.0, 1, X),
    rep_array(0.0, 1, X),
    fishing_parameter_size_step,
    I
  );
}
