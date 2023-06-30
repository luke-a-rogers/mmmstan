functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=2> N; // Number of model steps (months/quarters/years) in the study
  int<lower=2> D; // Maximum duration at large in model steps
  int<lower=1> L; // Number of size or sex classes
  int<lower=2> X; // Number of geographic regions
  int<lower=1> T; // Number of years in the study
  // Constants
  int<lower=1> K; // Number of steps (months/quarters/years) per year
  // Movement index arrays
  array[N] int<lower=1, upper=T> n_to_t; // Model step to year index
  // Tag data
  array[N - 1, D, L, X, X] int<lower=0> tags;
  // Movement index (will be paired with a matrix version)
  array[X, X] int<lower=0, upper=1> movement_index;
  // Movement step priors
  vector<lower=0, upper=1>[X] mu_movement_step_diag;
  vector<lower=0>[X] sd_movement_step_diag;
  // Fishing rate priors
  array[T] vector<lower=0>[X] mu_fishing_rate;
  real<lower=0> cv_fishing_rate;
  // Selectivity priors
  array[L - 1] vector<lower=0, upper=1>[X] mu_selectivity_short;
  vector<lower=0>[L > 1 ? 1 : 0] cv_selectivity;
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
  array[6] int simplex_dims = assemble_simplex_dims(movement_index);
  // Declare tags released
  array[N - 1, L] vector[X] tags_released = assemble_tags_released(tags);
  // Declare tags transpose (permute x and y)
  array[N - 1, D, L, X, X] int tags_transpose = assemble_tags_transpose(tags);
  // Declare movement possible values
  array[D] matrix[X, X] movement_possible = assemble_movement_possible(
    movement_matrix,
    D
  );
}

parameters {
  // Stepwise movement rate simplexes
  array[L, simplex_dims[1]] simplex[1] m1;
  array[L, simplex_dims[2]] simplex[2] m2;
  array[L, simplex_dims[3]] simplex[3] m3;
  array[L, simplex_dims[4]] simplex[4] m4;
  array[L, simplex_dims[5]] simplex[5] m5;
  array[L, simplex_dims[6]] simplex[6] m6;
  // Instantaneous stepwise rates
  array[T] vector<lower=0>[X] fishing_step;
  vector<lower=0>[X] natural_mortality_step;
  real<lower=0> ongoing_loss_step;
  // Fractional (per tag) stepwise rates
  vector<lower=0, upper=1>[X] reporting_step;
  real<lower=0, upper=1> initial_loss_step;
  // Selectivity (per fish)
  array[L - 1] vector<lower=0, upper=1>[X] selectivity_short;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Stepwise movement rate
  array[L] matrix<lower=0, upper=1>[X, X] movement_step;
  // Stepwise survival rate
  array[T, K, L] vector<lower=0, upper=1>[X] survival_step;
  // Stepwise transition rate
  array[N, L] matrix<lower=0, upper=1>[X, X] transition_step;
  // Stepwise observation rate
  array[N, L] vector<lower=0, upper=1>[X] observation_step;
  // Instantaneous annual rates
  array[T] vector<lower=0>[X] fishing_rate;
  vector<lower=0>[X] natural_mortality_rate = natural_mortality_step * K;
  real<lower=0> ongoing_loss_rate = ongoing_loss_step * K;
  // Fractional (per tag) rates
  vector<lower=0, upper=1>[X] reporting_rate = reporting_step;
  real<lower=0, upper=1> initial_loss_rate = initial_loss_step;
  // Selectivity
  array[L] vector<lower=0, upper=1>[X] selectivity;
  for (l in 1:L) {
    if (l == L) {
      selectivity[l] = rep_vector(1.0, X);
    } else {
      selectivity[l] = selectivity_short[l];
    }
  }
  //  // Fishing weight
  //  array[K] vector<lower=0, upper=1>[X] fishing_weight;
  // Assemble stepwise movement rates [L] [X, X]
  movement_step = assemble_movement_step(
    m1, m2, m3, m4, m5, m6,
    movement_index,
    L
  );
  //  // Assemble fishing weight [K][X]
  //  fishing_weight = assemble_fishing_weight(fishing_weight_transpose);
  // Assemble stepwise survival rate [T, K, L][X]
  survival_step = assemble_survival_step(
    fishing_step,
    //    fishing_weight,
    selectivity,
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
    selectivity,
    reporting_step,
    K, L
  );
  // Assemble fishing rate [T][X]
  fishing_rate = assemble_fishing_rate(fishing_step, K);
}

model {
  // Declare enumeration values
  array[N - 1, D, L] matrix[X, X] abundance;
  array[N - 1, D, L] matrix[X, X] predicted;
  array[C] int observed;
  array[C] real expected;
  // Initialize count
  int count = 0;
  // Populate released abundance
  for (n in 1:(N - 1)) { // Model step
    for (l in 1:L) { // Released size
      abundance[n, 1, l] = diag_matrix(
        tags_released[n, l] * (1 - initial_loss_step)
      );
    }
  }
  // Compute expected recoveries
  for (n in 1:(N - 1)) { // Model step
    for (d in 2:min(N - n + 1, D)) { // Duration at large
      for (l in 1:L) { // Released size
        // Propagate abundance
        abundance[n, d, l] = abundance[n, d - 1, l]
        * transition_step[n + d - 2, l]; // Previous step
        // Compute predicted
        predicted[n, d, l] = diag_post_multiply(
          abundance[n, d, l],
          observation_step[n + d - 1, l] // Current step
        );
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags_released[n, l, x] > 0) { // Were any tags released?
              if (movement_possible[d][x, y] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[n, d, l, y, x]; // Integer
                expected[count] = predicted[n, d, l, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End d
  } // End n
  // Movement step priors
  for (l in 1:L) {
    diagonal(movement_step[l]) ~ normal(
      mu_movement_step_diag,
      sd_movement_step_diag
    );
  }
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
  // Selectivity
  for (l in 1:L) {
    if (l < L) {
      selectivity_short[l] ~ normal(
        mu_selectivity_short[l],
        mu_selectivity_short[l] * cv_selectivity[1]
      );
    }
  }
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, sd_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);
}

generated quantities {
  // Annual movement rate
  array[L] matrix<lower=0, upper=1>[X, X] movement_rate;
  // Assemble movement rate [L][X, X]
  movement_rate = assemble_movement_rate(movement_step, K);
}
