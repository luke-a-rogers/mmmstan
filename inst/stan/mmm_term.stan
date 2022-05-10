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
  matrix<lower=0, upper=1>[X, X] mu_movement_step_mean;
  real<lower=0> mu_cv_random_walk;
  real<lower=0> mu_dispersion;
  // Prior standard deviations
  matrix<lower=0, upper=1>[X, X] sd_movement_step_mean; // Not used
  real<lower=0> sd_cv_random_walk;
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
}

parameters {
  // Movement simplexes
  array[simplex_dimensions[1]] simplex[1] a1;
  array[simplex_dimensions[2]] simplex[2] a2;
  array[simplex_dimensions[3]] simplex[3] a3;
  array[simplex_dimensions[4]] simplex[4] a4;
  array[simplex_dimensions[5]] simplex[5] a5;
  array[simplex_dimensions[6]] simplex[6] a6;
  // Random walk coefficient of variation
  real<lower=0> cv_random_walk;
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
//  matrix<lower=0, upper=1>[X, X] movement_product; // [X, X]
  // Assemble movement step [X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6,
    movement_index,
    I, D // Each set to one
  );
  // Assemble movement product [X, X]
//  movement_product = assemble_movement_product(
//    movement_step,
//    K
//  );
  // Assemble survival step
  survival_step = rep_array(rep_vector(0.85, X), N, S);
  // Assemble observed step
  observed_step = rep_array(rep_vector(0.05, X), N, S);
}

model {
  // Declare enumeration values
  array[N - 1, S, L] matrix[X, X] abundance;
  array[N - 1, S, L] matrix[X, X] predicted;
  array[C] int observed;
  array[C] real expected;
  // Declare index values
  int count;
  // Populate released abundance
  for (n in 1:(N - 1)) { // Model step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        abundance[n, s, 1] = diag_matrix(
          tags_released[n, s] * (1 - initial_loss_step)
        );
      }
    }
  }
  // Initialize count
  count = 0;
  // Compute expected recoveries
  for (n in 1:(N - 1)) { // Partial sum index range within released step
    for (s in 1:S) { // Released size
      for (l in 2:min(N - n + 1, L)) { // Liberty step
        // Propagate abundance
        abundance[n, s, l] = abundance[n, s, l - 1]
        * diag_post_multiply(
            movement_step[n_to_i[n], s_to_d[s]],
            survival_step[n + l - 2, s]
          );
        // Compute predicted
        predicted[n, s, l] = diag_post_multiply(
          abundance[n, s, l],
          observed_step[n + l - 1, s]
        );
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags_released[n, s, x] > 0) { // Were any tags released?
              if (movement_possible[l][x, y] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[n, s, l, y, x]; // Integer
                expected[count] = predicted[n, s, l, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End s
  } // End n
  // Random walk coefficient of variation prior
  cv_random_walk ~ normal(mu_cv_random_walk, sd_cv_random_walk);
  // Movement step priors
  for (i in 2:I) {
    for (y in 1:X) {
      for (x in 1:X) {
        if (movement_index[x, y] == 1) {
          movement_step[i, 1][x, y] ~ normal(
            movement_step[i - 1, 1][x, y],
            mu_movement_step_mean[x, y] * cv_random_walk
          );
        }
      }
    }
  }
  for (y in 1:X) {
    for (x in 1:X) {
      if (movement_index[x, y] == 1) {
        movement_step[1, 1][x, y] ~ normal(
          movement_step[I, 1][x, y],
          mu_movement_step_mean[x, y] * cv_random_walk
        );
      }
    }
  }

//  for (i in 2:I) {
//    to_vector(movement_step[i, 1]) ~ normal(
//      to_vector(movement_step[i - 1, 1]),
//      to_vector(mu_movement_step) * cv_random_walk + 1e-12
//    );
//  }
//  to_vector(movement_step[1, 1]) ~ normal(
//    to_vector(movement_step[I, 1]),
//    to_vector(mu_movement_step) * cv_random_walk + 1e-12
//  );
//  to_vector(movement_product) ~ normal(
//    to_vector(mu_movement_mean),
//    to_vector(sd_movement_mean)
//  );
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, cv_dispersion * mu_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);
}

generated quantities {
  // Assemble movement mean
  array[I] matrix[X, X] movement_term = assemble_movement_term(
    movement_step,
    K
  );
}
