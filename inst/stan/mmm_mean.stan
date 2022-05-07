functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=2> N; // Number of model steps (months) in the study
  int<lower=1> S; // Number of size classes
  int<lower=2> L; // Number of maximum model steps at liberty
  int<lower=2> X; // Number of geographic regions
  int<lower=1> T; // Number of times (years)
  int<lower=1> I; // Number of terms (seasons) per time (year)
  int<lower=1> J; // Number of steps per time (year)
  int<lower=1> K; // Number of steps per term (season)
  int<lower=1> P; // Number of movement step parameters in one [N, S] slice
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
  array[8] int simplex_dimensions = assemble_simplex_dimensions(
    movement_index,
    1 // Movement step [1, X, X]
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
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Fix some parameters
  real<lower=0, upper=1> initial_loss_step = 0.1;
  // Declare movement rates
  matrix<lower=0, upper=1>[X, X] movement_step;
  // Assemble movement step [X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6, a7, a8,
    movement_index,
    1 // [1, X, X]
  )[1]; // [X, X]
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
        abundance[n, s, 1, x, x] = tags[n, s, 1, x, x]
        * (1 - initial_loss_step);
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
        * movement_step; // Square matrix [X, X]
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
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, cv_dispersion * mu_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);
}

generated quantities {
  // Assemble movement mean
  matrix[X, X] movement_mean = matrix_power(movement_step, J);
  // Assemble movement study
  matrix[X, X] movement_total = matrix_power(movement_step, N);
}
