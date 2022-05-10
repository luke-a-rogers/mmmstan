functions {
  #include functions.stan
}

data {
  // Model form
  int<lower=0, upper=3> form; // 0: mean; 1: time; 2: term; 3: size
  // Model index limits
  int<lower=2> N; // Number of model steps (months/quarters/years) in the study
  int<lower=1> S; // Number of size or sex classes
  int<lower=2> L; // Number of maximum model steps at liberty
  int<lower=2> X; // Number of geographic regions
  // Movement index limits
  int<lower=1> I; // Number of movement steps
  int<lower=1> D; // Number of movement sizes
  // Fishing index limits
  int<lower=1> T; // Number of fishing steps (almost always years)
  int<lower=1> W; // Number of fishing weight steps (seasons per year)
  // Conversion index limits
  int<lower=1> J; // Factor to convert instantaneous stepwise rates to rates
  int<lower=1> K; // Matrix power to convert movement steps (to rates)
  int<lower=1> P; // Number of parameters in one [i,d][x,y] movement step slice
  // Movement index arrays
  array[N] int<lower=1> n_to_i; // Model step to movement step index
  array[S] int<lower=1> s_to_d; // Model size to movement size index
  // Fishing index arrays
  array[N] int<lower=1> n_to_t; // Model step to fishing step index
  array[N] int<lower=1> n_to_w; // Model step to fishing weight step index
  // Tag data
  array[N - 1, S, L, X, X] int<lower=0> tags;
  // Movement index (will be paired with a matrix version)
  array[X, X] int<lower=0, upper=1> movement_index;
  // Movement step mean priors
  matrix<lower=0, upper=1>[X, X] mu_movement_step_mean;
  matrix<lower=0, upper=1>[X, X] sd_movement_step_mean;
  // Fishing rate priors
//  array[T] vector<lower=0>[X] mu_fishing_rate;
//  real<lower=0> cv_fishing_rate;
  // Selectivity priors
//  vector<lower=0, upper=1>[S - 1] mu_selectivity;
//  real<lower=0> cv_selectivity;
  // Fishing weight priors
//  array[W] vector<lower=0, upper=1> mu_fishing_weight;
//  real<lower=0> cv_fishing_weight;
  // Natural mortality rate priors
//  vector<lower=0>[X] mu_mortality_rate;
//  real<lower=0> cv_mortality_rate;
  // Fractional (per tag) reporting rate priors
//  vector<lower=0>[X] mu_reporting_rate;
//  real<lower=0> cv_reporting_rate;
  // Fractional (per tag) initial loss rate priors
  real<lower=0, upper=1> mu_initial_loss_rate;
  real<lower=0> cv_initial_loss_rate;
  // Instantaneous ongoing loss rate priors
//  real<lower=0> mu_ongoing_loss_rate;
//  real<lower=0> cv_ongoing_loss_rate;
  // Dispersion priors
  real<lower=0> mu_dispersion;
  real<lower=0> sd_dispersion;
  // Tolerance values
  real<lower=0> tolerance_expected;
//  real<lower=0> tolerance_fishing;
}

transformed data {
  // Declare upper bound on number of observations
  int<lower=0> C = N * S * L * X * X;
  // Declare movement matrix version of movement index
  matrix[X, X] movement_matrix = assemble_movement_matrix(movement_index);
  // Declare simplex dimensions
  array[6] int simplex_dimensions = assemble_simplex_dimensions(
    movement_index,
    I, D
  );
  // Declare tags released
  array[N - 1, S] vector[X] tags_released = assemble_tags_released(tags);
  // Declare tags transpose (permute x and y)
  array[N - 1, S, L, X, X] int tags_transpose = assemble_tags_transpose(tags);
  // Declare movement possible values
  array[L] matrix[X, X] movement_possible = assemble_movement_possible(
    movement_matrix,
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
  // Movement mean stepwise rate
//  matrix<lower=0, upper=1>[X, X] movement_step_mean;
  // Instantaneous stepwise rates
//  array[T] vector<lower=0>[X] fishing_step;
//  vector<lower=0>[X] mortality_step;
//  real<lower=0> ongoing_loss_step;
  // Fractional (per tag) stepwise rates
//  vector<lower=0, upper=1>[X] reporting_step;
  real<lower=0, upper=1> initial_loss_step;
  // Selectivity (per fish)
//  vector<lower=0, upper=1>[S - 1] selectivity_short;
  // Seasonal fishing weights
//  array[X] simplex[W] fishing_weight_transpose;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Declare stepwise rates
  array[I, D] matrix<lower=0, upper=1>[X, X] movement_step;
  array[N, S] vector<lower=0, upper=1>[X] survival_step;
  array[N, S] vector<lower=0, upper=1>[X] observed_step;
  // Declare stepwise movement deviations
//  array[I, D] matrix<lower=-1, upper=1> movement_deviation;
  // Declare instantaneous rates
//  array[T] vector<lower=0>[X] fishing_rate;
//  vector<lower=0>[X] mortality_rate = mortality_step * J;
//  real<lower=0> ongoing_loss_rate = ongoing_loss_step * J;
  // Declare fractional (per tag) rates
//  vector<lower=0, upper=1>[X] reporting_rate = reporting_step;
//  real<lower=0, upper=1> initial_loss_rate = initial_loss_step;
  // Selectivity
//  vector<lower=0, upper=1>[S] selectivity = append_row(selectivity_short, 1.0);
  // Declare fishing weight
//  array[W] vector<lower=0, upper=1> fishing_weight;
  // Assemble movement step [I, D][X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6,
    movement_index,
    I, D
  );
  // Assemble movement deviation
//  movement_deviation = assemble_movement_deviation(
//    movement_step,
//    movement_step_mean
//  );
  // Assemble fishing rate
//  fishing_rate = assemble_fishing_rate(fishing_step, J);
  // Assemble fishing weight
//  fishing_weight = assemble_fishing_weight(fishing_weight_transpose);
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
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, sd_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);
}

generated quantities {
  // Declare movement mean
  matrix[form == 0 ? X : 0, form == 0 ? X : 0] movement_mean;
  // Populate movement mean
  if (form == 0) {
    movement_mean = matrix_power(movement_step[1, 1], K);
  }
}
