functions {
  #include functions.stan
}

data {
  // Model form
  int<lower=0, upper=1> model_time;
  int<lower=0, upper=1> model_term;
  int<lower=0, upper=1> model_size;
  // Index limits
  int<lower=2> N; // Number of model steps (months/quarters/years) in the study
  int<lower=2> D; // Maximum duration at large in model steps
  int<lower=1> T; // Number of years in the study
  int<lower=1> K; // Number of terms (quarters/trimesters/halfs) per year
  int<lower=2> L; // Number of size or sex classes
  int<lower=2> X; // Number of geographic regions
  // Constants
  int<lower=1> H; // Number of steps per term
  int<lower=1> J; // Number of steps per year
  // Movement index arrays
  array[N] int<lower=1, upper=T> n_to_t; // Model step to year index
  array[N] int<lower=1, upper=K> n_to_k; // Model step to term index
  // Tag data
  array[N - 1, L, D, X, X] int<lower=0> tags;
  // Movement index (will be paired with a matrix version)
  array[X, X] int<lower=0, upper=1> movement_index;
  // Fishing rate priors
  array[T] vector<lower=0>[X] mu_fishing_rate;
  real<lower=0> cv_fishing_rate;
  // Selectivity priors
//  vector<lower=0, upper=1>[L - 1] mu_selectivity;
//  real<lower=0> cv_selectivity;
  // Fishing weight priors
//  array[K] vector<lower=0, upper=1>[X] mu_fishing_weight;
//  real<lower=0> cv_fishing_weight;
  // Natural mortality rate priors
  vector<lower=0>[X] mu_mortality_rate;
  vector<lower=0>[X] sd_mortality_rate;
  // Fractional (per tag) reporting rate priors
  vector<lower=0, upper=1>[X] mu_reporting_rate;
  vector<lower=0>[X] sd_reporting_rate;
  // Fractional (per tag) initial loss rate priors
  real<lower=0, upper=1> mu_initial_loss_rate;
  real<lower=0> sd_initial_loss_rate;
  // Instantaneous ongoing loss rate priors
  real<lower=0> mu_ongoing_loss_rate;
  real<lower=0> sd_ongoing_loss_rate;
  // Autoregression priors
//  real<lower=0, upper=1> mu_autoregress;
//  real<lower=0, upper=1> sd_autoregress;
//  real<lower=0> mu_sigma;
//  real<lower=0> sd_sigma;
  // Dispersion priors
  real<lower=0> mu_dispersion;
  real<lower=0> sd_dispersion;
  // Tolerance values
  real<lower=0> tolerance_expected;
  real<lower=0> tolerance_fishing;
}

transformed data {
  // Maximum number of observations
  int<lower=0> C = N * L * D * X * X;
  // Matrix version of movement index
  matrix[X, X] movement_matrix = assemble_movement_matrix(movement_index);
  // Simplex dimensions
  array[6] int simplex_dimensions_mean = assemble_simplex_dimensions(
    movement_index,
    1, 1
  );
  array[6] int simplex_dimensions_time = assemble_simplex_dimensions(
    movement_index,
    model_time,
    T
  );
  array[6] int simplex_dimensions_term = assemble_simplex_dimensions(
    movement_index,
    model_term,
    K
  );
  array[6] int simplex_dimensions_size = assemble_simplex_dimensions(
    movement_index,
    model_size,
    L
  );
  // Declare tags released
  array[N - 1, L] vector[X] tags_released = assemble_tags_released(tags);
  // Declare tags transpose (permute x and y)
  array[N - 1, L, D, X, X] int tags_transpose = assemble_tags_transpose(tags);
  // Declare movement possible values
  array[D] matrix[X, X] movement_possible = assemble_movement_possible(
    movement_matrix,
    D
  );
}

parameters {
  // Stepwise movement mean simplexes
  array[simplex_dimensions_mean[1]] simplex[1] m1;
  array[simplex_dimensions_mean[2]] simplex[2] m2;
  array[simplex_dimensions_mean[3]] simplex[3] m3;
  array[simplex_dimensions_mean[4]] simplex[4] m4;
  array[simplex_dimensions_mean[5]] simplex[5] m5;
  array[simplex_dimensions_mean[6]] simplex[6] m6;
  // Stepwise movement time simplexes
  array[simplex_dimensions_time[1]] simplex[1] t1;
  array[simplex_dimensions_time[2]] simplex[2] t2;
  array[simplex_dimensions_time[3]] simplex[3] t3;
  array[simplex_dimensions_time[4]] simplex[4] t4;
  array[simplex_dimensions_time[5]] simplex[5] t5;
  array[simplex_dimensions_time[6]] simplex[6] t6;
  // Stepwise movement term simplexes
  array[simplex_dimensions_term[1]] simplex[1] k1;
  array[simplex_dimensions_term[2]] simplex[2] k2;
  array[simplex_dimensions_term[3]] simplex[3] k3;
  array[simplex_dimensions_term[4]] simplex[4] k4;
  array[simplex_dimensions_term[5]] simplex[5] k5;
  array[simplex_dimensions_term[6]] simplex[6] k6;
  // Stepwise movement size simplexes
  array[simplex_dimensions_size[1]] simplex[1] l1;
  array[simplex_dimensions_size[2]] simplex[2] l2;
  array[simplex_dimensions_size[3]] simplex[3] l3;
  array[simplex_dimensions_size[4]] simplex[4] l4;
  array[simplex_dimensions_size[5]] simplex[5] l5;
  array[simplex_dimensions_size[6]] simplex[6] l6;

  // Instantaneous stepwise rates
  array[T] vector<lower=0>[X] fishing_step;
  vector<lower=0>[X] mortality_step;
  real<lower=0> ongoing_loss_step;
  // Fractional (per tag) stepwise rates
  vector<lower=0, upper=1>[X] reporting_step;
  real<lower=0, upper=1> initial_loss_step;
//  // Selectivity (per fish)
//  vector<lower=0, upper=1>[L - 1] selectivity_short;
//  // Seasonal fishing weights
//  array[X] simplex[K] fishing_weight_transpose;
//  // Autoregression priors
//  vector<lower=0, upper=1>[form > 0 ? X : 0] autoregress;
//  vector<lower=0>[form > 0 ? X : 0] sigma;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Stepwise movement rate components
  martix<lower=0, upper=1>[X, X] movement_step_mean;
  array[T] matrix<lower=0, upper=1>[X, X] movement_tran_time;
  array[K] matrix<lower=0, upper=1>[X, X] movement_tran_term;
  array[L] matrix<lower=0, upper=1>[X, X] movement_tran_size;
  // Stepwise movement rate
  array[T, K, L] matrix<lower=0, upper=1>[X, X] movement_step;
  // Stepwise survival rate
  array[T, K, L] vector<lower=0, upper=1>[X] survival_step;
  // Stepwise transition rate
  array[N, L] matrix<lower=0, upper=1>[X, X] transition_step;
  // Stepwise observation rate
  array[N, L] vector<lower=0, upper=1>[X] observation_step;

  // Declare instantaneous rates
  array[T] vector<lower=0>[X] fishing_rate;
  vector<lower=0>[X] mortality_rate = mortality_step * J;
  real<lower=0> ongoing_loss_rate = ongoing_loss_step * J;
  // Declare fractional (per tag) rates
  vector<lower=0, upper=1>[X] reporting_rate = reporting_step;
  real<lower=0, upper=1> initial_loss_rate = initial_loss_step;
//  // Selectivity
//  vector<lower=0, upper=1>[S] selectivity = append_row(selectivity_short, 1.0);
//  // Declare fishing weight
//  array[K] vector<lower=0, upper=1> fishing_weight;

// TODO:
  // Assemble stepwise movement rate components
  movement_step_mean = assemble_movement_step_mean();
  movement_tran_time = assemble_movement_transformation();
  movement_tran_term = assemble_movement_transformation();
  movement_tran_size = assemble_movement_transformation();
  // Assemble stepwise movement rate [T, K, L][X, X]
  movement_step = assemble_movement_step();
  // Assemble stepwise survival rate [T, K, L][X]
  survival_step = assemble_survival_step();
  // Assemble stepwise transition rate [N, L][X, X]
  transition_step = assemble_transition_step();
  // Assemble stepwize observation rate [N, L][X]
  observation_step = assemble_observation_step();

  // Assemble fishing rate [T][X]
  fishing_rate = assemble_fishing_rate(fishing_step, J);
//  // Assemble fishing weight [K][X]
//  fishing_weight = assemble_fishing_weight(fishing_weight_transpose);

// TODO: Remove
  // Assemble movement step [I, D][X, X]
  movement_step = assemble_movement_step(
    a1, a2, a3, a4, a5, a6,
    movement_index,
    I, D
  );
  // Assemble survival step
  survival_step = assemble_survival_step(
    fishing_step,
    mortality_step,
    ongoing_loss_step,
    n_to_t,
    S
  );
  // Assemble observed step
  observed_step = assemble_observed_step(
    fishing_step,
    reporting_step,
    n_to_t,
    S
  );
}

model {
  // Declare enumeration values
  array[N - 1, L, D] matrix[X, X] abundance;
  array[N - 1, L, D] matrix[X, X] predicted;
  array[C] int observed;
  array[C] real expected;
  // Declare index values
  int count;
  // Populate released abundance
  for (n in 1:(N - 1)) { // Model step
    for (l in 1:L) { // Released size
      abundance[n, l, 1] = diag_matrix(
        tags_released[n, l] * (1 - initial_loss_step)
      );
    }
  }
  // Initialize count
  count = 0;
  // Compute expected recoveries
  for (n in 1:(N - 1)) { // Model step
    for (l in 1:L) { // Released size
      for (d in 2:min(N - n + 1, D)) { // Duration at large
        // Propagate abundance
        abundance[n, l, d] = abundance[n, l, d - 1]
        * transition_step[n + d - 2, l];

//        * diag_post_multiply(
//            movement_step[n_to_i[n], s_to_d[s]], // apparent index error
//            survival_step[n + l - 2, s]
//          );

        // Compute predicted
        predicted[n, l, d] = diag_post_multiply(
          abundance[n, l, d],
          observation_step[n + d - 1, l]
        );
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags_released[n, l, x] > 0) { // Were any tags released?
              if (movement_possible[d][x, y] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[n, l, d, y, x]; // Integer
                expected[count] = predicted[n, l, d, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End s
  } // End n
  // Stepwise movement rate priors

  // Fishing rate prior
  for (t in 1:T) {
    fishing_rate[t] ~ normal(
      mu_fishing_rate[t],
      mu_fishing_rate[t] * cv_fishing_rate + tolerance_fishing
    );
  }
  // Mortality rate prior
  mortality_rate ~ normal(mu_mortality_rate, sd_mortality_rate);
  // Reporting rate prior
  reporting_rate ~ normal(mu_reporting_rate, sd_reporting_rate);
  // Ongoing loss rate prior
  ongoing_loss_rate ~ normal(mu_ongoing_loss_rate, sd_ongoing_loss_rate);
  // Initial loss rate prior
  initial_loss_rate ~ normal(mu_initial_loss_rate, sd_initial_loss_rate);
  // Dispersion prior
  dispersion ~ normal(mu_dispersion, sd_dispersion);
  // Sampling statement (var = mu + mu^2 / dispersion)
  observed[1:count] ~ neg_binomial_2(expected[1:count], dispersion);
}

generated quantities {
  // Movement rate components
  matrix<lower=0, upper=1>[X, X] movement_mean;
  array[T] matrix<lower=0, upper=1>[X, X] movement_time;
  array[K] matrix<lower=0, upper=1>[X, X] movement_term;
  array[L] matrix<lower=0, upper=1>[X, X] movement_size;
  // Movement rate
  array[T, K, L] matrix<lower=0, upper=1>[X, X] movement_rate;

  // Populate movement mean
  // Populate movement time
  // Populate movement term
  // Populate movement size
  // Populate movement rate
}
