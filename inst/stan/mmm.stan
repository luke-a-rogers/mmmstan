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
  // Movement mean priors
  vector<lower=0, upper=1>[X] mu_movement_mean_diag;
  vector<lower=0>[X] sd_movement_mean_diag;
  // Fishing rate priors
  array[T] vector<lower=0>[X] mu_fishing_rate;
  real<lower=0> cv_fishing_rate;
  // Selectivity priors
  //  vector<lower=0, upper=1>[L - 1] mu_selectivity;
  //  real<lower=0> cv_selectivity;
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

// Continue from here: array of movement rate matrices [L][X, X]









// Old below here --------------------------------------------------------------



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
//  int<lower=1> W; // Number of fishing weight steps (seasons per year)
  // Conversion index limits
  int<lower=1> J; // Factor to convert instantaneous stepwise rates to rates
  int<lower=1> K; // Matrix power to convert movement steps (to rates)
  int<lower=1> P; // Number of parameters in one [i,d][x,y] movement step slice
  // Movement index arrays
  array[N] int<lower=1, upper=I> n_to_i; // Model step to movement step index
  array[S] int<lower=1, upper=D> s_to_d; // Model size to movement size index
  // Fishing index arrays
  array[N] int<lower=1, upper=T> n_to_t; // Model step to fishing step index
//  array[N] int<lower=1, upper=W> n_to_w; // Model step to fishing weight step index
  // Tag data
  array[N - 1, S, L, X, X] int<lower=0> tags;
  // Movement index (will be paired with a matrix version)
  array[X, X] int<lower=0, upper=1> movement_index;
  // Movement step mean priors
  matrix<lower=0, upper=1>[X, X] mu_movement_step_mean;
  matrix<lower=0, upper=1>[X, X] sd_movement_step_mean;
  // Fishing rate priors
  array[T] vector<lower=0>[X] mu_fishing_rate;
  real<lower=0> cv_fishing_rate;
  // Selectivity priors
//  vector<lower=0, upper=1>[S - 1] mu_selectivity;
//  real<lower=0> cv_selectivity;
  // Fishing weight priors
//  array[W] vector<lower=0, upper=1>[X] mu_fishing_weight;
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
  real<lower=0, upper=1> mu_autoregress;
  real<lower=0, upper=1> sd_autoregress;
  real<lower=0> mu_sigma;
  real<lower=0> sd_sigma;
  // Dispersion priors
  real<lower=0> mu_dispersion;
  real<lower=0> sd_dispersion;
  // Tolerance values
  real<lower=0> tolerance_expected;
  real<lower=0> tolerance_fishing;
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
//  vector<lower=0, upper=1>[X] movement_step_diag;
  // Instantaneous stepwise rates
  array[T] vector<lower=0>[X] fishing_step;
  vector<lower=0>[X] mortality_step;
  real<lower=0> ongoing_loss_step;
  // Fractional (per tag) stepwise rates
  vector<lower=0, upper=1>[X] reporting_step;
  real<lower=0, upper=1> initial_loss_step;
  // Selectivity (per fish)
//  vector<lower=0, upper=1>[S - 1] selectivity_short;
  // Seasonal fishing weights
//  array[X] simplex[W] fishing_weight_transpose;
  // Autoregression priors
  vector<lower=0, upper=1>[form > 0 ? X : 0] autoregress;
  vector<lower=0>[form > 0 ? X : 0] sigma;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;
}

transformed parameters {
  // Declare stepwise rates
  array[I, D] matrix<lower=0, upper=1>[X, X] movement_step;
  array[N, S] vector<lower=0, upper=1>[X] survival_step;
  array[N, S] vector<lower=0, upper=1>[X] observed_step;
  // Declare stepwise movement deviations
  array[I, D] vector<lower=-1, upper=1>[X] movement_deviation;
  // Declare instantaneous rates
  array[T] vector<lower=0>[X] fishing_rate;
  vector<lower=0>[X] mortality_rate = mortality_step * J;
  real<lower=0> ongoing_loss_rate = ongoing_loss_step * J;
  // Declare fractional (per tag) rates
  vector<lower=0, upper=1>[X] reporting_rate = reporting_step;
  real<lower=0, upper=1> initial_loss_rate = initial_loss_step;
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
  if (form == 0) {
    movement_deviation = rep_array(rep_vector (0.0, X), I, D);
  } else {
    movement_deviation = assemble_movement_deviation(
      movement_step,
      mu_movement_step_mean
//      movement_step_mean
    );
  }
  // Assemble fishing rate
  fishing_rate = assemble_fishing_rate(fishing_step, J);
  // Assemble fishing weight
//  fishing_weight = assemble_fishing_weight(fishing_weight_transpose);
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
        * diag_pre_multiply(
            survival_step[n + l - 2, s],
            movement_step[n_to_i[n + l - 2], s_to_d[s]]
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
  // Movement step diag
//  movement_step_diag ~ normal(
//    diagonal(mu_movement_step_mean),
//    diagonal(sd_movement_step_mean)
//  );
  // Movement deviation prior
  if (form == 1) {
    // Autoregression priors
    autoregress ~ normal(mu_autoregress, sd_autoregress);
    sigma ~ normal(mu_sigma, sd_sigma);
    // Initial value
    for (d in 1:D) {
      movement_deviation[1, d] ~ normal(0.0, sigma);
    }
    // Autoregressive values
    for (i in 2:I) {
      for(d in 1:D) {
        movement_deviation[i, d] ~ normal(
         movement_deviation[i - 1, d] .* autoregress,
         sigma
        );
      }
    }
  } else if (form == 2) {
    // Autoregression priors
    autoregress ~ normal(mu_autoregress, sd_autoregress);
    sigma ~ normal(mu_sigma, sd_sigma);
    // Initial value
    for (d in 1:D) {
      movement_deviation[1, d] ~ normal(
        movement_deviation[I, d] .* autoregress,
        sigma
      );
    }
    // Autoregressive values
    for (i in 2:I) {
      for (d in 1:D) {
        movement_deviation[i, d] ~ normal(
         movement_deviation[i - 1, d] .* autoregress,
         sigma
        );
      }
    }
  } else if (form == 3) {
    // Autoregression priors
    autoregress ~ normal(mu_autoregress, sd_autoregress);
    sigma ~ normal(mu_sigma, sd_sigma);
    // Autoregressive values
    for (i in 1:I) {
      for (d in 2:D) {
        movement_deviation[i, d] ~ normal(
         movement_deviation[i, d - 1] .* autoregress,
         sigma
        );
      }
    }
  }
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
  // Declare movement mean
  matrix[X, X] movement_mean = rep_matrix(0.0, X, X);
  array[I] matrix[X, X] movement_time = rep_array(rep_matrix(0.0, X, X), I);
  array[I] matrix[X, X] movement_term = rep_array(rep_matrix(0.0, X, X), I);
  array[D] matrix[X, X] movement_size = rep_array(rep_matrix(0.0, X, X), D);
  // Populate movement mean
  if (form == 0) {
    movement_mean = matrix_power(movement_step[1, 1], K);
  }
  // Populate movement time
  if (form == 1) {
    for (i in 1:I) {
      movement_time[i] = matrix_power(movement_step[i, 1], K);
    }
  }
  // Populate movement term
  if (form == 2) {
    for (i in 1:I) {
      movement_term[i] = matrix_power(movement_step[i, 1], K);
    }
  }
  // Populate movement size
  if (form == 3) {
    for (d in 1:D) {
      movement_size[d] = matrix_power(movement_step[1, d], K);
    }
  }
}
