functions {
  #include functions.stan
}

data {
  // Index limits
  int<lower=2> T; // Number of model steps (months/quarters/years) in the study
  int<lower=2> D; // Maximum duration at large in model steps
  int<lower=1> L; // Number of size or sex classes
  int<lower=2> S; // Number of geographic regions
  int<lower=1> Years; // Number of years in the study
  // Constants
  int<lower=1> K; // Number of steps (months/quarters/years) per year
  // Movement index arrays
  // array[T] int<lower=1, upper=Years> n_to_t; // Model step to year index
  // Tag data
  array[T - 1, D, L, S, S] int<lower=0> tags;
  // Movement index (will be paired with a matrix version)
  array[S, S] int<lower=0, upper=1> movement_index;
  // Movement step priors
  vector<lower=0, upper=1>[S] mu_movement_step_diag;
  vector<lower=0>[S] sd_movement_step_diag;
  // Fishing rate priors
  array[Years] vector<lower=0>[S] mu_fishing_rate;
  real<lower=0> cv_fishing_rate;
  // Selectivity priors
  array[L - 1] vector<lower=0, upper=1>[S] mu_selectivity_short;
  vector<lower=0>[L > 1 ? 1 : 0] cv_selectivity;
  // Fishing weight priors
  //  array[K] vector<lower=0, upper=1>[S] mu_fishing_weight;
  //  real<lower=0> cv_fishing_weight;
  // Natural mortality rate priors
  vector<lower=0>[S] mu_natural_mortality_rate;
  vector<lower=0>[S] sd_natural_mortality_rate;
  // Fractional (per tag) reporting rate priors
  vector<lower=0, upper=1>[S] mu_reporting_rate;
  vector<lower=0>[S] sd_reporting_rate;
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
  int<lower=0> C = T * D * L * S * S;
  // Matrix version of movement index
  matrix[S, S] movement_matrix = to_matrix(movement_index);
  // Simplex dimensions
  array[6] int simplex_dims = assemble_simplex_dims(movement_index);
  // Declare tags released
  array[T - 1, L] vector[S] tags_released = assemble_tags_released(tags);
  // Declare tags transpose (permute s0 and s)
  array[T - 1, D, L, S, S] int tags_transpose = assemble_tags_transpose(tags);
  // Declare movement possible values
  array[D] matrix[S, S] movement_possible = assemble_movement_possible(
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
  // Instantaneous yearly rates
  array[Years] vector<lower=0>[S] fishing_rate;
  vector<lower=0>[S] natural_mortality_rate;
  real<lower=0> ongoing_loss_rate;
  // Fractional (per tag) rates
  vector<lower=0, upper=1>[S] reporting_rate;
  real<lower=0, upper=1> initial_loss_rate;
  // Selectivity (per fish)
  array[L - 1] vector<lower=0, upper=1>[S] selectivity_short;
  // Negative binomial dispersion parameter
  real<lower=0> dispersion;

  // INFO: Updated above 2024-08-23
  // Instantaneous stepwise rates
  // array[Years] vector<lower=0>[S] fishing_step;
  // vector<lower=0>[S] natural_mortality_step;
  // real<lower=0> ongoing_loss_step;

  // Fractional (per tag) stepwise rates
  // vector<lower=0, upper=1>[S] reporting_step;
  // real<lower=0, upper=1> initial_loss_step;

  // Selectivity (per fish)
  // array[L - 1] vector<lower=0, upper=1>[S] selectivity_short;
  // Negative binomial dispersion parameter
  // real<lower=0> dispersion;
}

transformed parameters {
  // Stepwise movement rate
  array[L] matrix<lower=0, upper=1>[S, S] movement_step;
  // Selectivity per fish
  array[L] vector<lower=0, upper=1>[S] selectivity;
  // Stepwise selected weighted fishing rate
  array[T, L] vector<lower=0>[S] selected_weighted_fishing_step;
  // Stepwise natural mortality plus loss rate
  vector<lower=0>[S] natural_mortality_plus_loss_step;
  // Stepwise survive then move rate
  array[T, L] matrix<lower=0, upper=1>[S, S] survive_then_move_step;
  // Stepwise get caught and reported rate
  array[T, L] vector<lower=0, upper=1>[S] get_caught_and_reported_step;

  // INFO: Updated above 2024-08-23
  // Stepwise survival rate
  // array[T, L] vector<lower=0, upper=1>[S] survival_step;

  // INFO: Updated above 2024-08-23
  // Stepwise transition rate
  // array[T, L] matrix<lower=0, upper=1>[S, S] transition_step;

  // // Stepwise observation rate
  // array[T, L] vector<lower=0, upper=1>[S] observation_step;

  //  // Fishing weight
  //  array[K] vector<lower=0, upper=1>[S] fishing_weight;

  // INFO: Updated above 2024-08-23
  // Instantaneous annual rates
  // array[Years] vector<lower=0>[S] fishing_rate;
  // vector<lower=0>[S] natural_mortality_rate = natural_mortality_step * K;
  // real<lower=0> ongoing_loss_rate = ongoing_loss_step * K;
  // Fractional (per tag) rates
  // vector<lower=0, upper=1>[S] reporting_rate = reporting_step;
  // real<lower=0, upper=1> initial_loss_rate = initial_loss_step;
  // Selectivity
  // array[L] vector<lower=0, upper=1>[S] selectivity;
  // for (l in 1:L) {
  //   if (l == L) {
  //     selectivity[l] = rep_vector(1.0, S);
  //   } else {
  //     selectivity[l] = selectivity_short[l];
  //   }
  // }

  // Assemble stepwise movement rates [L][S, S]
  movement_step = assemble_movement_step(
    m1, m2, m3, m4, m5, m6,
    movement_index,
    L
  );
  //  // Assemble fishing weight [K][S]
  //  fishing_weight = assemble_fishing_weight(fishing_weight_transpose);

  // Assemble selectivity [L][S]
  selectivity = assemble_selectivity(selectivity_short, L, S);

  // Assemble selected weighted fishing step [T, L][S]
  selected_weighted_fishing_step = assemble_selected_weighted_fishing_step(
    fishing_rate,
    // fishing_weight,
    selectivity,
    T, K
  );

  // Assemble natural mortality plus loss step [S]
  natural_mortality_plus_loss_step = (1/(1.0 * K))
  * (natural_mortality_rate + ongoing_loss_rate);

  // Assemble survive then move step [T, L][S, S]
  survive_then_move_step = assemble_survive_then_move_step(
    movement_step,
    selected_weighted_fishing_step,
    natural_mortality_plus_loss_step
  );

  // Assemble get caught and reported step [T, L][S]
  get_caught_and_reported_step = assemble_get_caught_and_reported_step(
    selected_weighted_fishing_step,
    reporting_rate
  );

  // INFO: Updated above 2024-08-23
  // // Assemble stepwise survive then move [T, L][S, S]
  // survive_then_move_step = assemble_survive_then_move_step(
  //   movement_step,
  //   fishing_rate,
  //   // fishing_weight,
  //   selectivity,
  //   natural_mortality_rate,
  //   ongoing_loss_rate,
  //   T, K
  // );

  // INFO: Updated above 2024-08-23
  // // Assemble stepwise survival rate [T, L][S]
  // survival_step = assemble_survival_step(
  //   fishing_rate,
  //   // fishing_weight,
  //   selectivity,
  //   natural_mortality_rate,
  //   ongoing_loss_rate,
  //   T, K
  // );

  // INFO: Updated above 2024-08-23
  // // Assemble stepwise survival rate [Years, K, L][S]
  // survival_step = assemble_survival_step(
  //   fishing_step,
  //   //    fishing_weight,
  //   selectivity,
  //   natural_mortality_step,
  //   ongoing_loss_step,
  //   K, L
  // );
  // // Assemble stepwise transition rate [T, L][S, S]
  // transition_step = assemble_transition_step(
  //   movement_step,
  //   survival_step
  // );

  // INFO: Updated above 2024-08-23
  // // Assemble stepwize observation rate [T, L][S]
  // observation_step = assemble_observation_step(
  //   fishing_rate,
  //   // fishing_weight,
  //   selectivity,
  //   reporting_rate,
  //   T, K
  // );

  // INFO: Updated above 2024-08-23
  // // Assemble fishing rate [Years][S]
  // fishing_rate = assemble_fishing_rate(fishing_step, K);
}

model {
  // Declare enumeration values
  array[T - 1, D, L] matrix[S, S] abundance;
  array[T - 1, D, L] matrix[S, S] predicted;
  array[C] int observed;
  array[C] real expected;
  // Initialize count
  int count = 0;
  // Populate released abundance
  for (t in 1:(T - 1)) { // Model step
    for (l in 1:L) { // Released size
      abundance[t, 1, l] = diag_matrix(
        tags_released[t, l] * (1 - initial_loss_rate)
      );
    }
  }
  // Compute expected recoveries
  for (t in 1:(T - 1)) { // Model step
    for (d in 2:min(T - t + 1, D)) { // Duration at large
      for (l in 1:L) { // Released size
        // Propagate abundance
        abundance[t, d, l] = abundance[t, d - 1, l]
        // * diag_pre_multiply(survival_step[t + d - 2, l], movement_step[l]);
        * survive_then_move_step[t + d - 2, l];

        // INFO: Updated above 2024-08-03
        // abundance[t, d, l] = abundance[t, d - 1, l]
        // * transition_step[t + d - 2, l]; // Previous step

        // Compute predicted
        predicted[t, d, l] = diag_post_multiply(
          abundance[t, d, l],
          // observation_step[t + d - 1, l] // Current step
          get_caught_and_reported_step[t + d - 1, l]
        );


        // Compute vectors
        for (s in 1:S) { // Current region
          for (s0 in 1:S) { // Released region
            if (tags_released[t, l, s0] > 0) { // Were any tags released?
              if (movement_possible[d][s0, s] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[t, d, l, s, s0]; // Integer
                expected[count] = predicted[t, d, l, s0, s]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End s0
        } // End s
      } // End l
    } // End d
  } // End t
  // Movement step priors
  for (l in 1:L) {
    diagonal(movement_step[l]) ~ normal(
      mu_movement_step_diag,
      sd_movement_step_diag
    );
  }
  // Fishing rate prior
  for (year in 1:Years) {
    fishing_rate[year] ~ normal(
      mu_fishing_rate[year],
      mu_fishing_rate[year] * cv_fishing_rate + tolerance_fishing
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
  array[L] matrix<lower=0, upper=1>[S, S] movement_rate;
  // Assemble movement rate [L][S, S]
  movement_rate = assemble_movement_rate(movement_step, K);
}
