
functions {
  #include functions.stan
  #include parameters.stan
}

data {
  // Index limits
  int<lower=2> A; // Number of release and recovery areas
  int<lower=1> G_released; // Number of release groups
  int<lower=1> G_harvest; // Number of harvest rate groups
  int<lower=1> T_released; // Number of release time steps
  int<lower=2> T_liberty; // Maximum number of time steps at liberty
  int<lower=2> T_study; // T_liberty + (T_released - 1) or T_liberty + 1
  int<lower=1> T_movement; // Number of movement time steps
  int<lower=1> T_harvest; // Number of harvest rate time steps
  int<lower=1> T_reporting; // Number of tag reporting rate time steps
  // Constants
  int<lower=1> T_year; // Number of time steps per year
  // Input rates
  real<lower=0> m; // Natural mortality rate
  real<lower=0, upper=1> u; // Initial tag retention rate (proportion)
  real<lower=0> v; // Tag loss rate
  // Reporting rate
  real<lower=0, upper=1> w[T_reporting, A]; // Tag reporting rate array
  // Tag data
  int<lower=0> x[T_released, A, G_released]; // Released array
  int<lower=0> y[T_released, A, G_released, T_liberty, A]; // Recovered array
  // Movement index array
  int<lower=0, upper=1> z[A, A]; // Movement index array
  // Index vectors
  int<lower=1> p_time_index[T_study]; // Movement rate time step index
  int<lower=1> h_time_index[T_study]; // Harvest rate time step index
  int<lower=1> w_time_index[T_study]; // Reporting rate time step index
  int<lower=1> h_group_index[G_released]; // Harvest rate group index
  // Option constants
  int<lower=0, upper=1> random_walk; // Include random walk on retention rates
  // Prior parameters
  real<lower=0> h_prior_mean[G_harvest, T_harvest, A];
  real<lower=0> h_prior_sd[G_harvest, T_harvest, A];
  real<lower=0> phi_prior_mean;
  real<lower=0> phi_prior_sd;
  real<lower=0> sigma_prior_mean[(random_walk == 1 && T_movement > 1) ? A : 0];
  real<lower=0> sigma_prior_sd[(random_walk == 1 && T_movement > 1) ? A : 0];
  // Fudge constants
  real<lower=0> p_fudge;
  real<lower=0> y_fudge;
}

transformed data {
  // Transformed prior parameters
  real<lower=0> h_alpha[G_harvest, T_harvest, A];
  real<lower=0> h_beta[G_harvest, T_harvest, A];
  real<lower=0> phi_alpha;
  real<lower=0> phi_beta;
  real<lower=0> sigma_alpha[(random_walk == 1 && T_movement > 1) ? A : 0];
  real<lower=0> sigma_beta[(random_walk == 1 && T_movement > 1) ? A : 0];
  real<lower=0> pars[2] = rep_array(0.0, 2);
  // Simplex dimensions
  int simplex_dims[6] = simplex_dimensions(A, G_released, T_movement, z);
  // Transformed indexes
  real vs = v / T_year;
  real ms = m / T_year;
  int Y = 0; // Number of recovery observations
  for (rt in 1:T_released) {
    for (ra in 1:A) {
      for (rg in 1:G_released) {
        if (x[rt, ra, rg] > 0) {
          Y += (min(T_liberty, T_study - rt) - 1) * A; // Realized liberty
        }
      }
    }
  }
  // Populate harvest prior parameters
  for (hg in 1:G_harvest) {
    for (ht in 1:T_harvest) {
      for (ca in 1:A) {
        pars = beta_pars(h_prior_mean[hg, ht, ca], h_prior_sd[hg, ht, ca]);
        h_alpha[hg, ht, ca] = pars[1];
        h_beta[hg, ht, ca] = pars[2];
      }
    }
  }
  // Populate negative dispersion prior parameters
  pars = gamma_pars(phi_prior_mean, phi_prior_sd);
  phi_alpha = pars[1];
  phi_beta = pars[2];
  // Populate random walk standard deviation parameters
  if (random_walk == 1 && T_movement > 1) {
    for (ca in 1:A) {
      pars = gamma_pars(sigma_prior_mean[ca], sigma_prior_sd[ca]);
      sigma_alpha[ca] = pars[1];
      sigma_beta[ca] = pars[2];
    }
  }
}

parameters {
  // Harvest rate
  real<lower=0, upper=1> h[G_harvest, T_harvest, A];
  // Movement simplexes
  simplex[1] s1[simplex_dims[1]]; // Not used
  simplex[2] s2[simplex_dims[2]];
  simplex[3] s3[simplex_dims[3]];
  simplex[4] s4[simplex_dims[4]];
  simplex[5] s5[simplex_dims[5]];
  simplex[6] s6[simplex_dims[6]];
  // Random walk standard deviation
  real<lower=0> sigma[(random_walk == 1 && T_movement > 1) ? A : 0];
  // Negative binomial dispersion var = mu + mu^2 / phi
  real<lower=0> phi;
}

transformed parameters {
  // Initialize array version of movement simplexes
  real p1[simplex_dims[1], 1];
  real p2[simplex_dims[2], 2];
  real p3[simplex_dims[3], 3];
  real p4[simplex_dims[4], 4];
  real p5[simplex_dims[5], 5];
  real p6[simplex_dims[6], 6];
  // Initialize fishing mortality rates
  real f[G_harvest, T_harvest, A];
  real fs[G_harvest, T_harvest, A];
  // Populate array version of movement simplexes
  for (i in 1:simplex_dims[1]) {for (j in 1:1) {p1[i, j] = s1[i, j]; } }
  for (i in 1:simplex_dims[2]) {for (j in 1:2) {p2[i, j] = s2[i, j]; } }
  for (i in 1:simplex_dims[3]) {for (j in 1:3) {p3[i, j] = s3[i, j]; } }
  for (i in 1:simplex_dims[4]) {for (j in 1:4) {p4[i, j] = s4[i, j]; } }
  for (i in 1:simplex_dims[5]) {for (j in 1:5) {p5[i, j] = s5[i, j]; } }
  for (i in 1:simplex_dims[6]) {for (j in 1:6) {p6[i, j] = s6[i, j]; } }
  // Populate fishing mortality rates
  for (hg in 1:G_harvest) {
    for (ht in 1:T_harvest) {
      for (ca in 1:A) {
        f[hg, ht, ca] = -log(1 - h[hg, ht, ca]);
        fs[hg, ht, ca] = f[hg, ht, ca] / T_year;
      }
    }
  }
}

model {
  // Initialize values
  real ps[G_released, T_movement, A, A]; // [ , , ca, pa] Stepwise rates
  real n[T_released, A, G_released, T_liberty, A]; // Predicted abundance
  real ss[G_released, T_study, A]; // Stepwise survival rates
  int y_obs[Y]; // Recoveries
  real y_hat[Y]; // Predicted recoveries
  int y_ind; // Recovery index counter
  real n_sub[T_released, A];
  int E = 0; // Number of realized steps at libery
  n = rep_array(rep_array(0, T_liberty, A), T_released, A, G_released);
  ss = rep_array(0, G_released, T_study, A);
  y_obs = rep_array(0, Y);
  y_hat = rep_array(0, Y);
  y_ind = 1;

  // Create stepwise movement rates
  ps = p_step(A, G_released, T_movement, z, p1, p2, p3, p4, p5, p6, p_fudge);

  // Compute survival
  for (rg in 1:G_released) {
    for (ct in 1:T_study) {
      for (ca in 1:A) {
        ss[rg, ct, ca] = exp(
          -fs[h_group_index[rg], h_time_index[ct], ca]
          - ms
          - vs);
      }
    }
  }

  // Populate initial abundances
  for (rt in 1:T_released) {
    for (ra in 1:A) {
      for (rg in 1:G_released) {
        n[rt, ra, rg, 1, ra] = u * x[rt, ra, rg];
      }
    }
  }

  // Compute predicted recoveries
  for (rt in 1:T_released) {
    for (ra in 1:A) {
      for (rg in 1:G_released) {
        if (x[rt, ra, rg] > 0) {
          E = min(T_liberty, T_study - rt);
          for (lt in 2:E) { // Populate abundance array n
            for (ca in 1:A) {
              for (pa in 1:A) {
                n[rt, ra, rg, lt, ca] += n[rt, ra, rg, lt - 1, pa]
                * ss[rg, rt + lt - 2, pa]
                * ps[rg, p_time_index[rt + lt - 2], ca, pa];
              } // End for pa
            } // End for ca
          } // End for lt
          for (lt in 2:E) { // Compute recoveries and predicted recoveries
            for (ca in 1:A) {
              y_obs[y_ind] = y[rt, ra, rg, lt, ca];
              y_hat[y_ind] = n[rt, ra, rg, lt, ca]
              * (1 - exp(-fs[h_group_index[rg], h_time_index[rt + lt - 1], ca]))
              * w[w_time_index[rt + lt - 1], ca]
              + y_fudge;
              y_ind += 1;
            } // End for ca
          } // End for lt
        } // End if
      } // End for rg
    } // End for ra
  } // End for rt

  // Dispersion prior
  phi ~ gamma(phi_alpha, phi_beta);

  // Harvest rate priors
  for (hg in 1:G_harvest) {
    for (ht in 1:T_harvest) {
      for (ca in 1:A) {
        h[hg, ht, ca] ~ beta(h_alpha[hg, ht, ca], h_beta[hg, ht, ca]);
      }
    }
  }

  // Random walk priors
  if (random_walk == 1 && T_movement > 1) {
    sigma ~ gamma(sigma_alpha, sigma_beta);
  }

  // Random walk on retention rates (self-movement rates)
  if (random_walk == 1 && T_movement > 1) {
    for (rg in 1:G_released) {
      for (mt in 2:T_movement) {
        for (ca in 1:A) {
          ps[rg, mt, ca, ca] ~ normal(ps[rg, mt - 1, ca, ca], sigma[ca]);
        }
      }
    }
  }

  // Sampling statement
  // y_obs ~ poisson(y_hat); // 11.9 sec simplest model
  y_obs ~ neg_binomial_2(y_hat, phi);
}

generated quantities {
  // Initialize
  real ps[G_released, T_movement, A, A]; // Stepwise version [ , , ca, pa]
  matrix[A, A] ps_matrix[G_released, T_movement]; // Stepwise intermediary
  matrix[A, A] p_matrix[G_released, T_movement]; // Annual intermediary
  real p[A, A, T_movement, G_released]; // Annual user version [pa, ca, , ,]

  // Create stepwise movement rates
  ps = p_step(A, G_released, T_movement, z, p1, p2, p3, p4, p5, p6, p_fudge);

  // Populate annual movement rate array
  for (rg in 1:G_released) {
    for (mt in 1:T_movement) {
      // Populate p_matrix
      for (pa in 1:A) {
        for (ca in 1:A) {
          ps_matrix[rg, mt, pa, ca] = ps[rg, mt, ca, pa];
        }
      }
      // Populate annual p_matrix
      p_matrix[rg, mt] = matrix_power(ps_matrix[rg, mt], T_year);
      // Populate annual array p
      for (pa in 1:A) {
        for (ca in 1:A) {
          p[pa, ca, mt, rg] = p_matrix[rg, mt, pa, ca];
        }
      }
    }
  }
}
