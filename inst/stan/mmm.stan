
functions {
  #include create_simplex_dimensions.stan
  #include create_p_step.stan
}

data {
  // Index limits
  int<lower=2> A; // Number of release and recovery areas
  int<lower=1> G; // Number of release groups
  int<lower=2> L; // Maximum number of time steps at liberty
  int<lower=1> T; // Number of release time steps
  int<lower=1> H; // Number of harvest rate time steps
  int<lower=2> S; // Number of study time steps L + (T - 1) or L + 1
  int<lower=1> P; // Number of movement time steps
  int<lower=1> Q; // Number of harvest rate groups
  int<lower=1> W; // Number of tag reporting rate time steps
  // Constants
  int<lower=1> Y; // Number of time steps per year
  // Tag data
  int<lower=0> x[T, A, G, L, A]; // Tag array
  // Reporting rate
  real<lower=0, upper=1> w[W, A]; // Tag reporting rate array
  // Movement index array
  int<lower=0, upper=1> z[A, A]; // Movement index array
  // Input rates
  real<lower=0, upper=1> u; // Initial tag retention rate (proportion)
  real<lower=0> v; // Tag loss rate
  real<lower=0> m; // Natural mortality rate
  // Index vectors
  int<lower=1> h_index[S]; // Harvest rate time step index
  int<lower=1> p_index[S]; // Movement rate time step index
  int<lower=1> q_index[G]; // Harvest rate group index
  int<lower=1> w_index[S]; // Reporting rate time step index
  // Prior parameters
  real<lower=0> h_alpha[Q, H, A];
  real<lower=0> h_beta[Q, H, A];
  // Fudge constants
  real<lower=0> p_fudge;
  real<lower=0> y_fudge;
}

transformed data {
  // Initialize
  int simplex_dimensions[6] = create_simplex_dimensions(P, A, G, z);
  // Transformed indexes
  real v_step = v / Y;
  real m_step = m / Y;
  int R = 0; // Number of recovery observations
  for (mt in 1:T) {
    for (ma in 1:A) {
      for (mg in 1:G) {
        if (x[mt, ma, mg, 1, ma] > 0) {
          R += (min(L, S - mt) - 1) * A; // Realized liberty after release
        }
      }
    }
  }
}

parameters {
  // Harvest rate
  real<lower=0, upper=1> h[Q, H, A];
  // Movement simplexes
  simplex[1] s1[simplex_dimensions[1]]; // Not used
  simplex[2] s2[simplex_dimensions[2]];
  simplex[3] s3[simplex_dimensions[3]];
  simplex[4] s4[simplex_dimensions[4]];
  simplex[5] s5[simplex_dimensions[5]];
  simplex[6] s6[simplex_dimensions[6]];
  // Negative binomial dispersion var = mu + mu^2 / phi
  real<lower=0.1,upper=10> phi;
}

transformed parameters {
  // Initialize array version of movement simplexes
  real p1[simplex_dimensions[1], 1];
  real p2[simplex_dimensions[2], 2];
  real p3[simplex_dimensions[3], 3];
  real p4[simplex_dimensions[4], 4];
  real p5[simplex_dimensions[5], 5];
  real p6[simplex_dimensions[6], 6];
  // Initialize fishing mortality rates
  real f[Q, H, A];
  real f_step[Q, H, A];
  // Populate array version of movement simplexes
  for (i in 1:simplex_dimensions[1]) {for (j in 1:1) {p1[i, j] = s1[i, j]; } }
  for (i in 1:simplex_dimensions[2]) {for (j in 1:2) {p2[i, j] = s2[i, j]; } }
  for (i in 1:simplex_dimensions[3]) {for (j in 1:3) {p3[i, j] = s3[i, j]; } }
  for (i in 1:simplex_dimensions[4]) {for (j in 1:4) {p4[i, j] = s4[i, j]; } }
  for (i in 1:simplex_dimensions[5]) {for (j in 1:5) {p5[i, j] = s5[i, j]; } }
  for (i in 1:simplex_dimensions[6]) {for (j in 1:6) {p6[i, j] = s6[i, j]; } }
  // Populate fishing mortality rates
  for (cg in 1:Q) {
    for (ct in 1:H) {
      for (ca in 1:A) {
        f[cg, ct, ca] = -log(1 - h[cg, ct, ca]);
        f_step[cg, ct, ca] = f[cg, ct, ca] / Y;
      }
    }
  }
}

model {
	// Initialize values
	real p_step[G, P, A, A]; // [ , , ca, pa] Movement rates
  real n[T, A, G, L, A]; // Predicted abundance
	real s_step[G, S, A]; // Survival rate
	int y_obs[R]; // Recoveries
	real y_hat[R]; // Predicted recoveries
	int y_ind; // Recovery index counter
	real n_sub[T, A];
	int E = 0; // Number of realized steps at libery
	n = rep_array(rep_array(0, L, A), T, A, G);
	s_step = rep_array(0, G, S, A);
	y_obs = rep_array(0, R);
	y_hat = rep_array(0, R);
	y_ind = 1;

	// Create stepwise movement rates
	p_step = create_p_step(p_fudge, P, A, G, p1, p2, p3, p4, p5, p6, z);

	// Compute survival
	for (mg in 1:G) {
		for (ct in 1:S) {
			for (ca in 1:A) {
				s_step[mg, ct, ca] = exp(
				  -f_step[q_index[mg], h_index[ct], ca]
				  - m_step
				  - v_step);
			}
		}
	}

	// Populate initial abundances
	for (mt in 1:T) {
	  for (ma in 1:A) {
	    for (mg in 1:G) {
	      n[mt, ma, mg, 1, ma] = u * x[mt, ma, mg, 1, ma];
	    }
	  }
	}

	// Compute predicted recoveries
	for (mt in 1:T) {
		for (ma in 1:A) {
			for (mg in 1:G) {
			  if (x[mt, ma, mg, 1, ma] > 0) {
			    E = min(L, S - mt);
  				for (cl in 2:E) { // Populate abundance array n
	  				for (ca in 1:A) {
		  				for (pa in 1:A) {
			  				n[mt, ma, mg, cl, ca] += n[mt, ma, mg, cl - 1, pa]
		  					* s_step[mg, mt + cl - 2, pa]
		  					* p_step[mg, p_index[mt + cl - 2], ca, pa];
	  					} // End for pa
  					} // End for ca
  				} // End for cl
  				for (cl in 2:E) { // Compute recoveries and predicted recoveries
  					for (ca in 1:A) {
  					  y_obs[y_ind] = x[mt, ma, mg, cl, ca];
  						y_hat[y_ind] = n[mt, ma, mg, cl, ca]
  						* (1 - exp(-f_step[q_index[mg], h_index[mt + cl - 1], ca]))
  						* w[w_index[mt + cl - 1], ca]
  						+ y_fudge;
  						y_ind += 1;
  					} // End for ra
  				} // End for cl
			  } // End if
			} // End for mg
		} // End for ma
	} // End for mt

	// Priors
	for (cg in 1:Q) {
	  for (ct in 1:H) {
	    for (ca in 1:A) {
        h[cg, ct, ca] ~ beta(h_alpha[cg, ct, ca], h_beta[cg, ct, ca]);
	    }
	  }
	}

	// Sampling statement
	// y_obs ~ poisson(y_hat); // 11.9 sec simplest model
	y_obs ~ neg_binomial_2(y_hat, phi);
}

generated quantities {
  // Initialize
  real p_step[G, P, A, A]; // Stepwise model version [ , , ca, pa]
  matrix[A, A] p_step_matrix[G, P]; // Stepwise intermediary
  matrix[A, A] p_matrix[G, P]; // Annual intermediary
  real p[A, A, P, G]; // Annual user version [pa, ca, , ,]

  // Create stepwise movement rates
	p_step = create_p_step(p_fudge, P, A, G, p1, p2, p3, p4, p5, p6, z);

  // Populate annual movement rate array
  for (mg in 1:G) {
    for (ct in 1:P) {
      // Populate p_matrix
      for (pa in 1:A) {
        for (ca in 1:A) {
          p_step_matrix[mg, ct, pa, ca] = p_step[mg, ct, ca, pa];
        }
      }
      // Populate annual p_matrix
      p_matrix[mg, ct] = matrix_power(p_step_matrix[mg, ct], Y);
      // Populate annual array p
      for (pa in 1:A) {
        for (ca in 1:A) {
          p[pa, ca, ct, mg] = p_matrix[mg, ct, pa, ca];
        }
      }
    }
  }
}
