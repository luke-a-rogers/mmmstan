
data {
  // Index limits
  int<lower=2> A; // Number of release and recovery areas
  int<lower=1> G; // Number of release groups
  int<lower=2> L; // Maximum number of time steps at liberty
  int<lower=1> T; // Number of release time steps
  int<lower=1> H; // Number of harvest rate time steps
  int<lower=2> I; // Number of study time steps (T + L - 1)
  // Constants
  int<lower=1> Y; // Number of time steps per year
  // Tag data
  int<lower=0> x[T, A, G, L, A]; // Tag array
  // Input rates
  real<lower=0, upper=1> u; // Initial tag retention rate (proportion)
  real<lower=0> v; // Tag loss rate
  real<lower=0> m; // Natural mortality rate
  // Index vectors
  int<lower=1> h_index[I];
  // Prior parameters
  real<lower=0> h_alpha[A];
  real<lower=0> h_beta[A];
  // Fudge constants
  real<lower=0> p_fudge;
  real<lower=0> y_fudge;
}

transformed data {
  // Transformed indexes
  real v_step = v / Y;
  real m_step = m / Y;
  int YT = 0; // TODO: Rename (Number of recovery observations)
  for (mt in 1:T) {
		for (ma in 1:A) {
			for (mg in 1:G) {
			  if (x[mt, ma, mg, 1, ma] > 0) {
			    YT += L * A;
			  }
			}
		}
  }
}

parameters {
  // Harvest rate
  real<lower=0, upper=1> h[H, A];
	// Length class 1
	simplex[2] p11;
	simplex[3] p12;
	simplex[2] p13;
	// // Length class 2
	// simplex[2] p21;
	// simplex[3] p22;
	// simplex[2] p23;
	// // Length class 3
	// simplex[2] p31;
	// simplex[3] p32;
	// simplex[2] p33;
}

transformed parameters {
  real f[H, A]; // Fishing mortality
  real f_step[H, A];
  for (ct in 1:H) {
    for (ca in 1:A) {
      f[ct, ca] = -log(1 - h[ct, ca]);
      f_step[ct, ca] = f[ct, ca] / Y;
    }
  }
}

// transformed parameters {
//   real p[G, A, A]; // [mg, ca, pa] Movement rates
// 	p = rep_array(1e-12, G, A, A);
//   // Length class 1
// 	p[1, 1, 1] = p11[1];
// 	p[1, 2, 1] = p11[2];
// 	p[1, 1, 2] = p12[1];
// 	p[1, 2, 2] = p12[2];
// 	p[1, 3, 2] = p12[3];
// 	p[1, 2, 3] = p13[1];
// 	p[1, 3, 3] = p13[2];
// 	// Length class 1
// 	p[2, 1, 1] = p21[1];
// 	p[2, 2, 1] = p21[2];
// 	p[2, 1, 2] = p22[1];
// 	p[2, 2, 2] = p22[2];
// 	p[2, 3, 2] = p22[3];
// 	p[2, 2, 3] = p23[1];
// 	p[2, 3, 3] = p23[2];
// 	// Length class 1
// 	p[3, 1, 1] = p31[1];
// 	p[3, 2, 1] = p31[2];
// 	p[3, 1, 2] = p32[1];
// 	p[3, 2, 2] = p32[2];
// 	p[3, 3, 2] = p32[3];
// 	p[3, 2, 3] = p33[1];
// 	p[3, 3, 3] = p33[2];
// }

model {
  real p_step[G, A, A]; // [mg, ca, pa] Movement rates
	// Initialize values
	real n[T, A, G, L, A]; // Predicted abundance
	real s_step[G, I, A]; // Survival rate
	int y_vec[YT]; // Recoveries
	real y_hat[YT]; // Predicted recoveries
	int y_ind; // Recovery index counter
	real n_sub[T, A];
	n = rep_array(rep_array(0, L, A), T, A, G);
	s_step = rep_array(0, G, I, A);
	y_vec = rep_array(0, YT);
	y_hat = rep_array(0, YT);
	y_ind = 1;

	// Slipped in here
	p_step = rep_array(p_fudge, G, A, A);
  // Length class 1
	p_step[1, 1, 1] = p11[1];
	p_step[1, 2, 1] = p11[2];
	p_step[1, 1, 2] = p12[1];
	p_step[1, 2, 2] = p12[2];
	p_step[1, 3, 2] = p12[3];
	p_step[1, 2, 3] = p13[1];
	p_step[1, 3, 3] = p13[2];

	// Compute survival
	for (mg in 1:G) {
		for (ct in 1:I) {
			for (ca in 1:A) {
				s_step[mg, ct, ca] = exp(-f_step[h_index[ct], ca] - m_step - v_step);
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
  				for (cl in 2:L) { // Populate abundance array n
	  				for (ca in 1:A) {
		  				for (pa in 1:A) {
			  				n[mt, ma, mg, cl, ca] += n[mt, ma, mg, cl - 1, pa]
		  					* s_step[mg, mt + cl - 2, pa] * p_step[mg, ca, pa];
	  					} // End for pa
  					} // End for ca
  				} // End for cl
  				for (cl in 2:L) { // Compute recoveries and predicted recoveries
  					for (ca in 1:A) {
  					  y_vec[y_ind] = x[mt, ma, mg, cl, ca];
  						y_hat[y_ind] = n[mt, ma, mg, cl, ca]
  						* (1 - exp(-f_step[h_index[mt + cl - 1], ca])) + y_fudge;
  						y_ind += 1;
  					} // End for ra
  				} // End for cl
			  } // End if
			} // End for mg
		} // End for ma
	} // End for mt

	// Priors
	for (ct in 1:H) {
    h[ct] ~ beta(h_alpha, h_beta);
	}

	// Sampling statement
	y_vec ~ poisson(y_hat);
}

generated quantities {
  real p_step[G, A, A]; // [mg, ca, pa] Movement rates
  // Annual movement rates
  matrix[A, A] p_step_matrix[G];
  matrix[A, A] p_matrix[G];
  real p[A, A, G];

  // Slipped in here
	p_step = rep_array(p_fudge, G, A, A);
  // Length class 1
	p_step[1, 1, 1] = p11[1];
	p_step[1, 2, 1] = p11[2];
	p_step[1, 1, 2] = p12[1];
	p_step[1, 2, 2] = p12[2];
	p_step[1, 3, 2] = p12[3];
	p_step[1, 2, 3] = p13[1];
	p_step[1, 3, 3] = p13[2];

  for (mg in 1:G) {
    // Populate p_matrix
    for (pa in 1:A) {
      for (ca in 1:A) {
        p_step_matrix[mg, pa, ca] = p_step[mg, ca, pa];
      }
    }
    // Populate p_matrix_annual
    p_matrix[mg] = matrix_power(p_step_matrix[mg], Y);
    // Populate p_annual
    for (pa in 1:A) {
      for (ca in 1:A) {
        p[pa, ca, mg] = p_matrix[mg, pa, ca];
      }
    }
  }
}
