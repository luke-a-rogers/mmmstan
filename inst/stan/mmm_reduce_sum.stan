
functions {
  // Partial sum function for within-chain parallel threading
  real partial_sum_lpmf(
    int[] release_steps_slice,
    int start,
    int end,
    int T,
    int A,
    int G,
    int L,
    real u,
    real[] f,
    real[,,] p,
    real[,,] s,
    int[,,,,] x
  ) {
    // Instantiate objects
    int span = end - start + 1;
    int Y = span * A * G * L * A;
    int y = 1;
    int y_obs[Y];
    real y_hat[Y];
    real n[span, A, G, L, A]; // Predicted abundance

    // Initialize objects
    n = rep_array(rep_array(0, L, A), span, A, G);
    y_obs = rep_array(0, Y);
    y_hat = rep_array(0, Y);

  	// Populate initial abundances
  	for (mt in start:end) {
  	  for (ma in 1:A) {
  	    for (mg in 1:G) {
  	      n[mt - start + 1, ma, mg, 1, ma] = u * x[mt, ma, mg, 1, ma];
  	    }
  	  }
  	}
	  // Compute predicted recoveries
	  for (mt in start:end) {
	  	for (ma in 1:A) {
		  	for (mg in 1:G) {
  			  for (cl in 2:L) { // Populate abundance array n
	  				for (ca in 1:A) {
		  		  	for (pa in 1:A) {
			  		  	n[mt - start + 1, ma, mg, cl, ca]
			  		  	+= n[mt - start + 1, ma, mg, cl - 1, pa]
		  				  * s[mg, mt + cl - 2, pa] * p[mg, ca, pa];
	  				  } // End for pa
  				  } // End for ca
  			  } // End for cl
  			  for (cl in 2:L) { // Compute recoveries and predicted recoveries
  				  for (ca in 1:A) {
  				    y_obs[y] = x[mt, ma, mg, cl, ca];
  					  y_hat[y] = n[mt - start + 1, ma, mg, cl, ca]
  					  * (1 - exp(-f[ca])) + 1e-12;
  					  y += 1;
  				  } // End for ra
  			  } // End for cl
			  } // End for mg
		  } // End for ma
	  } // End for mt
	  // Return partial sum
    return poisson_lupmf(y_obs | y_hat);
  }
}

data {
  // Index limits
  int<lower=2> A; // Number of release and recovery areas
  int<lower=1> G; // Number of release groups
  int<lower=2> L; // Maximum number of time steps at liberty
  int<lower=1> T; // Number of release time steps
  // Constants
  int<lower=1> Y; // Number of time steps per year
  // Tag data
  int<lower=0> x[T, A, G, L, A]; // Tag array
  // Input rates
  real<lower=0, upper=1> u; // Initial tag retention rate (proportion)
  real<lower=0> v; // Tag loss rate
  real<lower=0> m; // Natural mortality rate
  // Prior parameters
  real<lower=0> h_alpha[A];
  real<lower=0> h_beta[A];
}

transformed data {
  // Transformed indexes
  int ST = T + (L - 1); // TODO: Rename (Number of survival time steps)
}

parameters {
  // Harvest rate
  real<lower=0, upper=1> h[A];
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
  real f[A]; // Fishing mortality
  for (ca in 1:A) {
    f[ca] = -log(1 - h[ca]);
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
  real p[G, A, A]; // [mg, ca, pa] Movement rates
	// Initialize values
	real s[G, ST, A]; // Survival rate
	int grainsize = 1;
	int release_steps[T];
	s = rep_array(0, G, ST, A);

	// Slipped in here
	p = rep_array(1e-12, G, A, A);
  // Length class 1
	p[1, 1, 1] = p11[1];
	p[1, 2, 1] = p11[2];
	p[1, 1, 2] = p12[1];
	p[1, 2, 2] = p12[2];
	p[1, 3, 2] = p12[3];
	p[1, 2, 3] = p13[1];
	p[1, 3, 3] = p13[2];

	// Populate release steps
	for (mt in 1:T) {
	  release_steps[mt] = mt;
	}

	// Compute survival
	for (mg in 1:G) {
		for (ct in 1:ST) {
			for (ca in 1:A) {
				s[mg, ct, ca] = exp(-f[ca] - m - v);
			}
		}
	}

	// Priors
  h ~ beta(h_alpha, h_beta);

	// Sampling statement
	// y_vec ~ poisson(y_hat);
	// Likelihood statement using reduce_sum()
	target += reduce_sum(partial_sum_lupmf, release_steps, grainsize,
	  T, A, G, L, u, f, p, s, x);
}

generated quantities {
  real p[G, A, A]; // [mg, ca, pa] Movement rates
  // Annual movement rates
  matrix[A, A] p_matrix[G];
  matrix[A, A] p_matrix_annual[G];
  real p_annual[A, A, G];

  // Slipped in here
	p = rep_array(1e-12, G, A, A);
  // Length class 1
	p[1, 1, 1] = p11[1];
	p[1, 2, 1] = p11[2];
	p[1, 1, 2] = p12[1];
	p[1, 2, 2] = p12[2];
	p[1, 3, 2] = p12[3];
	p[1, 2, 3] = p13[1];
	p[1, 3, 3] = p13[2];

  for (mg in 1:G) {
    // Populate p_matrix
    for (pa in 1:A) {
      for (ca in 1:A) {
        p_matrix[mg, pa, ca] = p[mg, ca, pa];
      }
    }
    // Populate p_matrix_annual
    p_matrix_annual[mg] = matrix_power(p_matrix[mg], Y);
    // Populate p_annual
    for (pa in 1:A) {
      for (ca in 1:A) {
        p_annual[pa, ca, mg] = p_matrix_annual[mg, pa, ca];
      }
    }
  }
}
