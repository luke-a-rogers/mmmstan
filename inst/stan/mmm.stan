
functions {
  // Matrix power (Remove upon release RStan 2.3)
  matrix mat_power(matrix a, int n);
  matrix mat_power(matrix a, int n) {
    if (n == 0)
      return diag_matrix(rep_vector(1, rows(a)));
    else if (n == 1)
      return a;
    else
      return a *  mat_power(a, n - 1);
  } // End mat_power
}

data {
  // Index limits
	int A; // Number of release and recovery areas
	int G; // Number of release groups
	int L; // Maximum number of time steps at liberty
	int T; // Number of release time steps

	// Tag data & movement index
	int x[T, A, G, L, A]; // Tag array

	// Input rates
	real d; // Tags attached after release (proportion)
	real f; // Fishing mortality rate
	real h; // Tag loss rate
	real m; // Natural mortality rate
}

transformed data {
  int ST = T + (L - 1);
  int Y = 0;
  for (mt in 1:T) {
		for (ma in 1:A) {
			for (mg in 1:G) {
			  if (x[mt, ma, mg, 1, ma] > 0) {
			    Y += L * A;
			  }
			}
		}
  }
}

parameters {
	// real<lower=0> f; // Fishing mortality rate
	// Length class 1
	simplex[2] p11;
	simplex[3] p12;
	simplex[2] p13;
	// Length class 2
	simplex[2] p21;
	simplex[3] p22;
	simplex[2] p23;
	// Length class 3
	simplex[2] p31;
	simplex[3] p32;
	simplex[2] p33;
}

transformed parameters {
  real p[G, A, A]; // [mg, ca, pa] Movement rates
	p = rep_array(1e-12, G, A, A);
  // Length class 1
	p[1, 1, 1] = p11[1];
	p[1, 2, 1] = p11[2];
	p[1, 1, 2] = p12[1];
	p[1, 2, 2] = p12[2];
	p[1, 3, 2] = p12[3];
	p[1, 2, 3] = p13[1];
	p[1, 3, 3] = p13[2];
	// Length class 1
	p[2, 1, 1] = p21[1];
	p[2, 2, 1] = p21[2];
	p[2, 1, 2] = p22[1];
	p[2, 2, 2] = p22[2];
	p[2, 3, 2] = p22[3];
	p[2, 2, 3] = p23[1];
	p[2, 3, 3] = p23[2];
	// Length class 1
	p[3, 1, 1] = p31[1];
	p[3, 2, 1] = p31[2];
	p[3, 1, 2] = p32[1];
	p[3, 2, 2] = p32[2];
	p[3, 3, 2] = p32[3];
	p[3, 2, 3] = p33[1];
	p[3, 3, 3] = p33[2];
}

model {
	// Initialize values
	real n[T, A, G, L, A]; // Predicted abundance
	real s[G, ST, A]; // Survival rate
	int y_vec[Y]; // Recoveries
	real y_hat[Y]; // Predicted recoveries
	int y_ind; // Recovery index counter
	real n_sub[T, A];
	n = rep_array(rep_array(0, L, A), T, A, G);
	s = rep_array(0, G, ST, A);
	y_vec = rep_array(0, Y);
	y_hat = rep_array(0, Y);
	y_ind = 1;

	// Compute survival
	for (mg in 1:G) {
		for (ct in 1:ST) {
			for (ca in 1:A) {
				s[mg, ct, ca] = exp(-f - m - h);
			}
		}
	}

	// Populate initial abundances
	for (mt in 1:T) {
	  for (ma in 1:A) {
	    for (mg in 1:G) {
	      n[mt, ma, mg, 1, ma] = d * x[mt, ma, mg, 1, ma];
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
		  					* s[mg, mt + cl - 2, pa] * p[mg, ca, pa];
	  					} // End for pa
  					} // End for ca
  				} // End for cl
  				for (cl in 2:L) { // Compute recoveries and predicted recoveries
  					for (ca in 1:A) {
  					  y_vec[y_ind] = x[mt, ma, mg, cl, ca];
  						y_hat[y_ind] = n[mt, ma, mg, cl, ca] * (1 - exp(-f));
  						y_ind += 1;
  					} // End for ra
  				} // End for cl
			  } // End if
			} // End for mg
		} // End for ma
	} // End for mt

	// Sampling statement
	y_vec ~ poisson(y_hat);
}

generated quantities {
  // Annual movement rates
  matrix[A, A] p_matrix[G];
  matrix[A, A] p_matrix_annual[G];
  real p_annual[A, A, G];
  for (mg in 1:G) {
    // Populate p_matrix
    for (pa in 1:A) {
      for (ca in 1:A) {
        p_matrix[mg, pa, ca] = p[mg, ca, pa];
      }
    }
    // Populate p_matrix_annual
    p_matrix_annual[mg] = mat_power(p_matrix[mg],  12);
    // Populate p_annual
    for (pa in 1:A) {
      for (ca in 1:A) {
        p_annual[pa, ca, mg] = p_matrix_annual[mg, pa, ca];
      }
    }
  }
}
