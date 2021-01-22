
data {
	// Index limits
	int A; // Number of release and recovery areas
	int G; // Number of release groups
	int L; // Maximum number of time steps at liberty
	int T; // Number of release time steps

	// Tag data & movement index
	int x[T, A, G]; // Tag releases
	int y[T, A, G, L, A]; // Tag recoveries

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
			  if (x[mt, ma, mg] > 0) {
			    Y += L * A;
			  }
			}
		}
  }
}

parameters {
	// real<lower=0> f; // Fishing mortality rate
	// Movement parameters
	simplex[3] p1;
	simplex[3] p2;
	simplex[3] p3;
}

transformed parameters {
  real p[A, A]; // [ca, pa]
  for (ca in 1:A) {
    p[ca, 1] = p1[ca];
    p[ca, 2] = p2[ca];
    p[ca, 3] = p3[ca];
  }
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
	      n[mt, ma, mg, 1, ma] = d * x[mt, ma, mg];
	    }
	  }
	}

	// Compute predicted recoveries
	for (mt in 1:T) {
		for (ma in 1:A) {
			for (mg in 1:G) {
			  if (x[mt, ma, mg] > 0) {
  				for (cl in 2:L) { // Populate abundance array n
	  				for (ca in 1:A) {
		  				for (pa in 1:A) {
			  				n[mt, ma, mg, cl, ca] += n[mt, ma, mg, cl - 1, pa]
		  					* s[mg, mt + cl - 2, pa] * p[ca, pa];
	  					} // End for pa
  					} // End for ca
  				} // End for cl
  				for (cl in 2:L) { // Compute recoveries and predicted recoveries
  					for (ca in 1:A) {
  					  y_vec[y_ind] = y[mt, ma, mg, cl, ca];
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
