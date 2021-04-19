
// Partial sum function for within-chain parallel threading
real partial_sum_lpmf(
  int[] release_steps_slice,
  int start,
  int end,
  int S,
  int A,
  int G,
  int L,
  real u,
  real y_phi,
  real y_fudge,
  int[] h_index,
  int[] p_index,
  int[] q_index,
  int[] w_index,
  real[,] w,
  real[,,] f_step,
  real[,,] s_step,
  real[,,,] p_step,
  int[,,,,] x
) {
  // Instantiate objects
  int N = end - start + 1; // Number of release steps for this slice
  int obs_count[N] = create_obs_count(start, end, S, A, G, L, x);
  int R = sum(obs_count);
  int E = 0; // Number of realized steps at libery
  int y_ind = 1;
  int y_obs[R];
  real y_hat[R];
  real n[N, A, G, L, A]; // Predicted abundance

  // Initialize objects
  n = rep_array(rep_array(0, L, A), N, A, G);
  y_obs = rep_array(0, R);
  y_hat = rep_array(0, R);

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
	  	  if (x[mt, ma, mg, 1, ma] > 0) {
			    E = min(L, S - mt);
  				for (cl in 2:E) { // Populate abundance array n
    				for (ca in 1:A) {
	    		  	for (pa in 1:A) {
		  		  	  n[mt - start + 1, ma, mg, cl, ca]
		  		  	  += n[mt - start + 1, ma, mg, cl - 1, pa]
	  				    * s_step[mg, mt + cl - 2, pa]
		  				  * p_step[mg, p_index[mt + cl - 2], ca, pa];
    				  } // End for pa
  				  } // End for ca
  			  } // End for cl
  				for (cl in 2:E) { // Compute recoveries and predicted recoveries
				    for (ca in 1:A) {
				      y_obs[y_ind] = x[mt, ma, mg, cl, ca];
  						y_hat[y_ind] = n[mt - start + 1, ma, mg, cl, ca]
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

  // Return partial sum
  // return poisson_lupmf(y_obs | y_hat);
  return neg_binomial_2_lupmf(y_obs | y_hat, y_phi);
}
