
// Partial sum function for within-chain parallel threading
real partial_sum_lpmf(
  int[] release_steps_slice,
  int start,
  int end,
  int A,
  int G_released,
  int T_liberty,
  int T_study,
  int[,,] x,
  int[,,,,] y,
  real[,] w,
  real u,
  int[] p_time_index,
  int[] h_time_index,
  int[] h_group_index,
  int[] w_time_index,
  real y_fudge,
  real[,,] fs,
  real[,,] ss,
  real[,,,] ps,
  real phi
) {
  // Instantiate objects
  int N = end - start + 1; // Number of release steps for this slice
  int obs[N] = obs_count(start, end, A, G_released, T_liberty, T_study, x);
  int Y = sum(obs);
  int E = 0; // Number of realized steps at libery
  int y_ind = 1;
  int y_obs[Y];
  real y_hat[Y];
  real n[N, A, G_released, T_liberty, A]; // Predicted abundance

  // Initialize objects
  n = rep_array(rep_array(0, T_liberty, A), N, A, G_released);
  y_obs = rep_array(0, Y);
  y_hat = rep_array(0, Y);

	// Populate initial abundances
	for (rt in start:end) {
	  for (ra in 1:A) {
	    for (rg in 1:G_released) {
	      n[rt - start + 1, ra, rg, 1, ra] = u * x[rt, ra, rg];
	    }
	  }
	}

  // Compute predicted recoveries
  for (rt in start:end) {
  	for (ra in 1:A) {
	  	for (rg in 1:G_released) {
	  	  if (x[rt, ra, rg] > 0) {
			    E = min(T_liberty, T_study - rt);
  				for (lt in 2:E) { // Populate abundance array n
    				for (ca in 1:A) {
	    		  	for (pa in 1:A) {
		  		  	  n[rt - start + 1, ra, rg, lt, ca]
		  		  	  += n[rt - start + 1, ra, rg, lt - 1, pa]
	  				    * ss[rg, rt + lt - 2, pa]
		  				  * ps[rg, p_time_index[rt + lt - 2], ca, pa];
    				  } // End for pa
  				  } // End for ca
  			  } // End for lt
  				for (lt in 2:E) { // Compute recoveries and predicted recoveries
				    for (ca in 1:A) {
				      y_obs[y_ind] = y[rt, ra, rg, lt, ca];
  						y_hat[y_ind] = n[rt - start + 1, ra, rg, lt, ca]
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

  // Return partial sum
  // return poisson_lupmf(y_obs | y_hat);
  return neg_binomial_2_lupmf(y_obs | y_hat, phi);
}
