// Create stepwise movement rates p_step[ , , ca, pa]
real[,,,] create_p_step(
  int A,
  int G_released,
  int T_movement,
  int[,] z,
  real[,] p1,
  real[,] p2,
  real[,] p3,
  real[,] p4,
  real[,] p5,
  real[,] p6,
  real p_fudge) {
    // Initialize
    real p_step[G_released, T_movement, A, A] = rep_array(rep_array(p_fudge, A, A), G_released, T_movement);
    int sim_ind[6] = rep_array(0, 6); // Simplex index
    int row_sum; // z
    int col_ind; // z
    // Populate stepwise movement rates
    for (rg in 1:G_released) {
      for (mt in 1:T_movement) {
        for (pa in 1:A) {
          row_sum = 0;
          for (ca in 1:A) {
            row_sum += z[pa, ca];
          }
          if (row_sum > 0) {
            sim_ind[row_sum] += 1;
            col_ind = 0;
            for (ca in 1:A) {
              if (z[pa, ca] == 1) {
                col_ind += 1;
                if (row_sum == 1) {
                  p_step[rg, mt, ca, pa] = p1[sim_ind[row_sum], col_ind];
                } else if (row_sum == 2) {
                  p_step[rg, mt, ca, pa] = p2[sim_ind[row_sum], col_ind];
                } else if (row_sum == 3) {
                  p_step[rg, mt, ca, pa] = p3[sim_ind[row_sum], col_ind];
                } else if (row_sum == 4) {
                  p_step[rg, mt, ca, pa] = p4[sim_ind[row_sum], col_ind];
                } else if (row_sum == 5) {
                  p_step[rg, mt, ca, pa] = p5[sim_ind[row_sum], col_ind];
                } else if (row_sum == 6) {
                  p_step[rg, mt, ca, pa] = p6[sim_ind[row_sum], col_ind];
                }
              }
            }
          }
        }
      }
    }
    return p_step;
}
