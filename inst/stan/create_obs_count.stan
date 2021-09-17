// Create observation count
int[] create_obs_count(
  int start,
  int end,
  int A,
  int G_released,
  int T_liberty,
  int T_study,
  int[,,,,] x
) {
  // Initialize
  int N = end - start + 1;
  int obs_count[N] = rep_array(0, N);
  // Populate observation count
  for (rt in start:end) {
    for (ra in 1:A) {
      for (rg in 1:G_released) {
        if (x[rt, ra, rg, 1, ra] > 0) {
          obs_count[rt - start + 1] += (min(T_liberty, T_study - rt) - 1) * A;
        }
      }
    }
  }
  // Return observation count
  return obs_count;
}
