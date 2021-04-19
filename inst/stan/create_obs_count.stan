// Create observation count
int[] create_obs_count(
  int start,
  int end,
  int S,
  int A,
  int G,
  int L,
  int[,,,,] x
) {
  // Initialize
  int N = end - start + 1;
  int obs_count[N] = rep_array(0, N);
  // Populate observation count
  for (mt in start:end) {
    for (ma in 1:A) {
      for (mg in 1:G) {
        if (x[mt, ma, mg, 1, ma] > 0) {
          obs_count[mt - start + 1] += (min(L, S - mt) - 1) * A;
        }
      }
    }
  }
  // Return observation count
  return obs_count;
}
