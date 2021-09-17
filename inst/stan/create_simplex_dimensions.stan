// Create simplex dimensions
int[] create_simplex_dimensions(
  int A,
  int G_released,
  int T_movement,
  int[,] z) {
  int simplex_dimensions[6] = rep_array(0, 6);
  int row_sum;
  for (pa in 1:A) {
    row_sum = 0;
    for (ca in 1:A) {
      row_sum += z[pa, ca];
    }
    if (row_sum > 0) {
      simplex_dimensions[row_sum] += T_movement * G_released;
    }
  }
  return simplex_dimensions;
}
