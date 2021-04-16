// Create simplex dimensions
int[] create_simplex_dimensions(int P, int A, int G, int[,] z) {
  int simplex_dimensions[6] = rep_array(0, 6);
  int row_sum;
  for (pa in 1:A) {
    row_sum = 0;
    for (ca in 1:A) {
      row_sum += z[pa, ca];
    }
    if (row_sum > 0) {
      simplex_dimensions[row_sum] += P * G;
    }
  }
  return simplex_dimensions;
}
