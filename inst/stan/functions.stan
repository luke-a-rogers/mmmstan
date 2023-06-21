array[] int assemble_simplex_dims (array[,] int mindex) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[A] int simplex_dimensions = rep_array(0, A);
  int row_x_sum;
  // Populate simplex dimensions
  for (x in 1:X) {
    row_x_sum = sum(mindex[x]);
    if (row_x_sum > 0) {
      if (row_x_sum > A) {
        reject("row_x_sum: ", row_x_sum);
      }
      simplex_dimensions[row_x_sum] += 1;
    }
  }
  // Return value
  return simplex_dimensions;
}

array[,] vector assemble_tags_released(array[,,,,] int tags) {
  // Get dimensions
  int N = dims(tags)[1]; // Here N = model N - 1
  int L = dims(tags)[3];
  int X = dims(tags)[4];
  // Declare values
  array[N, L] vector[X] tags_released;
  // Populate tags released
  for (n in 1:N) { // Here N = model N - 1
    for (l in 1:L) { // Released size
      for (x in 1:X) { // Released region
        tags_released[n, l, x] = tags[n, 1, l, x, x] * 1.0;
      }
    }
  }
  // Return tags released
  return tags_released;
}

array [,,,,] int assemble_tags_transpose (array[,,,,] int tags) {
  // Get dimensions
  int N = dims(tags)[1]; // Here N = model N - 1
  int D = dims(tags)[2];
  int L = dims(tags)[3];
  int X = dims(tags)[4];
  // Declare values
  array[N, D, L, X, X] int tags_transpose;
  // Populate tag array
  for (n in 1:N) { // Here N = model N - 1
    for (d in 1:D) {
      for (l in 1:L) {
        for (x in 1:X) {
          for (y in 1:X) {
            tags_transpose[n, d, l, y, x] = tags[n, d, l, x, y];
          }
        }
      }
    }
  }
  // Return tag array
  return tags_transpose;
}

/**
* Assemble Matrices That Indicate Possible Movement at Each Duration Step
*
* @param movement_matrix, a matrix of dimension [X, X]
* @param D, an integer giving the maximum duration at large in steps
*
* @return an array of dimension [D] holding matrices of dimension [X, X]
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one model step.
*/
array[] matrix assemble_movement_possible (
  matrix mmatrix,
  int D
) {
  // Get dimensions
  int X = dims(mmatrix)[1];
  // Initialize values
  array[D] matrix[X, X] movement_possible;
  // Populate movement possible
  movement_possible[1] = rep_matrix(0.0, X, X);
  movement_possible[2] = mmatrix;
  // Iterate higher indexes
  for (d in 3:D) {
    if (min(movement_possible[d - 1]) > 0.0) {
      movement_possible[d] = movement_possible[d - 1];
    } else {
      movement_possible[d] = movement_possible[d - 1] * mmatrix;
    }
  }
  // Return movement possible
  return movement_possible;
}

array[] matrix assemble_movement_step (
  array[,] vector m1, // [L, ]
  array[,] vector m2, // [L, ]
  array[,] vector m3, // [L, ]
  array[,] vector m4, // [L, ]
  array[,] vector m5, // [L, ]
  array[,] vector m6, // [L, ]
  array[,] int mindex, // [X, X]
  int L
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[L] matrix[X, X] movement_step = rep_array(rep_matrix(0.0, X, X), L);
  array[L, A] int index = rep_array(0, L, A); // Simplex array row index
  int row_x_sum; // Movement index row sum
  int column; // Simplex index (column)
  // Populate movement step
  for (l in 1:L) {
    for (x in 1:X) {
      row_x_sum = sum(mindex[x]);
      if (row_x_sum > 0) {
        index[l, row_x_sum] += 1;
        column = 0;
        for (y in 1:X) {
          if (mindex[x, y] == 1) {
            column += 1;
            if (row_x_sum == 1) {
              movement_step[l, x, y] = m1[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 2) {
              movement_step[l, x, y] = m2[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 3) {
              movement_step[l, x, y] = m3[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 4) {
              movement_step[l, x, y] = m4[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 5) {
              movement_step[l, x, y] = m5[l, index[l, row_x_sum], column];
            } else if (row_x_sum == 6) {
              movement_step[l, x, y] = m6[l, index[l, row_x_sum], column];
            } else {
              reject("row_x_sum: ", row_x_sum);
            }
          }
        }
      }
    }
  }
  // Return movement step
  return movement_step;
}

array[,,] vector assemble_survival_step (
  array[] vector fishing_step,
  // array[] vector fishing_weight,
  vector selectivity,
  vector natural_mortality_step,
  real ongoing_loss_step,
  int K,
  int L
) {
  // Get dimensions
  int T = dims(fishing_step)[1];
  int X = dims(fishing_step)[2];
  // Initialize values
  array[T, K, L] vector[X] survival_step;
  // Populate survival step
  for (t in 1:T) {
    for (k in 1:K) {
      for (l in 1:L) {
        survival_step[t, k, l] = exp(
          -fishing_step[t] * selectivity[l] // .* fishing_weight[k] * selectivity[l]
          - natural_mortality_step
          - ongoing_loss_step
        );
      }
    }
  }
  // Return survival step
  return survival_step;
}

array[,] matrix assemble_transition_step (
  array[] matrix movement_step,
  array[,,] vector survival_step
) {
  // Get dimensions
  int T = dims(survival_step)[1];
  int K = dims(survival_step)[2];
  int L = dims(survival_step)[3];
  int X = dims(survival_step)[4];
  int N = T * K;
  // Declare values
  array[N, L] matrix[X, X] transition_step;
  int n = 1;
  // Populate transition step
  for (t in 1:T) {
    for (k in 1:K) {
      for (l in 1:L) {
        transition_step[n, l] = diag_pre_multiply( // A_n = A_{n-1}S_{n-1}\Gamma
          survival_step[t, k, l],
          movement_step[l]
        );
      }
      n += 1;
    }
  }
  // Return transition_step
  return transition_step;
}

array[,] vector assemble_observation_step (
  array[] vector fishing_step,
  // array[] vector fishing_weight,
  vector selectivity,
  vector reporting_step,
  int K,
  int L
) {
  // Get dimensions
  int T = dims(fishing_step)[1];
  int X = dims(fishing_step)[2];
  int N = T * K;
  // Declare values
  array[N, L] vector[X] observation_step;
  int n = 1;
  // Populate observation step
  for (t in 1:T) {
    for (k in 1:K) {
      for (l in 1:L) {
        observation_step[n, l] = reporting_step
        .* (1.0 - exp(-fishing_step[t] * selectivity[l])); // .* fishing_weight[k] * selectivity[l]
      }
      n += 1;
    }
  }
  // Return observation step
  return observation_step;
}

array[] vector assemble_fishing_weight (
  array[] vector fishing_weight_transpose
) {
  // Get dimensions
  int W = dims(fishing_weight_transpose)[2];
  int X = dims(fishing_weight_transpose)[1];
  // Declare fishing weight
  array[W] vector[X] fishing_weight;
  // Populate fishing weight
  for (w in 1:W) {
    for (x in 1:X) {
      fishing_weight[w, x] = fishing_weight_transpose[x, w];
    }
  }
  // Return fishing weight
  return fishing_weight;
}

array[] matrix assemble_movement_rate (
  array[] matrix movement_step,
  int K
) {
  // Get dimensions
  int L = dims(movement_step)[1];
  int X = dims(movement_step)[2];
  // Declare values
  array[L] matrix[X, X] movement_rate;
  // Populate movement rate
  for (l in 1:L) {
    movement_rate[l] = matrix_power(movement_step[l], K);
  }
  // Return movement rate
  return movement_rate;
}

array[] vector assemble_fishing_rate (
  array[] vector fishing_step,
  int K
) {
  // Get dimensions
  int T = dims(fishing_step)[1];
  int X = dims(fishing_step)[2];
  // Declare values
  array[T] vector[X] fishing_rate;
  // Populate fishing rate
  for (t in 1:T) {
    fishing_rate[t] = fishing_step[t] * K;
  }
  // Return fishing rate
  return fishing_rate;
}

array[] int index_n_to_r (int start, int end) {
  // Declare values
  array[end] int n_to_r;
  int r = 1;
  // Populate index array
  for (n in start:end) {
    n_to_r[n] = r;
    r += 1;
  }
  // Return array
  return n_to_r;
}

real partial_sum_lpmf (
  array[] int index,
  int start,
  int end,
  int N,
  int D,
  int L,
  int X,
  array[,,,,] int tags_transpose,
  array[,] vector tags_released,
  array[,] matrix transition_step,
  array[,] vector observation_step,
  array[] matrix movement_possible,
  real initial_loss_step,
  real tolerance_expected,
  real dispersion
) {
  // Declare index limits
  int R = (end - start + 1); // R stands in for N - 1
  int C = R * D * L * X * X;
  // Declare index arrays
  array[end] int n_to_r = index_n_to_r(start, end);
  // Declare enumeration values
  array[R, D, L] matrix[X, X] abundance;
  array[R, D, L] matrix[X, X] predicted;
  array[C] int observed;
  array[C] real expected;
  // Initialize count
  int count = 0;
  // Populate released abundance
  for (n in start:end) { // Model step
    for (l in 1:L) { // Released size
      abundance[n_to_r[n], 1, l] = diag_matrix(
        tags_released[n, l] * (1 - initial_loss_step)
      );
    }
  }
  // Compute expected recoveries
  for (n in start:end) { // Partial sum index range within released step
    for (d in 2:min(N - n + 1, D)) { // Duration at large
      for (l in 1:L) { // Released size
        // Propagate abundance
        abundance[n_to_r[n], d, l] = abundance[n_to_r[n], d - 1, l]
        * transition_step[n + d - 2, l]; // Previous step
        // Compute predicted
        predicted[n_to_r[n], d, l] = diag_post_multiply(
          abundance[n_to_r[n], d, l],
          observation_step[n + d - 1, l] // Current step
        );
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags_released[n, l, x] > 0) { // Were any tags released?
              if (movement_possible[d][x, y] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[n, d, l, y, x]; // Integer
                expected[count] = predicted[n_to_r[n], d, l, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End d
  } // End n
  // Return likelihood contribution
  return neg_binomial_2_lupmf(observed[1:count] | expected[1:count],dispersion);
}
