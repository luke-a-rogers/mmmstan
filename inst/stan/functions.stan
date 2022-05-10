matrix assemble_movement_matrix (array[,] int mindex) {
  // Get dimensions
  int X = dims(mindex)[1];
  // Declare movement matrix
  matrix[X, X] movement_matrix;
  matrix[X, X] movement_matrix_transpose;
  // Populate movement matrix transpose
  for (x in 1:X) {
    for (y in 1:X) {
      movement_matrix_transpose[y, x] = mindex[x, y] * 1.0;
    }
  }
  // Populate movement matrix
  movement_matrix = movement_matrix_transpose';
  // Return movement matrix
  return movement_matrix;
}

/**
* Assemble An Integer Array of Simplex Dimensions
*
* @param mindex, an array of dimension [X, X]
* @param I, positive integer number of movement steps
* @param D, positive integer number of movement sizes
*
* @return an integer array of dimension [A] giving the number of simplexes
* of each dimension 1:A
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one model step.
*/
array[] int assemble_simplex_dimensions (
  array[,] int mindex,
  int I,
  int D
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[A] int simplex_dimensions = rep_array(0, A);
  int row_sum;
  // Populate simplex dimensions
  for (x in 1:X) {
    row_sum = 0;
    for (y in 1:X) {
      row_sum += mindex[x, y];
    }
    if (row_sum > 0) {
      if (row_sum > A) {
        reject("row_sum: ", row_sum);
      }
      simplex_dimensions[row_sum] += I * D;
    }
  }
  // Return simplex dimensions
  return simplex_dimensions;
}

array[,] vector assemble_tags_released(array[,,,,] int tags) {
  // Get dimensions
  int N = dims(tags)[1]; // Here N = model N - 1
  int S = dims(tags)[2];
  int X = dims(tags)[4];
  // Declare values
  array[N, S] vector[X] tags_released;
  // Populate tags released
  for (n in 1:N) { // Here N = model N - 1
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        tags_released[n, s, x] = tags[n, s, 1, x, x] * 1.0;
      }
    }
  }
  // Return tags released
  return tags_released;
}

array [,,,,] int assemble_tags_transpose (array[,,,,] int tags) {
  // Get dimensions
  int N = dims(tags)[1]; // Here N = model N - 1
  int S = dims(tags)[2];
  int L = dims(tags)[3];
  int X = dims(tags)[4];
  // Declare values
  array[N, S, L, X, X] int tags_transpose;
  // Populate tag array
  for (n in 1:N) { // Here N = model N - 1
    for (s in 1:S) {
      for (l in 1:L) {
        for (x in 1:X) {
          for (y in 1:X) {
            tags_transpose[n, s, l, y, x] = tags[n, s, l, x, y];
          }
        }
      }
    }
  }
  // Return tag array
  return tags_transpose;
}

/**
* Assemble Matrices That Indicate Possible Movement at Each Liberty Step
*
* @param movement_matrix, a matrix of dimension [X, X]
* @param L, an integer giving the maximum steps at liberty
*
* @return an array of dimension [L] holding matrices of dimension [X, X]
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one model step.
*/
array[] matrix assemble_movement_possible (
  matrix movement_matrix,
  int L
) {
  // Get dimensions
  int X = dims(movement_matrix)[1];
  // Initialize values
  array[L] matrix[X, X] movement_possible;
  // Populate movement possible
  movement_possible[1] = rep_matrix(0.0, X, X);
  movement_possible[2] = movement_matrix;
  // Iterate higher indexes
  for (l in 3:L) {
    if (min(movement_possible[l - 1]) > 0.0) {
      movement_possible[l] = movement_possible[l - 1];
    } else {
      movement_possible[l] = movement_possible[l - 1]
      * movement_matrix;
    }
  }
  // Return movement possible
  return movement_possible;
}

/**
* Assemble An Array of Stepwise Movement Rate Matrices
*
* @param a1, ..., a8, arrays of stepwise movement rate simplexes
* @param mindex, an array of dimension [X, X]
* @param I, positive integer number of movement steps
* @param D, positive integer number of movement sizes
*
* @return an array of dimension [I, D] holding stepwise movement rate
* matrices [X, X]
*
* The arguments a1, ..., a8 are arrays of stepwise movement rate simplexes
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one model step.
*/
array[,] matrix assemble_movement_step (
  array[] vector a1,
  array[] vector a2,
  array[] vector a3,
  array[] vector a4,
  array[] vector a5,
  array[] vector a6,
  array[,] int mindex,
  int I,
  int D
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[I, D] matrix[X, X] movement_step = rep_array(rep_matrix(0.0,X,X),I,D);
  array[A] int index = rep_array(0, A); // Simplex array row index
  int row_sum; // Movement index row sum
  int column; // Simplex index (column)
  // Populate
  for (i in 1:I) {
    for (d in 1:D) {
      for (x in 1:X) {
        row_sum = sum(mindex[x, 1:X]);
        if (row_sum > 0) {
          index[row_sum] += 1;
          column = 0;
          for (y in 1:X) {
            if (mindex[x, y] == 1) {
              column += 1;
              if (row_sum == 1) {
                movement_step[i, d, x, y] = a1[index[row_sum], column];
              } else if (row_sum == 2) {
                movement_step[i, d, x, y] = a2[index[row_sum], column];
              } else if (row_sum == 3) {
                movement_step[i, d, x, y] = a3[index[row_sum], column];
              } else if (row_sum == 4) {
                movement_step[i, d, x, y] = a4[index[row_sum], column];
              } else if (row_sum == 5) {
                movement_step[i, d, x, y] = a5[index[row_sum], column];
              } else if (row_sum == 6) {
                movement_step[i, d, x, y] = a6[index[row_sum], column];
              } else {
                reject("row_sum: ", row_sum);
              }
            }
          }
        }
      }
    }
  }
  // Return movement term
  return movement_step;
}

array[,] matrix assemble_movement_deviation (
  array[,] matrix mstep,
  matrix mmean
) {
  // Get dimensions
  int I = dims(mstep)[1];
  int D = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Initialize values
  array[I, D] matrix[X, X] movement_deviation;
  // Populate movement deviation
  for (i in 1:I) {
    for (d in 1:D) {
      movement_deviation[i, d] = mstep[i, d] - mmean;
    }
  }
  // Return movement deviation
  return movement_deviation;
}

array[] vector assemble_fishing_weight (array[] vector fweight) {
  // Get dimensions
  int W = dims(fweight)[2];
  int X = dims(fweight)[1];
  // Declare fishing weight
  array[W] vector[X] fishing_weight;
  // Populate fishing weight
  for (w in 1:W) {
    for (x in 1:X) {
      fishing_weight[w, x] = fweight[x, w];
    }
  }
  // Return fishing weight
  return fishing_weight;
}

array[] vector assemble_fishing_rate (
  array[] vector fstep,
  int J
) {
  // Get dimensions
  int T = dims(fstep)[1];
  int X = dims(fstep)[2];
  // Declare values
  array[T] vector[X] fishing_rate;
  // Populate fishing rate
  for (t in 1:T) {
    fishing_rate[t] = fstep[t] * J;
  }
  // Return fishing rate
  return fishing_rate;
}


// New current above here ------------------------------------------------------



array[] int assemble_n_to_r (int start, int end) {
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
  int S,
  int L,
  int X,
  array[] int n_to_i,
  array[] int s_to_d,
  array[,,,,] int tags_transpose,
  array[,] vector tags_released,
  array[,] matrix movement_step,
  array[,] vector survival_step,
  array[,] vector observed_step,
  array[] matrix movement_possible,
  real initial_loss_step,
  real tolerance_expected,
  real dispersion
) {
  // Declare index limits
  int R = (end - start + 1); // R stands in for N - 1
  int C = R * S * L * X * X;
  // Declare index arrays
  array[end] int n_to_r = assemble_n_to_r(start, end);
  // Declare enumeration values
  array[R, S, L] matrix[X, X] abundance;
  array[R, S, L] matrix[X, X] predicted;
  array[C] int observed;
  array[C] real expected;
  // Declare index values
  int count;
  // Populate released abundance
  for (n in start:end) { // Model step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        abundance[n_to_r[n], s, 1] = diag_matrix(
          tags_released[n, s] * (1 - initial_loss_step)
        );
      }
    }
  }
  // Initialize count
  count = 0;
  // Compute expected recoveries
  for (n in start:end) { // Partial sum index range within released step
    for (s in 1:S) { // Released size
      for (l in 2:min(N - n + 1, L)) { // Liberty step
        // Propagate abundance
        abundance[n_to_r[n], s, l] = abundance[n_to_r[n], s, l - 1]
        * diag_post_multiply(
            movement_step[n_to_i[n], s_to_d[s]],
            survival_step[n + l - 2, s]
          );
        // Compute predicted
        predicted[n_to_r[n], s, l] = diag_post_multiply(
          abundance[n_to_r[n], s, l],
          observed_step[n + l - 1, s]
        );
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags_released[n, s, x] > 0) { // Were any tags released?
              if (movement_possible[l][x, y] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags_transpose[n, s, l, y, x]; // Integer
                expected[count] = predicted[n_to_r[n], s, l, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End s
  } // End n
  return neg_binomial_2_lupmf(observed[1:count] | expected[1:count],dispersion);
}

array[] matrix assemble_movement_term (
  array[,] matrix mstep,
  int K
) {
  // Get dimensions
  int I = dims(mstep)[1];
  int D = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Declare values
  array[I] matrix[X, X] movement_term;
  // Populate movement term
  for (i in 1:I) {
    movement_term[i] = matrix_power(mstep[i, 1], K);
  }
  // Return movement term
  return movement_term;
}



// Current above here ----------------------------------------------------------




array[,] matrix assemble_movement_rate (
  array[,] matrix mstep,
  int mpower
) {
  // Get dimensions
  int I = dims(mstep)[1];
  int S = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Initialize values
  array[I, S] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0,X,X),I,S);
  // Populate movement rate
  for (i in 1:I) {
    for (s in 1:S) {
      movement_rate[i, s] = matrix_power(mstep[i, s], mpower);
    }
  }
  // Return movement rate
  return movement_rate;
}




/**
* Assemble a movement mean matrix
*
* @param mstep, an array [N, S] of matrices [X, X]
* @param mpower, matrix power
*
* @return a matrix of dimension [X, X]
*/
matrix assemble_movement_mean (
  array[,] matrix mstep,
  int mpower
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int S = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Initialize values
  matrix[X, X] movement_mean = rep_matrix(0.0, X, X);
  matrix[X, X] movement_sum = rep_matrix(0.0, X, X);
  real movement_count = (N * 1.0) * (S * 1.0);
  // Populate movement sum
  for (n in 1:N) {
    for (s in 1:S) {
      for (y in 1:X) {
        for (x in 1:X) {
          movement_sum[x, y] += mstep[n, s, x, y];
        }
      }
    }
  }
  // Populate movement mean
  movement_mean = movement_sum / movement_count;
  // Possibly convert stepwise to rate
  movement_mean = matrix_power(movement_mean, mpower);
  // Return movement mean
  return movement_mean;
}

array[] matrix assemble_movement_time_rate (
  array[,] matrix mstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int T = ntime;
  int I = nterm;
  int S = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Initialize values
  array[T] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), T);
  matrix[X, X] movement_time;
  matrix[X, X] movement_mean;
  int n = 1;
  // Populate movement rate
  for (t in 1:T) {
    movement_time = diag_matrix(rep_vector(1.0, X));
    for (i in 1:I) {
      movement_mean = rep_matrix(0.0, X, X);
      for (s in 1:S) {
        movement_mean += mstep[n, s] / (S * 1.0);
      }
      movement_time = movement_time * movement_mean;
      n += 1;
    }
    movement_rate[t] = movement_time;
  }
  // Return movement rate
  return movement_rate;
}

array[] matrix assemble_movement_term_rate (
  array[,] matrix mstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int T = ntime;
  int I = nterm;
  int S = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Initialize values
  array[I] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), I);
  array[I] matrix[X, X] movement_sum = rep_array(rep_matrix(0.0, X, X), I);
  int n = 1;
  // Populate movement sum
  for (t in 1:T) {
    for (i in 1:I) {
      for (s in 1:S) {
        movement_rate[i] = movement_rate[i] + mstep[n, s] / (T * S * 1.0);
      }
      n += 1;
    }
  }
  // Return stepwise movement rate
  return movement_rate;
}

array[] matrix assemble_movement_size_rate (
  array[,] matrix mstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int T = ntime;
  int I = nterm;
  int S = dims(mstep)[2];
  int X = dims(mstep)[3];
  // Initialize values
  array[S] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), S);
  array[S] matrix[X, X] movement_prod;
  int n = 1;
  // Populate movement sum
  for (t in 1:T) {
    movement_prod = rep_array(diag_matrix(rep_vector(1.0, X)), S);
    for (i in 1:I) {
      for (s in 1:S) {
        movement_prod[s] = movement_prod[s] * mstep[n, s];
      }
      n += 1;
    }
    for (s in 1:S) {
      movement_rate[s] = movement_rate[s] + movement_prod[s] / (T * 1.0);
    }
  }
  // Return movement rate
  return movement_rate;
}

array[] matrix assemble_movement_time_deviation (
  array[] matrix mtime,
  matrix mmean
) {
  // Get dimensions
  int T = dims(mtime)[1];
  int X = dims(mtime)[2];
  // Declare values
  array[T] matrix[X, X] movement_deviation = rep_array(rep_matrix(0.0,X,X),T);
  // Populate movement deviation
  for (t in 1:T) {
    movement_deviation[t] = mtime[t] - mmean;
  }
  // Return movement
  return movement_deviation;
}

array[] matrix assemble_movement_term_deviation (
  array[] matrix mterm,
  matrix mmean
) {
  // Get dimensions
  int I = dims(mterm)[1];
  int X = dims(mterm)[2];
  // Declare values
  array[I] matrix[X, X] movement_deviation = rep_array(rep_matrix(0.0,X,X),I);
  // Populate movement deviation
  for (i in 1:I) {
    movement_deviation[i] = mterm[i] - mmean;
  }
  // Return movement
  return movement_deviation;
}

array[] matrix assemble_movement_size_deviation (
  array[] matrix msize,
  matrix mmean
) {
  // Get dimensions
  int S = dims(msize)[1];
  int X = dims(mmean)[2];
  // Declare values
  array[S] matrix[X, X] movement_deviation = rep_array(rep_matrix(0.0,X,X),S);
  // Populate movement deviation
  for (s in 1:S) {
    movement_deviation[s] = msize[s] - mmean;
  }
  // Return movement
  return movement_deviation;
}



array[] matrix assemble_fishing_term (array[] vector fsimp) {
  // Get dimensions
  int X = dims(fsimp)[1];
  int I = dims(fsimp)[2];
  // Declare array
  array[I] matrix[X, X] fishing_weight = rep_array(
    diag_matrix(rep_vector(1.0, X)),
    I
  );
  // Populate fishing weight
  for (i in 1:I) {
    for (x in 1:X) {
      fishing_weight[i, x, x] = fsimp[x, i] * I;
    }
  }
  // Return fishing weight
  return fishing_weight;
}

/**
* Assemble a survival step array [N] of diagonal matrices [X, X]
*
* @param fstep, an array [T] of vector [X] stepwise fishing rates
* @param fterm, an array [I] of matrix [X, X] stepwise fishing term weights
* @param sstep, a vector [X] of selectivity fractions
* @param mstep, a vector [X] of stepwise natural mortality rates
* @param ostep, a real ongoing stepwise tag loss rate
*
* @return an array [N, S] of matrices [X, X]
*/
array[,] matrix assemble_survival_step (
  array [] vector fstep,
  array [] matrix fterm,
  vector sstep,
  vector mstep,
  real ostep
) {
  // Get dimensions
  int T = dims(fstep)[1];
  int I = dims(fterm)[1];
  int S = dims(sstep)[1];
  int X = dims(fstep)[2];
  int N = T * I;
  // Initialize values
  array[N, S] matrix[X, X] survival_step = rep_array(rep_matrix(0.0, X, X),N,S);
  int n = 1;
  // Populate survival step
  for (t in 1:T) {
    for (i in 1:I) {
      for (s in 1:S) {
        survival_step[n, s] = diag_matrix(exp(-fterm[i] * fstep[t] * sstep[s]))
        * diag_matrix(exp(-mstep))
        * exp(-ostep);
      }
      n += 1;
    }
  }
  // Return survival step
  return survival_step;
}

/**
* Assemble an observation step array [N] of diagonal matrices [X, X]
*
* @param fstep, an array [T] of vector [X] stepwise fishing rates
* @param fterm, an array [I] of matrix [X, X] stepwise fishing term weights
* @param sstep, a vector [X] of selectivity fractions
* @param rstep, a vector [X] of tag reporting steps
*
* @return an array [N, S] of matrices [X, X]
*/
array[,] matrix assemble_observed_step (
  array [] vector fstep,
  array [] matrix fterm,
  vector sstep,
  vector rstep
) {
  // Get dimensions
  int T = dims(fstep)[1];
  int I = dims(fterm)[1];
  int S = dims(sstep)[1];
  int X = dims(fstep)[2];
  int N = T * I;
  // Initialize values
  array[N, S] matrix[X, X] observed_step = rep_array(rep_matrix(0.0, X, X),N,S);
  int n = 1;
  // Populate reporting step
  for (t in 1:T) {
    for (i in 1:I) {
      for (s in 1:S) {
        observed_step[n, s] = diag_matrix(1 - exp(-fterm[i] * fstep[t] * sstep[s]))
        * diag_matrix(rstep);
      }
      n += 1;
    }
  }
  // Return observed step
  return observed_step;
}



real partial_sum_lpmf (
  array[] int index,
  int start,
  int end,
  int X,
  int T,
  int I,
  int S,
  int N,
  int L,
  int C,
  array[,,,,] int tags,
  array[,] matrix movement_step,
  array[,] matrix survival_step,
  array[,] matrix observed_step,
  array[] matrix movement_possible_transpose,
  real initial_loss_step,
  real tolerance_expected,
  real dispersion
) {
  // Declare enumeration values
  array[N-1,S,L] matrix[X,X] abundance = rep_array(rep_matrix(0.0,X,X),N-1,S,L);
  array[N-1,S,L] matrix[X,X] predicted = rep_array(rep_matrix(0.0,X,X),N-1,S,L);
  array[C] int observed = rep_array(0, C);
  array[C] real expected = rep_array(0.0, C);
  // Declare index values
  int count;
  // Populate released abundance
  for (n in start:end) { // Partial sum index range within released step
    for (s in 1:S) { // Released size
      for (x in 1:X) { // Released region
        abundance[n, s, 1, x, x] = tags[n, s, 1, x, x] // Integer scalar
        * (1 - initial_loss_step); // Real scalar
      }
    }
  }
  // Initialize count
  count = 0;
  // Compute expected recoveries
  for (n in start:end) { // Partial sum index range within released step
    for (s in 1:S) { // Released size
      for (l in 2:min(N - n + 1, L)) { // Liberty step
        // Propagate abundance
        abundance[n, s, l] = abundance[n, s, l - 1]
        * survival_step[n + l - 2, s] // Previous step; diagonal matrix [X, X]
        * movement_step[n + l - 2, s]; // Prevous step; square matrix [X, X]
        // Compute predicted
        predicted[n, s, l] = abundance[n, s, l] // Square matrix [X, X]
        * observed_step[n + l - 1, s]; // Current step; diagonal matrix [X, X]
        // Compute vectors
        for (y in 1:X) { // Current region
          for (x in 1:X) { // Released region
            if (tags[n, s, 1, x, x] > 0) { // Were any tags released?
              if (movement_possible_transpose[l][y, x] > 0) {
                // Increment observation count
                count += 1;
                // Populate observed and expected values
                observed[count] = tags[n, s, l, x, y]; // Integer
                expected[count] = predicted[n, s, l, x, y]
                + tolerance_expected; // Real
              } // End if
            } // End if
          } // End x
        } // End y
      } // End l
    } // End s
  } // End n

  // Return partial sum
  return neg_binomial_2_lupmf(observed[1:count] | expected[1:count],dispersion);
}


