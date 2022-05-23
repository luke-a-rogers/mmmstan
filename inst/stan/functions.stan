array[] int simplex_dims (array[,] int mindex) {
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

array[] int simplex_dims (
  array[,] int mindex,
  int use_tran,
  int B
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[A] int simplex_dimensions = rep_array(0, A);
  int row_x_sum;
  // Conditional
  if (use_tran == 1) {
    // Populate simplex dimensions
    for (x in 1:X) {
      row_x_sum = sum(mindex[x]);
      if (row_x_sum > 0) {
        if (row_x_sum > A) {
          reject("row_x_sum: ", row_x_sum);
        }
        simplex_dimensions[row_x_sum] += B * row_x_sum;
      }
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
* Assemble Matrices That Indicate Possible Movement at Each Liberty Step
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

matrix assemble_movement_step_mean (
  array[] vector m1,
  array[] vector m2,
  array[] vector m3,
  array[] vector m4,
  array[] vector m5,
  array[] vector m6,
  array[,] int mindex
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  matrix[X, X] movement_step_mean = rep_matrix(0.0, X, X);
  array[A] int index = rep_array(0, A); // Simplex array row index
  int row_x_sum; // Movement index row sum
  int column; // Simplex index (column)
  // Populate
  for (x in 1:X) {
    row_x_sum = sum(mindex[x]);
    if (row_x_sum > 0) {
      index[row_x_sum] += 1;
      column = 0;
      for (y in 1:X) {
        if (mindex[x, y] == 1) {
          column += 1;
          if (row_x_sum == 1) {
            movement_step_mean[x, y] = m1[index[row_x_sum], column];
          } else if (row_x_sum == 2) {
            movement_step_mean[x, y] = m2[index[row_x_sum], column];
          } else if (row_x_sum == 3) {
            movement_step_mean[x, y] = m3[index[row_x_sum], column];
          } else if (row_x_sum == 4) {
            movement_step_mean[x, y] = m4[index[row_x_sum], column];
          } else if (row_x_sum == 5) {
            movement_step_mean[x, y] = m5[index[row_x_sum], column];
          } else if (row_x_sum == 6) {
            movement_step_mean[x, y] = m6[index[row_x_sum], column];
          } else {
            reject("row_x_sum: ", row_x_sum);
          }
        }
      }
    }
  }
  // Return movement step mean
  return movement_step_mean;
}

array[,] matrix assemble_movement_transformation (
  array[] vector s1,
  array[] vector s2,
  array[] vector s3,
  array[] vector s4,
  array[] vector s5,
  array[] vector s6,
  array[,] int mindex,
  int use_tran,
  int B
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int A = 6;
  // Declare values
  array[B, X] matrix[X, X] movement_tran = rep_array(identity_matrix(X), B, X);
  array[A] int index = rep_array(0, A); // Simplex array row index
  int row_w_sum; // Movement index row sum
  int column; // Simplex index (column)
  // Conditional
  if (use_tran == 1) {
    // Populate
    for (b in 1:B) { // One of T, K, L
      for (w in 1:X) { // Row of movement step mean to transform
        row_w_sum = sum(mindex[w]);
        if (row_w_sum > 0) {
          for (x in 1:X) { // The row of the transformation matrix
            if (mindex[w, x] == 1) {
              index[row_w_sum] += 1;
              column = 0;
              for (y in 1:X) { // The column of the transformation matrix
                if (mindex[w, y] == 1) {
                  column += 1;
                  if (row_w_sum == 1) {
                    movement_tran[b, w, x, y] = s1[index[row_w_sum], column];
                  } else if (row_w_sum == 2) {
                    movement_tran[b, w, x, y] = s2[index[row_w_sum], column];
                  } else if (row_w_sum == 3) {
                    movement_tran[b, w, x, y] = s3[index[row_w_sum], column];
                  } else if (row_w_sum == 4) {
                    movement_tran[b, w, x, y] = s4[index[row_w_sum], column];
                  } else if (row_w_sum == 5) {
                    movement_tran[b, w, x, y] = s5[index[row_w_sum], column];
                  } else if (row_w_sum == 6) {
                    movement_tran[b, w, x, y] = s6[index[row_w_sum], column];
                  } else {
                    reject("row_w_sum: ", row_w_sum);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // Return movement tran
  return movement_tran;
}

array[,,] matrix assemble_movement_step (
  matrix mmean,
  array[,] matrix mtime,
  array[,] matrix mterm,
  array[,] matrix msize,
  int model_time,
  int model_term,
  int model_size
) {
  // Get dimensions
  int T = dims(mtime)[1];
  int K = dims(mterm)[1];
  int L = dims(msize)[1];
  int X = dims(mmean)[1];
  // Declare movement step
  array[T, K, L] matrix[X, X] movement_step;
  // Populate movement step
  if (!model_time && !model_term && !model_size) {
    movement_step = rep_array(mmean, T, K, L);
  } else if (!model_time && !model_term && model_size) {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x] * msize[l, x];
          }
        }
      }
    }
  } else if (!model_time && model_term && !model_size) {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x] * mterm[k, x];
          }
        }
      }
    }
  } else if (model_time && !model_term && !model_size) {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x] * mtime[t, x];
          }
        }
      }
    }
  } else if (!model_time && model_term && model_size) {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x]
            * mterm[k, x]
            * msize[l, x];
          }
        }
      }
    }
  } else if (model_time && model_term && !model_size) {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x]
            * mtime[t, x]
            * mterm[k, x];
          }
        }
      }
    }
  } else if (model_time && !model_term && model_size) {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x]
            * mtime[t, x]
            * msize[l, x];
          }
        }
      }
    }
  } else {
    for (t in 1:T) { // Year
      for (k in 1:K) { // Term
        for (l in 1:L) { // Size
          for (x in 1:X) { // Row within matrix movement step mean
            movement_step[t, k, l][x] = mmean[x]
            * mtime[t, x]
            * mterm[k, x]
            * msize[l, x];
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
  // vector selectivity,
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
          -fishing_step[t] // .* fishing_weight[k] * selectivity[l]
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
  array[,,] matrix movement_step,
  array[,,] vector survival_step
) {
  // Get dimensions
  int T = dims(movement_step)[1];
  int K = dims(movement_step)[2];
  int L = dims(movement_step)[3];
  int X = dims(movement_step)[4];
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
          movement_step[t, k, l]
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
  // vector selectivity,
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
        .* (1.0 - exp(-fishing_step[t])); // .* fishing_weight[k] * selectivity[l]
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

array[] vector assemble_fishing_rate (
  array[] vector fishing_step,
  int J
) {
  // Get dimensions
  int T = dims(fishing_step)[1];
  int X = dims(fishing_step)[2];
  // Declare values
  array[T] vector[X] fishing_rate;
  // Populate fishing rate
  for (t in 1:T) {
    fishing_rate[t] = fishing_step[t] * J;
  }
  // Return fishing rate
  return fishing_rate;
}

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

array[] matrix assemble_movement_rate (
  matrix movement_step_mean,
  array[,] matrix movement_tran,
  int use_tran,
  int R
) {
  // Get dimensions
  int B = dims(movement_tran)[1];
  int X = dims(movement_step_mean)[1];
  // Declare value
  array[B] matrix[X, X] movement_step;
  array[B] matrix[X, X] movement_rate;
  // Conditional
  if (use_tran) {
    // Populate value
    for (b in 1:B) {
      for (x in 1:X) {
        movement_step[b][x] = movement_step_mean[x] * movement_tran[b, x];
      }
      movement_rate[b] = matrix_power(movement_step[b], R);
    }
  } else {
    movement_rate = rep_array(matrix_power(movement_step_mean, R), B);
  }
  // Return value
  return movement_rate;
}


// New current above here ------------------------------------------------------





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







/**

array[,] matrix assemble_observed_step_old (
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

array[,] matrix assemble_survival_step_old (
  array [] vector fstep,
  array [] matrix fterm, // fweight
  vector sstep, // selectivity
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

*/
