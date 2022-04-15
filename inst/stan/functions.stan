array[] assemble_simplex_dimensions (
  array[,] mindex,
  int ntime,
  int nterm,
  int nsize
) {
  // Get dimensions
  int X = dims(mindex)[1];
  // Declare values
  array[6] int simplex_dimensions = rep_array(0, 6);
  int row_sum;
  // Populate simplex dimensions
  for (x in 1:X) {
    row_sum = 0;
    for (y in 1:X) {
      row_sum += mindex[x, y];
    }
    if (row_sum > 0) {
      simplex_dimensions[row_sum] += ntime * nterm * nsize;
    }
  }
  // Return simplex dimensions
  return simplex_dimensions;
}

array[,] matrix assemble_movement_step(
  array[] vector a1,
  array[] vector a2,
  array[] vector a3,
  array[] vector a4,
  array[] vector a5,
  array[] vector a6,
  array[,] mindex,
  int ntime,
  int nterm,
  int nsize,
  real tolm
) {
  // Get dimensions
  int N = ntime * nterm;
  int T = ntime;
  int I = nterm;
  int S = nsize;
  // Initialize movement step
  array[N, S] matrix[X, X] movement_step = rep_array(rep_matrix(tolm,X,X),N,S);
  // Initialize indexes
  array[6] int index = rep_array(0, 6); // Simplex array row index
  int row_sum; // Movement index row sum
  int column; // Simplex index (column)
  // Populate movement step
  for (n in 1:N) {
    for (s in 1:S) {
      for (x in 1:X) {
        row_sum = sum(mindex[x, 1:X]);
        if (row_sum > 0) {
          index[row_sum] += 1;
          column = 0;
          for (y in 1:X) {
            if (mindex[x, y] == 1) {
              column += 1;
              if (row_sum == 1) {
                movement_step[n, s, x, y] = a1[index[row_sum], column];
              } else if (row_sum == 2) {
                movement_step[n, s, x, y] = a2[index[row_sum], column];
              } else if (row_sum == 3) {
                movement_step[n, s, x, y] = a3[index[row_sum], column];
              } else if (row_sum == 4) {
                movement_step[n, s, x, y] = a4[index[row_sum], column];
              } else if (row_sum == 5) {
                movement_step[n, s, x, y] = a5[index[row_sum], column];
              } else if (row_sum == 6) {
                movement_step[n, s, x, y] = a6[index[row_sum], column];
              } else {
                reject("row_sum: ", row_sum);
              }
            }
          }
        }
      }
    }
  }
  // Return movement step
  return movement_step;
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
  array[] real mstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int T = ntime;
  int I = nterm;
  int S = dims(mstep[2]);
  int X = dims(mstep[3]);
  // Initialize values
  array[T] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), T);
  matrix[X, X] movement_time;
  matrix[X, X] movement_mean;
  int n = 1;
  // Populate movement rate
  for (t in 1:T) {
    movement_time = diag_matrix(1.0, X, X);
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
  array[] real mstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int T = ntime;
  int I = nterm;
  int S = dims(mstep[2]);
  int X = dims(mstep[3]);
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
  array[] real mstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(mstep)[1];
  int T = ntime;
  int I = nterm;
  int S = dims(mstep[2]);
  int X = dims(mstep[3]);
  // Initialize values
  array[S] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), S);
  array[S] matrix[X, X] movement_prod;
  int n = 1;
  // Populate movement sum
  for (t in 1:T) {
    movement_prod = rep_array(diag_matrix(1.0, X, X), S)
    for (i in 1:I) {
      for (s in 1:S) {
        movement_prod[s] = movement_prod[s] * mstep[n, s];
      }
      n += 1;
    }
    for (s in 1:S) {
      movement_rate[s] = movement_sum[s] + movement_prod[s] / (T * 1.0);
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
  int X = dims(mterm)[2];
  // Declare values
  array[S] matrix[X, X] movement_deviation = rep_array(rep_matrix(0.0,X,X),S);
  // Populate movement deviation
  for (s in 1:S) {
    movement_deviation[s] = msize[s] - mmean;
  }
  // Return movement
  return movement_deviation;
}








/**
* Assemble a movement term rate array
*
* @param pmean, an array of dimension [P]
* @param pterm, an array of dimension [I, P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [I] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The second argument holds the parameter deviations
* associated with movement rate deviations from the mean for term, pterm.
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year) and acts as the matrix power for multiplying
* movement matrices [X, X].
*/
array[] matrix assemble_movement_term_rate (
  array[] real pmean,
  array[,] real pterm,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int I = dims(pterm)[1];
  int X = dims(mindex)[1];
  // Initialize arrays
  array[I, 1] matrix[X, X] movement_step = assemble_movement_step(
    pmean,
    rep_array(0.0, 1, P),
    pterm,
    rep_array(0.0, 1, P),
    mindex
  );
  array[I] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), I);
  // Populate movement rates
  for (i in 1:I) {
    movement_rate[i] = matrix_power(movement_step[i, 1], nterm);
  }
  // Return movement rate
  return movement_rate;
}

/**
* Assemble a movement size rate array
*
* @param pmean, an array of dimension [P]
* @param psize, an array of dimension [S, P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [S] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The second argument holds the parameter deviations
* associated with movement rate deviations from the mean for size, psize.
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year) and acts as the matrix power for multiplying
* movement matrices [X, X].
*/
array[] matrix assemble_movement_size_rate (
  array[] real pmean,
  array[,] real psize,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int S = dims(psize)[1];
  int X = dims(mindex)[1];
  // Initialize arrays
  array[1, S] matrix[X, X] movement_step = assemble_movement_step(
    pmean,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    psize,
    mindex
  );
  array[S] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), S);
  // Populate movement rates
  for (s in 1:S) {
    movement_rate[s] = matrix_power(movement_step[1, s], nterm);
  }
  // Return movement rate
  return movement_rate;
}








/**
* Inverse multinomial logit
*
* @param params, an array of dimension [U]
*
* @return an array of dimension [V] that sums to one
*/
array[] real inv_multi_logit (array[] real params) {
  // Get dimension
  int U = size(params); // Number of parameters (in a given mindex row)
  int V = U + 1; // Number of resultant movement rates
  // Declare scalars
  real sum_exp_params = 0.0;
  real sum_movements_u = 0.0;
  // Declare simplex
  array[V] real movements = rep_array(0.0, V);
  // Populate sum
  for (u in 1:U) {
    sum_exp_params += exp(params[u]);
  }
  // Populate movement rates
  for (u in 1:U) {
    movements[u] = exp(params[u]) / (1.0 + sum_exp_params);
    sum_movements_u += movements[u];
  }
  movements[V] = 1.0 - sum_movements_u;
  // Check movement rates
  if (sum(movements) > 1) {
    print("params: ", params);
    print("movements: ", movements);
    reject("sum(movements): ", sum(movements));
  }
  if (sum(movements) < 0) {
    print("params: ", params);
    print("movements: ", movements);
    reject("sum(movements): ", sum(movements));
  }
  if (max(movements) > 1) {
    print("params: ", params);
    print("movements: ", movements);
    reject("max(movements): ", max(movements));
  }
  if (min(movements) < 0) {
    print("params: ", params);
    print("movements: ", movements);
    reject("min(movements): ", min(movements));
  }
  // Return movement rates
  return movements;
}

/**
* Inverse multi-logit
*
* @param arg, an array of dimension [P]
*
* @return an array of dimension [P + 1] that sums to one
*/
array[] real inverse_multi_logit (array[] real arg) {
  // Get dimensions
  int P = size(arg);
  // Instantiate sums
  real sum_exp_arg = 0;
  real sum_value = 0;
  // Instantiate value
  array[P + 1] real value = rep_array(0.0, P + 1);
  // Check dimensions
  if (P < 1) {reject("P must not be < 1; found P = ", P);}
  // Get sum of exponential
  for (p in 1:P) {
    sum_exp_arg += exp(arg[p]);
  }
  // Populate value
  for (p in 1:P) {
    value[p] = exp(arg[p]) / (1 + sum_exp_arg);
    sum_value += value[p];
  }
  value[P + 1] = 1 - sum_value;
  // Return value
  return value;
}

/**
* Assemble a movement parameter array 'slab'
*
* @param pmean, an array of dimension [P]
* @param ptime, an array of dimension [T, P]
* @param pterm, an array of dimension [I, P]
* @param psize, an array of dimension [S, P]
*
* @return an array of dimension [N, S, P]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The remaining arguments hold the parameter deviations
* associated with movement rate deviations from the mean for time (ptime),
* term (pterm), and size (psize).
*
* Parameter slabs associated with movement rate means only,
* or movement rates plus certain types of deviations can be
* assembled by passing rep_array(0.0, 1, P) for arguments corresponding
* to the deviations to exclude.
*
* Note that N = T * I, which is the number of study steps and one more than
* the number of released steps, N - 1.
*/
array[,,] real assemble_movement_parameter_slab (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  // Set N and instantiate index
  int N = T * I;
  int n;
  // Instantiate slabs
  array[N, S, P] real parameter_slab = rep_array(0.0, N, S, P);
  // Populate parameter slab
  n = 1;
  for (t in 1:T) {
    for (i in 1:I) {
      for (s in 1:S) {
        for (p in 1:P) {
          parameter_slab[n, s, p] = pmean[p]
          + ptime[t, p]
          + pterm[i, p]
          + psize[s, p];
        }
      }
      n += 1;
    }
  }
  // Return parameter slab
  return parameter_slab;
}

/**
* Assemble a movement array 'slab'
*
* @param pmean, an array of dimension [P]
* @param ptime, an array of dimension [T, P]
* @param pterm, an array of dimension [I, P]
* @param psize, an array of dimension [S, P]
* @param mindex, an array of dimension [X, X]
*
* @return an array of dimension [N, S, M]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The remaining arguments hold the parameter deviations
* associated with movement rate deviations from the mean for time (ptime),
* term (pterm), and size (psize).
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* Movement slabs associated with movement rate means only,
* or movement rates plus certain types of deviations can be
* assembled by passing rep_array(0.0, 1, P) for arguments corresponding
* to the deviations to exclude.
*
* Note that N = T * I, which is the number of study steps and one more than
* the number of released steps, N - 1.
*/
array[,,] real assemble_movement_slab (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] int mindex
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  int X = dims(mindex)[1];
  // Set dimensions
  int M = P + X;
  int N = T * I;
  // Instantiate index ranges
  array[2] int p;
  array[2] int m;
  // Initialize slabs
  array[N, S, P] real parameter_slab = assemble_movement_parameter_slab(
    pmean,
    ptime,
    pterm,
    psize
  );
  array[N, S, M] real movement_slab = rep_array(0.0, N, S, M);
  // Populate movement slab
  for (n in 1:N) { // N steps
    for (s in 1:S) { // S size classes
      p[1] = 1;
      m[1] = 1;
      for (x in 1:X) { // X sets of parameters to convert
        // Set the upper limit of the index ranges
        p[2] = p[1] + sum(mindex[x]) - 2;
        m[2] = m[1] + sum(mindex[x]) - 1;
        // Assign values to the movement slab
        movement_slab[n, s, m[1]:m[2]] = inv_multi_logit(
          parameter_slab[n, s, p[1]:p[2]]
        );
        // Update the lower limit of the index ranges
        p[1] = p[2] + 1;
        m[1] = m[2] + 1;
      }
    }
  }
  // Return movement slab
  return movement_slab;
}

/**
* Assemble a movement array at timescale 'step'
*
* @param pmean, an array of dimension [P]
* @param ptime, an array of dimension [T, P]
* @param pterm, an array of dimension [I, P]
* @param psize, an array of dimension [S, P]
* @param mindex, an array of dimension [X, X]
*
* @return an array of dimension [N, S] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The remaining arguments hold the parameter deviations
* associated with movement rate deviations from the mean for time (ptime),
* term (pterm), and size (psize).
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* Movement associated with movement means only, or movement plus certain
* types of deviations can be assembled by passing rep_array(0.0, 1, P)
* for arguments corresponding to the deviations to exclude.
*
* Note that N = T * I, which is the number of study steps and one more than
* the number of released steps, N - 1.
*/
array[,] matrix assemble_movement_step (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] int mindex
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  int X = dims(mindex)[1];
  // Set dimensions
  int M = P + X;
  int N = T * I;
  // Initialize arrays
  array[N, S, M] real movement_slab = assemble_movement_slab(
    pmean,
    ptime,
    pterm,
    psize,
    mindex
  );
  array[N, S] matrix[X, X] movement_step = rep_array(rep_matrix(0.0,X,X),N,S);
  // Instantiate matrix
  matrix[X, X] transpose_step;
  // Instantiate count
  int m;
  // Populate movement step
  for (n in 1:N) {
    for (s in 1:S) {
      m = 1;
      // Accommodate column-major indexing for matrices
      transpose_step = rep_matrix(0.0, X, X);
      for (x in 1:X) {
        for (y in 1:X) {
          if (mindex[x, y] == 1) {
            transpose_step[y, x] = movement_slab[n, s, m];
            m += 1;
          }
        }
      }
      movement_step[n, s] = transpose_step'; // Transpose of transpose
    }
  }
  // Return movement step
  return movement_step;
}



/**
* Assemble a movement time deviation array
*
* @param pmean, an array of dimension [P]
* @param ptime, an array of dimension [T, P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [T] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The second argument holds the parameter deviations
* associated with movement rate deviations from the mean for time, ptime.
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year) and acts as the matrix power for multiplying
* movement matrices [X, X].
*/
array[] matrix assemble_movement_time_deviation (
  array[] real pmean,
  array[,] real ptime,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int T = dims(ptime)[1];
  int X = dims(mindex)[1];
  // Initialize movement rates
  matrix[X, X] movement_mean = assemble_movement_mean_rate(
    pmean,
    mindex,
    nterm
  );
  array[T] matrix[X, X] movement_type = assemble_movement_time_rate(
    pmean,
    ptime,
    mindex,
    nterm
  );
  array[T] matrix[X, X] movement_deviation;
  // Populate movement deviation
  for (t in 1:T) {
    for (y in 1:X) {
      for (x in 1:X) {
        movement_deviation[t,x,y] = movement_type[t,x,y] - movement_mean[x,y];
      }
    }
  }
  // Return movement deviation
  return movement_deviation;
}

/**
* Assemble a movement term deviation array
*
* @param pmean, an array of dimension [P]
* @param pterm, an array of dimension [I, P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [I] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The second argument holds the parameter deviations
* associated with movement rate deviations from the mean for term, pterm.
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year) and acts as the matrix power for multiplying
* movement matrices [X, X].
*/
array[] matrix assemble_movement_term_deviation (
  array[] real pmean,
  array[,] real pterm,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int I = dims(pterm)[1];
  int X = dims(mindex)[1];
  // Initialize movement rates
  matrix[X, X] movement_mean = assemble_movement_mean_rate(
    pmean,
    mindex,
    nterm
  );
  array[I] matrix[X, X] movement_type = assemble_movement_term_rate(
    pmean,
    pterm,
    mindex,
    nterm
  );
  array[I] matrix[X, X] movement_deviation;
  // Populate movement deviation
  for (i in 1:I) {
    for (y in 1:X) {
      for (x in 1:X) {
        movement_deviation[i,x,y] = movement_type[i,x,y] - movement_mean[x,y];
      }
    }
  }
  // Return movement deviation
  return movement_deviation;
}

/**
* Assemble a movement size deviation array
*
* @param pmean, an array of dimension [P]
* @param psize, an array of dimension [S, P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [S] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The second argument holds the parameter deviations
* associated with movement rate deviations from the mean for size, psize.
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year) and acts as the matrix power for multiplying
* movement matrices [X, X].
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year) and acts as the matrix power for multiplying
* movement matrices [X, X].
*/
array[] matrix assemble_movement_size_deviation (
  array[] real pmean,
  array[,] real psize,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int S = dims(psize)[1];
  int X = dims(mindex)[1];
  // Initialize movement rates
  matrix[X, X] movement_mean = assemble_movement_mean_rate(
    pmean,
    mindex,
    nterm
  );
  array[S] matrix[X, X] movement_type = assemble_movement_size_rate(
    pmean,
    psize,
    mindex,
    nterm
  );
  array[S] matrix[X, X] movement_deviation;
  // Populate movement deviation
  for (s in 1:S) {
    for (y in 1:X) {
      for (x in 1:X) {
        movement_deviation[s,x,y] = movement_type[s,x,y] - movement_mean[x,y];
      }
    }
  }
  // Return movement deviation
  return movement_deviation;
}

/**
* Assemble Matrices That Indicate Possible Movement at Each Liberty Step
*
* @param mindex, an array of dimension [X, X]
* @param nliberty, an integer giving the maximum steps at liberty
*
* @return an array of dimension [L] holding matrices of dimension [X, X]
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* The argument nliberty is the maximum number of steps at liberty
*/
array[] matrix assemble_movement_possible_transpose (
  array[,] int mindex,
  int nliberty
) {
  // Get dimensions
  int X = dims(mindex)[1];
  int L = nliberty;
  // Initialize values
  array[L] matrix[X, X] movement_possible_transpose;
  matrix[X, X] movement_transpose;
  // Populate movement transpose
  for (x in 1:X) {
    for (y in 1:X) {
      movement_transpose[y, x] = mindex[x, y] * 1.0;
    }
  }
  // Populate movement possible transpose
  movement_possible_transpose[1] = rep_matrix(0.0, X, X);
  movement_possible_transpose[2] = movement_transpose;
  // Iterate higher indexes
  for (l in 3:L) {
    if (min(movement_possible_transpose[l - 1]) > 0.0) {
      movement_possible_transpose[l] = movement_possible_transpose[l - 1];
    } else {
      movement_possible_transpose[l] = movement_possible_transpose[l - 1]
      * movement_transpose;
    }
  }
  // Return movement possible transpose
  return movement_possible_transpose;
}

/**
* Assemble a fishing rate array
*
* @param fstep, an array [N] of vectors [X]
* @param ntime, an integer giving the number of years in the study
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array [T] of vectors [X]
*
* The argument fstep is an array of vectors of stepwise fishing rates.
*
* The arguement ntime is the number of times (years) in the study.
*
* The argument nterm is the number of terms (often quarters) per time
* (usually year).
*/
array[] vector assemble_fishing_rate (
  array[] vector fstep,
  int ntime,
  int nterm
) {
  // Get dimensions
  int N = dims(fstep)[1];
  int X = dims(fstep)[2];
  int T = ntime;
  int I = nterm;
  // Initialize values
  array[T] vector[X] fishing_rate = rep_array(rep_vector(0.0, X), T);
  // Initialize index count
  int count = 1;
  // Populate fishing rate
  for (t in 1:T) {
    for (i in 1:I){
      fishing_rate[t] = fishing_rate[t] + fstep[count];
      count += 1;
    }
  }
  // Return fishing rate
  return fishing_rate;
}

/**
* Assemble a survival step array [N] of diagonal matrices [X, X]
*
* @param fstep, an array [N] of vector [X] stepwise fishing rates
* @param mstep, a vector [X] of stepwise natural mortality rates
* @param ostep, a real ongoing stepwise tag loss rate
*
* @return an array [N] of matrices [X, X]
*/
array[] matrix assemble_survival_step (
  array [] vector fstep,
  vector mstep,
  real ostep
) {
  // Get dimensions
  int N = dims(fstep)[1];
  int X = dims(fstep)[2];
  // Initialize values
  array[N] matrix[X, X] survival_step = rep_array(rep_matrix(0.0, X, X), N);
  // Populate survival step
  for (n in 1:N) {
    survival_step[n] = diag_matrix(exp(-fstep[n]))
    * diag_matrix(exp(-mstep))
    * exp(-ostep);
  }
  // Return survival step
  return survival_step;
}

/**
* Assemble an observation step array [N] of diagonal matrices [X, X]
*
* @param fstep, an array [N] of vector [X] stepwise fishing rates
* @param rstep, a vector [X] of tag reporting steps
*
* @return an array [N] of matrices [X, X]
*/
array[] matrix assemble_observal_step (
  array [] vector fstep,
  vector rstep
) {
  // Get dimensions
  int N = dims(fstep)[1];
  int X = dims(fstep)[2];
  // Initialize values
  array[N] matrix[X, X] observal_step = rep_array(rep_matrix(0.0, X, X), N);
  // Populate reporting step
  for (n in 1:N) {
    observal_step[n] = diag_matrix(1 - exp(-fstep[n])) * diag_matrix(rstep);
  }
  // Return observal step
  return observal_step;
}

/**
* Assemble an fishing term deviation array [T, I] of vectors [X]
*
* @param fstep, an array [N] of vector [X] stepwise fishing rates
* @param frate, an array [N] of vector [X] fishing rates
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array [T, I] of vectors [X]
*/
array[,] vector assemble_fishing_term_deviation(
  array[] vector fstep,
  array[] vector frate,
  int nterm
) {
  // Get dimensions
  int N = dims(fstep)[1];
  int X = dims(fstep)[2];
  int T = dims(frate)[1];
  int I = nterm;
  // Initialize values
  array[T, I] vector[X] fishing_term_deviation = rep_array(rep_vector(0.0,X),T,I);
  // Initialize index count
  int count = 1;
  // Populate fishing term deviation
  for (t in 1:T) {
    for (i in 1:I) {
      fishing_term_deviation[t, i] = frate[t] / (I * 1.0) - fstep[count];
      count += 1;
    }
  }
  // Return fishing term deviation
  return fishing_term_deviation;
}
