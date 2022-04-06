/**
* Inverse logit
*
* @param arg, a real value
*
* @return a value in [0, 1]
*/
real inverse_logit(real arg) {
  real value;
  if (is_inf(arg)) {
    if (arg > 0) {
      value = 1;
    } else {
      value = 0;
    }
  } else {
    value = exp(arg) / (1 + exp(arg));
  }
  return value;
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
  array[P + 1] real value;
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
        movement_slab[n, s, m[1]:m[2]] = inverse_multi_logit(
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
* Assemble a movement mean rate matrix
*
* @param pmean, an array of dimension [P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return a matrix of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates.
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*/
matrix assemble_movement_mean_rate (
  array[] real pmean,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int X = dims(mindex)[1];
  // Initialize arrays
  array[1, 1] matrix[X, X] movement_step = assemble_movement_step(
    pmean,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    mindex
  );
  // Initialize matrix
  matrix[X, X] movement_rate = matrix_power(movement_step[1, 1], nterm);
  // Return movement rate
  return movement_rate;
}

/**
* Assemble a movement time rate array
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
*/
array[] matrix assemble_movement_time_rate (
  array[] real pmean,
  array[,] real ptime,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int T = dims(ptime)[1];
  int X = dims(mindex)[1];
  // Initialize arrays
  array[T, 1] matrix[X, X] movement_step = assemble_movement_step(
    pmean,
    ptime,
    rep_array(0.0, 1, P),
    rep_array(0.0, 1, P),
    mindex
  );
  array[T] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0, X, X), T);
  // Populate movement rates
  for (t in 1:T) {
    movement_rate[t] = matrix_power(movement_step[t, 1], nterm);
  }
  // Return movement rate
  return movement_rate;
}

/**
* Assemble a movement term rate array
*
* @param pmean, an array of dimension [P]
* @param pterm, an array of dimension [I, P]
* @param mindex, an array of dimension [X, X]
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
*/
array[] matrix assemble_movement_term_rate (
  array[] real pmean,
  array[,] real pterm,
  array[,] int mindex
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
    movement_rate[i] = matrix_power(movement_step[i, 1], 1);
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
  // Set dimensions
  int N = T * I;
  int Z = N * S;
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
* Assemble a movement (full) rate array
*
* @param pmean, an array of dimension [P]
* @param ptime, an array of dimension [T, P]
* @param pterm, an array of dimension [I, P]
* @param psize, an array of dimension [S, P]
* @param mindex, an array of dimension [X, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [N, S] holding matrices of dimension [X, X]
*
* The first argument pmean holds the parameters associated with the mean
* movement rates. The next three arguments hold the parameter deviations
* associated with movement rate deviations from the mean for time (ptime),
* term (pterm), and size (psize).
*
* The argument mindex is a square integer array of zeros and ones indicating
* that movement is permitted (one) or not permitted (zero) between a given
* source region (row) and destination region (column) in one step.
*
* Note that N = T * I, which is the number of study steps and one more than
* the number of released steps, N - 1.
*/
array[,] matrix assemble_movement_full_rate (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] int mindex,
  int nterm
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  int X = dims(mindex)[1];
  // Set dimensions
  int N = T * I;
  // Initialize arrays
  array[N, S] matrix[X, X] movement_step = assemble_movement_step(
    pmean,
    ptime,
    pterm,
    psize,
    mindex
  );
  array[N, S] matrix[X, X] movement_rate = rep_array(rep_matrix(0.0,X,X),N,S);
  // Populate movement rates
  for (n in 1:N) {
    for (s in 1:S) {
      movement_rate[n, s] = matrix_power(movement_step[n, s], nterm);
    }
  }
  // Return movement rate
  return movement_rate;
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
array[] matrix assemble_movement_deviation (
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
*/
array[] matrix assemble_movement_term_deviation (
  array[] real pmean,
  array[,] real pterm,
  array[,] int mindex
) {
  // Get dimensions
  int P = dims(pmean)[1];
  int I = dims(pterm)[1];
  int X = dims(mindex)[1];
  // Initialize movement rates
  matrix[X, X] movement_mean = assemble_movement_mean_rate(
    pmean,
    mindex,
    1 // nterm = 1 for movement term rates
  );
  array[I] matrix[X, X] movement_type = assemble_movement_term_rate(
    pmean,
    pterm,
    mindex
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
*/
array[] matrix assemble_movement_deviation (
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
* Assemble a fishing parameter array 'slab'
*
* @param pmean, an array of dimension [X]
* @param ptime, an array of dimension [T, X]
* @param pterm, an array of dimension [I, X]
* @param psize, an array of dimension [S, X]
*
* @return an array of dimension [N, S, X]
*
* The first argument fmean holds the parameters associated with the mean
* fishing rates. The remaining arguments hold the parameter deviations
* associated with fishing rate deviations from the mean for time (ptime),
* term (pterm), and size (psize).
*
* Parameter slabs associated with fishing rate means only,
* or fishing rates plus certain types of deviations can be
* assembled by passing rep_array(0.0, 1, X) for arguments corresponging
* to the deviations to exclude.
*
* Note that N = T * I, which is the number of study steps and one more than
* the number of released steps, N - 1.
*/
array[,,] real assemble_fishing_parameter_slab (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize
) {
  // Get dimensions
  int X = dims(pmean)[1];
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
* Assemble a fishing array at timescale 'step'
*
* @param pmean, an array of dimension [X]
* @param ptime, an array of dimension [T, X]
* @param pterm, an array of dimension [I, X]
* @param psize, an array of dimension [S, X]
*
* @return an array of dimension [N, S, X]
*
* The first argument pmean holds the parameters associated with the mean
* fishing rates. The remaining arguments hold the parameter deviations
* associated with fishing rate deviations from the mean for time (ptime),
* term (pterm), and size (psize).
*
* Fishing associated with fishing rate means only, or fishing
* rates plus certain types of deviations can be assembled by setting the
* arguments for deviations to exclude to rep_array(0.0, 1, X).
*
* Note that N = T * I, which is the number of study steps and one more than
* the number of released steps, N - 1.
*/
array[,,] real assemble_fishing_step (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize
) {
  // Get dimensions
  int X = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  // Set dimensions
  int N = T * I;
  // Initialize arrays
  array[N, S, X] real parameter_slab = assemble_fishing_parameter_slab(
    pmean,
    ptime,
    pterm,
    psize
  );
  array[N, S, X] real fishing_step = rep_array(0.0, N, S, X);
  // Populate fishing step
  for (n in 1:N) {
    for (s in 1:S) {
      for (x in 1:X) {
        // Fix this:
        // Maps values in parameter_slab (-Inf, Inf) to values in (-Inf, 0)
        // for fishing_step by log(exp(a) / (1 + exp(a))) = a - log(1 + exp(a))
        fishing_step[n, s, x] = parameter_slab[n, s, x]
        - log(1 + exp(parameter_slab[n, s, x]));
      }
    }
  }
  // Return fishing step
  return fishing_step;
}





// Old below here --------------------------------------------------------------


/**
* Assemble a fishing rate array
*
* @param pmean, an array of dimension [X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [X]
*
* The first argument pmean holds the parameters associated with the mean
* fishing mortality rates.
*
*/
array[] real assemble_fishing_rate (
  array[] real pmean,
  int nterm
) {
  // Get dimensions
  int X = dims(pmean)[1];
  // Initialize arrays
  array[1, 1, X] real fishing_step = assemble_fishing_step(
    pmean,
    rep_array(0.0, 1, X),
    rep_array(0.0, 1, X),
    rep_array(0.0, 1, X)
  );
  // Initialize array
  array[X] real fishing_rate = pow(fishing_step[1, 1], nterm);
  // Return fishing rate
  return fishing_rate;
}

/**
* Assemble a fishing rate array
*
* @param pmean, an array of dimension [X]
* @param ptime, an array of dimension [T, X]
* @param pterm, an array of dimension [I, X]
* @param psize, an array of dimension [S, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [Z, X] holding fishing mortality rate
*
* The first argument pmean holds the parameters associated with the mean
* fishing rates. The next three arguments hold the parameter deviations
* associated with fishing rate deviations from the mean for time (ptime),
* term (pterm), and size (psize). Two of these three should be placeholders
* rep_array(0.0, 1, X).
*
*/
array[,] real assemble_fishing_rate (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  int nterm
) {
  // Get dimensions
  int X = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  // Set dimensions
  int N = T * I;
  int Z = N * S;
  // Initialize arrays
  array[N, S, X] real fishing_step = assemble_fishing_step(
    pmean,
    ptime,
    pterm,
    psize
  );
  array[Z, X] real fishing_rate = rep_array(0.0, Z, X);
  // Initialize count
  int z;
  // Check dimensions
  if (X < 1) {reject("X must not be < 1; found X = ", X);}
  if (T < 1) {reject("T must not be < 1; found T = ", T);}
  if (I < 1) {reject("I must not be < 1; found I = ", I);}
  if (S < 1) {reject("S must not be < 1; found S = ", S);}
  // if (Z > max(T, max(I, S))) {reject("T * I * S must not be > max(T, I, S)");}
  // Populate fishing rates
  z = 1;
  for (n in 1:N) {
    for (s in 1:S) {
      fishing_rate[z] = pow(fishing_step[N, S], nterm);
      z += 1;
    }
  }
  // Return fishing rate
  return fishing_rate;
}

/**
* Assemble a fishing rate deviation array
*
* @param pmean, an array of dimension [X]
* @param ptime, an array of dimension [T, X]
* @param pterm, an array of dimension [I, X]
* @param psize, an array of dimension [S, X]
* @param nterm, an integer giving the number of seasons (terms) per year (times)
*
* @return an array of dimension [Z, X]
*
* The first argument pmean holds the parameters associated with the mean
* fishing rates. The next three arguments hold the parameter deviations
* associated with fishing rate deviations from the mean for time (ptime),
* term (pterm), and size (psize). Two of these three should be placeholders
* rep_array(0.0, 1, X).
*
*/
array[,] real assemble_fishing_deviation (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  int nterm
) {
  // Get dimensions
  int X = dims(pmean)[1];
  int T = dims(ptime)[1];
  int I = dims(pterm)[1];
  int S = dims(psize)[1];
  // Set dimensions
  int Z = T * I * S;
  // Initialize movement rates
  array[X] real fishing_mean = assemble_fishing_rate(
    pmean,
    nterm
  );
  array[Z, X] real fishing_type = assemble_fishing_rate(
    pmean,
    ptime,
    pterm,
    psize,
    nterm
  );
  array[Z, X] real fishing_deviation;
  // Check dimensions
  if (X < 1) {reject("X must not be < 1; found X = ", X);}
  if (T < 1) {reject("T must not be < 1; found T = ", T);}
  if (I < 1) {reject("I must not be < 1; found I = ", I);}
  if (S < 1) {reject("S must not be < 1; found S = ", S);}
  // if (Z > max(T, max(I, S))) {reject("T * I * S must not be > max(T, I, S)");}
  // Populate movement deviation
  for (z in 1:Z) {
    for (x in 1:X) {
      fishing_deviation[z, x] = fishing_type[z, x] - fishing_mean[x];
    }
  }
  // Return movement deviation
  return fishing_deviation;
}
