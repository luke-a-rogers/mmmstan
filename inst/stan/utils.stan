/**
* Forward declarations. Allows function definitions in alphabetical-ish order.
*/
array[,,] real assemble_parameter_slab (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize
);
array[,,] real assemble_movement_slab (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] real mindex
);
array[,] matrix assemble_movement_step (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] real mindex
);
array[,] matrix assemble_movement_rate (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] real mindex
);
array[] real inverse_multi_logit (array[] real arg);

/**
* Assemble a parameter array 'slab'
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
* Parameter slabs associated with movement rate means only, or movement
* rates plus certain types of deviations can be assembled by setting the
* arguments for deviations to exclude to rep(0.0, 1, P). Note that N = T * I.
*/
array[,,] real assemble_parameter_slab (
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
  // Check dimensions
  if (P < 1) {reject("P must not be < 1; found P = ", P);}
  if (T < 1) {reject("T must not be < 1; found T = ", T);}
  if (I < 1) {reject("I must not be < 1; found I = ", I);}
  if (S < 1) {reject("S must not be < 1; found S = ", S);}
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
* Movement slabs associated with movement rate means only, or movement
* rates plus certain types of deviations can be assembled by setting the
* arguments for deviations to exclude to rep(0.0, 1, P). Note that N = T * I.
*/
array[,,] real assemble_movement_slab (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] real mindex
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
  array[N, S, P] real parameter_slab = assemble_parameter_slab(
    pmean,
    ptime,
    pterm,
    psize
  );
  array[N, S, M] real movement_slab = rep_array(0.0, N, S, M);
  // Check dimensions
  if (P < 1) {reject("P must not be < 1; found P = ", P);}
  if (X < 1) {reject("X must not be < 1; found X = ", X);}
  if (N < 1) {reject("N must not be < 1; found N = ", N);}
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
* Movement associated with movement rate means only, or movement
* rates plus certain types of deviations can be assembled by setting the
* arguments for deviations to exclude to rep(0.0, 1, P). Note that N = T * I.
*/
array[,] matrix assemble_movement_step (
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] real mindex
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
  // Check dimensions
  if (P < 1) {reject("P must not be < 1; found P = ", P);}
  if (X < 1) {reject("X must not be < 1; found X = ", X);}
  if (N < 1) {reject("N must not be < 1; found N = ", N);}
  // Populate movement step
  for (n in 1:N) {
    for (s in 1:S) {
      m = 1;
      transpose_step = rep_matrix(0.0, X, X);
      for (x in 1:X) {
        for (y in 1:X) {
          if (mindex[x, y] == 1) {
            transpose_step[y, x] = movement_slab[n, s, m];
            m += 1;
          }
        }
      }
      movement_step[n, s] = transpose_step';
    }
  }
  // Return movement step
  return movement_step;
}

// Check above and check more thoroughly below here

/**
* Assemble a movement rate array
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
* Movement rates associated with movement rate means only, or movement
* rates plus certain types of deviations can be assembled by setting the
* arguments for deviations to exclude to rep(0.0, 1, P). Note that N = T * I.
*/
array[,] matrix assemble_movement_rate ( // May need to separate into two (or more) functions: one or two for rate, one for deviations
  array[] real pmean,
  array[,] real ptime,
  array[,] real pterm,
  array[,] real psize,
  array[,] real mindex
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
      movement_rate[n, s] = matrix_power(movement_step[n, s], I); // TODO: may need to replace I
    }
  }
  // Return movement rate
  return movement_rate;
}

/**
* Inverse multi-logit
*
* @param arg, an array of one dimension
*
* @return an array of one dimension that sums to one
*/
array[] real inverse_multi_logit (array[] real arg) {
  // Get dimensions
  int P = size(arg);
  // Instantiate sums
  real sum_exp_arg = 0;
  real sum_values = 0;
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
    sum_values += value[p];
  }
  value[P + 1] = 1 - sum_values;
  // Return value
  return value;
}
