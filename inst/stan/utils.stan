// Movement slab
real[,,] assemble_movement_slab (
  real[] pmean,
  real[,] ptime,
  real[,] pterm,
  real[,] psize
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
  array[N, S, P] real parameter_slab = assemble_parameter_slab(
    pmean,
    ptime,
    pterm,
    psize
  );
  array[N, S, P + 1] real movement_slab = rep_array(0.0, N, S, P + 1);
  // Check dimensions
  if (P < 1) {reject("P must not be < 1; found P = ", P);}
  if (T < 1) {reject("T must not be < 1; found T = ", T);}
  if (I < 1) {reject("I must not be < 1; found I = ", I);}
  if (S < 1) {reject("S must not be < 1; found S = ", S);}
  // Populate movement slab
  n = 1;
  for (t in 1:T) {
    for (i in 1:I) {
      for (s in 1:S) {
        movement_slab[n, s] = inverse_multi_logit(parameter_slab[n, s]);
      }
      n += 1;
    }
  }
  // Return movement slab
  return movement_slab;
}

// Parameter slab
real[,,] assemble_parameter_slab (
  real[] pmean,
  real[,] ptime,
  real[,] pterm,
  real[,] psize
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

// Inverse multi-logit
real[] inverse_multi_logit (real[] arg) {
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
