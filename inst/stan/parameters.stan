// Beta parameters
real[] beta_pars (real mu, real sigma) {
  real pars[2] = rep_array(0, 2);
  real nu = (mu * (1 - mu) / (sigma * sigma)) - 1;
  real alpha = mu * nu;
  real beta = (1 - mu) * nu;
  pars[1] = alpha;
  pars[2] = beta;
  return pars;
}

// Gamma parameters
real[] gamma_pars (real mu, real sigma) {
  real pars[2] = rep_array(0, 2);
  real alpha = mu^2 / sigma^2;
  real beta = mu / sigma^2;
  pars[1] = alpha;
  pars[2] = beta;
  return pars;
}
