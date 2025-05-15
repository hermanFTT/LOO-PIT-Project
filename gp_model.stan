

functions {
vector diagSPD_EQ(real alpha, real rho, real L, int M) {
  return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
}
vector diagSPD_Matern32(real alpha, real rho, real L, int M) {
   return 2*alpha * (sqrt(3)/rho)^1.5 * inv((sqrt(3)/rho)^2 + ((pi()/2/L) * linspaced_vector(M, 1, M))^2);
}
vector diagSPD_periodic(real alpha, real rho, int M) {
  real a = 1/rho^2;
  vector[M] q = exp(log(alpha) + 0.5 * (log(2) - a + to_vector(log_modified_bessel_first_kind(linspaced_int_array(M, 1, M), a))));
  return append_row(q,q);
}
matrix PHI(int N, int M, real L, vector x) {
  return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
}
matrix PHI_periodic(int N, int M, real w0, vector x) {
  matrix[N,M] mw0x = diag_post_multiply(rep_matrix(w0*x, M), linspaced_vector(M, 1, M));
  return append_col(cos(mw0x), sin(mw0x));
}

}
data {
  int<lower=1> N;      // number of observations
  vector[N] x;         // univariate covariate
  vector[N] y;         // target variable
        
  real<lower=0> c_f;   // factor c to determine the boundary value L
  int<lower=1> M_f;    // number of basis functions
  real<lower=0> lengthscale_f; // lengthscale of f
  real<lower=0> sigma_f;       // scale of f
}
transformed data {
  // Normalize data
  real xmean = mean(x);
  real ymean = mean(y);
  real xsd = sd(x);
  real ysd = sd(y);
  vector[N] xn = (x - xmean)/xsd;
  vector[N] yn = (y - ymean)/ysd;
  // Basis functions for f
  real L_f = c_f*max(xn);
  matrix[N,M_f] PHI_f = PHI(N, M_f, L_f, xn);
}
parameters {
  real intercept_f;
  vector[M_f] beta_f;          // the basis functions coefficients
  real<lower=0> sigman;         // noise sigma
}
model {
  // spectral densities for f and g
  vector[M_f] diagSPD_f = diagSPD_EQ(sigma_f, lengthscale_f, L_f, M_f);
  // priors
  intercept_f ~ normal(0, 0.1);
  beta_f ~ normal(0, 1);
  lengthscale_f ~ normal(0, 1);
  sigma_f ~ normal(0, 1);
  sigman ~ normal(0, 1);
  // model
  yn ~ normal(intercept_f + PHI_f * (diagSPD_f .* beta_f), sigman);
}
generated quantities {
  vector[N] f;
  vector[N] log_lik;
  real sigma = sigman*ysd;
  {
    // spectral densities
    vector[M_f] diagSPD_f = diagSPD_EQ(sigma_f, lengthscale_f, L_f, M_f);
    // function scaled back to the original scale
    f = (intercept_f + PHI_f * (diagSPD_f .* beta_f))*ysd + ymean;
    for (n in 1:N) {
      log_lik[n] = normal_lpdf(yn[n] | intercept_f + PHI_f[n] * (diagSPD_f .* beta_f),
          sigman);
    }
  }
}
