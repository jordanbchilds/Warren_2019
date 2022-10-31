data {
  int D; //number of dimensions
  int K; //number of mixtures
  
  vector[D] mean_mu1;
  vector[D] mean_mu2;
  cov_matrix[D] var_mu1;
  cov_matrix[D] var_mu2;
  int<lower=D> n_1;
  int<lower=D> n_2;
  cov_matrix[D] S_1;
  cov_matrix[D] S_2;
  
  real<lower=0> g_alpha;
  real<lower=0> h_alpha;
  real<lower=0> g_beta;
  real<lower=0> h_beta;
}
parameters {
  vector[D] mu1_;
  vector[D] mu2_;
  cov_matrix[D] V1_; //precision of component one
  cov_matrix[D] V2_;
  real<lower=0> alpha_;
  real<lower=0> beta_;
  real<lower=0, upper=1> probctrl_; //'like-ctrl' patient proportion
}
model {
  mu1_ ~ multi_normal(mean_mu1, var_mu1);
  mu2_ ~ multi_normal(mean_mu2, var_mu2);
  V1_ ~ inv_wishart(n_1, S_1);
  V2_ ~ inv_wishart(n_2, S_2);
  alpha_ ~ gamma(g_alpha, h_alpha);
  beta_ ~ gamma(g_beta, h_beta);
  probctrl_ ~ beta(alpha_, beta_);
}
generated quantities{
  matrix[D,K] comp;
  vector[D] pred;
  int<lower=0, upper=1> z;
  vector[D] mu1;
  vector[D] mu2;
  cov_matrix[D] V1; //precision of component one
  cov_matrix[D] V2;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0, upper=1> probctrl;
  
  mu1 = multi_normal_rng(mean_mu1, var_mu1);
  mu2 = multi_normal_rng(mean_mu2, var_mu2);
  V1 = inv_wishart_rng(n_1, S_1);
  V2 = inv_wishart_rng(n_2, S_2);
  alpha = gamma_rng(g_alpha, h_alpha);
  beta = gamma_rng(g_beta, h_beta);
  probctrl = beta_rng(alpha, beta);
  z = bernoulli_rng(probctrl);
  
  if(z){
    pred = multi_normal_rng(mu1, V1);
  } else {
    pred = multi_normal_rng(mu2, V2);
  }
  comp[,1] = multi_normal_rng(mu1, V1);
  comp[,2] = multi_normal_rng(mu2, V2);
}

