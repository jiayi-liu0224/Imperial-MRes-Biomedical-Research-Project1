//hierarchical_model

data {
  int<lower = 0> N_t; // number of traits
  int<lower = 0> N_c; // number of cells
  real<lower=-1,upper=1> z0; //the root value
  matrix[N_c,N_c] C;  // variance-covariance matrix of the tree
  matrix<lower=-1,upper=1>[N_c,N_t] Y;   //trait value
}

transformed data {
  vector[N_c] mu_vector;
  mu_vector=rep_vector(z0, N_c);
}
parameters {
  vector[N_t] eta; 
  real<lower=0> mu;
  real<lower=0> tau;
}

transformed parameters {
  vector<lower=0>[N_t] sig2; //the mutation rate
  sig2 = mu + tau*eta;
}


model {
  // Hyperpriors
  mu~normal(0.001,0.05); 
  tau~normal(0,0.05);

  for (n in 1:N_t){
  //Set the prior
  eta[n]~normal(0,1);
  }
  for (n in 1:N_t){
  Y[:,n] ~ multi_normal(mu_vector,sig2[n]*C);
  }

}
