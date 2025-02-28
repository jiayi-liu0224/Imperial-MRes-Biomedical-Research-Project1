//Non_hierarchical_model

data {
  int<lower = 0> N_t; // number of traits
  int<lower = 0> N_c; // number of cells
  real z0; //the root value
  cov_matrix[N_c] C;  // variance-covariance matrix of the tree
  matrix[N_c,N_t] Y;   //trait value
}


parameters {
  vector<lower=0>[N_t] sig2; //the mutation rate
}

model {
  //Set the prior
  sig2~normal(0.001,0.1);
  for (n in 1:N_t){
  Y[:,n] ~ multi_normal(rep_vector(z0,N_c),sig2[n]*C);
  }
  
}
