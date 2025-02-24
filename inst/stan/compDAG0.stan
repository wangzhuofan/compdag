data {
  int<lower = 0> n;
  int<lower = 0> p1;
  matrix[p1,n] y1;
}

parameters {
  real mu1;
}
model {
  mu1~gamma(1,1);
  for(i in 1:n){
    target += dirichlet_lpdf( y1[,i] |rep_vector(mu1,p1));
  }
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = dirichlet_lpdf(y1[,i] |rep_vector(mu1,p1));
}
