data {
  int<lower = 0> n;
  int<lower = 0> p1;
  int<lower = 0> p2;
  matrix[p1,n] y1;
  matrix[p2,n] y2;
  real mu2;
}

parameters {
  real pi1;
  simplex[p2] m1[p1];
}
transformed parameters {
  matrix[p2,p1] mt1;
  matrix[p2,n] mx1;
  matrix[p2,n] mx;
  for(i in 1:p1){
    mt1[,i] = m1[i];
  }
  mx1 = pi1*mt1*y1;
  mx = mx1;
}
model {
  pi1~gamma(1,1);
  for(j in 1:p1) m1[j]~dirichlet(rep_vector(1.0/p2,p2));
  for(i in 1:n){
    target += dirichlet_lpdf( y2[,i] |(mx[,i]+mu2));
  }
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = dirichlet_lpdf(y2[,i] |(mx[,i]+mu2));
}
