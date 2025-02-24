data {
  int<lower = 0> n;
  int<lower = 0> p1;
  int<lower = 0> p2;
  int<lower = 0> p3;
  matrix[p1,n] y1;
  matrix[p2,n] y2;
  matrix[p3,n] y3;
  real mu3;
}

parameters {
  real pi1;real pi2;
  simplex[p3] m1[p1];
  simplex[p3] m2[p2];
}
transformed parameters {
  matrix[p3,p1] mt1;matrix[p3,p2] mt2;
  matrix[p3,n] mx1;matrix[p3,n] mx2;
  matrix[p3,n] mx;
  for(i in 1:p1){
    mt1[,i] = m1[i];
  }
  for(i in 1:p2){
    mt2[,i] = m2[i];
  }
  mx1 = pi1*mt1*y1;
  mx2 = pi2*mt2*y2;
  mx = mx1+mx2;
}
model {
  pi1~gamma(1,1);
  pi2~gamma(1,1);
  for(j in 1:p1) m1[j]~dirichlet(rep_vector(1.0/p3,p3));
  for(j in 1:p2) m2[j]~dirichlet(rep_vector(1.0/p3,p3));
  for(i in 1:n){
    target += dirichlet_lpdf( y3[,i] |(mx[,i]+mu3));
  }
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = dirichlet_lpdf(y3[,i] |(mx[,i]+mu3));
}
