data {
  int<lower = 0> n;
  int<lower = 0> p1;
  int<lower = 0> p2;
  int<lower = 0> p3;
  int<lower = 0> p4;
  int<lower = 0> p5;
  matrix[p1,n] y1;
  matrix[p2,n] y2;
  matrix[p3,n] y3;
  matrix[p4,n] y4;
  matrix[p5,n] y5;
  real mu5;
}

parameters {
  real pi1;real pi2;real pi3;real pi4;
  simplex[p5] m1[p1];
  simplex[p5] m2[p2];
  simplex[p5] m3[p3];
  simplex[p5] m4[p4];
}
transformed parameters {
  matrix[p5,p1] mt1;matrix[p5,p2] mt2;matrix[p5,p3] mt3;matrix[p5,p4] mt4;
  matrix[p5,n] mx1;matrix[p5,n] mx2;matrix[p5,n] mx3;matrix[p5,n] mx4;
  matrix[p5,n] mx;
  for(i in 1:p1){
    mt1[,i] = m1[i];
  }
  for(i in 1:p2){
    mt2[,i] = m2[i];
  }
  for(i in 1:p3){
    mt3[,i] = m3[i];
  }
  for(i in 1:p4){
    mt4[,i] = m4[i];
  }
  mx1 = pi1*mt1*y1;
  mx2 = pi2*mt2*y2;
  mx3 = pi3*mt3*y3;
  mx4 = pi4*mt4*y4;
  mx = mx1+mx2+mx3+mx4;
}
model {
  pi1~gamma(1,1);
  pi2~gamma(1,1);
  pi3~gamma(1,1);
  pi4~gamma(1,1);
  for(j in 1:p1) m1[j]~dirichlet(rep_vector(1.0/p5,p5));
  for(j in 1:p2) m2[j]~dirichlet(rep_vector(1.0/p5,p5));
  for(j in 1:p3) m3[j]~dirichlet(rep_vector(1.0/p5,p5));
  for(j in 1:p4) m4[j]~dirichlet(rep_vector(1.0/p5,p5));
  for(i in 1:n){
    target += dirichlet_lpdf( y5[,i] |(mx[,i]+mu5));
  }
}
generated quantities {
  vector[n] log_lik;
  for (i in 1:n) log_lik[i] = dirichlet_lpdf(y5[,i] |(mx[,i]+mu5));
}
