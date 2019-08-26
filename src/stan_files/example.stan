//
//  example stan
//

data {
  //existing data
  int<lower = 1>  N;
  real            Y[N];
}


parameters {
  real          theta;
}

model {
  //prior
  theta ~ normal(0, 1000);

  //likelihood
  Y ~ normal(theta, 1);
}
