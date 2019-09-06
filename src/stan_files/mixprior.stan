// Mixture prior for platform design
// To take into account time effect
functions {
  real icdf(real p, int nsmp, vector smp) {
    real rst;
    int  low = 1;
    int  up  = nsmp;
    int  mid; 

    while (up > low + 1) {
      mid = (low + up)/2;
      if (smp[mid] > p) {
        up = mid;
      } else {
        low = mid;
      }
    }

    /* while (low < nsmp) { */
    /*   up = min(low + 500, nsmp); */
    /*   if (smp[up] < p) { */
    /*     low = up;  */
    /*   } else { */
    /*     break; */
    /*   } */
    /* } */

    /* l = low; */
    /* for (i in low:up) { */
    /*   if (smp[i] > p) */
    /*     break; */
    /*   l = l+1; */
    /* } */

    rst = (1.0*low)/nsmp;
    return(rst);
  }

  real icdf2(real p, int nsmp, vector smp) {
    real rst;
    real   l;
    vector[nsmp + 1] vv;
    vv  = append_row(p, smp);
    l   = rank(vv, 1);
    rst = l/nsmp;
    return(rst);
  }

}

data {
  int<lower = 1>  N; //total number of patients
  int<lower = 0>  Y; //number of responsders

  int<lower = 1>  NSMP;  //number of prior samples 
  vector[NSMP]    SMP;   //prior samples 
}

parameters {
  real<lower = 0, upper = 1> p;
}

transformed parameters {
  real theta = icdf(p, NSMP, SMP);
}

model {
  // log-likelihood
  target += uniform_lpdf(p | 0, 1);
  target += binomial_lpmf(Y | N, theta); 
}

generated quantities {
}
