data {
  int nt;
  array[nt] real t;
  array[nt] real yobs;
}

transformed data{
  int n;
  n=1;
}

parameters {
  array[n+1] real p ;
  real<lower=0> sigma;
}

transformed parameters {
	real beta0 = p[1] ;
	real beta1 = p[2] ;
}

model {
  
  array[nt] real y ;
  for (i in 1:nt) {
      y[i]=0;
          for (j in 1:n+1){
               y[i] = y[i] + p[j]*pow(t[i],j-1); 
              }
   }
           

  
  yobs ~ normal(y, sigma);
}

generated quantities {

  vector[nt] log_lik;
  array[nt] real y;


  for (i in 1:nt) {
      y[i]=0;
          for (j in 1:n+1){
               y[i] = y[i] + p[j]*pow(t[i],j-1); 
              }
   }

  for (i in 1:nt){
  log_lik[i] = normal_lpdf(yobs[i] | y[i], sigma);
  }
  
  
  array[nt] real y_rep ;
 
  for (i in 1:nt){
  	y_rep[i] = normal_rng(y[i], sigma);
  }
  
}
