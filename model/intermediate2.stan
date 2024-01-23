
data {

	int <lower=0> nData ; //number of data points
	int nObs ; //number of mesurement points

	array[nData] real TIME ;
	array[nData] int EVID ;
	array[nData] int II ;
	array[nData] int SS ;
	array[nData] real AMT ;
	array[nData] int ID ;
	array[nData] real RATE ;
	array[nData] int ADDL ;
	array[nData] int CMT ;

	vector [nObs] cObs ;
	
	int <lower=0> nSub ;
	array[nSub] int ID_list ;
	array[nSub] int starts ; //delimiter for each individual data 
	array[nSub] int ends ; //delimiter for each individual data 

	array[nSub] int nDataSub ; 
	array[nSub] int nObsSub ;
	
	array[nObs] int iObs ; //indices of mesurements (observations)
	
	real <lower=0> CL_pop_meanPrior ;
	real <lower=0> V1_pop_meanPrior ;
	real <lower=0> ka_pop_meanPrior ;
}


transformed data {

	int nCmt = 2 ; // 1 compartment model = 2 cmt
	int nTheta = 3 ; //Number of parameters 
	int K = 3 ; //Number of PK parameters with IIV

	vector[nObs] logcObs = log(cObs) ;

}


parameters {

	real <lower=0, upper=100> CL_pop ;
	real <lower=0, upper=150> V1_pop ;
	real <lower=0, upper=50> ka_pop ; 

	// IIV variablity
	cholesky_factor_corr[K] L  ;
	vector <lower=0> [K] omega ; 
	matrix [K, nSub] etaStd;

	// Residual variability
	real <lower=0> sigma ;

}

transformed parameters {

	vector <lower=0> [K] THETA_pop ; 
	THETA_pop = to_vector({CL_pop, V1_pop, ka_pop}) ;

	//Normalized params
	matrix<lower = 0>[nSub, K] THETA_norm ;
	THETA_norm = (rep_matrix(THETA_pop, nSub) .* exp(diag_pre_multiply(omega, L * etaStd)))' ;

	// Individual params
	//array[nSub] vector <lower=0> [nTheta] THETA_ind ; 
	matrix<lower=0>[nSub, nTheta] THETA_ind ;
	//For now there is no covariate so THETA_ind = THETA_norm 
	THETA_ind = THETA_norm ;

	matrix<lower = 0>[nData, nCmt] X ; 
  	vector<lower = 0>[nData] cHat;
	vector[nObs] cHatObs ;
	
	for (i in 1:nSub) {

		X[starts[i]:ends[i], ] = pmx_solve_onecpt(TIME[starts[i]:ends[i]], AMT[starts[i]:ends[i]], 												RATE[starts[i]:ends[i]], II[starts[i]:ends[i]],
											EVID[starts[i]:ends[i]], CMT[starts[i]:ends[i]],
											ADDL[starts[i]:ends[i]], SS[starts[i]:ends[i]], 
											to_array_1d(THETA_ind[i, ]))';

		cHat[starts[i]:ends[i]] = col(X[starts[i]:ends[i]], 2) ./ THETA_ind[i, 2] ;
	}

	cHatObs = cHat[iObs] ;
	vector[nObs] logcHatObs = log(cHatObs) ;
	
	
}


model {

	// PK parameters
	CL_pop ~ normal(CL_pop_meanPrior, 5) ;
    V1_pop ~ normal(V1_pop_meanPrior, 5);
    ka_pop ~ normal(ka_pop_meanPrior, 0.5);

	// IIV variability
	L ~ lkj_corr_cholesky(1);

	to_vector(etaStd) ~ std_normal() ;
  	omega ~ normal(0, 0.1);

	// Residual variability
    sigma ~ normal(0, 1);

	// Likelihood
    fmax(machine_precision(), logcObs) ~ normal(fmax(machine_precision(), logcHatObs), sigma); 

}


generated quantities {

	vector[nObs] log_lik ;
	matrix[K, nSub] etaStdPred ;
	matrix<lower = 0>[nSub, K] THETAPred_norm ;
	matrix<lower = 0>[nSub, K] THETAPred_ind ;

	matrix[nData, 2] XPred ;
	vector[nData] cHatPred ;
	vector[nData] cPred_pop ; // population predictions
	vector[nData] cPred_ind ; // individual predictions

	for(i in 1:nSub) {
    	for(j in 1:K) {
      		etaStdPred[j, i] = normal_rng(0, 1);
    	}
	}

	THETAPred_norm = (rep_matrix(THETA_pop, nSub) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';

	//For now there is no covariate so THETA_ind = THETA_norm 
	THETAPred_ind = THETAPred_norm ;

	for (i in 1:nSub) {

		XPred[starts[i]:ends[i], ] = pmx_solve_onecpt(TIME[starts[i]:ends[i]], AMT[starts[i]:ends[i]],
														RATE[starts[i]:ends[i]], II[starts[i]:ends[i]],
														EVID[starts[i]:ends[i]], CMT[starts[i]:ends[i]],
														ADDL[starts[i]:ends[i]], SS[starts[i]:ends[i]], 
														to_array_1d(THETAPred_ind[i, ]))';

		cHatPred[starts[i]:ends[i]] = col(XPred[starts[i]:ends[i]], 2) ./ THETAPred_ind[i, 2] ;                
	}
	

	for (t in 1:nData) {

		if (cHat[t] < 0.0001) {cPred_ind[t] = 0 ;}
		else {cPred_ind[t] = lognormal_rng(log(cHat[t]), sigma);}

		if (cHatPred[t] < 0.0001) {cPred_pop[t] = 0 ;}
		else {cPred_pop[t] = lognormal_rng(log(cHatPred[t]), sigma) ;}
	}

	for (i_obs in 1:nObs) {
	
		log_lik[i_obs] = normal_lpdf(cObs[i_obs] | cHatObs[i_obs], sigma);
	}


}

