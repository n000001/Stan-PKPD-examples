functions {

	vector expMatSolve(real t0, real t, vector init, 
						real CL, real Q, real V1, real V2, real ka, real ke0) {
 
		matrix[4, 4] K ;
		K = rep_matrix(0, 4, 4) ; 
		K[1, 1] = - ka ;
		K[2, 1] = ka ;
		K[2, 2] = - (CL/V1 + Q/V1) ;
		K[2, 3] = Q/V2 ;
		K[3, 2] = Q/V1 ;
		K[3, 3] = - Q/V2 ;
		K[4, 2] = ke0 ;
		K[4, 4] = -ke0 ;
		
		
		vector[4] x ;
		x = matrix_exp((t - t0) * K) * init ;
	   
	   	return x ;}


	matrix solvetwoCptModel(array[] real time, 
							array[] real amt, array[] int cmt, array[] int evid, 
							real CL, real Q, real V1, real V2, real ka, real ke0) {

  		int nt = size(time) ;
		vector[4] init = rep_vector(0, 4) ;
		matrix[nt, 4] result;

		for (i in 1:nt){
		
			init = expMatSolve(time[max(1, i - 1)], time[i], init, CL, Q, V1, V2, ka, ke0) ;

			if (evid[i] == 1) init[cmt[i]] += amt[i] ; // dose

			result[i, ] = init' ;
		}

		return result ;}

}


data{

	//General data
	int<lower = 1> nt;
	int< lower = 1 > N ;
	array[N] real time ;
	array[N] int cmt ;
	array[N] int evid ;
	array[N] int addl ;
	array[N] int ii ;
	array[N] int ss ;
	array[N] real amt ;
	array[N] real dose ;

	//Observations
	int <lower = 1> nObsPK ;
	array[nObsPK] int <lower = 1> iObsPK ;
	vector <lower=0> [nObsPK] cObs ;
	vector <lower=0> [nObsPK] effObs ;

	//Individuals related data
	int <lower=1> nSubjects;
	vector <lower=0> [nSubjects] weight ;
	array[nSubjects] int <lower=0> start;
	array[nSubjects] int <lower=0> end;

}


transformed data{

	int nTheta = 7 ; //number of params
	int nCmt = 4 ; //3components model
	int nIIV = 5 ; //number of individual related params
	int K = nIIV ;


	// Chol. factor param 

	int eta = 1 ; //choose value
	

	// IIV params variability prior params

	vector[K] mu_omega = rep_vector(0.25, K) ; //choose value
	matrix[K,K] cov_omega = diag_matrix(rep_vector(0.25, K)) ;
	

	// Adapt data to multi-dosing schedule 

	array[N] real new_amt = dose ;
	array[N] int new_evid = evid ;
	array[N] int new_cmt = cmt ;

	real <lower=0> eps_dose = 0.05 ; //tolerance for dose administration time
	// there is a problem in the case where time vector does not have values close to dosing times
	// then it will be necessary to create a new time vector 

	for (t in 1:N) {

		if ((evid[t] == 1) && (dose[t]*addl[t]*ii[t] != 0)){

			int t0 = t ;
			int T = t ;
			int ad = 0 ;
			while ((T < N-1) && (ad < addl[t])) {
				T += 1 ;
				if (abs( (time[T]-time[t0]) - ii[t] ) <= eps_dose) {
					new_amt[T] += dose[t] ;
					new_evid[T] = 1 ;
					new_cmt[T] = 1 ;
					ad += 1 ;
					t0 = T ;
				}
			}
		}
	}



}

parameters{


	// Population params

	// with individual variations
	real <lower=0> CL_pop ;
	real <lower=0> Q_pop ;
	real <lower=0> VC_pop ;
	real <lower=0> VP_pop ;
	real <lower=0> ke0_pop ;

	// without individual variations
	real <lower=0> ka_pop ; 
	real <lower=0> EC50_pop ;


	// IIV variablity

	cholesky_factor_corr[K] L;
	vector <lower=0> [K] omega ; 
	matrix[K, nSubjects] etaStd;
	

	// Residual var param

	real <lower=0> sigma_c ;
	real <lower=0> sigma_eff ;

}


transformed parameters {

	//params that have individual variations

	vector <lower=0> [K] THETA_pop ; 
	THETA_pop = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ke0_pop}) ;


	//weight-normalized params

	matrix<lower = 0>[nSubjects, K] THETA_norm ;
	THETA_norm = (rep_matrix(THETA_pop, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStd)))' ;


	// Individual params

	array[nSubjects] vector <lower=0> [nTheta] THETA_ind ; 



	// Solving the 2 CPTs model:

	matrix [N, nCmt] X ; 
	vector [N] cHat ;
	vector [nObsPK] cHatObs ;

	vector <lower=0> [N] c_e ;
	vector <lower=0> [N] effHat ;
	vector[nObsPK] effHatObs ;

	for (i in 1:nSubjects) {

	
		THETA_ind[i][1] = THETA_norm[i, 1] .* (weight[i] / 70)^0.75; 
		THETA_ind[i][2] = THETA_norm[i, 2] .* (weight[i] / 70)^0.75;
		THETA_ind[i][3] = THETA_norm[i, 3] .* (weight[i] / 70); 
		THETA_ind[i][4] = THETA_norm[i, 4] .* (weight[i] / 70); 
		THETA_ind[i][5] = ka_pop ;
		THETA_ind[i][6] = THETA_norm[i, 5] ; 
		THETA_ind[i][7] = EC50_pop ;

		
		X[start[i]:end[i], ] = solvetwoCptModel(time[start[i]:end[i]],
										new_amt[start[i]:end[i]],
										new_cmt[start[i]:end[i]],
										new_evid[start[i]:end[i]],
										THETA_ind[i][1], 
										THETA_ind[i][2],
										THETA_ind[i][3], 
										THETA_ind[i][4], 
										THETA_ind[i][5],
										THETA_ind[i][6]) ; 


		cHat[start[i]:end[i]] = X[start[i]:end[i], 2] ./ THETA_ind[i][3] ;

		c_e[start[i]:end[i]] = X[start[i]:end[i], 4] ./ THETA_ind[i][3] ;
		effHat[start[i]:end[i]] = 100 * c_e[start[i]:end[i]] ./ (THETA_ind[i][7] + c_e[start[i]:end[i]]) ;
		
	}

	cHatObs = cHat[iObsPK] ;
	effHatObs = effHat[iObsPK] ;
		
}


model {

	// Population params

	log(THETA_pop[1]) ~ normal(log(10), 1) ;
	log(THETA_pop[2]) ~ normal(log(10), 1) ;
	log(THETA_pop[3]) ~ normal(log(30), 1) ;
	log(THETA_pop[4]) ~ normal(log(100), 1) ;
	log(ka_pop) ~  normal(log(2), 1) ;

	//effect
	//log(THETA_pop[5]) ~ normal(log(2), 1) ;
	//log(EC50_pop) ~ normal(log(150), 1) ;


	// Covariance matrix related params

	L ~ lkj_corr_cholesky(1);
  	to_vector(etaStd) ~ normal(0, 1);
  	omega ~ cauchy(0, 0.5);


	// Residual var.

	sigma_c ~ normal(0, 0.1);
	sigma_eff ~ normal(0, 0.1);

	// Likelihood 
	log(cObs) ~ normal(log(fmax(machine_precision(), cHatObs)), sigma_c) ;
	log(effObs) ~ normal(log(fmax(machine_precision(), effHatObs)), sigma_eff) ;
}


generated quantities {

	matrix[K, nSubjects] etaStdPred ;
	matrix<lower = 0>[nSubjects, K] THETAPred_norm ;
	matrix<lower = 0>[nSubjects, K] THETAPred_ind ;

	matrix[N, nCmt] XPred ;
	vector[N] cHatPred ;
	vector[N] c_ePred ;
	vector[N] effHatPred ;

	vector[N] cPred_pop ; // population predictions
	vector[N] cPred_ind ; // individula predictions
	vector[N] effPred_pop ; 
	vector[N] effPred_ind ; 


	for(i in 1:nSubjects) {
    	for(j in 1:K) {
      		etaStdPred[j, i] = normal_rng(0, 1);
    	}
  	}
	THETAPred_norm = (rep_matrix(THETA_pop, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';


	for(i in 1:nSubjects) {

		// Solving the 2 CPTs model

		THETAPred_ind[i][1] = THETAPred_norm[i, 1] .* (weight[i] / 70)^0.75; 
		THETAPred_ind[i][2] = THETAPred_norm[i, 2] .* (weight[i] / 70)^0.75;
		THETAPred_ind[i][3] = THETAPred_norm[i, 3] .* (weight[i] / 70); 
		THETAPred_ind[i][4] = THETAPred_norm[i, 4] .* (weight[i] / 70); 
		THETAPred_ind[i][5] = THETAPred_norm[i, 5];
		
		XPred[start[i]:end[i], ] = solvetwoCptModel(time[start[i]:end[i]],
													new_amt[start[i]:end[i]],
													new_cmt[start[i]:end[i]],
													new_evid[start[i]:end[i]],
													THETAPred_ind[i][1], 
													THETAPred_ind[i][2],
													THETAPred_ind[i][3], 
													THETAPred_ind[i][4], 
													ka_pop,
													THETAPred_ind[i][5]) ; 

		//concentration
		cHatPred[start[i]:end[i]] = XPred[start[i]:end[i], 2] ./ THETAPred_ind[i][3] ;

		//effect
		c_ePred[start[i]:end[i]] = XPred[start[i]:end[i], 4] ./ THETAPred_ind[i][3] ;
		effHatPred[start[i]:end[i]] = 100 * c_ePred[start[i]:end[i]] ./ (EC50_pop + c_ePred[start[i]:end[i]]) ;
	}

	for (t in 1:N) {

		if(time[t] == 0) {
			cPred_ind[t] = 0 ;
			cPred_pop[t] = 0 ;
			effPred_ind[t] = 0 ;
			effPred_pop[t] = 0 ;
		}
		else {

			cPred_ind[t] = lognormal_rng(log(fmax(machine_precision(), cHat[t])), sigma_c) ;
			cPred_pop[t] = lognormal_rng(log(fmax(machine_precision(), cHatPred[t])), sigma_c) ;

			effPred_ind[t] = lognormal_rng(log(fmax(machine_precision(), effHat[t])), sigma_eff) ;
			effPred_pop[t] = lognormal_rng(log(fmax(machine_precision(), effHatPred[t])), sigma_eff) ;
			


		}
	}



}







