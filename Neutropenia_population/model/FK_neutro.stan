functions {



	// Matrix exponential based solution of the two CPTs ODE

	vector expMatSolve(real t0, real t, vector init, 
						real CL, real Q, real V1, real V2, real ka) {
 
		matrix[3, 3] K ;
		K = rep_matrix(0, 3, 3) ; 
		K[1, 1] = - ka ;
		K[2, 1] = ka ;
		K[2, 2] = - (CL/V1 + Q/V1) ;
		K[2, 3] = Q/V2 ;
		K[3, 2] = Q/V1 ;
		K[3, 3] = - Q/V2 ;
		
		vector[3] x ;
		x = matrix_exp((t - t0) * K) * init ;
	   
	   	return x ; }



	// Analytical solution of the two CPTs ODE

	matrix solvetwoCptModel(array[] real time, 
							array[] real amt, array[] int cmt, array[] int evid, 
							real CL, real Q, real V1, real V2, real ka) {

  		int nt = size(time) ;
		vector[3] init = rep_vector(0, 3) ;
		matrix[nt, 3] result;

		for (i in 1:nt){
		
			init = expMatSolve(time[max(1, i - 1)], time[i], init, CL, Q, V1, V2, ka) ;

			if (evid[i] == 1) init[cmt[i]] += amt[i] ; // dose

			result[i, ] = init' ;
		}

		return result ; }



	// ODE system for the FK model

	array[] real NeutODE(real t, array[] real x, array[] real params, 
							array[] real rdummy, array[] int idummy) {
		
			real CL = params[1];
			real Q = params[2];
			real V1 = params[3];
			real V2 = params[4];
			real ka = params[5];
			real mtt = params[6];	
			real circ0 = params[7];
			real gamma = params[8];
			real alpha = params[9];

			real conc = params[10]; // analytical solution
		
			real ktr = 4.0 / mtt ;
			real prol = x[1] + circ0 ;
			real transit1 = x[2] + circ0 ;
			real transit2 = x[3] + circ0 ;
			real transit3 = x[4] + circ0 ;
			real circ = fmax(machine_precision(), x[5] + circ0) ;
			real EDrug = fmin(1, alpha * conc) ;
			 
			// FK ODE system 
			array[5] real dxdt ;

			// PD
			dxdt[1] = ktr * prol * ((1 - EDrug) * ((circ0/circ)^gamma) - 1) ;
			dxdt[2] = ktr * (prol - transit1) ;
			dxdt[3] = ktr * (transit1 - transit2) ;
			dxdt[4] = ktr * (transit2 - transit3) ;
			dxdt[5] = ktr * (transit3 - circ) ;
			
			return dxdt ; }



	// Optimization based solution of the FK model

	matrix solvetwoCptNeut(array[] real time, array[] real amt, array[] int cmt, array[] int evid, 
							array[] real params, vector conc,
							data array[] real rdummy, data array[] int idummy) {
		
			int nCmt = 5;
			int nt = size(time) ;
			real t0 = time[1] ;
			array[size(params) + 1] real new_params ; // with cHat 
			array[nCmt] real init ; // must be array bc of integrate_ode_rk45
			init = rep_array(0, nCmt);
			matrix [nt, nCmt] result ;

			for (i in 1:nt) {
				
				array[1] real c = {conc[i]}; // analytical solution for cHat[i]
				new_params = append_array(params, c) ;
				
				if (t0 == time[i]) {;}
				else {init = to_array_1d(integrate_ode_rk45(NeutODE, 
															init, 
															t0, 
															time[i:i], 
															new_params, 
															rdummy,
															idummy)) ;}
															// time[i:i] bc need array type

				for (j in 1:nCmt) result[i, j] = init[j] ;
				t0 = time[i] ;
			}
			
		return result ; }



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

	real mttPrior;
	real mttPriorCV;
	real circ0Prior;
	real circ0PriorCV;
	real gammaPrior ;
	real gammaPriorCV;
	real alphaPrior;
	real alphaPriorCV;


	//Observations

	//PK 
	int<lower = 1> nObsPK ;
	array[nObsPK] int <lower = 1> iObsPK ;
	vector<lower=0>[nObsPK] cObs ;
	// PD
	int<lower = 1> nObsPD ;
	array[nObsPD] int <lower = 1> iObsPD ;
	vector[nObsPD] neutObs ;


	//Individual data

	int<lower=1> nSubjects;
	vector <lower=0> [nSubjects] weight ;
	array[nSubjects] int <lower=0> start;
	array[nSubjects] int <lower=0> end;


	//for easier access
	array[nObsPK] real <lower=0> timePK ; 
	array[nObsPD] real <lower=0> timePD ;
	array[nSubjects] int <lower=0> startPK ;
	array[nSubjects] int <lower=0> endPK ;
	array[nSubjects] int <lower=0> startPD ;
	array[nSubjects] int <lower=0> endPD ;

}


transformed data{


	int nThetaPK = 5 ; 
	int nThetaPD = 4 ;
	int nTheta = 9 ;

	int nCmtPK = 3 ; //two compartments model
	int nCmtPD = 5 ; //FK compartments

	int nIIV = 4 + 3; //number of individual related params
	int K = nIIV ;


	// Chol. factor param 

	int eta = 1 ; //choose value
	

	// Adapt data to multi-dosing schedule 

	array[N] real new_amt = dose ;
	array[N] int new_evid = evid ;
	array[N] int new_cmt = cmt ;
	real <lower=0> eps_dose = 0.05 ; //tolerance for dose administration time
	// there will be a problem in the case where time vector does not have values close to dosing times
	// then it will be necessary to create a new time vector else it will not be accurate

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


	array[0] real rdummy ; // needed for ODE function
	array[0] int idummy ; // needed for ODE function


}



parameters{


	// Population params

	// with individal variations
	real <lower=0> CL_pop ;
	real <lower=0> Q_pop ;
	real <lower=0> VC_pop ;
	real <lower=0> VP_pop ;
	real <lower=0> mtt_pop ;
	real <lower=0> circ0_pop ;
	real <lower=0> alpha_pop ;
	
	// without individual variations
	real <lower=0> ka_pop ;
	real <lower=0> gamma_pop ; 
	


	// IIV variablity

	cholesky_factor_corr[K] L;
	vector <lower=0> [K] omega ; 
	matrix[K, nSubjects] etaStd;
	

	// Residual var

	real <lower=0> sigma ;
	real <lower=0> sigmaNeut ;

}


transformed parameters {

	vector <lower=0> [K] THETA_pop ; 
	THETA_pop = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, mtt_pop, circ0_pop, alpha_pop});

	
	matrix<lower = 0> [nSubjects, K] THETA_norm ; //weight-normalized params
	array[nSubjects, K] real <lower=0> THETA_ind ; //individual params
	
	for (i in 1:nSubjects) {
		THETA_norm = (rep_matrix(THETA_pop, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStd)))' ;
	}


	matrix [N, nCmtPK] X ; 
	vector [N] cHat ;
	vector [nObsPK] cHatObs ;

	array[nSubjects, nTheta] real <lower=0> params_ind ;//all params
	matrix [N, nCmtPD] Y ;
	vector [N] neutHat ;
	vector [nObsPD] neutHatObs ;


	for (i in 1:nSubjects) {


		// Solving the 2 CPTs model:

		THETA_ind[i][1] = THETA_norm[i, 1] .* (weight[i] / 70)^0.75; 
		THETA_ind[i][2] = THETA_norm[i, 2] .* (weight[i] / 70)^0.75;
		THETA_ind[i][3] = THETA_norm[i, 3] .* (weight[i] / 70); 
		THETA_ind[i][4] = THETA_norm[i, 4] .* (weight[i] / 70); 
		
		X[start[i]:end[i], ] = solvetwoCptModel(time[start[i]:end[i]],
										new_amt[start[i]:end[i]],
										new_cmt[start[i]:end[i]],
										new_evid[start[i]:end[i]],
										THETA_ind[i][1], 
										THETA_ind[i][2],
										THETA_ind[i][3], 
										THETA_ind[i][4], 
										ka_pop) ; 

		cHat[start[i]:end[i]] = X[start[i]:end[i], 2] ./ THETA_ind[i][3] ;


		// Solving the Fk model

		THETA_ind[i][5] = THETA_norm[i, 5] ; //mtt
		THETA_ind[i][6] = THETA_norm[i, 6] ;//circ0
		THETA_ind[i][7] = THETA_norm[i, 7] ;//alpha

		params_ind[i] = {THETA_ind[i][1], 
						 THETA_ind[i][2], 
						 THETA_ind[i][3], 
						 THETA_ind[i][4], 
						 ka_pop, 
						 THETA_ind[i][5],
						 THETA_ind[i][6],
						 gamma_pop,
						 THETA_ind[i][7]} ;


		
		Y[start[i]:end[i], ] = solvetwoCptNeut(time[start[i]:end[i]], 
											 new_amt[start[i]:end[i]], 
											 new_cmt[start[i]:end[i]], 
											 new_evid[start[i]:end[i]], 
											 params_ind[i], 
											 cHat[start[i]:end[i]], 
											 rdummy, idummy) ;
		
	  	neutHat[start[i]:end[i]] = Y[start[i]:end[i], 5] + THETA_ind[i][6] ;

		
	}

	cHatObs = cHat[iObsPK] ; 
	neutHatObs = neutHat[iObsPD] ;
	
	
}
	
	

model {

	// Population params priors

	log(THETA_pop[1]) ~ normal(log(10), 1) ;
	log(THETA_pop[2]) ~ normal(log(10), 1) ;
	log(THETA_pop[3]) ~ normal(log(30), 1) ;
	log(THETA_pop[4]) ~ normal(log(100),1) ;
	log(ka_pop) ~  normal(log(2), 0.5) ;

	log(THETA_pop[5]) ~ normal(log(mttPrior), mttPriorCV) ;
	log(THETA_pop[6]) ~ normal(log(circ0Prior), circ0PriorCV) ;
	log(THETA_pop[7]) ~ normal(log(alphaPrior), alphaPriorCV) ;
	log(gamma_pop) ~ normal(log(gammaPrior), gammaPriorCV) ;


	// Covariance matrix related params

	L ~ lkj_corr_cholesky(1);
  	to_vector(etaStd) ~ normal(0, 1);
  	omega ~ cauchy(0, 0.5);


	// Residual var.

	sigma ~ normal(0, 0.1);
	sigmaNeut ~ normal(0, 0.1) ;


	// Likelihood 

	log(cObs) ~ normal(fmax(machine_precision(), log(cHatObs)), sigma) ;
	log(neutObs) ~ normal(fmax(machine_precision(), log(neutHatObs)), sigmaNeut);


}


generated quantities {

	vector[N] cPred ;
	vector[N] neutPred ;

	for (t in 1:N) {

		if(time[t] == 0) cPred[t] = 0 ;
		else cPred[t] = lognormal_rng(log(fmax(machine_precision(), cHat[t])), sigma) ;

		neutPred[t] = lognormal_rng(log(fmax(machine_precision(), neutHat[t])), sigmaNeut) ;
	}













}









