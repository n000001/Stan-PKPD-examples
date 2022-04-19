functions {

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
		x = matrix_exp((t - t0) * K) * to_vector(init) ;
	   
	   	return x ;
	   	}


	matrix solvetwoCptModel(array[] real time, 
							array[] real amt, array[] int cmt, array[] int evid, 
							real CL, real Q, real V1, real V2, real ka) {

  		int nt = size(time) ;
		vector[3] init = rep_vector(0, 3) ;
		matrix[nt, 3] result;

		for (i in 1:nt){
		
			init = expMatSolve(time[max(1, i - 1)], time[i], init, CL, Q, V1, V2, ka) ;

			// at first iteration t = t0 so init will not change (i.e. (0, 0, 0))
			// then we add the first dose (if it is given at this time)

			if (evid[i] == 1) init[cmt[i]] += amt[i] ; // dose

			result[i, ] = init' ;
		}

		return result;
		}


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
		real conc = params[10];
    
		real ktr = 4.0 / mtt ;
		real prol = x[1] + circ0 ;
		real transit1 = x[2] + circ0 ;
		real transit2 = x[3] + circ0 ;
		real transit3 = x[4] + circ0 ;
		real circ = fmax(machine_precision(), x[5] + circ0) ;
		real EDrug = fmin(1, alpha * conc) ;
		 
		// ODE system 
		array[5] real dxdt ;

    	// PD
		dxdt[1] = ktr * prol * ((1 - EDrug) * ((circ0/circ)^gamma) - 1) ;
		dxdt[2] = ktr * (prol - transit1) ;
		dxdt[3] = ktr * (transit1 - transit2) ;
		dxdt[4] = ktr * (transit2 - transit3) ;
		dxdt[5] = ktr * (transit3 - circ) ;
		
		return dxdt ;
		}
		
		
	matrix solvetwoCptNeut(array[] real time, 
							array[] real amt, array[] int cmt, array[] int evid, 
							array[] real params, vector conc,
							data array[] real rdummy, data array[] int idummy) {
	
		int nCmt = 5;
		int nt = size(time) ;

		array[size(params) + 1] real new_params ;

		array[nCmt] real init ; // must be array bc of integrate_ode_rk45
		init = rep_array(0, nCmt);
		real t0 = time[1] ;
		matrix [nt, nCmt] result ;

		for (i in 1:nt) {
			
			array[1] real c = {conc[i]};
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

			if (evid[i] == 1 && cmt[i] >= 4) {init[cmt[i]-3] += amt[i] ;} // dose

			for (j in 1:nCmt) result[i, j] = init[j] ;
			t0 = time[i] ;
		}
		
	return result ;
	}
			
}


data {

	int nt;
	array[nt] real time;
	int nObsPD ;
	int nObsPK ;
	array[nObsPD] int iObsPD;
	array[nObsPK] int iObsPK;
	
	array[nt] real amt;
	array[nt] int cmt;
	array[nt] int evid;
	array[nObsPK] real cObs;
	array[nObsPD]  real neutObs;
	
	real circ0Prior;
	real circ0PriorCV;
	real CLPrior;
	real CLPriorCV;
	real gammaPrior ;
	real gammaPriorCV;
	real alphaPrior;
	real alphaPriorCV;
	real kaPrior ;
	real kaPriorCV;
	real mttPrior;
	real mttPriorCV;
	real QPrior;
	real QPriorCV;
	real V1Prior;
	real V1PriorCV;
	real V2Prior;
	real V2PriorCV;
	
}


transformed data {

  array[0] real rdummy ; // need for ODE function
  array[0] int idummy ; // need for ODE function

}


parameters {

	// PK parameters
	real < lower = 0 > CL ;
	real < lower = 0 > Q ;
	real < lower = 0 > V1 ;
	real < lower = 0 > V2 ;
	real < lower = 0 > ka ;
	real < lower = 0 > sigma ;

	// PD parameters
	real < lower = 0 > mtt ;	
	real < lower = 0 > circ0 ;
	real < lower = 0 > alpha ;
	real < lower = 0 > gamma ;
	real < lower = 0 > sigmaNeut ; 

}


transformed parameters {
	
  	matrix[nt, 3] y  ;
  	y = solvetwoCptModel(time, amt, cmt, evid, CL, Q, V1, V2, ka)  ;
  	vector[nt] cHat = y[, 2] ./ V1 ;
  	vector[nObsPK] cHatObs = cHat[iObsPK] ;

	array[9] real < lower = 0 > params = {CL, Q, V1, V2, ka, mtt, circ0, gamma, alpha} ;
	matrix[nt, 5] x ;
  	x = solvetwoCptNeut(time, amt, cmt, evid, params, cHat, rdummy, idummy) ;	
  	vector [nt] neutHat = x[, 5] + circ0 ;
	vector [nObsPD] neutHatObs = neutHat[iObsPD] ;

}


model {

	// PK priors
	log(CL) ~ normal(log(CLPrior), CLPriorCV) ;
	log(Q) ~ normal(log(QPrior), QPriorCV) ;
	log(V1) ~ normal(log(V1Prior), V1PriorCV) ;
	log(V2) ~ normal(log(V2Prior), V2PriorCV) ;
	log(ka) ~ normal(log(kaPrior), kaPriorCV) ;
	sigma ~ normal(0, 1) ;
	
	// PK likelihood
	log(cObs) ~ normal(log(cHatObs), sigma) ;

	// PD priors
	log(mtt) ~ normal(log(mttPrior), mttPriorCV) ;
	log(circ0) ~ normal(log(circ0Prior), circ0PriorCV) ;
	log(alpha) ~ normal(log(alphaPrior), alphaPriorCV) ;
	log(gamma) ~ normal(log(gammaPrior), gammaPriorCV) ;
	sigmaNeut ~ normal(0 , 1) ;
	
	// PD likelihood
	log(neutObs) ~ normal(log(fmax(machine_precision(), neutHatObs)), sigmaNeut) ;

}

generated quantities {

	vector[nt] cPred ;
 	vector[nt] neutPred ;

	for (i in 1:nt) {
    	if(time[i] == 0) {
			cPred[i] = 0 ;
    	} 
		else {
			cPred[i] = lognormal_rng(log(fmax(machine_precision(), cHat[i])), sigma) ;
    	}	
    	neutPred[i] = lognormal_rng(log(fmax(machine_precision(), neutHat[i])), sigmaNeut) ;
  }

}



