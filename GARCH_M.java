package TimeSeriesAnalysis;

import Mathematics.MatrixOperations;

public class GARCH_M extends GARCH{

	public GARCH_M(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag);
	}
	
	//public GARCH_M(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag,String usedDistribution) {
	//	super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag, usedDistribution);
	//}

	
	public static double [][] calc_volatilies_from_GARCH_M(){
		
		double [][] h_t = new double [n_usedObservations][1];
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = volaPars[0][0];
		}
		
		double [][] estValues = calc_est_values_from_GARCH();
		
		double h_t_prev = 0.0;
	
		for(int t=0; t<n_usedObservations; t++){
			h_t_prev += Math.pow(observed_variables[startIdx+t][0]-estValues[t][0],2.0);
		}
    
		h_t_prev /= n_usedObservations;
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4residuals);
				
		for(int m=0; m<lag4residuals; m++){
			double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m+1), n_usedObservations, (lag4observedVariables+1));
			
			for(int t=0; t<n_usedObservations; t++){
				double u_t = 0.0;
				for(int c=0; c<(lag4observedVariables+1); c++){
					u_t += higher_lagged_vars[t][c]*arPars[c][0];
				}
				
				int prevCounter = t-m-1;
				
				if(prevCounter<0){
					u_t = conv_lagged_vars[t][m+1]-u_t-arPars[(lag4observedVariables+2)][0]*Math.sqrt(h_t_prev);
				}else{
					u_t = conv_lagged_vars[t][m+1]-u_t-arPars[(lag4observedVariables+2)][0]*Math.sqrt(h_t[t-m-1][0]);
				}
				
				u_t = Math.pow(u_t,2.0);
				h_t[t][0] += maPars[m][0]*u_t;		
			}			
		}
		
		for(int t=0; t<n_usedObservations; t++){
			for(int p=0; p<lag4volatility; p++){
				int prevCounter = t-p-1;
				if(prevCounter<0){
					h_t[t][0] += volaPars[1+p][0]*h_t_prev;
				}else{
					h_t[t][0] += volaPars[1+p][0]*h_t[t-p-1][0];
				}
					
			}			
		}
				
		return h_t;
		
	}
	
	
	public static void set_GARCH_M_pars_from_vec(double [] par_vec){
		
		arPars   = new double [lag4observedVariables+2][1];
		volaPars = new double [lag4volatility+1][1];
		maPars   = new double [lag4residuals][1];
		
		int idx = 0;
		
		for(int i=0; i<(lag4observedVariables+2); i++){
			arPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			volaPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			maPars[i][0] = par_vec[idx];
			idx++;
		}
				
	}
	
	
	public static double [][] calc_est_values_from_GARCH_M(double [][] h_t){
		
		double [][] prev_obs_values = MatrixOperations.get_matrix_from_vec(lagged_variables.get(0), n_usedObservations, (lag4observedVariables+1));
		
		int n_dates = prev_obs_values.length;
	    int n_vars = prev_obs_values[0].length;	    
		double [][] X = new double [n_dates][n_vars];
		
		for(int i=0; i<n_dates; i++){
			for(int j=0; j<n_vars; j++){
				if(j<(n_vars-1)){
					X[i][j] = prev_obs_values[i][j];
				}else{
					X[i][j] = Math.sqrt(h_t[i][0]);
				}
				
			}
		}
		
		double [][] est_vars = MatrixOperations.multiplication(X, arPars);
		
		return est_vars;
		
	}
	
	
	public static boolean check_GARCH_M_par_restrictions(){
		
		//TODO: Function not implemented yet!
		
	}
	
	
	public static double calc_log_likelihood_4_GARCH_M(){
		
		boolean truePars = check_GARCH_M_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH_M();	
		double [][] est_vars = calc_est_values_from_GARCH_M(h_t);
		
		for(int t=0; t<n_usedObservations; t++){
			
			if(h_t[t][0]<=0.0){
				h_t[t][0] = 1e-10;
			}
			
			logLik += -0.5*(Math.log(h_t[t][0]) + Math.pow((observed_variables[startIdx+t][0]-est_vars[t][0]), 2.0)/h_t[t][0]);
		}
		
		logLik += (-1.0)*n_usedObservations/2.0*Math.log(2.0*Math.PI);
			
		if(Double.isInfinite(logLik) == true){
			logLik = -1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = -1e+100;
		}
		
		System.out.println(-logLik);
		
		return -logLik;
		
	}
	
	
	
	
}
