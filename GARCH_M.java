package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.NewtonMethod;
import Optimization.SimulatedAnnealing;
import Regression.LinearRegression;

public class GARCH_M extends GARCH{

	public static String garch_m_type;
	
	public GARCH_M(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag);
		n_arPars = lag4observedVariables+2;
		garch_m_type = "linear";
	}
	
	//public GARCH_M(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag,String usedDistribution) {
	//	super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag, usedDistribution);
	//}

	
	public static double [][] calc_volatilies_from_GARCH_M(){
		
		double [][] h_t = new double [n_usedObservations][1];
	
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = 1e-100;
		}
		
		double [][] estValues = calc_est_values_from_GARCH_M(h_t);
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = volaPars[0][0];
		}
		
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
					u_t = conv_lagged_vars[t][m+1]-u_t-arPars[(lag4observedVariables+1)][0]*h_t_prev;
				}else{
					u_t = conv_lagged_vars[t][m+1]-u_t-arPars[(lag4observedVariables+1)][0]*h_t[t-m-1][0];
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
	    int n_vars = prev_obs_values[0].length+1;	    
		double [][] X = new double [n_dates][n_vars];
		
		for(int i=0; i<n_dates; i++){
			for(int j=0; j<n_vars; j++){
				if(j<(n_vars-1)){
					X[i][j] = prev_obs_values[i][j];
				}else{
					
					if(garch_m_type == "linear"){
						X[i][j] = h_t[i][0];
					}
					
					if(garch_m_type == "log"){
						X[i][j] = Math.log(h_t[i][0]);
					}
					
					if(garch_m_type == "sqrt"){
						X[i][j] = Math.sqrt(h_t[i][0]);
					}
					
				}			
			}
		}
		
		double [][] est_vars = MatrixOperations.multiplication(X, arPars);
		
		return est_vars;
		
	}
	
	
	public static boolean check_GARCH_M_par_restrictions(){
		
		boolean restrictions_satisfied = check_GARCH_par_restrictions();
		
		return restrictions_satisfied;
		
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
	
	
	public static double opti_log_likelihood_4_GARCH_M(double [] pars, double [] further_ars){
		
		set_GARCH_M_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_M();
		
		return logLik;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH_M(){
		
    	double [] start_value = get_start_values_4_est_GARCH_M();
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
    		ArrayList<List<Double>> limits = get_par_limits_4_MLE(start_value);
    		int n_pars = limits.get(0).size();        	
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
        	
    		for(int i=0; i<n_pars; i++){
    			lower_values[i] = limits.get(0).get(i);
    		    upper_values[i] = limits.get(1).get(i);
    		}
    		
        	if(optimizer == "DEoptim"){
        		DifferentialEvolution optim = new DifferentialEvolution(GARCH_M::opti_log_likelihood_4_GARCH_M, 200);
        		optim.set_convergence_criterion(1e-02);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	set_GARCH_M_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(GARCH_M::opti_log_likelihood_4_GARCH_M, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	set_GARCH_M_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    	
    	if(optimizer == "Newton"){
    		NewtonMethod optim = new NewtonMethod(GARCH_M::opti_log_likelihood_4_GARCH_M, 100000);  	
        	optim.set_convergence_criterion(1e-08);
        	optim.do_Newton_Optimization(start_value);
        	
        	set_GARCH_M_pars_from_vec(optim.get_optimal_candidate());        	
        	logLikelihood = (-1.0)*optim.get_optimal_value();
    	}
    			
    	boolean truePars = check_GARCH_M_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("M-GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	volatilities = calc_volatilies_from_GARCH_M();
    	fittedValues = calc_est_values_from_GARCH_M(volatilities);
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_start_values_4_est_GARCH_M(){
		
		int n_pars = lag4observedVariables+lag4volatility+lag4residuals+3;
		
		int parIdx = 0;
		
		double [] start_values = new double [n_pars];
		
		//First OLS regression w.r.t. observed variables x_t = a + b_1*x_t-1 + ... + b_p*x_t-p + u_t
    	double [][] y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);
    	double [][] X = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4observedVariables);    	
    	
    	LinearRegression obj_lm = new LinearRegression(y, X, false);    	
    	obj_lm.do_parameter_estimation();
		 	
    	for(int i=0; i<(lag4observedVariables+1); i++){
            start_values[parIdx] = (obj_lm.get_est_parameters())[i][0];
            parIdx++;
        }
    	
    	double [][] residuals = obj_lm.get_residuals();
    	int nDates = residuals.length;
    		 
    	y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(residuals, (startIdx+lag4residuals), (nDates-1), 0, 0);
    	X = get_lagged_Y_for_lag(residuals, (startIdx+lag4residuals), (nDates-1), lag4residuals); 
    		
    	obj_lm = new LinearRegression(y, X, false);    	
    	obj_lm.do_parameter_estimation();

    	double [][] estMaPars = obj_lm.get_est_parameters();   	
    	double [][] h_t = MatrixOperations.multiplication(X, estMaPars);    	
    	double [][] H_lagged = get_lagged_Y_for_lag(h_t, (startIdx+lag4volatility), (h_t.length-1), lag4volatility);
    	
    	int n=0;
    	
    	if(H_lagged.length>X.length){
    		n=X.length;
    	}else{
    		n=H_lagged.length;
    	}
    	
    	double [][] X_new = new double [n][(lag4volatility+lag4residuals+1)];
    	
    	for(int t=0; t<n; t++){
    		X_new[t][0] = 1.0;
    		for(int p=0; p<lag4volatility; p++){
    			X_new[t][(1+p)] = H_lagged[t][p+1];
    		}
    		for(int p=0; p<lag4residuals; p++){
    			X_new[t][(lag4volatility+1+p)] = X[t][p+1];
    		}
    	}
    	
    	h_t = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(h_t, 0, (n-1), 0, 0);
    	obj_lm = new LinearRegression(h_t, X_new, false);    	
    	obj_lm.do_parameter_estimation();
    	
    	for(int i=0; i<(lag4volatility+lag4residuals+1); i++){
    		start_values[parIdx] = (obj_lm.get_est_parameters())[i][0];
    		parIdx++;
    	}
		
		return start_values;
		
	}
	
	
	public static String [] get_valid_garch_m_types(){
		
		String [] valid_types = {"linear", "log", "sqrt"};
		
		return valid_types;
		
	}
	
	
	public static void set_garch_m_types(String type){
		
		String [] valid_types = get_valid_garch_m_types();
		int [] idx = Utilities.Utilities.get_idx(valid_types, type);
		
		if(idx[0] == -1){
			throw new RuntimeException(type + " is not a valid type for M-GARCH specification.");
		}
		
		garch_m_type = type;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"Germany"};
    	
		InputDataManager inputData = new InputDataManager();		
		inputData.fileReader(file, true, true, true);
    	
    	int nData = inputData.numberOfRows-1;
    	String [] rownames = new String [nData];
    	for(int i=0; i<nData;i++){
    		rownames[i] = Integer.toString(i+1);
    	}
		
		inputData.selectLoadedData(rownames, colnames);
		double [][] obsData = inputData.selectedDblFileData;

		int obsLag  = 1;
		int volaLag = 1;
		int maLag   = 1;
		int start_idx = obsLag+maLag;
		int end_idx = obsData.length;
		
		GARCH_M obj_arch = new GARCH_M(obsData, start_idx, end_idx, obsLag, volaLag, maLag);
		obj_arch.set_garch_m_types("linear");	
		
		//obj_arch.set_GARCH_optimizer("SANN");
		obj_arch.do_MLE_4_GARCH_M();
		
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
	
}
