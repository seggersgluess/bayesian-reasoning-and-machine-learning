package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.SimulatedAnnealing;

public class GARCH_T extends GARCH{

	static double [][] alpha_plus;
	static double [][] alpha_minus; 
	
	public GARCH_T(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag);
	}

	
	public static double [][] calc_volatilies_from_GARCH_T(){
		
		double [][] h_t = new double [n_usedObservations][1];
		
		double [][] estValues = calc_est_values_from_GARCH();
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = volaPars[0][0];
		}
		
		double h_t_prev = 0.0;
	
		for(int t=0; t<n_usedObservations; t++){
			h_t_prev += Math.pow(observed_variables[startIdx+t][0]-estValues[t][0],2.0);
		}
    
		h_t_prev /= n_usedObservations;
		h_t_prev = Math.sqrt(h_t_prev);
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4residuals);
				
		for(int m=0; m<lag4residuals; m++){
			double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m+1), n_usedObservations, (lag4observedVariables+1));
			
			for(int t=0; t<n_usedObservations; t++){
				double u_t = 0.0;
				for(int c=0; c<(lag4observedVariables+1); c++){
					u_t += higher_lagged_vars[t][c]*arPars[c][0];
				}
				
				u_t = conv_lagged_vars[t][m+1]-u_t;
				
				double u_t_plus  = 0.0;
				double u_t_minus = 0.0;
				if(u_t > 0.0){
					u_t_plus = u_t;
				}else{
					u_t_minus = u_t;
				}
				
				h_t[t][0] += (alpha_plus[m][0]*u_t_plus + alpha_minus[m][0]*u_t_minus)*u_t;		
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
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = Math.pow(h_t[t][0],2.0);
		}
		
		return h_t;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH_T(){
		
    	double [] start_value = get_start_values_4_est_GARCH_T();
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
    		ArrayList<List<Double>> limits = get_par_limits_4_MLE(start_value);
    		int n_pars = limits.get(0).size();        	
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
        	int n_GARCH_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
        	
    		for(int i=0; i<n_pars; i++){
    			lower_values[i] = limits.get(0).get(i);
    		    upper_values[i] = limits.get(1).get(i);
    		    if(i>=n_GARCH_pars){
    		    	lower_values[i] = -3.0;
    		    	upper_values[i] = 3.0;
    		    }
    		}
    		
        	if(optimizer == "DEoptim"){
        		DifferentialEvolution optim = new DifferentialEvolution(GARCH_T::opti_log_likelihood_4_GARCH_T, 200);
        		optim.set_convergence_criterion(1e-02);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	set_GARCH_T_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(GARCH_T::opti_log_likelihood_4_GARCH_T, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	set_GARCH_T_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    	  			
    	boolean truePars = check_GARCH_T_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("T-GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	volatilities = calc_volatilies_from_GARCH_T();
    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    		
	}
	
	
	public static double opti_log_likelihood_4_GARCH_T(double [] pars, double [] further_ars){
		
		set_GARCH_T_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_T();
		
		return logLik;
		
	}
	
	
	public static double calc_log_likelihood_4_GARCH_T(){
		
		boolean truePars = check_GARCH_T_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH_T();	
		double [][] est_vars = calc_est_values_from_GARCH();
		
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
	
	
	public static void set_GARCH_T_pars_from_vec(double [] par_vec){
		
		arPars      = new double [lag4observedVariables+1][1];
		volaPars    = new double [lag4volatility+1][1];
		alpha_plus  = new double [lag4residuals][1];
		alpha_minus = new double [lag4residuals][1];
		
		int idx = 0;
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			arPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			volaPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			alpha_plus[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			alpha_minus[i][0] = par_vec[idx];
			idx++;
		}
		
	}
	
	
	public static boolean check_GARCH_T_par_restrictions(){
		
		boolean restrictions_satisfied = true;
		
		if(volaPars[0][0] <= 0.0){
			restrictions_satisfied = false;
		}
		
		for(int i=0; i<lag4volatility; i++){
			if(volaPars[i+1][0] < 0.0){
				restrictions_satisfied = false;
				break;
			}
		}
			
		return restrictions_satisfied;
		
	}
	
	
	public static double [] get_start_values_4_est_GARCH_T(){
		
		//Note: Start values for asymmetry parameters alpha_plus & alpha_minus are set 0.0		
		int n_GARCH_pars = lag4observedVariables+lag4volatility+2; //Note: GARCH without the maPars (alpha_plus/alpha_minus used instead).
		int n_pars = n_GARCH_pars+2*lag4residuals;
		
		double [] start_values = new double [n_pars];		
		double [] start_values_4_GARCH = get_start_values_4_est_GARCH();
		
		for(int i=0; i<n_GARCH_pars; i++){
			start_values[i] = start_values_4_GARCH[i];
		}
		
		return start_values;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"France"};
    	
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
		int maLag   = 2;
		int start_idx = obsLag+maLag;
		int end_idx = obsData.length;
		
		GARCH_T obj_arch = new GARCH_T(obsData, start_idx, end_idx, obsLag, volaLag, maLag);
		obj_arch.do_MLE_4_GARCH_T();
		
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		
		System.out.println("");
		MatrixOperations.print_matrix(alpha_plus);
		MatrixOperations.print_matrix(alpha_minus);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
}
