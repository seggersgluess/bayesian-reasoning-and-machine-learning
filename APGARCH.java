package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.SimulatedAnnealing;

public class APGARCH extends GARCH{

	static double [][] asymmetryPars;
	static double asymmetryPower;
	
	public APGARCH(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag);
	}

	
	public static double [][] calc_volatilies_from_GARCH_AP(){
		
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
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4residuals);
				
		for(int m=0; m<lag4residuals; m++){
			double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m+1), n_usedObservations, (lag4observedVariables+1));
			
			for(int t=0; t<n_usedObservations; t++){
				double u_t = 0.0;
				for(int c=0; c<(lag4observedVariables+1); c++){
					u_t += higher_lagged_vars[t][c]*arPars[c][0];
				}
				
				u_t = conv_lagged_vars[t][m+1]-u_t;
				u_t = Math.abs(u_t)-asymmetryPars[m][0]*u_t;
				u_t = Math.pow(u_t, 2.0/asymmetryPower);				
				h_t[t][0] += maPars[m][0]*u_t;		
			}			
		}
		
		for(int t=0; t<n_usedObservations; t++){
			for(int p=0; p<lag4volatility; p++){
				int prevCounter = t-p-1;
				if(prevCounter<0){
					h_t[t][0] += volaPars[1+p][0]*Math.pow(h_t_prev,1.0/asymmetryPower);
				}else{
					h_t[t][0] += volaPars[1+p][0]*h_t[t-p-1][0];
				}
					
			}			
		}
				
		return h_t;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH_A(){
		
    	double [] start_value = get_start_values_4_est_GARCH_AP();
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
    		ArrayList<List<Double>> limits = get_par_limits_4_MLE(start_value);
    		int n_pars = limits.get(0).size();        	
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
        	int n_GARCH_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
        	
    		for(int i=0; i<n_pars; i++){
    			lower_values[i] = limits.get(0).get(i);
    		    upper_values[i] = limits.get(1).get(i);
    		    if(i>=n_GARCH_pars && i<(n_pars-1)){
    		    	lower_values[i] = -1.0;
    		    	upper_values[i] = 1.0;
    		    }
    		    if(i==(n_pars-1)){
    		    	lower_values[i] = 0.0;
    		    	upper_values[i] = 5.0;
    		    }
    		}
    		
        	if(optimizer == "DEoptim"){
        		DifferentialEvolution optim = new DifferentialEvolution(APGARCH::opti_log_likelihood_4_GARCH_AP, 200);
        		optim.set_convergence_criterion(1e-02);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	set_GARCH_AP_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(APGARCH::opti_log_likelihood_4_GARCH_AP, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	set_GARCH_AP_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    	  			
    	boolean truePars = check_GARCH_AP_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("Asymmetric-Power-GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	volatilities = calc_volatilies_from_GARCH_AP();
    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    		
	}
	
	
	public static boolean check_GARCH_AP_par_restrictions(){
		
		boolean truePars = check_GARCH_par_restrictions();
		
		if(truePars == false){
			return truePars;
		}
		
		for(int i=0; i<lag4residuals; i++){
			if(Math.abs(asymmetryPars[i][0])>1.0){
				truePars = false;
				return truePars;
			}
		}
		
		if(asymmetryPower<0.0){
			truePars = false;
		}
		
		return truePars;
		
	}
	
	
	public static double calc_log_likelihood_4_GARCH_AP(){
		
		boolean truePars = check_GARCH_AP_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH_AP();	
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
	
	
	public static double opti_log_likelihood_4_GARCH_AP(double [] pars, double [] further_ars){
		
		set_GARCH_AP_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_AP();
		
		return logLik;
		
	}
	
	
	public static void set_GARCH_AP_pars_from_vec(double [] par_vec){
		
		arPars   = new double [lag4observedVariables+1][1];
		volaPars = new double [lag4volatility+1][1];
		maPars   = new double [lag4residuals][1];
		asymmetryPars = new double [lag4residuals][1];
			
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
			maPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			asymmetryPars[i][0] = par_vec[idx];
			idx++;
		}
		
		asymmetryPower = par_vec[idx];
		
	}
	
	
	public static double [] get_start_values_4_est_GARCH_AP(){
		
		//Note: Start values for asymmetry parameters are set 0.0		
		int n_GARCH_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
		int n_pars = n_GARCH_pars+lag4residuals+1;
		
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
		
		APGARCH obj_arch = new APGARCH(obsData, start_idx, end_idx, obsLag, volaLag, maLag);
		obj_arch.do_MLE_4_GARCH_A();
		
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);
		
		System.out.println("");
		MatrixOperations.print_matrix(asymmetryPars);
		System.out.println(asymmetryPower);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
}
