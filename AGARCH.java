package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.NumDeriv;
import Optimization.SimulatedAnnealing;

public class AGARCH extends GARCH{

	static double [][] asymmetryPars;
	static double [][] standErrors_asymmetryPars;
	static double [][] tValues_asymmetryPars;
	
	public AGARCH(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag);
		GARCH_model = "A-GARCH";
	}

	
	public static double [][] calc_volatilies_from_GARCH_A(){
		
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
				u_t -= asymmetryPars[m][0];
				u_t = Math.pow(u_t, 2.0);				
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
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH_A(){
		
    	double [] start_value = get_start_values_4_est_GARCH_A();
    	double [] optimal_value = new double [start_value.length];
    	
    	n_modelPars = start_value.length;
    	
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
        		DifferentialEvolution optim = new DifferentialEvolution(AGARCH::opti_log_likelihood_4_GARCH_A, 200);
        		optim.set_convergence_criterion(convergence_criterion);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_GARCH_A_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(AGARCH::opti_log_likelihood_4_GARCH_A, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_GARCH_A_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}   
        	
    	}
    	  			
    	boolean truePars = check_GARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("Asymmetric-GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	if(convergence == false){
    		System.out.println("MLE has not converged for " + GARCH_model);
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	volatilities = calc_volatilies_from_GARCH_A();
    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    		
    	hessian = NumDeriv.hessian(AGARCH::opti_log_likelihood_4_GARCH_A, optimal_value, null);
    	estParCovariance = MatrixOperations.inverse(hessian);
    	set_GARCH_A_standard_errors_and_t_values_of_est_pars();
    	  
    	calc_information_criteria();
    	
	}
	
	
	public static double calc_log_likelihood_4_GARCH_A(){
		
		boolean truePars = check_GARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH_A();	
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
		
		return -logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_GARCH_A(double [] pars, double [] further_ars){
		
		set_GARCH_A_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_A();
		
		return logLik;
		
	}
	
	
	public static void set_GARCH_A_standard_errors_and_t_values_of_est_pars(){
		
		standErrors_arPars   = new double [lag4observedVariables+1][1];
		standErrors_volaPars = new double [lag4volatility+1][1];
		standErrors_maPars   = new double [lag4residuals][1];
		standErrors_asymmetryPars = new double [lag4residuals][1];
		
		tValues_arPars   = new double [lag4observedVariables+1][1];
		tValues_volaPars = new double [lag4volatility+1][1];
		tValues_maPars   = new double [lag4residuals][1];
		tValues_asymmetryPars = new double [lag4residuals][1];
		
		int idx = 0;
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			standErrors_arPars[i][0] = Math.sqrt(estParCovariance[idx][idx]);
			tValues_arPars[i][0] = arPars[i][0]/standErrors_arPars[i][0];
			idx++;
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			standErrors_volaPars[i][0] = Math.sqrt(estParCovariance[idx][idx]);
			tValues_volaPars[i][0] = volaPars[i][0]/standErrors_volaPars[i][0];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			standErrors_maPars[i][0] = Math.sqrt(estParCovariance[idx][idx]);
			tValues_maPars[i][0] = maPars[i][0]/standErrors_maPars[i][0];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			standErrors_asymmetryPars[i][0] = Math.sqrt(estParCovariance[idx][idx]);
			tValues_asymmetryPars[i][0] = asymmetryPars[i][0]/standErrors_asymmetryPars[i][0];
			idx++;
		}
		
	}
	
	
	public static void set_GARCH_A_pars_from_vec(double [] par_vec){
		
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
		
	}
	
	
	public static double [] get_start_values_4_est_GARCH_A(){
		
		//Note: Start values for asymmetry parameters are set 0.0		
		int n_GARCH_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
		int n_pars = n_GARCH_pars+lag4residuals;
		
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
		
		AGARCH obj_arch = new AGARCH(obsData, start_idx, end_idx, obsLag, volaLag, maLag);
		obj_arch.do_MLE_4_GARCH_A();
		
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);	
		MatrixOperations.print_matrix(asymmetryPars);
		System.out.println("");
		System.out.println("SE´s of parameters:");
		MatrixOperations.print_matrix(standErrors_arPars);
		MatrixOperations.print_matrix(standErrors_volaPars);
		MatrixOperations.print_matrix(standErrors_maPars);	
		MatrixOperations.print_matrix(standErrors_asymmetryPars);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
	
}
