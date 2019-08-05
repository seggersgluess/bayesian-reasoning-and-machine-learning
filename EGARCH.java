package TimeSeriesAnalysis;

import DataManagement.InputDataManager;
import Mathematics.GammaFunction;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.NumDeriv;
import Optimization.SimulatedAnnealing;

public class EGARCH extends GARCH{

	static double asymmetry_parameter; 
	static double standError_asymmetry_parameter;
	static double tValue_asymmetry_parameter;
	
	public EGARCH(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag, "GED");
		GARCH_model = "E-GARCH";
	}
	
	
	public static double [][] calc_volatilies_from_EGARCH(){
			
		double [][] h_t = new double [n_usedObservations][1];
		double [][] exp_h_t = new double [n_usedObservations][1];
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = volaPars[0][0];
		}
		
		double [][] estValues = calc_est_values_from_GARCH();
		
		double h_t_prev = 0.0;
		
		for(int t=0; t<n_usedObservations; t++){
			h_t_prev += Math.pow(observed_variables[startIdx+t][0]-estValues[t][0],2.0);
		}
    
		h_t_prev /= n_usedObservations;
		h_t_prev = Math.log(h_t_prev);	
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4residuals);
				
		for(int t=0; t<n_usedObservations; t++){	
			
			for(int m=0; m<lag4residuals; m++){
					
				double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m+1), n_usedObservations, (lag4observedVariables+1));					
				double u_t = 0.0;
				
				for(int c=0; c<(lag4observedVariables+1); c++){
					u_t += higher_lagged_vars[t][c]*arPars[c][0];
				}
				
				int prevCounter = t-m-1;
				
				if(prevCounter<0){
					u_t = (conv_lagged_vars[t][m+1]-u_t)/Math.sqrt(Math.exp(h_t_prev));
				}else{
					u_t = (conv_lagged_vars[t][m+1]-u_t)/Math.sqrt(Math.exp(h_t[t-m-1][0]));
				}
				
				u_t = Math.abs(u_t)-error_mean+asymmetry_parameter*u_t;				
				h_t[t][0] += maPars[m][0]*u_t;
	
			}		
			
			for(int p=0; p<lag4volatility; p++){
				int prevCounter = t-p-1;
				if(prevCounter<0){
					h_t[t][0] += volaPars[1+p][0]*h_t_prev;
				}else{
					h_t[t][0] += volaPars[1+p][0]*h_t[t-p-1][0];
				}
					
			}
			
			exp_h_t[t][0] = Math.exp(h_t[t][0]);
			
		}
			
		return exp_h_t;
		
	}
	
	
	public static void set_EGARCH_standard_errors_and_t_values_of_est_pars(){
		
		standErrors_arPars   = new double [lag4observedVariables+1][1];
		standErrors_volaPars = new double [lag4volatility+1][1];
		standErrors_maPars   = new double [lag4residuals][1];
		
		tValues_arPars   = new double [lag4observedVariables+1][1];
		tValues_volaPars = new double [lag4volatility+1][1];
		tValues_maPars   = new double [lag4residuals][1];
		
		int idx = 0;
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			standErrors_arPars[i][0] = estParCovariance[idx][idx];
			tValues_arPars[i][0] = arPars[i][0]/standErrors_arPars[i][0];
			idx++;
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			standErrors_volaPars[i][0] = estParCovariance[idx][idx];
			tValues_volaPars[i][0] = volaPars[i][0]/standErrors_volaPars[i][0];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			standErrors_maPars[i][0] = estParCovariance[idx][idx];
			tValues_maPars[i][0] = maPars[i][0]/standErrors_maPars[i][0];
			idx++;
		}
		
		standError_asymmetry_parameter = estParCovariance[idx][idx];
		tValue_asymmetry_parameter     = asymmetry_parameter/standError_asymmetry_parameter;
		
		standError_tail_parameter      = estParCovariance[idx+1][idx+1];
		tValue_tail_parameter          = tail_parameter/standError_tail_parameter ;
		
	}
	
	
	public static void set_EGARCH_pars_from_vec(double [] par_vec){
		
		arPars   = new double [lag4observedVariables+1][1];
		volaPars = new double [lag4volatility+1][1];
		maPars   = new double [lag4residuals][1];
		
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
		
		asymmetry_parameter = par_vec[idx];
		tail_parameter      = par_vec[idx+1];
			
	}
	
	
	public static boolean check_EGARCH_par_restrictions(){
		
		boolean restrictions_satisfied = true;
		
		if(tail_parameter<=0.0){
			restrictions_satisfied = false;
		}
		
		return restrictions_satisfied;
		
	}
	
	
	public static double calc_log_likelihood_4_EGARCH(){
		
		boolean truePars = check_EGARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		calc_mean_of_generalized_error_dist_4_GARCH();
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_EGARCH();		
		double [][] est_vars = calc_est_values_from_GARCH();
		
		for(int t=0; t<n_usedObservations; t++){
			
			if(h_t[t][0]<=0.0){
				h_t[t][0] = 1e-10;
			}
			
			double residuals = observed_variables[startIdx+t][0]-est_vars[t][0];
			residuals = Math.abs(residuals/(Math.sqrt(Math.exp(h_t[t][0]))*lambda));
			
			logLik += -0.5*(Math.pow(residuals, tail_parameter)+h_t[t][0]);
		
		}
		
		logLik += n_usedObservations*(Math.log(tail_parameter/lambda)-(1.0+1.0/tail_parameter)*Math.log(2.0)-Math.log(GammaFunction.gamma(1.0/tail_parameter)));
			
		if(Double.isInfinite(logLik) == true){
			logLik = -1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = -1e+100;
		}
		
		return -logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_EGARCH(double [] pars, double [] further_ars){
		
		set_EGARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_EGARCH();
		
		return logLik;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_EGARCH(){
		
    	double [] start_value = get_start_values_4_est_GARCH();
    	double [] optimal_value = new double [start_value.length];
    	
    	n_modelPars = start_value.length;
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
    		double boundary_par = 0.8;
    		
        	int n_pars = start_value.length+1;
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
    		
        	for(int i=0; i<(lag4observedVariables+1);i++){
        		if(start_value[i]<0.0){
        			upper_values[i] = start_value[i]*(1.0-boundary_par);
        			lower_values[i] = start_value[i]*(1.0+boundary_par);
        		}
        		if(start_value[i]>0.0){
        			upper_values[i] = start_value[i]*(1.0+boundary_par);
        			lower_values[i] = start_value[i]*(1.0-boundary_par);
        		}
        		if(start_value[i]==0.0){
        			upper_values[i] = boundary_par;
        			lower_values[i] =-boundary_par;
        		}       		
        	}       	
        	
        	int idx = lag4observedVariables+1;
        	
        	for(int i=0; i<(lag4volatility+lag4residuals+1); i++){     		       		
        		if(start_value[idx]<=0.0){
        			upper_values[idx] = boundary_par;
        		}else{
        			upper_values[idx] = start_value[idx]*(1.0+boundary_par);
        		}  
        		idx++;
        	}
        	
        	//Asymmetry parameter
        	upper_values[idx] = 3.0;
        	lower_values[idx] = -3.0;
        	
        	//Tail parameter
        	upper_values[idx+1] = 10.0;
        	lower_values[idx+1] = 0.01;
        	
        	if(optimizer == "DEoptim"){
        		DifferentialEvolution optim = new DifferentialEvolution(EGARCH::opti_log_likelihood_4_EGARCH, 200);
        		optim.set_convergence_criterion(convergence_criterion);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_EGARCH_pars_from_vec(optimal_value);          	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(EGARCH::opti_log_likelihood_4_EGARCH, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_EGARCH_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}   
        	
    	}
    			
    	boolean truePars = check_EGARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("EGARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	if(convergence == false){
    		System.out.println("MLE has not converged for " + GARCH_model);
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    	volatilities = calc_volatilies_from_EGARCH();
    	
    	hessian = NumDeriv.hessian(EGARCH::opti_log_likelihood_4_EGARCH, optimal_value, null);
    	estParCovariance = MatrixOperations.inverse(hessian);
    	set_EGARCH_standard_errors_and_t_values_of_est_pars();
    	
    	calc_information_criteria();
    	
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"UK"};
    	
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
		
		EGARCH obj_arch = new EGARCH(obsData, start_idx, end_idx, obsLag, volaLag, maLag);
			
		obj_arch.do_MLE_4_EGARCH();
		
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);	
		System.out.println(tail_parameter);
		System.out.println(asymmetry_parameter);
		System.out.println("");
		System.out.println("SE´s of parameters:");
		MatrixOperations.print_matrix(standErrors_arPars);
		MatrixOperations.print_matrix(standErrors_volaPars);
		MatrixOperations.print_matrix(standErrors_maPars);	
		System.out.println(standError_tail_parameter);
		System.out.println(standError_asymmetry_parameter);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
}
