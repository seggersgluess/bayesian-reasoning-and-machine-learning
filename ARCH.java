package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

import DataManagement.InputDataManager;
import Mathematics.GammaFunction;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.NumDeriv;
import Optimization.SimulatedAnnealing;
import Regression.LinearRegression;

public class ARCH extends GARCH_ParentClass{
	
	//constructor
	public ARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag){
		
		if(start_idx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(end_idx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(start_idx >= end_idx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(obsLag < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		if(volaLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4volatility = volaLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		int n_lagged_vars    = lag4volatility+1;
		lagged_variables = new ArrayList<List<Double>>(n_lagged_vars);
		
		for(int m=0; m<n_lagged_vars; m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		distribution = "Normal";
		
		GARCH_model = "ARCH";
		
	}
	
	
	public ARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, String usedDistribution){
		
		if(start_idx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(end_idx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(start_idx >= end_idx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(obsLag < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		if(volaLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4volatility = volaLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		int n_lagged_vars    = lag4volatility+1;
		lagged_variables = new ArrayList<List<Double>>(n_lagged_vars);
		
		for(int m=0; m<n_lagged_vars; m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		distribution = usedDistribution;
		
		GARCH_model = "ARCH";
			
	}
	
	
	public static double [][] calc_est_values_from_ARCH(){
		
		double [][]conv_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(0), n_usedObservations, (lag4observedVariables+1));
		double [][] est_vars = MatrixOperations.multiplication(conv_lagged_vars, arPars);
		
		return est_vars;
		
	}
	
	
	public static double [][] calc_volatilies_from_ARCH(){
		
		double [][] h_t = new double [n_usedObservations][1];
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = volaPars[0][0];
		}
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4volatility);
				
		for(int m=0; m<lag4volatility; m++){
			double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m+1), n_usedObservations, (lag4observedVariables+1));
			
			for(int r=0; r<n_usedObservations; r++){
				double u_t = 0.0;
				for(int c=0; c<(lag4observedVariables+1); c++){
					u_t += higher_lagged_vars[r][c]*arPars[c][0];
				}
				u_t = conv_lagged_vars[r][m+1]-u_t;
				u_t = Math.pow(u_t,2.0);
				h_t[r][0] += volaPars[m+1][0]*u_t;
			}
			
		}
		
		return h_t;
		
	}
	
	
	public static double calc_log_likelihood_4_ARCH(){
		
		boolean truePars = check_ARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_ARCH();		
		double [][] est_vars = calc_est_values_from_ARCH();
		
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
	
	
	public static double calc_log_likelihood_4_ARCH_with_Student_Dist(){
		
		boolean truePars = check_ARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_ARCH();		
		double [][] est_vars = calc_est_values_from_ARCH();
		
		for(int t=0; t<n_usedObservations; t++){
			
			if(h_t[t][0]<=0.0){
				h_t[t][0] = 1e-10;
			}
			
			logLik += -0.5*(Math.log(h_t[t][0]) + (student_df+1.0)*Math.log(1.0+Math.pow((observed_variables[startIdx+t][0]-est_vars[t][0]), 2.0)/(h_t[t][0]*(student_df-2.0))));
		}
		
		logLik += n_usedObservations*Math.log((GammaFunction.gamma((student_df+1.0))/2.0)/(Math.sqrt(Math.PI)*GammaFunction.gamma(student_df/2.0))*Math.pow(student_df-2.0, -0.5));
		
		if(Double.isInfinite(logLik) == true){
			logLik = -1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = -1e+100;
		}
		
		return -logLik;
		
	}
	
	
	public static boolean check_ARCH_par_restrictions(){
		
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
		
		if(distribution == "Student"){
			if(student_df <= 2.0){
				restrictions_satisfied = false;
			}
		}
		
		return restrictions_satisfied;
		
	}
	
	
	public static double opti_log_likelihood_4_ARCH(double [] pars, double [] further_ars){
		
		set_ARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_ARCH();
		
		return logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_ARCH_with_Student_Dist(double [] pars, double [] further_ars){
		
		set_ARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_ARCH_with_Student_Dist();
		
		return logLik;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_ARCH(){
		
    	BiFunction<double[], double[], Double> target_function = null;
    	
    	if(distribution == "Normal"){
    		target_function = ARCH::opti_log_likelihood_4_ARCH;
    	}
    	
    	if(distribution == "Student"){
    		target_function = ARCH::opti_log_likelihood_4_ARCH_with_Student_Dist;
    	}
    	
    	double [] start_value = get_start_values_4_est_ARCH();
    	double [] optimal_value = new double [start_value.length];
    	
    	n_modelPars = start_value.length;
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
        	int n_pars = start_value.length;
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
        	
        	double boundary_par = 3.0;
        	
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
        	
        	for(int i=0; i<(lag4volatility+1); i++){     		
        		int idx = lag4observedVariables+1+i;
        		if(start_value[idx]<=0.0){
        			upper_values[idx] = 0.2;
        		}else{
        			upper_values[idx] = start_value[idx]*1.2;
        		}      		
        	}
        	
        	if(distribution == "Student"){
        		lower_values[n_pars-1] = 2.0;
        		upper_values[n_pars-1] = n_usedObservations-1;
        	}
        	
        	if(optimizer == "DEoptim"){
        		DifferentialEvolution optim = new DifferentialEvolution(target_function, 500);
        		optim.set_convergence_criterion(convergence_criterion);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_ARCH_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(target_function, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_ARCH_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}   
        	
    	}
    			
    	boolean truePars = check_ARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	if(convergence == false){
    		System.out.println("MLE has not converged for " + GARCH_model);
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	fittedValues = calc_est_values_from_ARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    	volatilities = calc_volatilies_from_ARCH();
    	
    	hessian = NumDeriv.hessian(target_function, optimal_value, null);
    	estParCovariance = MatrixOperations.inverse(hessian);
    	set_ARCH_standard_errors_and_t_values_of_est_pars();
    	
    	calc_information_criteria();
    	
	}
	
	
	public static void set_ARCH_standard_errors_and_t_values_of_est_pars(){
		
		standErrors_arPars   = new double [lag4observedVariables+1][1];
		standErrors_volaPars = new double [lag4volatility+1][1];
		tValues_arPars       = new double [lag4observedVariables+1][1];
		tValues_volaPars     = new double [lag4volatility+1][1];
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			standErrors_arPars[i][0] = Math.sqrt(estParCovariance[i][i]);
			tValues_arPars[i][0] = arPars[i][0]/standErrors_arPars[i][0];
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			standErrors_volaPars[i][0] = Math.sqrt(estParCovariance[(lag4observedVariables+1+i)][(lag4observedVariables+1+i)]);
			tValues_volaPars[i][0] = volaPars[i][0]/standErrors_volaPars[i][0];
		}
		
		if(distribution == "Student"){
			standError_student_df = Math.sqrt(estParCovariance[(lag4observedVariables+lag4volatility+2)][(lag4observedVariables+lag4volatility+2)]);
		    tValue_student_df = student_df/standError_student_df;
		}
		
	}
	
	
	public static void set_ARCH_pars_from_vec(double [] par_vec){
		
		arPars = new double [lag4observedVariables+1][1];
		volaPars = new double [lag4volatility+1][1];
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			arPars[i][0] = par_vec[i];
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			volaPars[i][0] = par_vec[(lag4observedVariables+1+i)];
		}
		
		if(distribution == "Student"){
			student_df = par_vec[(lag4observedVariables+lag4volatility+2)];
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_start_values_4_est_ARCH(){
		
		int n_pars = lag4observedVariables + lag4volatility + 2;
		
		if(distribution == "Student"){
			n_pars++;
		}
		
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
    	int nResiduals = residuals.length;
    	
        for(int i=0; i<nResiduals; i++){
        	residuals[i][0] = Math.pow(residuals[i][0], 2.0);
        }
    	
    	y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(residuals, (startIdx+lag4volatility), (nResiduals-1), 0, 0);
    	X = get_lagged_Y_for_lag(residuals, (startIdx+lag4volatility), (nResiduals-1), lag4volatility); 
    		
    	obj_lm = new LinearRegression(y, X, false);    	
    	obj_lm.do_parameter_estimation();

    	for(int i=0; i<(lag4volatility+1); i++){
    		start_values[parIdx] = (obj_lm.get_est_parameters())[i][0];
    		parIdx++;
    	}
		    
    	if(distribution == "Student"){
    		start_values[parIdx] = 3.0;
    	}
    	
		return start_values;
		
	}
	

	public static double [][] get_ARCH_fitted_values(){
		return fittedValues;
	}
	
	
	public static double [][] get_ARCH_residuals(){
		return residuals;
	}
	
	
	public static double [][] get_ARCH_volatilities(){
		return volatilities;
	}
	
	
	public static void get_ARCH_forecast_4_Gaussian(int n_steps){
		
		if(n_steps <= 0){
			System.out.println("Invalid number of steps supplied for " + GARCH_model + " volatility forecasting.");
		}
		
		forecastingSteps = n_steps;
		volatilityForecast = new double [n_steps][1];
		
		for(int s=0; s<n_steps; s++){
			volatilityForecast[s][0] += volaPars[0][0];
			int idx=0;
			for(int i=0; i<lag4volatility; i++){
				int histStep = s-i;
				if(histStep <=0){
					int lagIdx = n_usedObservations-i-1+idx;
					double u_t = residuals[lagIdx][0];
					u_t = Math.pow(u_t, 2.0);
					volatilityForecast[s][0] += volaPars[(1+i)][0]*u_t;
				}else{
					volatilityForecast[s][0] += volaPars[(1+i)][0]*volatilityForecast[(s-i-1)][0];
					idx++;
				}				
			}
		}
		
	}
	
	
	public static void get_ARCH_forecast_4_Student_Dist(int n_steps){
		
		if(n_steps <= 0){
			System.out.println("Invalid number of steps supplied for " + GARCH_model + " volatility forecasting.");
		}
		
		forecastingSteps = n_steps;
		volatilityForecast = new double [n_steps][1];
		
		for(int s=0; s<n_steps; s++){
			volatilityForecast[s][0] += volaPars[0][0];
			int idx=0;
			for(int i=0; i<lag4volatility; i++){
				int histStep = s-i;
				if(histStep <=0){
					int lagIdx = n_usedObservations-i-1+idx;
					double u_t = residuals[lagIdx][0];
					u_t = Math.pow(u_t, 2.0);
					volatilityForecast[s][0] += volaPars[(1+i)][0]*u_t;
				}else{
					volatilityForecast[s][0] += (student_df/(student_df-1))*volaPars[(1+i)][0]*volatilityForecast[(s-i-1)][0];
					idx++;
				}			
			}
		}
		
	}
	
	
	public static void get_ARCH_forecast(int n_steps){
		
		if(distribution == "Normal"){
			get_ARCH_forecast_4_Gaussian(n_steps);
		}
	
		if(distribution == "Student"){
			get_ARCH_forecast_4_Student_Dist(n_steps);
		}
		
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

		int obsLag = 2;
		int volaLag = 1;
		int start_idx = obsLag+volaLag;
		int end_idx = obsData.length;
		
		ARCH obj_arch = new ARCH(obsData, start_idx, end_idx, obsLag, volaLag);
		//obj_arch.set_ARCH_optimizer("SANN");		
		obj_arch.do_MLE_4_ARCH();
		
		obj_arch.get_ARCH_forecast(12);
		
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(arPars);
		System.out.println("");
		MatrixOperations.print_matrix(volaPars);
		System.out.println("");
		System.out.println("SE´s of parameters:");
		MatrixOperations.print_matrix(standErrors_arPars);
		System.out.println("");
		MatrixOperations.print_matrix(standErrors_volaPars);
		System.out.println("");
		System.out.println("");
		
		MatrixOperations.print_matrix(volatilityForecast);
		
		//System.out.println(student_df);
		
		plot_GARCH_estimates();
		//plot_GARCH_volatility_forecast();
		//plot_histogram_of_residuals(100);
		
	}
	
	
}
