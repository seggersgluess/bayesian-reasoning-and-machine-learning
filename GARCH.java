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

public class GARCH extends GARCH_ParentClass{

	//Case of Generalized Error Distribution
	static double tail_parameter;
	static double standError_tail_parameter;
	static double tValue_tail_parameter;
	
	static double lambda;
	static double error_mean;
	
	
	//constructor
	public GARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag){
		
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
		
		if(resLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4volatility = volaLag;
		lag4residuals = resLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		int n_lagged_vars    = lag4residuals+1;
		lagged_variables = new ArrayList<List<Double>>(n_lagged_vars);
		
		for(int m=0; m<n_lagged_vars; m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		n_arPars = lag4observedVariables+1;
		n_heteroPars = lag4volatility+lag4residuals+1;
		
		distribution = "Normal";
		
		GARCH_model = "GARCH";
		
	}
	
	
	public GARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag, String usedDistribution){
		
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
		
		if(resLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4volatility = volaLag;
		lag4residuals = resLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		int n_lagged_vars    = lag4residuals+1;
		lagged_variables = new ArrayList<List<Double>>(n_lagged_vars);
		
		for(int m=0; m<n_lagged_vars; m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		n_arPars = lag4observedVariables+1;
		n_heteroPars = lag4volatility+lag4residuals+1;
		
		distribution = usedDistribution;
		
		GARCH_model = "GARCH";
			
	}
	
	
	public static double [][] calc_est_values_from_GARCH(){
		
		double [][] conv_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(0), n_usedObservations, (lag4observedVariables+1));
		double [][] est_vars = MatrixOperations.multiplication(conv_lagged_vars, arPars);
		
		return est_vars;
		
	}
	
	
	public static double [][] calc_volatilies_from_GARCH(){
		
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
				u_t = conv_lagged_vars[t][m+1]-u_t;
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
	
	
	public static double calc_log_likelihood_4_GARCH(){
		
		boolean truePars = check_GARCH_par_restrictions();
 
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH();		
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
	
	
	public static double calc_log_likelihood_4_GARCH_with_Student_Dist(){
		
		boolean truePars = check_GARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH();		
		double [][] est_vars = calc_est_values_from_GARCH();
		
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
	
	
	public static double calc_log_likelihood_4_GARCH_with_Cauchy_Dist(){
		
		boolean truePars = check_GARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH();		
		double [][] est_vars = calc_est_values_from_GARCH();
		
		for(int t=0; t<n_usedObservations; t++){
			
			if(h_t[t][0]<=0.0){
				h_t[t][0] = 1e-10;
			}
			
			logLik += Math.log(h_t[t][0]/(Math.pow((observed_variables[startIdx+t][0]-est_vars[t][0]), 2.0)+Math.pow(h_t[t][0],2.0)));
		}
		
		logLik += n_usedObservations*Math.log(1.0/Math.PI);
		
		if(Double.isInfinite(logLik) == true){
			logLik = -1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = -1e+100;
		}
		
		return -logLik;
		
	}
	
	
	public static double calc_log_likelihood_4_GARCH_with_GED(){
		
		boolean truePars = check_GARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		calc_mean_of_generalized_error_dist_4_GARCH();
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH();		
		double [][] est_vars = calc_est_values_from_GARCH();
		
		for(int t=0; t<n_usedObservations; t++){
			
			if(h_t[t][0]<=0.0){
				h_t[t][0] = 1e-10;
			}
			
			double residuals = observed_variables[startIdx+t][0]-est_vars[t][0];
			residuals = Math.abs(residuals/(Math.sqrt(h_t[t][0])*lambda));
			
			logLik += -0.5*(Math.pow(residuals, tail_parameter)+Math.log(h_t[t][0]));
		
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
	
	
	public static void calc_mean_of_generalized_error_dist_4_GARCH(){
			
		lambda = (Math.pow(2.0, -2.0/tail_parameter)*GammaFunction.gamma(1.0/tail_parameter))/GammaFunction.gamma(3.0/tail_parameter);
		lambda = Math.sqrt(lambda);
		
		error_mean = lambda*Math.pow(2.0, 1.0/tail_parameter)*GammaFunction.gamma(2.0/tail_parameter)/GammaFunction.gamma(1.0/tail_parameter);
		
	}
	
	
	public static boolean check_GARCH_par_restrictions(){
		
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
		
		for(int i=0; i<lag4residuals; i++){
			if(maPars[i][0] < 0.0){
				restrictions_satisfied = false;
				break;
			}
		}
		
		if(distribution == "Student"){
			if(student_df <= 2.0){
				restrictions_satisfied = false;
			}
		}
		
		if(distribution == "GED"){
			if(tail_parameter <= 0.0){
				restrictions_satisfied = false;
			}
		}
		
		return restrictions_satisfied;
		
	}
	
	
	public static double opti_log_likelihood_4_GARCH(double [] pars, double [] further_ars){
		
		set_GARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH();
		
		return logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_GARCH_with_Student_Dist(double [] pars, double [] further_ars){
		
		set_GARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_with_Student_Dist();
		
		return logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_GARCH_with_Cauchy_Dist(double [] pars, double [] further_ars){
		
		set_GARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_with_Cauchy_Dist();
		
		return logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_GARCH_with_GED(double [] pars, double [] further_ars){
		
		set_GARCH_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_with_GED();
		
		return logLik;
		
	}
	
	
	public static ArrayList<List<Double>> get_par_limits_4_MLE(double [] start_value){
		
		double boundary_par = 3.0;
		
    	int n_pars = start_value.length;
    	
    	double [] lower_values = new double [n_pars];
    	double [] upper_values = new double [n_pars];
		
    	for(int i=0; i<(n_arPars);i++){
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
    	
    	int idx = n_arPars;
    	
    	for(int i=0; i<n_heteroPars; i++){     		       		
    		if(start_value[idx]<=0.0){
    			upper_values[idx] = boundary_par;
    		}else{
    			upper_values[idx] = start_value[idx]*(1.0+boundary_par);
    		}  
    		idx++;
    	}
    	
    	if(distribution == "Student"){
    		lower_values[n_pars-1] = 2.0;
    		upper_values[n_pars-1] = n_usedObservations-1;
    	}
    	
    	if(distribution == "GED"){
    		lower_values[n_pars-1] = 0.01;
    		upper_values[n_pars-1] = n_usedObservations-1;
    	}
		
    	ArrayList<List<Double>> limits = new ArrayList<List<Double>>(2);
    	List<Double> lower_limits = new ArrayList<Double>();
    	for(int i=0; i<lower_values.length; i++){
    		lower_limits.add(lower_values[i]);
    	}
 
    	List<Double> upper_limits = new ArrayList<Double>();
    	for(int i=0; i<upper_values.length; i++){
    		upper_limits.add(upper_values[i]);
    	}
    	
    	limits.add(lower_limits);
    	limits.add(upper_limits);
    	
    	return limits;
    	
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH(){
		
    	BiFunction<double[], double[], Double> target_function = null;
    	
    	if(distribution == "Normal"){
    		target_function = GARCH::opti_log_likelihood_4_GARCH;
    	}
    	
    	if(distribution == "Student"){
    		target_function = GARCH::opti_log_likelihood_4_GARCH_with_Student_Dist;
    	}
    	
    	if(distribution == "Cauchy"){
    		target_function = GARCH::opti_log_likelihood_4_GARCH_with_Cauchy_Dist;
    	}
    	
    	if(distribution == "GED"){
    		target_function = GARCH::opti_log_likelihood_4_GARCH_with_GED;
    	}
    	
    	double [] start_value = get_start_values_4_est_GARCH();
    	double [] optimal_value = new double [start_value.length];
    	
    	n_modelPars = start_value.length;
    	
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
        		DifferentialEvolution optim = new DifferentialEvolution(target_function, 500);
        		optim.set_convergence_criterion(convergence_criterion);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_GARCH_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(target_function, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	optimal_value = optim.get_optimal_candidate();
            	set_GARCH_pars_from_vec(optimal_value);        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
            	convergence = optim.get_convergence_info();
        	}   
        	
    	}
    			
    	boolean truePars = check_GARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	if(convergence == false){
    		System.out.println("MLE has not converged for " + GARCH_model);
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    	volatilities = calc_volatilies_from_GARCH();
    	
    	hessian = NumDeriv.hessian(target_function, optimal_value, null);
    	estParCovariance = MatrixOperations.inverse(hessian);
    	set_GARCH_standard_errors_and_t_values_of_est_pars();
    	
    	calc_information_criteria();
    	
	}
	
	
	public static void set_GARCH_standard_errors_and_t_values_of_est_pars(){
		
		standErrors_arPars   = new double [n_arPars][1];
		standErrors_volaPars = new double [lag4volatility+1][1];
		standErrors_maPars   = new double [lag4residuals][1];
		
		tValues_arPars   = new double [n_arPars][1];
		tValues_volaPars = new double [lag4volatility+1][1];
		tValues_maPars   = new double [lag4residuals][1];
		
		int idx=0;
		
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
		
		if(distribution == "Student"){
			standError_student_df = Math.sqrt(estParCovariance[idx][idx]);
			tValue_student_df = student_df/standError_student_df;
			idx++;
		}
		
		if(distribution == "GED"){
			standError_tail_parameter = Math.sqrt(estParCovariance[idx][idx]);
			tValue_tail_parameter = tail_parameter/standError_tail_parameter;
			idx++;
		}
		
	}
	
	
	public static void set_GARCH_pars_from_vec(double [] par_vec){
		
		arPars   = new double [n_arPars][1];
		volaPars = new double [lag4volatility+1][1];
		maPars   = new double [lag4residuals][1];
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			arPars[i][0] = par_vec[i];
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			volaPars[i][0] = par_vec[(lag4observedVariables+1+i)];
		}
		
		for(int i=0; i<lag4residuals; i++){
			maPars[i][0] = par_vec[(lag4observedVariables+lag4volatility+2)+i];
		}
		
		if(distribution == "Student"){
			student_df = par_vec[(lag4observedVariables+lag4volatility+lag4residuals+2)];
		}
		
		if(distribution == "GED"){
			tail_parameter = par_vec[(lag4observedVariables+lag4volatility+lag4residuals+2)];
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_start_values_4_est_GARCH(){
		
		int n_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
		
		if(distribution == "Student"){
			n_pars++;
		}
		
		if(distribution == "GED"){
			n_pars++;
		}
		
		int parIdx = 0;
		
		double [] start_values = new double [n_pars];
		
		//First OLS regression w.r.t. observed variables x_t = a + b_1*x_t-1 + ... + b_p*x_t-p + u_t
    	double [][] y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);
    	double [][] X = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4observedVariables);    	
    	
    	LinearRegression obj_lm = new LinearRegression(y, X, false);
    	
    	obj_lm.do_parameter_estimation();
		
    	for(int i=0; i<n_arPars; i++){
    		start_values[parIdx] = (obj_lm.get_est_parameters())[i][0];
    		parIdx++;
    	}
    	
    	double [][] residuals = obj_lm.get_residuals();
    	int nResiduals = residuals.length;
    	
        for(int i=0; i<nResiduals; i++){
        	residuals[i][0] = Math.pow(residuals[i][0], 2.0);
        }
    	
    	y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(residuals, (startIdx+lag4residuals), (nResiduals-1), 0, 0);
    	X = get_lagged_Y_for_lag(residuals, (startIdx+lag4residuals), (nResiduals-1), lag4residuals); 
    		
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
    	
    	double [][] X_new = new double [n][n_heteroPars];
    	
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
    	
    	for(int i=0; i<n_heteroPars; i++){
    		start_values[parIdx] = (obj_lm.get_est_parameters())[i][0];	
    		parIdx++;
    	}
		      		
    	if(distribution == "Student"){
    		start_values[parIdx] = 3.0;
    	}
    	
    	if(distribution == "GED"){
    		start_values[parIdx] = 2.0;
    	}
    	
		return start_values;
		
	}
	
	
	public static void get_GARCH_forecast_4_Gaussian(int n_steps){
		
		if(n_steps <= 0){
			System.out.println("Invalid number of steps supplied for " + GARCH_model + " volatility forecasting.");
		}
		
		forecastingSteps = n_steps;
		volatilityForecast = new double [n_steps][1];
		
		for(int s=0; s<n_steps; s++){
			volatilityForecast[s][0] += volaPars[0][0];
			int idx=0;
			for(int i=0; i<lag4residuals; i++){
				int histStep = s-i;
				if(histStep <=0){
					int lagIdx = n_usedObservations-i-1+idx;
					double u_t = residuals[lagIdx][0];
					u_t = Math.pow(u_t, 2.0);
					volatilityForecast[s][0] += maPars[i][0]*u_t;
				}else{
					volatilityForecast[s][0] += maPars[i][0]*volatilityForecast[(s-i-1)][0];
					idx++;
				}				
			}
			
			idx=0;
			for(int i=0; i<lag4volatility; i++){				
				int histStep = s-i;
				if(histStep <=0){
					int lagIdx = n_usedObservations-i-1+idx;
					double h_t = volatilities[lagIdx][0];
					volatilityForecast[s][0] += volaPars[(1+i)][0]*h_t;
				}else{
					volatilityForecast[s][0] += volaPars[(1+i)][0]*volatilityForecast[(s-i-1)][0];
					idx++;
				}
			}
			
		}
		
	}
	
	
	public static void get_GARCH_forecast_4_Student_Dist(int n_steps){
		
		if(n_steps <= 0){
			System.out.println("Invalid number of steps supplied for " + GARCH_model + " volatility forecasting.");
		}
		
		forecastingSteps = n_steps;
		volatilityForecast = new double [n_steps][1];
		
		for(int s=0; s<n_steps; s++){
			volatilityForecast[s][0] += volaPars[0][0];
			int idx=0;
			for(int i=0; i<lag4residuals; i++){
				int histStep = s-i;
				if(histStep <=0){
					int lagIdx = n_usedObservations-i-1+idx;
					double u_t = residuals[lagIdx][0];
					u_t = Math.pow(u_t, 2.0);
					volatilityForecast[s][0] += maPars[i][0]*u_t;
				}else{
					volatilityForecast[s][0] += (student_df/(student_df-1))*maPars[i][0]*volatilityForecast[(s-i-1)][0];
					idx++;
				}				
			}
			
			idx=0;
			for(int i=0; i<lag4volatility; i++){				
				int histStep = s-i;
				if(histStep <=0){
					int lagIdx = n_usedObservations-i-1+idx;
					double h_t = volatilities[lagIdx][0];
					volatilityForecast[s][0] += volaPars[(1+i)][0]*h_t;
				}else{
					volatilityForecast[s][0] += volaPars[(1+i)][0]*volatilityForecast[(s-i-1)][0];
					idx++;
				}
			}
			
		}
		
	}

	
	public static void get_GARCH_forecast(int n_steps){
		
		if(distribution == "Normal"){
			get_GARCH_forecast_4_Gaussian(n_steps);
		}
	
		if(distribution == "Student"){
			get_GARCH_forecast_4_Student_Dist(n_steps);
		}
		
	}
	
	
	public static void set_arPars(double [][] ar_pars){
		arPars = ar_pars;
	}
	
	
	public static void set_volaPars(double [][] vola_pars){
		volaPars = vola_pars;
	}
	
	
	public static void set_maPars(double [][] ma_pars){
		maPars = ma_pars;
	}
	
	
	public static double [][] get_GARCH_fitted_values(){
		return fittedValues;
	}
	
	
	public static double [][] get_GARCH_residuals(){
		return residuals;
	}
	
	
	public static double [][] get_GARCH_volatilities(){
		return volatilities;
	}
	
	
	public static void set_GARCH_optimizer(String opti_algorithm){
		optimizer = opti_algorithm;
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
		int maLag   = 1;
		int start_idx = obsLag+maLag;
		int end_idx = obsData.length;
		
		GARCH obj_arch = new GARCH(obsData, start_idx, end_idx, obsLag, volaLag, maLag, "Cauchy");
		obj_arch.set_convergence_criterion(1e-08);
		//obj_arch.set_GARCH_optimizer("SANN");
		obj_arch.do_MLE_4_GARCH();
		//obj_arch.get_GARCH_forecast(100);
		
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);		
		System.out.println("");
		System.out.println("SE�s of parameters:");
		MatrixOperations.print_matrix(standErrors_arPars);
		MatrixOperations.print_matrix(standErrors_volaPars);
		MatrixOperations.print_matrix(standErrors_maPars);	
		
		System.out.println(student_df);
		System.out.println(standError_student_df);
		
		System.out.println(logLikelihood);
		
		//MatrixOperations.print_matrix(estParCovariance);
		
		plot_GARCH_estimates();
		//plot_GARCH_volatility_forecast();
		
	}
	
}
