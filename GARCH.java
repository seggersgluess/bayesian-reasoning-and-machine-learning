package TimeSeriesAnalysis;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.GammaFunction;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.NewtonMethod;
import Optimization.SimulatedAnnealing;
import Regression.LinearRegression;

public class GARCH {

	static double [][] observed_variables;
	
	static int n_observations;
	static int n_usedObservations;
	static int n_variables; 
	
	static int startIdx;
	static int endIdx;
	
	static int lag4observedVariables;
	static int lag4volatility;
	static int lag4residuals;
	
	static double [][] arPars;
	static double [][] volaPars;
	static double [][] maPars;
	
	static double [][] fittedValues;
	static double [][] residuals;
	static double [][] volatilities;
	
	static String distribution;
	static double student_df;
	
	static double logLikelihood;
	
	static ArrayList<List<Double>> lagged_variables;
	
	static String optimizer = "DEoptim";
	
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
		
		distribution = "Normal";
		
		//Check input consistency!!!
		
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
		
		distribution = usedDistribution;
		
		//Check input consistency!!!
		
	}
	
	
	public static double [][] calc_est_values_from_GARCH(){
		
		double [][]conv_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(0), n_usedObservations, (lag4observedVariables+1));
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
			int prevCounter = lag4volatility-t;
		
			if(prevCounter > 0){
				for(int p=0; p<prevCounter; p++){
					int p_backward = lag4volatility-p-1;
					h_t[t][0] += volaPars[1+p_backward][0]*h_t_prev;
				}
				
				for(int p=0; p<(lag4volatility-prevCounter); p++){
					h_t[t][0] += volaPars[1+p][0]*h_t[t-p-1][0];
				}
			}else{
				for(int p=0; p<lag4volatility; p++){
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
			logLik = 1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = 1e+100;
		}
		
		System.out.println(-logLik);
		
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
			logLik = 1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = 1e+100;
		}
		
		System.out.println(-logLik);
		
		return -logLik;
		
	}
	
	
	public static boolean check_GARCH_par_restrictions(){
		
		boolean restrictions_satisfied = true;
		
		if(volaPars[0][0] <= 0.0){
			restrictions_satisfied = false;
		}
		
		for(int i=0; i<lag4volatility; i++){
			if(volaPars[i][0] < 0.0){
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
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH(){
		
    	BiFunction<double[], double[], Double> target_function = null;
    	
    	if(distribution == "Normal"){
    		target_function = GARCH::opti_log_likelihood_4_GARCH;
    	}
    	
    	if(distribution == "Student"){
    		target_function = GARCH::opti_log_likelihood_4_GARCH_with_Student_Dist;
    	}
    	
    	double [] start_value = get_start_values_4_est_GARCH();
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
        	int n_pars = start_value.length;
        	
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
    		
        	for(int i=0; i<(lag4observedVariables+1);i++){
        		if(start_value[i]<0.0){
        			upper_values[i] = start_value[i]*0.8;
        			lower_values[i] = start_value[i]*1.2;
        		}
        		if(start_value[i]>0.0){
        			upper_values[i] = start_value[i]*1.2;
        			lower_values[i] = start_value[i]*0.8;
        		}
        		if(start_value[i]==0.0){
        			upper_values[i] = 0.2;
        			lower_values[i] =-0.2;
        		}       		
        	}       	
        	
        	for(int i=0; i<(lag4volatility+lag4residuals+1); i++){     		
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
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	set_GARCH_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(target_function, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	set_GARCH_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    	
    	if(optimizer == "Newton"){
    		NewtonMethod optim = new NewtonMethod(target_function, 100000);  	
        	optim.set_convergence_criterion(1e-08);
        	optim.do_Newton_Optimization(start_value);
        	
        	set_GARCH_pars_from_vec(optim.get_optimal_candidate());        	
        	logLikelihood = (-1.0)*optim.get_optimal_value();
    	}
    			
    	boolean truePars = check_GARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    	volatilities = calc_volatilies_from_GARCH();
    	
	}
	
	
	public static void set_GARCH_pars_from_vec(double [] par_vec){
		
		arPars   = new double [lag4observedVariables+1][1];
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
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_start_values_4_est_GARCH(){
		
		int n_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
		
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
		    
    	if(distribution == "Student"){
    		start_values[parIdx] = 3.0;
    	}
    	
		return start_values;
		
	}
	
	
	//sets TxP matrix X_bar = [Y_t-1-m,...,Y_t-p-m] w.r.t. additional lag m beside lag p
	public static List<Double> get_lagged_Y(int m){
		
		List<Double> lagged_Y = new ArrayList<Double>((lag4observedVariables+1));
		
		int T = n_usedObservations;
				
		for(int p=0; p<(lag4observedVariables+1); p++){	
			for(int t=0; t<T; t++){	
				if(p==0){
					lagged_Y.add(1.0);
				}else{
					int curIdx = startIdx+t;
					int lagIdx = curIdx-p-m;
					lagged_Y.add(observed_variables[lagIdx][0]);
				}						
			}
		}
			
		return lagged_Y;
		
	}
	
	
	//sets TxP matrix X_bar = [Y_t-1,...,Y_t-p] for specific lag p
	public static double [][] get_lagged_Y_for_lag(double [][] obs_data, int start_idx, int end_idx, int lag){
		
		int T = end_idx-start_idx+1;
		
		double [][] lagged_Y = new double [T][(lag+1)];
		
		for(int t=0; t<T; t++){
			lagged_Y[t][0] = 1.0;
			int curIdx = start_idx+t;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				lagged_Y[t][(p+1)] = obs_data[lagIdx][0];			
			}
		}
			
		return lagged_Y;
		
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
	public static void plot_GARCH_estimates(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(2);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

	 	//MatrixOperations.print_matrix(volatilities);
	 	
	 	obj_graph.plotLines(xAxis, obs_data, true);	
	 	obj_graph.plotLines(xAxis, fittedValues, false, Color.RED);
	 	obj_graph.plotLines(xAxis, volatilities, true);
	 	
	 	String [] title  = {"Observed vs. fitted Variable", "GARCH("+lag4volatility+","+lag4residuals+") impl. Volatility"};
	 	String [] yLabel = title;
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
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

		int obsLag  = 2;
		int volaLag = 2;
		int maLag   = 2;
		int start_idx = obsLag+maLag;
		int end_idx = obsData.length;
		
		GARCH obj_arch = new GARCH(obsData, start_idx, end_idx, obsLag, volaLag, maLag);
			
		obj_arch.set_GARCH_optimizer("SANN");
		obj_arch.do_MLE_4_GARCH();
		
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);
		
		System.out.println(student_df);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
	
}