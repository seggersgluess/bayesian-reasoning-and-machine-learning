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

public class ARCH {

	static double [][] observed_variables;
	
	static int n_observations;
	static int n_usedObservations;
	static int n_variables; 
	
	static int startIdx;
	static int endIdx;
	
	static int lag4observedVariables;
	static int lag4heteroscedasticity;
	
	static double [][] arPars;
	static double [][] heteroPars;
	
	static double [][] fittedValues;
	static double [][] residuals;
	static double [][] volatilities;
	
	static String distribution;
	static double student_df;
	
	static double logLikelihood;
	
	static ArrayList<List<Double>> lagged_variables;
	
	static String optimizer = "DEoptim";
	
	//constructor
	public ARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int heteroscedasticityLag){
		
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
		
		if(heteroscedasticityLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4heteroscedasticity = heteroscedasticityLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		int n_lagged_vars    = lag4heteroscedasticity+1;
		lagged_variables = new ArrayList<List<Double>>(n_lagged_vars);
		
		for(int m=0; m<n_lagged_vars; m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		distribution = "Normal";
		
		//Check input consistency!!!
		
	}
	
	
	public ARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int heteroscedasticityLag, String usedDistribution){
		
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
		
		if(heteroscedasticityLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4heteroscedasticity = heteroscedasticityLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		int n_lagged_vars    = lag4heteroscedasticity+1;
		lagged_variables = new ArrayList<List<Double>>(n_lagged_vars);
		
		for(int m=0; m<n_lagged_vars; m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		distribution = usedDistribution;
		
		//Check input consistency!!!
		
	}
	
	
	public static double [][] calc_est_values_from_ARCH(){
		
		double [][]conv_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(0), n_usedObservations, (lag4observedVariables+1));
		double [][] est_vars = MatrixOperations.multiplication(conv_lagged_vars, arPars);
		
		return est_vars;
		
	}
	
	
	public static double [][] calc_volatilies_from_ARCH(){
		
		double [][] h_t = new double [n_usedObservations][1];
		
		for(int t=0; t<n_usedObservations; t++){
			h_t[t][0] = heteroPars[0][0];
		}
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag4heteroscedasticity);
				
		for(int m=0; m<lag4heteroscedasticity; m++){
			double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m+1), n_usedObservations, (lag4observedVariables+1));
			
			for(int r=0; r<n_usedObservations; r++){
				double u_t = 0.0;
				for(int c=0; c<(lag4observedVariables+1); c++){
					u_t += higher_lagged_vars[r][c]*arPars[c][0];
				}
				u_t = conv_lagged_vars[r][m+1]-u_t;
				u_t = Math.pow(u_t,2.0);
				h_t[r][0] += heteroPars[m+1][0]*u_t;
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
		
		System.out.println(-logLik);
		
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
		
		System.out.println(-logLik);
		
		return -logLik;
		
	}
	
	
	public static boolean check_ARCH_par_restrictions(){
		
		boolean restrictions_satisfied = true;
		
		if(heteroPars[0][0] <= 0.0){
			restrictions_satisfied = false;
		}
		
		for(int i=0; i<lag4heteroscedasticity; i++){
			if(heteroPars[i+1][0] < 0.0){
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
        	
        	for(int i=0; i<(lag4heteroscedasticity+1); i++){     		
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
            	set_ARCH_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(target_function, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	set_ARCH_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    	
    	if(optimizer == "Newton"){
    		NewtonMethod optim = new NewtonMethod(target_function, 100000);  	
        	optim.set_convergence_criterion(1e-08);
        	optim.do_Newton_Optimization(start_value);
        	
        	set_ARCH_pars_from_vec(optim.get_optimal_candidate());        	
        	logLikelihood = (-1.0)*optim.get_optimal_value();
    	}
    			
    	boolean truePars = check_ARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	fittedValues = calc_est_values_from_ARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    	volatilities = calc_volatilies_from_ARCH();
    	
	}
	
	
	public static void set_ARCH_pars_from_vec(double [] par_vec){
		
		arPars = new double [lag4observedVariables+1][1];
		heteroPars = new double [lag4heteroscedasticity+1][1];
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			arPars[i][0] = par_vec[i];
		}
		
		for(int i=0; i<(lag4heteroscedasticity+1); i++){
			heteroPars[i][0] = par_vec[(lag4observedVariables+1+i)];
		}
		
		if(distribution == "Student"){
			student_df = par_vec[(lag4observedVariables+lag4heteroscedasticity+2)];
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_start_values_4_est_ARCH(){
		
		int n_pars = lag4observedVariables + lag4heteroscedasticity + 2;
		
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
    	
    	y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(residuals, (startIdx+lag4heteroscedasticity), (nResiduals-1), 0, 0);
    	X = get_lagged_Y_for_lag(residuals, (startIdx+lag4heteroscedasticity), (nResiduals-1), lag4heteroscedasticity); 
    		
    	obj_lm = new LinearRegression(y, X, false);    	
    	obj_lm.do_parameter_estimation();

    	for(int i=0; i<(lag4heteroscedasticity+1); i++){
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
	
	
	public static double [][] get_ARCH_fitted_values(){
		return fittedValues;
	}
	
	
	public static double [][] get_ARCH_residuals(){
		return residuals;
	}
	
	
	public static double [][] get_ARCH_volatilities(){
		return volatilities;
	}
	
	
	public static void set_ARCH_optimizer(String opti_algorithm){
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
	 	
	 	String [] title  = {"Observed vs. fitted Variable", "(" + distribution + ") ARCH impl. Volatility"};
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

		int obsLag = 1;
		int heteroscedasticityLag = 1;
		int start_idx = obsLag+heteroscedasticityLag;
		int end_idx = obsData.length;
		
		ARCH obj_arch = new ARCH(obsData, start_idx, end_idx, obsLag, heteroscedasticityLag, "Student");
		obj_arch.set_ARCH_optimizer("SANN");		
		obj_arch.do_MLE_4_ARCH();
		
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(heteroPars);
		
		System.out.println(student_df);
		
		plot_GARCH_estimates();
		
	}
	
	
}
