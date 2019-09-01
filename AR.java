package TimeSeriesAnalysis;

import DataManagement.InputDataManager;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Regression.LinearRegression;


public class AR extends TimeSeriesAnalysis{

	static double [][] lagged_variables;
	
	static String method;
	
	static double mean;
	
	static double intercept;
	
	static double standError_intercept;
	
	static double tValue_intercept;
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     * @param obsLag integer lag order (max order if automatic selection is used)
     * @param usedMethod string method for selecting estimation procedure for AR
     */
	public AR(double [][] obs_variables, int start_idx, int end_idx, int obsLag, String usedMethod){
		super(obs_variables, start_idx, end_idx, obsLag);
		model = "AR";
		String [] validMethods = get_valid_ar_est_methods();
		
		int [] validIdx = Utilities.Utilities.get_idx(validMethods, usedMethod);
		if(validIdx[0] == -1){
			System.out.println(usedMethod + " is not a valid estimation method for AR. Use Yule-Walker method.");
			method = "Yule-Walker";
		}else{
			method = usedMethod;
		}
		
	    demean_series_4_ar();
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public AR(double [][] obs_variables, int start_idx, int end_idx, int obsLag) {
		super(obs_variables, start_idx, end_idx, obsLag);	
		model = "AR";
		method = "Yule-Walker";
		demean_series_4_ar();	
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     */
	public AR(double [][] obs_variables, int start_idx, int end_idx){
		super(obs_variables, start_idx, end_idx);
		model = "AR";
		method = "Yule-Walker";
	    demean_series_4_ar();					
	}
	
	
	/**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param obsLag integer lag order (max order if automatic selection is used)
     * @param usedMethod string method for selecting estimation procedure for AR
     */
	public AR(double [][] obs_variables, int obsLag, String usedMethod){
		super(obs_variables, obsLag);			
		model = "AR";
		String [] validMethods = get_valid_ar_est_methods();
		
		int [] validIdx = Utilities.Utilities.get_idx(validMethods, usedMethod);
		if(validIdx[0] == -1){
			System.out.println(usedMethod + " is not a valid estimation method for AR. Use Yule-Walker method.");
			method = "Yule-Walker";
		}else{
			method = usedMethod;
		}
		
	    demean_series_4_ar();
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public AR(double [][] obs_variables, int obsLag){			
		super(obs_variables, obsLag);	
		model = "AR";
		method = "Yule-Walker";
	    demean_series_4_ar();
					
	}
	
	public static void do_automatic_ar_model_estimation_and_selection(String usedIC){
		
		ic = usedIC;		
		int orgLag = lag;
		ic4modelSelection = new double [orgLag];
		
		for(int i=0; i<orgLag; i++){
			lag = i+1;
			do_ar_estimation();
			
			if(ic == "aic"){
				ic4modelSelection[i] = aic;
			}
			
			if(ic == "bic"){
				ic4modelSelection[i] = bic;
			}
			
		}
		
		double min_ic  = Utilities.Utilities.getMin(ic4modelSelection);		
		int [] min_idxs = Utilities.Utilities.get_idx(ic4modelSelection, min_ic);
		int optLag = min_idxs[0]+1;
			
		System.out.println("Lag order " + optLag + " identified by minimizing " + ic + ".");
		
		lag = optLag;
		do_ar_estimation();
		
	}
	
	
	public static void do_ar_estimation(){
		
		if(method == "Yule-Walker"){
			do_yule_walker_ar_estimation();
		}
		
		if(method == "OLS"){
			do_ols_ar_estimation();
		}
	
		calc_information_criteria();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_yule_walker_ar_estimation(){
		
		lagged_variables = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag);	
		
		double [][] obsData = MatrixOperations.get_sub_matrix_between_row_idxs(observed_variables, (startIdx-lag), endIdx);
		
		AutoCorrelation ac = new AutoCorrelation(obsData);
		
		ac.calculate_acf(lag);
		
		parameters = new double [lag][1];
		
		double [][] corr_matrix = new double [lag][lag];
		double [][] cov_matrix  = new double [lag][lag];
		double [][] corr_vector = new double [lag][1];
		for(int j=0; j<lag; j++){
			corr_vector[j][0] = ac.autoCorrFunc[j+1][0];
			int idx = 1;
	 		for(int k=(j+1); k<lag; k++){
	 			cov_matrix[k][j] = ac.autoCovFunc[idx][0];
	 			cov_matrix[j][k] = cov_matrix[k][j];
	 			corr_matrix[k][j] = ac.autoCorrFunc[idx][0];	
	 			corr_matrix[j][k] = corr_matrix[k][j];
	 			idx++;
			}
            cov_matrix[j][j] = ac.autoCovFunc[0][0];
	 		corr_matrix[j][j] = 1.0;
		}

		cov_matrix = MatrixOperations.inverse(cov_matrix);
		double [][] pars = MatrixOperations.multiplication(MatrixOperations.inverse(corr_matrix),corr_vector);
		
		for(int i=0; i<lag; i++){
			parameters[i][0] = pars[i][0];
		}
				
		sigma = ac.autoCovFunc[0][0]*(1.0-MatrixOperations.multiplication(MatrixOperations.transpose(corr_vector), pars)[0][0]);
    	sigma *= (endIdx+1.0)/(endIdx-lag+1.0);
		
    	calc_ar_fitted_values_and_residuals();
    	
    	cov_matrix = MatrixOperations.scalar_multiplication(sigma/endIdx, cov_matrix);  		
    	standErrors_parameters = GeneralMath.sqrt(MatrixOperations.get_diagonal_from_matrix(cov_matrix));
		 
    	calculate_t_statistics_4_ar_est_pars();
    	  	
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_ols_ar_estimation(){
	
		lagged_variables = get_lagged_Y_for_lag(observed_variables, startIdx, endIdx, lag);	
		
		double [][] unlagged_variables = MatrixOperations.get_sub_matrix_between_row_idxs(observed_variables, startIdx, endIdx);
		
    	LinearRegression obj_lm = new LinearRegression(unlagged_variables, lagged_variables, true);
    	
    	obj_lm.do_parameter_estimation();
		
    	parameters = obj_lm.get_est_parameters();
    	standErrors_parameters = obj_lm.get_parameter_errors();
    	
    	for(int i=0; i<lag; i++){
    		standErrors_parameters[i][0] *=Math.sqrt((n_usedObservations-(lag+1.0))/n_usedObservations);
    	}
    	
    	intercept = obj_lm.get_est_constant();
    	standError_intercept = obj_lm.get_constant_error();
    	standError_intercept*= Math.sqrt((n_usedObservations-(lag+1.0))/n_usedObservations);
    	 	
    	residuals = obj_lm.get_residuals();
    	fittedValues = obj_lm.get_fitted_values();
    	  	
    	double ss_res        = MatrixOperations.multiplication(MatrixOperations.transpose(residuals), residuals)[0][0];
		sigma                = ss_res/(n_usedObservations);
    	
		calculate_t_statistics_4_ar_est_pars();
		
    	double sd = Math.sqrt(sigma);
    	standardized_residuals = new double [n_usedObservations][1];
    	
    	for(int t=0; t<n_usedObservations; t++){
    		standardized_residuals[t][0] = residuals[t][0]/sd;
    	}
    	
	}
	
	
	public static void calc_ar_fitted_values_and_residuals(){
		
		fittedValues = new double [n_usedObservations][1];
		residuals = new double [n_usedObservations][1];
		standardized_residuals = new double [n_usedObservations][1];
		
		double sd = Math.sqrt(sigma);
		
		for(int t=0; t<n_usedObservations; t++){
			int curIdx = startIdx+t;
			for(int p=0; p<lag; p++){
				fittedValues[t][0] += parameters[p][0]*lagged_variables[t][p];
			}	
			residuals[t][0] = observed_variables[curIdx][0]-fittedValues[t][0];
			standardized_residuals[t][0] = residuals[t][0]/sd;
		}
		
	}
	
	
	public static void calc_information_criteria(){
				
		calc_aic_4_ar();
		calc_bic_4_ar();
		
	}
	
	
	public static void calc_aic_4_ar(){
		aic = Math.log(sigma) + 2.0*lag/n_usedObservations;
	}
	
	
	public static void calc_bic_4_ar(){
		bic = Math.log(sigma) + lag*Math.log(n_usedObservations)/n_usedObservations;
	}
	
	
	// returns the t-values for the estimated parameters
	public static void calculate_t_statistics_4_ar_est_pars(){
		
		tValues_parameters = new double [lag][1];
		
		for(int i=0; i<lag; i++){			
			tValues_parameters[i][0] = parameters[i][0]/standErrors_parameters[i][0];			
		}
			
		if(method == "OLS"){
			tValue_intercept = intercept/standError_intercept;
		}
		
	}
	
	
	public static double get_ar_est_incercept(){
		return intercept;		
	}
	
	
	public static double get_ar_intercept_stand_errors(){
		return standError_intercept;
	}
	
	
	public static double get_ar_intercept_t_value(){
		return tValue_intercept;
	}
	
	
	public static String get_ar_used_method(){
		return method;
	}
	
	
	public static String [] get_valid_ar_est_methods(){
		
		String [] valid_methods = {"Yule-Walker",
				                   "OLS"};	
		return valid_methods;
	}
	
		
	public static void demean_series_4_ar(){
		
		mean = GeneralMath.mean(observed_variables);
		
		for(int t=0; t<n_observations; t++){
			observed_variables[t][0] -= mean;
		}
		
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

		int lag = 5;
		//int start_idx = lag;
		//int end_idx = obsData.length;
		String method = "OLS";
		
		AR obj_ar = new AR(obsData,lag,method);		
		//obj_ar.do_yule_walker_ar_estimation();
		//obj_ar.do_ols_ar_estimation();
			
		obj_ar.do_ar_estimation();
		//obj_ar.do_automatic_ar_model_estimation_and_selection("bic");
		
		obj_ar.plot_time_series_and_fitted_values_with_residuals();
		//obj_ar.plot_time_series_and_fitted_values();
		//obj_ar.plot_residuals();		
		
		//obj_ar.plot_ic_of_model_selection();
		//obj_ar.create_qq_plot();
		//obj_ar.plot_standardized_residuals();
		//obj_ar.plot_acf_of_residuals();
		//obj_ar.plot_pacf_of_residuals();
		//obj_ar.plot_histogram_of_residuals();
		
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(parameters);
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(standErrors_parameters);
		//System.out.println(intercept);
		//System.out.println(standError_intercept);
		//MatrixOperations.print_matrix(residuals);
		
		//calc_information_criteria();
		
	}
	
}
