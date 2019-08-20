package TimeSeriesAnalysis;

import java.awt.Color;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Graphics.HistGraphics;
import Graphics.QQ_Plot;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Regression.LinearRegression;

public class AR {

	static double [][] observed_variables;
	
	static int n_observations;
	static int n_usedObservations;
	static int n_variables; 
	
	static int startIdx;
	static int endIdx;
	
	static int lag;
	
	static double [][] lagged_variables;
	
	static String model = "AR";
	static String method;
	
	static double mean;
	
	static double [][] parameters;
	static double intercept;
	
	static double [][] standErrors_parameters;
	static double standError_intercept;
	
	static double [][] tValues_parameters;
	static double tValue_intercept;
		
	static double [][] fittedValues;
	static double [][] residuals;
	static double [][] standardized_residuals;
	
	static double sigma;
		
	static double aic;
	static double bic;
	
	static String ic;
	static double [] ic4modelSelection;
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     * @param obsLag integer lag order (max order if automatic selection is used)
     * @param usedMethod string method for selecting estimation procedure for AR
     */
	public AR(double [][] obs_variables, int start_idx, int end_idx, int obsLag, String usedMethod){
		
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
		
		if((start_idx-obsLag)<0){
			throw new RuntimeException("Supplied lag order is too long with respect to supplied start idx.");
		}
		
		String [] validMethods = get_valid_ar_est_methods();
		
		int [] validIdx = Utilities.Utilities.get_idx(validMethods, usedMethod);
		if(validIdx[0] == -1){
			System.out.println(usedMethod + " is not a valid estimation method for AR. Use Yule-Walker method.");
			method = "Yule-Walker";
		}else{
			method = usedMethod;
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		
		n_usedObservations = endIdx-startIdx+1;
		
		lag = obsLag;

	    demean_series();
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public AR(double [][] obs_variables, int start_idx, int end_idx, int obsLag){
		
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
		
		if((start_idx-obsLag)<0){
			throw new RuntimeException("Supplied lag order is too long with respect to supplied start idx.");
		}
		
		method = "Yule-Walker";
			
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		
		n_usedObservations = endIdx-startIdx+1;
		
		lag = obsLag;

	    demean_series();
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     */
	public AR(double [][] obs_variables, int start_idx, int end_idx){
		
		if(start_idx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(end_idx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(start_idx >= end_idx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		lag = start_idx;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		
		n_usedObservations = endIdx-startIdx+1;

	    demean_series();
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param obsLag integer lag order (max order if automatic selection is used)
     * @param usedMethod string method for selecting estimation procedure for AR
     */
	public AR(double [][] obs_variables, int obsLag, String usedMethod){
			
		if(obsLag < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		String [] validMethods = get_valid_ar_est_methods();
		
		int [] validIdx = Utilities.Utilities.get_idx(validMethods, usedMethod);
		if(validIdx[0] == -1){
			System.out.println(usedMethod + " is not a valid estimation method for AR. Use Yule-Walker method.");
			method = "Yule-Walker";
		}else{
			method = usedMethod;
		}
		
		lag = obsLag;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = lag;
		endIdx                 = obs_variables.length-1;
		
		n_usedObservations = endIdx-startIdx+1;

	    demean_series();
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public AR(double [][] obs_variables, int obsLag){
			
		if(obsLag < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		method = "Yule-Walker";
		
		lag = obsLag;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = lag;
		endIdx                 = obs_variables.length-1;
		
		n_usedObservations = endIdx-startIdx+1;

	    demean_series();
					
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
	
	
	public static double [][] get_ar_est_parameters(){
		return parameters;
	}
	
	
	public static double get_ar_est_incercept(){
		return intercept;		
	}
	
	
	public static double [][] get_ar_parameter_stand_errors(){
		return standErrors_parameters;
	}
	
	
	public static double get_ar_intercept_stand_errors(){
		return standError_intercept;
	}
	
	
	public static double [][] get_ar_parameter_t_values(){
		return tValues_parameters;
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
	
	
	@SuppressWarnings("static-access")
	public static void plot_time_series(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

	 	obj_graph.plotLines(xAxis, obs_data, true);	
	 	
	 	String [] title  = {"Observed Variable"};
	 	String [] subTitle = {"(Mean Adjusted)"};
	 	String [] yLabel = title;
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setSubTitle1(subTitle, null, "11");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_time_series_and_fitted_values(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

	 	obj_graph.plotPoints(xAxis, obs_data, true, Color.BLUE);	
	 	obj_graph.plotLines(xAxis, fittedValues, false, Color.RED);
	 	 	
	 	String [] title  = {"Observed vs. Fitted Variables"};
	 	String [] subTitle = {"(Mean Adjusted)"};
	 	String [] yLabel = title;
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setSubTitle1(subTitle, null, "11");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_time_series_and_fitted_values_with_residuals(){
		
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

	 	obj_graph.plotPoints(xAxis, obs_data, true, Color.BLUE);	
	 	obj_graph.plotLines(xAxis, fittedValues, false, Color.RED);
	 	 	
	 	obj_graph.plotLines(xAxis, residuals, true);
	 	
	 	String title1  = "Observed vs. Fitted Variables";
	 	String subTitle1 = "(Mean Adjusted)";
	 	
	 	String title2  = "Residuals";
	 	String subTitle2 = "";
	 	
	 	String [] title = {title1, title2};
	 	String [] subTitle = {subTitle1, subTitle2};
	 	
	 	String [] yLabel = title;
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setSubTitle1(subTitle, null, "11");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_residuals(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	obj_graph.plotLines(xAxis, residuals, true);	
	 	
	 	String [] title  = {"Residuals"};
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
	public static void plot_standardized_residuals(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	obj_graph.plotLines(xAxis, standardized_residuals, true);	
	 	
	 	String [] title  = {"Standardized Residuals"};
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
	public static void plot_acf_of_residuals(){
		
		int lag4acf = (int)Math.log10(n_usedObservations)*10;
		if(lag4acf > (n_usedObservations-1)){
			lag4acf = n_usedObservations-1;
		}
		AutoCorrelation obj_acf = new AutoCorrelation(residuals);
		obj_acf.calculate_acf(lag4acf);
		obj_acf.plot_acf("Autocorrelation");
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_pacf_of_residuals(){
		
		int lag4acf = (int)Math.log10(n_usedObservations)*10;
		if(lag4acf > (n_usedObservations-1)){
			lag4acf = n_usedObservations-1;
		}
		AutoCorrelation obj_acf = new AutoCorrelation(residuals);
		obj_acf.calculate_pacf(lag4acf);
		obj_acf.plot_pacf();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_ic_of_model_selection(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

        double [][] ic_values = get_ic_values_of_model_selection();
        int n_ic_evals = ic_values.length;
        
	 	double [][] xAxis = new double [n_ic_evals][1];
	 	
	 	for(int i=0; i<n_ic_evals; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	obj_graph.plotLines(xAxis, ic_values, true, Color.RED);	
	 	
	 	String [] title  = {ic.toUpperCase()};
	 	String [] subTitle = {"(" + model + " Model Selection)"};
	 	String [] yLabel = title;
	 	String [] xLabel = {"Lag Order"};
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setSubTitle1(subTitle, null, "11");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setXLabel(xLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void create_qq_plot(){
		
		QQ_Plot qq_plot = new QQ_Plot();
		qq_plot.set_width(1000);
		qq_plot.set_heigth(600);
		qq_plot.create_qq_plot(standardized_residuals);
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_histogram_of_residuals(){
		
		int nBins = (int) n_usedObservations/4;
		
		HistGraphics obj_graph = new HistGraphics();
		
	 	String [] titles = {"Histogram"};
	 	String [] subTitles = {"Residuals"};
	 	String [] yLabels = {""};
	 	String [] xLabels = {"Residuals"};
	 	
	 	obj_graph.plotHistogram(residuals,true,true);
	 	
	 	obj_graph.setNumberOfPlotColums(1);
	 	obj_graph.setNumberOfPlotRows(1);
	 	
	 	obj_graph.setGraphWidth(1000);
	 	obj_graph.setGraphHeight(600);
	 	
	 	obj_graph.set_numberOfBins(nBins);
	 	
	 	obj_graph.setTitle(titles, null, null);
	 	obj_graph.setSubTitle1(subTitles, null, null);
	 	obj_graph.setYLabel(yLabels, null, null);
	 	obj_graph.setXLabel(xLabels, null, null);
	 	obj_graph.setNumberOfDigits4XAxis(2);
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold",11);
	 	obj_graph.setFontOfYAxisUnits("bold",11);
	 	obj_graph.freqHist(false);	
	 	obj_graph.setPDFLineWidth(2);
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_histogram_of_residuals(int nBins){
		
		HistGraphics obj_graph = new HistGraphics();
		
	 	String [] titles = {"Histogram"};
	 	String [] subTitles = {"Residuals"};
	 	String [] yLabels = {""};
	 	String [] xLabels = {"Residuals"};
	 	
	 	obj_graph.plotHistogram(residuals,true,true);
	 	
	 	obj_graph.setNumberOfPlotColums(1);
	 	obj_graph.setNumberOfPlotRows(1);
	 	
	 	obj_graph.setGraphWidth(1000);
	 	obj_graph.setGraphHeight(600);
	 	
	 	obj_graph.set_numberOfBins(nBins);
	 	
	 	obj_graph.setTitle(titles, null, null);
	 	obj_graph.setSubTitle1(subTitles, null, null);
	 	obj_graph.setYLabel(yLabels, null, null);
	 	obj_graph.setXLabel(xLabels, null, null);
	 	obj_graph.setNumberOfDigits4XAxis(2);
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold",11);
	 	obj_graph.setFontOfYAxisUnits("bold",11);
	 	obj_graph.freqHist(false);	
	 	obj_graph.setPDFLineWidth(2);
	 	obj_graph.plot();
		
	}
	
	
	public static String get_ic_4_model_secection(){
		return ic;
	}
	
	
	public static double [][] get_ic_values_of_model_selection(){
		
		int n = ic4modelSelection.length;
		double [][] ic_values = new double [n][1];
		for(int i=0; i<n; i++){
			ic_values[i][0] = ic4modelSelection[i];
		}
		
		return ic_values;
	}
	
	
	public static double [][] get_residuals(){
		return residuals;
	}
	
	
	public static double [][] get_fitted_values(){
		return fittedValues;
	}
	
	
	public static void demean_series(){
		
		mean = GeneralMath.mean(observed_variables);
		
		for(int t=0; t<n_observations; t++){
			observed_variables[t][0] -= mean;
		}
		
	}
	
	
	//sets TxP matrix X_bar = [Y_t-1,...,Y_t-p] for specific lag p
	public static double [][] get_lagged_Y_for_lag(double [][] obs_data, int start_idx, int end_idx, int lag){
			
		int T = end_idx-start_idx+1;
			
		double [][] lagged_Y = new double [T][lag];
			
		for(int t=0; t<T; t++){
			int curIdx = start_idx+t;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				lagged_Y[t][p] = obs_data[lagIdx][0];			
			}
		}
				
		return lagged_Y;
			
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

		int lag = 11;
		//int start_idx = lag;
		//int end_idx = obsData.length;
		String method = "Yule-Walker";
		
		AR obj_ar = new AR(obsData,lag,method);		
		//obj_ar.do_yule_walker_ar_estimation();
		//obj_ar.do_ols_ar_estimation();
			
		//obj_ar.do_ar_estimation();
		obj_ar.do_automatic_ar_model_estimation_and_selection("bic");
		
		//obj_ar.plot_time_series_and_fitted_values_with_residuals();
		//obj_ar.plot_time_series_and_fitted_values();
		//obj_ar.plot_residuals();		
		
		//obj_ar.plot_ic_of_model_selection();
		obj_ar.create_qq_plot();
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
