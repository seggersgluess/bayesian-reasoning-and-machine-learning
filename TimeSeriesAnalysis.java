package TimeSeriesAnalysis;

import java.awt.Color;

import Graphics.GenGraphics;
import Graphics.HistGraphics;
import Graphics.QQ_Plot;
import Mathematics.MatrixOperations;

public class TimeSeriesAnalysis {

	static String model;
	
    static double [][] observed_variables;
	
	static int n_observations;
	static int n_usedObservations;
	static int n_variables; 
	
	static int startIdx;
	static int endIdx;
	
	static int lag;
		
	static double [][] parameters;
	static double [][] standErrors_parameters;
	static double [][] tValues_parameters;
	
	static double [][] fittedValues;
	static double [][] residuals;
	static double [][] standardized_residuals;
	
	static double sigma;
	
	static double log_likelihood;
	
	static double aic;
	static double bic;
	
	static String ic;
	static double [] ic4modelSelection;
	
	/**
     * Constructor for Moving Average Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public TimeSeriesAnalysis(double [][] obs_variables, int start_idx, int end_idx, int obsLag){
		
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
				
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		
		n_usedObservations = endIdx-startIdx+1;
		
		lag = obsLag;
					
	}
	
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     */
	public TimeSeriesAnalysis(double [][] obs_variables, int start_idx, int end_idx){
		
		if(start_idx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(end_idx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(start_idx >= end_idx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(start_idx == 0) {
			start_idx = 1;
		}
		
		lag = start_idx;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		
		n_usedObservations = endIdx-startIdx+1;
					
	}
		
	
    /**
     * Constructor for Autoregressive Model class
     * @param obs_variables double vector containing data series
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public TimeSeriesAnalysis(double [][] obs_variables, int obsLag){
			
		if(obsLag < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		lag = obsLag;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = lag;
		endIdx                 = obs_variables.length-1;
		
		n_usedObservations = endIdx-startIdx+1;
					
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
		
        double [][] ic_values = get_ic_values_of_model_selection();
        
        if(ic_values != null) {
        	
        	int n_ic_evals = ic_values.length;
    		
            GenGraphics obj_graph = new GenGraphics();
    		       
            obj_graph.setNumberOfPlotColums(1);
            obj_graph.setNumberOfPlotRows(1);
    	 	
            obj_graph.setGraphWidth(1000);
            obj_graph.setGraphHeight(600);
     
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
        }else {
        	System.out.println("No automatic model selection done for " + model + " model.");
        }
        	
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
		
		if(ic4modelSelection == null) {
			return null;
		}
		
		int n = ic4modelSelection.length;
		double [][] ic_values = new double [n][1];
		for(int i=0; i<n; i++){
			ic_values[i][0] = ic4modelSelection[i];
		}
		
		return ic_values;
	}
	
	public static double get_est_sigma() {
		return sigma;
	}
	
	
	public static double [][] get_est_parameters(){
		return parameters;
	}
	
	
	public static double [][] get_parameter_stand_errors(){
		return standErrors_parameters;
	}
	
	
	public static double [][] get_parameter_t_values(){
		return tValues_parameters;
	}
	
	
	public static double [][] get_residuals(){
		return residuals;
	}
	
	
	public static double [][] get_fitted_values(){
		return fittedValues;
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
	
	
}
