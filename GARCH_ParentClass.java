package TimeSeriesAnalysis;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import Graphics.GenGraphics;
import Graphics.HistGraphics;
import Mathematics.InformationCriteria;
import Mathematics.MatrixOperations;

public class GARCH_ParentClass {

	static double [][] observed_variables;
	
	static int n_observations;
	static int n_usedObservations;
	static int n_variables; 
	
	static int startIdx;
	static int endIdx;
	
	static String GARCH_model;
	
	static int lag4observedVariables;
	static int lag4volatility;
	static int lag4residuals;
	
	static double [][] arPars;
	static double [][] volaPars;
	static double [][] maPars;
	
	static double [][] standErrors_arPars;
	static double [][] standErrors_volaPars;
	static double [][] standErrors_maPars;
	
	static double [][] tValues_arPars;
	static double [][] tValues_volaPars;
	static double [][] tValues_maPars;
	
	static int n_modelPars;
	static int n_arPars;
	static int n_heteroPars;
	
	static double [][] fittedValues;
	static double [][] residuals;
	static double [][] volatilities;
	
	static String distribution;
	static double student_df;
	static double standError_student_df;
	static double tValue_student_df;
	
	static double logLikelihood;
	
	static ArrayList<List<Double>> lagged_variables;
	
	static String optimizer = "DEoptim";
	
	static boolean convergence;
	static double convergence_criterion = 1e-06;
	
	static double [][] hessian;
	static double [][] estParCovariance;
	
	//information criteria
	static double aic;
	static double bic;
	static double hq;
	
	static int forecastingSteps;
	static double [][] volatilityForecast;
	
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

	 	obj_graph.plotLines(xAxis, obs_data, true);	
	 	obj_graph.plotLines(xAxis, fittedValues, false, Color.RED);
	 	obj_graph.plotLines(xAxis, volatilities, true);
	 	
	 	String title1 = "Observed vs. fitted Variable";
	 	String title2 = "("+ distribution +") " +GARCH_model +"("+lag4volatility+","+ lag4residuals + ") impl. Volatility";
	 	if(GARCH_model == "ARCH"){
	 		title2 = "("+ distribution +") " +GARCH_model +"("+lag4volatility+") impl. Volatility";
	 	}
	 	String [] title  = {title1, title2};
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
	public static void plot_GARCH_volatility_forecast(){
		
		if(volatilityForecast == null){
			throw new RuntimeException("Volatility forecasting for " + GARCH_model + " not yet done.");
		}
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

        int n_values4plot = n_usedObservations+forecastingSteps;
        
	 	double [][] xAxis = new double [n_values4plot][1];
	 	double [][] volas = new double [n_values4plot][1];
	 	
	 	for(int i=0; i<n_values4plot; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		volas[i][0] = volatilities[i][0];
	 	}
	 	
	 	for(int i=0; i<forecastingSteps; i++){
	 		int idx = n_usedObservations+i;
	 		volas[idx][0] = volatilityForecast[i][0];
	 	}

	 	obj_graph.plotLines(xAxis, volas, true);
	 	
	 	String title = "("+ distribution +") " +GARCH_model +"("+lag4volatility+","+ lag4residuals + ") impl. Volatility";
	 	if(GARCH_model == "ARCH"){
	 		title = "("+ distribution +") " +GARCH_model +"("+lag4volatility+") Volatility";
	 	}
	 	
	 	String [] mainTitle  = {title};
	 	String [] yLabel = mainTitle;
	 	String [] subTitle = {"Realized & Forecasted"}; 
	 	obj_graph.setTitle(mainTitle, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setSubTitle1(subTitle, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_histogram_of_residuals(int n_bins){
		
		HistGraphics obj_hist = new HistGraphics();
		 	
	 	//double [] min = {-3.0, -5.0, -15.0, -5.0, -5.0, -10.0};
	 	//double [] max = {5.0, 5.0, 25.0, 1.0, 5.0, 10.0};
	 	
	 	String [] titles = {"Residuals"};
	 	
	 	String subTitle = "("+ distribution +") " +GARCH_model +"("+lag4volatility+","+ lag4residuals + ") impl. Volatility";
	 	if(GARCH_model == "ARCH"){
	 		subTitle = "("+ distribution +") " +GARCH_model +"("+lag4volatility+") Volatility";
	 	}
	 	
	 	String [] subTitles = {subTitle};
	 	String [] yLabels = {""};
	 	String [] xLabels = {"Residuals"};
	 	
	 	obj_hist.plotHistogram(residuals,true,true);

	 	//plotHistogram(y1,false,true);
	 	obj_hist.setNumberOfPlotColums(1);
	 	obj_hist.setNumberOfPlotRows(1);
	 	obj_hist.setGraphWidth(1000);
	 	obj_hist.setGraphHeight(600);
	 	
	 	obj_hist.set_numberOfBins(n_bins);
	 	obj_hist.setTitle(titles, null, null);
	 	obj_hist.setSubTitle1(subTitles, null, null);
	 	obj_hist.setYLabel(yLabels, null, null);
	 	obj_hist.setXLabel(xLabels, null, null);
	 	
	 	obj_hist.setNumberOfDigits4XAxis(2);
	 	obj_hist.setNumberOfDigits4YAxis(2);
	 	obj_hist.setFontOfXAxisUnits("bold",11);
	 	obj_hist.setFontOfYAxisUnits("bold",11);
	 	//set_max_x_value(max);
	 	//set_min_x_value(min);
	 	obj_hist.freqHist(false);	
	 	obj_hist.setPDFLineWidth(2);
	 	//noLinesAroundBars();
	 	obj_hist.plot();
		
	}
		
	
	@SuppressWarnings("static-access")
	public static void calc_information_criteria(){
		
		InformationCriteria ic = new InformationCriteria();
		
		double likelihood = Math.exp(logLikelihood);
		
		aic = ic.aic(likelihood, n_modelPars);
		bic = ic.bic(likelihood, n_modelPars, n_usedObservations);
		hq  = ic.hq(likelihood, n_modelPars, n_usedObservations);
		
	}
	
	
	public static double [][] get_fitted_values(){
		return fittedValues;
	}
	
	
	public static double [][] get_residuals(){
		return residuals;
	}
	
	
	public static double [][] get_volatilities(){
		return volatilities;	
	}
	
	
	public static void set_convergence_criterion(double conv_criterion){
		convergence_criterion = conv_criterion;
	}
	
	
	
}
