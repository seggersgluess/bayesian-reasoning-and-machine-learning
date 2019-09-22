package TimeSerieFilter;

import java.awt.Color;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class TrigonometricRegression {

	static double [][] observed_variables;
	static int T;
	
	static double [][] trend;
	static double [][] cycle;
	
	static double upper_freq;
	static double lower_freq;
	
	static double [] freq_range;
	
	public TrigonometricRegression(double [][] obs_vars, double lower, double upper, boolean conventional_type) {
		
		if(upper < 2.0) {
			throw new RuntimeException("Upper parameter has to be >/= 2.0.");
		}
		
		T = obs_vars.length;
		
		if(lower > T) {
			throw new RuntimeException("Lower parameter has to be </= T.");
		}
		
		observed_variables = obs_vars;
			
		upper_freq = Math.floor(T/lower); 
		lower_freq = Math.ceil(T/upper);
		
		int n_freqs = (int)upper_freq-(int)lower_freq+1;		
		freq_range = new double [n_freqs];
		
		freq_range[0]         = (int) lower_freq;
		freq_range[n_freqs-1] = (int) upper_freq;
		
		int n_add_freq = n_freqs-1;
		int idx = 1;
		for(int i=0; i<(n_add_freq); i++) {			
			freq_range[idx] = (int) lower_freq+idx;
			idx++;
		}
		
		if(conventional_type == true) {
			calc_conv_tr_filter();
		}else {
			calc_mod_tr_filter();
		}
			
	}
	
	
	public TrigonometricRegression(double [][] obs_vars, double lower, double upper) {
		
		if(upper < 2.0) {
			throw new RuntimeException("Upper parameter has to be >/= 2.0.");
		}
		
		T = obs_vars.length;
		
		if(lower > T) {
			throw new RuntimeException("Lower parameter has to be </= T.");
		}
		
		observed_variables = obs_vars;
			
		upper_freq = Math.floor(T/lower); 
		lower_freq = Math.ceil(T/upper);
		
		int n_freqs = (int)upper_freq-(int)lower_freq+1;		
		freq_range = new double [n_freqs];
		
		freq_range[0]         = (int) lower_freq;
		freq_range[n_freqs-1] = (int) upper_freq;
		
		int n_add_freq = n_freqs-1;
		int idx = 1;
		for(int i=0; i<(n_add_freq); i++) {			
			freq_range[idx] = (int) lower_freq+idx;
			idx++;
		}
		
		calc_conv_tr_filter();
			
	}
	
	
	public static void calc_mod_tr_filter() {
		
		cycle = new double [T][1];
		int n_freqs = freq_range.length;
		
		if(upper_freq != (int)T/2.0) {
			for(int t=0; t<T; t++) {
				int time1 = t+1;
				int time2 = time1-T;
				for(int l=time2; l<time1; l++) {
					int obsTimeIdx = time1-1-l;
					double element = 0.0;
					for(int j=0; j<n_freqs; j++) {
						double freq = 2.0*Math.PI*freq_range[j]/T;
						element += Math.cos(freq*l);
					}
					element *=2.0/T*observed_variables[obsTimeIdx][0];
					cycle[t][0] += element;
				}	
			}
		}else {
			for(int t=0; t<T; t++) {
				int time1 = t+1;
				int time2 = time1-T;
				for(int l=time2; l<time1; l++) {
					int obsTimeIdx = time1-1-l;
					double element = 0.0;
					for(int j=0; j<n_freqs; j++) {
						double freq = 2.0*Math.PI*freq_range[j]/T;						
						element += Math.cos(freq*l);						
					}
					element *= 2.0/T;
					element += 1.0/T*Math.cos(Math.PI*(time1-l))*Math.cos(Math.PI*time1);
					element *= observed_variables[obsTimeIdx][0];
					cycle[t][0] += element; 
				}	
			}
		}
		
		trend = MatrixOperations.substract(observed_variables, cycle);
		
	}
	
	
	public static void calc_conv_tr_filter() {
		
		double [][] coeffs = calc_tr_alpha_and_beta_coefficients();
		
		cycle = new double [T][1];
		
		int n_freqs = freq_range.length;
		
		for(int t=0; t<T; t++) {
			for(int i=0; i<n_freqs; i++) {
				double freq = 2.0*Math.PI*freq_range[i]/T;
				double time = t+1;
				cycle[t][0] += coeffs[i][0]*Math.cos(freq*time)+coeffs[i][1]*Math.sin(freq*time);
			}
		}
		
		trend = MatrixOperations.substract(observed_variables, cycle);
		
	}
	
	
	public static double [][] calc_tr_alpha_and_beta_coefficients(){
		
		int n_freqs = freq_range.length;
		double [][] coeffs = new double [n_freqs][2];
		
		for(int i=0; i<n_freqs; i++) {
			double freq = 2.0*Math.PI*freq_range[i]/T;			
			if(freq_range[i] != (int)T/2) {
				for(int t=0; t<T; t++) {
					double time = t+1;
					coeffs[i][0] += Math.cos(freq*time)*observed_variables[t][0];
					coeffs[i][1] += Math.sin(freq*time)*observed_variables[t][0];
				}
				coeffs[i][0] *= 2.0/T;
				coeffs[i][1] *= 2.0/T;
			}else {
				for(int t=0; t<T; t++) {
					double time = t+1;
					coeffs[i][0] += Math.cos(Math.PI*time)*observed_variables[t][0];
					coeffs[i][1] += observed_variables[t][0];//Math.sin(Math.PI*time)*observed_variables[t][0];
				}
				coeffs[i][0] *= 1.0/T;
				coeffs[i][1] *= 1.0/T;
			}
		}
		
		return coeffs;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_tr_results(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(2);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

        int n = observed_variables.length;
        
	 	double [][] xAxis = new double [n][1];
	 	
	 	for(int i=0; i<n; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	obj_graph.plotLines(xAxis, observed_variables, true, Color.BLUE);	
	 	obj_graph.plotLines(xAxis, trend, false, Color.RED);
	 	 	
	 	obj_graph.plotLines(xAxis, cycle, true);
	 	
	 	String title1  = "Trigonometric Regression Filter";
	 	String subTitle1 = "Trend Component";
	 	
	 	String title2  = "Trigonometric Regression Filter";
	 	String subTitle2 = "Cyclical Component";
	 	
	 	String [] title = {title1, title2};
	 	String [] subTitle = {subTitle1, subTitle2};
	 	
	 	String [] yLabel = {"Observed/TR-Trend", "TR-Cycle"};
	 		
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
	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Filter/Unemployment.txt";
    	String [] colnames = {"Unemployment"};
    	
		InputDataManager inputData = new InputDataManager();		
		inputData.fileReader(file, true, true, true);
    	
    	int nData = inputData.numberOfRows-1;
    	String [] rownames = new String [nData];
    	for(int i=0; i<nData;i++){
    		rownames[i] = Integer.toString(i+1);
    	}
		
		inputData.selectLoadedData(rownames, colnames);
		double [][] obsData = inputData.selectedDblFileData;
		
		TrigonometricRegression trFilter = new TrigonometricRegression(obsData,2.0,32.0,false);
		
		trFilter.plot_tr_results();
		
		MatrixOperations.print_matrix(cycle);
		
	}
	
	
}
