package TimeSerieFilter;

import java.awt.Color;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class BaxterKing {

	static double [][] observed_variables;
	static double [][] used_observed_variables;
	
	static double [][] trend;
	static double [][] cycle;
	
	static double upper_freq;
	static double lower_freq;
	static double order;
	
	public BaxterKing(double [][] obs_vars, double lower, double upper, double used_order) {
		
		observed_variables = obs_vars;
		
		upper_freq = 2.0*Math.PI/lower;
		lower_freq = 2.0*Math.PI/upper;
		order = used_order;
		
		int n = (int)(observed_variables.length-2.0*order);
		
		cycle = new double [n][1];
		used_observed_variables = new double [n][1];
		
		int int_order = (int) order;
		int n_iterations = 2*int_order+1;
		
		double [] b = new double [n_iterations];
		double sum_b = 0.0;
		int idx = -int_order;
		for(int k=0; k<n_iterations; k++) {
			if(idx==0) {
				b[k] = (upper_freq-lower_freq)/Math.PI;
			}else {
				b[k] = (Math.sin(upper_freq*(double)idx)-Math.sin(lower_freq*(double)idx))/(Math.PI*(double)idx);
			}	
			sum_b += b[k];
			idx++;
		}
				
		double theta = -sum_b/(2.0*order+1.0);
		int idx2 = int_order;	
		
		for(int t=0; t<n; t++) {
			idx = -(int_order);
			for(int k=0; k<n_iterations; k++) {			
				cycle[t][0] += (b[k]+theta)*observed_variables[idx+t+idx2][0];
				idx++;
			}
			used_observed_variables[t][0] = observed_variables[t+idx2][0];	
		}
		
		trend = MatrixOperations.substract(used_observed_variables, cycle);
		
	}
	
	
	public static double [][] get_bk_trend(){
		return trend;
	}
	
	
	public static double [][] get_bk_cycle(){
		return cycle;
	}
	
	
	public static double get_bk_lower_frequency() {
		return lower_freq;
	}
	
	
	public static double get_bk_upper_frequency() {
		return upper_freq;
	}
	
	
	public static double get_bk_order() {
		return order;
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_bk_results(){
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(2);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

        int n = used_observed_variables.length;
        
	 	double [][] xAxis = new double [n][1];
	 	
	 	for(int i=0; i<n; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	obj_graph.plotLines(xAxis, used_observed_variables, true, Color.BLUE);	
	 	obj_graph.plotLines(xAxis, trend, false, Color.RED);
	 	 	
	 	obj_graph.plotLines(xAxis, cycle, true);
	 	
	 	String title1  = "Baxter-King Filter";
	 	String subTitle1 = "Trend Component";
	 	
	 	String title2  = "Baxter-King Filter";
	 	String subTitle2 = "Cyclical Component";
	 	
	 	String [] title = {title1, title2};
	 	String [] subTitle = {subTitle1, subTitle2};
	 	
	 	String [] yLabel = {"Observed/BK-Trend", "BK-Cycle"};
	 		
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
		
		BaxterKing bkFilter = new BaxterKing(obsData, 2.0, 32.0, 7.0);
		
		bkFilter.plot_bk_results();
		
	}
	
}
