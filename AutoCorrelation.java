package TimeSeriesAnalysis;

import java.awt.Color;
import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class AutoCorrelation {

	static double [][] obsData;
	static int n_obsData;
	
	static double [][] autoCovFunc;
	static double [][] autoCorrFunc;
	
	static double [][] partialAutoCorrFunc;
	
	static double lconf;
	static double upconf;
	
	static double [][] pacf;
	static int maxLag;
	
	public AutoCorrelation(double [][] observedData){
		obsData = observedData;
		n_obsData = obsData.length;
		
		lconf = -1.96/Math.sqrt(n_obsData);
		upconf = 1.96/Math.sqrt(n_obsData);
	}
	
	
	public static void calculate_acf1(int lag){
		
		if(lag<=0){
			throw new RuntimeException("Invalid lag supplied for acf calculation.");
		}
		
		maxLag = lag;
		int usedLags = maxLag+1;
		autoCorrFunc = new double [usedLags][1];
		autoCovFunc = new double [usedLags][1];
		
		double mean = GeneralMath.mean(obsData);
		double [][] adj = new double [n_obsData][1];
		for(int t=0; t<n_obsData; t++){
			adj[t][0] = obsData[t][0]-mean;
		}
		
	    double orgSD = GeneralMath.sd(obsData);
		
	    autoCovFunc[0][0] = Math.pow(orgSD, 2.0);
	    autoCorrFunc[0][0] = 1.0;
	    
		for(int i=1; i<usedLags; i++){
			double [][] newData = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(obsData, i, (n_obsData-1), 0, 0);
			mean = GeneralMath.mean(newData);
			double newSD = 0.0;
			for(int t=0; t<(n_obsData-i); t++){
				double newAdj = (newData[t][0]-mean);
				autoCovFunc[i][0] += adj[t][0]*newAdj;
				newSD += Math.pow(newAdj, 2.0);
			}
			newSD/=(n_obsData);
			newSD = Math.sqrt(newSD);
			autoCovFunc[i][0]/= n_obsData;
			autoCorrFunc[i][0] = autoCovFunc[i][0]/(orgSD*newSD);
		}
		
	}
	
	
	public static void calculate_acf(int lag){
		
		if(lag<=0){
			throw new RuntimeException("Invalid lag supplied for acf calculation.");
		}
		
		maxLag = lag;
		int usedLags = maxLag+1;
		autoCorrFunc = new double [usedLags][1];
		autoCovFunc = new double [usedLags][1];
		
		double mean = GeneralMath.mean(obsData);
		double [][] adj = new double [n_obsData][1];
		for(int t=0; t<n_obsData; t++){
			adj[t][0] = obsData[t][0]-mean;
		}
		
	    double variance = GeneralMath.variance(obsData);
	    
	    autoCovFunc[0][0]  = variance;
	    autoCorrFunc[0][0] = 1.0;
	    
		for(int i=1; i<usedLags; i++){
			double [][] newData = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(obsData, i, (n_obsData-1), 0, 0);
			for(int t=0; t<(n_obsData-i); t++){
				double newAdj = (newData[t][0]-mean);
				autoCovFunc[i][0] += adj[t][0]*newAdj;
			}
			autoCovFunc[i][0]/= (n_obsData-1);
			autoCorrFunc[i][0] = autoCovFunc[i][0]/variance;
		}
		
	}
	
	
	public static void calculate_pacf(int lag){
		
	    calculate_acf(lag);
			
		partialAutoCorrFunc = new double [lag][1];
		
		for(int i=0; i<lag; i++){			
			double [][] corr_matrix = new double [i+1][i+1];
			double [][] corr_vector = new double [i+1][1];
			for(int j=0; j<(i+1); j++){
				corr_vector[j][0] = autoCorrFunc[j+1][0];
				int idx = 1;
	 			for(int k=(j+1); k<(i+1); k++){
	 				corr_matrix[k][j] = autoCorrFunc[idx][0];	
	 				corr_matrix[j][k] = corr_matrix[k][j];
	 				idx++;
				}

	 			corr_matrix[j][j] = 1.0;
			}

			double [][] pars = MatrixOperations.multiplication(MatrixOperations.inverse(corr_matrix),corr_vector);
			
			partialAutoCorrFunc[i][0] = pars[i][0];
			
		}
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_acf(String type){
   		
		String [] validTypes = {"Autocovariance", "Autocorrelation"};
		int [] idx = Utilities.get_idx(validTypes,type);
		
		if(idx[0]==-1){
			throw new RuntimeException(type + " is not a valid.");
		}
		
		int usedLags = maxLag+1;
	           	
	 	double [][] x = new double [usedLags][1];	 	
	 	double [][] y1 = new double [usedLags][1];	 	 	
	 	
	 	double [][] xBand = new double [usedLags][2];
	 	double [][] confBand = new double [usedLags][2];
	 	
	 	for(int i =0; i<usedLags ; i++) {
	 		x[i][0] = i;
	 		xBand[i][0] = i;
	 		xBand[i][1] = i;
	 		confBand[i][0] = upconf;
	 		confBand[i][1] = lconf;
	 	}
	 	  
	 	GenGraphics obj_graph = new GenGraphics();
	 	
	 	obj_graph.setNumberOfPlotColums(1);
	 	obj_graph.setNumberOfPlotRows(1);
	 	obj_graph.setGraphWidth(1000);
	 	obj_graph.setGraphHeight(600);
	 	
	 	String [] title = {""};
	 	
	 	if(type == "Autocorrelation"){
	 		obj_graph.plotPoints(x, autoCorrFunc,false, Color.BLUE);	 		
	 		title[0] = "Autocorrelation";
	 	}
	 	
	 	if(type == "Autocovariance"){
	 		obj_graph.plotPoints(x, autoCovFunc,false, Color.BLUE);
	 		title[0] = "Autocovariance";
	 	}
	 	
	 	obj_graph.set_point_width(10);
	 	obj_graph.plotLines(x, y1,false, new Color(200,200,200,200));
	 	obj_graph.drawErrorBars(true);
	 	obj_graph.setColor4ErrorBars(new Color(128,128,128, 200));
	 	
	 	if(type == "Autocorrelation"){
	 		obj_graph.plotShadedArea(xBand, confBand, new Color(100,149,237,200));
	 	}
	 	
	 	String [] yLabel   = title;
	 	String [] xLabel   = {"Lag"};
	 	
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setXLabel(xLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("plain", 10);
	 	obj_graph.setFontOfYAxisUnits("plain", 10);
	 	
	 	obj_graph.plot();
   		
   	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_pacf(){
   			
	 	double [][] x = new double [maxLag][1];	 	
	 	double [][] y1 = new double [maxLag][1];	 	 	
	 	
	 	double [][] xBand = new double [maxLag][2];
	 	double [][] confBand = new double [maxLag][2];
	 	
	 	for(int i =0; i<maxLag ; i++) {
	 		x[i][0] = i+1;
	 		xBand[i][0] = i+1;
	 		xBand[i][1] = i+1;
	 		confBand[i][0] = upconf;
	 		confBand[i][1] = lconf;
	 	}
	 	  
	 	GenGraphics obj_graph = new GenGraphics();
	 	
	 	obj_graph.setNumberOfPlotColums(1);
	 	obj_graph.setNumberOfPlotRows(1);
	 	obj_graph.setGraphWidth(1000);
	 	obj_graph.setGraphHeight(600);
	 	
	 	String [] title = {""};
	 	
	 	obj_graph.plotPoints(x, partialAutoCorrFunc,false, Color.BLUE);	 		
	 	title[0] = "Partial Autocorrelation";
	 	
	 	obj_graph.set_point_width(10);
	 	obj_graph.plotLines(x, y1,false, new Color(200,200,200,200));
	 	obj_graph.drawErrorBars(true);
	 	obj_graph.setColor4ErrorBars(new Color(128,128,128, 200));
	 	

	 	obj_graph.plotShadedArea(xBand, confBand, new Color(100,149,237,200));
	 	
	 	String [] yLabel   = title;
	 	String [] xLabel   = {"Lag"};
	 	
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setXLabel(xLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("plain", 10);
	 	obj_graph.setFontOfYAxisUnits("plain", 10);
	 	
	 	obj_graph.plot();
   		
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
		double [][] observedData = inputData.selectedDblFileData;
		
		AutoCorrelation obj_acf = new AutoCorrelation(observedData);
		
		//obj_acf.calculate_acf(10);
		
		obj_acf.calculate_pacf(25);
		
		MatrixOperations.print_matrix(obj_acf.partialAutoCorrFunc);
		
		obj_acf.plot_pacf();
		//obj_acf.plot_acf("Autocorrelation");
		
	}
	
}
