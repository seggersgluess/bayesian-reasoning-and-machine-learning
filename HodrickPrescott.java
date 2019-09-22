package TimeSerieFilter;

import java.awt.Color;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class HodrickPrescott {

	static double [][] observed_variables;
	
	static double [][] trend;
	static double [][] cycle;
	
	static double penalty;
	
	
	public HodrickPrescott(double [][] obs_vars, double lambda_penalty_par) {
		
		if(obs_vars.length<5) {
			throw new RuntimeException("Not enough dates supplied for applying HP filter.");
		}
		
		observed_variables = obs_vars;
		penalty = lambda_penalty_par;
		
		double [][] hp_matrix = calc_hp_matrix();
		
		trend = MatrixOperations.multiplication(MatrixOperations.inverse(hp_matrix),observed_variables);
		cycle = MatrixOperations.substract(observed_variables, trend);
		
	}
	
	
	public static double [][] calc_hp_matrix(){
		
		int n = observed_variables.length;
		double [][] hp_matrix = new double [n][n];
	
		double first_row_coeff_1 = 1.0 + penalty;
		double first_row_coeff_2 = -2.0*penalty;
		double first_row_coeff_3 = penalty;
		
		double sec_row_coeff_1 = -2.0*penalty;
		double sec_row_coeff_2 = 1.0+5.0*penalty;
		double sec_row_coeff_3 = -4.0*penalty;
		double sec_row_coeff_4 = penalty;
		
		double middle_row_coeff_1 = penalty;
		double middle_row_coeff_2 = -4.0*penalty;
		double middle_row_coeff_3 = 1.0+6.0*penalty;
		double middle_row_coeff_4 = -4.0*penalty;
		double middle_row_coeff_5 = penalty;
		
		hp_matrix[0][0] = first_row_coeff_1;
		hp_matrix[0][1] = first_row_coeff_2;
		hp_matrix[0][2] = first_row_coeff_3;
		
		hp_matrix[(n-1)][(n-1)] = first_row_coeff_1;
		hp_matrix[(n-1)][(n-2)] = first_row_coeff_2;
		hp_matrix[(n-1)][(n-3)] = first_row_coeff_3;
		
		hp_matrix[1][0] = sec_row_coeff_1;
		hp_matrix[1][1] = sec_row_coeff_2;
		hp_matrix[1][2] = sec_row_coeff_3;
		hp_matrix[1][3] = sec_row_coeff_4;
		
		hp_matrix[(n-2)][(n-1)] = sec_row_coeff_1;
		hp_matrix[(n-2)][(n-2)] = sec_row_coeff_2;
		hp_matrix[(n-2)][(n-3)] = sec_row_coeff_3;
		hp_matrix[(n-2)][(n-4)] = sec_row_coeff_4;
		
		int nRows = n-4;
		int rowIdx = 2;
		int colIdx = 0;
	
		for(int i=0; i<nRows; i++) {			
			hp_matrix[rowIdx][colIdx]   = middle_row_coeff_1;
			hp_matrix[rowIdx][colIdx+1] = middle_row_coeff_2;
			hp_matrix[rowIdx][colIdx+2] = middle_row_coeff_3;
			hp_matrix[rowIdx][colIdx+3] = middle_row_coeff_4;
			hp_matrix[rowIdx][colIdx+4] = middle_row_coeff_5;
			
			colIdx++;
			rowIdx++;			
		}
			
		return hp_matrix;
		
	}
	
	
	public static double [][] get_hp_trend(){
		return trend;
	}
	
	
	public static double [][] get_hp_cycle(){
		return cycle;
	}
	
	
	public static double get_hp_penalty(){
		return penalty;
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_hp_results(){
		
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
	 	
	 	String title1  = "Hodrick-Prescott Filter";
	 	String subTitle1 = "Trend Component";
	 	
	 	String title2  = "Hodrick-Prescott Filter";
	 	String subTitle2 = "Cyclical Component";
	 	
	 	String [] title = {title1, title2};
	 	String [] subTitle = {subTitle1, subTitle2};
	 	
	 	String [] yLabel = {"Observed/HP-Trend", "HP-Cycle"};
	 		
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
		
		HodrickPrescott hpFilter = new HodrickPrescott(obsData,1600.0);
		
		hpFilter.plot_hp_results();
		
	}
	
}
