package Kernels;

import java.awt.Color;
import java.util.HashMap;
import java.util.Set;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;
import Regression.LinearRegression;

public class TestKernelRidgeRegression {

public static void test1KernelRidgeRegression() {
		
    	double[][] data = null;
    	
    	try {
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
    		data = inputData.selectedDblFileData;
    		
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
		
    	int lag = 20;
    	int n_usedObs = data.length-lag;
    	double [][] laggedData = new double [n_usedObs][lag];
    	double [][] unlaggedData = new double [n_usedObs][1];
    	
    	int startIdx = lag-1;
    	
        for(int i=0; i<n_usedObs; i++) {
        	unlaggedData[i][0] = data[startIdx+i+1][0];
        	for(int l=0; l<lag; l++) {
        		laggedData[i][l] = data[startIdx-l+i][0];
        	}	
        }
    	
        double [][] y = MatrixOperations.get_sub_matrix_between_column_idxs(unlaggedData, 0, 0);
    	double [][] X = MatrixOperations.get_sub_matrix_between_column_idxs(laggedData, 1, lag-1); 
    	
    	LinearRegression obj_lm = new LinearRegression(y, X, true);
    	long startTime = System.currentTimeMillis();
    	obj_lm.do_parameter_estimation();
    	double [][] reg_est_values = obj_lm.get_fitted_values();
    	long endTime = System.currentTimeMillis();
    	System.out.println("Finished AR-regression after " + ((endTime-startTime)/1000.0) + " secs.");
    	
    	//--- Start kernel ridge regressions ---
    	startTime = System.currentTimeMillis();
    	
    	double penalty = 0.1;
    	
    	HashMap<String, HashMap<String,Object[][]>> kernel_ridge_est_res = new HashMap<String, HashMap<String,Object[][]>>();
    	
    	HashMap<String,Object[][]> par_str = new HashMap<String,Object[][]>();	
    	HashMap<String,Object[][]> list_fitted_values = new HashMap<String,Object[][]>();
    	HashMap<String,Object[][]> list_residuals = new HashMap<String,Object[][]>();
    	
    	KernelRidgeRegression obj_ridge = new KernelRidgeRegression(y, X, penalty, true, "linear");
    	obj_ridge.fit();
    	obj_ridge.predict(X);
    	
    	double [][] pred_Values = obj_ridge.get_predictedValues();   
    	    	
		String [][] pars = new String [1][1];
		pars[0][0] = "Kernel: LINEAR, Penalty: " + penalty;
			
		par_str.put("LINEAR", pars);
    	
    	double [][] prim_fittedValues = obj_ridge.get_fitted_values();   	 	
    	double [][] prim_residuals = obj_ridge.get_residuals();
    	
    	int n_obs = prim_fittedValues.length;
    	Double [][] fittedValues = new Double [n_obs][1];
    	Double [][] residuals = new Double [n_obs][1]; 
    	
    	for(int j=0; j<n_obs; j++) {
    		fittedValues[j][0] = prim_fittedValues[j][0];
    		residuals[j][0]    = prim_residuals[j][0];
    	}
    	
    	list_fitted_values.put("LINEAR", fittedValues);
    	list_residuals.put("LINEAR",residuals);
		
    	//Check prediction function for kernel = LINEAR
    	checkPredictionsFromPredFunction(prim_fittedValues, pred_Values);
		
    	String [] kernels = obj_ridge.get_valid_kernels();
    	int n_kernels = kernels.length;
    	
    	for(int i=0; i<n_kernels; i++) {	
    		String k = kernels[i].toUpperCase();
    		System.out.println("---Do " + k + " kernel ridge fitting---");
    		if(kernels[i].contentEquals("rbf") == true) {
    			
    			double bandwidth = 2.0;
    			pars = new String [1][1];
    			pars[0][0] = "Penalty: " + penalty + ", Bandwidth: " + bandwidth;
    				
    			par_str.put(k, pars);
    			
        		obj_ridge = new KernelRidgeRegression(y, X, penalty, true, "rbf");
        		obj_ridge.fit(bandwidth);
        		obj_ridge.predict(X, bandwidth);
            	pred_Values = obj_ridge.get_predictedValues(); 
        	}
        	
        	if(kernels[i].contentEquals("sigmoid") == true) {
        		
        		double gamma = 0.001;
        		double alpha = 0.001;
        		
        		pars = new String [1][1];
    			pars[0][0] = "Penalty: " + penalty + ", Alpha: " + alpha + ", Gamma: " + gamma;
    				
    			par_str.put(k, pars);
        		
    			obj_ridge = new KernelRidgeRegression(y, X, penalty, true, "sigmoid");
    			obj_ridge.fit(gamma, alpha);
            	obj_ridge.predict(X, gamma, alpha);
            	pred_Values = obj_ridge.get_predictedValues(); 
        	}
        	
        	if(kernels[i].contentEquals("polynomial") == true) {
        		
        		double gamma = 0.1;
        		double alpha = 0.1;
        		double degree = 4.0;
        		
        		pars = new String [1][1];
    			pars[0][0] = "Penalty: " + penalty + ", Alpha: " + alpha + ", Gamma: " + gamma + ", Degree: " + degree;
    				
    			par_str.put(k, pars);
        		
    			obj_ridge = new KernelRidgeRegression(y, X, penalty, true, "polynomial");
        		obj_ridge.fit(degree, gamma, alpha);
            	obj_ridge.predict(X, degree, gamma, alpha);
            	pred_Values = obj_ridge.get_predictedValues(); 
        	}
        	
        	if(kernels[i].contentEquals("se") == true) {
        		
        		pars = new String [1][1];
    			pars[0][0] = "Penalty: " + penalty + ", Sigma: I";
    				
    			par_str.put(k, pars);
        		
    			obj_ridge = new KernelRidgeRegression(y, X, penalty, true, "se");
        		//--- Use identity matrix as covariance matrix ---
        		double [][] Sigma = MatrixOperations.identity(X[0].length);
        		obj_ridge.fit(Sigma);
            	obj_ridge.predict(X, Sigma);
            	pred_Values = obj_ridge.get_predictedValues();
        	}
        	
        	prim_fittedValues = obj_ridge.get_fitted_values();
        	prim_residuals = obj_ridge.get_residuals();
        	
        	fittedValues = new Double [n_obs][1];
        	residuals = new Double [n_obs][1]; 
        	
        	for(int j=0; j<n_obs; j++) {
        		fittedValues[j][0] = prim_fittedValues[j][0];
        		residuals[j][0]    = prim_residuals[j][0];
        	}
        	
        	list_fitted_values.put(k, fittedValues);
        	list_residuals.put(k,residuals);
        	
        	//Check prediction function for respective kernel
        	checkPredictionsFromPredFunction(prim_fittedValues, pred_Values);
        	
    	}
    	
    	kernel_ridge_est_res.put("Parameters", par_str);
    	kernel_ridge_est_res.put("FittedValues", list_fitted_values);
    	kernel_ridge_est_res.put("Residuals", list_residuals);
    	
    	endTime = System.currentTimeMillis();
    	    	
    	System.out.println("Finished " + n_kernels + " AR kernel ridge regression after " + ((endTime-startTime)/1000.0) + " secs.");
    	
    	HashMap<String, double[][]> linear_est_res = new HashMap<String, double [][]>(3);
    	linear_est_res.put("Observed", y);
    	linear_est_res.put("Linear", reg_est_values);
    	linear_est_res.put("Residuals_linear", obj_lm.get_residuals());

    	
    	plot4KernelRidgeRegressionTest1(kernel_ridge_est_res, linear_est_res);
	}


	public static void checkPredictionsFromPredFunction(double [][] fitFunctionValues, double [][] predFunctionValues) {
		int n = fitFunctionValues.length;
		for(int i=0; i<n; i++) {
			double d = fitFunctionValues[i][0]-predFunctionValues[i][0];
			if(d != 0.0) {
				throw new RuntimeException("Test failed. Difference not 0!");
			}
		}
	}


	@SuppressWarnings("static-access")
	public static void plot4KernelRidgeRegressionTest1(HashMap<String, HashMap<String,Object[][]>> kernel_ridge_est_res, HashMap<String, double [][]> linear_est_res) {
		
		Set<String> kernelSet = kernel_ridge_est_res.get("Parameters").keySet();
		int n_kernels = kernelSet.size();
		
		String [] kernels = new String [n_kernels];
		int idx = 0;
		for (String kernel : kernelSet) {
            kernels[idx++] = kernel; 
		}
		
		double [][] y_obs          = linear_est_res.get("Observed");
		double [][] y_linear       = linear_est_res.get("Linear");
		double [][] residuals_linear = linear_est_res.get("Residuals_linear");
	
		int n_obs = y_obs.length;
		
		double [][] x = new double [n_obs][1];
		
		for(int i=0; i<n_obs; i++) {
			x[i][0] = i+1;
		}
		
		GenGraphics graph = new GenGraphics();
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setNumberOfPlotColums(2);
	 	graph.setNumberOfPlotRows(n_kernels);
	 	
	 	graph.setGraphWidth(650);
	 	graph.setGraphHeight(600);
	 	
	 	String [] titles    = new String [2*n_kernels];
	 	String [] subTitles = new String [2*n_kernels];
	 	String [] xLabel    = new String [2*n_kernels];
	 	String [] yLabel    = new String [2*n_kernels];
	 	
	 	idx = 0;
	 	
	 	for(int k=0; k<n_kernels; k++) {
	 		
	 		String kernel = kernels[k];
	 		
	 		titles[idx] = kernel + " Kernel implied Fitted Values";
	 		subTitles[idx] = (String) kernel_ridge_est_res.get("Parameters").get(kernel)[0][0];
	 		xLabel[idx] = "Observation No.";
	 		yLabel[idx] = "Values";
	 		idx++;
	 		titles[idx] = kernel + " Kernel implied Residuals";
	 		subTitles[idx] = (String) kernel_ridge_est_res.get("Parameters").get(kernel)[0][0];
	 		xLabel[idx] = "Observation No.";
	 		yLabel[idx] = "Residuals";
	 		idx++;
	 		
	 		double [][] kernel_fit = new double [n_obs][1];
	 		double [][] kernel_res = new double [n_obs][1];
	 		
	 		for(int i=0; i<n_obs; i++) {
	 			kernel_fit[i][0] = (double) kernel_ridge_est_res.get("FittedValues").get(kernel)[i][0];
	 			kernel_res[i][0] = (double) kernel_ridge_est_res.get("Residuals").get(kernel)[i][0];
	 		}
	 		
		 	graph.plotLines(x, y_linear ,true, Color.RED);
		 	graph.plotLines(x, kernel_fit ,false, Color.BLUE);
		 	graph.plotPoints(x, y_obs, false, Color.BLACK);
		 	graph.set_point_width(3);
		 	
		 	graph.plotLines(x, residuals_linear ,true, Color.RED);
		 	graph.plotLines(x, kernel_res ,false, Color.BLUE);
	 		
	 	}
	 		 	 	
	 	graph.setTitle(titles, null, "10");
	 	graph.setSubTitle1(subTitles, null, "9");
	 	graph.setYLabel(yLabel, null, "7");
	 	graph.setXLabel(xLabel, null, "7");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(2);
	 	graph.setFontOfXAxisUnits("plain", 7);
	 	graph.setFontOfYAxisUnits("plain", 7);
	 	graph.setNumberOfXDivisions(5);
	 	
	 	graph.plot();
		
	}


	//Calculation of Gram matrix from specific kernel
	public static void test2KernelRidgeRegression() {
	
		double [][] X = new double [3][3];
		
		double counter = 0.0;
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				counter++;
				X[i][j] = counter;
			}
		}
		
		MatrixOperations.print_matrix(X);
		
		KernelFunctions kf = new KernelFunctions("linear");
		double [][] g = kf.calc_gram_matrix(X);
		
		System.out.println("Gram matrix:");
		MatrixOperations.print_matrix(g);
		
	}

	
	public static void main(String[] args) {
		test1KernelRidgeRegression();
		//test2KernelRidgeRegression();
	}

}
