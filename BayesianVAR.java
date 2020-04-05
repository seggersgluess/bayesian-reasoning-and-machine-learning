package TimeSeriesAnalysis;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;

import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class BayesianVAR {

	double [][] observed_variables;
	
	int n_observations;
	int n_usedObservations;
	int n_variables; 
	
	int lag;
	int startIdx;
	int endIdx;
	
	String [][] variable_names;
	
	double [][] const_vec;
	ArrayList<double [][]> ar_matrices;
	double [][] sigma_matrix;
	
	
	double lambda = 0.0;
	double epsilon = 1e-06;
	double [][] deltas;
	
	double [][] fittedValues;
	double [][] residuals;
	
	public BayesianVAR(double [][] obs_variables, int start_idx, int end_idx, int lag_number, double lambda){
		
		if(start_idx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(end_idx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(start_idx >= end_idx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(lag_number < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		this.lambda = lambda;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx = start_idx;
		endIdx   = end_idx-1;
		lag      = lag_number;
	
		n_usedObservations = endIdx-startIdx+1;
		
		check_inputs_for_constructor();

	}
	
	
	public BayesianVAR(double [][] obs_variables, int lag_number, double lambda){
		
		if(lag_number < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		this.lambda = lambda;
		
		lag = lag_number;
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = lag;
		endIdx                 = obs_variables.length-1;
		
		n_usedObservations = endIdx-startIdx+1;
					
	}
	
	
	public void estimate_bayesianVAR() {
		
		set_initial_sigmas_from_ar();
		
		HashMap<String, double [][]> data_struct = calc_data_structure_with_dummies();
		
		double [][] X = data_struct.get("X");
		double [][] Y = data_struct.get("Y");
		
		int k = n_variables*lag+1;
		
		double [][] X_trans = MatrixOperations.transpose(X);
		double [][] X_inv = new double [k][k];
		
		if(k<100) {
			X_inv = MatrixOperations.inverse(MatrixOperations.multiplication(X_trans, X));
		}else {
			X_inv = MatrixOperations.inverse_fast(MatrixOperations.multiplication_fast(X_trans, X));
		}
		
		double [][] B = MatrixOperations.multiplication(X_inv, X_trans);
		B = MatrixOperations.multiplication(B, Y);
		double [][] res = MatrixOperations.substract(Y, MatrixOperations.multiplication(X,B));
		sigma_matrix = MatrixOperations.multiplication(MatrixOperations.transpose(res), res);
		
		X = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(X, 0, n_usedObservations-1, 0, X[0].length-1);
		Y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(Y, 0, n_usedObservations-1, 0, Y[0].length-1);
		fittedValues = MatrixOperations.multiplication(X, B);
		residuals = MatrixOperations.substract(Y,fittedValues);
	}
	
	
	@SuppressWarnings("static-access")
	public void set_initial_sigmas_from_ar() {
		sigma_matrix = new double [n_variables][n_variables];
		for(int i=0; i<n_variables; i++) {
			double [][] series = MatrixOperations.get_column_vec_from_matrix(observed_variables, i);
			AR obj_ar = new AR(series,startIdx, endIdx,lag,"OLS");		
			obj_ar.do_ar_estimation();
			sigma_matrix[i][i] = obj_ar.get_est_sigma();
		}
	}
	
	
	public HashMap<String, double [][]> calc_data_structure_with_dummies() {
		
		HashMap<String, double [][]> data_struct = new HashMap<String, double [][]>();
		
		int k = n_variables*lag+1;
		int n_rows = n_usedObservations+k+n_variables;
		
		double [][] Y = new double [n_rows][n_variables];
		double [][] X = new double [n_rows][k];
		
		double [][] J_p = new double [lag][lag];
		double [][] scaled_sigma_diag = new double [n_variables][n_variables];
		
		//Observed variables Y,X
		for(int t=0; t<n_usedObservations; t++) {
			int curIdx = startIdx+t;
			int idx = 0;
			for(int i=0; i<n_variables; i++) {
				Y[t][i] = observed_variables[curIdx][i];		
			}
			for(int p=0; p<lag; p++) {
				int lagIdx = curIdx-p-1;
				for(int i=0; i<n_variables; i++) {
					X[t][idx] = observed_variables[lagIdx][i];
					idx++;
				}
				J_p[p][p] = p+1.0;
			}			
		}
		
		//Dummy variables Y_d, X_d
		int startIdx1 = n_usedObservations;
		int startIdx2 = startIdx1+n_variables*lag;
		
		for(int i=0; i<n_variables; i++) {
			double sigma = Math.sqrt(sigma_matrix[i][i]);
			Y[startIdx1+i][i] = sigma*deltas[i][0]/lambda;
			Y[startIdx2+i][i] = sigma;		
			scaled_sigma_diag[i][i] = sigma/lambda;
		}
		
		double [][] kron_prod = MatrixOperations.kronecker(J_p, scaled_sigma_diag);
		
		int startRowIdx = n_usedObservations;
		int endRowIdx = startRowIdx+lag*n_variables-1;
		int startColIdx = 0;
		int endColIdx = lag*n_variables-1;
		
		X = MatrixOperations.set_sub_matrix_to_matrix(X, kron_prod, startRowIdx, endRowIdx, startColIdx, endColIdx);
		X[n_rows-1][k-1] = epsilon;
		
		data_struct.put("Y", Y);
		data_struct.put("X", X);
		
		return data_struct;
	}
	
	
	public void set_est_pars_from_B(double [][] B) {
		
		ar_matrices = new ArrayList<double [][]>();
		
		int startColIdx = 0;
		int endColIdx = n_variables-1;
		
		int startRowIdx = 0;
		int endRowIdx = 0;
		
		for(int p=0; p<lag; p++) {
			endRowIdx = startRowIdx+n_variables-1;
			double [][] ar_matrix = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(B, startRowIdx, endRowIdx, startColIdx, endColIdx);
			startRowIdx = endRowIdx+1;
			ar_matrices.add(ar_matrix);
		}
		
		const_vec = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(B, startRowIdx+1, startRowIdx+1, startColIdx, endColIdx);		
	}
	
	
	public void check_inputs_for_constructor(){
		
		if(observed_variables.length == 0){
			throw new RuntimeException("No input data loaded yet.");		
		}
		
		int totalNumberOfUsedObs = n_usedObservations+lag;
		
		if(totalNumberOfUsedObs > n_observations){
			throw new RuntimeException("Not enough data for estimating model.");
		}
		
	}
	
	
	public double [][] get_sigma_matrix() {
		if(sigma_matrix == null) {
			throw new RuntimeException("Bayesian VAR not yet estimated. No covariance matrix found.");
		}
		return sigma_matrix;
	}
	
	
	public double [][] get_const_vec() {
		if(const_vec == null) {
			throw new RuntimeException("Bayesian VAR not yet estimated. No vector of constants found.");
		}
		return const_vec;
	}
	
	
	//lag = 1,2,... (t-1, t-2,...)
	public double [][] get_ar_coeff_matrix(int lag) {
		if(ar_matrices == null) {
			throw new RuntimeException("Bayesian VAR not yet estimated. No AR-matrix found.");
		}
		if(lag<=0) {
			throw new RuntimeException(lag + " is not a vaild lag. Only integers > 0 allowed.");
		}
		if(this.lag<=lag) {
			throw new RuntimeException(lag + " is not a vaild lag. Lag of estimation is: " + this.lag);
		}
		
		return ar_matrices.get(lag-1);
	}
	
	
	public void set_deltas(double [][] deltas) {
		this.deltas = deltas;
	}
	
	
	public void set_variable_names(String [][] labels) {
		if(labels.length != n_variables) {
			System.out.println("Invalid number of variable names supplied.");
		}else {
			variable_names = labels;
		}
	}
	
	
	@SuppressWarnings("static-access")
	public void plot_time_series_and_fitted_values() {
		
		if(fittedValues == null) {
			throw new RuntimeException("No fitted data found for bayesian VAR. First fit the bayesian VAR.");
		}
		
        GenGraphics obj_graph = new GenGraphics();
		       
        obj_graph.setNumberOfPlotColums(8);
        obj_graph.setNumberOfPlotRows(6);
	 	
        obj_graph.setGraphWidth(1100);
        obj_graph.setGraphHeight(700);

	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	String [] titles = new String [n_variables];
	 	String [] yLabels = new String [n_variables];
	 	
	 	for(int i=0; i<n_variables; i++) {
		 	double [][] obs_data    = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, i, i);
		 	double [][] fitted_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(fittedValues, 0, n_usedObservations-1, i, i);	
		 	obj_graph.plotLines(xAxis, fitted_data, true, Color.RED);
		 	obj_graph.plotPoints(xAxis, obs_data, false, Color.BLUE);
		 	if(variable_names == null) {
		 		titles[i]  = "Observed vs. Fitted Variables";
		 	}else {
		 		titles[i] = variable_names[i][0];
		 	}
		 	yLabels[i] = "";
	 	}
	 	
	 	obj_graph.setTitle(titles, null, "9");
	 	obj_graph.setYLabel(yLabels, null, "9");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 8);
	 	obj_graph.setFontOfYAxisUnits("bold", 8);
	 	obj_graph.setNumberOfXDivisions(4);
	 	obj_graph.setNumberOfYDivisions(3);
	 	obj_graph.set_point_width(2);
	 	obj_graph.set_line_widht(2);
	 	
	 	obj_graph.plot();
		
	}
	
	
}
