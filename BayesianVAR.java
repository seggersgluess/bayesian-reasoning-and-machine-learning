package TimeSeriesAnalysis;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Distributions.InvWishart;
import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.Cholesky;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

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
	
	//Posterior infos
	int df = 0;
	double [][] inv_XX;
	
	//Impulse Response
	ArrayList<double [][]> wold_ma_rep;
	int [] impulse_vars;
	int [] response_vars;
	
	ArrayList<ArrayList<List<Double>>> orthogonal_ir;
	int ir_steps = 0;
	double [][] s_shocks;
	int n_draws = 10;
	
	ArrayList<Double> confidence = new ArrayList<Double>();
	
	ArrayList<ArrayList<ArrayList<double [][]>>> upper_conf_orthogonal_ir;
	ArrayList<ArrayList<ArrayList<double [][]>>> lower_conf_orthogonal_ir;
	
	
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
		confidence.add(0.9);
		
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
		confidence.add(0.9);			
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
		set_est_pars_from_B(B);
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
		
		df = n_rows+2-k;
		
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
			X[t][k-1] = 1.0;
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
	
	
	private void set_est_pars_from_B(double [][] B) {
		
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
		
		const_vec = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(B, startRowIdx, startRowIdx, startColIdx, endColIdx);		
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
	
	
	//lag = 1,2,... (t-1, t-2,...) -> Starting with p=1!
	public double [][] get_ar_matrix(int lag) {
		if(ar_matrices == null) {
			throw new RuntimeException("Bayesian VAR not yet estimated. No AR-matrix found.");
		}
		if(lag<=0) {
			throw new RuntimeException(lag + " is not a vaild lag. Only integers > 0 allowed.");
		}
		if(this.lag<lag) {
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
	
	
	public void calc_Wold_MA_representation(int n_steps){
		
		n_steps++;
		wold_ma_rep = new ArrayList<double [][]>(n_steps);
		
		double [][] psi = MatrixOperations.identity(n_variables);
		double [][] psi_new = new double [n_variables][n_variables];
		
		wold_ma_rep.add(psi);
		
		for(int i=1; i<n_steps; i++){		
			psi_new = new double [n_variables][n_variables];
			for(int j=0; j<i; j++){					
				if(j>(lag-1)){
					break;					
				}else{
					int idx = i-j-1;
					double [][] ar_matrix = get_ar_matrix(j+1);
					psi = wold_ma_rep.get(idx);
					psi_new = MatrixOperations.add(psi_new, MatrixOperations.multiplication(ar_matrix, psi));
				}	
			}			
			wold_ma_rep.add(psi_new);
		}		
	}
	
	
	public double [][] get_wold_ma_rep_matrix(int step){
		
		if(step > wold_ma_rep.size()){
			throw new RuntimeException("Supplied step " + step + " not allowed.");
		}		
		double [][] wold_ma_matrix = wold_ma_rep.get(step);
		
		return wold_ma_matrix;
	}
	
	
	//TODO: specification of s_shock: double [n_variables][impulse.length]-> In columns spec. of shock from impulse variable
	public void calc_orthogonal_ir(int [] impulse, int [] response, int n_steps, double [][] s_shocks){
		
		if(s_shocks != null) {
			if(s_shocks.length != n_variables) {
				throw new RuntimeException("Invalid number of supplied shocks. Number of rows has to be " + n_variables);
			}
			if(s_shocks[0].length != impulse.length) {
				throw new RuntimeException("Invalid number of supplied shocks. Supply for every impulse variable a shock.");
			}
			
			for(int i=0; i<n_variables; i++) {
				for(int j=0; j<impulse.length; j++) {
					if(s_shocks[i][j] != 0.0) {
						if(i!=impulse[j]) {
							throw new RuntimeException("Shock supplied for variable " + i + " instead for impulse variable " + impulse[j]);
						}
					}
				}
			}		
			this.s_shocks = s_shocks;
		}
		
		ir_steps = n_steps;
		
		calc_Wold_MA_representation(n_steps+1);
		
		double [][] P = Cholesky.decompose(get_sigma_matrix());
		
		int n_impulses = impulse.length;
		int n_responses = response.length;
		
		impulse_vars = impulse;
		response_vars = response;
		
		orthogonal_ir = new ArrayList<ArrayList<List<Double>>>(n_impulses);
		
		for(int i=0; i<n_impulses; i++){
			ArrayList<List<Double>> responseList = new ArrayList<List<Double>>(n_responses);
			double [][] ortho_shock = new double [n_variables][1];
			if(s_shocks == null) {
				ortho_shock = MatrixOperations.get_column_vec_from_matrix(P, impulse[i]);
			}else {
				ortho_shock = s_shocks;
			}			
			for(int j=0; j<n_responses; j++){
				ArrayList<Double> responses = new ArrayList<Double>(n_steps);			
				for(int t=0; t<n_steps; t++){
					double [][] wold_ma_coeff = get_wold_ma_rep_matrix(t);
					wold_ma_coeff = MatrixOperations.multiplication(wold_ma_coeff, ortho_shock);
					responses.add(wold_ma_coeff[response[j]][0]);
				}
				responseList.add(responses);
			}	
			orthogonal_ir.add(responseList);
		}		
	}	
	
	
	public double [][] get_orthogonal_ir(int impulse, int response){
		
		int [] impIdx = Utilities.get_idx(impulse_vars, impulse);
		
		if(impIdx[0] == -1){
			throw new RuntimeException("Invalid impulse variable. Calculate IR function for this variable.");
		}
		
		int [] respIdx = Utilities.get_idx(response_vars, response);
		
		if(respIdx[0] == -1){
			throw new RuntimeException("Invalid response variable. Calculate IR function for this variable.");
		}
		
		int n_steps = orthogonal_ir.get(impIdx[0]).get(respIdx[0]).size();
		
		double [][] ir = new double [n_steps][1];
		
		for(int t=0; t<n_steps; t++){
			ir[t][0] = orthogonal_ir.get(impIdx[0]).get(respIdx[0]).get(t);
		}
		
		return ir;		
	}
	
	
	public void calc_ir_confidence_intervals() {
		
		if(confidence == null) {
			confidence.add(0.9);
		}
		
		upper_conf_orthogonal_ir = ir_conf_bands_struct();
		lower_conf_orthogonal_ir = ir_conf_bands_struct();
		
		InvWishart invW = new InvWishart(sigma_matrix, df);
		double [][] B_vec = MatrixOperations.vecAs2dArray(get_B_matrix());
		
		if(inv_XX == null) {
			//TODO: generate inv_XX!
		}
		
		//TODO: Check if this is constant!
		ArrayList<double[][]> org_ar_matrices = ar_matrices;
		double [][] org_const_vec = const_vec;
		ArrayList<double [][]> org_wold_ma_rep = wold_ma_rep;
		ArrayList<ArrayList<List<Double>>> org_orthogonal_ir = orthogonal_ir;
		ArrayList<ArrayList<double [][]>> sim_ir = sim_ir_struct();
		
		int n_imp_vars = impulse_vars.length;
		int n_resp_vars = response_vars.length;
		
		for(int i=0; i<n_draws; i++) {
			double [][] psi = invW.sample();
			double [][] parsCov = MatrixOperations.kronecker(psi, inv_XX);
			NormalDistribution normDist = new NormalDistribution(B_vec, parsCov);
			double [][] B = normDist.sample();
			set_est_pars_from_B(B);
			calc_orthogonal_ir(impulse_vars, response_vars, ir_steps, s_shocks);
			for(int j=0; j<n_imp_vars; j++) {
				for(int k=0; k<n_resp_vars; k++) {
					for(int s=0; s<ir_steps; s++) {
						//TODO: Check if values are set to sim_ir!
						sim_ir.get(j).get(k)[i][s] = orthogonal_ir.get(j).get(k).get(s);
					}					
				}
			}
		}
		
		int n_conf_levels = confidence.size();
		
		for(int i=0; i<n_imp_vars; i++) {
			for(int j=0; j<n_resp_vars; j++) {
				for(int s=0; s<ir_steps; s++) {
					double [] bootIRcoeff = MatrixOperations.get_column_from_matrix(sim_ir.get(i).get(j), s);
					for(int q=0; q<n_conf_levels; q++) {
						double uConf = confidence.get(q);
						//TODO: Check if values are set to conf bands struct!
						upper_conf_orthogonal_ir.get(i).get(j).get(q)[0][s] = GeneralMath.quantile(bootIRcoeff,uConf);
						double lConf = 1.0-confidence.get(q);
						lower_conf_orthogonal_ir.get(i).get(j).get(q)[0][s] = GeneralMath.quantile(bootIRcoeff,lConf);
					}
				}			
			}
		}
			
		ar_matrices = org_ar_matrices;
		const_vec = org_const_vec;
		wold_ma_rep = org_wold_ma_rep;
		orthogonal_ir = org_orthogonal_ir;	
	}
	
	
	public ArrayList<double [][]> get_orthogonal_ir_confidence(int impulse, int response){
		
		if(lower_conf_orthogonal_ir == null || upper_conf_orthogonal_ir == null) {
			return null;
		}
		
		int [] impIdx = Utilities.get_idx(impulse_vars, impulse);
		
		if(impIdx[0] == -1){
			throw new RuntimeException("Invalid impulse variable. Calculate IR function for this variable.");
		}
		
		int [] respIdx = Utilities.get_idx(response_vars, response);
		
		if(respIdx[0] == -1){
			throw new RuntimeException("Invalid response variable. Calculate IR function for this variable.");
		}
		
		int n_steps = orthogonal_ir.get(impIdx[0]).get(respIdx[0]).size();
		
		ArrayList<double [][]> ir_conf_bands = new ArrayList<double [][]>();
		int n_conf_levels = confidence.size();
		for(int q=0; q<n_conf_levels; q++) {
			double [][] conf_bands = new double [n_steps][2];
			for(int t=0; t<n_steps; t++){
				conf_bands[t][0] = lower_conf_orthogonal_ir.get(impIdx[0]).get(respIdx[0]).get(q)[0][t];
				conf_bands[t][1] = upper_conf_orthogonal_ir.get(impIdx[0]).get(respIdx[0]).get(q)[0][t];
			}
			ir_conf_bands.add(conf_bands);
		}
		
		return ir_conf_bands;
	}
	
	
	private ArrayList<ArrayList<double [][]>> sim_ir_struct() {
		
		ArrayList<ArrayList<double [][]>> sim_ir_struct = new ArrayList<ArrayList<double [][]>>();
		
		int n_imp_vars = impulse_vars.length;
		int n_resp_vars = response_vars.length;
		
		for(int i=0; i<n_imp_vars; i++) {
			for(int j=0; j<n_resp_vars; j++) {
				ArrayList<double [][]> resp_struct = new ArrayList<double [][]>();
				resp_struct.add(new double [n_draws][ir_steps]);
				sim_ir_struct.add(resp_struct);
			}		
		}
		return sim_ir_struct;
	}
	
	
	private ArrayList<ArrayList<ArrayList<double [][]>>> ir_conf_bands_struct() {
		
		ArrayList<ArrayList<ArrayList<double [][]>>> ir_conf_bands_struct = new ArrayList<ArrayList<ArrayList<double [][]>>>();
		
		int n_imp_vars = impulse_vars.length;
		int n_resp_vars = response_vars.length;
		int n_conf_levels = confidence.size();
		
		for(int i=0; i<n_imp_vars; i++) {
			for(int j=0; j<n_resp_vars; j++) {
				ArrayList<ArrayList<double [][]>> resp_struct = new ArrayList<ArrayList<double [][]>>();
				for(int q=0; q<n_conf_levels; q++) {
					ArrayList<double [][]> conf_struct = new ArrayList<double [][]>();
					conf_struct.add(new double [1][ir_steps]);
					resp_struct.add(conf_struct);
				}
				ir_conf_bands_struct.add(resp_struct);
			}		
		}
		return ir_conf_bands_struct;
	}
	
	
	private double [][] get_B_matrix() {
		double [][] B = new double [n_variables*lag+1][n_variables];
		
		int startColIdx = 0;
		int endColIdx = n_variables-1;		
		int startRowIdx = 0;
		int endRowIdx = 0;
		
		for(int p=0; p<lag; p++) {
			endRowIdx = startRowIdx+n_variables-1;
			B = MatrixOperations.set_sub_matrix_to_matrix(B, ar_matrices.get(p+1),startRowIdx, endRowIdx, startColIdx, endColIdx);
			startRowIdx = endRowIdx+1;
		}
		
		B = MatrixOperations.set_sub_matrix_to_matrix(B, const_vec, startRowIdx, startRowIdx, startColIdx, endColIdx);		

		return B;
	}
	
	
	public int [] get_idxs_of_variables(String [] variables) {
		
		if(variable_names == null) {
			throw new RuntimeException("No variable names found. Estimate first the bayesian VAR.");
		}
		
		int n_vars = variables.length;
		
		int [] idxs = new int [n_vars];
		for(int i=0; i<n_vars; i++) {
			int [] idx = Utilities.get_idx(variable_names, variables[i]);
			if(idx[0] == -1) {
				throw new RuntimeException(variables[i] + " not found in supplied data set.");
			}
			idxs[i] = idx[0];
		}
		return idxs;
	}
	
	
	//Deletes large matrix inv_XX necessary for inference from posterior (and IR)
	public void freeMemory() {
		inv_XX = null;
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
	
	
	@SuppressWarnings("static-access")
	public void plotImpulseResponse(int impulse, int response, int n_steps){
   		
		double [][] irData = new double [n_steps][1];
		
		double [][] ir = get_orthogonal_ir(impulse, response);
		
		if(ir.length<n_steps) {
			throw new RuntimeException("Invalid number of steps. Only" + ir.length + " steps calculated.");
		}
		
		ArrayList<double[][]> confBands = get_orthogonal_ir_confidence(impulse, response);
		int n_conf_bands = confidence.size();
		
        for(int i=0; i<n_steps; i++){
        	irData[i][0] = ir[i][0];
        }
		        
        GenGraphics obj_graph = new GenGraphics();
		        
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_steps][1];
	 	double [][] xAxis4IR = new double [n_steps][2];
	 	
	 	for(int i=0; i<n_steps; i++){
	 		xAxis[i][0] = i+1;
	 		xAxis4IR[i][0] = i+1;
	 		xAxis4IR[i][1] = i+1;
	 	}
	 	
	 	String [] title = new String [1];
	 	String [] yLabel = new String [1];
	 	  
	 	obj_graph.plotLines(xAxis,irData,true);
		
	 	if(confBands != null) {
		 	//TODO: Check this!
		 	for(int i=0; i<n_conf_bands; i++) {
		 		obj_graph.plotLines(xAxis4IR,confBands.get(i),false,Color.GRAY);
		 	}
	 	}
	
	 	List<Color> lineColor = new ArrayList<Color>();
	 	lineColor.add(Color.RED);
	 	lineColor.add(Color.BLACK);
	 	lineColor.add(Color.BLACK);
	 	
	 	title[0]  = "Impulse: " + impulse + " Response:" + response;
	 	yLabel[0] = title[0];
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);

	 	obj_graph.setLineColor(lineColor);
	 	
	 	obj_graph.plot();
   		
   	}
	
	
}
