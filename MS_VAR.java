package HiddenMarkovModels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Regression.LinearRegression;
import Utilities.Utilities;

public class MS_VAR extends HMM{

	//constructor
	public MS_VAR(int startIdx4Obs, int endIdx4Obs, int lagNumber, int stateNumber, String ms_type){
		
		if(startIdx4Obs < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(endIdx4Obs < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(startIdx4Obs >= endIdx4Obs){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(lagNumber < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		if(stateNumber < 0){
			throw new RuntimeException("No valid number of states supplied.");
		}
		
		startIdx = startIdx4Obs;
		endIdx   = endIdx4Obs;
		lag      = lagNumber;
		n_states = stateNumber;
	
		n_usedObservations = endIdx-startIdx+1;
		
		check_inputs_for_constructor();
		set_MS_VAR_type(ms_type);
		set_default_values_MS_VAR();
		
	}
	
	
	public static void check_inputs_for_constructor(){
		
		if(observed_variables.length == 0){
			throw new RuntimeException("No input data loaded yet.");		
		}
		
		int totalNumberOfUsedObs = n_usedObservations+lag;
		
		if(totalNumberOfUsedObs > n_observations){
			throw new RuntimeException("Not enough data for estimating model.");
		}
		
	}
	
	
	//Updates transition matrix P (see e.g. Droumaguet 2012 pp. 12)
	public static void update_transition_matrix_4_MS_VAR(){
		
		int nRows = (int) Math.pow(n_states, 2.0);
		int usedDates = n_usedObservations-1;
		
		double [][] filterProbs = get_filtered_probs();
		double [][] smoothProbs = get_smoothed_probs();
		double [][] transposedTrans = MatrixOperations.transpose(transMatrix);
		
		double [][] xi_1 = new double [nRows][1];
		double [][] xi_2 = new double [nRows][1];
		double [][] e    = MatrixOperations.unit_vector(n_states);
		
		for(int t=0; t<usedDates; t++){
			
			double [][] selectedFilteredProbs = MatrixOperations.get_row_vec_from_matrix(filterProbs, t);
			double [][] updatedFilteredProbs  = MatrixOperations.multiplication(transposedTrans, selectedFilteredProbs);
			double [][] selectedSmoothedProbs = MatrixOperations.get_row_vec_from_matrix(smoothProbs, t+1);
			
			double [][] addTerm = MatrixOperations.divideVecByElement(selectedSmoothedProbs, updatedFilteredProbs);
			addTerm = MatrixOperations.kronecker(addTerm, selectedFilteredProbs);
			
			xi_2 = MatrixOperations.add(xi_2, addTerm);
			
		}
		
		xi_2 = MatrixOperations.multiplyVecByElement(MatrixOperations.vecAs2dArray(transMatrix), xi_2);
		
		xi_1 = MatrixOperations.kronecker(MatrixOperations.transpose(e),MatrixOperations.identity(n_states));
		xi_1 = MatrixOperations.multiplication(xi_1, xi_2);
		
		double [][] p_update = MatrixOperations.divideVecByElement(xi_2, MatrixOperations.kronecker(e,xi_1));
		
		transMatrix = MatrixOperations.get_matrix_from_vec(p_update, n_states, n_states);
		
	}
	
	
	public static double calc_log_likelihood_4_MS_VAR(){
		
		double logLik = 0.0;
		double [][] transposedTrans = MatrixOperations.transpose(transMatrix);
		
		for(int t=0; t<n_usedObservations; t++){
			
			//Xi_t|t
			double [][] filterProbs = get_filtered_probs(t);
			//Xi_t|t-1
			filterProbs = MatrixOperations.multiplication(transposedTrans, filterProbs);
			
			double [][] obsProbs     = MatrixOperations.get_column_vec_from_matrix(obs_model_probs, t);
			obsProbs                 = MatrixOperations.transpose(obsProbs);
			
			logLik += MatrixOperations.multiplication(obsProbs, filterProbs)[0][0];
			
		}
		
		logLik = Math.log(logLik);
		
		return logLik;
		
	}
	
	
	public static boolean my_switching(){
		
		boolean is_mySwitching = false;
		
		int [] idx = Utilities.get_idx(get_my_switching_models(), ms_var_type);
		
		if(idx[0]!=-1){
			is_mySwitching = true;
		}
		
		return is_mySwitching;
		
	}
	
	
	public static boolean ar_switching(){
		
		boolean is_mySwitching = false;
		
		int [] idx = Utilities.get_idx(get_ar_switching_models(), ms_var_type);
		
		if(idx[0]!=-1){
			is_mySwitching = true;
		}
		
		return is_mySwitching;
		
	}
	
	
	public static boolean sigma_switching(){
		
		boolean is_mySwitching = false;
		
		int [] idx = Utilities.get_idx(get_sigma_switching_models(), ms_var_type);
		
		if(idx[0]!=-1){
			is_mySwitching = true;
		}
		
		return is_mySwitching;
		
	}
	
	
	public static void calc_obs_model_probs_4_MS_VAR(){
		
		boolean mySwitch    = my_switching();
		boolean arSwitch    = ar_switching();
		boolean sigmaSwitch = sigma_switching();
		
		obs_model_probs = new double [n_states][n_usedObservations];
		
		double [][] my_vec      = new double [n_variables][1];
		double [][] arMatrix    = new double [n_variables][1];
		double [][] SigmaMatrix = new double [n_variables][1];
		
		for(int m=0; m<n_states; m++){
			
			if(m==0){
				my_vec      = get_my(m);
				SigmaMatrix = get_sigma_matrix(m);
			}else{
				if(mySwitch==true){
					my_vec = get_my(m);
				}
				if(sigmaSwitch==true){
					SigmaMatrix = get_sigma_matrix(m);
				}
			}
						
			for(int t=0; t<n_usedObservations; t++){
				
				int curIdx = startIdx+t;
				
				double [][] curObs  = MatrixOperations.get_row_vec_from_matrix(observed_variables, curIdx);
				double [][] meanVec = my_vec;
						
				for(int p=0; p<lag; p++){
					
					int lagIdx = curIdx-(p+1);
					double [][] laggedObs = MatrixOperations.get_row_vec_from_matrix(observed_variables, lagIdx);
					
					if(arSwitch==true){
						arMatrix = get_ar_matrix(m,p);
					}else{
						arMatrix = get_ar_matrix(0,p);
					}
													
					double [][] arTerm = MatrixOperations.multiplication(arMatrix, laggedObs);
					meanVec = MatrixOperations.add(meanVec, arTerm);
				}
					
				NormalDistribution normal = new NormalDistribution(meanVec, SigmaMatrix);
				obs_model_probs[m][t] = normal.get_multivariateNormalPDF(curObs);
				
			}
				
		}
		
	}
	
	
	public static void check_convergence_4_MS_VAR(double curLogLik, double prevLogLik){
		
		double convMetric = Math.abs((curLogLik-prevLogLik)/prevLogLik);
		
		if(convMetric <= convergence_criterion){
			convergence_reached = true;
		}
		
	}
	
	
	public static void expectation_maximization_4_MSI_VAR(){
		
		MSI_VAR ms_obj = new MSI_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
		
		//MatrixOperations.print_matrix(get_sigma_matrix(0));
		
		calc_obs_model_probs_4_MS_VAR();	
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();
					
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			//System.out.println("");
			//MatrixOperations.print_matrix(get_sigma_matrix(0));
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}
	
	
	public static void expectation_maximization_4_MSH_VAR(){
		
		MSH_VAR ms_obj = new MSH_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
			
		calc_obs_model_probs_4_MS_VAR();
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();
					
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}
	
	
	public static void expectation_maximization_4_MSA_VAR(){
		
		MSA_VAR ms_obj = new MSA_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
		
		calc_obs_model_probs_4_MS_VAR();
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();
								
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}
	
	
	public static void expectation_maximization_4_MSAH_VAR(){
		
		MSAH_VAR ms_obj = new MSAH_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
		
		calc_obs_model_probs_4_MS_VAR();
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();			
			
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}
	
	
	public static void expectation_maximization_4_MSIA_VAR(){
		
		MSIA_VAR ms_obj = new MSIA_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
			
		calc_obs_model_probs_4_MS_VAR();
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();
				
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}
	
	
	public static void expectation_maximization_4_MSIAH_VAR(){
		
		MSIAH_VAR ms_obj = new MSIAH_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
				
		calc_obs_model_probs_4_MS_VAR();
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();
							
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}

	
	public static void expectation_maximization_4_MSIH_VAR(){
		
		MSIH_VAR ms_obj = new MSIH_VAR();
		ms_obj.pre_calulation();
		ms_obj.initialize();
				
		calc_obs_model_probs_4_MS_VAR();
		
		double curLogLik  = Double.MIN_VALUE;
		double prevLogLik = Double.MIN_VALUE;	
		
		for(int i=0; i<max_iterations; i++){
						
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
	
			ms_obj.calc_and_set_parameters();
									
			update_transition_matrix_4_MS_VAR();
												
			calc_obs_model_probs_4_MS_VAR();
			
			curLogLik = calc_log_likelihood_4_MS_VAR();
			
			check_convergence_4_MS_VAR(curLogLik, prevLogLik);
			
			if(convergence_reached == true){
				System.out.println("EM algorithm for MS-VAR has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;
				break;
			}
			
			prevLogLik = curLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = curLogLik;
		
	}
	

	@SuppressWarnings("static-access")
	public static void plotInputSeries(){
   		
        GenGraphics obj_graph = new GenGraphics();
		
        int defaultCols = 3;
        int nCols = 0;
        int nRows = 1;       
        int counter = 0;
        
        for(int i=0; i<n_variables; i++){
        	if(i<=(defaultCols-1)){
        		nCols += 1;
        	}
        	if(counter == defaultCols){
        		nRows += 1;
        		counter = 0;
        	}
        	counter++;
        }
        
        obj_graph.setNumberOfPlotColums(nCols);
        obj_graph.setNumberOfPlotRows(nRows);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

	 	double [][] xAxis = new double [n_observations][1];
	 	
	 	for(int i=0; i<n_observations; i++){
	 		xAxis[i][0] = i +1;
	 	}
	 	
	 	String []   title = new String [n_observations];
	 	String []   yLabel = new String [n_observations];
	 	
	 	for(int i=0; i<n_variables; i++){	 		 		
	 		double [][] data = MatrixOperations.get_column_vec_from_matrix(observed_variables, i);
	 		obj_graph.plotLines(xAxis,data,true);	 		
	 		title[i]  = inputData.selected_colnames[i];
	 		yLabel[i] = title[i];
	 	}
	 		 	
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
   		
   	}
	
	
	//returns List of (i) const vector, (ii) (vectorized) p-dimensional AR matrices and (iii) (vectorized) Sigma matrix 
	@SuppressWarnings("static-access")
	public ArrayList<ArrayList<List<Double>>> est_conventional_VAR(){
		
		ArrayList<ArrayList<List<Double>>> est_res = new ArrayList<ArrayList<List<Double>>>(3);
		ArrayList<List<Double>> listOfARMatrices   = new ArrayList<List<Double>>(lag);
		ArrayList<List<Double>> listOfConstVec     = new ArrayList<List<Double>>(1);
		ArrayList<List<Double>> listOfSigmaMatrix  = new ArrayList<List<Double>>(1);
		
    	double [][] X = get_lagged_Y();    	
    	
    	int nARpars = n_variables*lag;
    	
		double [][] constant = new double [n_variables][1];				
    	double [][] arMatrices = new double [n_variables][nARpars];
    	double [][] residuals  = new double [n_usedObservations][n_variables];
    	
    	for(int i=0; i<n_variables; i++){
    		
    		double [][] y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, i, i);
      	        	
        	LinearRegression obj_lm = new LinearRegression(y, X, true);
        	
        	obj_lm.do_parameter_estimation();
    		
        	double [][] lmRes = obj_lm.get_residuals();
        	
        	for(int j=0; j<n_usedObservations; j++){
        		residuals[j][i]=lmRes[j][0];
        	}
        		
    		constant[i][0] = obj_lm.get_est_constant();
        	double [][] pars = obj_lm.get_est_parameters();
        	
        	for(int j=0; j<nARpars; j++){
        		//A_1,...,A_p
        		arMatrices[i][j] = pars[j][0];
        	}
        	
    	}
    	
    	double [][] SigmaMatrix = GeneralMath.cov(residuals);
  	
    	//Const vector
    	listOfConstVec.add(MatrixOperations.vecAsList(constant));
    	
    	//AR matrices A_1,...,A_p
    	int startColIdx = 0;
    	int endColIdx   = n_variables-1;

		for(int p=0; p<lag; p++){
			
			double [][] ar_matrix = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(arMatrices, 0, (n_variables-1), startColIdx, endColIdx);
			
			List<Double> arList = MatrixOperations.vecAsList(ar_matrix);
			listOfARMatrices.add(arList);
			startColIdx = endColIdx+1;
			endColIdx   = startColIdx+n_variables-1;
		}
		
		//Sigma matrix
    	listOfSigmaMatrix.add(MatrixOperations.vecAsList(SigmaMatrix));
    	
    	est_res.add(0,listOfConstVec);
    	est_res.add(1,listOfARMatrices);
    	est_res.add(2,listOfSigmaMatrix);
    	
    	return est_res;
    	
	}
	
	
	public static void set_default_values_MS_VAR(){
		
		max_iterations = 1000;
		convergence_criterion = 1e-06;
		
	}
	
	
	public static void calc_MS_VAR(){
		
		if(ms_var_type != null){
			
			if(ms_var_type == "MSI"){
				expectation_maximization_4_MSI_VAR();
			}
			
			if(ms_var_type == "MSH"){
				expectation_maximization_4_MSH_VAR();
			}
			
			if(ms_var_type == "MSA"){
				expectation_maximization_4_MSA_VAR();
			}
			
			if(ms_var_type == "MSAH"){
				expectation_maximization_4_MSAH_VAR();
			}
			
			if(ms_var_type == "MSIA"){
				expectation_maximization_4_MSIA_VAR();
			}
			
			if(ms_var_type == "MSIAH"){
				expectation_maximization_4_MSIAH_VAR();
			}
			
			if(ms_var_type == "MSIH"){
				expectation_maximization_4_MSIH_VAR();
			}
			
		}else{
			System.out.println("No MS-VAR type supplied for calculation.");
		}
		
	}
	
	
	public static double [][] get_MS_VAR_fitted_values(){
		
		double [][] fitted_values = new double [n_usedObservations][n_variables];
		
		boolean mySwitch = my_switching();
		boolean arSwitch = ar_switching();
			
		double [][] probs    = get_filtered_probs();
		
		int idx = startIdx;
		
		for(int t=0; t<n_usedObservations; t++){
			
			double [][] constant = new double [n_variables][1];
			double [][] arTerm   = new double [n_variables][1];
			
			if(mySwitch == true){
				for(int i=0; i<n_states; i++){
					double prob = probs[t][i];
					constant = MatrixOperations.add(MatrixOperations.scalar_multiplication(prob,get_my(i)),constant);
				}
			}else{
				constant = get_my(0);
			}
			
			if(arSwitch == true){
				for(int i=0; i<lag; i++){	
					double [][] laggedVar = MatrixOperations.get_row_vec_from_matrix(observed_variables, idx-(i+1));
					for(int j=0; j<n_states; j++){	
						double prob = probs[t][j];
						double [][] arMatrix = get_ar_matrix(j,i);						
						arTerm = MatrixOperations.add(MatrixOperations.scalar_multiplication(prob,MatrixOperations.multiplication(arMatrix, laggedVar)),arTerm);
					}		
				}
			}else{
				for(int i=0; i<lag; i++){
					double [][] laggedVar = MatrixOperations.get_row_vec_from_matrix(observed_variables, idx-(i+1));
					double [][] arMatrix = get_ar_matrix(0,i);						
					arTerm = MatrixOperations.add(MatrixOperations.multiplication(arMatrix, laggedVar),arTerm);
				}
			}

			for(int i=0; i<n_variables; i++){
				fitted_values[t][i]= constant[i][0] + arTerm[i][0];
			}
				
			idx++;
		}
		
		return fitted_values;
		
	}
	
	
	public static void print_MS_VAR_parameters(){
		
		boolean mySwitch = my_switching();
		boolean arSwitch = ar_switching();
		boolean sigmaSwitch = sigma_switching();
			
		if(mySwitch == true){
			for(int i=0; i<n_states; i++){
				System.out.println("");
				System.out.println("State "+(i+1)+ ": Constant");
				MatrixOperations.print_matrix(get_my(i));
			}
		}else{
			System.out.println("");
			System.out.println("Constant");
			MatrixOperations.print_matrix(get_my(0));
		}
		
		if(arSwitch == true){
			for(int i=0; i<n_states; i++){				
				for(int j=0; j<lag; j++){
					System.out.println("");
					System.out.println("State "+(i+1)+ ": AR["+(j+1)+"] matrix");
					MatrixOperations.print_matrix(get_ar_matrix(i,j));
				}		
			}
		}else{
			for(int j=0; j<lag; j++){
				System.out.println("");
				System.out.println("AR["+(j+1)+"] matrix");
				MatrixOperations.print_matrix(get_ar_matrix(0,j));
			}
		}
		
		if(sigmaSwitch == true){
			for(int i=0; i<n_states; i++){
				System.out.println("");
				System.out.println("State "+(i+1)+ ": Sigma matrix");
				MatrixOperations.print_matrix(get_sigma_matrix(i));
			}
		}else{
			System.out.println("");
			System.out.println("Sigma matrix");
			MatrixOperations.print_matrix(get_sigma_matrix(0));
		}
		
		System.out.println("");
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotInputSeriesWithErrorBars(){
   		
        GenGraphics obj_graph = new GenGraphics();
		
        int defaultCols = 3;
        int nCols = 0;
        int nRows = 1;       
        int counter = 0;
        
        for(int i=0; i<n_variables; i++){
        	if(i<=(defaultCols-1)){
        		nCols += 1;
        	}
        	if(counter == defaultCols){
        		nRows += 1;
        		counter = 0;
        	}
        	counter++;
        }
        
        obj_graph.setNumberOfPlotColums(nCols);
        obj_graph.setNumberOfPlotRows(nRows);	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);

        double [][] usedData   = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, (n_variables-1));
        double [][] fittedData = get_MS_VAR_fitted_values();
        
	 	double [][] xAxis = new double [n_usedObservations][1];
	 	
	 	for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
	 	
	 	String [] title = new String [n_variables];
	 	String [] yLabel = new String [n_variables];
	 	
	 	for(int i=0; i<n_variables; i++){	 		 		

	 		double [][] obsVar = MatrixOperations.get_column_vec_from_matrix(usedData, i);
	 		double [][] fittedVar = MatrixOperations.get_column_vec_from_matrix(fittedData, i);
	 		
	 		obj_graph.plotLines(xAxis,fittedVar,true);
	 		obj_graph.plotPoints(xAxis,obsVar,false);	    
		    obj_graph.drawErrorBars(true);
	 		 		
	 		title[i]  = inputData.selected_colnames[i];
	 		yLabel[i] = title[i];
	 		
	 	}
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
   		
   	}
	
	
    //test client
    @SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"France"};//{"France","USA","UK","Germany"};
    	int numberOfLags   = 1;
    	int numberOfStates = 2;
    	String msType = "MSI";
    	
    	read_input_data(file, true, true);
	    
    	int nData = inputData.numberOfRows-1;
    	String [] rownames = new String [nData];
    	for(int i=0; i<nData;i++){
    		rownames[i] = Integer.toString(i+1);
    	}
    	
    	select_input_data(rownames, colnames);
    	
    	//Plot input series:
    	//plotInputSeries();
    	
    	MS_VAR obj_test = new MS_VAR(numberOfLags, (nData-numberOfLags), numberOfLags, numberOfStates, msType);
    	
    	obj_test.calc_MS_VAR();
    	
    	
    	plotInputSeriesWithErrorBars();
    	
    	//Plot results
		double [][] f_probs = get_smoothed_probs();
		f_probs = MatrixOperations.get_column_vec_from_matrix(f_probs, 0);
		double [][] x = new double [f_probs.length][1];
		double [][] obs = new double [f_probs.length][1];
		String [] ylabel = {"1"};
		
		for(int j=0; j<f_probs.length; j++){
			x[j][0]=j;
			obs[j][0]=observed_variables[startIdx+j][0];
		}
		
		GenGraphics gr = new GenGraphics();
		gr.setNumberOfPlotColums(2);
		gr.setNumberOfPlotRows(1);		 	
		gr.setGraphWidth(1000);
		gr.setGraphHeight(500);
		gr.plotLines(x,f_probs,true,Color.RED);
		gr.plotLines(x,obs,true,Color.BLACK);
		gr.setYLabel(ylabel, null, "10");
		gr.setNumberOfDigits4XAxis(0);
		gr.plot();
		
		print_MS_VAR_parameters();	
		
		System.out.println(log_likelihood);
    	
    }
	
	
}
