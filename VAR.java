package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Regression.LinearRegression;

public class VAR {

	double [][] observed_variables;
	
	int n_observations;
	int n_usedObservations;
	int n_variables; 
	
	int lag;
	int startIdx;
	int endIdx;
	
	ArrayList<List<Double>> const_vec;
	ArrayList<List<Double>> ar_matrices;
	ArrayList<List<Double>> sigma_matrix;
	
	double [][] residuals;
	
	//constructor
	public VAR(double [][] observed_variables, int startIdx, int endIdx, int lagNumber){
		
		if(startIdx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(endIdx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(startIdx >= endIdx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(lagNumber < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		this.observed_variables = observed_variables;
		this.n_observations = observed_variables.length;
		this.n_variables    = observed_variables[0].length;
		
		this.startIdx = startIdx;
		this.endIdx   = endIdx;
		this.lag      = lagNumber;
	
		n_usedObservations = endIdx-startIdx+1;
		
		check_inputs_for_constructor();

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
	
	
	//calculates List of (i) const vector, (ii) (vectorized) p-dimensional AR matrices and (iii) (vectorized) Sigma matrix 
	@SuppressWarnings("static-access")
	public void est_conventional_VAR(){
			
		ar_matrices   = new ArrayList<List<Double>>(lag);
		const_vec     = new ArrayList<List<Double>>(1);
		sigma_matrix  = new ArrayList<List<Double>>(1);
			
	    double [][] X = get_lagged_Y();    	
	    	
	    int nARpars = n_variables*lag;
	    	
		double [][] constant = new double [n_variables][1];				
	    double [][] arMatrices = new double [n_variables][nARpars];
	    residuals  = new double [n_usedObservations][n_variables];
	    	
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
	    const_vec.add(MatrixOperations.vecAsList(constant));
	    	
	    //AR matrices A_1,...,A_p
	    int startColIdx = 0;
	    int endColIdx   = n_variables-1;

		for(int p=0; p<lag; p++){
				
			double [][] ar_matrix = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(arMatrices, 0, (n_variables-1), startColIdx, endColIdx);
				
			List<Double> arList = MatrixOperations.vecAsList(ar_matrix);
			ar_matrices.add(arList);
		    startColIdx = endColIdx+1;
			endColIdx   = startColIdx+n_variables-1;
		}
			
		//Sigma matrix
	    sigma_matrix.add(MatrixOperations.vecAsList(SigmaMatrix));
	    		    	
	}
	
	
	//returns List of (i) const vector, (ii) (vectorized) p-dimensional AR matrices and (iii) (vectorized) Sigma matrix 
	public ArrayList<ArrayList<List<Double>>> get_VAR_est_parameters(){
		
		ArrayList<ArrayList<List<Double>>> est_res = new ArrayList<ArrayList<List<Double>>>(3);
	    est_res.add(0,const_vec);
	    est_res.add(1,ar_matrices);
	    est_res.add(2,sigma_matrix);
	    
	    return est_res;
		
	}
	
	
	//returns matrix of VAR(p) residuals
	public double [][] get_VAR_residuals(){
		
		return residuals;
		
	}
	
	
	//sets TxKP matrix X_bar = [Y_-1,...,Y_-p]
	public double [][] get_lagged_Y(){
		
		int T = n_usedObservations;
		
		double [][] lagged_Y = new double [T][lag*n_variables];
		
		for(int t=0; t<T; t++){
			int curIdx = startIdx+t;
			int idx     = 0;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				for(int k=0; k<n_variables; k++){
					lagged_Y[t][idx] = observed_variables[lagIdx][k];
				    idx++;
				}
			}
		}
			
		return lagged_Y;
		
	}
	
}
