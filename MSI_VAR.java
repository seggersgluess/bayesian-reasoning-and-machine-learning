package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import TimeSeriesAnalysis.VAR;
import Utilities.Utilities;

public class MSI_VAR extends HMM{
	
	double [][] smoothed_probs;
	
	double [][] Y;
	double [][] X_bar;
	double [][] Xi_hat;
	double [][] xi_hat;
	double [][] xi_hat_diag;
	double [][] Z;
	double [][] B;
	double [][] U;
		
	
	//checks validity of supplied input parameter
	public void checkInputParameter(){
		
		if(startIdx > n_observations){
			throw new RuntimeException("Invalid start index supplied.");
		}
		
		if(endIdx > n_observations){
			throw new RuntimeException("Invalid end index supplied.");
		}
		
		if((startIdx-lag) < 0){
			throw new RuntimeException("Invalid lag number supplied.");
		}
			
	}
	
	
	//sets TxK matrix Y=[y_1,...,y_T] of observed variables
	public void set_Y(){
		
		int T = n_usedObservations;
		
		Y = new double [T][n_variables];
		
		for(int t=0; t<T; t++){
			int idx = startIdx+t;
			for(int k=0; k<n_variables; k++){
				Y[t][k] = observed_variables[idx][k];
			}	
		}
			
	}
	
	
	//sets TxKP matrix X_bar = [Y_-1,...,Y_-p]
	public void set_X_bar(){
		
		int T = n_usedObservations;
		
		X_bar = new double [T][lag*n_variables];
		
		for(int t=0; t<T; t++){
			int curIdx = startIdx+t;
			int idx    = 0;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				for(int k=0; k<n_variables; k++){
					X_bar[t][idx] = observed_variables[lagIdx][k];
				    idx++;
				}
			}
		}
			
	}
	
	
	//sets smoothed_probs (bold) Xi_hat
	public void set_smoothed_probs(){
	
		smoothed_probs = get_smoothed_probs();
		
	}
	
	
	//sets MTx1 xi_hat vector and MTxMT diagonal matrix of xi_hat on its diagonal
	public void set_smoothed_probs_vec(){
		
		int nRows = n_states*n_usedObservations;
		
		xi_hat = new double [nRows][1];
		
		int idx = 0;
		
		for(int r=0; r<n_states; r++){
			for(int t=0; t<n_usedObservations; t++){
				xi_hat[idx][0] = smoothed_probs[t][r];
				idx++;
			}
			
		}
		
		xi_hat_diag = MatrixOperations.diagonal(xi_hat);
		
	}
	
	
	//sets MxM matrix Xi_hat
	public void set_diag_matrix_summed_smoothed_probs(){
		
		double [][] e = new double [1][n_usedObservations];
		
		for(int i=0; i<n_usedObservations; i++){
			e[0][i] = 1.0;
		}
		
		Xi_hat = MatrixOperations.transpose(MatrixOperations.multiplication(e, smoothed_probs));
		Xi_hat = MatrixOperations.diagonal(Xi_hat);
		
	}
	
	
	//sets MTx(M+KP) matrix Z
	public void set_Z(double [][] X_bar){
		
		double [][] diag = MatrixOperations.identity(n_states);
		double [][] e_T  = MatrixOperations.unit_vector(n_usedObservations);
		double [][] e_M  = MatrixOperations.unit_vector(n_states);
		
		double [][] firstPart = MatrixOperations.kronecker(diag, e_T);
		double [][] secPart   = MatrixOperations.kronecker(e_M, X_bar);
		
		int nRows = n_states*n_usedObservations;
		int nCols = n_states + n_variables*lag;
		
		Z = new double [nRows][nCols];
		
		for(int r=0; r<nRows; r++){
			
			for(int c=0; c<n_states; c++){
				Z[r][c] = firstPart[r][c];
			}
			
			int idx = 0;
			
			for(int c=n_states; c<nCols; c++){
				Z[r][c] = secPart[r][idx];
				idx++;
			}
			
		}
	
	}
	
	
	//sets (M+Kp)xK matrix B^T
	public void set_B(){
		
		//First term of B
		int nRows = n_states+n_variables*lag;
		int nCols = nRows;
		
		double [][] B_1 = new double [nRows][nCols];
		
		//Upper left
		for(int r=0; r<n_states; r++){
			for(int c=0; c<n_states; c++){
				B_1[r][c] = Xi_hat[r][c];
			}
		}
		
		double [][] smoothdProbsTrans = MatrixOperations.transpose(smoothed_probs);
		double [][] X_bar_Trans       = MatrixOperations.transpose(X_bar);
		
		double [][] X = MatrixOperations.multiplication(smoothdProbsTrans, X_bar);
		
		//Upper right
		for(int r=0; r<n_states; r++){			
		    int idx = 0;		    
			for(int c=n_states; c<nCols; c++){
				B_1[r][c] = X[r][idx];
				idx++;
			}			
		}
		
		//Lower left
		X = MatrixOperations.multiplication(X_bar_Trans, smoothed_probs);
		
		int idx = 0;
		
		for(int r=n_states; r<nRows; r++){				
			for(int c=0; c<n_states; c++){
				B_1[r][c] = X[idx][c];
			}		
			idx++;			
		}
		
		//Upper right
		X = MatrixOperations.multiplication(X_bar_Trans, X_bar);
		
		idx = 0;
		
		for(int r=n_states; r<nRows; r++){	
			int idx2 = 0;
			for(int c=n_states; c<nCols; c++){
				B_1[r][c] = X[idx][idx2];
				idx2++;
			}		
			idx++;			
		}
		
		B_1 = MatrixOperations.inverse(B_1);
		
		//Second term of B
		nCols = n_variables;
		
		double [][] B_2 = new double [nRows][nCols];
		
		X = MatrixOperations.multiplication(smoothdProbsTrans, Y);
		
		for(int r=0; r<n_states; r++){
			for(int c=0; c<n_variables; c++){
				B_2[r][c] = X[r][c];
			}
		}
		
		X = MatrixOperations.multiplication(X_bar_Trans, Y);
		
		idx = 0;
		
		for(int r=n_states; r<nRows; r++){
			for(int c=0; c<n_variables; c++){
				B_2[r][c] = X[idx][c];
			}
			idx++;
		}
		
		B = MatrixOperations.multiplication(B_1, B_2);
		
	}
	
	
	public void set_U(){
		
		double [][] e = MatrixOperations.unit_vector(n_states);
		
		double [][] term1 = MatrixOperations.kronecker(e, Y);
		double [][] term2 = MatrixOperations.multiplication(Z,B);
				
		U = MatrixOperations.substract(term1, term2);
		
	}
	

	//returns (scalar) intercept my(s) w.r.t. selected state from parameter matrix B
	public double [][] get_my_of_MSI_4_state(int state){
			
		//B = [my(1)^T,...,my(M)^T,A_1^T,...,A_p^T]
		double [][] my = MatrixOperations.get_row_vec_from_matrix(B, state);
					
		return my;
			
	}
		
		
	//returns KxK matrix A(s)_p w.r.t. selected lag p
	public double [][] get_ar_matrix_of_MSI(int p){
			
		//B = [my(1),...,my(M),A_1^T,...,A_p^T]			
		double [][] ar_matrix = new double [n_variables][n_variables];
		int startIdx = n_states+n_variables*p;
			
		for(int r=0; r<n_variables; r++){
			for(int c=0; c<n_variables; c++){
				//A_p^T --> transpose here!
				ar_matrix[c][r] = B[startIdx+r][c];
			}
		}

		return ar_matrix;
			
	}
	
	
	public void set_my(){
		
		my = new ArrayList<List<Double>>();
		
		for(int m=0; m<n_states; m++){
			
			double [][] myVec = get_my_of_MSI_4_state(m);
			
			//Transform to List
			List<Double> myList = MatrixOperations.vecAsList(myVec);			
			my.add(myList);
		}
			
	}
	
	
	public void set_ARMatrix(){
		
		ARMatrices = new ArrayList<ArrayList<List<Double>>>();
		
		ArrayList<List<Double>> listOfARMatrices = new ArrayList<List<Double>>(lag); 
		for(int p=0; p<lag; p++){
			double [][] ar_matrix = get_ar_matrix_of_MSI(p);
			List<Double> arList = MatrixOperations.vecAsList(ar_matrix);
			listOfARMatrices.add(arList);
		}
		
		ARMatrices.add(listOfARMatrices);
			
	}
	
	
	public void set_Sigma(){
		
		Sigma = new ArrayList<List<Double>>();
		
		double T = 1.0/GeneralMath.sum(xi_hat);
		
		double [][] SigmaMatrix = new double [n_variables][n_variables];
		
		SigmaMatrix = MatrixOperations.multiplication(MatrixOperations.transpose(U), xi_hat_diag);
		SigmaMatrix = MatrixOperations.multiplication(SigmaMatrix, U);
		SigmaMatrix = MatrixOperations.scalar_multiplication(T, SigmaMatrix);
		
		Sigma.add(MatrixOperations.vecAsList(SigmaMatrix));
		
	}
	
	
	//Initialize with conventional VAR(p)-model
	public void initialize(){
		
		double [] summedSeries = new double [n_usedObservations];
		
		for(int i=0; i<n_usedObservations; i++){
			for(int j=0; j<n_variables; j++){
				summedSeries[i] += observed_variables[startIdx+i][j];
			}
		}
		
		int [] sortedIdxs = Utilities.get_idxs_for_sorted_vec(summedSeries);
		
		int splitLength = Math.round(n_usedObservations/n_states);
		int idx = 0;
		
		for(int m=0; m<n_states; m++){
			
			if(m==(n_states-1)){
				splitLength = n_usedObservations-idx;
			}
			
			double [][] sortedSeries = new double [splitLength][n_variables];
			
			for(int i=0; i<splitLength; i++){
				for(int j=0; j<n_variables; j++){
					sortedSeries[i][j] = observed_variables[(startIdx+sortedIdxs[idx])][j];	
				}
				idx++;
			}
					
			int startIdx4VAR = lag;
			int endIdx4VAR   = sortedSeries.length-1;
			VAR obj_VAR = new VAR(sortedSeries, startIdx4VAR, endIdx4VAR, lag);
			
			obj_VAR.est_conventional_VAR();
			ArrayList<ArrayList<List<Double>>> var_est_res = obj_VAR.get_VAR_est_parameters();
			
	    	my.add(var_est_res.get(0).get(0));
					
		}
				
		VAR obj_VAR = new VAR(observed_variables, startIdx, endIdx, lag);
		
		obj_VAR.est_conventional_VAR();
		ArrayList<ArrayList<List<Double>>> var_est_res = obj_VAR.get_VAR_est_parameters();
		 
		ARMatrices.add(var_est_res.get(1));
		Sigma.add(var_est_res.get(2).get(0));
			
		initialize_trans_probs();
		
	}
	
	
	//--- do all routines for EM algorithm ---
	
	//pre calculation of input
	public void pre_calulation(){
			
		set_Y();
		set_X_bar();
		set_Z(X_bar);
		
	}
	
	
	//parameter calculation and setting
	public void calc_and_set_parameters(){
		
		//Smoothed Probs as vector
		set_smoothed_probs();		
		//Xi_hat
		set_diag_matrix_summed_smoothed_probs();
		//xi_hat
		set_smoothed_probs_vec();
		
		set_B();
		set_U();
		
		set_my();
		set_ARMatrix();
		set_Sigma();
		
	}
	
}
