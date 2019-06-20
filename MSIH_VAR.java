package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import TimeSeriesAnalysis.VAR;
import Utilities.Utilities;

public class MSIH_VAR extends HMM{
		
	double [][] smoothed_probs;		
	double [][] y;
	double [][] Y;
	double [][] B;
	
	
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
	
	
	//sets TKx1 vector y=[y_1^T,...,y_T^T] of observed variables
	public void set_y(){
		
		int T = n_usedObservations;
		
		y = new double [T*n_variables][1];
		
		int idx1 = 0;		
		for(int t=0; t<T; t++){
			int idx2 = startIdx+t;
			for(int k=0; k<n_variables; k++){
				y[idx1][0] = observed_variables[idx2][k];
				idx1++;
			}	
		}
			
	}
	
	
	//sets smoothed_probs (bold) Xi_hat
	public void set_smoothed_probs(){
		
		smoothed_probs = get_smoothed_probs();
			
	}
	
	
	//sets TxKP matrix X_bar = [e_T x jotta_m,Y_-1,...,Y_-p]
	public double [][] get_X_bar_m(int state){
		
		int T = n_usedObservations;
		
		double [][] X_bar = new double [T][n_states+lag*n_variables];
		
		for(int t=0; t<T; t++){
			int curIdx = startIdx+t;
			X_bar[t][state] = 1.0;
			int idx = n_states;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				for(int k=0; k<n_variables; k++){
					X_bar[t][idx] = observed_variables[lagIdx][k];
				    idx++;
				}
			}
		}
		
		return X_bar;
		
	}
	
	
	public double [][] get_Xi_hat_m(int state){
		
		if(state > (n_states-1)){
			throw new RuntimeException("No valid state number supplied.");
		}
			
		double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, state);					
		double [][] Xi_hat_m = MatrixOperations.diagonal(smoothedProbs4State);
		
		return Xi_hat_m;
		
	}
	
	
	//sets (M+Kp)xK matrix B
	public void set_B(){
				
		double [][] sumTerm1    = new double [n_variables*(n_states+n_variables*lag)][n_variables*(n_states+n_variables*lag)];
		double [][] sumTerm2    = new double [n_variables*(n_states+n_variables*lag)][n_variables*n_usedObservations];
		
		for(int m=0; m<n_states; m++){
			
			double [][] X_bar_m = get_X_bar_m(m);
			double [][] X_bar_m_trans = MatrixOperations.transpose(X_bar_m);
			double [][] xi_hat_m    = get_Xi_hat_m(m);
						
			double [][] sigmaMatrix = MatrixOperations.inverse(get_sigma_matrix(m));
			
			sumTerm1 = MatrixOperations.add(sumTerm1,MatrixOperations.kronecker(MatrixOperations.multiplication(MatrixOperations.multiplication(X_bar_m_trans, xi_hat_m), X_bar_m),sigmaMatrix));
			sumTerm2 = MatrixOperations.add(sumTerm2, MatrixOperations.kronecker(MatrixOperations.multiplication(X_bar_m_trans,xi_hat_m),sigmaMatrix));
			
		}
		
		sumTerm1 = MatrixOperations.inverse(sumTerm1);
		
		B = MatrixOperations.multiplication(sumTerm1, sumTerm2);
		B = MatrixOperations.multiplication(B,y);		
		B = MatrixOperations.get_matrix_from_vec(B, n_variables, (n_states+n_variables*lag));
		//Set B^T for further calculations!
		B = MatrixOperations.transpose(B);
		
	}
	
	
	//set TxK matrix U of residuals
	public double [][] get_U_m(int state){
		
		double [][] X_bar_m = get_X_bar_m(state);
		//U = Y - X_bar_m * B^T
		double [][] U_m = MatrixOperations.substract(Y,MatrixOperations.multiplication(X_bar_m, B));

		return U_m;
		
	}
		
	
	//returns Kx1 intercept vector my(s) w.r.t. selected state from parameter matrix B
	public double [][] get_my_of_MSIH_4_state(int state){
			
		//B^T = [my(1)^T,...,my(M)^T,A_1^T,...,A_p^T]
		double [][] my = MatrixOperations.get_row_vec_from_matrix(B, state);
					
		return my;
			
	}
	
	
	//returns KxK matrix A(s)_p w.r.t.lag p from parameter matrix B
	public double [][] get_ar_matrix_of_MSIH(int p){
			
		//B^T = [my(1)^T,...,my(M)^T,A_1^T,...,A_p^T]			
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
			
			double [][] myVec = get_my_of_MSIH_4_state(m);
			
			//Transform to List
			List<Double> myList = MatrixOperations.vecAsList(myVec);			
			my.add(myList);
		}
			
	}
	
	
	public void set_ARMatrix(){
		
		ARMatrices = new ArrayList<ArrayList<List<Double>>>();
		
		ArrayList<List<Double>> listOfARMatrices = new ArrayList<List<Double>>(lag); 
		for(int p=0; p<lag; p++){
			double [][] ar_matrix = get_ar_matrix_of_MSIH(p);
			List<Double> arList = MatrixOperations.vecAsList(ar_matrix);
			listOfARMatrices.add(arList);
		}
		
		ARMatrices.add(listOfARMatrices);
			
	}
	
	
	//sets KxK sigma matrices dependent of states m=1,2,...,M
	public void set_Sigma(){
			
		Sigma = new ArrayList<List<Double>>();
		
		for(int m=0; m<n_states; m++){
				
			List<Double> sigmaVec = new ArrayList<Double>((int)Math.pow(n_variables, 2.0));
			double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, m);					
			double T_m = 1.0/GeneralMath.sum(smoothedProbs4State);
				
			double [][] xi_hat_m = get_Xi_hat_m(m);
			double [][] U_m = get_U_m(m);
			
			double [][] sigmaMatrix = MatrixOperations.multiplication(MatrixOperations.transpose(U_m), xi_hat_m);
			sigmaMatrix = MatrixOperations.multiplication(sigmaMatrix, U_m);
			sigmaMatrix = MatrixOperations.scalar_multiplication(T_m, sigmaMatrix);
				
			for(int i=0; i<n_variables; i++){				
				for(int j=0; j<n_variables; j++){				
					sigmaVec.add(sigmaMatrix[j][i]);								
				}			
			}
				
			Sigma.add(sigmaVec);
				
		}
				 	
	}
		
	
	//Initialize with conventional VAR(p)-model
	public void initialize(){
			
		//Estimate conventional VAR(p)
		VAR obj_VAR = new VAR(observed_variables, startIdx, endIdx, lag);
		
		obj_VAR.est_conventional_VAR();
		ArrayList<ArrayList<List<Double>>> var_est_res = obj_VAR.get_VAR_est_parameters();
				
	    //--- fill parameter list ---    		    
	    ARMatrices.add(var_est_res.get(1));
	    	    
	    double [][] varResiduals = obj_VAR.get_VAR_residuals();
	    double [] summedAbsResidualFluct = new double [n_usedObservations];
	    
	    for(int i=0; i<n_usedObservations; i++){
	    	for(int j=0; j<n_variables; j++){
	    		summedAbsResidualFluct[i] += Math.abs(varResiduals[i][j]);
	    	}
	    }
	    
		int [] sortedIdxs = Utilities.get_idxs_for_sorted_vec(summedAbsResidualFluct);
		
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
			obj_VAR = new VAR(sortedSeries, startIdx4VAR, endIdx4VAR, lag);
			
			obj_VAR.est_conventional_VAR();
			var_est_res = obj_VAR.get_VAR_est_parameters();
			
			my.add(var_est_res.get(0).get(0));
			Sigma.add(var_est_res.get(2).get(0));				
		}
	
	    initialize_trans_probs();
			
	}

	//--- do all routines for EM algorithm ---
	
	//pre calculation of input
	public void pre_calulation(){
			
		set_Y();
		set_y();
		
	}
	
	
	//parameter calculation and setting
	public void calc_and_set_parameters(){
		
		set_smoothed_probs();
		set_B();
		set_my();
		set_ARMatrix();
		set_Sigma();
		
	}
	
}
