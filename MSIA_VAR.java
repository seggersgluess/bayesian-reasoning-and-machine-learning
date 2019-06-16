package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.MatrixOperations;

public class MSIA_VAR extends HMM{

	double [][] smoothed_probs;
	
	double [][] Y;
	double [][] X_bar;
	double [][] Xi_hat;
	
	//List of M vectorized parameter matrices B_m
	ArrayList<List<Double>> B;
	
	double [][] U;
	
		
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
	
	
	//sets TxKP matrix X_bar = [1,Y_-1,...,Y_-p]
	public void set_X_bar(){
		
		int T = n_usedObservations;
		
		X_bar = new double [T][1+lag*n_variables];
		
		for(int t=0; t<T; t++){
			int curIdx = startIdx+t;
			X_bar[t][0] = 1.0;
			int idx    = 1;
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
	
	
	//returns Tx1 vector of smoothed probs for state m
	public double [][] get_Xi_hat_m(int state){
		
		if(state > (n_states-1)){
			throw new RuntimeException("No valid state number supplied.");
		}
			
		double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, state);					
		double [][] Xi_hat_m = MatrixOperations.diagonal(smoothedProbs4State);
		
		return Xi_hat_m;
		
	}
	
	
	//calculates vectorized smoothed probs and fills them into a diagonal matrix
	public void set_Xi_hat(){
		
		Xi_hat = new double [n_states*n_usedObservations][1];
		int idx = 0;
		for(int m=0; m<n_states; m++){
			for(int t=0; t<n_usedObservations; t++){
				Xi_hat[idx][0] = smoothed_probs[t][m];
				idx++;
			}
		}
		
		Xi_hat = MatrixOperations.diagonal(Xi_hat);
		
	}
	
	
	//calculates vectorized (Kp+1)xK matrices B_m and fills them into an ArrayList
	public void set_B(){
		
		B = new ArrayList<List<Double>>();
		
		double [][] X_bar_trans = MatrixOperations.transpose(X_bar);
		
		int nRows = n_variables*lag+1;
		int nCols = n_variables;
		
		for(int m=0; m<n_states; m++){
			
			List<Double> B_m_vec = new ArrayList<Double>();
			
			double [][] Xi_hat_m = get_Xi_hat_m(m);
			double [][] B_m = MatrixOperations.multiplication(MatrixOperations.multiplication(X_bar_trans, Xi_hat_m),X_bar);
			B_m = MatrixOperations.inverse(B_m);
			B_m = MatrixOperations.multiplication(B_m, MatrixOperations.multiplication(MatrixOperations.multiplication(X_bar_trans, Xi_hat_m),Y));
			
			for(int i=0; i<nCols; i++){				
				for(int j=0; j<nRows; j++){				
					B_m_vec.add(B_m[j][i]);								
				}			
			}
			
			B.add(B_m_vec);
			
		}
		
	}
	
	
	//returns (Kp+1)xK matrix of parameters dependent on state m
	public double [][] get_B(int state){
		
		if(state>(n_states-1)){
			throw new RuntimeException("Invalid state number supplied.");
		}
		
		double [][] B_m = new double [(n_variables*lag+1)][n_variables];
		int nRows = n_variables*lag+1;
		int nCols = n_variables;
		int idx=0;
		
		for(int i=0; i<nCols; i++){				
			for(int j=0; j<nRows; j++){				
				B_m[j][i] = B.get(state).get(idx);	
				idx++;
			}			
		}
		
		return B_m;
		
	}
	
	
	//sets MTxK matrix of residuals
	public void set_U(){
		
		double [][] B_m_block_matrix = new double [n_states*(n_variables*lag+1)][n_variables];
		
		int nRows = n_variables*lag+1;
		int nCols = n_variables;
		
		for(int m=0; m<n_states; m++){	
			int idx=0;
			int c = nRows*m;
			for(int i=0; i<nCols; i++){				
				for(int j=0; j<nRows; j++){				
					B_m_block_matrix[c+j][i] = B.get(m).get(idx);	
					idx++;
				}			
			}
			
		}
				
		double [][] firstTerm = MatrixOperations.kronecker(MatrixOperations.unit_vector(n_states), Y);
		double [][] secTerm   = MatrixOperations.kronecker(MatrixOperations.identity(n_states), X_bar);
		secTerm = MatrixOperations.multiplication(secTerm, B_m_block_matrix);
		
		U =  MatrixOperations.substract(firstTerm, secTerm);
		
	}
	
	
	//returns Kx1 intercept vector my(s) w.r.t. selected state from parameter matrix B
	public double [][] get_my_of_MSIA_4_state(int state){
		
		double [][] B_m = get_B(state);
		double [][] my = MatrixOperations.get_row_vec_from_matrix(B_m, 0);	
		
		return my;
		
	}
	
	
	//returns KxK matrix A(s)_p w.r.t. selected state and lag p from parameter matrix B
	public double [][] get_ar_matrix_of_MSIA_4_state(int state, int p){
				
		//B_m = [my^T,A_1^T,...,A_p^T]
		double [][] B_m = get_B(state);
		
		double [][] ar_matrix = new double [n_variables][n_variables];
		int startIdx = 1+p*n_variables;
		
		for(int r=0; r<n_variables; r++){
			for(int c=0; c<n_variables; c++){
				//A_p^T --> transpose here!
				ar_matrix[c][r] = B_m[startIdx+r][c];
			}
		}
	
		return ar_matrix;
		
	}
	
	
	public void set_my(){
		
		my = new ArrayList<List<Double>>();
		
		for(int m=0; m<n_states; m++){
			
			double [][] myVec = get_my_of_MSIA_4_state(m);
			
			//Transform to List
			List<Double> myList = MatrixOperations.vecAsList(myVec);			
			my.add(myList);
		}
			
	}

	
	public void set_ARMatrix(){
		
		ARMatrices = new ArrayList<ArrayList<List<Double>>>();
		
		for(int m=0; m<n_states; m++){
			ArrayList<List<Double>> listOfARMatrices = new ArrayList<List<Double>>(lag); 
			for(int p=0; p<lag; p++){
				double [][] ar_matrix = get_ar_matrix_of_MSIA_4_state(m,p);
				List<Double> arList = MatrixOperations.vecAsList(ar_matrix);
				listOfARMatrices.add(arList);
			}
			ARMatrices.add(listOfARMatrices);
		}
			
	}
	
	
	//sets KxK covariance matrix
	public void set_Sigma(){
		
		Sigma = new ArrayList<List<Double>>();
		
		double T = 1.0/MatrixOperations.trace(Xi_hat);
		
		double [][] SigmaMatrix = new double [n_variables][n_variables];
		
		SigmaMatrix = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(U), Xi_hat),U);
		SigmaMatrix = MatrixOperations.scalar_multiplication(T, SigmaMatrix);
		
		Sigma.add(MatrixOperations.vecAsList(SigmaMatrix));
		
	}
	
	
	//Initialize with conventional VAR(p)-model
	public void initialize(){
			
		MS_VAR obj_ms = new MS_VAR(startIdx, endIdx, lag, n_states, ms_var_type);
			
		ArrayList<ArrayList<List<Double>>> var_est_res = obj_ms.est_conventional_VAR();
			
	    //--- fill parameter list ---    		    
	    //state dependent const vectors my(1),...,m(M) and AR matrices A(1),...,A(M)
	    for(int m=0; m<n_states; m++){		 
	    	my.add(var_est_res.get(0).get(0));
	    	ARMatrices.add(var_est_res.get(1));	    	
	    }
	  
	    //Sigma
	    Sigma.add(var_est_res.get(2).get(0));
	    
	    initialize_trans_probs();
			
	}

	
	//--- do all routines for EM algorithm ---
	
	//pre calculation of input
	public void pre_calulation(){
			
		set_Y();
		set_X_bar();
		
	}
	
	
	//parameter calculation and setting
	public void calc_and_set_parameters(){
		
		set_smoothed_probs();
		set_Xi_hat();
		set_B();
		set_U();
		set_my();
		set_ARMatrix();
		set_Sigma();
		
	}
	
}
