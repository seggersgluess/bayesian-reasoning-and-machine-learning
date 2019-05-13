package HiddenMarkovModels;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class MSA_VAR extends HMM{

	//T=endIdx-startIdx+1
	int n_usedObservations; 
	//t=1
	int startIdx;
	//t=T
	int endIdx;
	//p
	int lag;
	
	double [][] smoothed_probs;
	
	double [][] X_bar;
	double [][] y;
	double [][] Y;
	
	//Matrix including all state dependent AR matrices
	double [][] B;
	
	double [][] U;
	double [][] SigmaMatrix;
	
	
	//constructor
	public MSA_VAR(int startIdx, int endIdx, int lag, int stateNumber){
		
		this.startIdx = startIdx;
		this.endIdx   = endIdx;
		this.lag      = lag;
		n_states      = stateNumber;
		
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
	
	
	//sets TxKP matrix X_bar = [Y_-1,...,Y_-p]
	public void set_X_bar(){
		
		int T = n_usedObservations;
		
		X_bar = new double [T][lag*n_variables];
		
		for(int t=0; t<T; t++){
			int curIdx = startIdx+t;
			int idx     = 0;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				for(int k=0; k<n_variables; k++){
					X_bar[t][idx] = observed_variables[lagIdx][k];
				    idx++;
				}
			}
		}
			
	}
	
	
	//returns diagonal matrix with smoothed probs of state m on its diagonal
	public double [][] get_Xi_hat_m(int state){
		
		if(state > (n_states-1)){
			throw new RuntimeException("No valid state number supplied.");
		}
			
		double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, state);					
		double [][] Xi_hat_m = MatrixOperations.diagonal(smoothedProbs4State);
		
		return Xi_hat_m;
		
	}
	
	
	//returns Tx(1+MKp) matrix [e_T,0,...,X_bar,...,0] with ones in first column and X_bar at m-th position
	public double [][] get_X_bar_m(int state){
		
		double [][] X_bar_m = new double [n_usedObservations][(1+n_states*n_variables*lag)];
		
		int startIdx = 1+n_states*n_variables*state;
		
		for(int m=0; m<n_states; m++){
			
			for(int i=0; i<n_usedObservations; i++){
				for(int j=0; j<(n_variables*lag); j++){
					X_bar_m[i][startIdx+j] = X_bar[i][j];
				}
			}
			
		}
	
		return X_bar_m;
		
	}
	
	
	//sets (1+MKP)xK matrix with intercept and state dependent AR matrices A_1, A_2,...,A_p
	public void set_B(){
		
		double [][] I_k = MatrixOperations.identity(n_variables);
		
		double [][] term1 = new double [(1+n_states*n_variables*lag)][(1+n_states*n_variables*lag)];
		double [][] term2 = new double [(1+n_states*n_variables*lag)][n_usedObservations];
		
		for(int m=0; m<n_states; m++){
			
			double [][] X_bar = get_X_bar_m(m);
			double [][] X_bar_trans = MatrixOperations.transpose(X_bar);
			double [][] Xi_hat_m = get_Xi_hat_m(m);
			term1 = MatrixOperations.add(term1, MatrixOperations.multiplication(MatrixOperations.multiplication(X_bar_trans, Xi_hat_m),X_bar));
			
			term2 = MatrixOperations.add(term2, MatrixOperations.multiplication(X_bar_trans, Xi_hat_m));
			
		}
		
		term1 = MatrixOperations.inverse(term1);
		B = MatrixOperations.kronecker(MatrixOperations.multiplication(term1, term2), I_k);
		B = MatrixOperations.multiplication(B, y);
		B = MatrixOperations.get_matrix_from_vec(B, (1+n_states*n_variables*lag), n_variables);
		
	}
	
	
	//calculates and sets KxK covariance matrix
	public void set_SigmaMatrix(){
		
		SigmaMatrix = new double [n_variables][n_variables];
		
		double T = 0.0;
		
		for(int m=0; m<n_states; m++){
			
			double [][] X_bar_m = get_X_bar_m(m);
			double [][] U_m = MatrixOperations.substract(Y, MatrixOperations.multiplication(X_bar_m, B));
			double [][] U_m_trans = MatrixOperations.transpose(U_m);
			
			SigmaMatrix = MatrixOperations.add(SigmaMatrix, MatrixOperations.multiplication(MatrixOperations.multiplication(U_m_trans, X_bar_m),U_m));
			
			T = T+ GeneralMath.sum(MatrixOperations.get_column_from_matrix(smoothed_probs, m));
			
		}
		
		T=1.0/T;
		
		SigmaMatrix = MatrixOperations.scalar_multiplication(T, SigmaMatrix);
		
	}
	

	//returns KxK covariance matrix
	public double [][] get_SigmaMatrix(){
		
		return SigmaMatrix;
		
	}
	
}
