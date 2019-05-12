package HiddenMarkovModels;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class MSI_VAR extends HMM{

	//T=endIdx-startIdx+1
	int n_usedObservations; 
	//t=1
	int startIdx;
	//t=T
	int endIdx;
	//p
	int lag;
	
	double [][] smoothed_probs;
	
	double [][] Y;
	double [][] X_bar;
	double [][] Xi_hat;
	double [][] xi_hat;
	double [][] Z;
	double [][] B;
	double [][] U;
	double [][] SigmaMatrix;
	
	
	//constructor
	public MSI_VAR(int startIdx, int endIdx, int lag, int stateNumber){
		
		this.startIdx = startIdx;
		this.endIdx   = endIdx;
		this.lag      = lag;
		n_states      = stateNumber;
		
	}
	
	
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
	
	
	//sets MTx1 xi_hat vector
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
		
	}
	
	
	//sets MxM matrix Xi_hat
	public void set_diag_matrix_summed_smoothed_probs(double [][] smoothed_probs){
		
		double [][] e = new double [1][n_usedObservations];
		
		for(int i=0; i<n_usedObservations; i++){
			e[0][i] = 1.0;
		}
		
		Xi_hat = MatrixOperations.multiplication(e, smoothed_probs);
				
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
			}
			
		}
		
	}
	
	
	//sets (M+Kp)xK matrix B
	public void set_B(){
		
		//First term of B
		int nRows = 2*n_states+n_variables*lag;
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
		double [][] term2 = MatrixOperations.multiplication(Z,MatrixOperations.transpose(B));
		
		U = MatrixOperations.substract(term1, term2);
		
	}
	
	
	public void set_SigmaMatrix(){
		
		double T = GeneralMath.sum(xi_hat);
		
		SigmaMatrix = MatrixOperations.multiplication(MatrixOperations.transpose(U), smoothed_probs);
		SigmaMatrix = MatrixOperations.multiplication(SigmaMatrix, U);
		SigmaMatrix = MatrixOperations.scalar_multiplication(T, SigmaMatrix);
		
	}
	
}
