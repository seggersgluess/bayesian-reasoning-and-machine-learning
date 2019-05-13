package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.MatrixOperations;

public class MSIAH_VAR extends HMM{

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
		
	//List of M vectorized parameter matrices B_m
	ArrayList<List<Double>> B;
		
    //List of M vectorized covariance matrices Sigma
	ArrayList<List<Double>> SigmaMatrices;;
		
		
	//constructor
	public MSIAH_VAR(int startIdx, int endIdx, int lag, int stateNumber){
			
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
		
	
	//calculates vectorized (Kp+1)xK matrices B_m and fills them into an ArrayList
	public void set_B(){
		
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

	
	//returns TxK matrix of residuals dependent on state m
	public double [][] get_U_m(int state){
		
		if(state>(n_states-1)){
			throw new RuntimeException("No valid state number supplied.");
		}
		
		double [][] B_m = get_B(state);		
		double [][] U = MatrixOperations.substract(Y,MatrixOperations.multiplication(X_bar, B_m));
		
		return U;
		
	}
	
	
	//sets KxK sigma matrices dependent of states m=1,2,...,M
	public void set_SigmaMatrices(){
			
		for(int m=0; m<n_states; m++){
				
			List<Double> sigmaVec = new ArrayList<Double>((int)Math.pow(n_variables, 2.0));				
				
			double [][] xi_hat_m = get_Xi_hat_m(m);
			double T_m = 1.0/MatrixOperations.trace(xi_hat_m);
			double [][] U_m = get_U_m(m);
			
			double [][] sigmaMatrix = MatrixOperations.multiplication(MatrixOperations.transpose(U_m), xi_hat_m);
			sigmaMatrix = MatrixOperations.multiplication(sigmaMatrix, U_m);
			sigmaMatrix = MatrixOperations.scalar_multiplication(T_m, sigmaMatrix);
				
			for(int i=0; i<n_variables; i++){				
				for(int j=0; j<n_variables; j++){				
					sigmaVec.add(sigmaMatrix[j][i]);								
				}			
			}
				
			SigmaMatrices.add(sigmaVec);
				
		}
				 	
	}
	
	
	//returns sigma matrix for state m from list of vectorized sigma matrices
	public double [][] get_sigma_matrix_4_state(int state){
			
		if(state > (n_states-1)){
			throw new RuntimeException("No valid state number supplied.");
		}
			
		double [][] sigmaMatrix = new double [n_variables][n_variables];
			
		int idx = 0;
			
		for(int i=0; i<n_variables; i++){				
			for(int j=0; j<n_variables; j++){				
				sigmaMatrix[j][i] = SigmaMatrices.get(state).get(idx);
				idx++;
			}			
		}
			
		return sigmaMatrix;
			
	}
	
}
