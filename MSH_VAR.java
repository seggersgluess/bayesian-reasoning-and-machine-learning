package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class MSH_VAR extends HMM{

	//T=endIdx-startIdx+1
	int n_usedObservations; 
	//t=1
	int startIdx;
	//t=T
	int endIdx;
	//p
	int lag;
		
	double [][] smoothed_probs;
		
	double [][] y;
	double [][] X_bar;
	double [][] B;
	double [][] U;
	ArrayList<List<Double>> SigmaMatrices;
	
	
	//constructor
	public MSH_VAR(int startIdx, int endIdx, int lag, int stateNumber){
		
		this.startIdx = startIdx;
		this.endIdx   = endIdx;
		this.lag      = lag;
		n_states      = stateNumber;
		
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
	
	
	public double [][] get_Xi_hat_m(int state){
		
		if(state > (n_states-1)){
			throw new RuntimeException("No valid state number supplied.");
		}
			
		double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, state);					
		double [][] Xi_hat_m = MatrixOperations.diagonal(smoothedProbs4State);
		
		return Xi_hat_m;
		
	}
	
	
	//sets smoothed_probs (bold) Xi_hat
	public void set_smoothed_probs(){
	
		smoothed_probs = get_smoothed_probs();
		
	}
	
	
	//sets (Kp)xK matrix B
	public void set_B(){
		
		double [][] X_bar_trans = MatrixOperations.transpose(X_bar);
		double [][] sumTerm1    = new double [(1+n_variables*lag)][(1+n_variables*lag)];
		double [][] sumTerm2    = new double [(1+n_variables*lag)][n_usedObservations];
		
		for(int m=0; m<n_states; m++){
			
			double [][] xi_hat_m    = get_Xi_hat_m(m);
			double [][] sigmaMatrix = MatrixOperations.inverse(get_sigma_matrix_4_state(m));
			
			sumTerm1 = MatrixOperations.add(sumTerm1,MatrixOperations.kronecker(MatrixOperations.multiplication(MatrixOperations.multiplication(X_bar_trans, xi_hat_m), X_bar),sigmaMatrix));
			sumTerm2 = MatrixOperations.add(sumTerm2, MatrixOperations.kronecker(MatrixOperations.multiplication(X_bar_trans,xi_hat_m),sigmaMatrix));
			
		}
		
		sumTerm1 = MatrixOperations.inverse(sumTerm1);
		
		B = MatrixOperations.multiplication(sumTerm1, sumTerm2);
		B = MatrixOperations.multiplication(B,y);
		
	}
	
	
	//set TxK matrix U of residuals
	public void set_U(){
		
		double [][] resVec = MatrixOperations.multiplication(MatrixOperations.kronecker(X_bar, MatrixOperations.identity(n_variables)),B);
		resVec = MatrixOperations.substract(y, U);
		
		U = new double [n_usedObservations][n_variables];
		
		int idx = 0;
		for(int t=0; t<n_usedObservations; t++){			
			for(int k=0; k<n_variables; k++){
				U[t][k] = resVec[idx][0];
				idx++;
			}
			
		}
		
	}
	
	
	//sets KxK sigma matrices dependent of states m=1,2,...,M
	public void set_SigmaMatrices(){
		
		for(int m=0; m<n_states; m++){
			
			List<Double> sigmaVec = new ArrayList<Double>((int)Math.pow(n_variables, 2.0));
			double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, m);					
			double T_m = GeneralMath.sum(smoothedProbs4State);
			
			double [][] xi_hat_m = get_Xi_hat_m(m);
			
			double [][] sigmaMatrix = MatrixOperations.multiplication(MatrixOperations.transpose(U), xi_hat_m);
			sigmaMatrix = MatrixOperations.multiplication(sigmaMatrix, U);
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
