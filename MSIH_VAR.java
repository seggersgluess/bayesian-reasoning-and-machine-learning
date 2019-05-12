package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class MSIH_VAR extends HMM{

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
	double [][] B;
	ArrayList<List<Double>> SigmaMatrices;
	
	
	//constructor
	public MSIH_VAR(int startIdx, int endIdx, int lag, int stateNumber){
		
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
	
	
	//sets (Kp)xK matrix B
	public void set_B(){
				
		double [][] sumTerm1    = new double [(n_states+n_variables*lag)][(n_states+n_variables*lag)];
		double [][] sumTerm2    = new double [(n_states+n_variables*lag)][n_usedObservations];
		
		for(int m=0; m<n_states; m++){
			
			double [][] X_bar_m = get_X_bar_m(m);
			double [][] X_bar_m_trans = MatrixOperations.transpose(X_bar_m);
			double [][] xi_hat_m    = get_Xi_hat_m(m);
			double [][] sigmaMatrix = MatrixOperations.inverse(get_sigma_matrix_4_state(m));
			
			sumTerm1 = MatrixOperations.add(sumTerm1,MatrixOperations.kronecker(MatrixOperations.multiplication(MatrixOperations.multiplication(X_bar_m_trans, xi_hat_m), X_bar_m),sigmaMatrix));
			sumTerm2 = MatrixOperations.add(sumTerm2, MatrixOperations.kronecker(MatrixOperations.multiplication(X_bar_m_trans,xi_hat_m),sigmaMatrix));
			
		}
		
		sumTerm1 = MatrixOperations.inverse(sumTerm1);
		
		B = MatrixOperations.multiplication(sumTerm1, sumTerm2);
		B = MatrixOperations.multiplication(B,y);
		
	}
	
	
	//set TxK matrix U of residuals
	public double [][] get_U_m(int state){
		
		double [][] X_bar_m = get_X_bar_m(state);
		double [][] resVec = MatrixOperations.substract(y,MatrixOperations.multiplication(X_bar_m, B));
		
		double [][] U_m = new double [n_usedObservations][n_variables];
		
		int idx = 0;
		for(int t=0; t<n_usedObservations; t++){			
			for(int k=0; k<n_variables; k++){
				U_m[t][k] = resVec[idx][0];
				idx++;
			}
			
		}
	
		return U_m;
		
	}
	
	
	//sets KxK sigma matrices dependent of states m=1,2,...,M
	public void set_SigmaMatrices(){
			
		for(int m=0; m<n_states; m++){
				
			List<Double> sigmaVec = new ArrayList<Double>((int)Math.pow(n_variables, 2.0));
			double [][] smoothedProbs4State = MatrixOperations.get_column_vec_from_matrix(smoothed_probs, m);					
			double T_m = GeneralMath.sum(smoothedProbs4State);
				
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
