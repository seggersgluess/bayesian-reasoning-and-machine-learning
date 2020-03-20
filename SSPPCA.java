package ComponentModels;

import java.util.HashMap;

import Mathematics.MatrixOperations;

public class SSPPCA extends SPPCA{

	double [][] X_1;
	double [][] X_2;
	
	int n_x1 = 0;
	int n_x2 = 0;
	
	double [][] X_1_scaled;
	double [][] X_2_scaled;
	
	double [][] labeled_factors;
	double [][] unlabeled_factors;
	
	double [][] labeled_factor_covariance;
	double [][] unlabeled_factor_covariance;
	

	public SSPPCA(double[][] X_1, double[][] Y, double [][] X_2, int n_factors, boolean center) {
		super(X_1, Y, n_factors, center);		
		//Clear data matrices
		X = null;
		X_scaled = null;
		
		if(X_1[0].length != X_2[0].length) {
			throw new RuntimeException("Inequal number of input features for labeled and unlabeled data X1 and X2 supplied.");
		}
		
		this.X_1 = X_1;
		this.X_2 = X_2;
		
		n_x1 = X_1.length;
		n_x2 = X_2.length;
		
		double [][] X_comb = MatrixOperations.rbind(X_1, X_2);
		my_x = get_sample_means(X_comb);
			
		X_1_scaled = new double [n_x1][n_variables];
		X_2_scaled = new double [n_x2][n_variables];
		
		for(int i=0; i<n_x1; i++) {
			for(int j=0; j<n_variables; j++) {
				X_1_scaled[i][j] = X_1[i][j]-my_x[j][0];
			}		
		}
		for(int i=0; i<n_x2; i++) {
			for(int j=0; j<n_variables; j++) {
				X_2_scaled[i][j] = X_2[i][j]-my_x[j][0];
			}		
		}	
	}

	
	public void em_4_SSPPCA() {
		
		HashMap<String, double [][]> pars = new HashMap<String, double [][]>();
		if(ext_init_values == null) {
			pars = get_init_pars_4_SSPPCA();
		}else {
			pars = ext_init_values;
		}
		
		W_x  = pars.get("W_x");
		W_y  = pars.get("W_y");
		sigma_x = pars.get("sigma_x")[0][0];
		sigma_y = pars.get("sigma_y")[0][0];
			
		double [][] Z1_trans = new double [n_x1][n_factors];
		double [][] Z2_trans = new double [n_x2][n_factors];
		double [][] M = new double [n_factors][n_factors];
		double [][] A_inv   = new double [n_factors][n_factors];
		double [][] diag_L = MatrixOperations.identity(n_factors);
		
		double normX = calc_sq_euclidian_X_4_SSPPCA();
		double normY = calc_sq_euclidian_Y_4_SSPPCA();
		
		for(int j=0; j<iterations; j++) {	
											
			double [][] W_x_trans = MatrixOperations.transpose(W_x); 
			double [][] W_y_trans = MatrixOperations.transpose(W_y);
			double [][] A = MatrixOperations.scalar_multiplication(1.0/sigma_x, MatrixOperations.multiplication(W_x_trans, W_x));
			A = MatrixOperations.add(A,MatrixOperations.scalar_multiplication(1.0/sigma_y, MatrixOperations.multiplication(W_y_trans, W_y)));
			A = MatrixOperations.add(A, diag_L);
				
			A_inv = MatrixOperations.inverse(A);
				
			double [][] W_x_new = new double [n_variables][n_factors];
			double [][] W1_x_new = new double [n_variables][n_factors];
			double [][] W2_x_new = new double [n_variables][n_factors];
			double [][] W_y_new = new double [n_y][n_factors];
			double [][] Cov_z1_i_sum = new double [n_factors][n_factors];
			double [][] Cov_z2_i_sum = new double [n_factors][n_factors];
			
			double [][] expected_z1_i_trans = new double [1][n_factors];
			double [][] expected_z2_i_trans = new double [1][n_factors];
						
			for(int i=0; i<n_x1; i++) {								
				double [][] x_i = MatrixOperations.get_row_vec_from_matrix(X_1_scaled, i);
				double [][] y_i = MatrixOperations.get_row_vec_from_matrix(Y_scaled, i);
				
				double [][] expected_z1_i = calc_expected_z_i_4_SPPCA(W_x_trans, W_y_trans, x_i, y_i);
					
				for(int k=0; k<n_factors; k++) {
					expected_z1_i_trans[0][k] = expected_z1_i[k][0];
					Z1_trans[k][i] = expected_z1_i[k][0];
				}
	
				double [][] expected_z1_i_z1_i_trans = calc_cov_z_i_4_SPPCA(A_inv, expected_z1_i, expected_z1_i_trans);
									 
				W1_x_new = MatrixOperations.add(W1_x_new,MatrixOperations.multiplication(x_i,expected_z1_i_trans));
				Cov_z1_i_sum = MatrixOperations.add(Cov_z1_i_sum, expected_z1_i_z1_i_trans);

				W_y_new = MatrixOperations.add(W_y_new,MatrixOperations.multiplication(y_i,expected_z1_i_trans));				
			}
			
			double [][] Psi = new double [n_factors][n_factors];
			for(int i=0; i<n_factors; i++) {
				Psi[i][i] = sigma_x;
			}
			
		    M = MatrixOperations.multiplication(W_x_trans, W_x);
			M = MatrixOperations.add(M, Psi);
			M = MatrixOperations.inverse(M);
			
			for(int i=0; i<n_x2; i++) {								
				double [][] x_i = MatrixOperations.get_row_vec_from_matrix(X_2_scaled, i);
				
				double [][] expected_z2_i = MatrixOperations.multiplication(MatrixOperations.multiplication(M, W_x_trans), x_i);
					
				for(int k=0; k<n_factors; k++) {
					expected_z2_i_trans[0][k] = expected_z2_i[k][0];
					Z2_trans[k][i] = expected_z2_i[k][0];
				}
	
				double [][] expected_z2_i_z2_i_trans = MatrixOperations.add(M, MatrixOperations.multiplication(expected_z2_i, expected_z2_i_trans));
									 
				W2_x_new = MatrixOperations.add(W2_x_new, MatrixOperations.multiplication(x_i,expected_z2_i_trans));
				Cov_z2_i_sum = MatrixOperations.add(Cov_z2_i_sum, expected_z2_i_z2_i_trans);			
			}
				
			double [][] Cov_z_i_sum = MatrixOperations.add(Cov_z1_i_sum, Cov_z2_i_sum);
			double [][] Cov_z_i_sum_inv = MatrixOperations.inverse(Cov_z_i_sum);
			
			double [][] W_x_sum = MatrixOperations.add(W1_x_new, W2_x_new);
			W_x_new = MatrixOperations.add(W_x_sum, Cov_z_i_sum_inv);
			double [][] W_x_new_trans = MatrixOperations.transpose(W_x_new);
			W_y_new = MatrixOperations.add(W_y_new, MatrixOperations.inverse(Cov_z1_i_sum));
			double [][] W_y_new_trans = MatrixOperations.transpose(W_y_new);	
			
			double sigma_x_new = 0.0;
			double sigma_y_new = 0.0;
			
			sigma_x_new += normX;
			sigma_y_new += normY;

			sigma_x_new += MatrixOperations.trace(MatrixOperations.multiplication(MatrixOperations.multiplication(W_x_new_trans, W_x_new),Cov_z_i_sum));
			sigma_x_new -= 2.0*MatrixOperations.trace(MatrixOperations.multiplication(W_x_new, MatrixOperations.transpose(W_x_sum)));
			sigma_x_new /= (n_variables*(n_x1+n_x2));
			
			sigma_y_new += MatrixOperations.trace(MatrixOperations.multiplication(MatrixOperations.multiplication(W_y_new_trans, W_y_new),Cov_z1_i_sum));
			sigma_y_new -= 2.0*MatrixOperations.trace(MatrixOperations.multiplication(Y_scaled, MatrixOperations.multiplication(W_y_new, Z1_trans)));
			sigma_y_new /= (n_y*n_x1);
				
			double [][] pars_prev = get_vectorized_pars(W_x, W_y, sigma_x, sigma_y);
			double [][] pars_new = get_vectorized_pars(W_x_new, W_y_new, sigma_x_new, sigma_y_new);
			
			int n_pars = pars_prev.length;
			double d = 0.0;
			for(int i=0; i<n_pars; i++) {
				d += Math.pow((pars_prev[i][0] - pars_new[i][0]),2.0);
			}
			
			//TODO: Does W_x not change when W_x_new is changed?
			W_x     = W_x_new;
			W_y     = W_y_new;
			sigma_x = sigma_x_new;
			sigma_y = sigma_y_new;
			
			n_iterations_done = j+1;
			if(d<=convergence_criterion){
				convergence = true;
				break;	
			}						
		}
		
		labeled_factors = MatrixOperations.transpose(Z1_trans);
		unlabeled_factors = MatrixOperations.transpose(Z2_trans);
		labeled_factor_covariance = A_inv;		
		unlabeled_factor_covariance = M;
	}
	
	
	public double calc_sq_euclidian_X_4_SSPPCA() {
		
		double norm = 0.0;
		
		for(int i=0; i<n_x1; i++) {
			for(int j=0; j<n_variables; j++) {
				norm += X_1_scaled[i][j]*X_1_scaled[i][j];
			}
		}
		for(int i=0; i<n_x2; i++) {
			for(int j=0; j<n_variables; j++) {
				norm += X_2_scaled[i][j]*X_2_scaled[i][j];
			}
		}
		return norm;
	}
	
	
	public double calc_sq_euclidian_Y_4_SSPPCA() {
		return calc_sq_euclidian_Y_4_SPPCA();
	}
	
	
	public HashMap<String, double [][]> get_init_pars_4_SSPPCA() {		
		HashMap<String, double [][]> init_pars = get_init_pars_4_SPPCA();
		return init_pars;
	}
	
}
