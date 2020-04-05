package ComponentModels;

import java.util.HashMap;

import Distributions.NormalDistribution;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class SPPCA extends ComponentModels{
	
	double [][] Y;
	double [][] Y_scaled;
	
	int n_y = 0;
	
	int n_factors;
	
	double [][] W_x;
	double [][] W_y;
	double [][] my_x;
	double [][] my_y;
	
	//sigma_x^2 & sigma_y^2
	double sigma_x;
	double sigma_y;
	
	double [][] factors;
	double [][] factor_covariance;
	
	int iterations = 500;
	int n_iterations_done = 0;
	
	boolean convergence = false;
	double convergence_criterion = 1e-06;
	
	double [][] rotated_Y;
	
	HashMap<String, double [][]> ext_init_values;
	
	
	public SPPCA(double[][] X, double [][] Y, int n_factors, boolean center) {
		//Centering the data is necessary.
		super(X, true, false);
			
		if(n_factors<1 || n_factors>n_variables) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_variables + " allowed.");
		}
		
		this.n_factors = n_factors;
		this.Y = Y;
		this.n_y = Y[0].length;
		
		my_x = MatrixOperations.transpose(center_pars);
		
		Y_scaled = new double [n_observations][n_y];
		my_y = get_sample_means(Y);
		for(int k=0; k<n_y; k++) {
			double [][] y = MatrixOperations.get_column_vec_from_matrix(Y, k);
			for(int i=0; i<n_observations; i++) {
				Y_scaled[i][k] = y[i][0]-my_y[k][0];
			}
		}
	}
	
	
	public void em_4_SPPCA() {
		
		HashMap<String, double [][]> pars = new HashMap<String, double [][]>();
		if(ext_init_values == null) {
			pars = get_init_pars_4_SPPCA();
		}else {
			pars = ext_init_values;
		}
		
		W_x  = pars.get("W_x");
		W_y  = pars.get("W_y");
		sigma_x = pars.get("sigma_x")[0][0];
		sigma_y = pars.get("sigma_y")[0][0];
			
		double [][] Z_trans = new double [n_factors][n_observations];
		double [][] A_inv   = new double [n_factors][n_factors];
		double [][] diag_L = MatrixOperations.identity(n_factors);
			
		double normX = calc_sq_euclidian_X_4_SPPCA();
		double normY = calc_sq_euclidian_Y_4_SPPCA();
		
		for(int j=0; j<iterations; j++) {	
											
			double [][] W_x_trans = MatrixOperations.transpose(W_x); 
			double [][] W_y_trans = MatrixOperations.transpose(W_y);
			double [][] A = MatrixOperations.scalar_multiplication(1.0/sigma_x, MatrixOperations.multiplication(W_x_trans, W_x));
			A = MatrixOperations.add(A,MatrixOperations.scalar_multiplication(1.0/sigma_y, MatrixOperations.multiplication(W_y_trans, W_y)));
			A = MatrixOperations.add(A, diag_L);
				
			A_inv = MatrixOperations.inverse(A);
				
			double [][] W_x_new = new double [n_variables][n_factors];
			double [][] W_y_new = new double [n_y][n_factors];
			double [][] Cov_z_i_sum = new double [n_factors][n_factors];
				
			double [][] expected_z_i_trans = new double [1][n_factors];
			
			for(int i=0; i<n_observations; i++) {								
				double [][] x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
				double [][] y_i = MatrixOperations.get_row_vec_from_matrix(Y_scaled, i);
				
				double [][] expected_z_i = calc_expected_z_i_4_SPPCA(W_x_trans, W_y_trans, x_i, y_i);
					
				for(int k=0; k<n_factors; k++) {
					expected_z_i_trans[0][k] = expected_z_i[k][0];
					Z_trans[k][i] = expected_z_i[k][0];
				}
	
				double [][] expected_z_i_z_i_trans = calc_cov_z_i_4_SPPCA(A_inv, expected_z_i, expected_z_i_trans);
									 
				W_x_new = MatrixOperations.add(W_x_new,MatrixOperations.multiplication(x_i,expected_z_i_trans));
				Cov_z_i_sum = MatrixOperations.add(Cov_z_i_sum, expected_z_i_z_i_trans);

				W_y_new = MatrixOperations.add(W_y_new,MatrixOperations.multiplication(y_i,expected_z_i_trans));				
			}
				
			double [][] Cov_z_i_sum_inv = MatrixOperations.inverse(Cov_z_i_sum);
			W_x_new = MatrixOperations.multiplication(W_x_new, Cov_z_i_sum_inv);
			double [][] W_x_new_trans = MatrixOperations.transpose(W_x_new);
			W_y_new = MatrixOperations.multiplication(W_y_new, Cov_z_i_sum_inv);
			double [][] W_y_new_trans = MatrixOperations.transpose(W_y_new);	
			
			double sigma_x_new = 0.0;
			double sigma_y_new = 0.0;
			
			sigma_x_new += normX;
			sigma_y_new += normY;
						
			sigma_x_new += MatrixOperations.trace(MatrixOperations.multiplication(MatrixOperations.multiplication(W_x_new_trans, W_x_new),Cov_z_i_sum));
			sigma_x_new -= 2.0*MatrixOperations.trace(MatrixOperations.multiplication(MatrixOperations.multiplication(X_scaled,W_x_new),Z_trans));
			sigma_x_new /= (n_variables*n_observations);
			
			sigma_y_new += MatrixOperations.trace(MatrixOperations.multiplication(MatrixOperations.multiplication(W_y_new_trans, W_y_new),Cov_z_i_sum));
			sigma_y_new -= 2.0*MatrixOperations.trace(MatrixOperations.multiplication(MatrixOperations.multiplication(Y_scaled,W_y_new),Z_trans));
			sigma_y_new /= (n_y*n_observations);
				
			double [][] pars_prev = get_vectorized_pars(W_x, W_y, sigma_x, sigma_y);
			double [][] pars_new = get_vectorized_pars(W_x_new, W_y_new, sigma_x_new, sigma_y_new);
			
			int n_pars = pars_prev.length;
			double d = 0.0;
			for(int i=0; i<n_pars; i++) {
				d += Math.pow((pars_prev[i][0] - pars_new[i][0]),2.0);
			}
			
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
		
		factors = MatrixOperations.transpose(Z_trans);
		factor_covariance = A_inv;		
	}
	
	
	public double [][] calc_expected_z_i_4_SPPCA(double [][] W_x_trans, double [][] W_y_trans, double [][] x_i, double [][] y_i) {
		
		double [][] expected_z_i = MatrixOperations.scalar_multiplication(1/sigma_x, MatrixOperations.multiplication(W_x_trans, x_i));
		expected_z_i = MatrixOperations.add(expected_z_i, MatrixOperations.scalar_multiplication(1/sigma_y, MatrixOperations.multiplication(W_y_trans, y_i)));
		
		return expected_z_i;
	}
	
	
	public double [][] calc_cov_z_i_4_SPPCA(double [][] A_inv, double [][] expected_z_i, double [][] expected_z_i_trans) {
		double [][] expected_z_i_z_i_trans = MatrixOperations.add(A_inv, MatrixOperations.multiplication(expected_z_i, expected_z_i_trans));		
		return expected_z_i_z_i_trans;
	}
	
	
	public double calc_sq_euclidian_Y_4_SPPCA() {
		
		double norm = 0.0;
		
		for(int i=0; i<n_observations; i++) {
			for(int j=0; j<n_y; j++) {
				norm += Y_scaled[i][j]*Y_scaled[i][j];
			}
		}
		return norm;
	}
	
	
	public double calc_sq_euclidian_X_4_SPPCA() {
		
		double norm = 0.0;
		
		for(int i=0; i<n_observations; i++) {
			for(int j=0; j<n_variables; j++) {
				norm += X_scaled[i][j]*X_scaled[i][j];
			}
		}
		return norm;
	}
	
	
	public double [][] get_vectorized_pars(double [][] W_x, double [][] W_y, double sigma_x, double sigma_y) {
		
		int n=n_variables*n_factors;
		n+= n_y*n_factors+2;
		
		double [][] par_vec = new double [n][1];
		int c = 0;
		for(int i=0; i<n_variables; i++) {
			for(int j=0; j<n_factors; j++) {
				par_vec[c][0] = W_x[i][j];
				c++;
			}
		}
		for(int i=0; i<n_y; i++) {
			for(int j=0; j<n_factors; j++) {
				par_vec[c][0] = W_y[i][j];
				c++;
			}
		}
		par_vec[c][0] = sigma_x;
		par_vec[(c+1)][0] = sigma_y;
		
		return par_vec;
	}
	
	
	public double [][] get_sample_means(double [][] X) {
		
		int n_vars = X[0].length;
		double [][] means = new double [n_vars][1];
		
		for(int i=0; i<n_vars; i++) {
			double [][] x_i = MatrixOperations.get_column_vec_from_matrix(X, i);
			means[i][0] = GeneralMath.mean(x_i);
		}
		return means;
	}
	
	
	public void do_SPPCA() {
		em_4_SPPCA();
		calc_rotated_X();	
		calc_rotated_Y();
	}
	
	
	public void calc_rotated_X() {
		
		if(W_x == null) {
			throw new RuntimeException("No rotation matrix W_x calculated yet. Do SPPCA estimation at first.");
		}
		rotated_X = MatrixOperations.multiplication(factors,MatrixOperations.transpose(W_x));
		for(int i=0; i<n_observations; i++) {
			for(int j=0; j<n_variables; j++) {
				rotated_X[i][j] += my_x[j][0];
			}
		}
	}
	
	
	public void calc_rotated_Y() {
		if(W_y == null) {
			throw new RuntimeException("No rotation matrix W_y calculated yet. Do SPPCA estimation at first.");
		}
		rotated_Y = MatrixOperations.multiplication(factors, MatrixOperations.transpose(W_y));
		for(int i=0; i<n_observations; i++) {
			for(int j=0; j<n_y; j++) {
				rotated_Y[i][j] += my_y[j][0];
			}
		}
	}
	
	
	public HashMap<String, double [][]> predict(double [][] X_new) {
		
		if(X_new[0].length != n_variables) {
			throw new RuntimeException("Invalid input data for prediction. Only " + n_variables + " input features allowed");
		}
		
		if(W_x == null) {
			throw new RuntimeException("SPPCA not trained yet. No prediction possible.");
		}
		
		HashMap<String, double [][]> pred_res = new HashMap<String, double [][]>();
		
		double [][] W_x_trans = MatrixOperations.transpose(W_x);
		double [][] Psi = new double [n_factors][n_factors];
		
		for(int i=0; i<n_factors; i++) {
			Psi[i][i] = sigma_x;
		}
		
		double [][] M = MatrixOperations.multiplication(W_x_trans, W_x);
		M = MatrixOperations.add(M, Psi);
		M = MatrixOperations.inverse(M);
		M = MatrixOperations.multiplication(M, W_x_trans);
		
		int n = X_new.length;
		
		double [][] Z_pred = new double [n][n_factors];
		
		for(int i=0; i<n; i++) {
			double [][] x = MatrixOperations.get_row_vec_from_matrix(X_new, i);
			double [][] mean_adj_x = MatrixOperations.substract(x, my_x);
			double [][] z = MatrixOperations.multiplication(M, mean_adj_x);
			for(int j=0; j<n_factors; j++) {
				Z_pred[i][j] = z[j][0];
			}
		}
		
		double [][] X_pred = MatrixOperations.multiplication(Z_pred, MatrixOperations.transpose(W_x));
		for(int i=0; i<n; i++) {
			for(int j=0; j<n_variables; j++) {
				X_pred[i][j] += my_x[j][0];
			}
		}
		
		double [][] Y_pred = MatrixOperations.multiplication(Z_pred, MatrixOperations.transpose(W_y));
		for(int i=0; i<n; i++) {
			for(int j=0; j<n_y; j++) {
				Y_pred[i][j] += my_y[j][0];
			}
		}
		
		pred_res.put("Z", Z_pred);
		pred_res.put("X",X_pred);
		pred_res.put("Y", Y_pred);
		
		return pred_res;
	}
	
	
	public double [][] get_factors() {
		if(factors == null) {
			throw new RuntimeException("No factors calculated yet. Do SPPCA estimation at first.");
		}
		return factors;
	}
	
	
	public double [][] get_W_x() {
		if(W_x == null) {
			throw new RuntimeException("No rotation matrix W_x calculated yet. Do SPPCA estimation at first.");
		}
		return W_x;
	}
	
	
	public double [][] get_W_y() {
		if(W_y == null) {
			throw new RuntimeException("No rotation matrix W_y calculated yet. Do SPPCA estimation at first.");
		}
		return W_y;
	}
	
	
	public double [][] get_rotated_Y() {
		if(rotated_Y == null) {
			throw new RuntimeException("No rotated y calculated yet. Do SPPCA estimation at first.");
		}
		return rotated_Y;
	}
	
	
	public double get_sigma_x() {
		if(W_x == null) {
			throw new RuntimeException("No noise level sigma_x calculated yet. Do SPPCA estimation at first.");
		}
		return sigma_x;
	}
	
	
	public double get_sigma_y() {
		if(W_y == null) {
			throw new RuntimeException("No noise level sigma_y calculated yet. Do SPPCA estimation at first.");
		}
		return sigma_y;
	}
	
	
	public double [][] get_my_x() {
		if(W_x == null) {
			throw new RuntimeException("No mean my_x calculated yet. Do SPPCA estimation at first.");
		}
		return my_x;
	}
	
	
	public double [][] get_my_y() {
		if(W_y == null) {
			throw new RuntimeException("No mean my_y calculated yet. Do SPPCA estimation at first.");
		}
		return my_y;
	}
	
	
	public int get_number_of_factors() {
		return n_factors;
	}
	
	
	public boolean isConverged() {
		return convergence;
	}
	
	
	public int get_numberOfIterationsDone() {
		return n_iterations_done;
	}
	
	
	public void set_iterations4EM(int iterations) {
		if(iterations < 1) {
			throw new RuntimeException("Invalid iterations for EM supplied.");
		}
		this.iterations = iterations;
	}
	
	
	public void set_convergence_criterion4EM(double criterion) {
		this.convergence_criterion = criterion;
	}
	
	
	public void set_external_init_pars(HashMap<String, double [][]> init_pars) {
		ext_init_values = init_pars;
	}
	
	
	public void set_external_init_pars(double [][] W_x, double [][] sigma_x, double [][] W_y, double [][] sigma_y) {
		
		ext_init_values = new HashMap<String, double [][]>();			
		ext_init_values.put("W_x", W_x);
		ext_init_values.put("sigma_x", sigma_x);		
		ext_init_values.put("W_y", W_y);
		ext_init_values.put("sigma_y", sigma_y);		
	}
	
	
	public HashMap<String, double [][]> get_init_pars_4_SPPCA() {
		
		HashMap<String, double [][]> init_pars = new HashMap<String, double [][]>();
		
		double [][] init_sigma_x = new double [n_variables][1];
		double [][] init_sigma_y = new double [n_y][1];
		
		init_sigma_x[0][0] = 1.0;
		init_sigma_y[0][0] = 1.0;
		
		double [][] init_W_x = new double [n_variables][n_factors];
		double [][] init_W_y = new double [n_y][n_factors];
		
		for(int i=0; i<n_variables; i++) {
			for(int j=0; j<n_factors; j++) {
				NormalDistribution normDist = new NormalDistribution(0.0, 1.0);
				init_W_x[i][j] = normDist.sample()[0][0];
			}
		}
		
		for(int i=0; i<n_y; i++) {
			for(int j=0; j<n_factors; j++) {
				NormalDistribution normDist = new NormalDistribution(0.0, 1.0);
				init_W_y[i][j] = normDist.sample()[0][0];
			}
		}
		
		init_pars.put("W_x", init_W_x);
		init_pars.put("sigma_x", init_sigma_x);
		init_pars.put("W_y", init_W_y);
		init_pars.put("sigma_y", init_sigma_y);
			
		return init_pars;
	}
	
}


