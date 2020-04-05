package ComponentModels;

import java.util.ArrayList;
import java.util.HashMap;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class PLS extends ComponentModels{

	double [][] Y;
	double [][] Y_scaled;
	
	int n_y = 0;
	
	double [][] center_pars_y;
	double [][] scale_pars_y;
	
	int n_factors = 1;
	
	double [][] W;
	double [][] W_x;
	double [][] W_y;
	
	double [][] reg_coeffs;
	
	double [][] factors_x;
	double [][] factors_y;
	
	double [][] rotated_Y;
	
	int iterations = 500;
	double convergence_criterion = 1e-06;
	
	
	public PLS(double[][] X, double [][] Y, int n_factors, boolean center, boolean scale) {
		super(X, center, scale);
		
		if(n_factors<1 || n_factors>n_variables) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_variables + " allowed.");
		}
		
		this.n_factors = n_factors;
		this.Y = Y;
		n_y = Y[0].length;
		
		if(center == true) {
			center_pars_y = new double [1][n_y];
			if(scale == true) {
				scale_pars_y = new double [1][n_y];
			}		
			for(int i=0; i<n_y; i++) {
				double [][] y_col = MatrixOperations.get_column_vec_from_matrix(Y, i);
				center_pars_y[0][i] = GeneralMath.mean(y_col);
				double scale_par = 1.0;
				if(scale == true) {
					scale_pars_y[0][i] = GeneralMath.sd(y_col);
					scale_par = scale_pars_y[0][i];
				}
				for(int j=0; j<n_observations; j++) {
					Y_scaled[j][i] = (Y[j][i] - center_pars[0][i])/scale_par;
				} 
			}		
		}
		
		if(scale == true && center == false) {
			scale_pars_y = new double [1][n_y]; 
			for(int i=0; i<n_y; i++) {
				double [][] y_col = MatrixOperations.get_column_vec_from_matrix(Y, i);
				scale_pars_y[0][i] = GeneralMath.sd(y_col);	
				for(int j=0; j<n_observations; j++) {
					Y_scaled[j][i] = Y[j][i]/scale_pars_y[0][i];
				} 
			}
		}
	}

	
	public void do_PLS() {
		do_nipals_4_pls();
		calc_rotated_X();
		calc_rotated_Y();
	}
	
	
	private void do_nipals_4_pls() {
		
		ArrayList<HashMap<String, double [][]>> pars = new ArrayList<HashMap<String, double [][]>>();
		
		//TODO: Check if copys are correct and do not overwrite X_scaled/ Y_scaled...
		double [][] X_mod   = X_scaled;
		double [][] Y_mod   = Y_scaled;
		
		for(int i=0; i<n_factors; i++) {
			
			HashMap<String, double [][]> factor_pars = new HashMap<String, double [][]>();
			
			double [][] X_trans = MatrixOperations.transpose(X_mod);
			double [][] Y_trans = MatrixOperations.transpose(Y_mod);
			
			double [][] w   = new double [n_variables][1];
			double [][] w_x = new double [n_variables][1];		
			double [][] w_y = new double [n_y][1];
			double [][] coeff = new double [1][1];
			double [][] scores_x = new double [n_observations][1];
			double [][] scores_y = MatrixOperations.get_column_vec_from_matrix(Y_mod, 0);
			
			double [][] w_old = new double [n_variables][1];
			double [][] w_y_old = new double [n_y][1];
			
			for(int j=0; j<iterations; j++) {
				
				w = MatrixOperations.multiplication(X_trans, scores_y);
				double norm = MatrixOperations.euclidian(w);
				w = MatrixOperations.scalar_multiplication(1.0/norm, w);
				
				scores_x = MatrixOperations.multiplication(X_mod, w);
				
				w_y = MatrixOperations.multiplication(Y_trans, scores_x);
				norm = MatrixOperations.euclidian(w_y);
				w_y = MatrixOperations.scalar_multiplication(1.0/norm, w_y);
				
				scores_y = MatrixOperations.multiplication(Y_mod, w_y);
				
				double d = calc_diff_norm_4_convergence_check(w, w_old, w_y, w_y_old);				
				if(d<convergence_criterion) {
					break;
				}else {
					w_old = w;
					w_y_old = w_y;
				}	
			}
			
			w_x = MatrixOperations.multiplication(X_trans, scores_x);
			double s = Math.pow(MatrixOperations.euclidian(scores_x),2.0);
			w_x = MatrixOperations.scalar_multiplication(1.0/s, w_x);	
			
			X_mod = MatrixOperations.substract(X_mod, MatrixOperations.multiplication(scores_x, MatrixOperations.transpose(w_x)));
			Y_mod = MatrixOperations.substract(Y_mod, MatrixOperations.multiplication(scores_y, MatrixOperations.transpose(w_y)));
			
			coeff = MatrixOperations.scalar_multiplication(s,MatrixOperations.multiplication(MatrixOperations.transpose(scores_x), scores_y));
			
			factor_pars.put("w", w_x);
			factor_pars.put("w_x", w_x);
			factor_pars.put("w_y", w_y);
			factor_pars.put("coeff", coeff);
			factor_pars.put("scores_x", scores_x);
			factor_pars.put("scores_y", scores_y);
					
			pars.add(factor_pars);
			
		}
		
		set_pars_2_obj(pars);		
	}
	
	
	private double calc_diff_norm_4_convergence_check(double [][] w_new, double [][] w_old, double [][] w_y_new, double [][] w_y_old) {
		
		double diff_norm = 0.0;
		
		int n_pars = n_variables + n_y;
		
		double [][] pars_new = new double [n_pars][1];
		double [][] pars_old = new double [n_pars][1];
		
		int c = 0;
		for(int i=0; i<n_variables; i++) {
			pars_new[c][0] = w_new[i][0];
			pars_old[c][0] = w_old[i][0];
			c++;
		}
		for(int i=0; i<n_y; i++) {
			pars_new[c][0] = w_y_new[i][0];
			pars_old[c][0] = w_y_old[i][0];
			c++;
		}
		
		double [][] diff = MatrixOperations.substract(pars_new, pars_old);
		diff_norm = MatrixOperations.euclidian(diff);
		
		return diff_norm;
	}
	
	
	
	private void set_pars_2_obj(ArrayList<HashMap<String, double [][]>> pars) {
		
		W   = new double [n_variables][n_factors];
		W_x = new double [n_variables][n_factors];
		W_y = new double [n_y][n_factors];
		
		reg_coeffs = new double [n_factors][1];
		
		factors_x = new double [n_observations][n_factors];
		factors_y = new double [n_observations][n_factors];
		
		for(int i=0; i<n_factors; i++) {
			for(int j=0; j<n_variables; j++) {
				W[j][i] = pars.get(i).get("w")[j][0];
			}
			for(int j=0; j<n_variables; j++) {
				W_x[j][i] = pars.get(i).get("w_x")[j][0];
			}
			for(int j=0; j<n_y; j++) {
				W_y[j][i] = pars.get(i).get("w_y")[j][0];
			}
			reg_coeffs[i][0] = pars.get(i).get("coeff")[0][0];
			for(int j=0; j<n_observations; j++) {
				factors_x[j][i] = pars.get(i).get("scores_x")[j][0];
				factors_y[j][i] = pars.get(i).get("scores_y")[j][0];
			}
		}
	}
	
	
	private void calc_rotated_X() {		
		if(W_x == null) {
			throw new RuntimeException("No partial least squares (PLS) done yet. Do at first the estimation.");
		}
		rotated_X = MatrixOperations.multiplication(factors_x, MatrixOperations.transpose(W_x));
		if(center == true) {
			for(int i=0; i<n_factors; i++) {
				double s = 1.0;
				if(scale == true) {
					s = scale_pars[0][i];
				}
				for(int j=0; j<n_observations; j++) {
					rotated_X[j][i] = rotated_X[j][i]*s+center_pars[0][i];
				}
			}
		}	
	}
	
	
	private void calc_rotated_Y() {		
		if(W_y == null) {
			throw new RuntimeException("No partial least squares (PLS) done yet. Do at first the estimation.");
		}
		rotated_Y = MatrixOperations.multiplication(factors_y, MatrixOperations.transpose(W_y));
		if(center == true) {
			for(int i=0; i<n_factors; i++) {
				double s = 1.0;
				if(scale == true) {
					s = scale_pars_y[0][i];
				}
				for(int j=0; j<n_observations; j++) {
					rotated_Y[j][i] = rotated_Y[j][i]*s+center_pars_y[0][i];
				}
			}
		}	
	}
	
	
	public HashMap<String, double [][]> predict(double [][] X_new) {
		
		int n = X_new.length;
		
		HashMap<String, double [][]> predictions = new HashMap<String, double [][]>();
				
		double [][] pred_factors_x = MatrixOperations.multiplication(X_new, W);
		double [][] pred_x = MatrixOperations.multiplication(pred_factors_x , MatrixOperations.transpose(W_x));
		double [][] pred_y = new double [n][n_y];
		
		for(int i=0; i<n; i++) {
			for(int j=0; j<n_factors; j++) {
				for(int k=0; k<n_y; k++) {
					pred_y[i][k] += reg_coeffs[j][0]*pred_factors_x[i][j]*W_y[k][j];
				}			 
			}
		}
		
		predictions.put("Y", pred_y);
		predictions.put("X", pred_x);
		predictions.put("Z", pred_factors_x);
		
		return predictions;
	}
	
	//TODO: getter/ setter ...
	
	
}
