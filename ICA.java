package ComponentModels;

import java.util.ArrayList;
import java.util.function.BiFunction;

import Mathematics.MatrixOperations;

public class ICA extends ComponentModels {

	int n_factors = 1;
	
	String prior = null;
	
	boolean whitening = true;
	double [][] whitening_matrix;
	
	double [][] mixing_matrix;
	double [][] factors;
	
	int iterations = 500;
	int iterationsDone = 0;
	boolean convergence = false;
	double convergence_criterion = 1e-04;
	
	ArrayList<BiFunction <Double, Double, Double>>  g;
	ArrayList<BiFunction <Double, Double, Double>>  g_dev;
	
	
	public ICA(double [][] X, int n_factors, boolean center, boolean whitening) {
		
		super(X, center, false);
		
		int n_vars = X[0].length;

		if(n_factors<1 || n_factors>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		this.n_factors = n_factors;
		this.whitening = whitening;
		
		//parent class method with: center = center, scale = false (only centering is done -> see super() constructor!).
		scale_and_center_data();	
		if(whitening == true) {
			do_whitening();
		}
	}
	
	
	public ICA(double [][] X, int n_factors, String prior, boolean center, boolean whitening) {
		
		super(X, center, false);
		
		int n_vars = X[0].length;

		if(n_factors<1 || n_factors>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		//Check valid prior
		boolean valid_prior = is_valid_prior(prior);
		if(valid_prior == false) {
			throw new RuntimeException(prior + " is not a valid prior for ICA.");
		}
		
		this.n_factors = n_factors;
		this.whitening = whitening;

		this.prior = prior.toLowerCase();
		set_priors();
				
		//parent class method with: center = center, scale = false (only centering is done).
		scale_and_center_data();
		if(whitening == true) {
			do_whitening();
		}
	}
	
	
	public void do_whitening() {
		
		double [][] cov = MatrixOperations.multiplication(MatrixOperations.transpose(X), X);		
		whitening_matrix = MatrixOperations.get_inv_square_root_of_matrix(cov);		
		X_scaled = MatrixOperations.multiplication(X, whitening_matrix);	
	}
	
	
	public void do_ICA() {
			
		double [][] W = get_init_mixing_matrix_W();
		double [][] C = MatrixOperations.multiplication(MatrixOperations.transpose(X_scaled), X_scaled);
		
		for(int i=0; i<iterations; i++) {
			double [][] betas = new double [n_factors][1];
			double [][] alphas = new double [n_factors][n_factors];
			double [][] A = new double [n_factors][n_factors];
			double [][] g_dev_sum = new double [n_factors][1];
			double [][] Z = MatrixOperations.multiplication(X_scaled, MatrixOperations.transpose(W));
			
			if(g == null) {
				
				this.g = new ArrayList<BiFunction <Double, Double, Double>>(n_factors);
				this.g_dev = new ArrayList<BiFunction <Double, Double, Double>>(n_factors);
				
				for(int k=0; k<n_factors; k++) {
					double [][] z_k = MatrixOperations.get_column_vec_from_matrix(Z, k);
					double sel_criterion = calc_criterion_4_prior_selection(z_k);
					if(sel_criterion > 0.0) {
						this.g.add(ICA::g_superGaussian);
						this.g_dev.add(ICA::g_dev_superGaussian);
					}else {
						this.g.add(ICA::g_subGaussian);
						this.g_dev.add(ICA::g_dev_subGaussian);
					}
				}
			}

			for(int j=0; j<n_observations; j++) {
				double [][] z_j = MatrixOperations.get_row_vec_from_matrix(Z, j);
				double [][] z_j_trans = new double [1][n_factors];
				double [][] g_vec = new double [n_factors][1];				
				for(int k=0; k<n_factors; k++) {
					z_j_trans[0][k] = z_j[k][0];
					g_vec[k][0] = g.get(k).apply(z_j[k][0],0.0);
					betas[k][0] -= z_j[k][0]*g_vec[k][0];
					g_dev_sum[k][0] += g_dev.get(k).apply(z_j[k][0], 0.0);
				}
				A = MatrixOperations.add(A, MatrixOperations.multiplication(g_vec, z_j_trans));
			}
			
			for(int k=0; k<n_factors; k++) {
				betas[k][0] /= n_observations;
				g_dev_sum[k][0] /= n_observations;
				alphas[k][k] = -1.0/(betas[k][0]+g_dev_sum[k][0]);
				for(int l=0; l<n_factors; l++) {
					A[k][l] /= n_observations;
					if(k==l) {
						A[k][l] = betas[k][0]+A[k][l];
					}else {
						A[k][l] = A[k][l];
					}
				}
			}
			A = MatrixOperations.multiplication(A, W);
			A = MatrixOperations.multiplication(alphas, A);
			double [][] W_new = new double [n_factors][n_variables];
			W_new = MatrixOperations.add(W, A);
				
			double [][] W_W_t = MatrixOperations.multiplication(W_new,C);
			W_W_t = MatrixOperations.multiplication(W_W_t, MatrixOperations.transpose(W_new));
			double [][] inv_sqr_W_W_t = MatrixOperations.get_inv_square_root_of_matrix(W_W_t);
			W_new = MatrixOperations.multiplication(inv_sqr_W_W_t, W_new);			
			
			check_convergence(W, W_new);	
						
			W = W_new;
			
			iterationsDone = i+1;
			if(convergence == true) {					
				break;
			}	
		}	
		
		mixing_matrix = MatrixOperations.transpose(W);
		factors = MatrixOperations.multiplication(X_scaled, mixing_matrix);
		calc_rotated_X();
	}
	
	
	public double [][] get_init_mixing_matrix_W() {
		double [][] W = new double [n_factors][n_variables];
		for(int i=0; i<n_factors; i++) {
			for(int j=0; j<n_variables; j++) {
				if(i==j) {
					W[i][i] = 1.0;
				}else {
					W[i][j] = 0.01;
				}			
			}
		}
		return W;
	}
	
	
	public double [][] get_factors() {
		if(factors == null) {
			System.out.println("ICA factors not calculated yet.");
		}
		return factors;
	}
	
	public void check_convergence(double [][] W_prev, double [][] W_new) {
		convergence = true;
		//Check if (w_j^Prev)^T w_j^New converges to 1.0 (orthonormality)
		for(int k=0; k<n_factors; k++) {
			double norm = 0.0;
			for(int j=0; j<n_variables; j++) {
				norm += W_prev[k][j]*W_new[k][j];
			}
			double absDiff = Math.abs(norm-1.0);
			if(absDiff>convergence_criterion) {
				convergence = false;
			}
		}
	}
	
	
	public void calc_rotated_X() {
		
		if(factors == null) {
			throw new RuntimeException("Factor Analysis not done yet. No factors found.");
		}
		
		rotated_X = new double [n_observations][n_variables];
		
		for(int i=0; i<n_observations; i++) {

			double [][] factor = MatrixOperations.get_row_vec_from_matrix(factors, i);
			double [][] x_ik = MatrixOperations.multiplication(mixing_matrix, factor);
			for(int j=0; j<n_variables; j++) {
					rotated_X[i][j] = x_ik[j][0];
			}	
		}	
	}
	
	
	public double [][] calc_uncentered_and_unwhitened_rotated_X() {
		
		double [][] inv_whitening_matrix = MatrixOperations.inverse(whitening_matrix);
		double [][] unadjusted_rotated_X = MatrixOperations.multiplication(rotated_X, inv_whitening_matrix);
		for(int j=0; j<n_variables; j++) {
			for(int i=0; i<n_observations; i++) {
				unadjusted_rotated_X[i][j] -=center_pars[0][j];
			}		
		}		
		return unadjusted_rotated_X;
	}
	
	
	public double calc_criterion_4_prior_selection(double [][] x) {
		
		int n_obs = x.length;
		double criterion = 0.0;
		
		for(int i=0; i<n_obs; i++) {
			double tanh = Math.tanh(x[i][0]);
			criterion += -tanh*x[i][0]+(1.0-tanh*tanh);
		}
		criterion /=n_obs;
		
		return criterion;
	}
	
	
	public String [] get_valid_priors() {
		String [] validPriors = {"cosh",
				                 "cube",
				                 "exp"
							    };
		return validPriors;
	}
	
	
	public boolean is_valid_prior(String s_prior) {
		boolean valid_prior = true;
		String [] validPriors = get_valid_priors();
		int [] idxs = Utilities.Utilities.get_idx(validPriors,s_prior);
		if(idxs[0] == -1) {
			valid_prior = false;
		}
		
		return valid_prior;
	}
	
	
	public void set_priors() {
		this.g = new ArrayList<BiFunction <Double, Double, Double>>(n_factors);
		this.g_dev = new ArrayList<BiFunction <Double, Double, Double>>(n_factors);
		if(prior == "cosh") {
			for(int k=0; k<n_factors; k++) {
				g.add(ICA::g_cosh);
				g_dev.add(ICA::g_dev_cosh);
			}
		}
		
		if(prior == "cube") {
			for(int k=0; k<n_factors; k++) {
				g.add(ICA::g_cube);
				g_dev.add(ICA::g_dev_cube);
			}
		}
		
		if(prior == "exp") {
			for(int k=0; k<n_factors; k++) {
				g.add(ICA::g_exp);
				g_dev.add(ICA::g_dev_exp);
			}
		}		
	}
	
	
	public static double g_superGaussian(double x, double furtherArg) {
		return -2.0*Math.tanh(x);
	}
	
	
	public static double g_dev_superGaussian(double x, double futherArg) {
		return Math.tanh(x)-x;
	}
	
	
	public static double g_subGaussian(double x, double furtherArg) {
		return -2.0*(1.0-Math.tanh(x)*Math.tanh(x));
	}
	
	
	public static double g_dev_subGaussian(double x, double furtherArg) {
		return -Math.tanh(x)*Math.tanh(x);
	}
	
	
	//TODO: Check function (Murphy 2012 p. 415)!
	public static double g_cosh(double x, double furtherArg) {
		double c = Math.PI/(2.0*Math.sqrt(3));
		return c*Math.sinh(x*c)/Math.cosh(x*c);
	}
	
	
	//TODO: Check derivation of g_cosh!
	public static double g_dev_cosh(double x, double furtherArg) {
		double c = Math.PI/(2.0*Math.sqrt(3));
		return 2.0*Math.pow(c,3.0)*Math.cosh(x*c)/Math.pow(Math.cosh(x*c),2.0);
	}
	
	
	public static double g_cube(double x, double furtherArg) {
		return x*x*x;
	}
	
	
	public static double g_dev_cube(double x, double furtherArg) {
		return 3.0*x*x;
	}
	
	
	public static double g_exp(double x, double furtherArg) {
		return (1.0-x*x)*Math.exp(-x*x/2.0);
	}
	
	
	public static double g_dev_exp(double x, double furtherArg) {
		return -4.0*x*Math.exp(-x*x/2.0);
	}
	
	
	public int get_done_iterations() {
		return iterationsDone;
	}
	
}
