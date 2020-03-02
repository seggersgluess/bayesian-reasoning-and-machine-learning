package Kernels;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class KernelRidgeRegression {

	double [][] explained_variable;
	double [][] explaining_variables;
	
	int n_observations;
	int n_explaining_variables;
	
	String kernel = "linear";
	double penalty = 0.0;
	
	double [][] fitted_explained_variable;
	double [][] residuals;
	
	//weights in dual (kernel) space
	double [][] weights;
	
	//primal weights (similar to normal Ridge Regression)
	double [][] primal_weights;
	
	double [][] prediction;
	
	boolean useConstant;
	
	double R_squared = -1.0;
	
	public KernelRidgeRegression(double [][] y, double [][] X, double penalty, boolean useConstant) {
		
		if(penalty < 0.0) {
			throw new RuntimeException("Negative complexity penalty not allowed.");
		}
		
		this.penalty               = penalty;
		this.useConstant            = useConstant;
		this.explained_variable    = y;
		this.n_observations         = y.length;
		
		if(useConstant == true) {
			this.explaining_variables = MatrixOperations.cbind(X,MatrixOperations.unit_vector(n_observations));
			this.n_explaining_variables = explaining_variables[0].length;
		}else {
			this.explaining_variables  = X;
			this.n_explaining_variables = explaining_variables[0].length;
		}
	}
	
	
	public KernelRidgeRegression(double [][] y, double [][] X, double penalty, boolean useConstant, String kernel) {
		
		if(penalty < 0.0) {
			throw new RuntimeException("Negative complexity penalty not allowed.");
		}
		
		this.penalty               = penalty;
		this.useConstant            = useConstant;
		//Check of valid kernel is done later in construction of object from class KernelFunctions.
		this.kernel = kernel;
		
		this.explained_variable    = y;
		this.n_observations         = y.length;
		
		if(useConstant == true) {
			this.explaining_variables = MatrixOperations.cbind(X,MatrixOperations.unit_vector(n_observations));
			this.n_explaining_variables = explaining_variables[0].length;
		}else {
			this.explaining_variables  = X;
			this.n_explaining_variables = explaining_variables[0].length;
		}
	}
	
	
	//Usage of linear kernel
	public void fit() {
		
		if(kernel.contentEquals("linear") == false) {
			throw new RuntimeException("Invalid kernel defined. Only linear kernel allowed.");
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gramMatrix = kf.calc_gram_matrix(explaining_variables); 
			
		calc_fit_results_from_gram_matrix(gramMatrix);

	}
	
	
	//Usage of rmf kernel
	public void fit(double sigma) {
		
		if(kernel.contentEquals("rbf") == false) {
			throw new RuntimeException("Invalid kernel defined. Only RMF kernel allowed.");
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gramMatrix = kf.calc_gram_matrix(explaining_variables, sigma); 
		
		calc_fit_results_from_gram_matrix(gramMatrix);

	}
	
	
	//Usage of sigmoid kernel
	public void fit(double gamma, double alpha) {
		
		if(kernel.contentEquals("sigmoid") == false) {
			throw new RuntimeException("Invalid kernel defined. Only sigmoid kernel allowed.");
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gramMatrix = kf.calc_gram_matrix(explaining_variables, gamma, alpha); 
		
		calc_fit_results_from_gram_matrix(gramMatrix);

	}
	
	
	//Usage of polynomial kernel
	public void fit(double degree, double gamma, double alpha) {
		
		if(kernel.contentEquals("polynomial") == false) {
			throw new RuntimeException("Invalid kernel defined. Only polynomial kernel allowed.");
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gramMatrix = kf.calc_gram_matrix(explaining_variables, degree, gamma, alpha); 
		
		calc_fit_results_from_gram_matrix(gramMatrix);

	}
	
	
	//Usage of squared exponential (se) kernel
	public void fit(double [][] Sigma) {
		
		if(kernel.contentEquals("se") == false) {
			throw new RuntimeException("Invalid kernel defined. Only se kernel allowed.");
		}
		
		double [][] Sigma_s = Sigma;
		
		if(useConstant == true) {
			int n = Sigma.length;
			Sigma_s = new double [n+1][n+1];		
			Sigma_s = MatrixOperations.set_sub_matrix_to_matrix(Sigma_s, Sigma, 1, n, 1, n);
			Sigma_s[0][0] = 1.0e-10;
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gramMatrix = kf.calc_gram_matrix(explaining_variables, Sigma_s); 
		
		calc_fit_results_from_gram_matrix(gramMatrix);

	}
	
	
	public void calc_fit_results_from_gram_matrix(double [][] gramMatrix) {
		
		weights = new double [n_observations][n_observations];
		
		for(int i=0; i<n_observations; i++) {
			for(int j=0; j<n_observations; j++) {
				if(i==j) {
					weights[i][j] = gramMatrix[i][j] + penalty;
				}else {
					weights[i][j] = gramMatrix[i][j];
				}
			}		
		}
		
		weights = MatrixOperations.inverse(weights);
		weights = MatrixOperations.multiplication(weights, explained_variable);	
		
		fitted_explained_variable = MatrixOperations.multiplication(gramMatrix,weights);
		residuals = MatrixOperations.substract(explained_variable, fitted_explained_variable);	
		
	}
	
	
	public void calc_est_pars() {
		if(weights == null) {
			throw new RuntimeException("Kernel ridge regression not fitted yet.");
		}
		primal_weights = MatrixOperations.multiplication(MatrixOperations.transpose(explaining_variables), weights);
	}
	
	
	//Prediction with linear kernel
	public void predict(double [][] X) {
		
		double [][] X_s = X;
		
		if(useConstant == true) {
			int n_obs = X.length;
			X_s = MatrixOperations.cbind(X, MatrixOperations.unit_vector(n_obs));
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gram_block = kf.calc_off_diagonal_block_gram_matrix(X_s, explaining_variables);
		prediction = MatrixOperations.multiplication(gram_block, weights);
		
	}
	
	
	//Prediction with rmf kernel
	public void predict(double [][] X, double sigma) {
		
		if(kernel.contentEquals("rbf") == false) {
			throw new RuntimeException("Invalid kernel defined. Only RMF kernel allowed.");
		}
		
		double [][] X_s = X;
		
		if(useConstant == true) {
			int n_obs = X.length;
			X_s = MatrixOperations.cbind(X, MatrixOperations.unit_vector(n_obs));
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gram_block = kf.calc_off_diagonal_block_gram_matrix(X_s, explaining_variables, sigma);
		prediction = MatrixOperations.multiplication(gram_block, weights);
		
	}
	
	
	//Prediction with sigmoid kernel
	public void predict(double [][] X, double gamma, double alpha) {
		
		if(kernel.contentEquals("sigmoid") == false) {
			throw new RuntimeException("Invalid kernel defined. Only sigmoid kernel allowed.");
		}
		
		double [][] X_s = X;
		
		if(useConstant == true) {
			int n_obs = X.length;
			X_s = MatrixOperations.cbind(X, MatrixOperations.unit_vector(n_obs));
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gram_block = kf.calc_off_diagonal_block_gram_matrix(X_s, explaining_variables, gamma, alpha);
		prediction = MatrixOperations.multiplication(gram_block, weights);
		
	}
	
	
	//Prediction with polynomial kernel
	public void predict(double [][] X, double degree, double gamma, double alpha) {
		
		if(kernel.contentEquals("polynomial") == false) {
			throw new RuntimeException("Invalid kernel defined. Only polynomial kernel allowed.");
		}
		
		double [][] X_s = X;
		
		if(useConstant == true) {
			int n_obs = X.length;
			X_s = MatrixOperations.cbind(X, MatrixOperations.unit_vector(n_obs));
		}
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gram_block = kf.calc_off_diagonal_block_gram_matrix(X_s, explaining_variables, degree, gamma, alpha);
		prediction = MatrixOperations.multiplication(gram_block, weights);
		
	}
	
	
	//Prediction with squared exponential (se) kernel (-> Sigma n_vars x n_vars covariance matrix)
	public void predict(double [][] X, double [][] Sigma) {
		
		if(kernel.contentEquals("se") == false) {
			throw new RuntimeException("Invalid kernel defined. Only squared exponential (se) kernel allowed.");
		}
		
		int n = Sigma.length;
		
		if(n != Sigma[0].length) {
			throw new RuntimeException("No square Sigma matrix supplied for squared exponential kernel ridge regression.");
		}
		
		double [][] Sigma_s = Sigma;
		
		double [][] X_s = X;
		
		if(useConstant == true) {
			int n_obs = X.length;
			X_s = MatrixOperations.cbind(X, MatrixOperations.unit_vector(n_obs));
			
			Sigma_s = new double [n+1][n+1];		
			Sigma_s = MatrixOperations.set_sub_matrix_to_matrix(Sigma_s, Sigma, 1, n, 1, n);
			Sigma_s[0][0] = 1.0e-10;
		}
				
		KernelFunctions kf = new KernelFunctions(kernel);
		
		double [][] gram_block = kf.calc_off_diagonal_block_gram_matrix(X_s, explaining_variables, Sigma_s);
		prediction = MatrixOperations.multiplication(gram_block, weights);
		
	}
	
	
	public double [][] get_fitted_values() {
		if(fitted_explained_variable == null) {
			throw new RuntimeException("No estimation with kernel regression done yet.");
		}
		return fitted_explained_variable;
	}
	
	
	public double [][] get_residuals() {
		if(residuals == null) {
			throw new RuntimeException("No estimation with kernel regression done yet.");
		}
		return residuals;
	}
	
	
	public double [][] get_kernel_weights() {
		if(weights == null) {
			throw new RuntimeException("No estimation with kernel regression done yet.");
		}
		return weights;
	}
	
	
	public double [][] get_est_parameters() {
		
		if(primal_weights == null) {
			calc_est_pars();
		}
		
		if(useConstant == false) {
			return primal_weights;
		}else {
			int n_pars = n_explaining_variables-1;
			return MatrixOperations.get_sub_matrix_between_row_and_col_idxs(primal_weights, 1, n_pars, 0, 0);
		}	
	}
	
	
	public double get_est_constant() {
		
		if(primal_weights == null) {
			calc_est_pars();
		}
		
		if(useConstant == false) {
			System.out.println("Kernel ridge regression fitted without constant.");
			return 0.0;
		}else {
			return primal_weights[0][0];
		}	
	}
	
	
	public double get_R_squared() {
		
		if(R_squared == -1.0) {
			if(fitted_explained_variable == null) {
				throw new RuntimeException("Ridge regression not fitted yet.");
			}
			double mean = GeneralMath.mean(explained_variable);
			R_squared = 0.0;
			for(int i=0; i<n_observations; i++) {
				double diff = fitted_explained_variable[i][0]-mean;
				R_squared = diff*diff;
				diff = explained_variable[i][0]-mean;
				R_squared /= diff*diff;			
			}
			R_squared = 1.0-R_squared;			
			return R_squared;
		}else {
			return R_squared;
		}	
	}
	
	
	public String get_used_kernel() {
		return kernel;
	}
	
	
	public double [][] get_predictedValues() {
		if(prediction == null) {
			throw new RuntimeException("No prediction from kernel ridge regression done yet.");
		}
		return prediction;
	}
	
	
	public String [] get_valid_kernels() {
		KernelFunctions kf = new KernelFunctions(kernel);
		String [] validKernels = kf.getValidKernels();
		return validKernels;
	}
	
}
