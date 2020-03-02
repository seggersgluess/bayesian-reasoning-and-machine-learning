package Kernels;

import Mathematics.DistanceMetrics;
import Mathematics.MatrixOperations;

public class KernelFunctions {

	public double [][] x1;
	public double [][] x2;
	
	public String kernel;
	
	
	public KernelFunctions(String kernel) {
		
		kernel = kernel.toLowerCase();
		String [] validKernels = getValidKernels();		
		int [] validIdx = Utilities.Utilities.get_idx(validKernels, kernel);
		if(validIdx[0] == -1) {
			throw new RuntimeException(kernel + " is not a valid kernel.");
		}
		this.kernel = kernel;
	}

	
	public KernelFunctions(double [][] x1, double [][] x2, String kernel) {
		
		kernel = kernel.toLowerCase();
		String [] validKernels = getValidKernels();		
		int [] validIdx = Utilities.Utilities.get_idx(validKernels, kernel);
		if(validIdx[0] == -1) {
			throw new RuntimeException(kernel + " is not a valid kernel.");
		}
		
		if(x1.length != x2.length) {
			throw new RuntimeException("Unequal length of supplied vector arguments.");
		}
		
		this.x1 = x1;
		this.x2 = x2;
		this.kernel = kernel;		
	}
	
	
	public double rbf_kernel(double sigma) {
		
		if(x1 == null || x2 == null) {
			throw new RuntimeException("No inputs x1 and x2 supplied yet.");
		}
		
		if(sigma<0.0) {
			throw new RuntimeException("Invalid sigma supplied for RBF kernel. Only values > 0.0 allowed");
		}
		
		DistanceMetrics distMetric = new DistanceMetrics(x1,x2);
		double d = distMetric.calcDistanceMetric("sqeuclidian");
		
		double k = Math.exp(-d/(2.0*sigma*sigma));
		
		return k;
	}
	
	
	public double squared_exp_kernel(double [][] Sigma) {
		
		if(x1 == null || x2 == null) {
			throw new RuntimeException("No inputs x1 and x2 supplied yet.");
		}
		
		int n = x1.length;
		
		if(Sigma.length != n || Sigma[0].length != n) {
			throw new RuntimeException("Invalid sigma matrix supplied for squared exponential kernel.");
		}
		
		double [][] d = new double [n][1];
		double [][] d_trans =  new double [1][n];
		
		for(int i=0; i<n; i++) {
			d[i][0] = x1[i][0]-x2[i][0];
			d_trans[0][i] = d[i][0];
		}
		
		double [][] sigmaInv = MatrixOperations.inverse(Sigma);
		
		double k = MatrixOperations.multiplication(MatrixOperations.multiplication(d_trans, sigmaInv),d)[0][0];
		k = Math.exp(-0.5*k);
		
		return k;
	}
	
	
	public double polynomial_kernel(double degree, double gamma, double alpha) {
		
		if(x1 == null || x2 == null) {
			throw new RuntimeException("No inputs x1 and x2 supplied yet.");
		}
		
		if(degree<1.0) {
			throw new RuntimeException("Invalid degree supplied for polynomial kernel. Only values > 0.0 valid.");
		}
		
		if(alpha<0.0) {
			throw new RuntimeException("Invalid alpha supplied for polynomial kernel. Only values > 0.0 valid.");
		}
		
		double k = MatrixOperations.multiplication(MatrixOperations.transpose(x1), x2)[0][0];
		k = Math.pow((gamma*k+alpha),degree);
		
		return k;
	}
	
	
	public double sigmoid_kernel(double gamma, double alpha) {
		
		if(x1 == null || x2 == null) {
			throw new RuntimeException("No inputs x1 and x2 supplied yet.");
		}
		
		if(alpha<0.0) {
			throw new RuntimeException("Invalid alpha supplied for sigmoid kernel. Only values > 0.0 valid.");
		}
		
		double k = MatrixOperations.multiplication(MatrixOperations.transpose(x1), x2)[0][0];
		k = Math.tanh(gamma*k+alpha);

		return k;
	}
	
	
	public double linear_kernel() {
		
		if(x1 == null || x2 == null) {
			throw new RuntimeException("No inputs x1 and x2 supplied yet.");
		}
		
		double k = MatrixOperations.multiplication(MatrixOperations.transpose(x1), x2)[0][0];
		
		return k;
	}
	
	
	public double calc_kernel() {
		
		double d = 0.0;
		
		if(kernel.contentEquals("linear")) {
			d = linear_kernel();
		}
		
		return d;
	}
	
	
	public double calc_kernel(double sigma) {
		
		double d = 0.0;
		
		if(kernel.contentEquals("rbf")) {
			d = rbf_kernel(sigma);
		}
		
		return d;	
	}
	
	
	public double calc_kernel(double gamma, double alpha) {
						
		double d = 0.0;
		
		if(kernel.contentEquals("sigmoid")) {
			d = sigmoid_kernel(gamma,alpha);
		}
	
		return d;
	}
	
	
	public double calc_kernel(double degree, double gamma, double alpha) {
		
		double d = 0.0;
		
		if(kernel.contentEquals("polynomial")) {
			d = polynomial_kernel(degree, gamma, alpha);
		}
	
		return d;
	}
	
	
	public double calc_kernel(double [][] Sigma) {
		
		double d = 0.0;
		
		if(kernel.contentEquals("se")) {
			d = squared_exp_kernel(Sigma);
		}
		
		return d;	
	}
	
	
	//Calculates Gram matrix for linear kernel
	public double [][] calc_gram_matrix(double [][] X) {
		
		int n_obs = X.length;
		
		double [][] gram_matrix = new double [n_obs][n_obs];
		
		for(int i=0; i<n_obs; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X, i);
			int idx = i+1;
			for(int j=0; j<idx; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X, j);
				double d = calc_kernel();
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates Gram matrix for rmf kernel
	public double [][] calc_gram_matrix(double [][] X, double sigma) {
		
		int n_obs = X.length;
		
		double [][] gram_matrix = new double [n_obs][n_obs];
		
		for(int i=0; i<n_obs; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X, i);
			int idx = i+1;
			for(int j=0; j<idx; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X, j);
				double d = calc_kernel(sigma);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	

	//Calculates Gram matrix for sigmoid kernel
	public double [][] calc_gram_matrix(double [][] X, double gamma, double alpha) {
		
		int n_obs = X.length;
		
		double [][] gram_matrix = new double [n_obs][n_obs];
		
		for(int i=0; i<n_obs; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X, i);
			int idx = i+1;
			for(int j=0; j<idx; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X, j);
				double d = calc_kernel(gamma, alpha);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}			
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates Gram matrix for polynomial kernel
	public double [][] calc_gram_matrix(double [][] X, double degree, double gamma, double alpha) {
		
		int n_obs = X.length;
		
		double [][] gram_matrix = new double [n_obs][n_obs];
		
		for(int i=0; i<n_obs; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X, i);
			int idx = i+1;
			for(int j=0; j<idx; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X, j);
				double d = calc_kernel(degree, gamma, alpha);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates Gram matrix for squared exponential (se) kernel
	public double [][] calc_gram_matrix(double [][] X, double [][] Sigma) {
		
		int n_obs = X.length;
		
		double [][] gram_matrix = new double [n_obs][n_obs];
		
		for(int i=0; i<n_obs; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X, i);
			int idx = i+1;
			for(int j=0; j<idx; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X, j);
				double d = calc_kernel(Sigma);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates n_1 x n_2 off diagonal block of Gram matrix with linear kernel  for supplied subsamples X_1 and X_2
	public double [][] calc_off_diagonal_block_gram_matrix(double [][] X_1, double [][] X_2) {
		
		int n1 = X_1.length;
		int n2 = X_2.length;
		
		if(X_1[0].length != X_2[0].length) {
			throw new RuntimeException("Unequal number of features in supplied subsamples X1 and X2.");
		}
		
		double [][] gram_matrix = new double [n1][n2];
		
		for(int i=0; i<n1; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X_1, i);
			for(int j=0; j<n2; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X_2, j);
				double d = calc_kernel();
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}			
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates n_1 x n_2 off diagonal block of Gram matrix with rmf kernel  for supplied subsamples X_1 and X_2
	public double [][] calc_off_diagonal_block_gram_matrix(double [][] X_1, double [][] X_2, double sigma) {
		
		int n1 = X_1.length;
		int n2 = X_2.length;
		
		if(X_1[0].length != X_2[0].length) {
			throw new RuntimeException("Unequal number of features in supplied subsamples X1 and X2.");
		}
		
		double [][] gram_matrix = new double [n1][n2];
		
		for(int i=0; i<n1; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X_1, i);
			for(int j=0; j<n2; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X_2, j);
				double d = calc_kernel(sigma);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates n_1 x n_2 off diagonal block of Gram matrix with sigmoid kernel  for supplied subsamples X_1 and X_2
	public double [][] calc_off_diagonal_block_gram_matrix(double [][] X_1, double [][] X_2, double gamma, double alpha) {
		
		int n1 = X_1.length;
		int n2 = X_2.length;
		
		if(X_1[0].length != X_2[0].length) {
			throw new RuntimeException("Unequal number of features in supplied subsamples X1 and X2.");
		}
		
		double [][] gram_matrix = new double [n1][n2];
		
		for(int i=0; i<n1; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X_1, i);
			for(int j=0; j<n2; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X_2, j);
				double d = calc_kernel(gamma, alpha);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}			
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates n_1 x n_2 off diagonal block of Gram matrix with polynomial kernel  for supplied subsamples X_1 and X_2
	public double [][] calc_off_diagonal_block_gram_matrix(double [][] X_1, double [][] X_2, double degree, double gamma, double alpha) {
		
		int n1 = X_1.length;
		int n2 = X_2.length;
		
		if(X_1[0].length != X_2[0].length) {
			throw new RuntimeException("Unequal number of features in supplied subsamples X1 and X2.");
		}
		
		double [][] gram_matrix = new double [n1][n2];
		
		for(int i=0; i<n1; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X_1, i);
			for(int j=0; j<n2; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X_2, j);
				double d = calc_kernel(degree, gamma, alpha);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	
	
	//Calculates n_1 x n_2 off diagonal block of Gram matrix with se kernel  for supplied subsamples X_1 and X_2
	public double [][] calc_off_diagonal_block_gram_matrix(double [][] X_1, double [][] X_2, double [][] Sigma) {
		
		int n1 = X_1.length;
		int n2 = X_2.length;
		
		if(X_1[0].length != X_2[0].length) {
			throw new RuntimeException("Unequal number of features in supplied subsamples X1 and X2.");
		}
		
		double [][] gram_matrix = new double [n1][n2];
		
		for(int i=0; i<n1; i++) {
			x1 = MatrixOperations.get_row_vec_from_matrix(X_1, i);
			for(int j=0; j<n2; j++) {
				x2 = MatrixOperations.get_row_vec_from_matrix(X_2, j);
				double d = calc_kernel(Sigma);
				gram_matrix[i][j] = d;
				if(i != j) {
					gram_matrix[j][i] = d;
				}				
			}
		}
		
		return gram_matrix;
	}
	
	
	public String [] getValidKernels() {
		String [] validKernels = new String [5];
		
		validKernels[0] = "rbf";
		validKernels[1] = "se";
		validKernels[2] = "polynomial";
		validKernels[3] = "sigmoid";
		validKernels[4] = "linear";
		
		
		return validKernels;
	}
	
}
