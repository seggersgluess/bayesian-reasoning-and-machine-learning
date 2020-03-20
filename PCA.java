package ComponentModels;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import Mathematics.MatrixOperations;

public class PCA extends ComponentModels{
	
	protected int n_pcs = 1;
	
	String method = "eigen";
	
	//Matrix of factor loadings W
	protected double [][] rotation_matrix;
	
	//Matrix of factors (or loadings) z
	protected double [][] factors;
		
	double [][] factor_variance;
	double [][] explained_variance_by_factors;
	double [][] cum_explained_variance_by_factors;
	
	//Covariance for eigendecompostion
	protected double [][] cov_matrix;
	
	public PCA(double [][] X, int n_pcs,  boolean scale) {
		
		super(X, scale);
		
		int n_vars = X[0].length;
		
		if(n_pcs<1 || n_pcs>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		this.n_pcs = n_pcs;

	}
	
	
	public PCA(double [][] X, int n_pcs, boolean center, boolean scale) {
		
		super(X, center, scale);
		
		int n_vars = X[0].length;

		if(n_pcs<1 || n_pcs>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		this.n_pcs = n_pcs;
	
	}
	
	
	public PCA(double [][] X, int n_pcs, String method, boolean center, boolean scale) {
		
		super(X, center, scale);
		
		int n_vars = X[0].length;
		
		if(n_vars > X.length) {
			throw new RuntimeException("Not enough observations supplied. (Number of features > number of observations)!");
		}
		
		if(n_pcs<1 || n_pcs>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		method = method.toLowerCase();
		String [] validMethods = get_valid_methods();
		int [] idxs = Utilities.Utilities.get_idx(validMethods, method);
		if(idxs[0] == -1) {
			throw new RuntimeException(method + " is not a valid method for PCA.");
		}
		
		this.n_pcs = n_pcs;
		this.method = method;
		
		scale_and_center_data();		
	}
	
	
	public void do_PCA() {
		
		if(method.contentEquals("eigen") == true) {
			//Here X is N times M (N: sample length, M number of variables)
			//-> so (XX^T) becomes here (X^T X)
			cov_matrix = MatrixOperations.multiplication(MatrixOperations.transpose(X_scaled), X_scaled);			
			do_PCA_by_EigenDec();
		}
		
		if(method.contentEquals("svd") == true) {
			do_PCA_by_SVD();
		}
		
		if(center == false && scale == false) {
			//free memory
			X_scaled = null;
		}		
	}
	
	
	public void do_PCA_by_SVD() {
		
		HashMap<String, double [][]> svd = MatrixOperations.get_fullSVD(X_scaled);
		double [][] right_singular_vectors = svd.get("V");
        double [][] singular_values        = svd.get("S");
		double [][] eigenValues            = new double [n_variables][1];
        
		//--- squared singular values ---
		for(int i=0; i<n_variables; i++) {
			eigenValues[i][0] = singular_values[i][0]*singular_values[i][0];
		}
		
		HashMap<String, List<Double>> sortedValues = Utilities.Utilities.get_sorted_elements_and_idxs_of_double_vector(eigenValues);
		List<Double> sortedEigenValues = sortedValues.get("SortedValues");
		List<Double> sortIdxs = sortedValues.get("Idxs");
		Collections.reverse(sortedEigenValues);
		Collections.reverse(sortIdxs);
		double [][] eigenVectors = right_singular_vectors;
		
		calc_pca_from_eigenvalues_and_eigenvectors(sortedEigenValues, eigenVectors, sortIdxs);		
	}
	

	public void do_PCA_by_EigenDec() {
		
		HashMap<String, double [][]> eigenDec = MatrixOperations.get_eigen_dec_4_symmetric_matrix(cov_matrix);
		
		double [][] eigenValues = MatrixOperations.get_diagonal_from_matrix(eigenDec.get("eigenvalues"));
		double [][] eigenVectors = eigenDec.get("eigenvectors");
		
		HashMap<String, List<Double>> sortedValues = Utilities.Utilities.get_sorted_elements_and_idxs_of_double_vector(eigenValues);
		List<Double> sortedEigenValues = sortedValues.get("SortedValues");
		List<Double> sortIdxs = sortedValues.get("Idxs");
		Collections.reverse(sortedEigenValues);
		Collections.reverse(sortIdxs);
		
		calc_pca_from_eigenvalues_and_eigenvectors(sortedEigenValues, eigenVectors, sortIdxs);
	}
	
	
	private void calc_pca_from_eigenvalues_and_eigenvectors(List<Double> sortedEigenValues, double [][] eigenVectors, List<Double> sortIdxs) {
		
		rotation_matrix = new double [n_variables][n_pcs];

		factor_variance = new double [1][n_pcs];
		explained_variance_by_factors = new double [1][n_pcs];
		cum_explained_variance_by_factors = new double [1][n_pcs];
		
		double cum_variance = 0.0;
		
		for(int i=0; i<n_variables; i++) {
			
			if(i<n_pcs) {
				int idx = sortIdxs.get(i).intValue();
				double [][] eigen_vec = MatrixOperations.get_column_vec_from_matrix(eigenVectors, idx);
				rotation_matrix = MatrixOperations.set_sub_matrix_to_matrix(rotation_matrix, eigen_vec, 0, (n_variables-1), i, i);
				//Variance with division by N-L where N is sample length and L is number of PCs
				factor_variance[0][i] = sortedEigenValues.get(i)/(n_observations-n_pcs);
			}	
			cum_variance += sortedEigenValues.get(i)/(n_observations-n_pcs);		
		}
				
		factors = MatrixOperations.multiplication(X_scaled, rotation_matrix);
		rotated_X = MatrixOperations.multiplication(factors, MatrixOperations.transpose(rotation_matrix));
		
		double cum_explained_variance = 0.0;
		
		for(int i=0; i<n_pcs; i++) {
			explained_variance_by_factors[0][i] = factor_variance[0][i]/cum_variance;
			cum_explained_variance += explained_variance_by_factors[0][i];
			cum_explained_variance_by_factors[0][i] = cum_explained_variance;
		}		
	}
	
	
	public String [] get_valid_methods() {
		String [] valid_methods = {"eigen",
								   "svd"};
		return valid_methods;
	}
	
	
	public double [][] get_factors() {
		if(factors == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			return factors;
		}
	}
	
	
	public double [][] get_rotation_matrix() {
		if(rotation_matrix == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			return rotation_matrix;
		}
	}
	
	
	public void clearCovMatrix() {
		cov_matrix = null;
	}


	public double [][] get_cov_matrix() {
		if(method == "eigen") {
			if(cov_matrix != null) {
				return cov_matrix;
			}else {
				System.out.println("No covariance matrix found.");
				return null;
			}
		}else {
			System.out.println("Method for PCA not eigendecomposition. No covariance found.");
			return null;
		}
	}
	
	
	public String get_method() {
		return method;
	}
	
	
	public int get_number_of_pcs() {
		return n_pcs;
	}
	

	public double [][] get_variance_of_factors() {
		if(factor_variance == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			return factor_variance;
		}
	}
	
	
	public double [][] get_sd_of_factors() {
		if(factor_variance == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			
			double [][] factor_sd = new double [1][n_pcs];
			for(int i=0; i<n_pcs; i++) {
				factor_sd[0][i] = Math.sqrt(factor_variance[0][i]);
			}
			return factor_sd;
		}
	}
	
	
	public double [][] get_explained_variance_of_factors() {
		if(explained_variance_by_factors == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			return explained_variance_by_factors;
		}
	}
	
	
	public double [][] get_cum_explained_variance_of_factors() {
		if(cum_explained_variance_by_factors == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			return cum_explained_variance_by_factors;
		}
	}
	
}
