package Kernels;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import ComponentModels.PCA;
import Mathematics.MatrixOperations;

public class KernelPCA extends PCA{

	String kernel = "linear";
	
	HashMap<String, double [][]> addPars4Kernel = new HashMap<String, double [][]>();
	
	public KernelPCA(double[][] X, int n_pcs, boolean scale) {
		super(X, n_pcs, scale);
	}

	
	public KernelPCA(double [][] X, int n_pcs, boolean center, boolean scale) {
		super(X, n_pcs, center, scale);
	}
	
	
	public KernelPCA(double [][] X, int n_pcs, String kernel, boolean center, boolean scale) {
		super(X, n_pcs, center, scale);
		this.kernel = kernel;
	}
	
	
	public KernelPCA(double [][] X, int n_pcs, String kernel, double sigma, boolean center, boolean scale) {
		super(X, n_pcs, center, scale);
		this.kernel = kernel;
		double [][] s = MatrixOperations.convNumberToArray(sigma);
		this.addPars4Kernel.put("sigma", s);
	}
	
	
	public KernelPCA(double [][] X, int n_pcs, String kernel, double gamma, double alpha, boolean center, boolean scale) {
		super(X, n_pcs, center, scale);
		this.kernel = kernel;
		double [][] g = MatrixOperations.convNumberToArray(gamma);
		double [][] a = MatrixOperations.convNumberToArray(alpha);
		this.addPars4Kernel.put("gamma", g);
		this.addPars4Kernel.put("alpha", a);
	}
	
	
	public KernelPCA(double [][] X, int n_pcs, String kernel, double degree, double gamma, double alpha, boolean center, boolean scale) {
		super(X, n_pcs, center, scale);
		this.kernel = kernel;
		double [][] d = MatrixOperations.convNumberToArray(degree);
		double [][] g = MatrixOperations.convNumberToArray(gamma);
		double [][] a = MatrixOperations.convNumberToArray(alpha);
		this.addPars4Kernel.put("degree", d);
		this.addPars4Kernel.put("gamma", g);
		this.addPars4Kernel.put("alpha", a);
	}
	
	
	public KernelPCA(double [][] X, int n_pcs, String kernel, double [][] Sigma, boolean center, boolean scale) {
		super(X, n_pcs, center, scale);
		this.kernel = kernel;
		this.addPars4Kernel.put("Sigma", Sigma);
	}
	
	
	public void do_kernelPCA() {
		
		KernelFunctions kf = new KernelFunctions(kernel);
		
		if(kernel.contentEquals("linear")) {
			cov_matrix = kf.calc_gram_matrix(X_scaled);
		}
		if(kernel.contentEquals("rbf")) {
			checkKernelPars();
			double sigma = addPars4Kernel.get("sigma")[0][0];
			cov_matrix = kf.calc_gram_matrix(X_scaled, sigma);
		}
		if(kernel.contentEquals("sigmoid")) {
			checkKernelPars();
			double gamma = addPars4Kernel.get("gamma")[0][0];
			double alpha = addPars4Kernel.get("alpha")[0][0];
			cov_matrix = kf.calc_gram_matrix(X_scaled, gamma, alpha);
		}
		if(kernel.contentEquals("polynomial")) {
			checkKernelPars();
			double degree = addPars4Kernel.get("degree")[0][0];
			double gamma = addPars4Kernel.get("gamma")[0][0];
			double alpha = addPars4Kernel.get("alpha")[0][0];
			cov_matrix = kf.calc_gram_matrix(X_scaled, degree, gamma, alpha);
		}
		if(kernel.contentEquals("se")) {
			checkKernelPars();
			double [][] Sigma = addPars4Kernel.get("Sigma");
			cov_matrix = kf.calc_gram_matrix(X_scaled, Sigma);
		}
		
		double scaleFactor = 1.0/((double) n_observations);
		double [][] I_N = MatrixOperations.get_matrix_with_equal_elments(scaleFactor, n_observations, n_observations);

		double [][] m_1 = MatrixOperations.multiplication(I_N, cov_matrix);
		double [][] m_2 = MatrixOperations.multiplication(cov_matrix, I_N);
		double [][] m_3 = MatrixOperations.multiplication(m_1, I_N);
		cov_matrix = MatrixOperations.substract(cov_matrix, m_1);
		cov_matrix = MatrixOperations.substract(cov_matrix, m_2);
		cov_matrix = MatrixOperations.add(cov_matrix, m_3);
		
		do_kernelPCA_by_EigenDec();
		
		if(center == false && scale == false) {
			//free memory
			X_scaled = null;
		}
	}
	
	
	public void do_kernelPCA_by_EigenDec() {
		
		HashMap<String, double [][]> eigenDec = MatrixOperations.get_eigen_dec_4_symmetric_matrix(cov_matrix);
		
		double [][] eigenValues = MatrixOperations.get_diagonal_from_matrix(eigenDec.get("eigenvalues"));
		double [][] eigenVectors = eigenDec.get("eigenvectors");
		
		HashMap<String, List<Double>> sortedValues = Utilities.Utilities.get_sorted_elements_and_idxs_of_double_vector(eigenValues);
		List<Double> sortedEigenValues = sortedValues.get("SortedValues");
		List<Double> sortIdxs = sortedValues.get("Idxs");
		Collections.reverse(sortedEigenValues);
		Collections.reverse(sortIdxs);
		
		calc_kernelPCA_from_eigenvalues_and_eigenvectors(sortedEigenValues, eigenVectors, sortIdxs);
	}
	
	
	private void calc_kernelPCA_from_eigenvalues_and_eigenvectors(List<Double> sortedEigenValues, double [][] eigenVectors, List<Double> sortIdxs) {
		
		int n = eigenVectors.length;
		rotation_matrix = new double [n][n_pcs];
		
		for(int i=0; i<n; i++) {		
			if(i<n_pcs) {
				int idx = sortIdxs.get(i).intValue();
				double [][] eigen_vec = MatrixOperations.get_column_vec_from_matrix(eigenVectors, idx);
				double eigen_vec_scaling = 1.0/Math.sqrt(sortedEigenValues.get(i));
				eigen_vec = MatrixOperations.scalar_multiplication(eigen_vec_scaling, eigen_vec);
				rotation_matrix = MatrixOperations.set_sub_matrix_to_matrix(rotation_matrix, eigen_vec, 0, (n-1), i, i);
			}	
		}
				
		factors = MatrixOperations.multiplication(cov_matrix, rotation_matrix);			
	}
	
	
	public void checkKernelPars() {
		ArrayList<double [][]> kernelPars4Check = new ArrayList<double [][]>(addPars4Kernel.values());
		int n=kernelPars4Check.size();
		
		if(n==0) {
			throw new RuntimeException("No parameters supplied for " + kernel + " kernel.");
		}
		
		for(int i=0; i<n; i++) {
			if(kernelPars4Check.get(i)==null) {
				throw new RuntimeException("No parameters supplied for " + kernel + " kernel.");
			}
		}
	}
	
	public String get_used_kernel() {
		return kernel;
	}
	
}
