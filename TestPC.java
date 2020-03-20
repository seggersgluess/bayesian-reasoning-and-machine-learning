package ComponentModels;

import DataManagement.InputDataManager;
import Kernels.KernelPCA;
import Mathematics.MatrixOperations;

public class TestPC {

	//Compare with prinComp in R (Attention: 1/N is used for Covariance XX^T in R)
	//(-> Eigenvalues are scaled by 1/N)
	
	public static void test1_PCA() {
		
		double [][] data = get_iris_data();
		
		//MatrixOperations.print_matrix(data);
		
		long startTime = System.currentTimeMillis();
		PCA pca = new PCA(data, 1, "eigen", true, true); 
		pca.do_PCA();
		long endTime = System.currentTimeMillis();
		System.out.println("PCA done after " + ((endTime-startTime)/1000.0) + " secs.");
		
		//MatrixOperations.print_matrix(pca.get_center_pars());
		
		double [][] x_rotated = pca.get_rotated_input();
		double [][] x_scaled  = pca.get_scaled_input_data();
		
		System.out.println("Diff. (org.) unreduced data vs. (approx.) reduced data");
		MatrixOperations.print_matrix(MatrixOperations.substract(x_scaled, x_rotated));
		
		System.out.println("");
		System.out.println("Rotation Matrix:");
		MatrixOperations.print_matrix(pca.get_rotation_matrix());
		System.out.println("");
		System.out.println("PCA Summary:");
		MatrixOperations.print_matrix(pca.get_sd_of_factors());
		MatrixOperations.print_matrix(pca.get_explained_variance_of_factors());
		MatrixOperations.print_matrix(pca.get_cum_explained_variance_of_factors());
		
	}
	

	public static void test2_kernelPCA() {
		
		double [][] data = get_iris_data();
		
		//MatrixOperations.print_matrix(data);
		
		long startTime = System.currentTimeMillis();
		KernelPCA pca = new KernelPCA(data, 3, "rbf", 1.0, true, false); 
		pca.do_kernelPCA();
		long endTime = System.currentTimeMillis();
		System.out.println("PCA done after " + ((endTime-startTime)/1000.0) + " secs.");
		
		System.out.println("");
		System.out.println("Rotation Matrix:");
		MatrixOperations.print_matrix(pca.get_rotation_matrix());
		System.out.println("");
	}
	
	
	public static double [][] get_iris_data() {
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/DecisionTrees/Classification/IrisData.txt";
		InputDataManager inputData = new InputDataManager();	
		try {
			inputData.fileReader(dirName, false, true, true);
		} catch (Exception e) {
			System.out.println("Can´t find data set for upload.");
			e.printStackTrace();
		}
		
		String [][] input = inputData.strFileData;
		
		int n_rows = input.length;
		int n_cols = input[0].length-1;
		double [][] data = new double [n_rows][n_cols];
		for(int i=0; i<n_cols; i++) {
			for(int j=0; j<n_rows; j++) {
				data[j][i] =  Double.parseDouble(input[j][i]);
			}			
		}
		return data;
	}
	
	
	public static void main(String[] args) throws Exception {		
		//test1_PCA();
		test2_kernelPCA();
	}
	
}
