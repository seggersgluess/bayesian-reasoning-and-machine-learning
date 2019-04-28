package MixtureModels;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class k_means {

	static InputDataManager inputData;
	
	static double [][] observed_variables;
	
	static int n_observations;
	static int n_variables;

	static int n_clusters;
	
	static double [][] my;
	static double [][] clusterIdxs4Observations;
	static ArrayList<ArrayList<List<Double>>> observedClusterData;
	
	static int max_iterations;
	static double convergence_criterion;
	static int convCount = 10;
	
	//k-means algorithm stats
	static boolean convergence_reached = false;
	static int n_iterations;
	
	
	public static void do_k_means_clustering(){
		
		int iterationCount = 0;
		
		clusterIdxs4Observations = new double [n_observations][1];
		
		//Initialize cluster centers my
		do_k_means_initizalization();
		
		double prevEuclidian      = Double.MAX_VALUE;
		double [][] prevClusterIdxs = new double [n_observations][1];
		
		for(int i=0; i<max_iterations; i++){
			
			//Asignment step
			for(int j=0; j<n_observations; j++){
				
				double [][] x    = MatrixOperations.get_row_vec_from_matrix(observed_variables, j);
				
				for(int c=0; c<n_clusters; c++){
					
					double [][] my_k = MatrixOperations.get_column_vec_from_matrix(my,c);					
					double curEuclidian = MatrixOperations.euclidian(MatrixOperations.substract(x, my_k));
					curEuclidian        = Math.pow(curEuclidian, 2.0);
					
					if(curEuclidian<prevEuclidian){
						prevEuclidian = curEuclidian;
						clusterIdxs4Observations[j][0] = c;
					}
					
				}
	
			}
			
			//Update step
			for(int c=0; c<n_clusters; c++){
				
				double [][] myUpdate = new double [n_variables][n_clusters];
				int counts = 0;
				
				for(int j=0; j<n_observations; j++){
					
					if(clusterIdxs4Observations[j][0] == c){
						
						double [][] x = MatrixOperations.get_row_vec_from_matrix(observed_variables, j);
						myUpdate      = MatrixOperations.add(myUpdate, x);
						counts++;
						
					}
					
				}
				
				if(counts != 0){
					myUpdate = MatrixOperations.scalar_multiplication(1/counts, myUpdate);
				}
				
				for(int j=0; j<n_variables; j++){
					my[j][c] = myUpdate[j][0];
				}				
				
			}
			
			double sumDiff = GeneralMath.sum(MatrixOperations.substract(clusterIdxs4Observations, prevClusterIdxs));
				
			if(sumDiff == convergence_criterion){
				
				if(iterationCount == convCount){
					System.out.println("k-means clustering has converged after " + i + " iterations.");
					convergence_reached = true;
					n_iterations = i;
					break;
				}

				iterationCount++;
				
			}
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		for(int c=0; c<n_clusters; c++){
			
			ArrayList<List<Double>> observedDateOfCluster = new ArrayList<List<Double>>();
			
			for(int i=0; i<n_observations; i++){
				
				if(clusterIdxs4Observations[i][0] == c){
					
					List<Double> obs = new ArrayList<Double>(n_observations);
					double [][] x = MatrixOperations.get_row_vec_from_matrix(observed_variables, i);
					
					for(int j=0; j<n_variables; j++){
						obs.add(x[j][0]);
					}
					
					observedDateOfCluster.add(obs);
					
				}

			}
			
			observedClusterData.add(observedDateOfCluster);
			
		}
			
	}
	
	
	public static void do_k_means_initialization(){
		
		//Forgy method
		int [] randomIdxs = Utilities.getRandomIntNumbers(0, n_observations-1, n_clusters);
		
		my = new double [n_variables][n_clusters];
		
		for(int c=0; c<n_clusters; c++){
			for(int i=0; i<n_variables; i++){
				my[i][c] = observed_variables[randomIdxs[c]][i];
			}
		}
		
	}
	
	
	public static double [][] get_observed_data_4_cluster(int clusterNumber){
		
		if(clusterNumber>(n_clusters-1)){
			throw new RuntimeException("Invalid cluster number supplied for getting observed cluster data.");
		}
		
		int n = observedClusterData.get(clusterNumber).size();
		
		double [][] obsClusterData = new double [n_variables][n];
		
		for(int i=0; i<n; i++){			
			for(int j=0; j<n_variables; j++){
				obsClusterData[j][i] = observedClusterData.get(clusterNumber).get(i).get(j);				
			}			
		}
		
		return obsClusterData;
		
	}	
	
	
	@SuppressWarnings("static-access")
	public static void read_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception{
		
		inputData = new InputDataManager();
		
		inputData.fileReader(fileName, false, hasRowNames, hasColNames);
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void select_input_data(String [] rownames, String [] colnames, String ref_class){
		
		inputData.selectLoadedData(rownames, colnames);
		
		observed_variables = inputData.selectedDblFileData;		
		n_observations     = observed_variables.length;
		n_variables        = observed_variables[0].length;
		
	}
	
	
	public static void remove_input_data_manager(){	
		
		inputData = null;		
		
	}
	
	
}
