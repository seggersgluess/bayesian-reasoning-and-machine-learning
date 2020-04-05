package Clustering;

import java.util.ArrayList;
import java.util.HashMap;

import Mathematics.DistanceMetrics;

public class Clustering {

	double [][] X;
	int n_observations = 0;
	int n_variables    = 0;
	
	ArrayList<HashMap<String, double [][]>> inputStructure;
	
	String metric = "euclidian";
	double e = 0.0;
	int minPts = 0;
	
	int n_clusters = 0;
	ArrayList<ArrayList<Integer>> clusterMemberIdxs;
	ArrayList<Integer> noiseIdxs;
	
	
	public Clustering(double [][] X, double epsilon, int minPoints) {
		
		if(epsilon<=0.0) {
			throw new RuntimeException("No valid epsilon supplied. Only positive values allowed.");
		}
		if(minPoints<=0) {
			throw new RuntimeException("No valid min. number of points supplied. Only positive values allowed.");
		}
		if(X == null){
			throw new RuntimeException("No input data supplied.");
		}
		
		this.X = X;
		this.n_observations = X.length;
		this.n_variables = X[0].length;
		this.e = epsilon;
		this.minPts = minPoints;
		
	}
	
	
	ArrayList<HashMap<String, double [][]>> get_neighborhood(HashMap<String, double [][]> point) {
		
		ArrayList<HashMap<String, double [][]>>  neighborhood = new ArrayList<HashMap<String, double [][]>>();
		
		double [][] x = point.get("features");
		for(int i=0; i<n_observations; i++) {
				double [][] y = inputStructure.get(i).get("features");
				double dist = new DistanceMetrics(x,y).calcDistanceMetric(metric);
				if(dist <= e) {
					neighborhood.add(inputStructure.get(i));				
				}	
		}
		return neighborhood;
	}
	
	
	public ArrayList<HashMap<String, double [][]>> get_noise_points() {
		
		ArrayList<HashMap<String, double [][]>> noise = new ArrayList<HashMap<String, double [][]>>();
		
		int n = noiseIdxs.size();
		if(n==0) {
			System.out.println("No noise points found.");
			noise = null;
		}else {
			for(int i=0; i<n; i++) {
				noise.add(inputStructure.get(noiseIdxs.get(i)));
			}
		}
		return noise;
	}
	
	
	//Returns features X for identified noise points
	public double [][] get_X_of_noise() {
		
		double [][] noise_X;
		
		int n = noiseIdxs.size();
		if(n==0) {
			System.out.println("No noise points found.");
			noise_X = null;
		}else {
			noise_X = new double [n][n_variables];
			for(int i=0; i<n; i++) {
				for(int j=0; j<n_variables; j++) {
					noise_X[i][j] = inputStructure.get(noiseIdxs.get(i)).get("features")[j][0];
				}
			}
		}
		return noise_X;
	}
	
	
	public ArrayList<HashMap<String, double [][]>> get_cluster_points(int clusterNumber) {
		
		if(clusterNumber>=n_clusters) {
			throw new RuntimeException("Invalid cluster number supplied.");
		}
		
		ArrayList<HashMap<String, double [][]>> cluster = new ArrayList<HashMap<String, double [][]>>();
		
		if(n_clusters==0) {
			System.out.println("No noise points found.");
			cluster = null;
		}else {
			int n_pointsInCluster = clusterMemberIdxs.get(clusterNumber).size();
			for(int i=0; i<n_pointsInCluster; i++) {
				cluster.add(inputStructure.get(clusterMemberIdxs.get(clusterNumber).get(i)));
			}
		}
		return cluster;
	}
	
	
	//Returns features X for identified points in clusterNumber (0,1,...n_clusters-1)
	public double [][] get_X_of_cluster(int clusterNumber) {
		
		if(clusterNumber>=n_clusters) {
			throw new RuntimeException("Invalid cluster number supplied.");
		}
		double [][] cluster_X;
		if(n_clusters==0) {
			System.out.println("No noise points found.");
			cluster_X = null;
		}else {
			int n_pointsInCluster = clusterMemberIdxs.get(clusterNumber).size();
			cluster_X = new double [n_pointsInCluster][n_variables];
			for(int i=0; i<n_pointsInCluster; i++) {
				for(int j=0; j<n_variables; j++) {
					cluster_X[i][j] = inputStructure.get(clusterMemberIdxs.get(clusterNumber).get(i)).get("features")[j][0];
				}
			}
		}
		return cluster_X;
	}
	
	
	public String [] get_valid_metrics() {
		String [] metrics = new String [5];
		metrics[0] = "euclidian";
		metrics[1] = "manhatten";
		metrics[2] = "chebyshev";
		metrics[3] = "cosine";
		metrics[4] = "sqeuclidian";
		return metrics;
	}
	
	
	public void set_metric(String metric) {
		metric = metric.toLowerCase();
		int [] idxs = Utilities.Utilities.get_idx(get_valid_metrics(), metric);
		if(idxs[0] == -1) {
			throw new RuntimeException(metric + " is not a valid metric for OPTICS clustering.");
		}
		this.metric = metric;
	}
	
	
	public String get_used_distance_metric() {
		return metric;
	}
	
	
	public int get_number_of_clusters() {
		return n_clusters;	
	}
	
}
