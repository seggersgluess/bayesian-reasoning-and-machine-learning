package Kernels;

import java.util.HashMap;

import Mathematics.DistanceMetrics;
import Mathematics.MatrixOperations;

public class KernelDensities {

	public double [][] X;
	public int nObservations = 0;
	public int n_explaining_variables = 0;
	
	public String kernel = "gaussian";
	public String metric = "euclidian";
	
	public double bandwidth = 0.0;
	
	public DistanceMetrics distMetric;
	
	public double [][] density;
	
	KernelDensities(double [][] X, double bandwidth) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
			
		this.X = X;
		this.nObservations = X.length;
		this.n_explaining_variables = X[0].length;
		this.bandwidth = bandwidth;
		this.distMetric = new DistanceMetrics(X,X);
		
		if(n_explaining_variables>1) {
			this.metric = "sqeuclidian";
		}
	}
	
	
	KernelDensities(double [][] X, double bandwidth, String kernel) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
			
		kernel = kernel.toLowerCase();
		String [] validKernels = getValidKernels();		
		int [] validIdx = Utilities.Utilities.get_idx(validKernels, kernel);
		if(validIdx[0] == -1) {
			throw new RuntimeException(kernel + " is not a valid kernel.");
		}
		
		this.X = X;
		this.nObservations = X.length;
		this.n_explaining_variables = X[0].length;
		this.bandwidth = bandwidth;
		this.kernel = kernel;
		this.distMetric = new DistanceMetrics(X,X);
		
		if(n_explaining_variables>1) {
			this.metric = "sqeuclidian";
		}
	}
	
	
	KernelDensities(double [][] X, double bandwidth, String kernel, String metric) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
		
		this.X = X;
		this.nObservations = X.length;
		this.n_explaining_variables = X[0].length;
		this.bandwidth = bandwidth;
		this.distMetric = new DistanceMetrics(X,X);
		
		metric.toLowerCase();
		String [] validMetrics = getValidDistanceMetrics();
		int [] validIdx = Utilities.Utilities.get_idx(validMetrics, metric);
		if(validIdx[0] == -1) {
			throw new RuntimeException(metric + " is not a valid distance metric.");
		}
		
		kernel = kernel.toLowerCase();
		String [] validKernels = getValidKernels();
		validIdx = Utilities.Utilities.get_idx(validKernels, kernel);
		if(validIdx[0] == -1) {
			throw new RuntimeException(kernel + " is not a valid kernel.");
		}
		
		this.kernel = kernel;
		
		if(n_explaining_variables>1 && metric.contentEquals("sqeuclidian") == false) {
			System.out.println("For multivariate input data X only squared euclidian is a valid metric.");
			this.metric = "sqeuclidian";
		}else {
			this.metric = metric;
		}	
	}
	
	
	public void calcKernelDensity() {
		
		density = new double [nObservations][1];
		
		int n_vars = X[0].length;
		
		for(int i=0; i<nObservations; i++) {
			double [][] x_i = MatrixOperations.get_row_vec_from_matrix(X, i);
			for(int j=0; j<nObservations; j++) {
				double [][] x_j = MatrixOperations.get_row_vec_from_matrix(X, j);
				distMetric.setVectors(x_i,x_j);
				double d = distMetric.calcDistanceMetric(metric);
				density[i][0] += calcKernel(d);			
			}
			density[i][0] /= nObservations;		
			if(n_vars>1) {
				density[i][0] -= 1.0/nObservations;
			}
		}	
				
	}
	
	
	public double [][] calcKernelDensity4NewObs(double [][] new_obs_X) {
		
		int n = new_obs_X.length;	
		int n_vars = new_obs_X[0].length;
		
		double [][] estDensities = new double [n][1];
				
		for(int i=0; i<n; i++) {
			double [][] new_x = MatrixOperations.get_row_vec_from_matrix(new_obs_X, i);
			for(int j=0; j<nObservations; j++) {
				double [][] x_j = MatrixOperations.get_row_vec_from_matrix(X, j);
				distMetric.setVectors(new_x,x_j);
				double d = distMetric.calcDistanceMetric(metric);
				estDensities[i][0] += calcKernel(d);			
			}
			estDensities[i][0] /= nObservations;
			if(n_vars>1) {
				estDensities[i][0] -= 1.0/nObservations;
			}
		}	
		
		return estDensities;
	}
	
	
	public HashMap<String, double [][]> calcKernelDensity4NewObsWithDensityComponents(double [][] new_obs_X) {
		
		int n = new_obs_X.length;
		int n_vars = new_obs_X[0].length;
		
		double [][] estDensities     = new double [n][1];
		double [][] compEstDensities = new double [n][nObservations];		
		
		for(int i=0; i<n; i++) {
			double [][] new_x = MatrixOperations.get_row_vec_from_matrix(new_obs_X, i);
			for(int j=0; j<nObservations; j++) {
				double [][] x_j = MatrixOperations.get_row_vec_from_matrix(X, j);
				distMetric.setVectors(new_x,x_j);
				double d = distMetric.calcDistanceMetric(metric);
				compEstDensities[i][j] = calcKernel(d);
				estDensities[i][0] += compEstDensities[i][j];			
			}
			estDensities[i][0] /= nObservations;
			if(n_vars>1) {
				estDensities[i][0] -= 1.0/nObservations;
			}
		}	
		
		HashMap<String, double [][]> estRes = new HashMap<String, double [][]>(2);
		estRes.put("densities", estDensities);
		estRes.put("components", compEstDensities);
		
		return estRes;
	}
	
	
	public double calcKernel(double distance) {
		
		double fittedDensity = 0.0;
		
		if(kernel.contentEquals("gaussian")) {
			double b_sq = bandwidth*bandwidth;
			fittedDensity = Math.exp(-distance*distance/(2*b_sq));
		}
		
		if(kernel.contentEquals("tophat")) {
			if(distance<bandwidth) {
				fittedDensity = 1.0;
			}
		}
		
		if(kernel.contentEquals("epanechnikov")) {
			double x = distance/bandwidth;
			if(Math.abs(x)<=1.0) {
				fittedDensity = 0.75*(1.0-x*x);
			}		
		}
		
		if(kernel.contentEquals("exponential")) {
			fittedDensity = Math.exp(-distance/bandwidth);
		}
		
		if(kernel.contentEquals("linear")) {
			if(distance<bandwidth) {
				fittedDensity = 1.0-distance/bandwidth;
			}		
		}
		
		if(kernel.contentEquals("cosine")) {
			if(distance<bandwidth) {
				double x = Math.PI*distance/2*bandwidth;
				fittedDensity = Math.cos(x);
			} 
		}
		
		return fittedDensity;
	}
	
	
	public void setAdditionalDistancePars(double [][] w) {
		
		if(metric.contentEquals("seuclidian")) {
			distMetric.setInput4StandardizedEuclidian(w);
		}
		
		if(metric.contentEquals("mahalanobis")) {
			distMetric.setInput4Mahalanobis(w);
		}
	}
	
	
	public void setAdditionalDistancePars(double p) {	
		distMetric.setInput4Minkowski(p);
	}
	
	
	public void setAdditionalDistancePars(double p, double [][] w) {		
		distMetric.setInput4WeightedMinkowski(p,w);
	}
	
	
	public String [] getValidDistanceMetrics() {
		String [] validMetrics = distMetric.getListOfDistanceMetrics();
		return validMetrics;
	}
	
	
	public String [] getValidKernels() {
		String [] validKernels = new String [6];
		
		validKernels[0] = "gaussian";
		validKernels[1] = "tophat";
		validKernels[2] = "epanechnikov";
		validKernels[3] = "exponential";
		validKernels[4] = "linear";
		validKernels[5] = "cosine";
		
		return validKernels;
	}
	
	
	public double [][] getFittedKernelDensities() {
		if(density == null) {
			System.out.println("No kernel densities fitted yet.");
			return null;
		}
		return density;
	}
	
}
