package Kernels;

import Mathematics.DistanceMetrics;
import Mathematics.MatrixOperations;

public class KernelDensities {

	public double [][] X;
	public int nObservations = 0;
	
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
		this.bandwidth = bandwidth;
		this.distMetric = new DistanceMetrics(X,X);
	}
	
	
	KernelDensities(double [][] X, double bandwidth, String kernel) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
				
		String [] validKernels = getValidKernels();
		int [] validIdx = Utilities.Utilities.get_idx(validKernels, kernel);
		if(validIdx[0] == -1) {
			throw new RuntimeException(kernel + " is not a valid kernel.");
		}
		
		this.X = X;
		this.nObservations = X.length;
		this.bandwidth = bandwidth;
		this.kernel = kernel;
		this.distMetric = new DistanceMetrics(X,X);
	}
	
	
	KernelDensities(double [][] X, double bandwidth, String kernel, String metric) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
		
		String [] validMetrics = getValidDistanceMetrics();
		int [] validIdx = Utilities.Utilities.get_idx(validMetrics, metric);
		if(validIdx[0] == -1) {
			throw new RuntimeException(metric + " is not a valid distance metric.");
		}
		
		String [] validKernels = getValidKernels();
		validIdx = Utilities.Utilities.get_idx(validKernels, kernel);
		if(validIdx[0] == -1) {
			throw new RuntimeException(kernel + " is not a valid kernel.");
		}
		
		this.X = X;
		this.nObservations = X.length;
		this.bandwidth = bandwidth;
		this.kernel = kernel;
		this.metric = metric;
		this.distMetric = new DistanceMetrics(X,X);
	}
	
	
	public void calcKernelDensity() {
		
		density = new double [nObservations][0];
		
		double summedDensities = 0.0;
		
		for(int i=0; i<nObservations; i++) {
			double [][] x_i = MatrixOperations.get_row_vec_from_matrix(X, i);
			for(int j=0; j<nObservations; j++) {
				double [][] x_j = MatrixOperations.get_row_vec_from_matrix(X, j);
				distMetric.setVectors(x_i,x_j);
				double d = distMetric.calcDistanceMetric(metric);
				density[i][0] += calcKernel(d);
			}
			summedDensities += density[i][0];
		}
		
		for(int i=0; i<nObservations; i++) {
			density[i][0] /= summedDensities;
		}		
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
			double b_sq = bandwidth*bandwidth;
			fittedDensity = (1.0-distance*distance/b_sq);
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
