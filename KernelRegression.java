package Kernels;

import java.util.HashMap;

import Mathematics.MatrixOperations;

public class KernelRegression {

	double [][] explained_variable;
	double [][] explaining_variables;
	
	int n_observations;
	int n_explaining_variables;
	
	String kernel = "gaussian";
	String metric = "euclidian";
	
	double bandwidth = 0.5;
	
	double [][] fitted_explained_variable;
	double [][] residuals;
	double [][] weights;
	
	double [][] prediction;
	double [][] weights4Prediction;
	
	
	public KernelRegression(double [][] y, double [][] X) {
		
		this.explained_variable    = y;
		this.explaining_variables  = X;
		
		this.n_explaining_variables = X[0].length;
		this.n_observations         = y.length;
	}
	
	
	public KernelRegression(double [][] y, double [][] X, double bandwidth) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
		
		this.bandwidth = bandwidth;
		this.explained_variable    = y;
		this.explaining_variables  = X;
		
		this.n_explaining_variables = X[0].length;
		this.n_observations         = y.length;
	}
	
	
	public KernelRegression(double [][] y, double [][] X, double bandwidth, String kernel) {
		
		if(bandwidth <= 0.0) {
			throw new RuntimeException("Invalid bandwidth. Only values larger than 0 valid.");
		}
			
		this.bandwidth = bandwidth;
		this.kernel = kernel;
		this.explained_variable    = y;
		this.explaining_variables  = X;
		
		this.n_explaining_variables = X[0].length;
		this.n_observations         = y.length;
	}
	
	
	public HashMap<String, double [][]> calcPredictionAndWeights(double [][] X) {
		
		int n = X.length;
		
		KernelDensities kde = new KernelDensities(explaining_variables, bandwidth, kernel, metric);
		HashMap<String, double [][]> estRes = kde.calcKernelDensity4NewObsWithDensityComponents(X);
		
		double [][] estDensities = estRes.get("densities");
		double [][] components   = estRes.get("components");
		
		double [][] predictions = new double [n][1];
		double [][] estWeights  = new double [n][n_observations];
		
		for(int i=0; i<n; i++) {
			for(int j=0; j<n_observations; j++) {
				estWeights[i][j] = components[i][j]/(n_observations*estDensities[i][0]);
				predictions[i][0] += estWeights[i][j]*explained_variable[j][0];
			}		
		}
		
		HashMap<String, double [][]> predRes = new HashMap<String, double [][]>();
		predRes.put("prediction", predictions);
		predRes.put("weights", estWeights);
		
		return predRes;	
	}
	
	
	public void predict(double [][] X) {
		
		HashMap<String, double [][]> predRes = calcPredictionAndWeights(X);
		
		prediction = predRes.get("prediction");
		weights4Prediction = predRes.get("weights");
	}
	
	
	public void inSamplePredictionKernelRegression() {
		
		HashMap<String, double [][]> predRes = calcPredictionAndWeights(explaining_variables);
		
		fitted_explained_variable = predRes.get("prediction");
		residuals = MatrixOperations.substract(explained_variable, fitted_explained_variable);
		weights = predRes.get("weights");		
	}
	
	
	public String [] getValidKernels() {
		
		double [][] dummyArray = new double [1][1];
		KernelDensities kd = new KernelDensities(dummyArray, 0.5);
		String [] validKernels = kd.getValidKernels();
		return validKernels;
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
	
	
	public double [][] get_weights() {
		if(weights == null) {
			throw new RuntimeException("No estimation with kernel regression done yet.");
		}
		return weights;
	}
	
	
	public double [][] get_predictedValues() {
		if(prediction == null) {
			throw new RuntimeException("No prediction from kernel regression done yet.");
		}
		return prediction;
	}
	
	
	public double [][] get_weightsOfPrediction() {
		if(weights4Prediction == null) {
			throw new RuntimeException("No prediction from kernel regression done yet.");
		}
		return weights4Prediction;
	}
	
	
	public String get_used_kernel() {
		return kernel;
	}
	
}
