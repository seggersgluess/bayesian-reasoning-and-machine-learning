package Distributions;

import Mathematics.GeneralMath;

public class Poisson {

	double [][] x;
	double lambda;
	
	Poisson(double x, double lambda) {
		
		if(x < 0.0) {
			throw new RuntimeException("Invalid x supplied. Only positive integer allowed.");
		}
		
		if(lambda <= 0.0) {
			throw new RuntimeException("Invalid mean supplied. Only positive values allowed.");
		}
		
		this.x = new double [1][1];
		this.x[0][0] = x;
		this.lambda = lambda;
	}
	
	
	Poisson(double [][] x, double lambda) {
		
		int n = x.length;
		
		for(int i=0; i<n; i++) {
			if(x[i][0] < 0.0) {
				throw new RuntimeException("Invalid x supplied. Only positive integer allowed.");
			}
		}
		
		if(lambda <= 0.0) {
			throw new RuntimeException("Invalid mean supplied. Only positive values allowed.");
		}
		
		this.x = x;
		this.lambda = lambda;
	}
	
	
	public double [][] pmf() {
			
		int n = x.length;
		double [][] probs = new double [n][1]; 
			
		for(int i=0; i<n; i++) {
			double x_fac = GeneralMath.factorial(x[i][0]);			
			probs[i][0] = Math.pow(lambda, x[i][0])*Math.exp(-lambda)/x_fac;
		}
		
		return probs; 
	}
	
	
	public double mean() {
		return lambda;
	}
	
	
	public double variance() {
		return lambda;
	}
	
}
