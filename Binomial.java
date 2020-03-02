package Distributions;

import Mathematics.GeneralMath;

public class Binomial {

	double n;
	double [][] x;
	double p;
	
	Binomial(double x, double n, double p) {
		
		this.x = new double [1][1];
		this.x[0][0] = x;
		
		this.n = n;
		this.p = p;
		
		checkInput();	
	}
	
	
	Binomial(double [][] x, double n, double p) {
		
		this.x = x;		
		this.n = n;
		this.p = p;
		
		checkInput();	
	}
	
	
	public double [][] pmf() {
		
		int n_vars = x.length;
		double [][] probs = new double [n_vars][1]; 
		
		double n_fac = GeneralMath.factorial(n);
		
		for(int i=0; i<n_vars; i++) {
			probs[i][0] = n_fac/(GeneralMath.factorial(n-x[i][0])*GeneralMath.factorial(x[i][0]));
			probs[i][0] *= Math.pow(p, x[i][0])*Math.pow((1.0-p), (n-x[i][0]));
		}
		
		return probs; 
	}
	
	
	public double mean() {
		return n*p;
	}
	
	
	public double variance() {
		return n*p*(1.0-p);
	}
	
	
	public void checkInput() {
		
		if(p<0.0 || p>1.0) {
			throw new RuntimeException("Invalid probability supplied. Only values between 0 and 1 allowed.");
		}
			
		if(n<=0.0) {
			throw new RuntimeException("Invalid n supplied. Only values > 0 allowed.");
		}
		
		int nVars = x.length;
		
		for(int i=0; i<nVars; i++) {
			if(x[i][0]<0.0 || x[i][0]>n) {
				throw new RuntimeException("Invalid x supplied. Only values between 0 and n allowed.");
			}
		}
			
	}
	
}
