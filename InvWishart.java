package Distributions;

import java.util.ArrayList;

import Mathematics.MatrixOperations;

public class InvWishart {

	double [][] S;
	int df = 0;
	int n = 0;
	
	public InvWishart(double [][] scale_matrix, int df) {
		
		if(scale_matrix.length != scale_matrix[0].length) {
			throw new RuntimeException("Supplied scale matrix is not a square matrix.");
		}
		
		if(df<=scale_matrix.length) {
			throw new RuntimeException("Supplied degrees of freedom have to be greater than " + scale_matrix.length);
		}
		
		this.S = scale_matrix;
		this.n = scale_matrix.length;
		this.df = df;		
	}
	
	
	public double [][] sample() {
		Wishart wDist = new Wishart(S,df);
		double [][] W = wDist.sample();
		W = MatrixOperations.inverse(W);
		return W;
	}
	
	
	public ArrayList<double [][]> sample(int n_samples) {
		ArrayList<double [][]> samples = new ArrayList<double [][]>(n_samples);
		for(int i=0; i<n_samples; i++) {
			samples.add(sample());
		}
		return samples;
	}
	
	
}
