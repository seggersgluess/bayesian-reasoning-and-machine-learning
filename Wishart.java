package Distributions;

import java.util.ArrayList;

import Mathematics.GammaFunction;
import Mathematics.LU;
import Mathematics.MatrixOperations;

public class Wishart {

	double [][] S;
	int df = 0;
	int n = 0;
	
	public Wishart(double [][] scale_matrix, int df) {
		
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
		
	
	public double get_WishartPDF(double [][] X) {
		
		double det_X = MatrixOperations.determinant(X);
		double det_S = MatrixOperations.determinant(S);
		double [][] S_inv = MatrixOperations.inverse(S);
		
		double mult_gamma = GammaFunction.gamma(df,n);
		
		double pdf = Math.pow(det_X, (double)((df-n-1.0)/2.0));
		pdf *= Math.exp(-1.0/2.0*MatrixOperations.trace(MatrixOperations.multiplication(X,S_inv)));
		pdf /= Math.pow(2.0, df*n/2.0)*Math.pow(det_S, df/2.0)*mult_gamma;
		
		return pdf;
	}
	
	
	public double [][] sample() {
		
		NormalDistribution normDist = new NormalDistribution(0.0, 1.0);
		double [][] W = new double [n][n];
		int n_mod = n-1;
		for(int i=0; i<n_mod; i++) {
			int idx = i+1;
			ChiSquared chi2Dist = new ChiSquared(df-idx+1);
			W[i][i] = Math.sqrt(chi2Dist.sample());		
			for(int j=idx; j<n; j++) {
				W[i][j] = normDist.sample()[0][0];
			}
		}
		ChiSquared chi2Dist = new ChiSquared(df-n+1);
		W[n-1][n-1] = Math.sqrt(chi2Dist.sample());
		W = MatrixOperations.multiplication(W, MatrixOperations.transpose(W));
		
		LU lu_dec = new LU(S);
		double [][] L = lu_dec.get_lower_triangular();
		W = MatrixOperations.multiplication(MatrixOperations.multiplication(L, W),MatrixOperations.transpose(L));
		
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
