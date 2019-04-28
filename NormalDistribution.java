package Distributions;

import java.util.*;

import Mathematics.Cholesky;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class NormalDistribution {

	static double [][] my;
	static double [][] sigma;
	
	static double [][] L = null;
	
	@SuppressWarnings("static-access")	
	public NormalDistribution(double [][] my, double [][] sigma){
		
		this.my   = my;
		this.sigma = sigma;
		
	}
	
	
	public static double [][] sample(){
		
		int n = my.length;
		
		if(L == null){
			
			L = Cholesky.decompose(sigma);
			
		}
			
		double [][] x = new double [n][1];
		
		Random r = new Random();
		
		for(int i=0; i<n; i++){
			
			x[i][0] = r.nextGaussian();
			
		}
		
		x = MatrixOperations.add(my, MatrixOperations.multiplication(L, x));
		
		return x;
		
	}
	
	
	public static double get_univariateNormalPDF(double x){
		
		return 1.0/(Math.sqrt(2.0*Math.PI*Math.pow(sigma[0][0],2.0)))*Math.exp(-Math.pow(x-my[0][0], 2.0)/(2.0*Math.pow(sigma[0][0], 2.0)));
			
	}
	
	
	public double get_multivariateNormalPDF(double [][] x){
		
		//ToDo´s: Determine determinant of Sigma!
		int n = x.length;
		
		double det = 1.0;
		double [][] sigma_inv = MatrixOperations.inverse(sigma);
		
		double [][] term = new double [n][1];
		double [][] termTrans = new double [1][n];
		
		for(int i=0; i<n; i++){			
			term[i][0] = x[i][0]-my[i][0];
			termTrans[0][i] = term[i][0];			
		}
		
		double density = MatrixOperations.multiplication(MatrixOperations.multiplication(termTrans, sigma_inv),term)[0][0];		
		density = Math.pow(2.0*Math.PI, my.length/2.0)*det*Math.exp(-1.0/2.0*density);
		
		return density;
		
	}
	
	
	public static double [][] sample(int n){
		
		int n_vars = my.length;
		
		double [][] X = new double [n_vars][n];
		
		for(int i=0; i<n; i++){
			
			if(i==0){
				
				X = sample();
				
			}else{
				
				X = MatrixOperations.cbind(X, sample());
				
			}
			
		}
		
		L = null;
		
		return X;
		
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	my = new double [2][1];
    	sigma = MatrixOperations.scalar_multiplication(2.0, MatrixOperations.identity(2));
    	
    	MatrixOperations.print_matrix(GeneralMath.mean_vec(MatrixOperations.transpose(sample(1000))));
    	
    }
	
	
}
