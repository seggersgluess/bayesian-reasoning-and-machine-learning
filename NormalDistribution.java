package Distributions;

import java.util.*;

import Mathematics.Cholesky;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class NormalDistribution {

	double [][] my;
	double [][] sigma;
	
	double [][] L = null;
	
	public NormalDistribution(double [][] my, double [][] sigma) {		
		this.my   = my;
		this.sigma = sigma;	
	}
	
	
	public NormalDistribution(double my, double sigma){
		
		double [][] myArray = new double [1][1];
		double [][] sigmaArray = new double [1][1];
		
		myArray[0][0] = my;
		sigmaArray[0][0] = sigma;
		
		this.my   = myArray;
		this.sigma = sigmaArray;		
	}
	
	
	public double [][] sample(){
		
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
		
	
	public double [][] sample(int n){
		
		int n_vars = my.length;
		
		double [][] X = new double [n_vars][n];
		
		if(n_vars == 1) {
			for(int i=0; i<n; i++){			
				if(i==0){				
					X = sample();				
				}else{
					X = MatrixOperations.rbind(X, sample());				
				}
				
			}
		}
		
		if(n_vars > 1) {
			for(int i=0; i<n; i++){			
				if(i==0){				
					X = sample();				
				}else{
					X = MatrixOperations.cbind(X, sample());				
				}			
			}
		}
		
		L = null;
		
		return X;		
	}
	
	
	public double get_univariateNormalPDF(double x){	
		return 1.0/(Math.sqrt(2.0*Math.PI*Math.pow(sigma[0][0],2.0)))*Math.exp(-Math.pow(x-my[0][0], 2.0)/(2.0*Math.pow(sigma[0][0], 2.0)));			
	}
	
	
	public double [][] get_univariateNormalPDF(double [][] x){
		
		int n=x.length;
		
		double [][] pdf = new double [n][1];
		
		for(int i=0; i<n; i++) {
			pdf[i][0] = get_univariateNormalPDF(x[i][0]);
		}
		
		return pdf;
	}
	
	
	//Zelen & Severo (1964) procedure for numerical approximation of the normal CDF
	public double get_univariateNormalCDF(double x){
		
		double [][] orgSigma = new double [1][1];
		double [][] orgMy = new double [1][1];
		
		orgSigma[0][0] = sigma[0][0];
		orgMy[0][0] = my[0][0];
		
		x = (x-my[0][0])/sigma[0][0];
		
		sigma[0][0] = 1.0;
		my[0][0] = 0.0;
		
		double b_0 = 0.2316419;
		double b_1 = 0.319381530;
		double b_2 = -0.356563782;
		double b_3 =  1.781477937;
		double b_4 = -1.821255978;
		double b_5 = 1.330274429;
		
		double t = 1.0/(1.0+b_0*x);
		double pdf = get_univariateNormalPDF(x);
		
		double cdf = (b_1*t+b_2*Math.pow(t, 2.0)+b_3*Math.pow(t, 3.0)+b_4*Math.pow(t, 4.0)+b_5*Math.pow(t, 5.0));
		cdf = pdf*cdf;
		cdf = 1.0-cdf;
		
		sigma[0][0] = orgSigma[0][0];
		my[0][0] = orgMy[0][0];
		
		return cdf;		
	}
	
	
	public double get_univariateNormalQuantile(double p){
	
		double q = p;
		if(p<0.5){
			q=1.0-p;
		}
		
		double z = -0.4115*((1.0-q)/q + Math.log((1.0-q)/q)-1.0);		
		
		z = z*sigma[0][0]+my[0][0];
		
		if(p<0.5){
			z=-z;
		}
		
		return z;		
	}	
		
	
	public double get_multivariateNormalPDF(double [][] x){
		
		int n = x.length;
		
		double det = MatrixOperations.determinant(sigma);
		double [][] sigma_inv = MatrixOperations.inverse(sigma);
		
		double [][] term = new double [n][1];
		double [][] termTrans = new double [1][n];
		
		for(int i=0; i<n; i++){			
			term[i][0] = x[i][0]-my[i][0];
			termTrans[0][i] = term[i][0];			
		}
		
		double density = MatrixOperations.multiplication(MatrixOperations.multiplication(termTrans, sigma_inv),term)[0][0];		
		density = Math.pow(2.0*Math.PI, -1.0*(my.length/2.0))*1.0/Math.sqrt(det)*Math.exp(-1.0/2.0*density);
		
		return density;	
	}
	
	
	public double [][] get_2d_confidence_ellipse(int nRandomNumbers, double [][] mean_vec, double [][] cov_matrix, double confidence) {
		
		HashMap<String, double[][]> eigenDec = MatrixOperations.get_eigen_dec_4_symmetric_matrix(cov_matrix);
		double [][] eigenvalues = eigenDec.get("eigenvalues");
		double [][] eigenvectors = eigenDec.get("eigenvectors");
		
		int idx = 0;
		if(eigenvalues[0][0]<eigenvalues[1][1]) {
			idx = 1;
		}
		
		if(confidence != 0.9 && confidence != 0.95) {
			throw new RuntimeException("Invalid confidence levels. Only 0.9 or 0.95 allowed.");
		}
		
		//confidence 0.9
		double conf = 4.605;
		if(confidence == 0.95) {
			//confidence 0.95
			conf = 5.991;
		}
		
		double majorLength = Math.sqrt(conf*eigenvalues[1][1]);
		double minorLength = Math.sqrt(conf*eigenvalues[0][0]);
		double angle       = Math.atan(eigenvectors[1][idx]/eigenvectors[0][idx]);
		
		int n = (int) (nRandomNumbers/2.0);
		double steps = majorLength/n;
		double q = minorLength/majorLength;
		
		double [][] ellipse = new double [4*n][2];
				
		for(int i=0;i<n; i++) {
			ellipse[i][0]     = -steps*(n-i);
			ellipse[i][1]     = -q*Math.sqrt((Math.pow(majorLength,2.0)-Math.pow(ellipse[i][0], 2.0)));
			ellipse[n+i][0]   = +steps*(i+1);
			ellipse[n+i][1]   = -q*Math.sqrt((Math.pow(majorLength,2.0)-Math.pow(ellipse[n+i][0], 2.0)));
			ellipse[2*n+i][0] = +steps*(n-i-1);
			ellipse[2*n+i][1] = +q*Math.sqrt((Math.pow(majorLength,2.0)-Math.pow(ellipse[2*n+i][0], 2.0)));
			ellipse[3*n+i][0] = -steps*(i+1);
			ellipse[3*n+i][1] = +q*Math.sqrt((Math.pow(majorLength,2.0)-Math.pow(ellipse[3*n+i][0], 2.0)));	
		}
		

		if(angle != 0.0) {
			double [][] rotation_matrix = GeneralMath.rotation_matrix(angle);
			ellipse = MatrixOperations.multiplication(ellipse, MatrixOperations.transpose(rotation_matrix));
		}
		
		n = ellipse.length;
		for(int i=0; i<n; i++) {
			ellipse[i][0] += mean_vec[0][0];
			ellipse[i][1] += mean_vec[1][0];
		}
		
		return ellipse;
	}
	
	
    // test client
	public static void main(String[] args) {
    	
    	//my = new double [2][1];
    	//sigma = MatrixOperations.scalar_multiplication(2.0, MatrixOperations.identity(2));
    	
    	//MatrixOperations.print_matrix(GeneralMath.mean_vec(MatrixOperations.transpose(sample(1000))));
    	
    	double [][] mean = new double [1][1];
    	double [][] sd = new double [1][1];
    	
    	mean[0][0] = 0.0;
    	sd[0][0] = 1.0;
    	
    	NormalDistribution nd = new NormalDistribution(mean, sd);
    	System.out.println(nd.get_univariateNormalCDF(-0.3));
    	System.out.println(nd.get_univariateNormalQuantile(0.05));
    	
    }
		
}
