package Mathematics;
import java.util.Arrays;
import java.util.List;


public class GeneralMath {

	public static double power(double x, int b){
		
	    double res =1;
	    for (int i = 0; i < b; i++) {
	        res *= x;
	    }
	    return res;
		
	}
	
	
	// calculates sum of a series x
	public static double sum(double [] x){
		
		double sum = 0;
		
		for(int i=0; i<x.length; i++){			
			sum = sum + x[i];			
		}
		
		return sum;
		
	}
	
	
	// calculates sum of a series x
	public static int sum(int [] x){
		
		int sum = 0;
		
		for(int i=0; i<x.length; i++){			
			sum = sum + x[i];			
		}
		
		return sum;
		
	}
	
	
	// calculates sum of a series x
	public static double sum(double [][] x){
		
		double sum = 0;
		
		for(int i=0; i<x.length; i++){			
			sum = sum + x[i][0];			
		}
		
		return sum;
		
	}
	
	
	// calculates sum of a series x
	public static int sum(List<Integer> x){
		
		int sum = 0;
		
		for(int i=0; i<x.size(); i++){			
			sum = sum + x.get(i);			
		}
		
		return sum;
		
	}
	
	
	// calculates sum of a series x
	public static double sumDblList(List<Double> x){
		
		double sum = 0.0;
		
		for(int i=0; i<x.size(); i++){			
			sum = sum + x.get(i);			
		}
		
		return sum;
		
	}
	
	
	public static double [][] cov(double [][] X){
		
		int n_obs = X.length;
		int n_variables = X[0].length;
		
		double [][] mean = new double [n_variables][1];
	    double [][] cov  = new double [n_variables][n_variables];
		
		for(int i=0; i<n_variables; i++){
			for(int j=0; j<n_obs; j++){
				mean[i][0] += X[j][i];
			}
			mean[i][0] = mean[i][0]/n_obs;
		}
		
		for(int i=0; i<n_variables; i++){
			for(int j=0; j<(i+1); j++){
				
				if(i==j){
					for(int k=0; k<n_obs; k++){
						cov[i][j] += Math.pow(X[k][i]-mean[i][0],2.0);
					}
					cov[i][j] = cov[i][j]/(n_obs-1.0);
				}
				
				if(i!=j){
					for(int k=0; k<n_obs; k++){
						cov[i][j] += (X[k][i]-mean[i][0])*(X[k][j]-mean[j][0]);
					}
					cov[i][j] = cov[i][j]/(n_obs-1);
					cov[j][i] = cov[i][j];
				}
				
			}
		}
		
		return cov;
		
	}
	
	
	// calculates factorial with x: x*(x-1)*(x-2)*...
	public static double factorial(double x){
		
		double fact = 1.0;
		
		for(int i=0; i<x; i++){			
			fact *= (i+1);			
		}
		
		return fact;
		
	}
	
	
	// calculates mean of a series x
	public static double mean(double [] x){
		
		int n_elements = x.length;
		double mean = 0.0;
		
		for(int i=0; i<n_elements; i++){			
			mean = mean + x[i];			
		}
		
		mean = mean/(double)n_elements;
		
		return mean;
		
	}
	
	
	// calculates mean of a series x
	public static double mean(List<Double> x){
		
		int n_elements = x.size();
		double mean = 0.0;
		
		for(int i=0; i<n_elements; i++){			
			mean = mean + x.get(i);			
		}
		
		mean = mean/(double)n_elements;
		
		return mean;
		
	}
	
	
	// calculates mean of a series x
	public static double mean(double [][] x){
		
		int n_elements = x.length;
		double mean = 0.0;
		
		for(int i=0; i<n_elements; i++){			
			mean += x[i][0];			
		}
		
		mean = mean/(double)n_elements;
		
		return mean;
		
	}
	
	
	// calculates mean of a series x
	public static double [][] mean_vec(double [][] X){
		
		int n_obs  = X.length;
		int n_vars = X[0].length;
				
		double [][] mean = new double [1][n_vars];
		
		for(int i=0; i<n_vars; i++){			
			for(int j=0; j<n_obs; j++){				
				mean[0][i] += X[j][i];				
			}
			
			mean[0][i] = mean[0][i]/(double)n_obs;
			
		}
		
		return mean;
		
	}
	
	
	// returns variance of a vector x
	public static double variance(double [] x){
		
		int n_elements = x.length;
		double mean = mean(x);
		double variance = 0;
		
		for(int i=0; i<n_elements; i++){			
			variance = variance + Math.pow((x[i]-mean),2.0);			
		}
		
		variance = variance/((double) n_elements - 1.0);
		
		return variance;
		
	}
	
	
	// returns variance of a vector x
	public static double variance(List<Double> x){
		
		int n_elements = x.size();
		double mean = mean(x);
		double variance = 0;
		
		for(int i=0; i<n_elements; i++){			
			variance = variance + Math.pow((x.get(i)-mean),2.0);			
		}
		
		variance = variance/((double) n_elements - 1.0);
		
		return variance;
		
	}
	
	
	// returns variance of a vector x
	public static double variance(double [][] x){
		
		int n_elements  = x.length;
		double mean     = mean(x);
		double variance = 0;
		
		for(int i=0; i<n_elements; i++){			
			variance = variance + Math.pow((x[i][0]-mean),2.0);			
		}
		
		variance = variance/((double) n_elements - 1.0);
		
		return variance;
		
	}
	
	
	// calculates mean of a series x
	public static double [][] variance_vec(double [][] X){
		
		int n_obs  = X.length;
		int n_vars = X[0].length;
		
		double [][] variance = new double [1][n_vars];
		
		for(int i=0; i<n_vars; i++){
			
			for(int j=0; j<n_obs; j++){				
				variance[0][i] = variance(MatrixOperations.get_column_from_matrix(X,i));				
			}
					
		}
		
		return variance;
		
	}
	
	
	// returns standard deviation of a vector x
	public static double sd(double [] x){
	
		double sd = Math.sqrt(variance(x));
		
		return sd;
		
	}
	
	
	// returns standard deviation of a vector x
	public static double sd(double [][] x){
	
		double sd = Math.sqrt(variance(x));
		
		return sd;
		
	}
		
	
	// returns standard deviation of a vector x
	public static double [][] sd_vec(double [][] X){
	
		double [][] sd = sqrt(MatrixOperations.transpose(variance_vec(X)));
		
		sd = MatrixOperations.transpose(sd);
		
		return sd;
		
	}
	
	
	// returns vector of square roots
	public static double [][] sqrt(double [][] x){
		
		int n_rows = x.length;
		
		double [][] sqrt_vec = new double [n_rows][1];
		
		for(int i=0; i<n_rows; i++){			
			sqrt_vec[i][0] = Math.sqrt(x[i][0]);			
		}
		
		return sqrt_vec;
		
	}
	
	
	public static double sincf(double x) {
		
		if(x==0.0) {
			return 1.0;
		}
		
		double p = Math.PI*x;
		return Math.sin(p)/p;
	}
	
	
	public static double quantile(double [] sample, double p){
		
		if(p>1.0 || p<0.0){		
			throw new RuntimeException("Invalid probability p supplied.");			
		}
		
		int n = sample.length;
		double quantile;
		double idx;
		
	    Arrays.sort(sample, 0, n);
			
	    idx = p*n;
	    if(Math.ceil(idx)-idx > 0.5){	        		
	    	idx = idx-1;        		
	    }
	    	
	    if(idx >= 0.0){        		
	        idx = Math.ceil(idx)-1.0; 		
	    }else{	
	        idx = 0.0;	        		
	    }
	    	
	    quantile = sample[(int)idx];
	        	
		return quantile;
		
	}
	
	
	public static double [][] rotation_matrix(double angle) {
		
		double [][] rotation_matrix = new double [2][2];
		
		rotation_matrix[0][0] = Math.cos(angle);
		rotation_matrix[1][1] = rotation_matrix[0][0];
		rotation_matrix[0][1] = -Math.sin(angle);
		rotation_matrix[1][0] = -rotation_matrix[0][1];
		
		return rotation_matrix;
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	double [][] x = {{1.0,2.0,3.0,4.0},{1.0,2.0,3.0,4.0}};
    	
    	double [] y = {1.0,2.0,3.0,4.0};
    	
    	MatrixOperations.print_matrix(sd_vec(MatrixOperations.transpose(x)));
    	
    	System.out.println(variance(y));
    	

    	double [][]A={{1.0,  10},
    				  {2.0, 203},
    				  {35.0,   3},
    				  {20.0,  40},
    				  {300.3,  39}};
    	
    	MatrixOperations.print_matrix(cov(A));
    	
    }
	
}
