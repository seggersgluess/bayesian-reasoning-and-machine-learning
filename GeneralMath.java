package Mathematics;
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
	
	
	// calculates factorial with x: x*(x-1)*(x-2)*...
	public static double factorial(double x){
		
		double fact = 1.0;
		
		for(int i=0; i<x; i++){			
			fact = fact*(i+1);			
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
			mean = mean + x[i][0];			
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
				mean[0][i] = mean[0][i] + X[j][i];				
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
				
	
    // test client
    public static void main(String[] args) {
    	
    	double [][] x = {{1.0,2.0,3.0,4.0},{1.0,2.0,3.0,4.0}};
    	
    	double [] y = {1.0,2.0,3.0,4.0};
    	
    	MatrixOperations.print_matrix(sd_vec(MatrixOperations.transpose(x)));
    	
    	System.out.println(variance(y));
    	
    }
	
}
