
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
	
	
	// calculates mean of a series x
	public static double mean(double [] x){
		
		int n_elements = x.length;
		double mean = 0.0;
		
		for(int i=0; i<n_elements; i++){
			
			mean = mean + x[i];
			
		}
		
		mean = mean/(double) n_elements;
		
		return mean;
		
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
    	
    	double [] a= {2.0, 3.0, 8.0};
    	
    	System.out.println(mean(a));
    	
    }
	
	
	
	
}
