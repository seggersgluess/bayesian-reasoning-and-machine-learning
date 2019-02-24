
import java.util.*;

public class NormalDistribution {

	static double [][] my;
	static double [][] sigma;
	
	static double [][] L = null;
	
	@SuppressWarnings("static-access")
	NormalDistribution(double [][] my, double [][] sigma){
		
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
