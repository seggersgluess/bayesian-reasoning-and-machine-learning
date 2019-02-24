import java.util.Random;

public class StudentDistribution {

	static double [][] my;
	static double [][] sigma;
	static double df;
	
	static double [][] L = null;
	
	@SuppressWarnings("static-access")
	StudentDistribution(double [][] my, double [][] sigma, double df){
		
		this.my    = my;
		this.sigma = sigma;
		this.df    = df;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [][] sample(){
		
		int n = my.length;
		
		if(L == null){
			
			L = Cholesky.decompose(sigma);
			
		}
			
		double [][] x = new double [n][1];
		
		Random r = new Random();
		GammaDistribution g = new GammaDistribution(df/2.0, 2.0*df);
		
		for(int i=0; i<n; i++){
			
			x[i][0] = r.nextGaussian();
			
		}
		
		x = MatrixOperations.multiplication(L, x);
		
		double rg = 1.0/Math.sqrt(g.sample(1)[0][0]);

		x = MatrixOperations.scalar_multiplication(rg, x);				
		x = MatrixOperations.add(my, x);
		
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
    	my[0][0] = 8.0;
    	my[1][0] = 12.0;
    	
    	sigma = MatrixOperations.scalar_multiplication(2.0, MatrixOperations.identity(2));
    	df = 5.0;
    	
    	System.out.println("Mean:");
    	
    	MatrixOperations.print_matrix(sample(1));
    	
    }
	
}
