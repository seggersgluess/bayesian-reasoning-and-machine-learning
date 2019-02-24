
import java.util.*;

public class GammaDistribution {
	
	static double shape;
	static double scale;
	static double rate;
	
	@SuppressWarnings("static-access")
	GammaDistribution(double shape, double scale){
		
		this.scale = scale;		
		this.shape = shape;
		
		rate = 1.0/scale;
		
	}
	
	
	public static double pdf(double x){
		
		double pdf = Math.pow(rate, shape)/GammaFunction.gamma(shape)*Math.pow(x, 1.0-shape)*Math.exp(-x*rate);
		
		return pdf;
		
	}
	
	
	// returns pdfs for n x 1 vector of gamma dist. random variables x
	public static double [][] pdf(double [][] x){
		
		int n = x.length;
		double [][] pdf = new double [n][1];
		
		for(int i=0; i<n; i++){
			
			pdf[i][0] = Math.pow(rate, shape)/GammaFunction.gamma(shape)*Math.pow(x[i][0], 1.0-shape)*Math.exp(-x[i][0]*rate);
			
		}
		
		return pdf;
		
	}
	
	
	// implements the Marsaglia and Tsang’s Method for generating gamma distributed random variables
	public static double [][] sample(int n){
		
		double[][] x     = new double [n][1];		
		double z         = 0;
		double u         = 0;
		double v         = 0;
		double threshold = 0;
		double modified_shape  = 0;
		
		if(shape >= 1.0){
			
			modified_shape = shape;
			
		}else{
			
			modified_shape = shape + 1.0;
		}
		
		double d = modified_shape - 1.0/3.0;
		double c = 1.0/Math.sqrt(9.0*d);
		
		Random r = new Random();
		
		for(int i=0; i<n; i++){
			
			boolean flag = false;
			
			while (flag == false) {
					
				z = r.nextGaussian();
				
				u = Math.random();
								
				v = Math.pow(1+c*z, 3.0);
							
				threshold = 0.5*Math.pow(z,2.0)+d-d*v + d*Math.log(v);		
					
				if((z <=-1.0/c && Math.log(u)>=threshold) || Double.isNaN(threshold) == true){
						
					flag = false;
						
				}else{
						
					flag = true;
						
				}
						
			}
							
			x[i][0] = d*v/rate;
							
			if(shape < 1.0){
			
				x[i][0] = x[i][0]*Math.pow(Math.random(),1.0/shape);
				
			}
			
		}

		return x;
		
	}
		
	
	public static double mean(){
		
		return shape/rate;
		
	}
	
	
	public static double variance(){
		
		return shape/Math.pow(rate, 2.0);
		
	}
	
	
	public static double sd(){
		
		return Math.sqrt(shape/Math.pow(rate, 2.0));
		
	}
	
	
	public static double skewness(){
		
		return 2.0/Math.sqrt(shape);
		
	}
	
	
	public static double excess_kurtosis(){
		
		return 6.0/shape;
		
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	shape = 0.5;
    	scale = 2.0;
    	rate  = 1.0/scale;
    	
    	//double [][] a = {{2.52220716}, {1.52526767}};
    	
    	System.out.println(mean());
    	
    	//MatrixOperations.print_matrix(pdf(a));
    	System.out.println(GeneralMath.mean(sample(1000)));
    	
    	
    }
	
}
