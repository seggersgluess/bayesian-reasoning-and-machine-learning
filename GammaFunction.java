package Mathematics;
public class GammaFunction {

	// implements simple routine by Lanczos
	public static double gamma(double z){
	
		double gamma = 0.0;
		
	    if(z < 0.5){	    	
	    	gamma = Math.PI/(Math.sin(Math.PI*z)*gamma(1.0-z));	  	
	    }else{	    	
			z = z-1.0;			
			double [][] c = {	
								{
									676.5203681218851,
			                     	-1259.1392167224028,
			                     	771.32342877765313,
			                     	-176.61502916214059,
			                     	12.507343278686905,
			                     	-0.13857109526572012,
			                     	9.9843695780195716e-06,
			                     	1.5056327351493116e-07		                     	
								}
							};
			
			c = MatrixOperations.transpose(c);
			
			int g = c.length;
			
			double [][] z_vec = new double [1][g];
			
			for(int i=0; i<g; i++){											
				z_vec[0][i] = 1.0/(z+i+1.0);				 	
			}
				
			double a_g = MatrixOperations.multiplication(z_vec, c)[0][0]+1.0;				
			gamma      = Math.sqrt(2*Math.PI)*Math.pow((z+g-0.5), z+0.5)*Math.exp(-1.0*(z+g-0.5))*a_g;    	
	    }

		return gamma;				
	}
	
	//Multivariate Gamma function 
	public static double gamma(double z, int p) {
		
		if(p<1) {
			throw new RuntimeException("Dimension of Gamma function has to be larger than 0");
		}
		
		double c = Math.pow(Math.PI,(double)(p*(p-1)/4.0));
		double gamma = 1.0;
		for(int i=0; i<p; i++) {
			double idx = i+1;
			gamma *= gamma(z+(1.0-idx)/2.0);
		}
		gamma *= c;
		
		return gamma;
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	System.out.println(gamma(0.11));
    	
    }
	
}
