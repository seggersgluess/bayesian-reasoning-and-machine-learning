
public class InvGammaDistribution {

	static double shape;
	static double scale;
	static double rate;
	
	@SuppressWarnings("static-access")
	InvGammaDistribution(double shape, double scale){
		
		if(shape <= 0.0){
			
			throw new RuntimeException("Unallowed shape parameter "+ shape +" less than 0.");
			
		}
		
		if(scale <= 0.0){
			
			throw new RuntimeException("Unallowed shape parameter "+ shape +" less than 0.");
			
		}
		
		this.scale = scale;		
		this.shape = shape;
		
		rate = 1.0/scale;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [][] sample(int n){
		
		GammaDistribution gammaDist = new GammaDistribution(shape, scale);
		
		double[][] x     = new double [n][1];
		
		x = gammaDist.sample(n);
		
		for(int i=0; i<n; i++){
			
			x[i][0] = 1.0/x[i][0];
			
		}
		
		return x;
		
	}
	
	
	public static double pdf(double x){
			
		double pdf = Math.pow(scale, shape)/(GammaFunction.gamma(shape))*Math.pow(x, -1.0*(shape+1.0))*Math.exp(-1.0*scale/x);
		
		return pdf;
		
	}
	
	
	// returns pdfs for n x 1 vector of inv gamma dist. random variables x
	public static double [][] pdf(double [][] x){
		
		int n = x.length;
		double [][] pdf = new double [n][1];
		
		for(int i=0; i<n; i++){
			
			pdf[i][0] = pdf(x[i][0]);
			
		}
		
		return pdf;
		
	}
	

	public static double mean(){
				
		return scale/(shape-1.0);
		
	}
	
	
	public static double variance(){
		
		if(shape <= 2.0){
			
			throw new RuntimeException("Secound moment can only be calculated for shape larger than 3.0.");
			
		}
		
		return Math.pow(scale, 2.0)/(Math.pow(shape-1.0, 2.0)*(shape-2.0));
		
	}
	
	
	public static double sd(){
		
		return Math.sqrt(variance());
		
	}
	
	
	public static double skewness(){
		
		if(shape <= 3.0){
			
			throw new RuntimeException("Skewness can only be calculated for shape larger than 3.0.");
			
		}
		
		return 4.0*Math.sqrt(shape-2.0)/(shape-2.0);
		
	}
	
	
	public static double excess_kurtosis(){
		
		if(shape <= 4.0){
			
			throw new RuntimeException("Excess kurtosis can only be calculated for shape larger than 3.0.");
			
		}
		
		return (30.0*shape-66.0)/((shape-3.0)*(shape-4.0));
		
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	shape = 2.0;
    	scale = 2.0;
    	
    	//double [][] a = {{2.0},{8.12}};
    	
    	double [][] sample = sample(100);
    	
    	MatrixOperations.print_matrix(sample);
    	
    	System.out.println("");
    	
    	System.out.println(GeneralMath.mean(sample));
    	
    	System.out.println("");
    	
    	System.out.println(mean());
    	
    }
	
}
