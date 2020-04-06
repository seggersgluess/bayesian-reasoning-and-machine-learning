package Distributions;

public class ChiSquared {

	int df = 0;
	
	public ChiSquared(int df) {
		if(df<=0) {
			throw new RuntimeException("Only degrees of freedom > 0 allowed for chi-squared distribution.");
		}
	}
	
	
	public double sample() {
		
		double r = 0.0;
		NormalDistribution normDist = new NormalDistribution(0.0,1.0);
		for(int i=0; i<df; i++) {
			r+=Math.pow(normDist.sample()[0][0], 2.0);
		}
		return r;
	}
	
	
	public double [][] sample(int n_samples) {
		
		double [][] r = new double [n_samples][1];
		NormalDistribution normDist = new NormalDistribution(0.0,1.0);
		for(int s=0; s<n_samples; s++) {
			for(int i=0; i<df; i++) {
				r[i][0] +=Math.pow(normDist.sample()[0][0], 2.0);
			}
		}

		return r;
	}
	
	
}
