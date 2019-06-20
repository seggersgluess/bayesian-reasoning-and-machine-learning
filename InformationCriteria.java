package Mathematics;

public class InformationCriteria {

    /**
     * Method for calculating Akaike information criterion
     * @param maxLikelihood maximal value of model´s likelihood
     * @param numberOfParameters number of estimated model parameters
     * @return Akaike information criterion
     */
	public static double aic(double maxLikelihood, int numberOfParameters){
		
		if(maxLikelihood <= 0.0){
			throw new RuntimeException("Invalid likelihood supplied.");
		}
		
		if(numberOfParameters <= 0){
			throw new RuntimeException("Invalid number of model parameters supplied.");
		}
		
		double aic = numberOfParameters*Math.log(maxLikelihood);
		
		return aic;
		
	}
	
	
    /**
     * Method for calculating Bayesian information criterion
     * @param maxLikelihood maximal value of model´s likelihood
     * @param numberOfParameters number of estimated model parameters
     * @param numberOfObservations number of observations (sample length)    
     * @return Bayesian information criterion
     */
	public static double bic(double maxLikelihood, int numberOfParameters, int numberOfObservations){
		
		if(maxLikelihood <= 0.0){
			throw new RuntimeException("Invalid likelihood supplied.");
		}
		
		if(numberOfParameters<= 0.0){
			throw new RuntimeException("Invalid likelihood supplied.");
		}
		
		if(numberOfObservations <= 0){
			throw new RuntimeException("Invalid number of observations supplied.");
		}
		
		double bic = Math.log(numberOfObservations)*numberOfParameters-2.0*Math.log(maxLikelihood);
		
		return bic;
		
	}
	
	
    /**
     * Method for calculating Hannan-Quinn information criterion
     * @param maxLikelihood maximal value of model´s likelihood
     * @param numberOfParameters number of estimated model parameters
     * @param numberOfObservations number of observations (sample length)    
     * @return Hannan-Quinn information criterion
     */
	public static double hq(double maxLikelihood, int numberOfParameters, int numberOfObservations){
		
		if(maxLikelihood <= 0.0){
			throw new RuntimeException("Invalid likelihood supplied.");
		}
		
		if(numberOfParameters<= 0.0){
			throw new RuntimeException("Invalid likelihood supplied.");
		}
		
		if(numberOfObservations <= 0){
			throw new RuntimeException("Invalid number of observations supplied.");
		}
		
		double bic = Math.log(Math.log(numberOfObservations))*2.0*numberOfParameters-2.0*Math.log(maxLikelihood);
		
		return bic;
		
	}
	
	
}
