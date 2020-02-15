package NaiveBayesClassifier;

import java.util.HashMap;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;

public class NBC {

	static double [][] explained_variable;
	static double [][] explaining_variables;	
	static String name_of_explained_variable;
	static String [] names_of_explaining_variables;
	
	static int n_observations;
	static int n_explaining_variables;
	
	static double [] classes;
	static int n_classes;
	
	static double [][] mod_explained_variable;
	
	static String distribution;
	
	static Object fittedNBC;
	
	double [][] pred_explained_variable;
	double [][] pred_probability;
	
	double accuracyRate;
	
	//Regularization for CategoricalNB & MultinomialNB
	double alpha = 1.0;
	
	public void fitNBC(String usedDistribution) {
		
		setDistribution(usedDistribution);
				
		convertExplainedVariable();
			
		if(distribution.contentEquals("GAUSSIAN")) {
			fittedNBC = new GaussianNB();
			((GaussianNB) fittedNBC).fit();
		}
		
		if(distribution.contentEquals("CATEGORICAL")) {
			fittedNBC = new CategoricalNB();
			((CategoricalNB) fittedNBC).fit();
		}
		
		if(distribution.contentEquals("MULTINOMIAL")) {
			fittedNBC = new MultinomialNB();
			((MultinomialNB) fittedNBC).fit();
		}
		
		if(distribution.contentEquals("BERNOULLI")) {			
			fittedNBC = new BernoulliNB();
			((BernoulliNB) fittedNBC).fit();
		}
		
	}
	
	
	public void fitNBC(String usedDistribution, double [][] y, double [][] X) {
	
		if(y.length != X.length) {
			throw new RuntimeException("Unequal length of y and X.");
		}
		
		n_observations = y.length;
		n_explaining_variables = X[0].length;
		explained_variable   = y;
		explaining_variables = X;
		
		fitNBC(usedDistribution);
		
	}
	
	
	public double [][] predictNB(double [][] X) {
		
		int n = X.length;
		
		double [][] y_hat = new double [n][1];
		
		if(distribution.contentEquals("GAUSSIAN")) {
			y_hat = ((GaussianNB) fittedNBC).predict(X).get("Prediction");
		}
		
		if(distribution.contentEquals("CATEGORICAL")) {
			y_hat = ((CategoricalNB) fittedNBC).predict(X).get("Prediction");
		}
		
		if(distribution.contentEquals("MULTINOMIAL")) {
			y_hat = ((MultinomialNB) fittedNBC).predict(X).get("Prediction");
		}
		
		if(distribution.contentEquals("BERNOULLI")) {
			y_hat = ((BernoulliNB) fittedNBC).predict(X).get("Prediction");
		}
		
		return y_hat;
	}
	
	
	public double [][] predictNB_probability(double [][] X) {
		
		int n = X.length;
		
		double [][] y_hat = new double [n][1];
		
		if(distribution.contentEquals("GAUSSIAN")) {
			y_hat = ((GaussianNB) fittedNBC).predict(X).get("Probability");
		}
		
		if(distribution.contentEquals("CATEGORICAL")) {
			y_hat = ((CategoricalNB) fittedNBC).predict(X).get("Probability");
		}
		
		if(distribution.contentEquals("MULTINOMIAL")) {
			y_hat = ((MultinomialNB) fittedNBC).predict(X).get("Probability");
		}
		
		if(distribution.contentEquals("BERNOULLI")) {
			y_hat = ((BernoulliNB) fittedNBC).predict(X).get("Probability");
		}
		
		return y_hat;
	}
	
	
	public HashMap<String, double [][]> predictNB_value_and_probability(double [][] X) {

		HashMap<String, double [][]> predRes = new HashMap<String, double [][]>();
		
		if(distribution.contentEquals("GAUSSIAN")) {
			predRes = ((GaussianNB) fittedNBC).predict(X);
		}
		
		if(distribution.contentEquals("CATEGORICAL")) {
			predRes = ((CategoricalNB) fittedNBC).predict(X);
		}
		
		if(distribution.contentEquals("MULTINOMIAL")) {
			predRes = ((MultinomialNB) fittedNBC).predict(X);
		}
		
		if(distribution.contentEquals("BERNOULLI")) {
			predRes = ((BernoulliNB) fittedNBC).predict(X);
		}
		
		return predRes;
	}
	
	
	public void inSamplePredictionNB() {	
		HashMap<String, double [][]> predRes = predictNB_value_and_probability(explaining_variables);
		pred_explained_variable = predRes.get("Prediction");	
		pred_probability = predRes.get("Probability");
		calcAccuracyRates();
	}
	
	
	public double logSumExp(double [] b) {
		
		double logSum = 0.0;
		double max = Utilities.Utilities.getMax(b);
		
		for(int c=0; c<n_classes; c++) {
			logSum += Math.exp(b[c]-max);
		}
		
		logSum = Math.log(logSum);
		logSum += max;
		
		return logSum;
	}
	
	
	public void calcAccuracyRates() {
		
		if(pred_explained_variable == null) {
			System.out.println("No insample prediction done yet. Make insamle prediction for NBC now.");
			inSamplePredictionNB();
		}
				
		double validClassCount = 0.0;
		for(int j=0; j<n_observations; j++) {
			if(pred_explained_variable[j][0] == explained_variable[j][0]) {
				validClassCount++;
			}
		}
		
		accuracyRate = validClassCount/n_observations;
						
	}
	
	
	public double calcAccuracyRates(double [][] predictedValues, double [][] trueValues) {
		
		if(predictedValues.length != trueValues.length) {
			throw new RuntimeException("Unequal length of predicted and true values.");
		}
				
		int n = predictedValues.length;
		double validClassCount = 0.0;
		for(int j=0; j<n; j++) {
			if(predictedValues[j][0] == trueValues[j][0]) {
				validClassCount++;
			}
		}
		
		return validClassCount/n;					
	}
	
	
	public void convertExplainedVariable() {
		
		if(explained_variable == null) {
			throw new RuntimeException("Explained variable not set yet.");
		}
		
		double [] y = MatrixOperations.get_column_from_matrix(explained_variable, 0);
        classes = Utilities.Utilities.get_unique_elements(y);
		
		n_classes = classes.length;

		mod_explained_variable = new double [n_observations][1];
		
		for(int i=0; i<n_observations; i++) {
			for(int c=0; c<n_classes; c++) {
				if(classes[c] == explained_variable[i][0]) {
					mod_explained_variable[i][0] = c;
				}
			}
		}
	
	}
	
	
	public boolean checkDistribution(String usedDistribution) {
		
		boolean valid = true;
		
		String [] validDist = getValidNBCDistribution();
		int [] idxs = Utilities.Utilities.get_idx(validDist, usedDistribution);
		
		if(idxs[0] == -1) {
			System.out.println(usedDistribution + " is not a valid algorithm for AdaBoost.");
			valid = false;
		}
		
		return valid;
	}
	
	
	public void setDistribution(String usedDistribution) {
		
		boolean validDist = checkDistribution(usedDistribution);
		if(validDist == false) {
			throw new RuntimeException(usedDistribution + " is not a valid Distribution for Naive Bayes Classifier.");
		}
		
		distribution = usedDistribution;		
	}
	
	
	public void setAlpha4Regularization(double alpha) {
		if(alpha < 0.0 || alpha >1.0) {
			throw new RuntimeException(alpha + " is not a valid alpha value for regularization only values between 0 and 1 allowed.");
		}
		this.alpha = alpha;
	}
	
	
	public String [] getValidNBCDistribution() {
		
		String [] validDist = {"GAUSSIAN",
							   "CATEGORICAL",
		                       "MULTINOMIAL",
		                       "BERNOULLI"};
		
		return validDist;
	}
	
	
	public double [][] get_explaining_variables() {
		return explaining_variables;
	}
	
	
	public double [][] get_explained_variable() {
		return explained_variable;
	}
	
	
	public double get_accuracyRate() {
		return accuracyRate;
	}
	
	
	public String get_usedNBC_distribution() {
		return distribution;
	}
	
	
	public double [][] get_inSamplePredProbability() {
		return pred_probability;
	}
	
	
	public void readAndSetInputParameters(String dirName, String name_explained_variable) {
		
		InputDataManager inputData = new InputDataManager();			
		try {
			inputData.fileReader(dirName, false, true, true);
		} catch (Exception e) {
			System.out.println("File with input data for naive Bayes cannot be uploaded.");
			e.printStackTrace();
		}
		
		String [] colNames = inputData.colnames;
		int [] idx = Utilities.Utilities.get_idx(colNames, name_explained_variable);
		if(idx[0] == -1) {
			throw new RuntimeException(name_explained_variable + " not found in loaded input data.");
		}
		
		n_explaining_variables = colNames.length-1;
		n_observations = inputData.numberOfRows-1;
		
		name_of_explained_variable = name_explained_variable;
		names_of_explaining_variables = new String [n_explaining_variables];
		
		String [] varNames = new String [n_explaining_variables+1];
		
		varNames[0] = name_explained_variable;
		int nameIdx = 0;
		for(int i=0; i<colNames.length; i++) {
			if(i != idx[0]) {
				varNames[nameIdx+1] = colNames[i];
				names_of_explaining_variables[nameIdx] = colNames[i];
				nameIdx++;
			}
		}
		
		explained_variable = new double [n_observations][1];
		explaining_variables = new double [n_observations][n_explaining_variables];
		
		String [] rowNames = new String[n_observations];
		for(int i=0; i<n_observations; i++){
    		rowNames[i] = Integer.toString(i+1);
    	}
		
		inputData.selectLoadedData(rowNames, varNames);
		
		for(int i=0; i<n_observations; i++) {
			explained_variable[i][0] = Double.parseDouble(inputData.selectedStrFileData[i][0]);
			for(int j=0; j<n_explaining_variables; j++) {
				explaining_variables[i][j] = Double.parseDouble(inputData.selectedStrFileData[i][j+1]);
			}
		}
		
	}
	
	
	public HashMap<String, Object> get_est_NBC_pars() {
		
		HashMap<String, Object> est_pars = new HashMap<String, Object>();
		
		if(distribution.contentEquals("GAUSSIAN")) {
			est_pars = ((GaussianNB) fittedNBC).get_est_pars();
		}
		
		if(distribution.contentEquals("CATEGORICAL")) {
			est_pars = ((CategoricalNB) fittedNBC).get_est_pars();
		}
		
		if(distribution.contentEquals("MULTINOMIAL")) {
			est_pars = ((MultinomialNB) fittedNBC).get_est_pars();
		}
		
		if(distribution.contentEquals("BERNOULLI")) {
			est_pars = ((BernoulliNB) fittedNBC).get_est_pars();
		}
		
		return est_pars;
	}
	
}
