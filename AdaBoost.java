package AdaptiveBasisModels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Mathematics.GeneralMath;
import Regression.LinearRegression;

public class AdaBoost extends CART {

	public static int iterations = 1;
	
	public static ArrayList<CART> treeClassifierInfos;
	public static List<Double> alphas;
	
	//Betas for AdaBoostRegression
	public static List<Double> betas;
	//Original y & X (before resampling)
	public static double [][] org_y;
	public static double [][] org_X;
	
	public static String algorithm;
	
	//InSample predictions of y_i after each boosting iteration i=1,2,...,iterations
	public static double [][] predicted_var;
	public static double [][] accuracyRates;
	
	
	//doAdaBoost without regLossType (for classification or regression with default LINEAR)
	public static void doAdaBoost(int treeDepth, String usedAlgorithm, int n_weakEstimators) {
		
		if(categorical_explained_var == true) {
			doAdaBoost(treeDepth, usedAlgorithm, n_weakEstimators, null);
		}else {
			doAdaBoost(treeDepth, usedAlgorithm, n_weakEstimators, "LINEAR");
		}
		
	}
	
	
	public static void doAdaBoost(int treeDepth, String usedAlgorithm, int n_weakEstimators, String regressionLossType) {
		
		if(treeDepth <= 0) {
			System.out.println("Invalid tree depth supplied.");
			return;
		}
		
		String [] validAlgorithms = getValidAlgorithms4AdaBoost();
		int [] idxs = Utilities.Utilities.get_idx(validAlgorithms, usedAlgorithm);
		if(idxs[0] == -1) {
			System.out.println(usedAlgorithm + " is not a valid algorithm for AdaBoost.");
			return;
		}
		
		if(n_weakEstimators <= 0) {
			System.out.println("Invalid number of weak estimators supplied.");
			return;
		}
		
		if(categorical_explained_var == true) {
			String [] classificationAlogrithms = {"SAMME", "SAMME.R"};
			idxs = Utilities.Utilities.get_idx(classificationAlogrithms, usedAlgorithm); 
			if(idxs[0] == -1) {
				System.out.println(usedAlgorithm + " is not a valid algorithm for classification in AdaBoost.");
				return;
			}
					
			if(usedAlgorithm.contentEquals("SAMME") == true) {
				doSAMME(treeDepth);
			}
			
			if(usedAlgorithm.contentEquals("SAMME.R") == true) {
				doSAMME_R(treeDepth);
			}
		}else {
			String [] classificationAlogrithms = {"REGRESSION"};
			idxs = Utilities.Utilities.get_idx(classificationAlogrithms, usedAlgorithm); 
			if(idxs[0] == -1) {
				System.out.println(usedAlgorithm + " is not a valid algorithm for regression in AdaBoost.");
				return;
			}
			doAdaBoostRegression(treeDepth,regressionLossType);
		}
		
	}
	
	
	//Conventional AdaBoost with decision stump
	public static void doSamme() {
		doSAMME(1);
	}
	
	
	//Stagewise Additive Modeling using Multi-class Exponential loss function (SAMME)
	@SuppressWarnings("static-access")
	public static void doSAMME(int treeDepth) {
		
		if(n_observations == 0) {
			System.out.println("No data for SAMME uploaded yet.");
			return;
		}
		
		algorithm = "SAMME";
		
		treeClassifierInfos = new ArrayList<CART>();
		
		weights = new ArrayList<Double>();
		
		double initWeight = 1.0/n_observations;
		
		for(int i=0; i<n_observations; i++) {
			weights.add(initWeight);
		}
		
		for(int i=0; i<iterations; i++) {			
			
			List<Double> updated_weights = new ArrayList<Double>();
			
			CART obj_CART = new CART();
			obj_CART.useCategoricalTree();
			obj_CART.set_CART_explained_variable(name_of_explained_variable);
			obj_CART.set_CART_inputData();
			
			obj_CART.set_minNumberOfElementsInKnot(1);
			obj_CART.set_maxNumberOfClassesInKnot(1);
			obj_CART.set_maxTreeDepth(treeDepth);
			
			obj_CART.fit_tree(false);			
			obj_CART.calc_least_square_regressor_weights();
			
			LinearRegression obj_linearReg = obj_CART.get_linearRegObject();
			double [][] fitted_explained_var = obj_linearReg.get_fitted_values();
			
			double error = 0.0;
			
			for(int j=0; j<n_observations; j++) {
				if(fitted_explained_var[j][0] != explained_variable[j][0]) {
					error += weights.get(j);
				}
			}
			
			double alpha = Math.log((1.0-error)/error) + Math.log(nClasses-1);
		    alphas.add(alpha);
		    double summed_weights = 0.0;
		    
		    for(int j=0; j<n_observations; j++) {
		    	double updated_weight = weights.get(j);
		    	if(fitted_explained_var[j][0] != explained_variable[j][0]) {
					updated_weight *= Math.exp(alpha);
				}	
		    	updated_weights.add(updated_weight);
		    	summed_weights += updated_weight;
		    }
		    
		    weights = new ArrayList<Double>();
		    for(int j=0; j<n_observations; j++) {
		    	weights.add(updated_weights.get(j)/summed_weights);
		    }
		    	
		    //TODO: Check if tree´s for different iterations are not the same
		    treeClassifierInfos.add(obj_CART);
		    
		}	
				
	}
	
	
	@SuppressWarnings("static-access")
	public static void doSAMME_R(int treeDepth) {
		
		if(n_observations == 0) {
			System.out.println("No data for SAMME.R uploaded yet.");
			return;
		}
		
		algorithm = "SAMME.R";
		
		treeClassifierInfos = new ArrayList<CART>();
		
		weights = new ArrayList<Double>();
		
		double initWeight = 1.0/n_observations;
		
		for(int i=0; i<n_observations; i++) {
			weights.add(initWeight);
		}
		
		for(int i=0; i<iterations; i++) {			
			
			List<Double> updated_weights = new ArrayList<Double>();
			double summed_weights = 0.0;
			
			CART obj_CART = new CART();
			obj_CART.useCategoricalTree();
			obj_CART.set_CART_explained_variable(name_of_explained_variable);
			obj_CART.set_CART_inputData();
			
			obj_CART.set_minNumberOfElementsInKnot(1);
			obj_CART.set_maxNumberOfClassesInKnot(1);
			obj_CART.set_maxTreeDepth(treeDepth);
			
			obj_CART.fit_tree(false);			
			
			double [] y = new double [nClasses];
			
			for(int c=0; c<nClasses; c++) {
				y[c] = -nClasses/(nClasses-1.0);
			}
			
			for(int j=0; j<n_observations; j++) {
				double [] y_i = y;
				int curClass = (int)explained_variable[j][0];
				y_i[curClass] = 1.0;
				
				double [][] x_i = new double [1][n_explaining_variables];
				for(int k=0; k<n_explaining_variables; k++) {
					x_i[0][k] = explaining_variables[j][k];
				}
				
				double [][] probPred = obj_CART.makeProbabilityPrediction(x_i);
				
				double scalarProd = 0.0;
				for(int k=0; k<nClasses; k++) {
					double prob = probPred[k][0];
					if(prob == 0.0) {
						prob = Double.MIN_NORMAL;
					}else {
						prob = Math.log(prob);
					}
					scalarProd += y_i[k]*prob;
				}
				
				scalarProd *= -(nClasses-1.0)/nClasses;
				scalarProd = Math.exp(scalarProd);
				
				double updated_weight = weights.get(i)*scalarProd;
				updated_weights.add(updated_weight);
				summed_weights += updated_weight;
				
			}

			weights = new ArrayList<Double>();
		    for(int j=0; j<n_observations; j++) {
		    	weights.add(updated_weights.get(j)/summed_weights);
		    }
		    	
		    //TODO: Check if tree´s for different iterations are not the same
		    treeClassifierInfos.add(obj_CART);
		    
		}	
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void doAdaBoostRegression(int treeDepth, String lossType) {
		
		if(n_observations == 0) {
			System.out.println("No data for AdaBoost Regression uploaded yet.");
			return;
		}
		
		String [] validLossTypes = get_validLossTypes4AdaBoostReg();
		int [] idxs = Utilities.Utilities.get_idx(validLossTypes, lossType);
		if(idxs[0] == -1) {
			System.out.println(lossType + " is not a valid loss type. Use type linear in the following.");
			lossType = "LINEAR";
		}
		
		algorithm = "REGRESSION";
		
		treeClassifierInfos = new ArrayList<CART>();
		
		weights = new ArrayList<Double>();
		
		double initWeight = 1.0/n_observations;
		
		for(int i=0; i<n_observations; i++) {
			weights.add(initWeight);
		}
		
		//TODO: Check if y & X do not change after resampling!
		org_y = explained_variable;
		org_X = explaining_variables;
		
		for(int i=0; i<iterations; i++) {			
			
			List<Double> updated_weights = new ArrayList<Double>();
			double summed_weights = 0.0;
			
			double [][] new_X = new double [n_observations][n_explaining_variables];
			double [][] new_y = new double [n_observations][1];
			
			//Sample new set {y,X}
			int [] sampledIdxs = Utilities.Utilities.getRandomIntNumbers(weights, n_observations);
			
			for(int j=0; j<n_observations; j++) {
				int idx = sampledIdxs[j];
				new_y[j][0] = org_y[idx][0];
				for(int k=0; k<n_explaining_variables; k++) {
					new_X[j][k] = org_X[idx][k];
				}
			}
			
			CART obj_CART = new CART();
			obj_CART.set_CART_explained_variable(name_of_explained_variable);
			obj_CART.set_CART_inputData();
			obj_CART.set_CART_explaining_variables_data(new_X);
			obj_CART.set_CART_explained_variable_data(new_y);
			
			obj_CART.set_minNumberOfElementsInKnot(1);
			obj_CART.set_maxNumberOfClassesInKnot(1);
			obj_CART.set_maxTreeDepth(treeDepth);
			
			obj_CART.fit_tree(false);			
					
			HashMap<String, double [][] >  lossInfos = loss4AdaBoostRegression(i, lossType);
			
			double [][] lossFunc = lossInfos.get("Loss");
			double meanLoss = lossInfos.get("MeanLoss")[0][0];
			double beta = meanLoss/(1.0-meanLoss);
			betas.add(beta);
						
			for(int j=0; j<n_observations; j++) {
				double weight = weights.get(j);
				double updated_weight = weight*Math.pow(beta, (1.0-lossFunc[j][0]));
				updated_weights.add(updated_weight);
				summed_weights += updated_weight;			
			}

			weights = new ArrayList<Double>();
		    for(int j=0; j<n_observations; j++) {
		    	weights.add(updated_weights.get(j)/summed_weights);
		    }
		    	
		    //TODO: Check if tree´s for different iterations are not the same
		    treeClassifierInfos.add(obj_CART);
		    
		}	
			
	}
	
	
	@SuppressWarnings("static-access")
	public static HashMap<String, double [][] > loss4AdaBoostRegression(int iteration, String lossType) {
		
		HashMap<String, double [][] > lossInfos = new HashMap<String, double [][] >();
		
		double [][] lossFunc = new double [n_observations][1];
		double maxLoss  = Double.MIN_VALUE;
		double [][] y = treeClassifierInfos.get(iteration).get_explained_variable();
		
		for(int i=0; i<n_observations; i++) {
			double [][] x_i = new double [1][n_explaining_variables];
			for(int k=0; k<n_explaining_variables; k++) {
				x_i[0][k] = explaining_variables[i][k];
			}
			
			double predValue = treeClassifierInfos.get(iteration).makePrediction(x_i);
			lossFunc[i][0] = Math.abs(predValue-y[i][0]); 
			
			if(lossFunc[i][0]>maxLoss) {
				maxLoss = lossFunc[i][0];
			}
			
		}
		
		double [][] meanLoss = new double [1][1];
		
		if(lossType == "LINEAR") {
			for(int i=0; i<n_observations; i++) {
				lossFunc[i][0] /= maxLoss;
				meanLoss[0][0] += weights.get(i)*lossFunc[i][0];
			}
		}
		
		if(lossType == "SQUARED") {
			for(int i=0; i<n_observations; i++) {
				lossFunc[i][0] = lossFunc[i][0]*lossFunc[i][0]/(maxLoss*maxLoss);
				meanLoss[0][0] += weights.get(i)*lossFunc[i][0];
			}
		}
		
		if(lossType == "EXPONENTIAL") {
			for(int i=0; i<n_observations; i++) {
				double lossValue = - lossFunc[i][0]/maxLoss;
				lossFunc[i][0] = 1.0-Math.exp(lossValue);
				meanLoss[0][0] += weights.get(i)*lossFunc[i][0];
			}
		}
		
		lossInfos.put("Loss", lossFunc);
		lossInfos.put("MeanLoss", meanLoss);
		
		return lossInfos;
		
	}
	
	
	public static String [] get_validLossTypes4AdaBoostReg() {
		String [] validLossTypes = {"LINEAR",
				                    "SQUARED",
				                    "EXPONENTIAL"};
		return validLossTypes;
	}
	
	
	public static void makeInSamplePredictionFromAdaBoost() {
		
		//AdaBoost classification
		if(categorical_explained_var == true) {
			predicted_var = new double [n_observations][iterations];
			double [][] pred_var = new double [1][iterations];
			
			for(int i=0; i<n_observations; i++) {
				double [][] x = new double [1][n_explaining_variables];
				for(int j=0; j<n_explaining_variables; j++) {
					x[0][j] = explaining_variables[i][j];
				}
				if(algorithm.contentEquals("SAMME") == true) {
					pred_var = make_predictionFromSAMME(x).get("ClassPredictions");		
				}
				if(algorithm.contentEquals("SAMME.R") == true) {
					pred_var = make_predictionFromSAMME_R(x).get("ClassPredictions");
				}
				for(int j=0; j<iterations; j++) {
					predicted_var[i][j] = pred_var[0][j];
				}
			}
			
			getAccuracyRates4AdaBoostClassification();
		}
		//AdaBoost regression
		if(categorical_explained_var == false) {
			predicted_var = new double [n_observations][1];
			double summedSquaredError = 0.0;
			double summedSquaredMeanDev = 0.0;
			double mean_y = GeneralMath.mean(org_y);
			
			for(int i=0; i<n_observations; i++) {
				double [][] x = new double [1][n_explaining_variables];
				for(int j=0; j<n_explaining_variables; j++) {
					x[0][j] = explaining_variables[i][j];
				}
				if(algorithm.contentEquals("REGRESSION") == true) {
					predicted_var[i][0] = make_predictionFromAdaBoostRegression(x);
				}
				summedSquaredError += Math.pow(predicted_var[i][0]-org_y[i][0],2.0);
				summedSquaredMeanDev += Math.pow(org_y[i][0]-mean_y, 2.0);
			}
			
			//R^2 (Measure of Determination)
			accuracyRates = new double [1][1];
			accuracyRates[0][0] = 1.0-summedSquaredError/summedSquaredMeanDev;			
		}
	
	}
	
	
	//returns accuracy rates for each of the i=1,2,...,m boosting iterations
	public static void getAccuracyRates4AdaBoostClassification() {
		
		if(predicted_var == null) {
			System.out.println("No insample prediction done yet. Make insamle prediction from AdaBoost now.");
			makeInSamplePredictionFromAdaBoost();
		}
		
		accuracyRates = new double [1][iterations];
		
		for(int i=0; i<iterations; i++) {
			double validClassCount = 0.0;
			for(int j=0; j<n_observations; j++) {
				if(predicted_var[j][i] == explained_variable[j][0]) {
					validClassCount++;
				}
			}
			accuracyRates[0][i] = validClassCount/n_observations;
		}
						
	}
	
	
	//Make prediction for y with new_x (has to be a 1xn_explaining_vars vector) from strong classifier for SAMME
	public static HashMap<String, double [][]> make_predictionFromSAMME(double [][] new_x) {
		
		if(treeClassifierInfos == null) {
			System.out.println("AdaBoost not trained yet. Can´t make prediction from new input.");
			return null;
		}
		
		if(algorithm.contentEquals("SAMME") == false) {
			System.out.println("Selected algorithm of AdaBoost is not SAMME. Cannot make prediction based on SAMME.");
			return null;
		}
		
		HashMap<String, double [][]> predResults = new HashMap<String, double [][]>();
		
		double classifierValue = Double.MIN_VALUE;
		double classPrediction = 0.0;
		int nIterations = alphas.size();
		
		double [][] predWeakClassifier = new double [nIterations][1];
		double [][] summedProbPredictions = new double [nClasses][1];
		
		for(int i=0; i<nIterations; i++) {
			predWeakClassifier[i][0] = treeClassifierInfos.get(i).makePrediction(new_x);	
			double [][] probPredictions = treeClassifierInfos.get(i).makeProbabilityPrediction(new_x);
			for(int c=0; c<nClasses; c++) {
				summedProbPredictions[c][0] += alphas.get(i)*probPredictions[c][0];
			}
		}
			
		double [][] classPredictions = new double [1][nIterations];
		double sumOfProbPreds = 0.0;
		for(int i=0; i<nIterations; i++) {
			double newValue = 0.0;
			for(int c=0; c<nClasses; c++) {
				if(predWeakClassifier[i][0] == c) {
					newValue += alphas.get(i);
				}	
				if(newValue>classifierValue) {
					classifierValue = newValue;
					classPrediction = c;
				}
				classPredictions[0][i] = classPrediction;
			}
		}
		
		for(int c=0; c<nClasses; c++) {
			sumOfProbPreds += summedProbPredictions[c][0];
		}
		
		for(int c=0; c<nClasses; c++) {
			summedProbPredictions[c][0] /= sumOfProbPreds;
		}
		
		double [][] classPred = new double [1][1];
		classPred[0][0] = classPrediction;
	
		predResults.put("FinalClassPrediction", classPred);
		predResults.put("ClassPredictions", classPredictions);
		predResults.put("ProbabilityPrediction", summedProbPredictions);
		
		return predResults;
		
	}
	
	
	//Make prediction for y with new_x (has to be a 1xn_explaining_vars vector) from strong classifier for SAMME.R
	public static HashMap<String, double [][]> make_predictionFromSAMME_R(double [][] new_x) {
			
		if(treeClassifierInfos == null) {
			System.out.println("AdaBoost not trained yet. Can´t make prediction from new input.");
			return null;
		}
		
		if(algorithm.contentEquals("SAMME.R") == false) {
			System.out.println("Selected algorithm of AdaBoost is not SAMME. Cannot make prediction based on SAMME.");
			return null;
		}
				
		HashMap<String, double [][]> predResults = new HashMap<String, double [][]>();
		
		double classifierValue = Double.MIN_VALUE;
		double classPrediction = 0.0;
		int nIterations = treeClassifierInfos.size();
		
		ArrayList<double [][]> probPredictions = new ArrayList<double [][]>();
		double [][] summedProbPredictions = new double [nClasses][1];
		ArrayList<Double> logSumsOfProbPredictions = new ArrayList<Double>();
		
		for(int i=0; i<nIterations; i++) {
			probPredictions.add(treeClassifierInfos.get(i).makeProbabilityPrediction(new_x));
			double logProbSum = 0.0;
			for(int c=0; c<nClasses; c++) {
				double prob = probPredictions.get(i)[c][0];
				summedProbPredictions[c][0] += prob;
				if(prob == 0.0) {
					logProbSum += Double.MIN_VALUE;
				}else {
					logProbSum += Math.log(prob);
				}
			}
			logSumsOfProbPredictions.add(logProbSum);
		}
		
		double [][] classPredictions = new double [1][nIterations];
		double sumOfProbPreds = 0.0;
		for(int i=0; i<nIterations; i++) {
			double newValue = 0.0;
			for(int c=0; c<nClasses; c++) {
				double prob = probPredictions.get(i)[c][0];
				if(prob == 0.0) {
					prob = Double.MIN_VALUE;
				}else {
					prob = Math.log(prob);
				}
				newValue += (nClasses-1.0)*(prob-1.0/nClasses*logSumsOfProbPredictions.get(i));
				if(newValue>classifierValue) {
					classifierValue = newValue;
					classPrediction = c;
				}
				classPredictions[0][i] = classPrediction;
			}
		}
		
		for(int c=0; c<nClasses; c++) {
			sumOfProbPreds += summedProbPredictions[c][0];
		}
		
		for(int c=0; c<nClasses; c++) {
			summedProbPredictions[c][0] /= sumOfProbPreds;
		}
		
		double [][] classPred = new double [1][1];
		classPred[0][0] = classPrediction;
	
		predResults.put("FinalClassPrediction", classPred);
		predResults.put("ClassPredictions", classPredictions);
		predResults.put("ProbabilityPrediction", summedProbPredictions);
		
		return predResults;
		
	}
	
	
	//Make prediction for y with new_x (has to be a 1xn_explaining_vars vector) from strong classifier of AdaBoostRegression algorithm
	public static double make_predictionFromAdaBoostRegression(double [][] new_x) {
		
		List<Double> predWeakClassifier = new ArrayList<Double>();
		double summedLogBetas = 0.0;
		for(int i=0; i<iterations; i++) {
			double beta = betas.get(i);
			if(beta == 0.0) {
				beta = Double.MIN_VALUE;
			}else {
				beta = Math.log(beta);
			}
			summedLogBetas += beta;
			predWeakClassifier.add(treeClassifierInfos.get(i).makePrediction(new_x));	
		}
		
		summedLogBetas *= 0.5;
		
		List<Double> sortedIdxs = Utilities.Utilities.get_sorted_elements_and_idxs_of_double_list(predWeakClassifier).get("Idxs");
		
		double partSummedLogBetas = 0.0;
		int predClassifier = -1;
		for(int i=0; i<iterations; i++) {
			int idx = (int) Math.round(sortedIdxs.get(i));
			double beta = betas.get(idx);
			if(beta == 0.0) {
				beta = Double.MIN_VALUE;
			}else {
				beta = Math.log(beta);
			}
			partSummedLogBetas += beta;
			if(partSummedLogBetas>=summedLogBetas) {
				predClassifier = idx;
			}
		}
		
		double prediction = predWeakClassifier.get(predClassifier);
		
		return prediction;
		
	}
	
	
	//Make prediction for y with new_x (has to be a 1xn_explaining_vars vector) from strong classifier of AdaBoost algorithms
	public static double makeAdaBoostPrediction(double [][] new_x) {
		
		String [] validAlgorithms = getValidAlgorithms4AdaBoost();
		int [] idxs = Utilities.Utilities.get_idx(validAlgorithms, algorithm);
		if(idxs[0] == -1) {
			System.out.println("No a valid algorithm set for AdaBoost. Cannot make a prediction.");
			return -1.0;
		}
		
		double prediction = -1.0;
		
		if(algorithm.contentEquals("SAMME") == true) {
			prediction = make_predictionFromSAMME(new_x).get("FinalClassPrediction")[0][0];
		}
		
		if(algorithm.contentEquals("SAMME.R") == true) {
			prediction = make_predictionFromSAMME_R(new_x).get("FinalClassPrediction")[0][0];
		}
		
		if(algorithm.contentEquals("REGRESSION") == true) {
			prediction = make_predictionFromAdaBoostRegression(new_x);
		}
			
		return prediction;
		
	}
	
	
	public static double [][] getAdaBoostPredictedClassProbabilities(double [][] new_x, boolean logProbs) {
		
		double [][] predProbs = new double [nClasses][1];
		
		String [] validAlgorithms = getValidAlgorithms4AdaBoost();
		int [] idxs = Utilities.Utilities.get_idx(validAlgorithms, algorithm);
		if(idxs[0] == -1) {
			System.out.println("No a valid algorithm set for AdaBoost. Cannot make a prediction.");
			return null;
		}
		
		if(algorithm.contentEquals("SAMME") == true) {
			predProbs = make_predictionFromSAMME(new_x).get("ProbabilityPrediction");
		}
		
		if(algorithm.contentEquals("SAMME.R") == true) {
			predProbs = make_predictionFromSAMME_R(new_x).get("ProbabilityPrediction");
		}
		
		if(logProbs == true) {
			for(int c=0; c<nClasses; c++) {
				double prob = predProbs[c][0];
				if(prob == 0.0) {
					prob = Double.MIN_VALUE;
				}else {
					prob = Math.log(prob);
				}
				predProbs[c][0] = prob;
			}
		}
		
		return predProbs;
		
	}
	
	
	public static CART getWeakClassifier(int classifierNumber) {
		if(treeClassifierInfos == null) {
			System.out.println("No classifier infos found.");
			return null;
		}
		if(classifierNumber>iterations) {
			System.out.println("Supplied classifier number " + classifierNumber + " not valid.");
			return null;
		}
		
		CART weakClassifier = treeClassifierInfos.get(classifierNumber);
		
		return weakClassifier;
		
	}
	
	
	public static String [] getValidAlgorithms4AdaBoost() {
		
		String [] validAlgorithms = {"SAMME",
				                     "SAMME.R",
				                     "REGRESSION"};
		
		return validAlgorithms;
		
	}
	
	
	public static void set_classes4AdaBoost(double [] inputClasses) {
		classes = inputClasses;
	}
	
	
	public static void read_AdaBoost_input_data(boolean classData, String fileName, boolean hasRowNames, boolean hasColNames) throws Exception {
		read_CART_input_data(classData, fileName, hasRowNames, hasColNames);
	}
	
	
	public static void set_AdaBoost_sampleName2InputData(String sampleName) {
		set_CART_sampleName(sampleName);
	}
	
	
	public static double [][] get_predicted_var() {
		if(predicted_var == null) {
			System.out.println("No insample prediction found yet.");
			return null;
		}
		return predicted_var;
	}
	
	
	//Returns for classification accuracy rate (after every boosting step) and for regression R^2 (after final boosting step)
	public static double [][] get_accuracy_rates() {
		if(accuracyRates == null) {
			System.out.println("No accuracy rates found for AdaBoost yet.");
			return null;
		}
		return accuracyRates;
	}
	
	
	public static void set_numberOfIterations4AdaBoost(int nIterations) {
		
		if(nIterations<0) {
			throw new RuntimeException("Invalid number of iterations.");
		}
		iterations = nIterations;
	}
		
}
