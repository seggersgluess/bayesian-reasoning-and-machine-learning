package AdaptiveBasisModels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import Regression.LinearRegression;

public class AdaBoost extends CART {

	public static int iterations = 1;
	
	public static ArrayList<CART> treeClassifierInfos;
	public static List<Double> alphas;
	
	public static String algorithm;
	
	//InSample predictions of y_i
	public static double [][] predicted_var;
	
	
	public static void doAdaBoost(int treeDepth, String usedAlgorithm, int n_weakEstimators) {
		
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
		
		if(usedAlgorithm.contentEquals("SAMME") == true) {
			doSAMME(treeDepth);
		}
		
		if(usedAlgorithm.contentEquals("SAMME.R") == true) {
			doSAMME_R(treeDepth);
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
	
	
	public static void makeInSamplePredictionFromAdaBoost() {
		
		predicted_var = new double [n_observations][1];
		
		for(int i=0; i<n_observations; i++) {
			double [][] x = new double [1][n_explaining_variables];
			for(int j=0; j<n_explaining_variables; j++) {
				x[0][j] = explaining_variables[i][j];
			}
			if(algorithm.contentEquals("SAMME") == true) {
				predicted_var[i][0] = make_predictionFromSAMME(x).get("ClassPrediction")[0][0];
			}
			if(algorithm.contentEquals("SAMME.R") == true) {
				predicted_var[i][0] = make_predictionFromSAMME_R(x).get("ClassPrediction")[0][0];
			}
		}
		
	}
	
	
	public static double getAccuracyRate4AdaBoost() {
		
		if(predicted_var == null) {
			System.out.println("No insample prediction done yet. Make insamle prediction from AdaBoost now.");
			makeInSamplePredictionFromAdaBoost();
		}
		
		double accuracyRate = 0.0;
		double validClassCount = 0.0;
		
		for(int i=0; i<n_observations; i++) {
			if(predicted_var[i][0] == explained_variable[i][0]) {
				validClassCount++;
			}
		}
		
		accuracyRate = validClassCount/n_observations;
		
		return accuracyRate;
		
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
		
		double prediction = Double.MIN_VALUE;
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
			
		double sumOfProbPreds = 0.0;
		for(int c=0; c<nClasses; c++) {
			double newPred = 0.0;
			for(int i=0; i<nIterations; i++) {
				if(predWeakClassifier[i][0] == c) {
					newPred += alphas.get(i);
				}	
			}
			if(newPred>prediction) {
				prediction = newPred;
			}
			sumOfProbPreds += summedProbPredictions[c][0];
		}
		
		for(int c=0; c<nClasses; c++) {
			summedProbPredictions[c][0] /= sumOfProbPreds;
		}
		
		double [][] classPred = new double [1][1];
		classPred[0][0] = prediction;
	
		predResults.put("ClassPrediction", classPred);
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
		
		double prediction = Double.MIN_VALUE;
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
		
		double sumOfProbPreds = 0.0;
		for(int c=0; c<nClasses; c++) {
			double newPred = 0.0;
			for(int i=0; i<nIterations; i++) {
				double prob = probPredictions.get(i)[c][0];
				if(prob == 0.0) {
					prob = Double.MIN_VALUE;
				}else {
					prob = Math.log(prob);
				}
				newPred += (nClasses-1.0)*(prob-1.0/nClasses*logSumsOfProbPredictions.get(i));
			}
			if(newPred>prediction) {
				prediction = newPred;
			}
			sumOfProbPreds += summedProbPredictions[c][0];
		}
		
		for(int c=0; c<nClasses; c++) {
			summedProbPredictions[c][0] /= sumOfProbPreds;
		}
		
		double [][] classPred = new double [1][1];
		classPred[0][0] = prediction;
	
		predResults.put("ClassPrediction", classPred);
		predResults.put("ProbabilityPrediction", summedProbPredictions);
		
		return predResults;
		
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
			prediction = make_predictionFromSAMME(new_x).get("ClassPrediction")[0][0];
		}
		
		if(algorithm.contentEquals("SAMME.R") == true) {
			prediction = make_predictionFromSAMME_R(new_x).get("ClassPrediction")[0][0];
		}
		
		return prediction;
		
	}
	
	
	public static double [][] getAdaBoostPredictedClassProbabilities(double [][] new_x, boolean log) {
		
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
		
		if(log == true) {
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
	
	
	public static String [] getValidAlgorithms4AdaBoost() {
		
		String [] validAlgorithms = {"SAMME",
				                     "SAMME.R"};
		
		return validAlgorithms;
		
	}
	
	
	public static void set_classes4AdaBoost(double [] inputClasses) {
		classes = inputClasses;
	}
	
	
	public static void read_AdaBoost_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception {
		read_CART_input_data(true, fileName, hasRowNames, hasColNames);
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
	
	
	public static void set_numberOfIterations4AdaBoost(int nIterations) {
		
		if(nIterations<0) {
			throw new RuntimeException("Invalid number of iterations.");
		}
		iterations = nIterations;
	}
		
}
