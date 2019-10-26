package AdaptiveBasisModels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class AdaBoost extends CART {

	public static int iterations = 1;
	
	public static ArrayList<HashMap<String,List<String>>> classifierInfos;
	public static List<Double> alphas;
	
	public static double [][] classifierValues;
	public static double [][] predicted_var;
	
	
	public static void doAdaBoost() {
		
		classifierValues = new double [n_observations][iterations];
		predicted_var = new double [n_observations][1];
		alphas = new ArrayList<Double>();
		
		List<Integer> rootIdxs = new ArrayList<Integer>(n_observations);
		double [][] weights = new double [1][n_observations];
		
		for(int i=0; i<n_observations; i++) {
			rootIdxs.add(i);
			//Start weights:
			weights[0][i] = 1.0/n_observations;
		}
		
		ArrayList<HashMap<String,List<String>>> stump_infos = get_decision_stump(rootIdxs, true);
		int nSplits = stump_infos.size();
		
		double [][] errorMatrix = new double [n_observations][nSplits];
		double [][] weightedErrorMatrix = new double [n_observations][nSplits]; 
		
		//binary class data {-1,+1}		
		for(int s=0; s<nSplits; s++) {			
			List<String> leftIdxs = stump_infos.get(s).get("LeftIdxs");
			for(int i=0; i<n_observations; i++) {
				double y_pred = 1.0;
				boolean isClass1 = leftIdxs.contains(Integer.toString(i));
				
				if(isClass1 == false) {
					y_pred = -1.0;
				}
				
				if(explained_variable[i][0] != y_pred) {
					errorMatrix[i][s] = 1.0;
				}				
			}
		}
		
		for(int i=0; i<iterations; i++) {
			
			weightedErrorMatrix = MatrixOperations.multiplication(weights, errorMatrix);
			
			double minIdx = Utilities.getMinWithIdx(weightedErrorMatrix).get("colIdx");
			
			classifierInfos.add(stump_infos.get((int) minIdx));
			
			double summedWeights = 0.0;
			
			for(int j=0; j<n_observations; j++) {
				summedWeights += weights[0][j];
			}
			
			double error = weightedErrorMatrix[0][(int) minIdx]/summedWeights;			
			double alpha = Math.log((1.0-error)/error);
			
			List<String> leftIdxs = stump_infos.get((int) minIdx).get("LeftIdxs");
			List<String> rightIdxs = stump_infos.get((int) minIdx).get("RightIdxs");
			
			int nLeftIdxs = leftIdxs.size();
			int nRightIdxs = rightIdxs.size();
			
			for(int j=0; j<nLeftIdxs; j++) {
				classifierValues[Integer.parseInt(leftIdxs.get(j))][i] = -1.0;
			}
			
			for(int j=0; j<nRightIdxs; j++) {
				classifierValues[Integer.parseInt(rightIdxs.get(j))][i] = 1.0;
			}
			
			for(int j=0; j<n_observations; j++) {
				weights[0][j] *= Math.exp(alpha*errorMatrix[j][(int) minIdx]);
				predicted_var[j][0] += alpha*classifierValues[j][i];
			}
			
			alphas.add(alpha);
			
		}
		
		for(int i=0; i<n_observations; i++) {
			if(predicted_var[i][0] < 0.0) {
				predicted_var[i][0] = -1.0;
			}
			
		}
		
	}
	
	
	//Make prediction for y with new_X (has to be a 1xn_explaining_vars vector)
	public static double make_predictionWithDecisionStumpInfos(double [][] new_X) {
		
		if(classifierInfos == null) {
			System.out.println("AdaBoost not trained yet. Can´t make prediction from new input.");
		}
		
		double prediction = 0.0;
		
		for(int i=0; i<iterations; i++) {
			
			String splitFeat = get_splittingFeature4Iteration(i);
			double threshold = get_threshold4Iteration(i);
			
			int idx = Utilities.get_idx(names_of_explaining_variables, splitFeat)[0];
			double inputValue4Feature = new_X[0][idx];
			
			if(inputValue4Feature <= threshold) {
				prediction -= alphas.get(i);
			}else {
				prediction += alphas.get(i);
			}
						
		}
		
		return prediction;
		
	}
	
	
	public static void read_AdaBoost_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception {
		read_CART_input_data(true, fileName, hasRowNames, hasColNames);
	}
	
	
	public static void set_AdaBoost_sampleName2InputData(String sampleName) {
		set_CART_sampleName(sampleName);
	}
	
	
	public static String get_splittingFeature4Iteration(int iteration) {
		String feature = classifierInfos.get(iteration).get("SplittingFeature").get(0);
		return feature;
	}
	
	
	public static double get_threshold4Iteration(int iteration) {
		double threshold = Double.parseDouble(classifierInfos.get(iteration).get("Threshold").get(0));
		return threshold;
	}
	
	
	public static double [][] get_predicted_var() {
		return predicted_var;
	}
	
	
	public static void set_numberOfIterations4AdaBoost(int nIterations) {
		
		if(nIterations<0) {
			throw new RuntimeException("Invalid number of iterations.");
		}
		iterations = nIterations;
	}
	
}
