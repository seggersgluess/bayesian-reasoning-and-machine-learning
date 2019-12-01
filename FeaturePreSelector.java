package ObjectDetection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class FeaturePreSelector {
	
	List<Double> f_scores;
	
	int nBestFeatures;
	List<Double> bestFeatureScores;
	List<Integer> bestIdxs;
	
	
	public FeaturePreSelector(List<Double> binaryLabels, ArrayList<List<Double>> features, int nBestFeatures) {
		calc_f_scores4BinaryClassification(binaryLabels,features);
		select_bestFeatures(nBestFeatures);
	}
	
	
	public void calc_f_scores4BinaryClassification(List<Double> binaryLabels, ArrayList<List<Double>> features) {
		
		if(binaryLabels == null) {
			throw new RuntimeException("No binary labels set for feature selection.");
		}
		
		if(features == null) {
			throw new RuntimeException("No features set for feature selection.");
		}
		
		List<Double> classes = Utilities.Utilities.get_unique_sorted_elements_of_double_list(binaryLabels);
		
		if(classes.size()!=2) {
			throw new RuntimeException("No binary data for y supplied.");
		}	
		
		f_scores = new ArrayList<Double>();
		
		int nSamples = binaryLabels.size();
		
		int nFeatures = features.get(0).size();
		
		for(int f=0; f<nFeatures; f++) {
			
			List<Integer> posIdxs = new ArrayList<Integer>();
			List<Integer> negIdxs = new ArrayList<Integer>();
			
			double mean = 0.0;
			double mean_pos = 0.0;
			double mean_neg = 0.0;
			
			for(int i=0; i<nSamples; i++) {
				double label = binaryLabels.get(i);
				double feature = features.get(i).get(f);
				if(label == classes.get(0)) {
					posIdxs.add(i);
					mean_pos += feature;
				}
				if(label == classes.get(1)) {
					negIdxs.add(i);
					mean_neg += feature;
				}	
				mean += feature;
			}
			
			int nPos = posIdxs.size();
			int nNeg = negIdxs.size();
			
			mean_pos /= nPos;
			mean_neg /= nNeg;
			mean /= nSamples;
			
			double pos_diff = mean_pos-mean;
		    double squared_pos_diff = pos_diff*pos_diff;
		    double neg_diff = mean_neg-mean;
		    double squared_neg_diff = neg_diff*neg_diff;
			double nominator = squared_pos_diff + squared_neg_diff;
			
			double var_pos = 0.0;
			double var_neg = 0.0;
			
			for(int i=0; i<nPos; i++) {
				int idx = posIdxs.get(i);
				double diff = features.get(idx).get(f)-mean_pos;
				double squared_diff = diff*diff;
				var_pos += squared_diff;
			}
			var_pos /= (nPos-1.0);
			
			for(int i=0; i<nNeg; i++) {
				int idx = negIdxs.get(i);
				double diff = features.get(idx).get(f)-mean_pos;
				double squared_diff = diff*diff;
				var_neg += squared_diff;
			}
			var_neg /= (nNeg-1.0);
			
			double denominator = var_pos + var_neg;
			
			double f_score = nominator/denominator;
			f_scores.add(f_score);
			
		}
		
	}
		
	
	public void select_bestFeatures(int nBestFeatures) {
		
		if(f_scores == null) {
			System.out.println("No F-scores calculated yet.");
		}
			
		this.nBestFeatures = nBestFeatures;
		
		HashMap<String, List<Double>> sortedScoreInfos = Utilities.Utilities.get_sorted_elements_and_idxs_of_double_list(f_scores);
		int nScores = sortedScoreInfos.get("SortedValues").size()-1;
		
		bestFeatureScores = new ArrayList<Double>();
		bestIdxs = new ArrayList<Integer>();
		
		for(int i=0; i<nBestFeatures; i++) {
			int idx = nScores-i;
			bestFeatureScores.add(sortedScoreInfos.get("SortedValues").get(idx));
			int featureIdx = (int) Math.round(sortedScoreInfos.get("Idxs").get(idx));
			bestIdxs.add(featureIdx);
		}
		
	}
		
	
	public List<Integer> get_selectedBest_f_scores_Idxs() {
		return bestIdxs;
	}
	
	
	public List<Double> get_selectedBest_f_scores() {
		return bestFeatureScores;
	}
	
	public List<Double> get_f_scores() {
		return f_scores;
	}
	
}
