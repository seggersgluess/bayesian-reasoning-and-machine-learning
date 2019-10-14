package AdaptiveBasisModels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import DataManagement.InputDataManager;
import Graphics.DecisionTreeGraphics;
import Regression.LinearRegression;

public class CART {

	static InputDataManager inputData;
	static String [] selectedRows;
	
	static double [][] explained_variable;
	static double [][] explaining_variables;	
	static String name_of_explained_variable;
	static String [] names_of_explaining_variables;
	
	static int n_observations;
	static int n_explaining_variables;
	
	static boolean categorical_explained_var;
	static double [] classes;
	static String [] classNames;
	static int nClasses;
	
	//User specified values to control tree growing
	static double minCostReduction = 0.0;
	static double maxNumberOfClassesInKnot = 0.0;
	static double minNumberOfElementsInKnot = 1.0;
	static int maxTreeDepth = 0;
	
	//Cost Measures
	static String costMeasure = "Gini";
	
	//Number of layers of the tree
	static int treeDepth;
	
	//Tree
	static HashMap <String, HashMap<String, HashMap<String,List<String>>>> tree;
	
	//CART model
	static boolean useConstant = true;
	static double [][] regressor_matrix;
	static LinearRegression obj_linearReg;
	
	//Graph control
	static boolean showLeafs = true;
	
	public static void fit_tree() {
		
		if(categorical_explained_var == false) {
			costMeasure = "Regression";
		}
		
		tree = new HashMap <String, HashMap<String, HashMap<String,List<String>>>>();
			
		boolean stop = false;
		
		//Get start idxs for root knot
		List<Integer> knotIdxs = new ArrayList<Integer>();
		List<String> strKnotIdxs = new ArrayList<String>();
		
		for(int i=0; i<n_observations; i++) {
			knotIdxs.add(i);
			strKnotIdxs.add(Integer.toString(i));
		}
			
		//--- Fill info for root knot ---
		HashMap<String,List<String>> knot = get_knot_structure();
		//Cost of knot
		String info = Double.toString(calc_cost(knotIdxs));	
		knot.get("Cost").add(info);
		//Number of elements in knot
		info = Integer.toString(n_observations);
		knot.get("NumberOfElements").add(info);	
		//Parent Knot
		knot.get("ParentKnot").add("Knot0");
		//Number of elements for classes in knot
		List<String> classDistribution = count_elements_of_classes(knotIdxs);
		knot.put("nElementsOfClasses",classDistribution);
		//Indices for elements in knot
		knot.put("IdxsOfElements", strKnotIdxs);
		//Threshold
		knot.get("Threshold").add("--");
		
		HashMap<String, HashMap<String,List<String>>> layer = new HashMap<String, HashMap<String,List<String>>>();
		layer.put("Knot1", knot);
		tree.put("Layer1", layer);		
		
		int layerCount = 1;
		
		while(stop == false) {
			
			String prevLayerLabel = "Layer" + layerCount;
			int nKnots = tree.get(prevLayerLabel).size();
			
			layerCount++;
			String layerLabel = "Layer" + layerCount;
			layer = new HashMap<String, HashMap<String,List<String>>>();
			
			int knotCount = 0;
			
			for(int k=0; k<nKnots; k++) {
				
				knotIdxs = new ArrayList<Integer>();
				String parentKnotLabel = "Knot" + (k+1);
				strKnotIdxs = tree.get(prevLayerLabel).get(parentKnotLabel).get("IdxsOfElements");
				int nElements = strKnotIdxs.size();
				for(int i=0; i<nElements; i++) {
					knotIdxs.add(Integer.parseInt(strKnotIdxs.get(i)));
				}
					
				boolean doSplit = true;
				Map<String,List<String>> splitInfos = get_split_infos_for_knot(knotIdxs);
				
				List<String> leftIdxs = splitInfos.get("IndicesLeftKnot");
				List<String> rightIdxs = splitInfos.get("IndicesRightKnot");
				
				int nLeftIdxs = leftIdxs.size();
				int nRightIdxs = rightIdxs.size();
				
				//1. Splitting criterion: Check minimum number of elements in knots
				if(minNumberOfElementsInKnot != 0.0) {
					if((nLeftIdxs < minNumberOfElementsInKnot) || (nRightIdxs  < minNumberOfElementsInKnot)) {
						doSplit = false;
					}
				}
				
				//2. Splitting criterion: Check cost reduction
				if(minCostReduction != 0.0) {
					double costReduction = Double.parseDouble(splitInfos.get("CostReduction").get(0));
					if(minCostReduction < costReduction) {
						doSplit = false;
					}
				}
				
				//3. Splitting criterion: Check homogenity of elements in knots
				if(maxNumberOfClassesInKnot != 0.0) {
								
					double [] leftExplainedVars = new double [nLeftIdxs];
				    double [] rightExplainedVars = new double [nRightIdxs];
				    
				    for(int i=0; i<nLeftIdxs; i++) {
				    	int idx = Integer.parseInt(leftIdxs.get(i));
				    	leftExplainedVars[i] = explained_variable[idx][0];
				    }
				    
				    for(int i=0; i<nRightIdxs; i++) {
				    	int idx = Integer.parseInt(rightIdxs.get(i));
				    	rightExplainedVars[i] = explained_variable[idx][0];
				    }
				    
					int nClassesLeftKnot  = count_classes_of_knot(leftExplainedVars);
					int nClassesRightKnot = count_classes_of_knot(rightExplainedVars);
					
					if((nClassesLeftKnot < maxNumberOfClassesInKnot) || (nClassesRightKnot < maxNumberOfClassesInKnot)) {
						doSplit = false;
					}
				}
				
				if(doSplit == true) {
					
					HashMap<String,List<String>> leftKnot = get_knot_structure();
					HashMap<String,List<String>> rightKnot = get_knot_structure();
					
					//SplittingFeature
					String splittingFeature = splitInfos.get("SplittingFeature").get(0);
					//tree.get(prevLayerLabel).get(parentKnotLabel).get("SplittingFeature").add(splittingFeature);
					leftKnot.get("SplittingFeature").add(splittingFeature);
					rightKnot.get("SplittingFeature").add(splittingFeature);
						
					//ParentKnot
					leftKnot.get("ParentKnot").add(parentKnotLabel);
					rightKnot.get("ParentKnot").add(parentKnotLabel);
					
					//Cost (absolute)
					leftKnot.get("Cost").add(splitInfos.get("Cost").get(1));
					rightKnot.get("Cost").add(splitInfos.get("Cost").get(2));
					
					//CostReduction (weighted by knots share of elements)
					leftKnot.get("CostReduction").add(splitInfos.get("CostReductionLeftKnot").get(0));
					rightKnot.get("CostReduction").add(splitInfos.get("CostReductionRightKnot").get(0));
					
					//Threshold
					leftKnot.get("Threshold").add(splitInfos.get("Threshold").get(0));
					rightKnot.get("Threshold").add(splitInfos.get("Threshold").get(0));
					
					//IdxsOfElements
					leftKnot.put("IdxsOfElements",leftIdxs);
					rightKnot.put("IdxsOfElements",rightIdxs);

					//NumberOfElements
					leftKnot.get("NumberOfElements").add(Integer.toString(nLeftIdxs));
					rightKnot.get("NumberOfElements").add(Integer.toString(nRightIdxs));
					
					//nElementsOfClasses	
					classDistribution = count_elements_of_classes_with_str_idxs(leftIdxs);
					leftKnot.put("nElementsOfClasses", classDistribution);
					classDistribution = count_elements_of_classes_with_str_idxs(rightIdxs);
					rightKnot.put("nElementsOfClasses", classDistribution);
					
					knotCount++;
					String knotLabel = "Knot"+knotCount;
					layer.put(knotLabel, leftKnot);
					knotCount++;
					knotLabel = "Knot"+knotCount;
					layer.put(knotLabel, rightKnot);
		
				}
				
			}
				
	        if(knotCount != 0) {
	        	tree.put(layerLabel, layer);
				//layerCount++;
	        }else {
	        	//Tree stops to grow (no further split)
	        	stop = true;
	        	layerCount--;
	        }
			
	        if(maxTreeDepth != 0) {
	        	if(layerCount == maxTreeDepth) {
	        		//Tree reached maximum depth
	        		stop = true;
	        	}
	        }
	        
		}
		
		treeDepth = layerCount;
		
	}
	
	
	public static HashMap<String, List<String>> get_knot_structure(){
		
		HashMap<String, List<String>> knot_struct = new HashMap<String, List<String>>();
		
		List<String> emptyList1 = new ArrayList<String>();		
		knot_struct.put("Cost", emptyList1);
		
		List<String> emptyList2 = new ArrayList<String>();
		knot_struct.put("SplittingFeature", emptyList2);
		
		List<String> emptyList3 = new ArrayList<String>();
		knot_struct.put("NumberOfElements", emptyList3);
		
		List<String> emptyList4 = new ArrayList<String>();
		knot_struct.put("nElementsOfClasses", emptyList4);
		
		List<String> emptyList5 = new ArrayList<String>();
		knot_struct.put("IdxsOfElements", emptyList5);
		
		List<String> emptyList6 = new ArrayList<String>();
		knot_struct.put("ParentKnot", emptyList6);
		
		List<String> emptyList7 = new ArrayList<String>();
		knot_struct.put("CostReduction", emptyList7);
		
		List<String> emptyList8 = new ArrayList<String>();
		knot_struct.put("Threshold", emptyList8);
		
		return knot_struct;
		
	}
	
	
    /**
     * Method for generating splitting infos
     * @param knotIdxs Indices of elements in current knot which has to be splitted 
     * @return ArrayList with 5 slots(
     * 1. Cost, 
     * 2. CostReduction,
     * 3. SplittingFeature,
     * 4. IndicesLeftKnot,
     * 5. IndicesRightKnot)
     */
	public static HashMap<String,List<String>> get_split_infos_for_knot(List<Integer> knotIdxs) {
		
		double minCost      = Double.MAX_VALUE;
		double minLeftCost  = Double.MAX_VALUE;
		double minRightCost = Double.MAX_VALUE;
		double overallCost  = Double.MAX_VALUE;
		double threshold    = 0.0;
		
		int knotSampleLength = knotIdxs.size();
		
		ArrayList<List<Double>> sorted_knot_sample = sort_sample_data(knotIdxs);
		
		List<Integer> optLeftIdxs = new ArrayList<Integer>();
		List<Integer> optRightIdxs = new ArrayList<Integer>();
		
		String selectedExplainingVariable = "";
		
		for(int j=0; j<n_explaining_variables; j++) {
			int n_sorted_elements = sorted_knot_sample.get(j).size();
			for(int k=0; k<n_sorted_elements; k++) {
				List<Integer> leftIdxs = new ArrayList<Integer>();
				List<Integer> rightIdxs = new ArrayList<Integer>();
				double t = sorted_knot_sample.get(j).get(k);
	            for(int i=0; i<knotSampleLength; i++) {	
	            	int sampleIdx = knotIdxs.get(i);
	            	if(explaining_variables[sampleIdx][j]<=t) {
	            		leftIdxs.add(sampleIdx);
	            	}else {
	            		rightIdxs.add(sampleIdx);
	            	}
	            }
	            
	            double leftCost = Double.MAX_VALUE;
	            double rightCost = Double.MAX_VALUE;
	            double cost = Double.MAX_VALUE;
	            
	            leftCost  = calc_cost(leftIdxs);
	            rightCost = calc_cost(rightIdxs);
	                     
	            cost = leftCost + rightCost;
	            
	            if(cost<minCost) {
	            	minCost = cost;
	            	minLeftCost = leftCost;
	            	minRightCost = rightCost;
	            	optLeftIdxs = leftIdxs;
	            	optRightIdxs = rightIdxs;
	            	selectedExplainingVariable = names_of_explaining_variables[j];      
	            	threshold = t;
	            }
	            	            
			}
			
		}
		
		overallCost = calc_cost(knotIdxs);
		
		HashMap<String,List<String>> split_infos = new HashMap<String, List<String>>();
		
    	//Cost
    	List<String> info = new ArrayList<String>(3);
    	info.add(Double.toString(minCost));
    	info.add(Double.toString(minLeftCost));
    	info.add(Double.toString(minRightCost));
    	split_infos.put("Cost",info);
    	
    	//Cost reduction
    	info = new ArrayList<String>(1);
    	double leftShare = (double) optLeftIdxs.size()/ (double) knotSampleLength;
    	double rightShare = (double) optRightIdxs.size()/(double) knotSampleLength;
    	double costReduction = overallCost-(leftShare*minLeftCost + rightShare*minRightCost);
    	info.add(Double.toString(costReduction));
    	split_infos.put("CostReduction",info);
    	
    	info = new ArrayList<String>(1);
    	info.add(Double.toString(leftShare*minLeftCost));
    	split_infos.put("CostReductionLeftKnot",info);
    	
    	info = new ArrayList<String>(1);
    	info.add(Double.toString(rightShare*minRightCost));
    	split_infos.put("CostReductionRightKnot",info);
    	
    	//Threshold
    	info = new ArrayList<String>(1);
    	info.add(Double.toString(threshold));
    	split_infos.put("Threshold",info);
    	
    	//SplittingFeature
    	info = new ArrayList<String>(1);
    	info.add(selectedExplainingVariable);
    	split_infos.put("SplittingFeature",info);
    	
    	//IndicesOfElementsInLeftKnot
    	int nIdxs = optLeftIdxs.size();
    	info = new ArrayList<String>(nIdxs);
    	for(int i=0; i<nIdxs; i++) {
    		info.add(Integer.toString(optLeftIdxs.get(i)));
    	}
    	split_infos.put("IndicesLeftKnot",info);
    	
    	//IndicesOfElementsInRightKnot
    	nIdxs = optRightIdxs.size();
    	info = new ArrayList<String>(nIdxs);
    	for(int i=0; i<nIdxs; i++) {
    		info.add(Integer.toString(optRightIdxs.get(i)));
    	}
    	split_infos.put("IndicesRightKnot",info);
    				
		return split_infos;
		
	}
	
	
	public static int count_classes_of_knot(double [] explained_vars_in_knot) {
		
		double [] unique_classes = Utilities.Utilities.get_unique_elements(explained_vars_in_knot);
		int n_classes = unique_classes.length;
		
		return n_classes;
		
	}
	
	
	public static List<String> count_elements_of_classes(List<Integer> idxs){
		
		int n = idxs.size();
		List<String> class_counts = new ArrayList<String>(nClasses);
		double [] knots_explained_vars = new double [n];
		
		for(int i=0; i<n; i++) {
			knots_explained_vars[i] = explained_variable[idxs.get(i)][0];
		}
		
		for(int c=0; c<nClasses; c++) {
			int [] idx = Utilities.Utilities.get_idx(knots_explained_vars, classes[c]);
			if(idx[0] == -1) {
				class_counts.add("0");
			}else {
				String count = Integer.toString(idx.length);
				class_counts.add(count);
			}
		}
		
		return class_counts;
		
	}
	
	
	public static List<String> count_elements_of_classes_with_str_idxs(List<String> idxs){
		
		int n = idxs.size();
		List<String> class_counts = new ArrayList<String>(nClasses);
		double [] knots_explained_vars = new double [n];
		
		for(int i=0; i<n; i++) {
			knots_explained_vars[i] = explained_variable[Integer.parseInt(idxs.get(i))][0];
		}
		
		for(int c=0; c<nClasses; c++) {
			int [] idx = Utilities.Utilities.get_idx(knots_explained_vars, classes[c]);
			if(idx[0] == -1) {
				class_counts.add("0");
			}else {
				String count = Integer.toString(idx.length);
				class_counts.add(count);
			}
		}
		
		return class_counts;
		
	}
	
	
	public static ArrayList<List<Double>> sort_sample_data(List<Integer> knotIdxs){
		
		int n_idxs = knotIdxs.size();
		
		double [][] knotSample = new double [n_idxs][n_explaining_variables];
		
		for(int i=0; i<n_explaining_variables; i++) {
			for(int j=0; j<n_idxs; j++) {
				int idx = knotIdxs.get(j);
				knotSample[j][i] = explaining_variables[idx][i];
			}
		}
		
		ArrayList<List<Double>> sorted_sample_data = new ArrayList<List<Double>>();
		
		for(int i=0; i<n_explaining_variables; i++) {
			List<Double> sorted_column = Utilities.Utilities.get_sorted_elements_of_matrix_column(knotSample, i);
			sorted_sample_data.add(sorted_column);
		}
		
		return sorted_sample_data;
		
	}
	
	
	public static double [][] calc_class_conditional_probs(List<Integer> idxs) {
		
		int n = idxs.size();
		
		double [][] probs = new double [nClasses][1];
		
		for(int c=0; c<nClasses; c++) {
			for(int i=0; i<n; i++) {
				if(explained_variable[idxs.get(i)][0] == classes[c]) {
					probs[c][0]++;
				}
			}
			probs[c][0] /= n;
		}
		
		return probs;
		
	}
	
	
	public static double calc_misclassification_rate(List<Integer> idxs) {
		
		double [][] probs = calc_class_conditional_probs(idxs);
		
		int maxIdx = 0;
	    double max = probs[0][0];
				
		for(int c=1; c<nClasses; c++) {
			if(probs[c][0]>max) {
				max = probs[c][0];
				maxIdx = c;
			}
		}
		
		double maxClass = classes[maxIdx];
		
		int n = idxs.size();
		double misclassification_rate = 0;
		
		for(int i=0; i<n; i++) {
			if(explained_variable[idxs.get(i)][0] == maxClass) {
				misclassification_rate++;
			}
		}
		
		misclassification_rate /=n;
		
		return misclassification_rate;
		
	}
	
	
	public static double calc_entropy(List<Integer> idxs) {
		
		double [][] class_cond_probs = calc_class_conditional_probs(idxs);
		
		double entropy = 0;
		
		for(int c=0; c<nClasses; c++) {
			entropy -= class_cond_probs[c][0]*Math.log(class_cond_probs[c][0]);
		}
		
		return entropy;
		
	}
	
	
	public static double calc_gini_index(List<Integer> idxs) {
		
		double [][] class_cond_probs = calc_class_conditional_probs(idxs);
		
		double gini = 0;
		
		for(int c=0; c<nClasses; c++) {
			gini += 1.0-Math.pow(class_cond_probs[c][0],2.0);
		}
		
		return gini;
		
	}
	
	
	public static double calc_regression_cost(List<Integer> idxs) {
		
		int n_idxs     = idxs.size();
		double mean    = 0.0;
		double regCost = 0.0;
		
		for(int i=0; i<n_idxs; i++) {
			mean += explained_variable[idxs.get(i)][0];
		}
		
		mean /= n_idxs;
		
		for(int i=0; i<n_idxs; i++) {
			double absDev = explained_variable[idxs.get(i)][0]-mean;
			regCost += Math.pow(absDev, 2.0);
		}
		
		return regCost;
		
	}
	
	
	public static double calc_cost(List<Integer> idxs) {
		
		double cost = 0.0;
		
		if(costMeasure == "Gini") {
			cost = calc_gini_index(idxs);
		}
		if(costMeasure == "Entropy") {
			cost = calc_entropy(idxs);
		}
		if(costMeasure == "MisclassificationRate") {
			cost = calc_misclassification_rate(idxs);
		}
		if(costMeasure == "Regression") {
			cost = calc_regression_cost(idxs);
		}
		
		return cost;
		
	}
		
	
	public static String [] get_valid_cost_measures() {
		
		String [] valid_measures = new String [4];
		
		valid_measures[0] = "Gini";
		valid_measures[0] = "Entropy";
		valid_measures[0] = "MisclassificationRate";
		valid_measures[0] = "Regression";
		
		return valid_measures;
	}
	
	
	public static HashMap <String, HashMap<String, HashMap<String,List<String>>>>  get_tree(){
		
		return tree;
		
	}
	
	public static HashMap<String, List<String>> get_knot_infos(int layerNumber, int knotNumber){
		
		if(layerNumber <= 0 || layerNumber > treeDepth) {
			throw new RuntimeException("Invalid layer number supplied.");
		}
		
		String layerLabel = "Layer"+ layerNumber;
		
		int nKnots = tree.get(layerLabel).size();
		
		if(knotNumber <= 0 || knotNumber > nKnots) {
			throw new RuntimeException("Invalid knot number supplied.");
		}
		
		String knotLabel  = "Knot" + knotNumber;
		
		HashMap<String, List<String>> knot_infos = tree.get(layerLabel).get(knotLabel);
		
		return knot_infos;
		
	}
	
	
	public static String get_knot_splitFeature(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		String splitFeature = knot_info.get("SplittingFeature").get(0);
		
		return splitFeature;
	
	}
	
	
	public static double get_knot_threshold(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		double threshold = Double.parseDouble(knot_info.get("Threshold").get(0));
		
		return threshold;
		
	}
	
	
	public static double get_knot_costReduction(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		double costReduction = Double.parseDouble(knot_info.get("CostReduction").get(0));
		
		return costReduction;
		
	}
	
	
	public static int get_knot_number_of_elements(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		int n = Integer.parseInt(knot_info.get("NumberOfElements").get(0));
		
		return n;
		
	}
	
	
	public static int [][] get_knot_class_distribution(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		List<String> distList = knot_info.get("nElementsOfClasses");
		
		int [][] classDist = new int [nClasses][1];
		
		for(int i=0; i<nClasses; i++) {
			classDist[i][0] = Integer.parseInt(distList.get(i));
		}
		
		return classDist;
		
	}
	
	
	public static double get_knot_cost(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		double cost = Double.parseDouble(knot_info.get("Cost").get(0));
		
		return cost;
		
	}
	
	
	public static List<Integer> get_knot_idxs(int layerNumber, int knotNumber){
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		List<String> strIdxs = knot_info.get("IdxsOfElements");
		int n = strIdxs.size();
		
		List<Integer> idxs = new ArrayList<Integer>(n);
		
		for(int i=0; i<n; i++) {
			idxs.add(Integer.parseInt(strIdxs.get(i)));
		}
		
		return idxs;
		
	}
	
	
	public static double [][] get_knot_explained_variables(int layerNumber, int knotNumber){
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		List<String> strIdxs = knot_info.get("IdxsOfElements");
		int n = strIdxs.size();
		
		double [][] explained_vars = new double [n][1];
		
		for(int i=0; i<n; i++) {
			int idx = Integer.parseInt(strIdxs.get(i));
			explained_vars[i][0] = explained_variable[idx][0];
		}
		
		return explained_vars;
		
	}
	
	
	public static double [][] get_knot_explaining_variables(int layerNumber, int knotNumber){
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		List<String> strIdxs = knot_info.get("IdxsOfElements");
		int n = strIdxs.size();
		
		double [][] explaining_vars = new double [n][n_explaining_variables];
		
		for(int i=0; i<n; i++) {
			int idx = Integer.parseInt(strIdxs.get(i));
			for(int j=0; j<n_explaining_variables; j++) {
				explaining_vars[i][j] = explaining_variables[idx][j];
			}
			
		}
		
		return explaining_variables;
		
	}
	
	
	public static String get_parent_knot(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		String parentKnot = knot_info.get("ParentKnot").get(0);
		
		return parentKnot;
		
	}

	
	public static void read_CART_categorical_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception{
		
		inputData = new InputDataManager();		
		inputData.fileReader(fileName, false, hasRowNames, hasColNames);
		categorical_explained_var = true;
		
	}
	
	
	public static void read_CART_numerical_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception{
		
		inputData = new InputDataManager();	
		inputData.fileReader(fileName, true, hasRowNames, hasColNames);
		categorical_explained_var = false;
		
	}
	
	
	public static void set_CART_explained_variable(String var_name) {
		name_of_explained_variable = var_name;
		if(inputData != null) {
			String [] colNames = inputData.colnames;
			int [] idx = Utilities.Utilities.get_idx(colNames, var_name);
			if(idx[0] == -1) {
				throw new RuntimeException(var_name + " not found in loaded input data.");
			}
			int n_vars = colNames.length;
			names_of_explaining_variables = new String [(n_vars-1)];
			for(int i=0; i<n_vars; i++) {
				if(colNames[i].contentEquals(var_name) == false) {
					names_of_explaining_variables[i] = colNames[i];
				}
			}
		}
	}
	
	
	public static void set_CART_explaining_variables(String [] var_names) {
		names_of_explaining_variables = var_names;
	}
	
	
	public static void set_CART_row_secection(String [] rowIdxs) {
		selectedRows = rowIdxs;
	}
	
	
	public static void set_CART_sampleName(String name) {
		inputData.setSampleName(name);
	}
	
	
	public static void set_CART_inputData() {
		
		n_explaining_variables = names_of_explaining_variables.length;
		
		String [] colNames = new String [1+n_explaining_variables];
		colNames[0] = name_of_explained_variable;
		for(int i=0; i<n_explaining_variables; i++) {
			colNames[i+1] = names_of_explaining_variables[i];
		}
		
		String [] rowNames = new String [n_observations];
		
		if(selectedRows != null) {
			n_observations = selectedRows.length;
			rowNames = selectedRows;
		}else {
			n_observations = inputData.numberOfRows-1;
			rowNames = new String [n_observations];
	    	for(int i=0; i<n_observations; i++){
	    		rowNames[i] = Integer.toString(i+1);
	    	}
		}
			
		explained_variable = new double [n_observations][1];
		explaining_variables = new double [n_observations][n_explaining_variables];
		
		if(categorical_explained_var == true) {
			
			String [] explained_var = new String [1];
			explained_var[0] = colNames[0];
			
			inputData.selectLoadedData(rowNames, explained_var);
			 
			classNames = Utilities.Utilities.get_unique_elements(inputData.selectedStrFileData);
			int nClasses = classNames.length;
			
			classes = new double [nClasses];
			
			for(int i=0; i<nClasses; i++) {
				classes[i] = i;
			}
					
			inputData.selectLoadedData(rowNames, colNames);			
			for(int i=0; i<n_observations; i++) {
				int [] idx = Utilities.Utilities.get_idx(classNames, inputData.selectedStrFileData[i][0]);
				explained_variable[i][0] = classes[idx[0]];
				for(int j=0; j<n_explaining_variables; j++) {
					explaining_variables[i][j] = Double.parseDouble(inputData.selectedStrFileData[i][j+1]);
				}
			}

		}else {
			inputData.selectLoadedData(rowNames, colNames);
			for(int i=0; i<n_observations; i++) {
				explained_variable[i][0] = inputData.selectedDblFileData[i][0];
				for(int j=0; j<n_explaining_variables; j++) {
					explaining_variables[i][j] = inputData.selectedDblFileData[i][j+1];
				}
			}					
		}
			
	}
	
	
	public static int getTreeDepth() {
		return treeDepth;
	}
	
	
	public static boolean isCategoricalTree() {
		return categorical_explained_var;
	}
	
	
	public static HashMap<String,List<String>> get_leafs_4_layers(){
		
		HashMap<String,List<String>> leafs = new HashMap<String,List<String>>();
		
		List<String> layerLeafs = new ArrayList<String>();
		
		if(treeDepth >= 2) {
			for(int i=1; i<treeDepth; i++) {
				
				int layerNumber = i+1;
				String layerLabel = "Layer"+layerNumber;
				layerLeafs = new ArrayList<String>();
				
				int prevLayerNumber = i;
				String prevLayerLabel = "Layer"+prevLayerNumber;
				HashMap<String, HashMap<String,List<String>>> layer = tree.get(prevLayerLabel);
				int n_prevKnots = layer.size();
				String [] prevLayerKnots = new String [n_prevKnots];
				for(int k=0; k<n_prevKnots; k++) {
					prevLayerKnots[k] = "Knot"+(k+1);
				}
									
				layer = tree.get(layerLabel);
				int nKnots = layer.size()/2;
				String [] parentKnots = new String [nKnots];
				int s = 0;
				for(int k=0; k<nKnots; k++) {
					parentKnots[k] = get_parent_knot(layerNumber, s+1);
					s += 2;
				}
				
				for(int k=0; k<n_prevKnots; k++) {
					String prevLayerKnot = prevLayerKnots[k];
					int [] idx = Utilities.Utilities.get_idx(parentKnots, prevLayerKnot);
					if(idx[0] == -1) {
						layerLeafs.add(prevLayerKnot);
					}
				}
				
				//if(layerLeafs.isEmpty() == false) {
				leafs.put(prevLayerLabel, layerLeafs);
				//}
																		
			}
		}
		
		//Leafs of the last layer
		String layerLabel = "Layer"+treeDepth;
		HashMap<String, HashMap<String,List<String>>> layer = tree.get(layerLabel);
		layerLeafs = new ArrayList<String>();
		int nKnots = layer.size();
		
		for(int k=0; k<nKnots; k++) {
			layerLeafs.add("Knot"+(k+1));
		}
		
		leafs.put(layerLabel, layerLeafs);
		
		return leafs;
		
	}
	
	
	public static HashMap<String, Integer> get_distOfSplittingFeatures() {
		
		HashMap<String, Integer> distOfSplitFeat = new HashMap<String, Integer>(n_explaining_variables);
		
		for(int i=0; i<n_explaining_variables; i++) {
			distOfSplitFeat.put(names_of_explaining_variables[i], 0);
		}
		
		for(int l=1; l<treeDepth; l++) {
			int layerNumber = (l+1);
			String layerLabel = "Layer"+layerNumber;
			int nKnots = tree.get(layerLabel).size();
			nKnots /=2;
			int s=0;
			for(int k=0; k<nKnots; k++) {
				int knotNumber = s+1;
				String splitFeat = get_knot_splitFeature(layerNumber,knotNumber);
				s += 2;
				int newNumber = distOfSplitFeat.get(splitFeat)+1;
				distOfSplitFeat.put(splitFeat, newNumber);
			}
		}
		
		return distOfSplitFeat;
		
	}
	
	
	public static void set_maxTreeDepth(int maxDepth) {
		maxTreeDepth = maxDepth;
	}
	
	
	public static void set_minNumberOfElementsInKnot(int minNumber) {
		minNumberOfElementsInKnot = minNumber;
	}
	
	
	public static HashMap<String, HashMap<String,Integer>> get_colPosOfLeafsInRegressionMatrix() {
		
		HashMap<String,List<String>> leafs = get_leafs_4_layers();
		int nLayers = leafs.size();		
		HashMap<String, HashMap<String,Integer>> colPosInRegMatrix = new HashMap<String, HashMap<String,Integer>>();
	    int idx = 0;
		for(int l=0; l<nLayers; l++) {
			String layerLabel = "Layer"+(l+1);
			int nLeafsInLayer = leafs.get(layerLabel).size();
			if(nLeafsInLayer != 0) {
				HashMap<String,Integer> colPos = new HashMap<String,Integer>();
				for(int i=0; i<nLeafsInLayer; i++) {
					String knotLabel = leafs.get(layerLabel).get(i);
					colPos.put(knotLabel, idx);					
					idx++;
				}
				colPosInRegMatrix.put(layerLabel, colPos);
			}
		}
		
		return colPosInRegMatrix;
		
	}
	
	
	public static void calc_regressor_matrix() {
		
		HashMap<String, HashMap<String,Integer>> colPosInRegMatrix = get_colPosOfLeafsInRegressionMatrix();
		
		int nLayersWithLeafs = colPosInRegMatrix.size();
		Object [] keys = colPosInRegMatrix.keySet().toArray();
		
		int nLeafs = 0;
		
		for(int l=0; l<nLayersWithLeafs; l++) {
			String layerLabel = (String) keys[l];
			nLeafs += colPosInRegMatrix.get(layerLabel).size();
		}
		
		regressor_matrix = new double [n_observations][nLeafs];
		
		for(int i=0; i<n_observations; i++) {

			for(int l=0; l<nLayersWithLeafs; l++) {
				String layerLabel = (String) keys[l];
				int layerNumber = Integer.parseInt(layerLabel.substring(5));
				Object [] leafNames = colPosInRegMatrix.get(layerLabel).keySet().toArray();	
				int nLeafNames = leafNames.length;
				for(int k=0; k<nLeafNames; k++) {
					String leafName = (String) leafNames[k];
					int knotNumber = Integer.parseInt(leafName.substring(4));
					List<Integer> leafIdxs = get_knot_idxs(layerNumber,knotNumber);
					int [] idxs = Utilities.Utilities.get_idx(leafIdxs, i);
					if(idxs[0] != -1) {
						int colIdx = colPosInRegMatrix.get(layerLabel).get(leafName);
						regressor_matrix[i][colIdx] = 1.0;
					}
				}
			}
					
		}
	
	}
	
	
	public static void calc_least_square_regressor_weights() {
		
		calc_regressor_matrix();
		
		if(regressor_matrix == null) {
			throw new RuntimeException("Regressor matrix for tree is not calculated yet.");
		}
			
    	obj_linearReg = new LinearRegression(explained_variable, regressor_matrix, useConstant);
    	
    	obj_linearReg.do_parameter_estimation();
		
    	System.out.println("");
    	
	}
	
	
	public static double [][] get_explained_variables() {
		return explained_variable;
	}
	
	
	public static String get_name_of_explained_variable() {
		return name_of_explained_variable;
	}
	
	
	public static String [] get_names_of_explaining_variables() {
		return names_of_explaining_variables;
	}
	
	
	public static LinearRegression get_linearRegObject() {
		if(obj_linearReg == null) {
			System.out.println("Linear regression not done yet for fitted tree.");
		}
		return obj_linearReg;
	}
	
	
	public static void set_dummy_tree() {
		
		HashMap <String, HashMap<String, HashMap<String,List<String>>>> dummyTree = new HashMap <String, HashMap<String, HashMap<String,List<String>>>>();
		
		//Layer no. 1 (Root Layer):
		HashMap<String, HashMap<String,List<String>>> layer = new HashMap<String, HashMap<String,List<String>>>();
		HashMap<String, List<String>> knot = get_knot_structure();
		knot.get("ParentKnot").add("Knot0");
		knot.get("Cost").add("100");
		knot.get("NumberOfElements").add("150");
		knot.get("nElementsOfClasses").add("50");
		knot.get("nElementsOfClasses").add("50");
		knot.get("nElementsOfClasses").add("50");
		layer.put("Knot1", knot);
		
		dummyTree.put("Layer1", layer);	
		
		//Layer no. 2 (Root Layer):
		layer = new HashMap<String, HashMap<String,List<String>>>();
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    knot.get("Cost").add("70");	    
	    knot.get("NumberOfElements").add("50");
	    knot.get("nElementsOfClasses").add("15");
	    knot.get("nElementsOfClasses").add("30");
	    knot.get("nElementsOfClasses").add("5");
	    knot.get("SplittingFeature").add("X1");
	    knot.get("Threshold").add("7.0");
	    knot.get("CostReduction").add("20.0");
	    layer.put("Knot1", knot);
		
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    knot.get("Cost").add("30");	    
	    knot.get("NumberOfElements").add("100");
	    knot.get("nElementsOfClasses").add("35");
	    knot.get("nElementsOfClasses").add("20");
	    knot.get("nElementsOfClasses").add("45");
	    knot.get("SplittingFeature").add("X1");
	    knot.get("Threshold").add("7.0");
	    knot.get("CostReduction").add("40.2");
	    layer.put("Knot2", knot);
	    
	    dummyTree.put("Layer2", layer);	
		
		//Layer no. 3 (Root Layer):
		layer = new HashMap<String, HashMap<String,List<String>>>();
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    layer.put("Knot1", knot);
		
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    layer.put("Knot2", knot);
	    
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot2");
	    layer.put("Knot3", knot);
	    
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot2");
	    layer.put("Knot4", knot);
	    
	    dummyTree.put("Layer3", layer);	
	    
		//Layer no. 4 (Root Layer):
		layer = new HashMap<String, HashMap<String,List<String>>>();
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    layer.put("Knot1", knot);
		
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    layer.put("Knot2", knot);
	    
	    knot = get_knot_structure(); 
	    knot.get("ParentKnot").add("Knot4");
	    knot.get("Cost").add("70");	    
	    knot.get("NumberOfElements").add("50");
	    knot.get("nElementsOfClasses").add("15");
	    knot.get("nElementsOfClasses").add("30");
	    knot.get("nElementsOfClasses").add("5");
	    knot.get("SplittingFeature").add("X1");
	    knot.get("Threshold").add("7.0");
	    knot.get("CostReduction").add("20.0");
	    layer.put("Knot3", knot);
	    	    	    
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot4");
	    knot.get("Cost").add("30");	    
	    knot.get("NumberOfElements").add("100");
	    knot.get("nElementsOfClasses").add("35");
	    knot.get("nElementsOfClasses").add("20");
	    knot.get("nElementsOfClasses").add("45");
	    knot.get("SplittingFeature").add("X1");
	    knot.get("Threshold").add("7.0");
	    knot.get("CostReduction").add("40.2");
	    layer.put("Knot4", knot);
	    
	    dummyTree.put("Layer4", layer);
	    
		//Layer no. 5 (Root Layer):
		layer = new HashMap<String, HashMap<String,List<String>>>();
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    layer.put("Knot1", knot);
		
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot1");
	    layer.put("Knot2", knot);
	    
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot2");
	    layer.put("Knot3", knot);
		
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot2");
	    layer.put("Knot4", knot);
	    
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot4");
	    layer.put("Knot5", knot);
		
	    knot = get_knot_structure();
	    knot.get("ParentKnot").add("Knot4");
	    knot.get("Cost").add("30");	    
	    knot.get("NumberOfElements").add("100");
	    knot.get("nElementsOfClasses").add("35");
	    knot.get("nElementsOfClasses").add("20");
	    knot.get("nElementsOfClasses").add("45");
	    knot.get("SplittingFeature").add("X1");
	    knot.get("Threshold").add("7.0");
	    knot.get("CostReduction").add("40.2");
	    layer.put("Knot6", knot);
	    
	    dummyTree.put("Layer5", layer);
	    
	    nClasses = 3;
		tree = dummyTree;
		treeDepth = tree.size();
				
	}

	
	public static void showLeafs(boolean show) {
		showLeafs = show;
	}
	
	
	@SuppressWarnings("static-access")
	public static CART clone_CART_obj_4_plots() {
				   
		CART cart_obj = new CART();
		cart_obj.tree                       = tree;
		cart_obj.classes                    = classes;
		cart_obj.classNames                 = classNames;
		cart_obj.nClasses                   = nClasses;
		cart_obj.treeDepth                  = treeDepth;
		cart_obj.explaining_variables       = explaining_variables;
		cart_obj.n_explaining_variables     = n_explaining_variables;
		cart_obj.explained_variable         = explained_variable;
		cart_obj.name_of_explained_variable = name_of_explained_variable;
		cart_obj.categorical_explained_var  = categorical_explained_var;
		cart_obj.obj_linearReg              = obj_linearReg;
		
		return cart_obj;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_tree() {
		
		//TODO: Delete this test method!
		//set_dummy_tree();
		   
		CART cart_obj = clone_CART_obj_4_plots();
		
		DecisionTreeGraphics obj_graph = new DecisionTreeGraphics();
		
		obj_graph.setCARTObject(cart_obj);
		obj_graph.setGraphWidth(800);
		obj_graph.setGraphHeight(600);
		obj_graph.showLeafs(showLeafs);
		String [] title    = {"Decision Tree"};
		obj_graph.setTitle(title, null, "15");
		if(inputData.sampleName != null) {
			String [] subTitle = {inputData.sampleName};
			obj_graph.setSubTitle1(subTitle, null, "12");
		}
		obj_graph.plotLayerGrid = true;
		obj_graph.setDecimalPlaces4InfoBox(2);
		obj_graph.plotDecisionTree();
        	
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotFittedExplainedVariable() {
		
		CART cart_obj = clone_CART_obj_4_plots();
		
		DecisionTreeGraphics obj_graph = new DecisionTreeGraphics();

		obj_graph.setCARTObject(cart_obj);
		obj_graph.setGraphWidth(600);
		obj_graph.setGraphHeight(600);
		
		obj_graph.plotFittedExplainedVariable();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotHistOfResiduals() {
		
		CART cart_obj = clone_CART_obj_4_plots();
		
		DecisionTreeGraphics obj_graph = new DecisionTreeGraphics();

		obj_graph.setCARTObject(cart_obj);
		obj_graph.setGraphWidth(600);
		obj_graph.setGraphHeight(600);
		
		obj_graph.plotHistOfResiduals();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotLeafWeights() {
		
		CART cart_obj = clone_CART_obj_4_plots();
		
		DecisionTreeGraphics obj_graph = new DecisionTreeGraphics();

		obj_graph.setCARTObject(cart_obj);
		obj_graph.setGraphWidth(600);
		obj_graph.setGraphHeight(600);
		
		obj_graph.plotLeafWeights();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_DistOfSplittingFeature() {
   
		CART cart_obj = clone_CART_obj_4_plots();
		
		DecisionTreeGraphics obj_graph = new DecisionTreeGraphics();
		
		obj_graph.setCARTObject(cart_obj);
		obj_graph.setGraphWidth(600);
		obj_graph.setGraphHeight(300);

		if(inputData.sampleName != null) {
			String [] subTitle = {inputData.sampleName};
			obj_graph.setSubTitle1(subTitle, null, "12");
		}
		
		obj_graph.plotLayerGrid = true;
		obj_graph.plotDistOfSplittingFeature();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
		
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/DecisionTrees/BostonHousePriceData.txt";
		
		CART obj_cart = new CART();
		
		obj_cart.read_CART_numerical_input_data(dirName, true, true);
		obj_cart.set_CART_explained_variable("MEDV");
		obj_cart.set_CART_sampleName("Boston Housing Data");
		obj_cart.set_CART_inputData();
		obj_cart.set_minNumberOfElementsInKnot(2);
		obj_cart.set_maxTreeDepth(10);
		obj_cart.fit_tree();
		
		calc_least_square_regressor_weights();
		
		
		//obj_cart.showLeafs(true);
		//obj_cart.plot_tree();
		//obj_cart.plot_DistOfSplittingFeature();
		
		obj_cart.plotFittedExplainedVariable();
		//obj_cart.plotHistOfResiduals();
		//obj_cart.plotLeafWeights();
		
		System.out.println("Finished.");
		
	}
	
}
