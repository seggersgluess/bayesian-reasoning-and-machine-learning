package AdaptiveBasisModels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import DataManagement.InputDataManager;

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
	static double minNumberOfElementsInKnot = 0.0;
	static int maxTreeDepth = 0;
	
	//Cost Measures
	static String costMeasure = "Gini";
	
	//Number of layers of the tree
	static int treeDepth;
	
	static HashMap <String, HashMap<String, HashMap<String,List<String>>>> tree;
	
	
	//Tree infos for knots
	//static ArrayList<List<Double>> Cost;
	//static ArrayList<List<Double>> CostReduction;
	//static ArrayList<List<String>> SplittingFeature;
	//static ArrayList<ArrayList<List<Integer>>> NumberOfElements4Class;
	//static ArrayList<ArrayList<List<Integer>>> IndicesOfElementsInKnot;
	
	
	
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
				String parentKnotLabel = "Knot" + k+1;
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
					
					String splittingFeature = splitInfos.get("SplittingFeature").get(0);
					tree.get(prevLayerLabel).get(parentKnotLabel).get("SplittingFeature").add(splittingFeature);
					
					HashMap<String,List<String>> leftKnot = get_knot_structure();
					HashMap<String,List<String>> rightKnot = get_knot_structure();
					
					//ParentKnot
					leftKnot.get("ParentKnot").add(parentKnotLabel);
					rightKnot.get("ParentKnot").add(parentKnotLabel);
					
					//Cost
					leftKnot.get("Cost").add(splitInfos.get("Cost").get(1));
					rightKnot.get("Cost").add(splitInfos.get("Cost").get(2));
					
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
				layerCount++;
	        }else {
	        	//Tree stops to grow (no further split)
	        	stop = true;
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
		
		List<String> emptyList = new ArrayList<String>();
		
		knot_struct.put("Cost", emptyList);
		knot_struct.put("SplittingFeature", emptyList);
		knot_struct.put("NumberOfElements", emptyList);
		knot_struct.put("nElementsOfClasses", emptyList);
		knot_struct.put("IdxsOfElements", emptyList);
		knot_struct.put("ParentKnot", emptyList);
		
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
		
		int knotSampleLength = knotIdxs.size();
		
		ArrayList<List<Double>> sorted_knot_sample = sort_sample_data(knotIdxs);
		
		List<Integer> optLeftIdxs = new ArrayList<Integer>();
		List<Integer> optRightIdxs = new ArrayList<Integer>();
		
		String selectedExplainingVariable = "";
		
		for(int j=0; j<n_explaining_variables; j++) {
			int n_sorted_elements = sorted_knot_sample.get(j).size();
			List<Integer> leftIdxs = new ArrayList<Integer>();
			List<Integer> rightIdxs = new ArrayList<Integer>();
			for(int k=0; k<n_sorted_elements; k++) {
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
    	double leftShare = optLeftIdxs.size()/knotSampleLength;
    	double rightShare = optRightIdxs.size()/knotSampleLength;
    	double costReduction = overallCost-(leftShare*minLeftCost + rightShare*minRightCost);
    	info.add(Double.toString(costReduction));
    	split_infos.put("CostReduction",info);
    	
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
	
	
	public static int get_knot_number_of_elements(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> knot_info = get_knot_infos(layerNumber, knotNumber);
		
		int splitFeature = Integer.parseInt(knot_info.get("NumberOfElements").get(0));
		
		return splitFeature;
		
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
	}
	
	
	public static void set_CART_explaining_variables(String [] var_names) {
		names_of_explaining_variables = var_names;
	}
	
	
	public static void set_CART_row_secection(String [] rowIdxs) {
		selectedRows = rowIdxs;
	}
	
	
	//TODO: Integrate this method into constructor of CART!
	public static void set_CART_inputData() {
		
		int n_explaining_variables = names_of_explaining_variables.length;
		
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
			
		inputData.selectLoadedData(rowNames, colNames);
		
		explained_variable = new double [n_observations][1];
		explaining_variables = new double [n_observations][n_explaining_variables];
		
		if(categorical_explained_var == true) {
			
			classNames = Utilities.Utilities.get_unique_elements(inputData.selectedStrFileData);
			int nClasses = classNames.length;
			
			classes = new double [nClasses];
			
			for(int i=0; i<nClasses; i++) {
				classes[i] = i;
			}
					
			for(int i=0; i<n_observations; i++) {
				int [] idx = Utilities.Utilities.get_idx(classNames, inputData.selectedStrFileData[i][0]);
				explained_variable[i][0] = classes[idx[0]];
				for(int j=0; j<n_explaining_variables; j++) {
					explaining_variables[i][j] = Double.parseDouble(inputData.selectedStrFileData[i][j+1]);
				}
			}

		}else {
			for(int i=0; i<n_observations; i++) {
				explained_variable[i][0] = inputData.selectedDblFileData[i][0];
				for(int j=0; j<n_explaining_variables; j++) {
					explaining_variables[i][j] = inputData.selectedDblFileData[i][j+1];
				}
			}					
		}
			
	}
	
	
	public static void plot_tree() {
		
		
		
	}
	
	
}
