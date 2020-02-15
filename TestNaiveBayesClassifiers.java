package NaiveBayesClassifier;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import TextMining.TextAnalyzer;

public class TestNaiveBayesClassifiers {

	//Test 1: Full data for training from loaded txt-file (Iris-Data)
	public static void test1GaussianNB() {
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/NaiveBayesClassifiers/IrisData.txt";
		String yName   = "species";
		NBC nbc = new NBC();
		nbc.readAndSetInputParameters(dirName,yName);
		nbc.fitNBC("GAUSSIAN");
		nbc.inSamplePredictionNB();
		
		double accuracy = nbc.get_accuracyRate();
		
		System.out.println("---- Test results for Gaussian NBC ----");
		System.out.println("In sample fit of Gaussian NBC with accuracy rate: " + accuracy*100.0 + "%");
		
		double []   testData = {5.7, 2.9, 4.2, 1.3};
		double [][] testVec  = MatrixOperations.transpose(MatrixOperations.convArrayToVec(testData));
		double [][] prediction = nbc.predictNB(testVec);
		
		System.out.println("Prediction from test vector: class=" + prediction[0][0] + " (Valid result: class=1.0)");
		
	}
	
	
	//Test 2: Sample split into training and test sub-samples (Iris-Data)
	public static void test2GaussianNB() {
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/NaiveBayesClassifiers/IrisData.txt";
		String yName   = "species";
		
		HashMap<String, double [][]> inputData = readInput4Test(dirName, yName);
		double [][] y = inputData.get("y");
		double [][] X = inputData.get("X");
		
		int n = y.length;
		int n_test = n/2;
		
		int [] trainingIdxs = Utilities.Utilities.getRandomIntNumbers(0, n-1, n_test);
		
		double [][] y_training = MatrixOperations.get_sub_matrix_4_row_idxs(y, trainingIdxs);
		double [][] X_training = MatrixOperations.get_sub_matrix_4_row_idxs(X, trainingIdxs);
		
		int [] testIdxs = Utilities.Utilities.getRandomIntNumbers(0, n-1, n_test);
		
		double [][] y_test = MatrixOperations.get_sub_matrix_4_row_idxs(y, testIdxs);
		double [][] X_test = MatrixOperations.get_sub_matrix_4_row_idxs(X, testIdxs);
		
		NBC nbc = new NBC();
		nbc.fitNBC("GAUSSIAN", y_training, X_training);
		nbc.inSamplePredictionNB();
		
		double accuracy = nbc.get_accuracyRate();
		
		System.out.println("---- Test results for Gaussian NBC ----");
		System.out.println("In sample fit of Gaussian NBC with accuracy rate: " + accuracy*100.0 + "%");
		
		double [][] prediction = nbc.predictNB(X_test);
		
		accuracy = nbc.calcAccuracyRates(prediction, y_test);
		System.out.println("Out of sample fit of Gaussian NBC with accuracy rate: " + accuracy*100.0 + "%");
		
	}
	
	
	public static void testCategoricalNB() {
		
	}
	
	
	public static void testBernoulliNB() {
		
		HashMap<Integer, HashMap<String, Integer>> fedPaperStats = readAndAnalyzeFederalPapers4MultinomialNBC();
		
		HashMap<String, List<Integer>> authorIdxs = getFedPaperIdxs4Authors(fedPaperStats);
		
		int nHamilton = authorIdxs.get("hamilton").size();
		int nMadison  = authorIdxs.get("madison").size();
		
		System.out.println("Papers by Hamilton:");
		System.out.println(authorIdxs.get("hamilton"));
		System.out.println("Papers by Madison:");
		System.out.println(authorIdxs.get("madison"));
		
		int nPapers4Training = nHamilton+nMadison;
		
		String [] keyWords = getKeyWords();
		int nKeyWords = keyWords.length;
		
		//X binary --> occurs word x=1 else x=0
		double [][] X = new double [nPapers4Training][nKeyWords];
		double [][] y = new double [nPapers4Training][1];
		
		for(int i=0; i<nHamilton; i++) {
			int idx = authorIdxs.get("hamilton").get(i);
			for(int w=0; w<nKeyWords; w++) {
				int wordCount = fedPaperStats.get(idx).get(keyWords[w]);
				if(wordCount > 0) {
					X[i][w] = 1.0;
				}
				y[i][0]  = 0.0;
			}	
		}
		
		for(int i=0; i<nMadison; i++) {
			int idx = authorIdxs.get("madison").get(i);
			for(int w=0; w<nKeyWords; w++) {
				int wordCount = fedPaperStats.get(idx).get(keyWords[w]);
				if(wordCount > 0) {
					X[i+nHamilton][w] = 1.0;
				}			 
				y[i+nHamilton][0]  = 1.0;
			}	
		}

		NBC nbc = new NBC();
		nbc.fitNBC("BERNOULLI", y, X);
		nbc.inSamplePredictionNB();
		 
		double accuracy = nbc.get_accuracyRate();
		
		System.out.println("---- Test results for Bernoulli NBC classification of Federalist Papers----");
		System.out.println("In sample fit of Bernoulli NBC with accuracy rate: " + accuracy*100.0 + "%");
		
		System.out.println("");
		System.out.println("In sample predicted probabilities of Bernoulli NBC:");
		MatrixOperations.print_matrix(nbc.get_inSamplePredProbability());
		
		int nUnknown = authorIdxs.get("unknown").size();
		
		//X binary --> occurs word x=1 else x=0
		double [][] X_unknown = new double [nUnknown][nKeyWords];
				
		for(int i=0; i<nUnknown; i++) {
			int idx = authorIdxs.get("unknown").get(i);
			for(int w=0; w<nKeyWords; w++) {
				int wordCount = fedPaperStats.get(idx).get(keyWords[w]);
				if(wordCount > 0) {
					X_unknown[i][w] = 1.0;
				}
			}	
		}
		
		HashMap<String, double[][]> pred_unknown = nbc.predictNB_value_and_probability(X_unknown);
		System.out.println("");
		System.out.println("Prediction of unknown Federalist Papers:");
		MatrixOperations.print_matrix(pred_unknown.get("Prediction"));
		MatrixOperations.print_matrix(pred_unknown.get("Probability"));
		
	}

	
	public static void test1MultinomialNB() {
		
		HashMap<Integer, HashMap<String, Integer>> fedPaperStats = readAndAnalyzeFederalPapers4MultinomialNBC();
		
		HashMap<String, List<Integer>> authorIdxs = getFedPaperIdxs4Authors(fedPaperStats);
		
		int nHamilton = authorIdxs.get("hamilton").size();
		int nMadison  = authorIdxs.get("madison").size();
		
		System.out.println("Papers by Hamilton:");
		System.out.println(authorIdxs.get("hamilton"));
		System.out.println("Papers by Madison:");
		System.out.println(authorIdxs.get("madison"));
		
		int nPapers4Training = nHamilton+nMadison;
		
		String [] keyWords = getKeyWords();
		int nKeyWords = keyWords.length;
		
		double [][] X = new double [nPapers4Training][nKeyWords];
		double [][] y = new double [nPapers4Training][1];
		
		for(int i=0; i<nHamilton; i++) {
			int idx = authorIdxs.get("hamilton").get(i);
			for(int w=0; w<nKeyWords; w++) {
				X[i][w] = fedPaperStats.get(idx).get(keyWords[w]);
				y[i][0]  = 0.0;
			}	
		}
		
		for(int i=0; i<nMadison; i++) {
			int idx = authorIdxs.get("madison").get(i);
			for(int w=0; w<nKeyWords; w++) {
				X[i+nHamilton][w] = fedPaperStats.get(idx).get(keyWords[w]);
				y[i+nHamilton][0]  = 1.0;
			}	
		}

		NBC nbc = new NBC();
		nbc.fitNBC("MULTINOMIAL", y, X);
		nbc.inSamplePredictionNB();
		 
		double accuracy = nbc.get_accuracyRate();
		
		System.out.println("---- Test results for Multinomial NBC classification of Federalist Papers ----");
		System.out.println("In sample fit of Multinomial NBC with accuracy rate: " + accuracy*100.0 + "%");
		
		int nUnknown = authorIdxs.get("unknown").size();
		
		double [][] X_unknown = new double [nUnknown][nKeyWords];
				
		for(int i=0; i<nUnknown; i++) {
			int idx = authorIdxs.get("unknown").get(i);
			for(int w=0; w<nKeyWords; w++) {
				X_unknown[i][w] = fedPaperStats.get(idx).get(keyWords[w]);
			}	
		}
		
		HashMap<String, double[][]> pred_unknown = nbc.predictNB_value_and_probability(X_unknown);
		System.out.println("");
		System.out.println("Prediction of unknown Federalist Papers:");
		MatrixOperations.print_matrix(pred_unknown.get("Prediction"));
		MatrixOperations.print_matrix(pred_unknown.get("Probability"));
		
	}
	
	
	//Reads 85 Federalist Papers from file and counts 10 key words
	public static HashMap<Integer, HashMap<String, Integer>> readAndAnalyzeFederalPapers4MultinomialNBC() {
		
		int nFedPapers = 85;
		
		HashMap<Integer, HashMap<String, Integer>> res = new HashMap<Integer, HashMap<String, Integer>>();
		
		String [] keyWords = getKeyWords();
		int nKeyWords = keyWords.length;
		
		String [] words = new String [nKeyWords+3]; 
		//Set first three words for identifying author
		words[0] = "hamilton";
		words[1] = "madison";
		words[2] = "jay";
		
		for(int i=0; i<nKeyWords; i++) {
			words[i+3] = keyWords[i];
		}
		
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/NaiveBayesClassifiers/FederalistPapers/";
    	
		for(int i=0; i<nFedPapers; i++) {
			
			String dirName4Paper = dirName + "Federalist_No" + Integer.toString(i+1) + ".txt";
			
			TextAnalyzer textAnalyzer = new TextAnalyzer();
			textAnalyzer.loadAndPrepareText(dirName4Paper);
			HashMap<String, Integer> wc = textAnalyzer.getWordCount(words);
			res.put((i+1), wc);
		}
		
		return res;	
	}
	
	
	public static String [] getKeyWords() {	
		String [] keyWords = {"although", "always", "commonly", "consequently", "considerable",
                              "enough", "there", "upon", "while", "whilst"};		
		return keyWords;
	}
	
	
	public static HashMap<String, List<Integer>> getFedPaperIdxs4Authors(HashMap<Integer, HashMap<String, Integer>> fedPaperStats) {
		
		HashMap<String, List<Integer>> authorIdxs = new HashMap<String, List<Integer>>(4);
		authorIdxs.put("madison", new ArrayList<Integer>());
		authorIdxs.put("hamilton", new ArrayList<Integer>());
		authorIdxs.put("jay", new ArrayList<Integer>());
		authorIdxs.put("unknown", new ArrayList<Integer>());
		
		int nFedPaper = fedPaperStats.size();
		
		for(int i=0; i<nFedPaper; i++) {
			int paperNo = i+1;
			int hamilton = fedPaperStats.get(paperNo).get("hamilton");
			int madison = fedPaperStats.get(paperNo).get("madison");
			int jay = fedPaperStats.get(paperNo).get("jay");
			if(madison != 0 && hamilton == 0){
				authorIdxs.get("madison").add(paperNo);
			}
			if(madison == 0 && hamilton != 0){
				authorIdxs.get("hamilton").add(paperNo);
			}
			if(madison != 0 && hamilton != 0){
				authorIdxs.get("unknown").add(paperNo);
			}
			if(jay !=0) {
				authorIdxs.get("jay").add(paperNo);
			}
		}
		
		return authorIdxs;
	}
		
	
	
	public static HashMap<String, double [][]> readInput4Test(String dirName, String name_explained_variable) {
		
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
		
		int n_explaining_variables = colNames.length-1;
		int n_observations = inputData.numberOfRows-1;
		
		String [] names_of_explaining_variables = new String [n_explaining_variables];
		
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
		
		double [][] explained_variable = new double [n_observations][1];
		double [][] explaining_variables = new double [n_observations][n_explaining_variables];
		
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
		
		HashMap<String, double [][]> loadedInput = new HashMap<String, double [][]>(2);
		loadedInput.put("y", explained_variable);
		loadedInput.put("X", explaining_variables);
		
		return loadedInput;
	}
	
	
	public static void main(String[] args) {
		//test2GaussianNB();
		//test1MultinomialNB();
		testBernoulliNB();
	}
	
}
