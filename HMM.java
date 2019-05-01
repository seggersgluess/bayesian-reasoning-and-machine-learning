package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;

public class HMM {

	static InputDataManager inputData;
	
	static double [][] observed_variables;
	
	static int n_observations;
	static int n_variables;

	static int n_states;
	
	//Results of forward-backward algorithm
	static ArrayList<List<Double>> filteredProbs;
	static ArrayList<List<Double>> smoothedProbs;
	static ArrayList<List<Double>> twoSliceMarginals;
	
	//Results of Viterbi algorithm
	static ArrayList<List<Double>> stateSequence;
	
	//Stats of the (Baum-Welch) EM algorithm
	static int max_iterations;
	
	static int n_iterations;
	static boolean convergence_reached = false;
	
	//Results of the (Baum-Welch) EM algorithm
	static double [][] initProbs;
	static double [][] transMatrix;
	
	//Gaussian observation model
	static ArrayList<List<Double>> my;
	static ArrayList<List<Double>> Sigma;
	
	
	static double [][] get_my(int stateNumber){
		
		if(stateNumber>n_states-1){
			throw new RuntimeException(stateNumber + " is not a valid number of states.");
		}
		
		double [][] myVec = new double [n_variables][1];
	
		for(int i=0; i<n_variables; i++){
			myVec[i][0] = my.get(stateNumber).get(i);
		}
		
		return myVec;
		
	}
	
	
	static double [][] get_sigma_matrix(int stateNumber){
		
		if(stateNumber>n_states-1){
			throw new RuntimeException(stateNumber + " is not a valid number of states.");
		}
		
		double [][] sigmaMatrix = new double [n_variables][n_variables];
		
		int idx = 0;
		
		for(int i=0; i<n_variables; i++){
			for(int j=0; j<n_variables; j++){
				sigmaMatrix[j][i] = my.get(stateNumber).get(idx);
				idx++;
			}
			
		}
		
		return sigmaMatrix;
		
	}
	
	
	static double [][] get_two_slice_marg_matrix(int observationNumber){
		
		if(observationNumber>n_observations-1){
			throw new RuntimeException(observationNumber + " is not in observed data sample.");
		}
		
		double [][] sliceMargMatrix = new double [n_states][n_states];
		
		int idx = 0;
		
		for(int i=0; i<n_states; i++){
			for(int j=0; j<n_states; j++){
				sliceMargMatrix[j][i] = twoSliceMarginals.get(observationNumber).get(idx);
				idx++;
			}
			
		}
		
		return sliceMargMatrix;
		
	}
	
	
	public static void expectation_maximization(){
		
		
		
	}
	
	
	public static void plot_probs(String [] probTypes){
		
		String [] allowed_types = get_allowed_prob_types();
		
	}
	
	
	public static String [] get_allowed_prob_types(){
	
		String [] allowed_types = new String [3];
		
		allowed_types[0] = "filtered";
		allowed_types[0] = "smoothed";
		allowed_types[0] = "Viterbi";
		
		return allowed_types;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void read_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception{
		
		inputData = new InputDataManager();
		
		inputData.fileReader(fileName, false, hasRowNames, hasColNames);
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void select_input_data(String [] rownames, String [] colnames, String ref_class){
		
		inputData.selectLoadedData(rownames, colnames);
		
		observed_variables = inputData.selectedDblFileData;		
		n_observations     = observed_variables.length;
		n_variables        = observed_variables[0].length;
		
	}
	
	
	public static void remove_input_data_manager(){	
		
		inputData = null;		
		
	}
	
	
}
