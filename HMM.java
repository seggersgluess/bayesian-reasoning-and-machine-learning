package HiddenMarkovModels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

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
	
	//n_states x n_observations
	static double [][] obs_model_probs;
	
	//Results of Viterbi algorithm
	//Most probable state sequence z_1, z_2, ..., z_T
	static int [][] stateSequence;
	static double [][] stateProbs;
	
	//Stats of the (Baum-Welch) EM algorithm
	static int max_iterations;
	static double convergence_criterion;

	static int n_iterations;
	static boolean convergence_reached = false;
	
	static double log_likelihood;
	
	//Results of the (Baum-Welch) EM algorithm
	static double [][] initProbs;
	static double [][] transMatrix;
	
	//Gaussian observation model
	static ArrayList<List<Double>> my;
	static ArrayList<List<Double>> Sigma;
	
	static List<String> probs4Plot;
	
	
	static void calcHMM(boolean calcViterbiSequence){
		
		expectation_maximization();
		
		if(calcViterbiSequence == true){
			runViterbiAlgorithm();
		}
		
		if(probs4Plot.size() == 0){
			plot_probs();
		}
		
		
	}
	
	
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
	
	
	public static double [][] get_filtered_probs(int observationNumber){
		
		if(observationNumber>n_observations-1){
			throw new RuntimeException(observationNumber + " is not in observed data sample.");
		}
	
		double [][] probs = new double [n_states][1];
		
		for(int i=0; i<n_states; i++){
			probs[i][0] = filteredProbs.get(observationNumber).get(i);
		}
	
		return probs;
		
	}
	
	
	public static double [][] get_filtered_probs(){

		double [][] probs = new double [n_observations][n_states];
		
		for(int t=0; t<n_observations; t++){
			for(int i=0; i<n_states; i++){
				probs[i][0] = filteredProbs.get(t).get(i);
			}
		}

		return probs;
		
	}
	
	
	public static double [][] get_smoothed_probs(){

		double [][] probs = new double [n_observations][n_states];
		
		for(int t=0; t<n_observations; t++){
			for(int i=0; i<n_states; i++){
				probs[i][0] = smoothedProbs.get(t).get(i);
			}
		}

		return probs;
		
	}
	
	
	public static double [][] get_Viterbi_probs(){
		
		return stateProbs;
		
	}
		
		
	public static void expectation_maximization(){
		
		obs_model_probs = new double [n_states][n_observations];
		
		double newLogLik  = 0.0;
		double prevLogLik = Double.MIN_VALUE;		
		double [][] n_1  = new double [n_states][1];
		double [][] n_j  = new double [n_states][1];
		double [][] margSliceSum = new double [n_states][n_states];
		double [][] my_exp       = new double [n_variables][1];
		double [][] sigma_exp    = new double [n_variables][n_variables];
		double [][] logLik_term3  = new double [n_states][1];
		
		my    = new ArrayList<List<Double>>(n_states);
		Sigma = new ArrayList<List<Double>>(n_states);
		
		for(int i=0; i<max_iterations; i++){
			
			newLogLik  = 0.0;
			
			ForwardBackwardAlgorithm fb = new ForwardBackwardAlgorithm();
			
			fb.runForwardAlgorithm();
			fb.runBackwardAlgorithm();
			
			for(int s=0; s<n_states; s++){
				
				my.add(s,new ArrayList<Double>(n_variables));
				Sigma.add(s,new ArrayList<Double>((int) Math.pow(n_variables,2.0)));
				
				n_1[s][0] = smoothedProbs.get(0).get(s);
				
				for(int t=1; t<n_observations; t++){
												
					double [][] x_t = MatrixOperations.get_row_vec_from_matrix(observed_variables, t);
					double [][] x_t_prod = MatrixOperations.multiplication(x_t, MatrixOperations.transpose(x_t));
					
					margSliceSum = MatrixOperations.add(margSliceSum,get_two_slice_marg_matrix(t));
					n_j[s][0]    = n_j[s][0] + smoothedProbs.get(t).get(s);
					my_exp       = MatrixOperations.add(my_exp,MatrixOperations.scalar_multiplication(smoothedProbs.get(t).get(s),x_t));
					sigma_exp    = MatrixOperations.add(sigma_exp,MatrixOperations.scalar_multiplication(smoothedProbs.get(t).get(s),x_t_prod));
								
					double [][] obs_model_probs = fb.get_obs_model_probs();
					
					logLik_term3[s][0] = logLik_term3[s][0]+smoothedProbs.get(t).get(s)*obs_model_probs[s][t];
						
				}
				
				int idx = 0;
				
				for(int j=0; j<n_variables; j++){
					my.get(s).add(j,my_exp[j][0]);
					for(int k=0; k<n_variables; k++){
						Sigma.get(s).add(idx,sigma_exp[k][j]);
						idx++;
					}						
				}
				
			}
				
			//E-Step
			for(int s=0; s<n_states; s++){
				
				double logLik_term2 = 0.0;
				
				for(int k=0; k<n_states; k++){
					logLik_term2 = margSliceSum[s][k]*Math.log(transMatrix[s][k]);
				}
				
				newLogLik = newLogLik + n_1[s][0]*Math.log(initProbs[s][0]) + logLik_term2 + logLik_term3[s][0];
			}
			
			//M-Step
			initProbs = n_1;
			
			for(int s=0; s<n_states; s++){
				
				List<Double> myUpdate    = new ArrayList<Double>(n_variables);
				List<Double> sigmaUpdate = new ArrayList<Double>((int) Math.pow(n_variables, 2.0));
				double sum = 0.0;
				
				for(int k=0; k<n_states; k++){
					sum = sum + margSliceSum[s][k];
				}
				
				for(int k=0; k<n_states; k++){
					transMatrix[s][k] = margSliceSum[s][k]/sum;	
				}
				
				double [][] myUpdateVec = new double [n_variables][0];
							
				for(int j=0; j<n_variables; j++){
					myUpdate.add(j,my.get(s).get(j)/n_j[s][0]);
					myUpdateVec[j][0] = myUpdate.get(j);
				}
				
				double [][] myUpdateProd = MatrixOperations.multiplication(myUpdateVec, MatrixOperations.transpose(myUpdateVec));
				myUpdateProd = MatrixOperations.vecAs2dArray(myUpdateProd);
				
				int idx = 0;
				
				for(int j=0; j<n_variables; j++){					
					for(int k=0; k<n_variables; k++){
						sigmaUpdate.add(idx,Sigma.get(s).get(idx)/n_j[s][0]-myUpdateProd[idx][0]);
						idx++;						
					}					
				}
				
				my.add(s,myUpdate);
				Sigma.add(s,sigmaUpdate);
				
			}

			if(Math.abs(prevLogLik-newLogLik)<=convergence_criterion){
				System.out.println("EM algorithm for HMM has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;					
				break;
			}
			
			prevLogLik = newLogLik;
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = newLogLik;
				
	}
	
	
	public static void runViterbiAlgorithm(){
		
		Viterbi v = new Viterbi();
		v.runViterbi();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_probs(){
		
		int nGraphs = probs4Plot.size();
		
		if(nGraphs == 0){
			return;
		}
		
		GenGraphics graphs  = new GenGraphics();
		
	 	graphs.setNumberOfPlotColums(nGraphs);
	 	graphs.setNumberOfPlotRows(1);
	 	
	 	graphs.setGraphWidth(500);
	 	graphs.setGraphHeight(300);
		
		double [][] data    = new double [n_observations][1];
		double [][] xValues = new double [n_observations][1];
		int probIdx = 0;
		
		String [] yLabels = new String [nGraphs];
		String [] xLabels = new String [nGraphs];
		
		List<Color> lineColor = new ArrayList<Color>();
		
		for(int i=0; i<n_observations; i++){
			xValues[i][0] = i+1;
		}
			
		for(int i=0; i<nGraphs; i++){
			
			if(probs4Plot.get(i)=="filtered"){
				for(int j=0; j<n_observations; j++){
					data[j][0] = filteredProbs.get(j).get(probIdx);
				}	
				
				yLabels[i] = "Filtered Prob. (State " + probIdx + ")";
			}
			
			if(probs4Plot.get(i)=="smoothed"){
				for(int j=0; j<n_observations; j++){
					data[j][0] = smoothedProbs.get(j).get(probIdx);
				}	
				
				yLabels[i] = "Smoothed Prob. (State " + probIdx + ")";
			}
			
			if(probs4Plot.get(i)=="Viterbi"){
				for(int j=0; j<n_observations; j++){
					data[j][0] = stateProbs[j][probIdx];
				}	
				
				yLabels[i] = "Viterbi Prob. (State " + probIdx + ")";
			}
			
			graphs.plotLines(xValues,data,true);
			
			lineColor.add(Color.BLACK);
			
		}
	 	      		
	 	String [] title    = get_prob_types_4_plot();
	
	 	graphs.setLineColor(lineColor);	 	
 	 	graphs.setTitle(title, "bold", "12");
 	 	graphs.setYLabel(yLabels, null, "10");
 	 	graphs.setXLabel(xLabels, null, "10");
 	 	graphs.setNumberOfDigits4XAxis(0);   
 	 	graphs.setNumberOfDigits4YAxis(2);
 	 	graphs.setFontOfXAxisUnits("plain", 10);
 	 	graphs.setFontOfYAxisUnits("plain", 10);
		
 	 	graphs.plot();
	 	
	}
	
	
	public static void set_probs4plot(String [] probTypes){
		
		int n_types = probTypes.length;
		probs4Plot  = new ArrayList<String>(n_types);
		String [] allowed_types = get_allowed_prob_types();
		
		for(int i=0; i<n_types; i++){
			int [] idx = Utilities.get_idx(allowed_types, probTypes[i]);
			if(idx[0] == -1){
				throw new RuntimeException(probTypes[i] + " are not valid probabilities for plot.");
			}
			
			probs4Plot.add(probTypes[i]);
				
		}
				
	}
	
	
	public static String [] get_allowed_prob_types(){
	
		String [] allowed_types = new String [3];
		
		allowed_types[0] = "filtered";
		allowed_types[0] = "smoothed";
		allowed_types[0] = "Viterbi";
		
		return allowed_types;
		
	}
	
	
	public static String [] get_prob_types_4_plot(){
		
		int nTypes = probs4Plot.size();
		
		String [] types = new String [nTypes];
		
		for(int i=0; i<nTypes; i++){
			types[i] = probs4Plot.get(i);
		}
		
		return types;
		
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
