package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import Optimization.L_BFGS;

public class ARCH {

	static double [][] observed_variables;
	
	static int n_observations;
	static int n_usedObservations;
	static int n_variables; 
	
	static int startIdx;
	static int endIdx;
	
	static int lag4observedVariables;
	static int lag4heteroscedasticity;
	
	static double [][] arPars;
	static double [][] heteroPars;
	
	static ArrayList<List<Double>> lagged_variables;
	
	//constructor
	public ARCH(double [][] obs_variables, int start_idx, int end_idx, int obsLag, int heteroscedasticityLag){
		
		if(start_idx < 0){
			throw new RuntimeException("No valid start index supplied.");
		}
		
		if(end_idx < 0){
			throw new RuntimeException("No valid end index supplied.");
		}
		
		if(start_idx >= end_idx){
			throw new RuntimeException("No valid start and end indices supplied.");
		}
		
		if(obsLag < 0){
			throw new RuntimeException("No valid number of lags supplied.");
		}
		
		if(heteroscedasticityLag < 0){
			throw new RuntimeException("No valid number of lags for heteroscedasticity term supplied.");
		}
		
		observed_variables = obs_variables;
		n_observations = observed_variables.length;
		n_variables    = observed_variables[0].length;
		
		startIdx               = start_idx;
		endIdx                 = end_idx-1;
		lag4observedVariables  = obsLag;
		lag4heteroscedasticity = heteroscedasticityLag;
		
		n_usedObservations = endIdx-startIdx+1;
			
		lagged_variables = new ArrayList<List<Double>>((lag4heteroscedasticity+1));
		
		for(int m=0; m<lagged_variables.size(); m++){			
			lagged_variables.add(get_lagged_Y(m));			
		}
		
		
		//Check input consistency!!!
		
	}
	
	
	public static double calc_log_likelihood_4_ARCH(){
		
		double logLik = 0.0;
		
		double [][] h_t = new double [n_usedObservations][1];
		
		double [][] conv_lagged_vars = get_lagged_Y_for_lag(lag4heteroscedasticity);
		
		for(int m=1; m<lag4heteroscedasticity; m++){
			double [][] higher_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(m), n_usedObservations, (lag4observedVariables+1));
			
			for(int r=0; r<n_usedObservations; r++){
				for(int c=0; c<(lag4observedVariables+1); c++){
					h_t[r][0] += higher_lagged_vars[r][c]*arPars[c][0];
				}
				h_t[r][0] = conv_lagged_vars[r][m]-h_t[r][0];
				h_t[r][0] = Math.pow(h_t[r][0],2.0);
			}
			
		}
		
		conv_lagged_vars = MatrixOperations.get_matrix_from_vec(lagged_variables.get(0), n_usedObservations, (lag4observedVariables+1));
		double [][] est_vars = MatrixOperations.multiplication(conv_lagged_vars, arPars);
		
		for(int t=0; t<n_usedObservations; t++){
			logLik += 0.5*(Math.log(h_t[t][0]) - Math.pow(observed_variables[startIdx+t][0]-est_vars[t][0], 2.0)/h_t[t][0]);
		}
		
		logLik += (-1.0)*n_usedObservations/2.0*Math.log(2.0*Math.PI);
		
		return logLik;
		
	}
	
	
	public static double opti_log_likelihood_4_ARCH(double [] pars, double [] further_ars){
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			arPars[i][0] = pars[i];
		}
		
		for(int i=0; i<(lag4heteroscedasticity+1); i++){
			heteroPars[i][0] = pars[(lag4observedVariables+1+i)];
		}
		
		double logLik = calc_log_likelihood_4_ARCH();
		
		return logLik;
		
	}
	
	
	public static void do_MLE_4_ARCH(){
		
    	double [] start_value = get_start_values_4_est_ARCH();
    	
    	L_BFGS optim = new L_BFGS(ARCH::opti_log_likelihood_4_ARCH, 100000);    	
    	optim.do_LBFGS_Optimization(start_value);
		
	}
	
	
	public static double [] get_start_values_4_est_ARCH(){
		
		int n_pars = lag4observedVariables + lag4heteroscedasticity + 2;
		
		double [] start_values = new double [n_pars];
		
		//Not implemented yet!
		
		//-> Use OLS for conv. AR estimation
		//-> Use OLS for squared lagged residuals
		
		return start_values;
		
	}
	
	
	//sets TxP matrix X_bar = [Y_t-1-m,...,Y_t-p-m] w.r.t. additional lag m beside lag p
	public static List<Double> get_lagged_Y(int m){
		
		List<Double> lagged_Y = new ArrayList<Double>((lag4observedVariables+1));
		
		int T = n_usedObservations;
		
		for(int t=0; t<T; t++){
			lagged_Y.add(1.0);
			int curIdx = startIdx+t;
			for(int p=0; p<lag4observedVariables; p++){
				int lagIdx = curIdx-p-1-m;
				lagged_Y.add(observed_variables[lagIdx][0]);			
			}
		}
			
		return lagged_Y;
		
	}
	
	
	//sets TxP matrix X_bar = [Y_t-1,...,Y_t-p] for specific lag p
	public static double [][] get_lagged_Y_for_lag(int lag){
		
		int T = n_usedObservations;
		
		double [][] lagged_Y = new double [T][(lag4observedVariables+1)];
		
		for(int t=0; t<T; t++){
			lagged_Y[t][0] = 1.0;
			int curIdx = startIdx+t;
			for(int p=0; p<lag; p++){
				int lagIdx = curIdx-p-1;
				lagged_Y[t][p] = observed_variables[lagIdx][0];			
			}
		}
			
		return lagged_Y;
		
	}
	

	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"Germany"};
    	int numberOfLags = 1;
    	
		InputDataManager inputData = new InputDataManager();		
		inputData.fileReader(file, true, true, true);
    	
    	int nData = inputData.numberOfRows-1;
    	String [] rownames = new String [nData];
    	for(int i=0; i<nData;i++){
    		rownames[i] = Integer.toString(i+1);
    	}
		
		inputData.selectLoadedData(rownames, colnames);
		double [][] obsData = inputData.selectedDblFileData;

		
		
		
	}
	
	
}
