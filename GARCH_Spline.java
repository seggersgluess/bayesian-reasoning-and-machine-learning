package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.SimulatedAnnealing;

public class GARCH_Spline extends GARCH{

	static int n_spline_knots;
	static int [] data_points_4_spline;
	static double [][] spline_pars;
	static double [][] spline_pars_limits;
	
	public GARCH_Spline(double[][] obs_variables, int start_idx, int end_idx, int obsLag, int volaLag, int resLag, int number_of_spline_knots) {
		super(obs_variables, start_idx, end_idx, obsLag, volaLag, resLag);
		n_spline_knots = number_of_spline_knots;
		calc_equi_dist_data_points_4_GARCH();
	}

	
	public static double [][] calc_volatilies_from_GARCH_Spline(){
		
		double [][] h_t_GARCH  = calc_volatilies_from_GARCH();
		
		double [][] h_t = new double [n_usedObservations][1];
		double [][] exp_spline = new double [n_usedObservations][1];
		
		for(int t=0; t<n_usedObservations; t++){
			for(int i=0; i<n_spline_knots; i++){
				exp_spline[t][0] += spline_pars[i][0]*Math.pow(((t+1)-data_points_4_spline[i]), 2.0);
			}
			h_t[t][0] = Math.exp(exp_spline[t][0])*h_t_GARCH[t][0];
		}
		
		return h_t;
		
	}
	
	
	public static double [][] get_exp_spline_4_GARCH(){
		
		double [][] exp_spline = new double [n_usedObservations][1];
		
		for(int t=0; t<n_usedObservations; t++){
			for(int i=0; i<n_spline_knots; i++){
				exp_spline[t][0] += spline_pars[i][0]*Math.pow(((t+1)-data_points_4_spline[i]), 2.0);
			}
			exp_spline[t][0] = Math.exp(exp_spline[t][0]);
		}
		
		return exp_spline;
		
	}
	
	
	public static void calc_equi_dist_data_points_4_GARCH(){
		
		int delta_t = (int) Math.ceil((double) n_usedObservations/(double) n_spline_knots);
		
		data_points_4_spline = new int [n_spline_knots];
		
		for(int i=0; i<(n_spline_knots-1); i++){
			data_points_4_spline[i] = delta_t*(i+1);		
		}
		
		data_points_4_spline[(n_spline_knots-1)] = n_usedObservations;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH_Spline(){
		
    	double [] start_value = get_start_values_4_est_GARCH_Spline();
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    		
    		ArrayList<List<Double>> limits = get_par_limits_4_MLE(start_value);
    		int n_pars = limits.get(0).size();        	
        	double [] lower_values = new double [n_pars];
        	double [] upper_values = new double [n_pars];
        	int n_GARCH_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
        	int idx =0;
        	
    		for(int i=0; i<n_pars; i++){
    			lower_values[i] = limits.get(0).get(i);
    		    upper_values[i] = limits.get(1).get(i);
    		    if(i>=n_GARCH_pars){
    		    	if(spline_pars_limits != null){   		    		
    		    		lower_values[i] = spline_pars_limits[idx][0];
    		    		upper_values[i] = spline_pars_limits[idx][1];
    		    		idx++;
    		    	}else{
    		    		lower_values[i] = -1e-08;
        		    	upper_values[i] = 1e-08;	
    		    	}   		    	
    		    }
    		}
    		
        	if(optimizer == "DEoptim"){
        		DifferentialEvolution optim = new DifferentialEvolution(GARCH_Spline::opti_log_likelihood_4_GARCH_Spline, 200);
        		optim.set_convergence_criterion(1e-02);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values, lower_values);
            	set_GARCH_Spline_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		SimulatedAnnealing optim = new SimulatedAnnealing(GARCH_Spline::opti_log_likelihood_4_GARCH_Spline, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values, lower_values);
            	set_GARCH_Spline_pars_from_vec(optim.get_optimal_candidate());        	
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    	  			
    	boolean truePars = check_GARCH_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("Spline-GARCH restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, 0, 0);

    	volatilities = calc_volatilies_from_GARCH_Spline();
    	fittedValues = calc_est_values_from_GARCH();
    	residuals = MatrixOperations.substract(obs_data, fittedValues);
    		
	}
	
	
	public static double opti_log_likelihood_4_GARCH_Spline(double [] pars, double [] further_ars){
		
		set_GARCH_Spline_pars_from_vec(pars);
		
		double logLik = calc_log_likelihood_4_GARCH_Spline();
		
		return logLik;
		
	}
	
	
	public static double calc_log_likelihood_4_GARCH_Spline(){
		
		boolean truePars = check_GARCH_par_restrictions();
		
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH_Spline();	
		double [][] est_vars = calc_est_values_from_GARCH();
		
		for(int t=0; t<n_usedObservations; t++){
			
			if(h_t[t][0]<=0.0){
				h_t[t][0] = 1e-10;
			}
			
			logLik += -0.5*(Math.log(h_t[t][0]) + Math.pow((observed_variables[startIdx+t][0]-est_vars[t][0]), 2.0)/h_t[t][0]);
		}
		
		logLik += (-1.0)*n_usedObservations/2.0*Math.log(2.0*Math.PI);
			
		if(Double.isInfinite(logLik) == true){
			logLik = -1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			logLik = -1e+100;
		}
		
		System.out.println(-logLik);
		
		return -logLik;
		
	}
	
	
	public static void set_GARCH_Spline_pars_from_vec(double [] par_vec){
		
		arPars      = new double [lag4observedVariables+1][1];
		volaPars    = new double [lag4volatility+1][1];
		maPars      = new double [lag4residuals][1];
		spline_pars = new double [n_spline_knots][1];
		
		int idx = 0;
		
		for(int i=0; i<(lag4observedVariables+1); i++){
			arPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<(lag4volatility+1); i++){
			volaPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<lag4residuals; i++){
			maPars[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<n_spline_knots; i++){
			spline_pars[i][0] = par_vec[idx];
			idx++;
		}
			
	}
	
	
	public static void set_spline_pars_limits(double [][] limits){
		
		if(limits.length != n_spline_knots){
			throw new RuntimeException("Invalid suppied limits for spline parameters.");
		}
		
		spline_pars_limits = limits;
		
	}
	
	
	public static double [] get_start_values_4_est_GARCH_Spline(){
		
		//Note: Start values of spline pars are set 0.0		
		int n_GARCH_pars = lag4observedVariables+lag4volatility+lag4residuals+2;
		int n_pars = n_GARCH_pars+n_spline_knots;
		
		double [] start_values = new double [n_pars];		
		double [] start_values_4_GARCH = get_start_values_4_est_GARCH();
		
		for(int i=0; i<n_GARCH_pars; i++){
			start_values[i] = start_values_4_GARCH[i];
		}
		
		return start_values;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    	
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"France"};
    	
		InputDataManager inputData = new InputDataManager();		
		inputData.fileReader(file, true, true, true);
    	
    	int nData = inputData.numberOfRows-1;
    	String [] rownames = new String [nData];
    	for(int i=0; i<nData;i++){
    		rownames[i] = Integer.toString(i+1);
    	}
		
		inputData.selectLoadedData(rownames, colnames);
		double [][] obsData = inputData.selectedDblFileData;

		int obsLag  = 1;
		int volaLag = 1;
		int maLag   = 1;
		int n_knots = 5;
		int start_idx = obsLag+maLag;
		int end_idx = obsData.length;
		
		double [][] limits = new double [n_knots][2];
		for(int i=0; i<n_knots; i++){
			limits[i][0] = -1e-06;
			limits[i][1] = 1e-06;
		}
		
		GARCH_Spline obj_arch = new GARCH_Spline(obsData, start_idx, end_idx, obsLag, volaLag, maLag, n_knots);
		obj_arch.set_spline_pars_limits(limits);
		obj_arch.do_MLE_4_GARCH_Spline();
		
		MatrixOperations.print_matrix(arPars);
		MatrixOperations.print_matrix(volaPars);
		MatrixOperations.print_matrix(maPars);
		
		System.out.println("");
		MatrixOperations.print_matrix(spline_pars);
		
		System.out.println(logLikelihood);
		
		plot_GARCH_estimates();
		
	}
	
}
