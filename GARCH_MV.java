package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;
import Optimization.DifferentialEvolution;
import Optimization.SimulatedAnnealing;
import Regression.LinearRegression;
import Utilities.Utilities;

public class GARCH_MV extends GARCH{

	static int n_observedVariables;
	
	static ArrayList<List<Double>> listOfVolaPars;
	static ArrayList<List<Double>> listOfMaPars;
	
	static List<Integer> listOfLag4Volatility;
	static List<Integer> listOfLag4Residuals;
	
	static int lag4StandDev;
	static int lag4Corr;
	static double [][] alpha_corr;
	static double [][] beta_corr;
	static ArrayList<List<Double>> R_t_list;
	
	//TODO: super() related to the parent class above GARCH class! Now it´s too dirty.
	public GARCH_MV(double[][] obs_variables, int start_idx, int end_idx, int [] volaLag, int [] resLag, int sdLag, int corrLag) {
		super(obs_variables, start_idx, end_idx, 0, volaLag[0], resLag[0]);
		
		n_observedVariables = obs_variables[0].length;
		
		set_volatility_lag_2_list(volaLag);
		set_ma_lag_2_list(resLag);
		
		lag4StandDev = sdLag;
		lag4Corr     = corrLag;
		
		n_heteroPars = get_number_of_GARCH_MV_pars();
		n_arPars = 0;
		
	}

	
	public static void set_volatility_lag_2_list(int [] volaLag){
		
		if(volaLag.length != n_observedVariables){
			throw new RuntimeException("Invalid number of volatility lags supplied for GARCH-MV.");
		}
		
		listOfLag4Volatility = new ArrayList<Integer>(n_observedVariables);
		
		for(int i=0; i<n_observedVariables; i++){
			listOfLag4Volatility.add(volaLag[i]);
		}
		
	}
	
	
	public static void set_ma_lag_2_list(int [] maLag){
		
		if(maLag.length != n_observedVariables){
			throw new RuntimeException("Invalid number of moving average lags supplied for GARCH-MV.");
		}
		
		listOfLag4Residuals = new ArrayList<Integer>(n_observedVariables);
		
		for(int i=0; i<n_observedVariables; i++){
			listOfLag4Residuals.add(maLag[i]);
		}
		
	}
	
	
	public static double calc_quasi_log_likelihood_4_GARCH_MV_first_step(){
		
		boolean truePars = check_GARCH_MV_first_step_par_restrictions();
		 
		if(truePars == false){			
			return 1e+100;
		}
		
		double logLik = 0.0;
		
		double [][] h_t = calc_volatilies_from_GARCH_MV();		
		
		for(int i=0; i<n_observedVariables; i++){
			
			for(int t=0; t<n_usedObservations; t++){
				
				if(h_t[t][0]<=0.0){
					return 1e+100;
				}
				
				logLik += -0.5*(Math.log(h_t[t][i]) + Math.pow((observed_variables[startIdx+t][i]),2.0)/h_t[t][i]);
			}
			
			if(Double.isInfinite(logLik) == true){
				return 1e+100;
			}   
			
			if(Double.isNaN(logLik) == true){
				return 1e+100;
			}
			
		}
		
		logLik += (-1.0)*n_usedObservations*n_observedVariables/2.0*Math.log(2.0*Math.PI);
		
		System.out.println(-logLik);
		
		return -logLik;
			
	}
	
	
	public static double opti_log_likelihood_4_GARCH_MV_first_step(double [] pars, double [] further_ars){
		
		set_GARCH_MV_pars_of_first_step_from_par_vec(pars);
		
		double logLik = calc_quasi_log_likelihood_4_GARCH_MV_first_step();
		
		return logLik;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_MLE_4_GARCH_MV(){
		    	
    	double [] start_value = get_start_values_4_first_step_est_GARCH_MV();
 	
    	set_start_values_4_second_step_est_GARCH_MV();
    	
    	if(optimizer == "DEoptim" || optimizer == "SANN"){
    	
    		ArrayList<List<Double>> limits = get_par_limits_4_MLE(start_value);
    		int n_pars_first_step = limits.get(0).size();   
    		int n_pars_sec_step   = lag4StandDev + lag4Corr;
        	double [] lower_values_first_step = new double [n_pars_first_step];
        	double [] upper_values_first_step = new double [n_pars_first_step];
        	double [] lower_values_sec_step   = new double [n_pars_sec_step];
        	double [] upper_values_sec_step   = new double [n_pars_sec_step];
        	
    		for(int i=0; i<n_pars_first_step; i++){
    			lower_values_first_step[i] = limits.get(0).get(i);
    		    upper_values_first_step[i] = limits.get(1).get(i);
    		}
    		
    		for(int i=0; i<n_pars_sec_step; i++){
    			lower_values_sec_step[i] = -1.0;
    			upper_values_sec_step[i] = 1.0;
    		}
    		
        	if(optimizer == "DEoptim"){
        		//First step Q-MLE:
        		DifferentialEvolution optim = new DifferentialEvolution(GARCH_MV::opti_log_likelihood_4_GARCH_MV_first_step, 500);
        		optim.set_convergence_criterion(convergence_criterion);
            	optim.set_number_of_function_eval(10);
            	optim.do_Differential_Evolution_Optimization(upper_values_first_step, lower_values_first_step);
            	set_GARCH_MV_pars_of_first_step_from_par_vec(optim.get_optimal_candidate());
            	volatilities = calc_volatilies_from_GARCH_MV();
            	//Second step Q-MLE:
            	optim = new DifferentialEvolution(GARCH_MV::opti_log_likelihood_4_GARCH_MV_second_step, 500);
            	optim.do_Differential_Evolution_Optimization(upper_values_sec_step, lower_values_sec_step);
            	set_corr_R_t();
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}
    		
        	if(optimizer == "SANN"){
        		//First step Q-MLE:
        		SimulatedAnnealing optim = new SimulatedAnnealing(GARCH_MV::opti_log_likelihood_4_GARCH_MV_first_step, 10000);
            	optim.do_Simulated_Annealing_Optimization(upper_values_first_step, lower_values_first_step);
            	set_GARCH_MV_pars_of_first_step_from_par_vec(optim.get_optimal_candidate());   
            	volatilities = calc_volatilies_from_GARCH_MV();
            	//Second step Q-MLE:
            	optim = new SimulatedAnnealing(GARCH_MV::opti_log_likelihood_4_GARCH_MV_second_step, 500);
            	optim.do_Simulated_Annealing_Optimization(upper_values_sec_step, lower_values_sec_step);
            	set_corr_R_t();
            	logLikelihood = (-1.0)*optim.get_optimal_value();
        	}   
        	
    	}
    			
    	boolean truePars = check_GARCH_MV_first_step_par_restrictions();
    	
    	if(truePars == false){
    		System.out.println("GARCH-MV restrictions for heteroscedasticity parameters violated.");
    	}
    	
    	
    	
	}
	

	public static boolean check_GARCH_MV_first_step_par_restrictions(){
		
		boolean restrictions_satisfied = true;
		
		double summedPars = 0.0;
		
		for(int k=0; k<n_observedVariables; k++){
			
			summedPars = 0.0;
			
			int lag4Vola = listOfLag4Volatility.get(k);
			int lag4Ma   = listOfLag4Residuals.get(k);
			
			if(listOfVolaPars.get(k).get(0) <= 0.0){
				restrictions_satisfied = false;
			}
			
			for(int i=0; i<lag4Vola; i++){
				summedPars += listOfVolaPars.get(k).get(i+1);
				if(listOfVolaPars.get(k).get(i+1) < 0.0){
					restrictions_satisfied = false;
					break;
				}
			}
			
			for(int i=0; i<lag4Ma; i++){
				summedPars += listOfMaPars.get(k).get(i);
				if(listOfMaPars.get(k).get(i) < 0.0){
					restrictions_satisfied = false;
					break;
				}
			}
				
		}
		
		if(summedPars >= 1.0){
			restrictions_satisfied = false;
		}
		
		return restrictions_satisfied;
		
	}
	
	
	public static void set_GARCH_MV_pars_of_first_step_from_par_vec(double [] par_vec){
		
		listOfVolaPars = new ArrayList<List<Double>>(n_observedVariables);
		listOfMaPars = new ArrayList<List<Double>>(n_observedVariables);
		
		int idx = 0;
		
		for(int i=0; i<n_observedVariables; i++){
			
			List<Double> pars = new ArrayList<Double>();
			int nVolaPars = listOfLag4Volatility.get(i)+1;
		    int nMaPars   = listOfLag4Residuals.get(i);
		    
			for(int j=0; j<nVolaPars; j++){
				pars.add(par_vec[idx]);
				idx++;
			}
			
			listOfVolaPars.add(pars);
			
			pars = new ArrayList<Double>();
			
			for(int j=0; j<nMaPars; j++){
				pars.add(par_vec[idx]);
				idx++;
			}
			
			listOfMaPars.add(pars);
			
		}
			
	}
	
	
	public static double [][] calc_volatilies_from_GARCH_MV(){
		
		double [][] h_t = new double [n_usedObservations][n_observedVariables];
		
		for(int i=0; i<n_observedVariables; i++){
			
			double [][] obsData = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, 0, endIdx, i, i);
			
			for(int t=0; t<n_usedObservations; t++){
				h_t[t][i] = listOfVolaPars.get(i).get(0);
			}

			double h_t_prev = 0.0;
		
			for(int t=0; t<n_usedObservations; t++){
				h_t_prev += Math.pow(obsData[startIdx+t][0],2.0);
			}
	    
			h_t_prev /= n_usedObservations;
			
			double [][] conv_lagged_vars = get_lagged_Y_for_lag(obsData, startIdx, endIdx, listOfLag4Residuals.get(i));
					
			for(int m=0; m<listOfLag4Residuals.get(i); m++){
				for(int t=0; t<n_usedObservations; t++){
					double u_t = conv_lagged_vars[t][m+1];
					u_t = Math.pow(u_t,2.0);
					h_t[t][i] += listOfMaPars.get(i).get(m)*u_t;		
				}			
			}
			
			for(int t=0; t<n_usedObservations; t++){
				for(int p=0; p<listOfLag4Volatility.get(i); p++){
					int prevCounter = t-p-1;
					if(prevCounter<0){
						h_t[t][i] += listOfVolaPars.get(i).get(1+p)*h_t_prev;
					}else{
						h_t[t][i] += listOfVolaPars.get(i).get(1+p)*h_t[t-p-1][0];
					}
						
				}			
			}
			
		}
			
		return h_t;
		
	}
	
	
	public static double [][] calc_standardized_obs_variables(){
		
		double [][] standardized_obs_vars = new double [n_usedObservations][n_observedVariables];
		
		for(int t=0; t<n_usedObservations; t++){
			for(int i=0; i<n_observedVariables; i++){
				standardized_obs_vars[t][i] = observed_variables[startIdx+t][i]/volatilities[t][i];
			}
		}
		
		return standardized_obs_vars;
		
	}
	
	
	public static double [][] calc_unconditional_cov(double [][] standardized_obs_vars){
		
		double [][] uncond_cov = new double [n_observedVariables][n_observedVariables];
		
		for(int t=0; t<n_usedObservations; t++){
			for(int i=0; i<n_observedVariables; i++){
				for(int j=0; j<n_observedVariables; j++){
					uncond_cov[i][j] += standardized_obs_vars[t][i]*standardized_obs_vars[t][j];
				}
			}
		}
		
		uncond_cov = MatrixOperations.scalar_multiplication(1.0/n_usedObservations, uncond_cov);
		
		return uncond_cov;
		
	}
	
	
	public static ArrayList<List<Double>> calc_Q_4_corr(double [][] uncond_cov, double [][] standardized_obs_vars){
		
		ArrayList<List<Double>> Q_list = new ArrayList<List<Double>>(n_usedObservations); 
		
		double firstTerm = 1.0-sum_Q_pars();
		
		double [][] Q_prev = new double [n_observedVariables][n_observedVariables]; 
		
		for(int i=0; i<n_observedVariables; i++){
			for(int j=0; j<n_observedVariables; j++){
				Q_prev[i][j] += standardized_obs_vars[0][i]*standardized_obs_vars[0][j];
			}
		}
		
		double [][] Q = new double [n_observedVariables][n_observedVariables];
		double [][] Q_matrix= MatrixOperations.scalar_multiplication(firstTerm, uncond_cov);
		
		for(int t=0; t<n_usedObservations; t++){
			
			for(int i=0; i<lag4StandDev; i++){
				int sdLag = t-i-1;
				if(sdLag<0){
					for(int j=0; j<n_observedVariables; j++){
						for(int k=0; k<n_observedVariables; k++){
							Q_matrix[j][k] += alpha_corr[i][0]*Q_prev[j][k];
						}
					}
				}else{
					for(int j=0; j<n_observedVariables; j++){
						for(int k=0; k<n_observedVariables; k++){
							Q_matrix[j][k] += alpha_corr[i][0]*standardized_obs_vars[sdLag][j]*standardized_obs_vars[sdLag][k];
						}
					}
				}
			}
			
			for(int i=0; i<lag4Corr; i++){
				int corrLag = t-i-1;
				if(corrLag<0){
					for(int j=0; j<n_observedVariables; j++){
						for(int k=0; k<n_observedVariables; k++){
							Q_matrix[j][k] += beta_corr[i][0]*Q_prev[j][k];
						}
					}
				}else{
					Q = MatrixOperations.get_matrix_from_vec(Q_list.get(corrLag), n_observedVariables, n_observedVariables);
					for(int j=0; j<n_observedVariables; j++){
						for(int k=0; k<n_observedVariables; k++){
							Q_matrix[j][k] += beta_corr[i][0]*Q[j][k];
						}
					}
				}
			}
			
			Q_list.add(MatrixOperations.vecAsList(Q_matrix));
			
		}
			
		return Q_list;
		
	}
	
	
	public static void set_start_values_4_second_step_est_GARCH_MV(){
		
		alpha_corr = new double [lag4StandDev][1];
		beta_corr  = new double [lag4Corr][1];
		
		alpha_corr[0][0] = 0.5;
		beta_corr[0][0] = 0.5;
		
	}
	
	
	//public static boolean check_GARCH_MV_second_step_par_restrictions(){
		
		//TODO: Not implemented yet!
		
	//}
	
	
	public static double calc_quasi_log_likelihood_4_GARCH_MV_second_step(){
		
		//TODO: check corr pars!
		//boolean truePars = check_GARCH_MV_second_step_par_restrictions();
		 
		//if(truePars == false){			
		//	return 1e+100;
		//}
		
		double logLik = 0.0;
		
		double [][] stand_obs_vars = calc_standardized_obs_variables();
		double [][] uncond_cov = calc_unconditional_cov(stand_obs_vars);
		ArrayList<List<Double>> Q_list = calc_Q_4_corr(uncond_cov, stand_obs_vars);
		double [][] Q_star = new double [n_observedVariables][n_observedVariables];
		
		for(int t=0; t<n_usedObservations; t++){
				
			double [][] Q_t = MatrixOperations.get_matrix_from_vec(Q_list.get(t), n_observedVariables, n_observedVariables);
				 
			for(int i=0; i<n_observedVariables; i++){
				Q_star[i][i] = 1.0/Math.sqrt(Q_t[i][i]);			
			}
			
			double [][] R_t = MatrixOperations.multiplication(MatrixOperations.multiplication(Q_star, Q_t),Q_star);
			
			//TODO: Error catching -> try inverse & determinant!
			double [][] R_t_inv = MatrixOperations.inverse(R_t);
			double R_t_det = MatrixOperations.determinant(R_t);
			R_t_det = Math.log(R_t_det);
			
			//TODO: Faster?
			double [][] stand_obs_vars_t = MatrixOperations.get_row_vec_from_matrix(stand_obs_vars, t);
						
			logLik += R_t_det+MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(stand_obs_vars_t), R_t_inv), stand_obs_vars_t)[0][0];
		
		}
		
		logLik/=-2.0;
		
		if(Double.isInfinite(logLik) == true){
			return 1e+100;
		}   
		
		if(Double.isNaN(logLik) == true){
			return 1e+100;
		}
		
		System.out.println(-logLik);
		
		return -logLik;
			
	}
	
	
	//public static double calc_total_log_likelihood_4_GARCH_MV(){
		
		
		
	//}
	
	
	public static double opti_log_likelihood_4_GARCH_MV_second_step(double [] pars, double [] further_ars){
		
		set_GARCH_MV_pars_of_second_step_from_par_vec(pars);
		
		double logLik = calc_quasi_log_likelihood_4_GARCH_MV_second_step();
		
		return logLik;
		
	}
	
	
	public static void set_corr_R_t(){
		
		R_t_list = new ArrayList<List<Double>>(n_usedObservations);
		
		double [][] stand_obs_vars = calc_standardized_obs_variables();
		double [][] uncond_cov = calc_unconditional_cov(stand_obs_vars);
		ArrayList<List<Double>> Q_list = calc_Q_4_corr(uncond_cov, stand_obs_vars);
		double [][] Q_star = new double [n_observedVariables][n_observedVariables];
		
		for(int t=0; t<n_usedObservations; t++){
				
			double [][] Q_t = MatrixOperations.get_matrix_from_vec(Q_list.get(t), n_observedVariables, n_observedVariables);
				 
			for(int i=0; i<n_observedVariables; i++){
				Q_star[i][i] = 1.0/Math.sqrt(Q_t[i][i]);			
			}
			
			double [][] R_t = MatrixOperations.multiplication(MatrixOperations.multiplication(Q_star, Q_t),Q_star);
			
			R_t_list.add(MatrixOperations.vecAsList(R_t));
			
		}	
		
	}
	
	
	public static double sum_Q_pars(){
		
		double summedPars = 0.0;
		
		summedPars += GeneralMath.sum(alpha_corr)+GeneralMath.sum(beta_corr);

		return summedPars;
		
	}
	
	
	public static void set_GARCH_MV_pars_of_second_step_from_par_vec(double [] par_vec){
		
		alpha_corr = new double [lag4StandDev][1];
		beta_corr  = new double [lag4Corr][1];
		
		int idx = 0;
		
		for(int i=0; i<lag4StandDev; i++){
			alpha_corr[i][0] = par_vec[idx];
			idx++;
		}
		
		for(int i=0; i<lag4Corr; i++){
			beta_corr[i][0] = par_vec[idx];
			idx++;
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_start_values_4_first_step_est_GARCH_MV(){
		
		double [] startValues = new double [n_heteroPars];
		
		int idx = 0;
		
		for(int k=0; k<n_observedVariables; k++){
			
			int lag4vola = listOfLag4Volatility.get(k);
			int lag4res = listOfLag4Residuals.get(k);			
			int n_volaPars = lag4vola+1;
			int n_maPars   = lag4res;
			int n_pars = n_volaPars + n_maPars;
			
			double [][] obsData = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, 0, endIdx, k, k);
			
			int nResiduals = obsData.length;
	    	  
	        for(int i=0; i<nResiduals; i++){
	        	obsData[i][0] = Math.pow(obsData[i][0], 2.0);
	        }
			
	    	double [][] y = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(obsData, (startIdx+lag4res), (nResiduals-1), 0, 0);
	    	double [][] X = get_lagged_Y_for_lag(obsData, (startIdx+lag4res), (nResiduals-1), lag4res); 
	    		
	    	LinearRegression obj_lm = new LinearRegression(y, X, false);    	
	    	obj_lm.do_parameter_estimation();

	    	double [][] estMaPars = obj_lm.get_est_parameters();   	
	    	double [][] h_t = MatrixOperations.multiplication(X, estMaPars);    	
	    	double [][] H_lagged = get_lagged_Y_for_lag(h_t, (startIdx+lag4vola), (h_t.length-1), lag4vola);
	    	
	    	int n=0;
	    	
	    	if(H_lagged.length>X.length){
	    		n=X.length;
	    	}else{
	    		n=H_lagged.length;
	    	}
	    	
	    	double [][] X_new = new double [n][n_pars];
	    	
	    	for(int t=0; t<n; t++){
	    		X_new[t][0] = 1.0;
	    		for(int p=0; p<lag4vola; p++){
	    			X_new[t][(1+p)] = H_lagged[t][p+1];
	    		}
	    		for(int p=0; p<lag4res; p++){
	    			X_new[t][(lag4vola+1+p)] = X[t][p+1];
	    		}
	    	}
	    	
	    	h_t = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(h_t, 0, (n-1), 0, 0);
	    	obj_lm = new LinearRegression(h_t, X_new, false);    	
	    	obj_lm.do_parameter_estimation();
	    	
	    	for(int i=0; i<n_pars; i++){
	    		startValues[idx] = (obj_lm.get_est_parameters())[i][0];	
	    		idx++;
	    	}
				
		}
		
		return startValues;
		
	}
	
	
	public static int get_number_of_GARCH_MV_pars(){
		
		int n_pars = 0;
		
		for(int i=0; i<n_observedVariables; i++){
			n_pars += listOfLag4Volatility.get(i)+1;
		}
		
		for(int i=0; i<n_observedVariables; i++){
			n_pars += listOfLag4Residuals.get(i);
		}
		
		return n_pars;
		
	}
	
	
	public static double [][] get_GARCH_volaPars(int numberOfVariable){
		
		int n_pars = listOfVolaPars.get(numberOfVariable).size();
		
		double [][] vola_pars = new double [n_pars][1];
		
		for(int i=0; i<n_pars; i++){
			vola_pars[i][0] = listOfVolaPars.get(numberOfVariable).get(i);
		}
		
		return vola_pars;
			
	}
	
	
	public static double [][] get_GARCH_maPars(int numberOfVariable){
		
		int n_pars = listOfMaPars.get(numberOfVariable).size();
		
		double [][] ma_pars = new double [n_pars][1];
		
		for(int i=0; i<n_pars; i++){
			ma_pars[i][0] = listOfMaPars.get(numberOfVariable).get(i);
		}
		
		return ma_pars;
				
	}
	
	
	public static double [][] get_GARCH_MV_corr_matrix(int timeStep){
		
		double [][] R_t = MatrixOperations.get_matrix_from_vec(R_t_list.get(timeStep), n_observedVariables, n_observedVariables);
		
		return R_t;
		
	}
	
	
	public static void do_mean_adjustment(){
		
		int n_dates      = observed_variables.length;	
		double [][] mean = GeneralMath.mean_vec(observed_variables);
		
		for(int i=0; i<n_observedVariables; i++){
			for(int t=0; t<n_dates; t++){
				observed_variables[t][i] -= mean[0][i];
			}	
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_GARCH_MV_estimates(){
		
        GenGraphics obj_graph = new GenGraphics();
		    
        int n_correlations = 0;
        for(int i=0; i<n_observedVariables; i++){
        	for(int j=(i+1); j<n_observedVariables; j++){
        		n_correlations++;
        	}
        }
        
        int addPlotLines = (int) Math.ceil((double)n_correlations/(double)n_observedVariables);
        
        obj_graph.setNumberOfPlotColums(n_observedVariables);
        obj_graph.setNumberOfPlotRows(2+addPlotLines);
	 	
        obj_graph.setGraphWidth(1000);
        obj_graph.setGraphHeight(600);
        
        String [] title = new String [2*n_observedVariables+n_correlations];
           
        double [][] xAxis = new double [n_usedObservations][1];
        double [][] vola  = new double [n_usedObservations][1];
        double [][] corr  = new double [n_usedObservations][1];
        
        for(int i=0; i<n_usedObservations; i++){
	 		xAxis[i][0] = i+1;
	 	}
        
        int idx = 0;
        
        for(int k=0; k<n_observedVariables; k++){
        	
        	title[idx]   = "Observed vs. fitted Variable";
        	    	 	
    	 	double [][] obs_data = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(observed_variables, startIdx, endIdx, k, k);

    	 	obj_graph.plotLines(xAxis, obs_data, true);	
    	 		 	
    	 	idx++;
        }
        
    	for(int k=0; k<n_observedVariables; k++){
    		
    		title[idx] = "("+ distribution +") MV-GARCH("+listOfLag4Volatility.get(k)+","+listOfLag4Residuals.get(k)+") impl. Volatility";
    		
        	for(int i=0; i<n_usedObservations; i++){
        		vola[i][0] = volatilities[i][k];
        	}
    		
        	obj_graph.plotLines(xAxis, vola, true);
        	
        	idx++;
        	
    	}
        	
    	for(int i=0; i<n_observedVariables; i++){
    		for(int j=(i+1); j<n_observedVariables; j++){
    			title[idx] = "Correlation Variables " + (i+1) + " and " + (j+1);
              	for(int t=0; t<n_usedObservations; t++){
              		double [][] corr_matrix = get_GARCH_MV_corr_matrix(t);
              		corr[t][0] = corr_matrix[i][j];
              	}
              	obj_graph.plotLines(xAxis, corr, true);
              	idx++;
    		}
        }
    	   	
	    String [] yLabel = title; 	 	
	
	 	obj_graph.setTitle(title, null, "10");
	 	obj_graph.setYLabel(yLabel, null, "8");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    		
    	//Load & select data input:
    	String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    	String [] colnames = {"Germany", "France", "UK", "USA"}; //{"Germany", "France", "UK", "USA"};
    	
		InputDataManager inputData = new InputDataManager();		
		inputData.fileReader(file, true, true, true);
    	
    	int nData = inputData.numberOfRows-1;
    	String [] rownames = new String [nData];
    	for(int i=0; i<nData;i++){
    		rownames[i] = Integer.toString(i+1);
    	}
		
		inputData.selectLoadedData(rownames, colnames);
		double [][] obsData = inputData.selectedDblFileData;

		int [] volaLag = {1,1,1,1};
		int [] maLag   = {1,1,1,1};
		int sdLag = 1;
		int corrLag = 1;
		
		int start_idx  = Utilities.getMax(maLag);
		int end_idx = obsData.length;
		
		GARCH_MV obj_garch_mv = new GARCH_MV(obsData, start_idx, end_idx, volaLag, maLag, sdLag, corrLag);
			
		//Using mean adjusted data
		obj_garch_mv.do_mean_adjustment();
		obj_garch_mv.do_MLE_4_GARCH_MV();
		
		//MatrixOperations.print_matrix(arPars);
		//MatrixOperations.print_matrix(volaPars);
		//MatrixOperations.print_matrix(maPars);
		System.out.println(logLikelihood);
		
		plot_GARCH_MV_estimates();
		
	}
	
}
