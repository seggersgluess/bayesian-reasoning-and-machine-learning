package TimeSeriesAnalysis;

import DataManagement.InputDataManager;
import Mathematics.GeneralMath;
import Mathematics.InformationCriteria;
import Mathematics.MatrixOperations;
import Optimization.L_BFGS;
import Optimization.NumDeriv;
import StateSpaceModels.KalmanFilter;

public class ARIMA extends TimeSeriesAnalysis{

	static int arLag;
	static int maLag;
	static int integration;
	
	static int n_stateVars;
	
	static KalmanFilter kf;
	
	static double [][] hessian;
	
	static int n_iterations = 1000;
	static boolean convergence = false;
	
	public ARIMA(double[][] obs_variables, int used_arLag, int used_maLag) {
		super(obs_variables, 0);
			
		if(used_arLag<0) {
			throw new RuntimeException("Invalid AR lag order supplied.");
		}
		
		if(used_maLag<0) {
			throw new RuntimeException("Invalid MA lag order supplied.");
		}
		
		arLag = used_arLag;
		maLag = used_maLag;
		
		model = "ARIMA";
		int maxLag = arLag;
		if(maLag>arLag) {
			maxLag = maLag;
		}
		
		startIdx               = maxLag;
		endIdx                 = observed_variables.length-1;
		
		n_usedObservations = endIdx-startIdx+1;
		
		if(maLag == 0) {
			n_stateVars = arLag;
		}else {
			n_stateVars = Math.max(arLag, (maLag+1));
		}
			
		demean_series_4_arima();
		
	}
	
	
	static double [][] get_arima_measurement_matrix(){
		
		double [][] m_matrix = new double [1][n_stateVars];
		m_matrix[0][0] = 1.0;
		
		if(maLag != 0) {
			for(int i=0; i<maLag; i++) {
				m_matrix[0][i+1] = parameters[(arLag+i)][0];
			}
		}
				
		return m_matrix;
			
	}
	
	
	static double [][] get_arima_measurement_cov(){
		
		double [][] m_cov = new double [1][1];
		
		return m_cov;
			
	}
	
	
	static double [][] get_arima_transition_matrix(){
		
		double [][] t_matrix = new double [n_stateVars][n_stateVars];
		
		//ARMA model
		if(arLag != 0 && maLag != 0) {
			for(int i=0; i<arLag; i++) {
				t_matrix[0][i] = parameters[i][0];
			}
			
			int n_identity = n_stateVars-1;
			for(int i=0; i<n_identity; i++) {
				t_matrix[i+1][i] = 1.0;
			}
			
		}
		
		//AR model
		if(arLag != 0 && maLag == 0) {
			for(int i=0; i<n_stateVars; i++) {
				t_matrix[0][i] = parameters[i][0];
			}

			int arLagMod = arLag-1;
			for(int i=0; i<arLagMod; i++) {
				t_matrix[i+1][i] = 1.0;
			}
			
		}
		
		//MA model
		if(arLag == 0 && maLag != 0) {
			for(int i=0; i<maLag; i++) {
				t_matrix[i+1][i] = 1.0;
			}			
		}
		
		return t_matrix;
			
	}
	
	
	static double [][] get_arima_transition_covariance(){
		
		double [][] cov = new double [n_stateVars][n_stateVars]; 
		
		cov[0][0] = sigma;
				
		return cov;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void calc_exact_log_likelihood_by_KF(){
			
		double [][] m_matrix = get_arima_measurement_matrix();
		double [][] m_cov    = new double [1][1];
		
		double [][] t_matrix = get_arima_transition_matrix();
		double [][] t_cov    = get_arima_transition_covariance();
					
		kf = new KalmanFilter(observed_variables, null, m_matrix, m_cov, null, t_matrix, t_cov);		
		kf.run_kalman_filter();
		
		log_likelihood = (-1.0)*kf.get_kf_log_likelihood();		
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_arima_estimation(){
		
		double [] init_values = get_init_values_4_MLE_of_arima();
		
		L_BFGS optim = new L_BFGS(ARIMA::calc_exact_log_likelihood_by_KF_from_par_vec, n_iterations);
		optim.set_convergence_criterion(1e-05);
		optim.do_LBFGS_Optimization(init_values);
		set_arima_pars_form_vec(optim.get_optimal_candidate());
		log_likelihood = (-1.0)*optim.get_optimal_value();
		calc_exact_log_likelihood_by_KF();
		log_likelihood *=(-1);
		calc_info_criteria_4_arima();
		calc_arima_fitted_values_and_residuals();
		
		calc_se_and_t_statistics_of_arima_from_log_lik();
		set_arima_pars_form_vec(optim.get_optimal_candidate());
		calc_exact_log_likelihood_by_KF();
		log_likelihood *=(-1);
		
		convergence = optim.get_convergence_info();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void do_arima_estimation(boolean doStatisticInference){
		
		double [] init_values = get_init_values_4_MLE_of_arima();
		
		L_BFGS optim = new L_BFGS(ARIMA::calc_exact_log_likelihood_by_KF_from_par_vec, n_iterations);
		optim.set_convergence_criterion(1e-05);
		optim.do_LBFGS_Optimization(init_values);
		set_arima_pars_form_vec(optim.get_optimal_candidate());
		log_likelihood = (-1.0)*optim.get_optimal_value();
		calc_exact_log_likelihood_by_KF();
		log_likelihood *=(-1);
		calc_info_criteria_4_arima();
		calc_arima_fitted_values_and_residuals();
		
		calc_se_and_t_statistics_of_arima_from_log_lik();
		set_arima_pars_form_vec(optim.get_optimal_candidate());
		calc_exact_log_likelihood_by_KF();
		log_likelihood *=(-1);
		
		convergence = optim.get_convergence_info();
		
		if(doStatisticInference = true) {
			calc_se_and_t_statistics_of_arima_from_log_lik();
			set_arima_pars_form_vec(optim.get_optimal_candidate());
			calc_exact_log_likelihood_by_KF();
			log_likelihood *=(-1);
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void calc_info_criteria_4_arima() {
		
		InformationCriteria ic = new InformationCriteria();
		
		//number of ma coeffs + number of ar_coeffs + sigma
		int nPars = parameters.length+1;
		
		aic = ic.aic(log_likelihood, nPars, true);
		bic = ic.bic(log_likelihood, nPars, n_usedObservations, true);
		
	}
	
	
	public static double calc_exact_log_likelihood_by_KF_from_par_vec(double [] par_vec, double [] further_ars){
		
		//par_vec = [ar_coeff_1,...,ar_coeff_p, ma_coeff_1,...,ma_coeff_q, sigma]		
		set_arima_pars_form_vec(par_vec);
		
		calc_exact_log_likelihood_by_KF();
		
		return log_likelihood;
		
	}
	
	
	public static void set_arima_pars_form_vec(double [] par_vec) {
		
		//par_vec = [ar_coeff_1,...,ar_coeff_p, ma_coeff_1,...,ma_coeff_q, sigma]
		int nPars = arLag + maLag;
		parameters = new double [nPars][1];
		
		int idx = 0;
		for(int i=0; i<nPars; i++){
			parameters[i][0] = par_vec[i];
			idx++;
		}
		
		sigma = par_vec[idx];
		
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] get_init_values_4_MLE_of_arima() {
		
		int nPars = arLag+maLag+1;
		
		double [] init_values = new double [nPars];
		
		int idx = 0;
		
		if(arLag != 0) {	
			//Initial variables from OLS estimate:
			double [][] obs_vars = observed_variables;
			int orgStartIdx = startIdx;
			int orgEndIdx   = endIdx;
					
			AR obj_ar = new AR(observed_variables, arLag);
			obj_ar.do_ols_ar_estimation();
						
			double [] initARpars = MatrixOperations.convVecToArray(obj_ar.get_est_parameters());
			
			observed_variables = obs_vars;
			startIdx = orgStartIdx;
			endIdx = orgEndIdx;
			model = "ARIMA";
			n_usedObservations = endIdx-startIdx+1;
			
			for(int i=0; i<arLag; i++){
				init_values[i] = initARpars[i];
				idx++;
			}
		}
	
		for(int i=0; i<maLag; i++) {
			init_values[idx+i] = 0.5/(i+1.0)-0.1;
		}
		
		init_values[(nPars-1)] = 1.0/((double) nPars);
		
		return init_values;
		
	}
	
	
	public static void calc_se_and_t_statistics_of_arima_from_log_lik(){
		
		int nPars = parameters.length+1;
		double [] optPars   = new double [nPars];
		
		for(int i=0; i<(nPars-1); i++) {
			optPars[i] = parameters[i][0];
		}
		
		optPars[(nPars-1)] = sigma;
		
		nPars--;
		
		standErrors_parameters  = new double [nPars][1];
		tValues_parameters     = new double [nPars][1];
		
		hessian = NumDeriv.hessian(ARIMA::calc_exact_log_likelihood_by_KF_from_par_vec, optPars, null);
		double [][] inv_hessian = MatrixOperations.inverse(hessian);
		
		for(int i=0; i<nPars; i++) {
			standErrors_parameters[i][0] = Math.sqrt(inv_hessian[i][i]);
			tValues_parameters[i][0] = parameters[i][0]/standErrors_parameters[i][0];
		}
				
	}

	
	@SuppressWarnings("static-access")
	public static void calc_arima_fitted_values_and_residuals(){
		fittedValues = new double [n_usedObservations][1];
		residuals = new double [n_usedObservations][1];
		standardized_residuals = new double [n_usedObservations][1];
		
		double sd = Math.sqrt(sigma);
		
		for(int t=0; t<n_usedObservations; t++){
			int curIdx = startIdx+t;
			residuals[t][0] = kf.get_kf_residual_vec(t)[0][0];
			fittedValues[t][0] = observed_variables[curIdx][0]-residuals[t][0];
			standardized_residuals[t][0] = residuals[t][0]/sd;
		}
	}
	
	
	public static void demean_series_4_arima(){
		
		double mean = GeneralMath.mean(observed_variables);
		
		for(int t=0; t<n_observations; t++){
			observed_variables[t][0] -= mean;
		}
		
	}
	
	
	public static double [][] get_hessian() {
		return hessian;
	}

	
	public static boolean convergence(){
		return convergence;
	}
			
	
	public static KalmanFilter get_kf() {
		return kf;
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

		int i_arLag = 1;
		int i_maLag = 1;
		//int start_idx = lag;
		//int end_idx = obsData.length;
		
		
		ARIMA obj_arima = new ARIMA(obsData,i_arLag,i_maLag);	
		
		obj_arima.do_arima_estimation();
		
		obj_arima.plot_time_series_and_fitted_values_with_residuals();
		//obj_arima.plot_time_series_and_fitted_values();
		//obj_arima.plot_residuals();		
		
		//obj_arima.plot_ic_of_model_selection();
		//obj_arima.create_qq_plot();
		//obj_arima.plot_standardized_residuals();
		//obj_arima.plot_acf_of_residuals();
		//obj_arima.plot_pacf_of_residuals();
		//obj_arima.plot_histogram_of_residuals();
		
		System.out.println("Parameter estimates:");
		MatrixOperations.print_matrix(parameters);
		System.out.println(sigma);
		//System.out.println(intercept);
		MatrixOperations.print_matrix(standErrors_parameters);
		System.out.println(log_likelihood);
		//MatrixOperations.print_matrix(residuals);
		
		//calc_information_criteria();
		
	}
	
	
}
