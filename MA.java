package TimeSeriesAnalysis;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;

public class MA extends ARIMA{
		
	/**
     * Constructor for Moving Average Model class
     * @param obs_variables double vector containing data series
     * @param start_idx integer index indicating start of the analysis
     * @param end_idx integer index indicating end of the analysis
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	//TODO: Define additional constructors!
	//public MA(double [][] obs_variables, int start_idx, int end_idx, int obsLag){
	//	super(obs_variables, start_idx, end_idx, obsLag);
	//	model = "MA";
	//    demean_series_4_ma();				
	//}
	
	
    /**
     * Constructor for Moving Average Model class
     * @param obs_variables double vector containing data series
     * @param obsLag integer lag order (max order if automatic selection is used)
     */
	public MA(double [][] obs_variables, int obsLag){
		super(obs_variables, 0, obsLag);
		model = "MA";				
	}
	
	
	public static void calc_info_criteria_4_ma() {
		calc_info_criteria_4_arima();
	}
	
	
	public static void do_ma_estimation() {
		do_arima_estimation();
	}
	
	
	public static void do_ma_estimation(boolean doStatisticInference) {
		do_arima_estimation(doStatisticInference);
	}
	
	
	public static void do_automatic_ma_model_estimation_and_selection(String usedIC){
		
		ic = usedIC;		
		int orgLag = lag;
		ic4modelSelection = new double [orgLag];
		ArrayList<List<Double>> parList = new ArrayList<List<Double>>(orgLag);
		
		for(int i=0; i<orgLag; i++){
			lag = i+1;
			do_ma_estimation(false);
			List<Double> estPars = MatrixOperations.vecAsList(parameters);
			estPars.add(sigma);
			parList.add(estPars);
			
			if(ic == "aic"){
				ic4modelSelection[i] = aic;
			}
			
			if(ic == "bic"){
				ic4modelSelection[i] = bic;
			}
			
		}
		
		double min_ic  = Utilities.Utilities.getMin(ic4modelSelection);		
		int [] min_idxs = Utilities.Utilities.get_idx(ic4modelSelection, min_ic);
		int optLag = min_idxs[0]+1;
			
		System.out.println("Lag order " + optLag + " identified by minimizing " + ic + ".");
		
		lag = optLag;
		
		double [] optPars = new double [optLag+1];
		
		for(int i=0; i<optLag+1; i++) {
			optPars[i] = parList.get(min_idxs[0]).get(i);
		}
		
		set_arima_pars_form_vec(optPars);
		calc_exact_log_likelihood_by_KF();
		log_likelihood *=(-1);
		calc_info_criteria_4_ma();
		calc_arima_fitted_values_and_residuals();
		calc_se_and_t_statistics_of_arima_from_log_lik();
		
		set_arima_pars_form_vec(optPars);
		calc_exact_log_likelihood_by_KF();
		log_likelihood *=(-1);
		
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

		int lag = 1;
		//int start_idx = lag;
		//int end_idx = obsData.length;
		
		
		MA obj_ma = new MA(obsData,lag);	
		
		obj_ma.do_ma_estimation();
		
		//obj_ma.do_automatic_ma_model_estimation_and_selection("aic");
		
		obj_ma.plot_time_series_and_fitted_values_with_residuals();
		//obj_ma.plot_time_series_and_fitted_values();
		//obj_ma.plot_residuals();		
		
		//obj_ma.plot_ic_of_model_selection();
		//obj_ma.create_qq_plot();
		//obj_ma.plot_standardized_residuals();
		//obj_ma.plot_acf_of_residuals();
		//obj_ma.plot_pacf_of_residuals();
		//obj_ma.plot_histogram_of_residuals();
		
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
