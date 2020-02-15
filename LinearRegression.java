package Regression;
import DataManagement.InputDataManager;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class LinearRegression {

	boolean use_constant;
	
	//Observations
	double [][] explained_variable;
	double [][] explaining_variables;
	
	int n_observations;
	int n_explaining_variables;
	
	//Estimated parameters & residuals
	double  constant;
	double  constant_error;
    double  constant_t_value;
	double  constant_p_value;
	
	double [][] parameters;
	double [][] parameter_errors;
	double [][] parameter_t_values;
	double [][] parameter_p_values;
	
	double [][] ss_pars;
	
	double [][] fitted_explained_variable;
	
	double [][] covariance;
	
	double ss_total;
	double ss_reg;
	
	double [][] residuals;
	double [][] standardized_residuals;
	double [][] studentized_residuals;
	double [][] external_studentized_residuals;
	double [][] cooks_distance;
	double [][] leverage;
	
	double ss_res;
	double sigma;
		
	double F;
	
	double R_squared;
	double adj_R_squared;
	
	
	public LinearRegression(double [][] y, double [][] X, boolean constant_usage){
		
		use_constant = constant_usage;
		
		explained_variable    = y;
		explaining_variables  = X;
		
		n_explaining_variables = X[0].length;
		n_observations         = y.length;
		
	}
	
	
	//Main routine for doing linear regression
	public void do_parameter_estimation(){
		
		double [][] regressor_matrix = explaining_variables;
	    double [][] J                = get_J_matrix();
		double [][] identity         = MatrixOperations.identity(n_observations);
	    
	    double [][] est_pars;
	    double [][] hat_matrix;
	    
	    double [][] term;
	    double [][] inv_term;
	    
	    double [][] errors;
	    
	    int par_idx_1 = 0;
	    int par_idx_2 = n_explaining_variables-1;
	    
		if(use_constant == true){
			
			double [][] unit_vector = MatrixOperations.unit_vector(n_observations);
			
			regressor_matrix        = MatrixOperations.cbind(regressor_matrix, unit_vector);
			
			 par_idx_1 = 1;
			 par_idx_2 = n_explaining_variables;
			
		}
				
        term                   = MatrixOperations.add(identity, MatrixOperations.scalar_multiplication(-1.0/n_observations, J));
		ss_total               = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(explained_variable),term), explained_variable)[0][0];
		
		inv_term               = MatrixOperations.inverse(MatrixOperations.multiplication(MatrixOperations.transpose(regressor_matrix),regressor_matrix));		
		term                   = MatrixOperations.multiplication(inv_term, MatrixOperations.transpose(regressor_matrix));
		
		est_pars               = MatrixOperations.multiplication(term, explained_variable);
		
		hat_matrix             = MatrixOperations.multiplication(regressor_matrix, term);
		
		term                   = MatrixOperations.add(hat_matrix, MatrixOperations.scalar_multiplication(-1.0/n_observations, J));		
		ss_reg                 = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(explained_variable),term), explained_variable)[0][0];
		
		parameters             = MatrixOperations.get_double_sub_vec(est_pars, par_idx_1, par_idx_2);
	
		fitted_explained_variable = MatrixOperations.multiplication(regressor_matrix, est_pars);
		
		residuals                 = MatrixOperations.add(explained_variable, MatrixOperations.scalar_multiplication(-1.0, fitted_explained_variable));		
		
		ss_res                    = MatrixOperations.multiplication(MatrixOperations.transpose(residuals), residuals)[0][0];
		sigma                     = ss_res/(n_observations-(n_explaining_variables+1.0));
			
		covariance                = MatrixOperations.scalar_multiplication(sigma,inv_term);
		errors                    = MatrixOperations.get_diagonal_from_matrix(covariance);
		
		parameter_errors          = GeneralMath.sqrt(MatrixOperations.get_double_sub_vec(errors, par_idx_1, par_idx_2));
		
		calculate_t_statistics_4_est_pars();
		
		if(use_constant == true){
			
			constant          = MatrixOperations.get_double_sub_vec(est_pars, 0, 0)[0][0];
			constant_error    = Math.sqrt(MatrixOperations.get_double_sub_vec(errors, 0, 0)[0][0]);	
			constant_t_value  = constant/constant_error;
			
		}
		
		F = (ss_reg/n_explaining_variables)/sigma;
		
		R_squared     = ss_reg/ss_total;
		adj_R_squared = 1.0-(n_observations-1.0)/(n_observations-(n_explaining_variables+1.0))*(1.0-R_squared);
		
		calculate_residual_measures(hat_matrix);
		
		calculate_sums_of_squares_4_est_pars(regressor_matrix);
				
		//print_estimation_res();
		//print_linear_model_stats();
		//print_covariance_matrix();
		//print_linear_model_diagnostics();
		
	}
	
	
	// returns matrix of ones
	public double [][] get_J_matrix(){
		
		double [][] J_matrix = new double [n_observations][n_observations];
		
		for(int i=0; i<n_observations; i++){
			
			for(int j=0; j<n_observations; j++){
				
				J_matrix[i][j] = 1.0;
				
			}
			
		}
		
		return J_matrix;
		
	}
	
	
	// returns the t-values for the estimated parameters
	public void calculate_t_statistics_4_est_pars(){
		
		parameter_t_values = new double [n_explaining_variables][1];
		
		for(int i=0; i<n_explaining_variables; i++){
			
			parameter_t_values[i][0] = parameters[i][0]/parameter_errors[i][0];
			
		}
			
	}
	
	
	// calculates studentized, external studentized and Cook´s distance
	public void calculate_residual_measures(double [][] hat_matrix){
		
		standardized_residuals         = new double [n_observations][1];
		studentized_residuals          = new double [n_observations][1];
		external_studentized_residuals = new double [n_observations][1];
		cooks_distance                 = new double [n_observations][1];
		leverage                       = new double [n_observations][1];
		
		double sqrt_sigma = Math.sqrt(sigma);
		
		for(int i=0; i<n_observations; i++){
			
			standardized_residuals[i][0]         = residuals[i][0]/sqrt_sigma;
			studentized_residuals[i][0]          = residuals[i][0]/Math.sqrt(sigma*(1.0-hat_matrix[i][i]));
			external_studentized_residuals[i][0] = residuals[i][0]*Math.pow((n_observations-n_explaining_variables)/(ss_res*(1.0-hat_matrix[i][i])-Math.pow(residuals[i][0], 2.0)), 0.5);
			cooks_distance[i][0]                 = Math.pow(studentized_residuals[i][0] , 2.0)/(n_explaining_variables+1.0)*(hat_matrix[i][i])/(1.0-hat_matrix[i][i]);
			
		}
		
		leverage = MatrixOperations.get_diagonal_from_matrix(hat_matrix);
		
	}
	
	
	// calculates sums of squares for the parameters beta_1, beta_2, ...
	public void calculate_sums_of_squares_4_est_pars(double [][] regressor_matrix){
		
		ss_pars = new double [n_explaining_variables][1];
		
		double [][] J = get_J_matrix();
		double [][] design_matrix;
		double [][]	term;
		
		double a = 0.0;
		
		if(use_constant == true){
			
			a = 1.0;
			
		}
		
		for(int i=1; i<n_explaining_variables+a; i++){
			
			design_matrix = MatrixOperations.get_sub_matrix_between_column_idxs(regressor_matrix, 0, i);
			
			term          = MatrixOperations.inverse(MatrixOperations.multiplication(MatrixOperations.transpose(design_matrix), design_matrix));
			
			design_matrix = MatrixOperations.multiplication(MatrixOperations.multiplication(design_matrix, term), MatrixOperations.transpose(design_matrix));
			
			term          = MatrixOperations.add(design_matrix, MatrixOperations.scalar_multiplication(-1.0/n_observations, J));
			
			ss_pars[i-1][0] = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(explained_variable), term), explained_variable)[0][0];
			
			if(i>1){
				
				ss_pars[i-1][0] = ss_pars[i-1][0]-ss_pars[(i-2)][0];
				
			}
			
		}
		
	}
	
		
	// prints estimated parameters, standard errors, t-values and p-values
	public void print_estimation_res(){
		
		int n_est_pars = n_explaining_variables;
		
		if(use_constant == true){
			
			n_est_pars = n_est_pars + 1; 
			
		}
		
		String [][] est_res = new String [4][n_est_pars+1];
		
		est_res[0][0] = "           ";
		est_res[1][0] = "Estimates  ";
	    est_res[2][0] = "Std. errors";
	    est_res[3][0] = "t-values   ";
		
	    int par_idx = 0;
	    
		for(int i=1; i<n_est_pars+1; i++){
			
			if(i==1 && use_constant == true){
				
				est_res[0][i] = "    Intercept";
				est_res[1][i] = Double.toString(constant);
				est_res[2][i] = Double.toString(constant_error);
				est_res[3][i] = Double.toString(constant_t_value);
				
			}else{
				
				est_res[0][i] = "        Parameter " + (i-1);
				est_res[1][i] = Double.toString(parameters[par_idx][0]);
				est_res[2][i] = Double.toString(parameter_errors[par_idx][0]);
				est_res[3][i] = Double.toString(parameter_t_values[par_idx][0]);
				
				par_idx = par_idx + 1;
				
			}
			
		}

		System.out.println("");
		System.out.println("Linear Model Estimation Results:");
		System.out.println("");
		
		MatrixOperations.print_matrix(est_res);
		
	}
	
	
	// prints estimated F-test, R^2, adjusted R^2, sigma
	public void print_linear_model_stats(){
			
		String [][] stats = new String [1][8];
		
		stats[0][0] = "F-Test:";
		stats[0][1] = Double.toString(F);
		stats[0][2] = "  R squared:";
		stats[0][3] = Double.toString(R_squared);
		stats[0][4] = "  adj. R squared:";
		stats[0][5] = Double.toString(adj_R_squared);
		stats[0][6] = "  Sigma:";
		stats[0][7] = Double.toString(sigma);
	    
		System.out.println("");
		MatrixOperations.print_matrix(stats);
		
		String [][] facts = new String [1][8];
		
		facts[0][0] = "Number of Observations:";
		facts[0][1] = Double.toString(n_observations);
		facts[0][2] = " Number of explaining Variables:";
		facts[0][3] = Double.toString(n_explaining_variables);
		facts[0][4] = " DF for Residuals:";
		facts[0][5] = Double.toString(n_observations-(n_explaining_variables+1.0));
		facts[0][6] = " DF for F-Test:";
		facts[0][7] = ""+ Double.toString(n_explaining_variables) + " and " + (n_observations-(n_explaining_variables+1.0));	
		
		MatrixOperations.print_matrix(facts);
		
	}
	
	
	// prints observed variable y, fitted variable, residuals and residual measures 
	public void print_linear_model_diagnostics(){
		
		String [][] diagnostics = new String [n_observations+1][8];
		
		diagnostics[0][0] = "Observed Values";
		diagnostics[0][1] = "Fitted Values  ";
		diagnostics[0][2] = "Residuals";
		diagnostics[0][3] = "Standardized Residuals";
		diagnostics[0][4] = "Studentized Residuals";
		diagnostics[0][5] = "Ext. Studentized Residuals";
		diagnostics[0][6] = "Leverage";
		diagnostics[0][7] = "Cook´s Distance";
	    
	    for(int i=0; i<n_observations; i++){
	    	
	    	diagnostics[i+1][0] = Double.toString(explained_variable[i][0]);
	    	diagnostics[i+1][1] = Double.toString(fitted_explained_variable[i][0]);
	    	diagnostics[i+1][2] = Double.toString(residuals[i][0]);
	    	diagnostics[i+1][3] = Double.toString(standardized_residuals[i][0]);
	    	diagnostics[i+1][4] = Double.toString(studentized_residuals[i][0]);
	    	diagnostics[i+1][5] = Double.toString(external_studentized_residuals[i][0]);
	    	diagnostics[i+1][6] = Double.toString(leverage[i][0]);
	    	diagnostics[i+1][7] = Double.toString(cooks_distance[i][0]);
	    	
	    }
	    
		System.out.println("");
		System.out.println("Linear Model Diagnostics:");
		System.out.println("");
	    
	    MatrixOperations.print_matrix(diagnostics); 
	    
	}
	
	
	public double get_est_constant(){
		
		return constant;
		
	}
	
	
	public double get_constant_error(){
		
		return constant_error;
		
	}
	
	
	public double get_constant_t_value(){
		
		return constant_t_value;
		
	}
	
	
	public double get_constant_p_value(){
		
		return constant_p_value;
		
	}
	
	
	public double [][] get_est_parameters(){
		
		return parameters;
		
	}
	
	
	public double [][] get_parameter_errors(){
		
		return parameter_errors;
		
	}
	
	
	public double [][] get_parameter_t_values(){
		
		return parameter_t_values;
		
	}
	
	
	public double [][] get_parameter_p_values(){
		
		return parameter_p_values;
		
	}
	
	
	public double get_sigma(){		
		return sigma;		
	}
	
	
	public double get_R_squared(){		
		return R_squared;		
	}
	
	
	public double get_adj_R_squared(){		
		return adj_R_squared;		
	}
	
	
	public double get_F_stat(){		
		return F;		
	}
	
	
	public double [][] get_residuals(){		
		return residuals;		
	}
	
	
	public double [][] get_fitted_values(){
		return fitted_explained_variable;
	}
	
	
	// prints observed variable y, fitted variable, residuals and residual measures 
	public void print_covariance_matrix(){
			
		int n_est_pars = n_explaining_variables;
		
		if(use_constant == true){
			
			n_est_pars = n_est_pars + 1; 
			
		}
		
		String [][] cov = new String [n_est_pars+1][n_est_pars+1];
		cov[0][0] = "         ";
		
		for(int i=0; i<n_est_pars+1; i++){
			
			for(int j=0; j<n_est_pars+1; j++){
				
				if(i==0){
					
					if(j == 1){
						
						cov[i][j] = "Intercept";
						
					}
					
					if(j > 1){
						
						cov[i][j] = "    Parameter " + (j-1);
						
					}
							
				}
				
				if(j==0){
					
					if(i == 1){
						
						cov[i][j] = "Intercept  ";
						
					}
					
					if(i > 1){
						
						cov[i][j] = "Parameter " + (i-1);
						
					}
							
				}
				
				if(i!=0 && j!=0){
					
					cov[i][j] = Double.toString(covariance[i-1][j-1]);
					
				}
					
			}
			
		}
		
		System.out.println("");
		System.out.println("Variance-Covariance-Matrix of estimated Parameters:");
		System.out.println("");
		
		MatrixOperations.print_matrix(cov);
		
	}
	
	
	public static void test1() {
		
    	double[][] data = null;
    	
    	try {
    		//data = ReadTextFile.readfile("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/RegTest.txt");
    		String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test_MS_Models/InterestRates.txt";
    		String [] colnames = {"Germany"};
        	
    		InputDataManager inputData = new InputDataManager();		
    		inputData.fileReader(file, true, true, true);
        	
        	int nData = inputData.numberOfRows-1;
        	String [] rownames = new String [nData];
        	for(int i=0; i<nData;i++){
        		rownames[i] = Integer.toString(i+1);
        	}
    		
    		inputData.selectLoadedData(rownames, colnames);
    		data = inputData.selectedDblFileData;
    		
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
		
    	int lag = 30;
    	int n_usedObs = data.length-lag;
    	double [][] laggedData = new double [n_usedObs][lag];
    	double [][] unlaggedData = new double [n_usedObs][1];
    	
    	int startIdx = lag-1;
    	
        for(int i=0; i<n_usedObs; i++) {
        	unlaggedData[i][0] = data[startIdx+i+1][0];
        	for(int l=0; l<lag; l++) {
        		laggedData[i][l] = data[startIdx-l+i][0];
        	}	
        }
    	
        double [][] y = MatrixOperations.get_sub_matrix_between_column_idxs(unlaggedData, 0, 0);
    	double [][] X = MatrixOperations.get_sub_matrix_between_column_idxs(laggedData, 1, lag-1); 
    	
    	LinearRegression obj_lm = new LinearRegression(y, X, true);
    	long startTime = System.currentTimeMillis();
    	obj_lm.do_parameter_estimation();
    	MatrixOperations.print_matrix(obj_lm.get_est_parameters());
    	long endTime = System.currentTimeMillis();
    	System.out.println("Finished AR-regression after " + ((endTime-startTime)/1000.0) + " secs.");
    	
	}
	
	
	public static void test2() {
		
    	double[][] data = null;
    	
    	try {
    		//data = ReadTextFile.readfile("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/RegTest.txt");
    		String file = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/RegTest.txt";
    		String [] colnames = {"2", "3", "4", "5", "6"};
        	
    		InputDataManager inputData = new InputDataManager();		
    		inputData.fileReader(file, true, true, true);
        	
        	int nData = inputData.numberOfRows-1;
        	String [] rownames = new String [nData];
        	for(int i=0; i<nData;i++){
        		rownames[i] = Integer.toString(i+1);
        	}
    		
    		inputData.selectLoadedData(rownames, colnames);
    		data = inputData.selectedDblFileData;
    		
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	
    	double [][] y = MatrixOperations.get_sub_matrix_between_column_idxs(data, 4, 4);
    	double [][] X = MatrixOperations.get_sub_matrix_between_column_idxs(data, 0, 3);    	
    	
    	LinearRegression obj_lm = new LinearRegression(y, X, true);
    	long startTime = System.currentTimeMillis();
    	obj_lm.do_parameter_estimation();
    	MatrixOperations.print_matrix(obj_lm.get_est_parameters());
    	long endTime = System.currentTimeMillis();
    	System.out.println("Finished regression after " + ((endTime-startTime)/1000.0) + " secs.");
	}
	
	
	// test client
    public static void main(String[] args) {
    	
    	test1();
    	
    }
	 	
}
