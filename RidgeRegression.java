

public class RidgeRegression {

	static boolean use_constant;
	
	double    penalty;
	double [] penalty_values;
	
	//Observations
	static double [][] explained_variable;
	static double [][] explaining_variables;
	
	static int n_observations;
	static int n_explaining_variables;
	
	static double  constant = 0.0;
	
	static double [][] parameters;
	static double [][] parameter_errors;
	static double [][] parameter_t_values;
	static double [][] variance_inflation;
	
	static double [][] unstandardized_parameters;
	static double [][] unstandardized_errors;
	
	// parameters for various penalties to generate ridge coefficient path
	static double [][] parameter_values_4_penalties;
	static double [][] variance_inflation_4_penalties;
	
	static double [][] covariance;
	
	static double [][] fitted_explained_variable;
	static double [][] unstand_fitted_explained_variable;
	
	static double [][] residuals;

	static double ss_total;
	static double ss_reg;
	static double ss_res;
	static double penalized_ss_res;
	static double sigma;
	
	static double R_squared;
	static double adj_R_squared;
	
	static double F;
	
	public RidgeRegression(double [][] y, double [][] X, double penalty, boolean constant_usage){
		
		use_constant = constant_usage;
		
		this.penalty          = penalty;
		explained_variable    = y;
		explaining_variables  = X;
		
		n_explaining_variables = X[0].length;
		n_observations         = y.length;
		
	}
	
		
	// main routine for doing linear regression
	public void do_parameter_estimation(){
		
		double [][] X = new double [n_observations][n_explaining_variables];
		double [][] y = new double [n_observations][1];
		
		double [][] J = get_J_matrix();
		
			
		X = getStandardizedValues(explaining_variables);
		y = getStandardizedValues(explained_variable);
			
		double [][] identity      = MatrixOperations.identity(n_observations);
		   		    
		double [][] penalty_term;
		double [][] inv_term;
		double [][] term; 
		double [][] XXt;
		double [][] hat_matrix;    
		
		term                      = MatrixOperations.add(identity, MatrixOperations.scalar_multiplication(-1.0/n_observations, J));
		ss_total                  = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(y),term), y)[0][0];
		
		identity                  = MatrixOperations.identity(n_explaining_variables);
		penalty_term              = MatrixOperations.scalar_multiplication(penalty, identity);
			
		XXt                       = MatrixOperations.multiplication(MatrixOperations.transpose(X),X);
		inv_term                  = MatrixOperations.add(penalty_term, XXt);
		inv_term                  = MatrixOperations.inverse(inv_term);		
		term                      = MatrixOperations.multiplication(inv_term, MatrixOperations.transpose(X));			
		parameters                = MatrixOperations.multiplication(term, y);
	
		fitted_explained_variable = MatrixOperations.multiplication(X, parameters);
			
		hat_matrix                = MatrixOperations.multiplication(X, term);
		term                      = MatrixOperations.add(hat_matrix, MatrixOperations.scalar_multiplication(-1.0/n_observations, J));
			
		residuals                 = MatrixOperations.add(y, MatrixOperations.scalar_multiplication(-1.0, fitted_explained_variable));			
		
		ss_reg                    = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(y),term), y)[0][0];
		ss_res                    = MatrixOperations.multiplication(MatrixOperations.transpose(residuals), residuals)[0][0];
		penalized_ss_res          = ss_res + penalty*Math.pow(MatrixOperations.euclidian(parameters),2.0);
		sigma                     = Math.sqrt(penalized_ss_res/(n_observations-n_explaining_variables));
		
		term                      = MatrixOperations.multiplication(inv_term,MatrixOperations.multiplication(XXt, inv_term));
		covariance                = MatrixOperations.scalar_multiplication(Math.pow(sigma,2.0),term);
		parameter_errors          = GeneralMath.sqrt(MatrixOperations.get_diagonal_from_matrix(covariance));
		
		variance_inflation        = MatrixOperations.scalar_multiplication((n_observations-1.0),MatrixOperations.get_diagonal_from_matrix(term));
		
		calculate_t_statistics_4_est_pars();
		
		F = (ss_reg/n_explaining_variables)/Math.pow(sigma,2.0);
		
		R_squared                 = ss_reg/ss_total;
		adj_R_squared             = 1.0-(n_observations-1.0)/(n_observations-(n_explaining_variables+1.0))*(1.0-R_squared);
		
		caluculate_unstandardized_parameters();
		
		if(use_constant == true){
			
			double [][] const_vec = MatrixOperations.scalar_multiplication(constant,MatrixOperations.unit_vector(n_observations));
			
			unstand_fitted_explained_variable = MatrixOperations.add(const_vec,MatrixOperations.multiplication(explaining_variables, parameters));
			
		}else{
			
			unstand_fitted_explained_variable = MatrixOperations.multiplication(explaining_variables, parameters);
			
		}
		
	}
	
	
	// returns ridge regression coefficient and variance inflation path
	public void generate_ridge_regression_path(double [] penalty_values){
		
		this.penalty_values = penalty_values;
		
		int n_penalties = penalty_values.length;
		
		parameter_values_4_penalties   = new double[n_explaining_variables][n_penalties];
		variance_inflation_4_penalties = new double[n_explaining_variables][n_penalties];
		
		double [][] X = new double [n_observations][n_explaining_variables];
		double [][] y = new double [n_observations][1];
			
		X = getStandardizedValues(explaining_variables);
		y = getStandardizedValues(explained_variable);
				    
		double [][] penalty_term;
		double [][] inv_term;
		double [][] term; 
		double [][] XXt;
  
		double [][] pars4penalty;
		double [][] variance_inf_4_penalty;
		
		double [][] identity      = MatrixOperations.identity(n_explaining_variables);
		
		for(int i=0; i<n_penalties; i++){
			
			penalty_term              = MatrixOperations.scalar_multiplication(penalty_values[i], identity);
			
			XXt                       = MatrixOperations.multiplication(MatrixOperations.transpose(X),X);
			inv_term                  = MatrixOperations.add(penalty_term, XXt);
			inv_term                  = MatrixOperations.inverse(inv_term);		
			term                      = MatrixOperations.multiplication(inv_term, MatrixOperations.transpose(X));			
			
			pars4penalty              = MatrixOperations.multiplication(term, y);
			
			term                      = MatrixOperations.multiplication(inv_term,MatrixOperations.multiplication(XXt, inv_term));

			variance_inf_4_penalty    = MatrixOperations.scalar_multiplication((n_observations-1.0),MatrixOperations.get_diagonal_from_matrix(term));
			
			for(int j=0; j<n_explaining_variables; j++){
				
				parameter_values_4_penalties[j][i]   = pars4penalty[j][0]; 
				variance_inflation_4_penalties[j][i] = variance_inf_4_penalty[j][0];
						
			}
				
		}
		
	}
	
	
	// returns matrix of ones
	public static double [][] get_J_matrix(){
		
		double [][] J_matrix = new double [n_observations][n_observations];
		
		for(int i=0; i<n_observations; i++){
			
			for(int j=0; j<n_observations; j++){
				
				J_matrix[i][j] = 1.0;
				
			}
			
		}
		
		return J_matrix;
		
	}
	
	
	// returns the t-values for the estimated parameters
	public static void calculate_t_statistics_4_est_pars(){
		
		parameter_t_values = new double [n_explaining_variables][1];
		
		for(int i=0; i<n_explaining_variables; i++){
			
			parameter_t_values[i][0] = parameters[i][0]/parameter_errors[i][0];
			
		}
			
	}
		
	
	// returns unstandardized parameters
	public static void caluculate_unstandardized_parameters(){
		
		unstandardized_parameters = new double [n_explaining_variables][1];
		unstandardized_errors     = new double [n_explaining_variables][1];
		
		for(int i=0; i<n_explaining_variables; i++){
			
			unstandardized_parameters[i][0] = parameters[i][0]*GeneralMath.sd(explained_variable)/GeneralMath.sd(MatrixOperations.get_column_from_matrix(explaining_variables,i));
			unstandardized_errors[i][0]     = parameter_errors[i][0]*GeneralMath.sd(explained_variable)/GeneralMath.sd(MatrixOperations.get_column_from_matrix(explaining_variables,i));
			
		}
		
		if(use_constant == true){
			
			constant = GeneralMath.mean(explained_variable)-MatrixOperations.multiplication(GeneralMath.mean_vec(explaining_variables), unstandardized_parameters)[0][0];

		}
		
	}
	
	
	// returns standardized values [x-mean(x)]/sd(x)
	public static double [][] getStandardizedValues(double [][] x){
		
		int nRows = x.length;
		int nCols = x[0].length;
		
		double mean;
		double sd;
		
		double [][] standardizedValues = new double [nRows][nCols];
		
		for(int i=0; i<nCols; i++){
			
			mean = GeneralMath.mean(MatrixOperations.get_column_from_matrix(x, i));
			sd   = GeneralMath.sd(MatrixOperations.get_column_from_matrix(x, i));
			
			for(int j=0; j<nRows; j++){
				
				standardizedValues[j][i] = (x[j][i] - mean)/sd;
				
			}
					
		}
		
		return standardizedValues;
		
	}
	
	
	// test client
    public static void main(String[] args) {
    	
    	double[][] data = null;
    	
    	try {
    		data = ReadTextFile.readfile("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/RegTest.txt");
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	
    	double [][] y = MatrixOperations.get_sub_matrix_between_column_idxs(data, 4, 4);
    	double [][] X = MatrixOperations.get_sub_matrix_between_column_idxs(data, 0, 3);    	
    	
    	RidgeRegression obj_ridge = new RidgeRegression(y, X, 0.17, true);
    		
    	obj_ridge.do_parameter_estimation();
    	
    	double [] penalties = {0.0, 0.17, 0.2, 0.3, 0.4, 0.5};
    	
    	obj_ridge.generate_ridge_regression_path(penalties);
    	
    	MatrixOperations.print_matrix(unstandardized_errors);

    	System.out.println(F);
    	
    }
	
	
}
