import java.util.Arrays;

public class LogisticRegression {

	static int max_iterations = 10000;
	
	static String optimizer = "L_BFGS";
	static boolean use_analytic_gradient = true;
	static double convergence_criterion = 1e-10;
	static double log_likelihood;
	static int number_of_iterations;
	static boolean convergence;
	
	static boolean use_L2_regularization = false;
	//Only inverse of V is stored under variable name V.
	static double [][] V;                      
	static double lambda;
	
	static boolean use_constant;
	
	static double [][] explained_variable;
	static double [][] explaining_variables;

	static int n_observations;
	static int n_explaining_variables;
	
	static int n_classes;
	static double [] classes;
	
	static double [][] Weights;
	static double [][] my;
	
	static double [][] Y;
	
	@SuppressWarnings("static-access")
	LogisticRegression(double [][] explained_variable, double [][] explaining_variables, boolean constant_usage){
		
		if(explained_variable.length != explaining_variables.length){
			
			throw new RuntimeException("Inconsistent number of observations between for supplied data.");
			
		}
				
		use_constant = constant_usage;
		
		n_observations = explained_variable.length;
		n_explaining_variables = explaining_variables[0].length;
		
		this.explained_variable   = explained_variable;
		
		if(use_constant == true){
			
			explaining_variables = MatrixOperations.cbind(explaining_variables, MatrixOperations.unit_vector(n_observations));
			
			n_explaining_variables = n_explaining_variables+1;
			
		}
		
		this.explaining_variables = MatrixOperations.transpose(explaining_variables);
	
		classes   = Arrays.stream(MatrixOperations.get_column_from_matrix(explained_variable,0)).distinct().toArray();
		
		n_classes = classes.length;
		
	}
	
	
	@SuppressWarnings("static-access")
	LogisticRegression(double [][] explained_variable, double [][] explaining_variables, boolean constant_usage, boolean use_L2_regularization){
		
		if(explained_variable.length != explaining_variables.length){
			
			throw new RuntimeException("Inconsistent number of observations between for supplied data.");
			
		}
		
		use_constant = constant_usage;
		
		n_observations = explained_variable.length;
		n_explaining_variables = explaining_variables[0].length;
		
		this.explained_variable   = explained_variable;
		
		if(use_constant == true){
			
			explaining_variables = MatrixOperations.cbind(explaining_variables, MatrixOperations.unit_vector(n_observations));
			
			n_explaining_variables = n_explaining_variables+1;
			
		}
		
		this.explaining_variables = MatrixOperations.transpose(explaining_variables);
	
		classes   = Arrays.stream(MatrixOperations.get_column_from_matrix(explained_variable,0)).distinct().toArray();
		
		n_classes = classes.length;
		
		this.use_L2_regularization = use_L2_regularization;
		
		if(use_L2_regularization == true){
			
			// V means inverse of the L2 reg. matrix.
			V = MatrixOperations.identity(n_explaining_variables);
			lambda = 1.0;
			
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public void do_logistic_regression(){
		
		Y = get_indicator_matrix(explained_variable);
			
		int n_start_values = (n_classes-1)*n_explaining_variables;
		
		double [] start_values = new double[n_start_values];
		
		for(int i=0; i<n_start_values; i++){
			
			start_values[i] = Math.random()+3.0;
			
		}
			
		
		
		if(optimizer == "L_BFGS"){
			
			L_BFGS optim;
			
			if(use_analytic_gradient == true){
				
				optim = new L_BFGS(LogisticRegression::logistic_reg_objective, LogisticRegression::get_logistic_reg_gradient, max_iterations);
				
			}else{
				
				optim = new L_BFGS(LogisticRegression::logistic_reg_objective, max_iterations);
				
			}
			
			optim.convergence_criterion = convergence_criterion;
			
			optim.do_LBFGS_Optimization(start_values);
			
			update_logistic_reg_parameters(optim.optimal_candidate);
			
			log_likelihood       = optim.optimal_value;			
			number_of_iterations = optim.number_of_iterations;
			convergence          = optim.convergence;
			
		}
		
		if(optimizer == "BFGS"){
			
			BFGS optim;
			
			if(use_analytic_gradient == true){
				
				optim = new BFGS(LogisticRegression::logistic_reg_objective, LogisticRegression::get_logistic_reg_gradient, max_iterations);
				
			}else{
				
				optim = new BFGS(LogisticRegression::logistic_reg_objective, max_iterations);
				
			}
			
			optim.convergence_criterion = convergence_criterion;
			
			optim.do_BFGS_Optimization(start_values);
			
			update_logistic_reg_parameters(optim.optimal_candidate);
			
			log_likelihood       = optim.optimal_value;			
			number_of_iterations = optim.number_of_iterations;
			convergence          = optim.convergence;
			
		}
		
		if(optimizer == "Levenberg_Marquardt"){
			
			LevenbergMarquardt optim;
			
			if(use_analytic_gradient == true){
				
				optim = new LevenbergMarquardt(LogisticRegression::logistic_reg_objective, LogisticRegression::get_logistic_reg_gradient, max_iterations);
				
			}else{
				
				optim = new LevenbergMarquardt(LogisticRegression::logistic_reg_objective, max_iterations);
				
			}
			
			optim.convergence_criterion = convergence_criterion;
			
			optim.do_Levenberg_Marquardt_Optimization(start_values);
			
			update_logistic_reg_parameters(optim.optimal_candidate);
			
			log_likelihood       = optim.optimal_value;			
			number_of_iterations = optim.number_of_iterations;
			convergence          = optim.convergence;
			
		}
		
		
		MatrixOperations.print_matrix(Weights);
		
		if(use_constant == true){
			
			n_explaining_variables = n_explaining_variables-1;
			explaining_variables = MatrixOperations.get_sub_matrix_between_column_idxs(MatrixOperations.transpose(explaining_variables), 1, n_explaining_variables);
			
		}
		
	}
	
	
	// returns the gradient of the logistic objective
	public static double [] get_logistic_reg_gradient(double [] weights, double [] further_args){
		
		update_logistic_reg_parameters(weights);
		
		double [][] grad = new double [(n_classes-1)*n_explaining_variables][1];
		
		double[][] weight_sum = null;
		
		double [][] D = MatrixOperations.substract(my,Y);
		
		for(int i=0; i<n_observations; i++){
			
			double [][] d = MatrixOperations.get_sub_matrix_between_column_idxs(D, i, i);
			
			double [][] x = MatrixOperations.get_sub_matrix_between_column_idxs(explaining_variables, i, i);
			double [][] k = MatrixOperations.kronecker(d, x);
 				
			grad = MatrixOperations.add(grad, k);
				
		}
			
		if(use_L2_regularization == true){
		
			weight_sum = new double [n_explaining_variables][1];
			
			for(int i=0; i<n_classes-1; i++){
			
				weight_sum = MatrixOperations.add(weight_sum, MatrixOperations.transpose(MatrixOperations.get_sub_matrix_between_row_idxs(Weights,i,i)));
			
			}
		
			grad = MatrixOperations.add(grad, MatrixOperations.multiplication(V, weight_sum));
			
		}
		
		double [] arrGrad = MatrixOperations.convVecToArray(grad);
		
		return arrGrad;
		
	}
	
	
	// returns the hessian of the logistic objective
	public static double [][] get_logistic_reg_hessian(double [] weights){
		
		update_logistic_reg_parameters(weights);
		
		double [][] hessian = new double [n_explaining_variables][n_classes];
		
		double [][] identity = null;
		
		for(int i=0; i<n_observations; i++){
			
			double [][] x = MatrixOperations.get_sub_matrix_between_column_idxs(explaining_variables, i, i);
						x = MatrixOperations.multiplication(x, MatrixOperations.transpose(x));
			
			double [][] m = MatrixOperations.get_sub_matrix_between_column_idxs(my, i, i);
			double [][] diag = MatrixOperations.diagonal(m);
						m = MatrixOperations.multiplication(m, MatrixOperations.transpose(m));
					
			double [][] M = MatrixOperations.substract(diag, m);
			
			hessian = MatrixOperations.add(hessian, MatrixOperations.kronecker(M, x));
			
		}
		
		if(use_L2_regularization == true){
			
			identity = MatrixOperations.identity((n_classes-1));			
			hessian  = MatrixOperations.add(hessian, MatrixOperations.kronecker(identity, V));
			
		}
		
		return hessian;
		
	}
	
	
	// returns the objective of the logistic regression
	public static double logistic_reg_objective(double [] weights, double [] further_args){
		
		update_logistic_reg_parameters(weights);
		
		double weight_sum  = 0.0;
		double [][] w_dist = null;
		
		double [][] WX = MatrixOperations.multiplication(Weights, explaining_variables);
			
		double log_likelihood = 0.0;
		
		double term_1   = 0.0;
		double exp_term = 0.0;
		
		for(int i=0; i<n_observations; i++){
			
			// weights for class C eq 0.0
			exp_term = 1.0;
			
			for(int j=0; j<n_classes-1;j++){
					
				exp_term = exp_term + Math.exp(WX[j][i]);
										
				term_1   = term_1 + Y[j][i]*WX[j][i];		
				
			}
			
			term_1 = term_1 - Math.log(exp_term);
				
		} 
		
		if(Double.isInfinite(log_likelihood) == true){
			
			log_likelihood = Double.MAX_VALUE;
			
		}
		
		//negative log likelihood
		log_likelihood = -term_1;
		
		if(use_L2_regularization == true){
			
			for(int i=0; i<n_classes-1; i++){
				
				w_dist = MatrixOperations.transpose(MatrixOperations.get_sub_matrix_between_row_idxs(Weights,i,i));
				w_dist = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(w_dist), V), w_dist);
				
				weight_sum = weight_sum + w_dist[0][0];
				
			}	
			
			log_likelihood = log_likelihood + lambda*weight_sum;
			
		}
		
		return log_likelihood;
		
	}
	
	
	// update weights, and my
	public static void update_logistic_reg_parameters(double [] weights){
		
		Weights = MatrixOperations.transpose(MatrixOperations.get_matrix_from_vec(weights, n_explaining_variables, n_classes-1));
	
		my      = soft_max(MatrixOperations.multiplication(Weights, explaining_variables));
				
	}
	

	// calculates soft max from a column vector x at index
	public static double [][] soft_max(double [][] X){
		
		int n_classes = X.length;
		int n_obs     = X[0].length;
		
		double [][] soft_max = new double [n_classes][n_obs];
		double exp;
		
		for(int i=0; i<n_obs; i++){
			
			double sum = 0;
			
			for(int j=0; j<n_classes; j++){
				
				exp = Math.exp(X[j][i]);
				sum = sum + exp;
				
				soft_max[j][i] = exp;
				
			}
			
			for(int j=0; j<n_classes; j++){
						
				soft_max[j][i] = soft_max[j][i]/sum;
				
			}
			
		}
		
		return soft_max;
		
	}
	
	
	// returns matrix Y with elements y(i,j) = 1 if y(i) = j else 0
	public static double [][] get_indicator_matrix(double [][] y){
		
		double [][] Y = new double [n_classes-1][n_observations];
		
		for(int i=0; i<n_observations; i++){
			
			for(int j=0; j<n_classes-1; j++){
				
				if(y[i][0] == classes[j]){
					
					Y[j][i] = 1.0;
					
				}
					
			}
			
		}
		
		return Y;
		
	}
	
	
	public static void set_L2_regularization_matrix(double [][] reg_matrix){
		
		int nrows = reg_matrix.length;
		int ncols = reg_matrix[0].length;
		
		if(nrows != ncols && nrows!= n_explaining_variables){
		
			throw new RuntimeException("No valid dimension of L2 regularization matrix. Check the matrix.");
			
		}
		
		V = MatrixOperations.inverse(reg_matrix);
		
	}
	
	
	public static void set_L2_lambda(double reg_lambda){
		
		lambda = reg_lambda;
		
	}
	
	
	// test client
	public static void main(String[] args){
		
    	double[][] data = null;
    	
    	try {
    		data = ReadTextFile.readfile("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/LogisticRegression/NewExample.txt");
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	
    	double [][] y = MatrixOperations.get_sub_matrix_between_column_idxs(data, 0, 0);
    	double [][] X = MatrixOperations.get_sub_matrix_between_column_idxs(data, 1, 2); 
	
    	LogisticRegression obj_logistic = new LogisticRegression(y, X, true, true);
	
    	obj_logistic.do_logistic_regression();
    	
	}
	
}
