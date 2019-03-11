
public class LogisticRegression {

	static int max_iterations = 10000;
	
	static String optimizer = "L_BFGS";
	static boolean use_analytic_gradient = true;
	static double convergence_criterion = 1e-10;
	static double log_likelihood;
	static double intercept_only_log_likelihood;
	static int number_of_iterations;
	static boolean convergence;
	
	static boolean use_L2_regularization = false;
	//Only inverse of diagonal V is stored under variable name V.
	static double [][] V;                      
	static double lambda;
	
	static boolean use_constant;
	
	static InputDataManager inputData;
	
	static double [][] explained_variable;
	static double [][] explaining_variables;

	static String [] names_of_explaining_variables;
	
	static int n_observations;
	static int n_explaining_variables;
	
	static String reference_class;
	static String [] str_classes;
	static double [] classes;
	
	static int n_classes;
	
	static double [][] Weights;
	static double [][] my;
	
	static double [][] Y;
	
	static double [][] standard_errors;
	static double [][] z_values;
	
	static double [][] rel_risk_Weights;
	
	//Pseudo R^2
	static double mcFadden;
	static double adj_mcFadden;
	static double cox_snell;
	static double cragg_uhler_nagelkerke;
	static double count;
	static double adj_count;
	
	
	LogisticRegression(boolean constant_usage){
			
		use_constant = constant_usage;
			
	}
	
	
	LogisticRegression(boolean constant_usage, boolean l2_regularization_usage){
		
		use_constant = constant_usage;
		use_L2_regularization =l2_regularization_usage;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void read_logistic_regression_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception{
		
		inputData = new InputDataManager();
		
		inputData.fileReader(fileName, false, hasRowNames, hasColNames);
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void select_logistic_regression_explained_variable(String [] rownames, String [] colnames, String ref_class){
		
		inputData.selectLoadedData(rownames, colnames);
		
		str_classes = Utilities.get_unique_elements(inputData.selectedStrFileData);
		n_classes   = str_classes.length;
		
		n_observations = inputData.selectedStrFileData.length;
		
		//Info: If no reference class is set, the last class will serve as reference class!
		if(ref_class != null){
			
			set_reference_class(ref_class);
			
		}
		
		convert_classes_into_numbers();
		
		for(int i=0; i<n_classes; i++){			
			classes[i] = (double) i;
		}
		
		Y = get_indicator_matrix(explained_variable);

		calc_log_lik_intercept_only();
		
		inputData.selectedStrFileData = null;
		
	}
	
	
	public static void set_reference_class(String ref_class){
		
		if(str_classes == null){
			throw new RuntimeException("Classes not set yet.");			
		}
		
		int idx = Utilities.get_idx(str_classes, ref_class)[0];
				
		if(idx == -1){
			throw new RuntimeException(ref_class + " not found as category in data for explained variable.");
		}
		
		reference_class = ref_class;
		
		//Rearrange classes with reference class at first position
		String [] new_str_classes = new String [n_classes];
		int counter = 0;
			
		idx = Utilities.get_idx(str_classes, reference_class)[0];
			
		new_str_classes[n_classes-1] = str_classes[idx];
			
		for(int i=0; i<n_classes; i++){
				
			if(i!=idx){
					
				new_str_classes[n_classes-counter-2] = str_classes[i];
				counter++;
					
			}
				
		}
			
		str_classes = new_str_classes;
				
	}
	
	
	@SuppressWarnings("static-access")
	public static void convert_classes_into_numbers(){
		
		classes = new double [n_classes];
		explained_variable = new double [n_observations][1];
		
		//Convert categories/classes into double numbers and set converted data as explained_variable
		for(int i=0; i<n_observations; i++){
			
			explained_variable[i][0] = (double) Utilities.get_idx(str_classes, inputData.selectedStrFileData[i][0])[0];
			
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void select_and_convert_logistic_regression_explaining_variables(String [] rownames, String [] colnames){
		
		inputData.selectLoadedData(rownames, colnames);
		
		int n_variables = inputData.selectedStrFileData[0].length;
		int n_obs       = inputData.selectedStrFileData.length;
		
		int n_category_variables = 0;
		int n_categories = 0;
	
		GenList categories = new GenList(n_variables);
		
		for(int i=0; i<n_variables; i++){
					
			try{
    			
				@SuppressWarnings("unused")
				double db     = Double.valueOf(inputData.selectedStrFileData[0][i]);
				double [] flag = null;
     			categories.add(flag);
    			
    		}catch(Exception e){			
	
    			int idx = 0;
    			
    			String [] unique_vec = new String [n_obs];		
    					
    			unique_vec[0] = inputData.selectedStrFileData[0][i];
    			
    			for(int j=1; j<n_obs; j++){
    				
    				int match = 0;
    				
    				for(int k=0; k<idx; k++){
    					
    					if(unique_vec[k].contentEquals(inputData.selectedStrFileData[j][i]) == true){
    						
    						match = 1;
    						break;
    						
    					}
    					
    				}
    				
    				if(match == 0){
    					
    					unique_vec[idx] = inputData.selectedStrFileData[j][i];
    					idx++;
    						
    				}
    				
    			}
    			
    			String [] unique_vec_2 = new String[idx];
    			
    			for(int j=0; j<idx; j++){
    				
    				unique_vec_2[j] = unique_vec[j];
    				
    			}
    			
    			categories.add(unique_vec_2);
    			
    			n_category_variables++;	
    			n_categories = n_categories + idx;
    			
    		}
				
		}
		
		int n = n_variables-n_category_variables+n_categories;
		int idx = 0;
		int const_idx = 0;
		boolean const_filled = false;
		
		if(use_constant == true){
			
			const_idx = 1;
			
		}
		
		//explaining_variables in transposed n_variables x n_obs form
		explaining_variables = new double [(n+const_idx)][n_obs];
		names_of_explaining_variables = new String [(n+const_idx)];
		
		if(use_constant == true){
			
			names_of_explaining_variables[0] = "Const";
			
		}
		
		for(int i=0; i<n_variables; i++){
			
			if(categories.get_string_array(i) == null){
				
				names_of_explaining_variables[(idx+const_idx)] = inputData.selected_colnames[i];
				
				for(int j=0; j<n_obs; j++){
					
					if(use_constant == true && const_filled == false){
						
						explaining_variables[0][j] = 1.0;
						
					}
					
					explaining_variables[(idx+const_idx)][j] = Double.valueOf(inputData.selectedStrFileData[j][i]);
					
				}
				
				const_filled = true;					
				
				idx++;
				
			}else{
				
				String [] variable_categories = categories.get_string_array(i);
				int n_variable_categories = variable_categories.length;
				
				for(int j=0; j<n_variable_categories; j++){
					
					names_of_explaining_variables[(idx+const_idx)] = variable_categories[j];
					
					for(int k=0; k<n_obs; k++){
						
						if(use_constant == true && const_filled == false){
							
							explaining_variables[0][k] = 1.0;
							
						}
						
						if(variable_categories[j].contentEquals(inputData.selectedStrFileData[k][i]) == true){
							
							explaining_variables[(idx+const_idx)][k] = 1.0;
							
						}	 
						
					}
					
					const_filled = true;
					
					idx++;
					
				}
					
			}
				
		}
			
		n_explaining_variables = explaining_variables.length;
		
		if(use_L2_regularization == true){
			
			// V means inverse of the L2 reg. matrix.
			V = MatrixOperations.identity(n_explaining_variables);
			lambda = 1.0;
			
		}
		
		inputData.selectedStrFileData = null;
		
	}
	
	
	public static void select_converted_logistic_regression_explaining_variables(String [] variable_names){
		
		if(explaining_variables == null){			
			throw new RuntimeException("No explaining variables selected yet.");
		}
		
		int [] variable_idxs = new int [variable_names.length];
		
		for(int i=0; i<variable_names.length; i++){
			
			variable_idxs[i] = Utilities.get_idx(names_of_explaining_variables, variable_names[i])[0];
			
			if(variable_idxs[i] == -1){
				throw new RuntimeException(variable_names[i] + " not found in selected data.");
			}
			
		}
		
		int const_idx = 0;
		
		if(use_constant == true){
			
			const_idx = 1;
			
		}
		
		int n_selected_vars = variable_idxs.length;
		
		double [][] selected_explaining_vars = new double [n_selected_vars+const_idx][n_observations];
		
		for(int i=0; i<n_observations; i++){
			
			if(use_constant == true){
				
				selected_explaining_vars[0][i] = 1.0;
				
			}
			
			for(int j=0; j<n_selected_vars; j++){
				
				selected_explaining_vars[j+const_idx][i] = explaining_variables[variable_idxs[j]][i];
				
			}
			
		}
		
		explaining_variables   = selected_explaining_vars;
		n_explaining_variables = explaining_variables.length;
		
		if(use_L2_regularization == true){
			
			// V means inverse of the L2 reg. matrix.
			V = MatrixOperations.identity(n_explaining_variables);
			lambda = 1.0;
			
		}
		
	}
	
	
	public static void remove_input_data_manager(){
		
		inputData = null;
		
	}
	
	
	public void do_logistic_regression(){
		
		//Check data consistency
		if(explained_variable == null){
			throw new RuntimeException("No explained variable selected for logistic regression.");
		}
		
		if(explaining_variables == null){
			throw new RuntimeException("No explaining variables selected for logistic regression.");
		}
		
		if(explained_variable.length != explaining_variables[0].length){
			throw new RuntimeException("Unequal number of observations for explained and explaining data for logistic regression supplied. Check your selected data.");		
		}
			
		int n_start_values = (n_classes-1)*n_explaining_variables;
		
		double [] start_values = new double[n_start_values];
		
		double []   opt_weights = null;
		
        opt_weights = do_MLE(start_values);
			
		calc_and_set_standard_errors(opt_weights);
		calc_z_values();	
		
		System.out.println("---------------------------------------------------------------------");
		
		System.out.println("Estimated weights:");
		MatrixOperations.print_matrix(Weights);
		
		System.out.println("Standard errors:");
		MatrixOperations.print_matrix(standard_errors);
		
		System.out.println("z-values:");
		MatrixOperations.print_matrix(z_values);
		
		System.out.println("---------------------------------------------------------------------");
		
		calc_and_set_pseudo_r_squred();
				
	}
	
	
	// returns the gradient of the logistic objective
	public static double [] get_logistic_reg_gradient(double [] weights, double [] further_args){
		
		update_logistic_reg_parameters(weights);
		
		double [][] grad = new double [(n_classes-1)*n_explaining_variables][1];
		
		double[][] l2_grad = null;
		
		double [][] D = MatrixOperations.substract(my,Y);
		
		for(int i=0; i<n_observations; i++){
			
			double [][] d = MatrixOperations.get_sub_matrix_between_column_idxs(D, i, i);
			
			double [][] x = MatrixOperations.get_sub_matrix_between_column_idxs(explaining_variables, i, i);
		
			double [][] k = MatrixOperations.kronecker(d, x);
			
			grad = MatrixOperations.add(grad, k);
				
		}

		if(use_L2_regularization == true){
		
			l2_grad = new double [n_explaining_variables*(n_classes-1)][1];
		
			int idx = 0;
			
			for(int i=0; i<n_classes-1; i++){
			
				for(int j=0; j<n_explaining_variables; j++){
									
					l2_grad[idx][0] = lambda*Weights[i][j]*V[j][j];
					
					idx = idx + 1;
					
				}
			
			}
		
			grad = MatrixOperations.add(grad, l2_grad);
			
		}
		
		double [] arrGrad = MatrixOperations.convVecToArray(grad);
		
		return arrGrad;
		
	}
	
	
	// returns the hessian of the logistic objective
	public static double [][] get_logistic_reg_hessian(double [] weights, double [] further_args){
		
		update_logistic_reg_parameters(weights);
		
		double [][] hessian = new double [n_explaining_variables*(n_classes-1)][n_explaining_variables*(n_classes-1)];
		
		double [][] l2_hessian = null;
		
		int idx = 0;
		
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
			
			l2_hessian = new double [n_explaining_variables*(n_classes-1)][n_explaining_variables*(n_classes-1)];		
			
			for(int i=0; i<n_explaining_variables*(n_classes-1); i++){
				
				if(idx == n_explaining_variables){
					
					idx = 0;
					
				}
				
				l2_hessian[i][i] = lambda*V[idx][idx];
				
				idx = idx+1;		
						
			}
			
			hessian  = MatrixOperations.add(hessian, l2_hessian);
			
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
			
			double sum = 1.0;
			
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
	
	
    //returns probability (odds) p(y=c|Wx) for observing y=c w.r.t. the estimated reg. model expressed as Wx. 
    //(x is column vector of explaining variables for which odds should be evaluated).
	public static double [][] get_odds(double [][] x){
		
		if(x.length != n_explaining_variables){			
			throw new RuntimeException("Invalid number of explaining variables supplied for calculating odds.");			
		}
		
		double [][] odds         = soft_max(MatrixOperations.multiplication(Weights, x));
	
		//Regarding probabilities for class 0 for which we have set w^T = 0	
		double [][] overall_odds = new double [n_classes][1];
		double      sum          = 0.0;
		
		for(int i=0; i<n_classes-1; i++){
				
			int idx = n_classes-2-i;
					
			overall_odds[(idx+1)][0] = odds[idx][0];
			sum                      = sum + odds[idx][0];
							
		}
		
		//probability of class 0
		overall_odds[0][0] = 1.0 - sum;
		
		return overall_odds;
		
	}
	
	
	//returns log odds log(p(y=c|Wx))
	public static double [][] get_log_odds(double [][] x){
			
		double [][] odds = get_odds(x);
		
		double [][] log_odds = new double [n_classes][1];
		
		for(int i=0; i<n_classes; i++){
			
			log_odds[i][0] = Math.log(odds[i][0]);
						
		}
		
		return log_odds;
		
	}
	
	
	public static double get_marginal_change_rate_of_log_odds(double [][] x, int class_idx, int x_idx){
		
		if(class_idx >= n_classes){			
			throw new RuntimeException("Invalid class index supplied. Check your index.");			
		}
		
		if(x_idx >= n_explaining_variables){			
			throw new RuntimeException("Invalid index for explaining variable supplied. Check your index.");			
		}
		
		double weight = 0.0;
		double marg_rates = 0.0;
		double denom = 0.0;
		
		double [][] odds = get_odds(x);
		
		double [] w_c = MatrixOperations.get_column_from_matrix(Weights, x_idx);
		
		for(int i=0; i<n_classes; i++){
			
			if(i==0){
				
				weight = 0.0;
				
			}else{
				
				weight =  w_c[(i-1)];
				
			}
			
			if(i != class_idx){
				
				marg_rates = marg_rates + weight*odds[i][0];
				denom      = denom + odds[i][0];
				
			}
			
		}
		
		if(class_idx == 0){
			
			marg_rates =  - marg_rates/denom;
			
		}else{
			
			marg_rates = w_c[class_idx-1] - marg_rates/denom;
			
		}
		
		return marg_rates;
		
	}
	
	
	public static double [] get_marginal_change_rate_of_log_odds_4_class_at_mean(int class_idx){
		
		double [][] mean = GeneralMath.mean_vec(explaining_variables);
	
		double [] marg_log_odd_changes = new double [n_explaining_variables];
		
    	double [][] x = new double [n_explaining_variables][1];
    	
    	for(int i=0; i<n_explaining_variables; i++){
    		
    		if(i==0){
    			
    			if(use_constant == true){
        			
        			x[i][0] = 1.0;
        			
        		}else{
        			
        			x[i][0] = mean[0][i];
        			
        		}
    			
    		}else{
    			
			if(use_constant == true){
        			
        			x[i][0] = mean[0][i-1];
        			
        		}else{
        			
        			x[i][0] = mean[0][i];
        			
        		}
    			
    		}
    							
    	}

    	for(int i=0; i<n_explaining_variables; i++){
    		
    		marg_log_odd_changes[i] = get_marginal_change_rate_of_log_odds(x, class_idx, i);
    		
    	}
    	
    	return marg_log_odd_changes;
		
	}
	
	
	public static void set_L2_regularization_matrix(double [][] reg_matrix){
		
		int nrows = reg_matrix.length;
		int ncols = reg_matrix[0].length;
		
		if(nrows != ncols && nrows!= n_explaining_variables){
		
			throw new RuntimeException("No valid dimension of L2 regularization matrix. Check the matrix.");
			
		}
		
		if(MatrixOperations.is_diagonal(reg_matrix) != true){
			
			throw new RuntimeException("Supplied regularization matrix is not diagonal. Check the matrix.");
			
		}
		
		V = MatrixOperations.inverse(reg_matrix);
		
	}
	
	
	public static void set_L2_lambda(double reg_lambda){
		
		lambda = reg_lambda;
		
	}
	
	
	public static double [] calc_standard_errors(double [] weights){
		
		double[][] hessian;
		
		//calculate standard errors 
		if(use_analytic_gradient == true){
			hessian = get_logistic_reg_hessian(weights, null);
		}else{
			hessian = NumDeriv.hessian(LogisticRegression::logistic_reg_objective, weights, null);		
		}
		
		double [][] inv_hessian = MatrixOperations.inverse(hessian);
		
		int n_args = inv_hessian.length;
		
		double [] se = new double [n_args];
	
		double [][] diag = MatrixOperations.get_diagonal_from_matrix(inv_hessian);
		
		for(int i=0; i<n_args; i++){
			
			se[i] = Math.sqrt(diag[i][0]);
			
		}
		
		return se;
			
	}
	
	
	public static void calc_and_set_standard_errors(double [] weights){
		
		double [] se = calc_standard_errors(weights);
		
		standard_errors = MatrixOperations.transpose(MatrixOperations.get_matrix_from_vec(se, n_explaining_variables, n_classes-1));
		
	}
	
	
	public static void calc_z_values(){
		
		z_values = new double [n_classes-1][n_explaining_variables];
		
		for(int i=0; i<n_classes-1; i++){
			
			for(int j=0; j<n_explaining_variables; j++){
				
				z_values[i][j] = Weights[i][j]/standard_errors[i][j];
				
			}
			
		}
		
	}
	
	
	public static void calc_and_set_relative_risk_ratios(){
		
		rel_risk_Weights         = new double [(n_classes-1)][n_explaining_variables];
		
		for(int i=0; i<n_classes-1; i++){
			
			for(int j=0; j<n_explaining_variables; j++){
				
				rel_risk_Weights[i][j] = Math.exp(Weights[i][j]);
				
			}
			
		}
					
	}
	
	
	@SuppressWarnings("static-access")
	public static double [] do_MLE(double [] start_values){
		
		double []   opt_weights = null;
		
		if(optimizer == "L_BFGS"){
			
			L_BFGS optim;
			
			if(use_analytic_gradient == true){
				
				optim = new L_BFGS(LogisticRegression::logistic_reg_objective, LogisticRegression::get_logistic_reg_gradient, max_iterations);
				
			}else{
				
				optim = new L_BFGS(LogisticRegression::logistic_reg_objective, max_iterations);
				
			}
			
			optim.convergence_criterion = convergence_criterion;
			
			optim.do_LBFGS_Optimization(start_values);			
			opt_weights = optim.optimal_candidate;
			
			update_logistic_reg_parameters(opt_weights);
			
			log_likelihood       = -optim.optimal_value;			
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
			
			opt_weights = optim.optimal_candidate;			
			update_logistic_reg_parameters(opt_weights);
			
			log_likelihood       = -optim.optimal_value;			
			number_of_iterations = optim.number_of_iterations;
			convergence          = optim.convergence;
			
		}
		
		if(optimizer == "Newton"){
			
			NewtonMethod optim;
			
			if(use_analytic_gradient == true){
				
				optim = new NewtonMethod(LogisticRegression::logistic_reg_objective, LogisticRegression::get_logistic_reg_gradient, LogisticRegression::get_logistic_reg_hessian, max_iterations);
				
			}else{
				
				optim = new NewtonMethod(LogisticRegression::logistic_reg_objective, max_iterations);
				
			}
			
			optim.set_version("ModNewton");
			
			optim.convergence_criterion = convergence_criterion;
			
			optim.do_Newton_Optimization(start_values);
			
			opt_weights = optim.optimal_candidate;			
			update_logistic_reg_parameters(opt_weights);
			
			log_likelihood       = -optim.optimal_value;			
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
			
			opt_weights = optim.optimal_candidate;			
			update_logistic_reg_parameters(opt_weights);
			
			log_likelihood       = -optim.optimal_value;			
			number_of_iterations = optim.number_of_iterations;
			convergence          = optim.convergence;
			
		}
		
		return opt_weights;
		
	}
	
	
	public static void set_optimizer(String selected_optimizer){
		
		int [] idx = Utilities.get_idx(get_allowed_optimizer(), selected_optimizer);
		
		if(idx[0] == -1){		
			throw new RuntimeException(selected_optimizer + " not an allowed optimizer for logistic regression.");			
		}
		
		optimizer = selected_optimizer;
		
	}
	
	
	public static String [] get_allowed_optimizer(){
		
		String [] optimizer = new String [4];
		
		optimizer[0] = "BFGS";
		optimizer[1] = "L_BFGS";
		optimizer[2] = "Levenberg_Marquardt";
		optimizer[3] = "Newton";
		
		return optimizer;		
				
	}
	
	
	public static void calc_log_lik_intercept_only(){
		
		boolean constant = false;
		
		if(use_constant == true){			
			constant = true;			
		}
		
		use_constant           = false;
		explaining_variables   = MatrixOperations.transpose(MatrixOperations.unit_vector(n_observations));
		n_explaining_variables = explaining_variables.length;
		
		if(use_L2_regularization == true){
			
			// V means inverse of the L2 reg. matrix.
			V = MatrixOperations.identity(n_explaining_variables);
			lambda = 1.0;
			
		}
		
		double [] start_values = new double [(n_classes-1)];
	    
		@SuppressWarnings("unused")
		double [] opt_weights  = do_MLE(start_values);
		
		intercept_only_log_likelihood = log_likelihood;
		
		use_constant = constant;
		
	}
	
	
	public static void calc_and_set_pseudo_r_squred(){
		
		double likelihood = Math.exp(log_likelihood);
		double intercept_only_likelihood = Math.exp(intercept_only_log_likelihood);
		
		//(log) likelihood based R^2
		mcFadden     = 1.0-log_likelihood/intercept_only_log_likelihood;
		adj_mcFadden = 1.0-(log_likelihood-n_explaining_variables)/intercept_only_log_likelihood; 
		cox_snell    = 1.0-Math.pow(intercept_only_likelihood/likelihood, 2.0/n_observations);
		cragg_uhler_nagelkerke = cox_snell/(1.0-Math.pow(intercept_only_likelihood, 2.0/n_observations));
		
		//count R^2
		double correct_count = count_correct_prediction_probs();
	    double max_freq      = count_most_frequent_class();
		count = correct_count/((n_classes-1)*n_observations);
		adj_count = (correct_count-max_freq)/((n_classes-1)*n_observations-max_freq);
		
	}
	
	
	public static double count_correct_prediction_probs(){
		
		double correct_pred_count = 0.0;

		for(int i=0; i<n_observations; i++){
			
			double obs_probs  = 0.0;
			double pred_probs = 0.0;
			
			for(int j=0; j<(n_classes-1); j++){
				
				obs_probs  = obs_probs + Y[j][i];
				pred_probs = pred_probs + my[j][i];
						
				if(Math.abs(Y[j][i]-my[j][i]) < 0.5){
						
					correct_pred_count = correct_pred_count + 1.0;
					
				}
							
			}
			
		}
			
		return correct_pred_count;
		
	}
	
	
	public static double count_most_frequent_class(){
		
		double [] counter = new double [n_classes];
		double max_count  = 0.0;
		double total = 0.0;
		
		for(int i=0; i<n_observations; i++){
			
			total = 0.0;
			
			for(int j=0; j<(n_classes-1.0); j++){
				
				total = total + Y[j][i];
				
				if(Y[j][i] == 1.0){
					
					counter[j] = counter[j] + 1.0;
					
				}
				
			}
			
			if(total == 0.0){
				
				counter[n_classes-1] = counter[n_classes-1] + 1.0;
				
			}
			
		}
		
		max_count = Utilities.getMax(counter);
		
		return max_count;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void unit_test_1() throws Exception{
		
		//Test with constant (true) and without l2 regularization (false)
		//See also: https://stats.idre.ucla.edu/stata/dae/multinomiallogistic-regression/
		
		//weights:
		double [][] w = {{5.218200299048554, 0.2913929856027784, -0.9826702373754551, -0.11360264235085257}, 
		                {2.8521862154952227, -0.5332909870086205, -1.1628320541752124, -0.057928409203357054}}; 
		//standard errors:
		double [][] se = {{1.163548973614939, 0.4763737074077733, 0.5955669294431913, 0.02221989917194932}, 
						  {1.1664391362907864, 0.443732131769415, 0.5142194604529419, 0.02141094600620421}};
		//z-values:
		double [][] z = {{4.484727688630533, 0.6116899003272395, -1.649974484470075, -5.112653368574506}, 
						 {2.445207921062235, -1.201830899380856, -2.261353650736887, -2.7055511319570487}}; 

		double intercept_log_lik = -204.09667423199426;
		
		String fileName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/LogisticRegression/Data.txt";
		
    	LogisticRegression obj_logistic = new LogisticRegression(true, false);
    	
    	obj_logistic.read_logistic_regression_input_data(fileName, true, true);
    	
    	String [] selected_y = {"prog"};
    	String    ref_class  = "academic";
     	
    	obj_logistic.select_logistic_regression_explained_variable(null, selected_y, ref_class);
    	
    	String [] selected_X = {"ses", "write"};
    	
    	obj_logistic.select_and_convert_logistic_regression_explaining_variables(null, selected_X);
    	
    	String [] selected_conv_X = {"middle", "high", "write"};
    	
    	obj_logistic.select_converted_logistic_regression_explaining_variables(selected_conv_X);
    	    	
    	obj_logistic.do_logistic_regression();
			
    	System.out.println("Delta est. weights:");
    	MatrixOperations.print_matrix(MatrixOperations.substract(Weights, w));
    	System.out.println("Delta standard errors:");
        MatrixOperations.print_matrix(MatrixOperations.substract(standard_errors, se));
    	System.out.println("Delta z-values:");
        MatrixOperations.print_matrix(MatrixOperations.substract(z_values, z));
    	System.out.println("Delta intercept log likelihood:");
    	System.out.println(intercept_log_lik-intercept_only_log_likelihood);
    	
	}
	
	
	// test client
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception{
				
		String fileName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/LogisticRegression/Data.txt";
		
    	LogisticRegression obj_logistic = new LogisticRegression(true, false);
    	
    	obj_logistic.read_logistic_regression_input_data(fileName, true, true);
    	
    	String [] selected_y = {"prog"};
    	String    ref_class  = "academic";
     	
    	obj_logistic.select_logistic_regression_explained_variable(null, selected_y, ref_class);
    	
    	String [] selected_X = {"ses", "write"};
    	
    	obj_logistic.select_and_convert_logistic_regression_explaining_variables(null, selected_X);
    	
    	String [] selected_conv_X = {"middle", "high", "write"};
    	
    	obj_logistic.select_converted_logistic_regression_explaining_variables(selected_conv_X);
    	    	
    	obj_logistic.do_logistic_regression();
    	
    	System.out.println(str_classes[0]);
    	MatrixOperations.print_vector(get_marginal_change_rate_of_log_odds_4_class_at_mean(0));
    	System.out.println(str_classes[1]);
    	MatrixOperations.print_vector(get_marginal_change_rate_of_log_odds_4_class_at_mean(1));
    	System.out.println(str_classes[2]);
    	MatrixOperations.print_vector(get_marginal_change_rate_of_log_odds_4_class_at_mean(2));
    	
	}
	
}
