
import java.util.Arrays; 

public class BayesianRegression {

	static String prior;

	static boolean use_constant;
	
	//Observations
	static double [][] explained_variable;
	static double [][] explaining_variables;
	
	static int n_observations;
	static int n_explaining_variables;
	
	static int n_draws;
	
	//Estimated parameters & residuals
	
	static double [] est_confidence = {0.025, 0.975};
	
	static double  constant;
	static double  constant_error;
	static double  constant_t_value;
	static double [] constant_confidence;
	
	static double [][] parameters;
	static double [][] parameter_errors;
	static double [][] parameter_t_values;
	static double [][] parameter_confidence;
	
	static double [][] parameter_samples;
	static double [][] sigma_samples;
	
	static double [][] ss_pars;
	
	static double [][] fitted_explained_variable;
	static double [][] residuals;
	
	static double [][] covariance;
	
	static double ss_total;
	static double ss_reg;
	
	static double ss_res;
	static double sigma;
		
	static double F;
	
	static double R_squared;
	static double adj_R_squared;
	
	static double init_shape;
	static double init_scale;
	static double [][] init_my;
	static double [][] init_sigma;
	
	static double updated_shape;
	static double updated_scale;
	static double [][] updated_my;
	static double [][] updated_sigma;
	
	//------------------------------
	// Predictive Bayesian
	//------------------------------
		
	static double [] confidence = {0.025, 0.975};
	
	static double [][] new_explaining_variables;
	static double [][] explained_variable_sample;
	static double [][] predictive_mean;
	static double [][] predictive_sigma;
	static double [][] predictive_confidence;
	
	
	@SuppressWarnings("static-access")
	public BayesianRegression(double [][] y, double [][] X, boolean constant_usage, int n_draws){
		
		use_constant = constant_usage;
		
		explained_variable    = y;
		explaining_variables  = X;
		
		n_explaining_variables = X[0].length;
		n_observations         = y.length;
		
		this.n_draws = n_draws;
		
	}
	
	
	public void set_initial_hyperparameters_4_NIG(double shape, double scale, double [][] my, double [][] sigma){
		
		int d = 0;
		
	    if(use_constant == true){
	    	
	    	d = 1;
	    	
	    }
		
    	if(my.length != n_explaining_variables+d){
    		
    		throw new RuntimeException("Incorrect dimension of initial parameter mean.");
    		
    	}
    	
    	if(sigma.length != n_explaining_variables+d || sigma[0].length != n_explaining_variables+d){
    		
    		throw new RuntimeException("Incorrect dimension of initial parameter sigma.");
    		
    	}
    	
		init_shape = shape;
		//scales are here defined as rates!
	    init_scale = 1.0/scale; 	  
		init_my    = my;
		init_sigma = sigma;
	     
	}
	
	
	@SuppressWarnings({ "static-access" })
	public void do_bayesian_regression(){
		
		double [][] regressor_matrix = get_regressor_matrix();
		
	    int par_idx_1 = 0;
	    int par_idx_2 = n_explaining_variables-1;
		
		if(use_constant == true){
				
			par_idx_1 = 1;
			par_idx_2 = n_explaining_variables;
			
		}
		
		parameters = new double [n_explaining_variables][1];
		
		prepare_hyperparameters();
		
		NormalDistribution normalDist     = new NormalDistribution(updated_my, updated_sigma);
		InvGammaDistribution invGammaDist = new InvGammaDistribution(updated_shape, updated_scale);
		
		for(int i=0; i<n_draws; i++){
			
			if(i==0){
										
				sigma_samples     = invGammaDist.sample(1);
						
				normalDist.sigma  = MatrixOperations.scalar_multiplication(sigma_samples[0][i], updated_sigma);						
				parameter_samples = normalDist.sample(1);
				
			}else{
							
				sigma_samples     = MatrixOperations.cbind(sigma_samples, invGammaDist.sample(1));
				
				normalDist.sigma  = MatrixOperations.scalar_multiplication(sigma_samples[0][0], updated_sigma);		
				parameter_samples = MatrixOperations.cbind(parameter_samples, normalDist.sample(1));
						
			}
	
		}
		
		parameter_samples   = MatrixOperations.transpose(parameter_samples);
		
		parameters           = MatrixOperations.transpose(GeneralMath.mean_vec(parameter_samples));
		parameter_errors     = MatrixOperations.transpose(GeneralMath.sd_vec(parameter_samples));	
		
		parameter_confidence      = calculate_confidence_intervals(parameter_samples, est_confidence);
		parameter_confidence      = MatrixOperations.transpose(parameter_confidence);
		
		fitted_explained_variable = MatrixOperations.multiplication(regressor_matrix, parameters);
		residuals                 = MatrixOperations.add(explained_variable, MatrixOperations.scalar_multiplication(-1.0, fitted_explained_variable));
		
		if(use_constant == true){
			
			constant            = MatrixOperations.get_double_sub_vec(parameters, 0, 0)[0][0];
			constant_error      = MatrixOperations.get_double_sub_vec(parameter_errors, 0, 0)[0][0];	
			constant_t_value    = constant/constant_error;
			constant_confidence = MatrixOperations.get_column_from_matrix(parameter_confidence , 0);
			
			parameters           = MatrixOperations.get_double_sub_vec(parameters, par_idx_1, par_idx_2);
			parameter_errors     = MatrixOperations.get_double_sub_vec(parameter_errors, par_idx_1, par_idx_2);
			parameter_confidence = MatrixOperations.get_sub_matrix_between_column_idxs(parameter_confidence, par_idx_1, par_idx_2);
			
		}
		
		calculate_t_statistics_4_est_pars(); 
		
		sigma              = GeneralMath.mean(MatrixOperations.transpose(sigma_samples));
		
		print_estimation_res();
		
	}
	
	
	// samples from the posterior predictive distribution
	@SuppressWarnings("static-access")
	public static void do_bayesian_prediction(){
		
		if(new_explaining_variables == null){
			
			throw new RuntimeException("No new inputs set for Bayesian prediction.");
			
		}
		
		if(new_explaining_variables[0].length != n_explaining_variables){
			
			throw new RuntimeException("Incorrect dimension of new input. Check your supplied data.");
			
		}
		
		if(updated_my == null || updated_sigma == null){
			
			throw new RuntimeException("Bayesian regression not done. Do at first regression.");
			
		}
		
		new_explaining_variables  = get_regressor_matrix();
		
		double[][] student_sigma = null;
		double df = 0;
		
		int n_new_obs             = new_explaining_variables.length;
		
		explained_variable_sample = new double [n_draws][n_new_obs];
		
		double [][] identity      = MatrixOperations.identity(n_new_obs);
		
		double [][] student_my    = MatrixOperations.multiplication(new_explaining_variables, updated_my);
	    
		student_sigma = MatrixOperations.multiplication(MatrixOperations.multiplication(new_explaining_variables, updated_sigma),MatrixOperations.transpose(new_explaining_variables));
        student_sigma = MatrixOperations.add(identity, student_sigma);
		
		if(prior == "CONJUGATE"){
						
            student_sigma = MatrixOperations.scalar_multiplication(updated_scale/updated_shape, student_sigma);
            df            = 2.0*updated_shape;
			
		}

		if(prior == "UNINFORMATIVE"){
			
			student_sigma = MatrixOperations.scalar_multiplication(sigma, student_sigma);
            df            = n_observations-n_explaining_variables;
			
		}
		
		StudentDistribution student = new StudentDistribution(student_my, student_sigma, df);
		
		for(int i=0; i<n_draws; i++){
			
			double [][] sample = student.sample();
			
			for(int j=0; j<n_new_obs; j++){
				
				explained_variable_sample[i][j] = sample[j][0];
				
			}
					
		}
		
		predictive_mean       = GeneralMath.mean_vec(explained_variable_sample);		
		predictive_sigma      = GeneralMath.sd_vec(explained_variable_sample);
		predictive_confidence = calculate_confidence_intervals(explained_variable_sample, confidence);
		
	}
	
	
	// setter for new inputs X
	public static void set_new_input_4_bayesian_prediction(double [][] X){
		
		new_explaining_variables = X;
		
	}
	
	
	//calculates confidence intervals
	public static double [][] calculate_confidence_intervals(double [][] sample, double [] confidence){
		
		if(Utilities.getMax(confidence)>1.0 || Utilities.getMin(confidence)<0.0){
			
			throw new RuntimeException("Invalid confidence levels supplied. Levels only in [0.0, 1.0] allowed.");
			
		}
		
		if(sample.length != n_draws){
			
			sample = MatrixOperations.transpose(sample);
			
		}

		int n_obs = sample[0].length;
		int n_confidences = confidence.length;		
		double [] column;
		double [][] confidence_intervals = new double [n_obs][n_confidences];
		double idx;
		
		for(int i=0; i<n_obs; i++){
			
			column = MatrixOperations.get_column_from_matrix(sample, i);
	        Arrays.sort(column, 0, n_draws-1);
			
	        for(int j=0; j<n_confidences; j++){
	        	
	        	idx = confidence[j]*n_draws;
	        	if(Math.ceil(idx)-idx > 0.5){
	        		
	        		idx = idx-1;
	        		
	        	}
	        	
	        	if(idx >= 0.0){
	        		
	        		idx = Math.ceil(idx)-1.0;
	        		
	        	}else{
	        		
	        		idx = 0.0;
	        		
	        	}
	        	
	        	confidence_intervals[i][j] = column[(int) idx];
	        	
	        }
	        
		}
			
		return confidence_intervals;
		
	}
	
	
	// setter for the confidence levels
	public static void set_confidence(double [] confidence_levels){
		
		if(Utilities.getMax(confidence_levels)>1.0 || Utilities.getMin(confidence_levels)<0.0){
			
			throw new RuntimeException("Invalid confidence levels supplied. Levels only in [0.0, 1.0] allowed.");
			
		}
		
		confidence = confidence_levels;
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void prepare_hyperparameters(){
				
		if(prior == "CONJUGATE"){
			
			if(init_my == null || init_sigma == null){
				
				System.out.println("No hyperparameters for priors supplied. Use default parameters.");
				
				set_default_prior_hyperparameters();
				
			}
			
			double [][] term;
			double [][] init_sigma_inv = MatrixOperations.inverse(init_sigma);
			
			double[][] regressor_matrix = get_regressor_matrix();
			
			term = MatrixOperations.multiplication(MatrixOperations.transpose(regressor_matrix),regressor_matrix);
			
			updated_sigma = MatrixOperations.inverse(MatrixOperations.add(term, init_sigma_inv));
			
			updated_my    = MatrixOperations.inverse(MatrixOperations.add(term, init_sigma_inv));
			updated_my    = MatrixOperations.multiplication(updated_sigma, MatrixOperations.add(MatrixOperations.multiplication(init_sigma_inv, init_my),MatrixOperations.multiplication(MatrixOperations.transpose(regressor_matrix), explained_variable)));
					
			updated_shape = init_shape + n_observations/2.0;
			
			updated_scale = MatrixOperations.multiplication(MatrixOperations.transpose(explained_variable), explained_variable)[0][0];
			updated_scale = updated_scale + MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(init_my), init_sigma_inv), init_my)[0][0];
			updated_scale = updated_scale - MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(updated_my), MatrixOperations.inverse(updated_sigma)), updated_my)[0][0];
			updated_scale = init_scale + 0.5*updated_scale;
			//scales are here defined as rates!
			updated_scale = 1.0/updated_scale;
			
		}
		
		if(prior == "UNINFORMATIVE"){
			
			LinearRegression linearReg = new LinearRegression(explained_variable, explaining_variables, use_constant);
			
			linearReg.do_parameter_estimation();
			
			updated_scale = linearReg.ss_res/2.0;
			//scales are here defined as rates!
			updated_scale = 1.0/updated_scale;
			
			updated_shape = (n_observations-n_explaining_variables)/2.0;
			
			updated_my  = linearReg.parameters;
			
			if(use_constant == true){
				
				double [][] regConst = new double [1][1];
				
				regConst[0][0] = linearReg.constant;
				updated_my     = MatrixOperations.transpose(MatrixOperations.cbind(MatrixOperations.transpose(updated_my),regConst));
				
			}
			
			updated_sigma   = MatrixOperations.scalar_multiplication(1.0/linearReg.sigma, linearReg.covariance);
	
		}
			
	}
	
	
	// setter for default hyperparameters of conjugate (unit information) priors
	public static void set_default_prior_hyperparameters(){
		
		int d = 0;
		
		if(use_constant == true){
			
			d = 1;
			
		}
		
		init_scale = 0.0;
		init_shape = 0.0;
		init_my    = new double [n_explaining_variables+d][1];
		
		double [][] regressor_matrix = get_regressor_matrix();
		
		init_sigma = MatrixOperations.scalar_multiplication(n_observations, MatrixOperations.inverse(MatrixOperations.multiplication(MatrixOperations.transpose(regressor_matrix), regressor_matrix)));	
		
	}
	
	
	// setter for prior type
	public static void set_prior_type(String used_prior){
		
		int [] valid_prior = Utilities.get_idx(get_valid_priors(), used_prior);
		
		if(valid_prior[0] == -1){
			
			throw new RuntimeException(used_prior + " is not a valid prior for Bayesian Regression.");
			
		}
				
		prior = used_prior;
		
	}
	
	
	// returns names of valid priors
	public static String [] get_valid_priors(){
		
		String [] valid_priors = {
									"CONJUGATE",
				                  	"UNINFORMATIVE"
								 };
		
		return valid_priors;
		
	}
	
	
	// returns regressor matrix
	public static double [][] get_regressor_matrix(){
		
		double [][] regressor_matrix = explaining_variables;
		
		if(use_constant == true){
			
			double [][] unit_vector = MatrixOperations.unit_vector(n_observations);
			
			regressor_matrix        = MatrixOperations.cbind(regressor_matrix, unit_vector);
		
		}
		
		return regressor_matrix;
		
	}
	
	
	// returns the t-values for the estimated parameters
	public static void calculate_t_statistics_4_est_pars(){
		
		parameter_t_values = new double [n_explaining_variables][1];
		
		for(int i=0; i<n_explaining_variables; i++){
			
			parameter_t_values[i][0] = parameters[i][0]/parameter_errors[i][0];
			
		}
			
	}
	
	
	// prints estimated parameters, standard errors, t-values and p-values
	public static void print_estimation_res(){
		
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
		System.out.println("Bayesian Linear Model Estimation Results:");
		System.out.println("");
		
		MatrixOperations.print_matrix(est_res);
		
	}
	
	
	// test client
    @SuppressWarnings("static-access")
	public static void main(String[] args) {
    	
    	double[][] data = null;
    	
    	try {
    		data = ReadTextFile.readfile("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/RegTest.txt");
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
    	
    	double [][] y = MatrixOperations.get_sub_matrix_between_column_idxs(data, 4, 4);
    	double [][] X = MatrixOperations.get_sub_matrix_between_column_idxs(data, 0, 3);    	
    	
    	BayesianRegression obj_blm = new BayesianRegression(y, X, false, 100);
    	
    	//obj_blm.set_prior_type("UNINFORMATIVE");
    	
    	obj_blm.set_prior_type("CONJUGATE");
    	
    	obj_blm.do_bayesian_regression();
    	
    	obj_blm.set_new_input_4_bayesian_prediction(X);
    	
    	System.out.println("");
    	
    	if(use_constant == true){
    		
    		System.out.println("Constant 95% Confidence");
    		MatrixOperations.print_vector(constant_confidence);
    		
    	}
    	
    	System.out.println("");
    	System.out.println("Parameter 95% Confidence");
    	MatrixOperations.print_matrix(parameter_confidence);
    	
    	obj_blm.do_bayesian_prediction();
    	
    	System.out.println("");
    	System.out.println("Predictive Mean");
    	MatrixOperations.print_matrix(predictive_mean);
    	
    	System.out.println("");
    	System.out.println("95% Prediction confidence");
    	MatrixOperations.print_matrix(obj_blm.predictive_confidence);
    		
    }
	
}
