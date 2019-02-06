import java.util.function.BiFunction;

public class UnconstrainedOptimizer {

	static double convergence_criterion;
	
	static double [] optimal_candidate;
	static double    optimal_value;
	static int       max_number_of_iterations;  //100000
	static int       number_of_iterations;
	static boolean   convergence;
	
	static double [][] candidate_trajectory;
	static double []   value_trajectory;
	
	static BiFunction <double [], double [], Double> f;
	
	public UnconstrainedOptimizer(int max_iterations){
			
		convergence_criterion    = 1e-15;
		
		optimal_candidate        = null;
		optimal_value            = 0;
		max_number_of_iterations = max_iterations;
		number_of_iterations     = 0;
		convergence              = false;
		
		candidate_trajectory     = null;
		value_trajectory         = null;
		
		f = null;
		
	}
	
	
	// computes target function
	public static double targetFunction(double [] opt_arg, double [] further_args){
		
		return f.apply(opt_arg, further_args);
		
	}
	
	
	// target function for golden section
	public static double get_golden_section_opt_res(double alpha, double [] further_args){
		
		int n_args = optimal_candidate.length;
		
		double [] prev_arg = MatrixOperations.get_double_sub_vec(further_args, 0, n_args-1);
		double [] direction = MatrixOperations.get_double_sub_vec(further_args, n_args, 2*n_args-1);
		double [] further_args_4_targ = MatrixOperations.get_double_sub_vec(further_args, 2*n_args, further_args.length-1);
			
		double [] args = new double [n_args];
		
		args    = MatrixOperations.add_vectors(prev_arg, MatrixOperations.scalar_vector_multiplication(alpha, direction));
		
		return targetFunction(args, further_args_4_targ);
		 
	}
	
	
}
