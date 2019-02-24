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
	
	static double [] further_args_f;
	
	static BiFunction <double [], double [], double []> grad;
	
	static double [] further_args_grad;
	
	static BiFunction <double [], double [], double [][]> hessian;
	
	static double [] further_args_hessian;
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, int max_iterations){
			
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
				
		this.f = objective_function;
		
		this.grad    = null;
		this.hessian = null;
		
		this.further_args_f       = null;
		this.further_args_grad    = null;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, int max_iterations){
			
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		
		this.grad    = null;
		this.hessian = null;
		
		this.further_args_f       = further_args;
		this.further_args_grad    = null;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = null;
	
		this.further_args_f       = null;
		this.further_args_grad    = null;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = null;
	
		this.further_args_f       = further_args;
		this.further_args_grad    = null;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, int max_iterations) {
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = null;
	
		this.further_args_f       = further_args;
		this.further_args_grad    = further_args_grad;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = null;
		this.further_args_grad    = null;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = further_args;
		this.further_args_grad    = null;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = further_args;
		this.further_args_grad    = further_args_grad;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = null;
		this.further_args_grad    = further_args_grad;
		this.further_args_hessian = null;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = further_args;
		this.further_args_grad    = further_args_grad;
		this.further_args_hessian = further_args_hessian;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, double [] further_args, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = further_args;
		this.further_args_grad    = null;
		this.further_args_hessian = further_args_hessian;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = null;
		this.further_args_grad    = further_args_grad;
		this.further_args_hessian = further_args_hessian;
		
	}
	
	
	@SuppressWarnings("static-access")
	public UnconstrainedOptimizer(BiFunction<double [], double [], Double> objective_function, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations){
		
		this.convergence_criterion    = 1e-15;
		
		this.optimal_candidate        = null;
		this.optimal_value            = 0;
		this.max_number_of_iterations = max_iterations;
		this.number_of_iterations     = 0;
		this.convergence              = false;
		
		this.candidate_trajectory     = null;
		this.value_trajectory         = null;
		
		this.f = objective_function;
		this.grad = grad;
		this.hessian = hessian;
		
		this.further_args_f       = null;
		this.further_args_grad    = null;
		this.further_args_hessian = further_args_hessian;
		
	}
	
	
	// computes target function
	public static double targetFunction(double [] opt_arg, double [] further_args){
		
		return f.apply(opt_arg, further_args);
		
	}
	
	
	// computes gradient
	public static double [] gradient(double [] opt_arg){
		
		if(grad != null){
			
			return grad.apply(opt_arg, further_args_grad);
			
		}else{
			
			throw new RuntimeException("No function for gradient supplied.");
			
		}
		
	}
	
	
	// computes hessian
	public static double [][] hessian(double [] opt_arg){
		
		if(hessian != null){
			
			return hessian.apply(opt_arg, further_args_hessian);
			
		}else{
			
			throw new RuntimeException("No function for Hessian supplied.");
			
		}
		
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
