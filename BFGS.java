import java.util.function.BiFunction;

public class BFGS extends UnconstrainedOptimizer{

	// constructor
	public BFGS(int max_iterations) {
		
		super(max_iterations);
	
	}
		
	
	// main (Quasi-Newton) optimization routine of the Broyden, Fletcher, Goldfarb and Shanno (BFGS) optimization routine
	public void do_BFGS_Optimization(double [] start_value, BiFunction<double [], double [], Double> g, double [] further_args){
	
		f = g;
		
		int n_args = start_value.length;
		
		optimal_candidate = new double [n_args];
		
		double alpha = 1.0;
		
		double [] grad_1 = new double [n_args];
		double [] grad_2 = new double [n_args];
		double [] delta_grads = new double [n_args];
		
		double [][] quasi_hessian = new double [n_args][n_args];
		double [][] inv_quasi_hessian = new double [n_args][n_args];
		
		double [] direction = new double [n_args];
		
		double [] args_1 =  new double [n_args];
		double [] args_2 = new double [n_args];
		double [] delta_args = new double [n_args];
		
		double [] further_args_4_gs;
		
		args_1 = start_value;
		grad_1 = NumDeriv.gradient(f, args_1, further_args);
		
		quasi_hessian = MatrixOperations.identity(n_args);
		inv_quasi_hessian = quasi_hessian;
		
		for(int i = 0; i < max_number_of_iterations; i++){
						
			if(i > 0){
				
				grad_2      = NumDeriv.gradient(f, args_2, further_args);
				
				delta_grads = MatrixOperations.add_vectors(grad_2, MatrixOperations.scalar_vector_multiplication(-1.0, grad_1));
				delta_args  = MatrixOperations.add_vectors(args_2, MatrixOperations.scalar_vector_multiplication(-1.0, args_1));
				
				quasi_hessian     = calculate_quasi_hessian(quasi_hessian, delta_args, delta_grads);
				inv_quasi_hessian = MatrixOperations.inverse(quasi_hessian);
				
				args_1 = args_2;
				grad_1 = grad_2;
				
			}
		
			direction = MatrixOperations.scalar_vector_multiplication(-1.0, MatrixOperations.multiplyMatrixWithVec(inv_quasi_hessian, grad_1));
							
			GoldenSection.lower_bound = 0.0;
			GoldenSection.upper_bound = 5.0;
				
			further_args_4_gs = MatrixOperations.combine_vectors(args_1, direction);
			
		    if(further_args != null){
		    	
		    	further_args_4_gs = MatrixOperations.combine_vectors(further_args_4_gs, further_args);
		    	
		    }
			
			alpha = GoldenSection.do_Golden_Section_Optimization(alpha, BFGS::get_golden_section_opt_res, further_args_4_gs);
						
			args_2    = MatrixOperations.add_vectors(args_1, MatrixOperations.scalar_vector_multiplication(alpha, direction));
																					
			if(Math.abs(targetFunction(args_1, further_args) - targetFunction(args_2, further_args)) < convergence_criterion){
				
				System.out.println("BFGS Optimization converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}
			
			number_of_iterations = i;
			
		}
		
		optimal_candidate = args_2;
		optimal_value     = targetFunction(args_2, further_args);
		
	}
	
	
	// returns quasi Hessian for BFGS
	public static double [][] calculate_quasi_hessian(double [][] prev_quasi_hessian, double [] delta_args, double [] delta_gradients){
		
		int n_args = delta_args.length;
		
		double [][] quasi_hessian = new double [n_args][n_args];
		double [][] term_2;
		double [][] term_3;
		double denom;
		
		double [][] vec_delta_args = MatrixOperations.convArrayToVec(delta_args);
		double [][] vec_delta_gradients = MatrixOperations.convArrayToVec(delta_gradients);
		
		term_2 = MatrixOperations.multiplication(vec_delta_gradients, MatrixOperations.transpose(vec_delta_gradients));
		denom  = 1.0/MatrixOperations.multiplication(MatrixOperations.transpose(vec_delta_gradients), vec_delta_args)[0][0];
		
		term_2 = MatrixOperations.scalar_multiplication(denom, term_2);
		
		term_3 = MatrixOperations.multiplication(MatrixOperations.multiplication(prev_quasi_hessian, vec_delta_args),MatrixOperations.transpose(MatrixOperations.multiplication(prev_quasi_hessian, vec_delta_args)));
		
		denom = (MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(vec_delta_args), prev_quasi_hessian), vec_delta_args))[0][0];
		denom = -1.0/denom;
		
		term_3 = MatrixOperations.scalar_multiplication(denom, term_3);
		
		quasi_hessian = MatrixOperations.add(MatrixOperations.add(prev_quasi_hessian,term_2),term_3);
		
		return quasi_hessian;
		
	}
	
		
	// test client
    public static void main(String[] args) {
    	
    	BFGS optim = new BFGS(100000);
    	
    	double [] start_value = {1.05, 1.1};
        double [] further_args = null;
        
        optim.do_BFGS_Optimization(start_value, TargetFunction::target_function, further_args);
          
    	//double [] further_args = {1.0, 100.0};
    	//optim.do_BFGS_Optimization(start_value, TargetFunction::target_function_with_further_args, further_args);
        
        
        MatrixOperations.print_vector(optimal_candidate);
        	
    }
    
}
