import java.util.function.BiFunction;

public class L_BFGS extends UnconstrainedOptimizer{

	// constructor
	public L_BFGS(BiFunction<double [], double [], Double> g, int max_iterations) {
		
		super(g, max_iterations);
	
	}
	
	
	public L_BFGS(BiFunction<double [], double [], Double> g, double [] further_args, int max_iterations) {
		
		super(g, further_args, max_iterations);
	
	}
	
	
	public L_BFGS(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		super(g, grad, max_iterations);
	
	}
	
	
	public L_BFGS(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		super(g, further_args, grad, max_iterations);
	
	}

	
	public L_BFGS(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, int max_iterations) {
		
		super(g, further_args, grad, further_args_grad, max_iterations);
	
	}
	
	
	// main (Quasi-Newton) optimization routine of the Broyden, Fletcher, Goldfarb and Shanno (BFGS) optimization routine
	public void do_LBFGS_Optimization(double [] start_value){
		
		int n_args = start_value.length;
		
		optimal_candidate = new double [n_args];
		
		double alpha = 1.0;
		
		double [] grad_1 = new double [n_args];
		double [] grad_2 = new double [n_args];
		double [] delta_grads = new double [n_args];
		
		double [][] inv_quasi_hessian = new double [n_args][n_args];
		double [] direction = new double [n_args];
		
		double [] args_1 =  new double [n_args];		
		double [] args_2 = new double [n_args];
		double [] delta_args = new double [n_args];
		
		double objective1;
		double objective2;
		
		double [] further_args_4_gs;
		
		args_1 = start_value;
		
		if(grad == null){
			
			grad_1 = NumDeriv.gradient(f, args_1, further_args_f);
			
		}else{
			
			grad_1 = gradient(args_1);
			
		}
			
		inv_quasi_hessian = MatrixOperations.identity(n_args);
		
		for(int i = 0; i < max_number_of_iterations; i++){
						
			if(i > 0){
				
				if(grad == null){
					
					grad_2 = NumDeriv.gradient(f, args_2, further_args_f);
					
				}else{
					
					grad_2 = gradient(args_2);
					
				}
				
				delta_grads = MatrixOperations.add_vectors(grad_2, MatrixOperations.scalar_vector_multiplication(-1.0, grad_1));
				delta_args  = MatrixOperations.add_vectors(args_2, MatrixOperations.scalar_vector_multiplication(-1.0, args_1));
				
				inv_quasi_hessian = calculate_inv_quasi_hessian(inv_quasi_hessian, delta_args, delta_grads);
				
				args_1 = args_2;
				grad_1 = grad_2;
				
			}
		
			//MatrixOperations.print_vector(grad_1);
			
			direction = MatrixOperations.scalar_vector_multiplication(-1.0, MatrixOperations.multiplyMatrixWithVec(inv_quasi_hessian, grad_1));
							
			GoldenSection.lower_bound = 0.0;
			GoldenSection.upper_bound = 5.0;
				
			further_args_4_gs = MatrixOperations.combine_vectors(args_1, direction);
			
		    if(further_args_f != null){
		    	
		    	further_args_4_gs = MatrixOperations.combine_vectors(further_args_4_gs, further_args_f);
		    	
		    }
			
			alpha = GoldenSection.do_Golden_Section_Optimization(alpha, L_BFGS::get_golden_section_opt_res, further_args_4_gs);
						
			args_2    = MatrixOperations.add_vectors(args_1, MatrixOperations.scalar_vector_multiplication(alpha, direction));
					
			objective1 = targetFunction(args_1, further_args_f);
			objective2 = targetFunction(args_2, further_args_f);
				
			System.out.println(objective2);
			
			if(Math.abs(objective1 - objective2) < convergence_criterion){
				
				System.out.println("L-BFGS Optimization converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}
				
			number_of_iterations = i;
			
		}
		
		optimal_candidate = args_2;
		optimal_value     = targetFunction(args_2, further_args_f);
			
	}
	
	
	// returns inverse of quasi Hessian
	public static double [][] calculate_inv_quasi_hessian(double [][] prev_inv_quasi_hessian, double [] delta_args, double [] delta_gradients){
		
		int n_args = delta_args.length;
		
		double [][] inv_quasi_hessian = new double [n_args][n_args];
		
		double [][] identity = MatrixOperations.identity(n_args);
		
		double [][] term_1;
		double [][] term_2;
		double [][] term_3;
		double denom;
		
		double [][] vec_delta_args = MatrixOperations.convArrayToVec(delta_args);
		double [][] vec_delta_gradients = MatrixOperations.convArrayToVec(delta_gradients);
		
		denom  = -1.0/MatrixOperations.multiplication(MatrixOperations.transpose(vec_delta_gradients), vec_delta_args)[0][0];
		
		term_1 = MatrixOperations.multiplication(vec_delta_args, MatrixOperations.transpose(vec_delta_gradients));		
		term_1 = MatrixOperations.scalar_multiplication(denom, term_1);
		term_1 = MatrixOperations.add(identity, term_1);
		
		
		term_2 = MatrixOperations.multiplication(vec_delta_gradients, MatrixOperations.transpose(vec_delta_args));	
		term_2 = MatrixOperations.scalar_multiplication(denom, term_2);
		term_2 = MatrixOperations.add(identity, term_2);
		
		term_3 = MatrixOperations.multiplication(vec_delta_args, MatrixOperations.transpose(vec_delta_args));		
		term_3 = MatrixOperations.scalar_multiplication((-1)*denom, term_3);
		
		inv_quasi_hessian = MatrixOperations.multiplication(term_1, prev_inv_quasi_hessian);
		inv_quasi_hessian = MatrixOperations.multiplication(inv_quasi_hessian, term_2);
		inv_quasi_hessian = MatrixOperations.add(inv_quasi_hessian,term_3);
		
		return inv_quasi_hessian;
		
	}
	
	
	// test client
    public static void main(String[] args) {
    	
    	double [] start_value = {3.05, -20.1};
    	double [] further_args = {1.0, 100.0};
    	
    	L_BFGS optim = new L_BFGS(TargetFunction::target_function_with_further_args, further_args, 100000);
    	//L_BFGS optim = new L_BFGS(TargetFunction::target_function, 100000);
    	
    	optim.do_LBFGS_Optimization(start_value);
        
        MatrixOperations.print_vector(optimal_candidate);
        	
    }
    
}
