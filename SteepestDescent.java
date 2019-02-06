import java.util.function.BiFunction;

public class SteepestDescent extends UnconstrainedOptimizer{

	
	// constructor
	public SteepestDescent(int max_iterations) {
		
		super(max_iterations);
	
	}
	
	
	// main optimization routine of the steepest descent optimization routine
	public void do_Steepest_Descent_Optimization(double [] start_value, BiFunction<double [], double [], Double> g, double [] further_args){
		
		f = g;
		
		int n_args = start_value.length;
		
		optimal_candidate = new double [n_args];
		
		double alpha = 1.0;
		
		double [] grad = new double [n_args];
		double [] direction = new double [n_args];
		
		double [] args_1 =  new double [n_args];
		double [] args_2 = new double [n_args];
		
		double [] further_args_4_gs;
		
		args_1 = start_value;
		
		for(int i = 0; i < max_number_of_iterations; i++){
				
			grad      = NumDeriv.gradient(f, args_1, further_args);
			direction = MatrixOperations.scalar_vector_multiplication(-1.0, grad);
							
			GoldenSection.lower_bound = 0.0;
		    GoldenSection.upper_bound = 5.0;
			
		    further_args_4_gs = MatrixOperations.combine_vectors(args_1, direction);
			
		    if(further_args != null){
		    	
		    	further_args_4_gs = MatrixOperations.combine_vectors(further_args_4_gs, further_args);
		    	
		    }
		    
		    alpha = GoldenSection.do_Golden_Section_Optimization(alpha, SteepestDescent::get_golden_section_opt_res, further_args_4_gs);
					
			args_2    = MatrixOperations.add_vectors(args_1, MatrixOperations.scalar_vector_multiplication(alpha, direction));
				
			if(Math.abs(targetFunction(args_1, further_args) - targetFunction(args_2, further_args)) < convergence_criterion){
				
				System.out.println("Steepest Descent Optimization converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}else{
				
				number_of_iterations = i;
				args_1 = args_2;
				
			}
					
		}
		
		optimal_candidate = args_2;
		optimal_value     = targetFunction(optimal_candidate, further_args);
		
		//return args_2;
		
	}
	

	// test client
    public static void main(String[] args) {
    		
    	SteepestDescent optim = new SteepestDescent(100000);
    	
        double [] start_value = {1.8, 2.1};
  
    	double [] further_args = {1.0, 100.0};
    	optim.do_Steepest_Descent_Optimization(start_value, TargetFunction::target_function_with_further_args, further_args);
              
        MatrixOperations.print_vector(optimal_candidate);
                	
    }

}
