import java.util.function.BiFunction;

public class ConjugateGradient extends UnconstrainedOptimizer{

	
	// constructor
	public ConjugateGradient(int max_iterations) {
		
		super(max_iterations);
	
	}


	// main optimization routine of the Fletcher Reeves conjugate gradient optimization routine
	public void do_Conjugate_Gradient_Optimization(double [] start_value, BiFunction<double [], double [], Double> g, double [] further_args){
		
		f = g;
		
		int n_args = start_value.length;
		
		optimal_candidate = new double [n_args];
		
		double alpha_1 = 1.0;
		double alpha_2 = 1.0;
		double lambda = 1.0;
		double [] grad_1 = new double [n_args];
		double [] grad_2 = new double [n_args];
		double [] direction = new double [n_args];
		double [] args_1 =  new double [n_args];
		double [] args_2 = new double [n_args];
		double [] args_3 = new double [n_args];
		
		double [] further_args_4_gs_1;
		double [] further_args_4_gs_2;
		
		args_1 = start_value;
		
		for(int i = 0; i < max_number_of_iterations; i++){
				
			grad_1    = NumDeriv.gradient(f, args_1, further_args);
			direction = MatrixOperations.scalar_vector_multiplication(-1.0, grad_1);
							
			GoldenSection.lower_bound = 0.0;
		    GoldenSection.upper_bound = 5.0;
					    	
		    further_args_4_gs_1 = MatrixOperations.combine_vectors(args_1, direction);
		    	
		    if(further_args != null){
		    	
		    	further_args_4_gs_1 = MatrixOperations.combine_vectors(further_args_4_gs_1, further_args);
		    	
		    }
		    	    
		    alpha_1   = GoldenSection.do_Golden_Section_Optimization(alpha_1, ConjugateGradient::get_golden_section_opt_res, further_args_4_gs_1);
			
			args_2    = MatrixOperations.add_vectors(args_1, MatrixOperations.scalar_vector_multiplication(alpha_1, direction));
			
			grad_2    = NumDeriv.gradient(f,  args_2, further_args);
			grad_2    = MatrixOperations.scalar_vector_multiplication(-1.0, grad_2);
			lambda    = Math.pow(MatrixOperations.euclidian(grad_2), 2.0)/Math.pow(MatrixOperations.euclidian(grad_1), 2.0);
			
			direction = MatrixOperations.add_vectors(grad_2, MatrixOperations.scalar_vector_multiplication(lambda, direction));

			further_args_4_gs_2 = MatrixOperations.combine_vectors(args_2, direction);
			
			if(further_args != null){
			    	
				further_args_4_gs_2 = MatrixOperations.combine_vectors(further_args_4_gs_2, further_args);
			    	
			}
			
			alpha_2 = GoldenSection.do_Golden_Section_Optimization(alpha_2, ConjugateGradient::get_golden_section_opt_res, further_args_4_gs_2);
	
			args_3    = MatrixOperations.add_vectors(args_2, MatrixOperations.scalar_vector_multiplication(alpha_2, direction));
				
			if(Math.abs(targetFunction(args_3, further_args) - targetFunction(args_2, further_args)) < convergence_criterion){
				
				System.out.println("Conjugate gradient Optimization converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}else{
				
				number_of_iterations = i;
				args_1 = args_2;
				
			}
					
		}
			
		optimal_candidate = args_3;
		optimal_value     = targetFunction(optimal_candidate, further_args);
			
	}
	
		
	// test client
    public static void main(String[] args) {
    	
    	ConjugateGradient cg1 = new ConjugateGradient(100000);
    	
    	double [] start_value = {3.05, -2.1};
        
    	double [] further_args = {1.0, 100.0};
    	cg1.do_Conjugate_Gradient_Optimization(start_value, TargetFunction::target_function_with_further_args, further_args);
        
        MatrixOperations.print_vector(optimal_candidate);
        System.out.println(number_of_iterations);
        System.out.println(convergence);
        System.out.println(convergence_criterion); 
        
        ConjugateGradient cg2 = new ConjugateGradient(100000);
        
        double [] sec_problem_further_args = null;
        cg2.do_Conjugate_Gradient_Optimization(start_value, TargetFunction::target_function, sec_problem_further_args); 
        
        MatrixOperations.print_vector(optimal_candidate);
        System.out.println(number_of_iterations);
        System.out.println(convergence);       
        System.out.println(convergence_criterion);
        
    }
    
}
