import java.util.function.BiFunction;

public class NewtonMethod extends UnconstrainedOptimizer{

	// constructor
	public NewtonMethod(int max_iterations) {
		
		super(max_iterations);
	
	}
	
	static String version = "Newton";
	
	
	// sets version of Newton algorithm
	public void set_version(String s_version){
		
		int [] valid_idxs = Utilities.get_idx(valid_versions(), s_version);
		
		if(valid_idxs[0] == -1.0){
			
			throw new RuntimeException("No valid Newton method supplied.");
			
		}else{
			
			version = s_version;
			
		}
		
	}
	
	
	// returns valid versions of Newton algorithm
	public static String [] valid_versions(){
		
		String [] versions = new String [2];
		
		versions[0] = "Newton";
		versions[1] = "ModNewton";
				
		return versions;
		
	}
	
	
	// main optimization routine of the Newton optimization routine
		public void do_Newton_Optimization(double [] start_value, BiFunction<double [], double [], Double> g, double [] further_args){
			
			f = g;
			
			int n_args = start_value.length;
			
			optimal_candidate = new double [n_args];
			
			double alpha = 1.0;
			
			double [] grad = new double [n_args];
			double [][] hessian = new double [n_args][n_args];
			double [][] inv_hessian = new double [n_args][n_args];
			
			double [] direction = new double [n_args];
			
			double [] args_1 =  new double [n_args];
			double [] args_2 = new double [n_args];
			
			double [] further_args_4_gs;
			
			args_1 = start_value;
			
			for(int i = 0; i < max_number_of_iterations; i++){
					
				grad      = NumDeriv.gradient(f, args_1, further_args);
				hessian   = NumDeriv.hessian(f, args_1, further_args);
				inv_hessian = MatrixOperations.inverse(hessian);
				
				direction = MatrixOperations.scalar_vector_multiplication(-1.0, MatrixOperations.multiplyMatrixWithVec(inv_hessian, grad));
				
				if(version == "ModNewton"){
					
					GoldenSection.lower_bound = 0.0;
					GoldenSection.upper_bound = 5.0;
					
					further_args_4_gs = MatrixOperations.combine_vectors(args_1, direction);
					
				    if(further_args != null){
				    	
				    	further_args_4_gs = MatrixOperations.combine_vectors(further_args_4_gs, further_args);
				    	
				    }
					
					alpha = GoldenSection.do_Golden_Section_Optimization(alpha, NewtonMethod::get_golden_section_opt_res, further_args_4_gs);
							
					args_2    = MatrixOperations.add_vectors(args_1, MatrixOperations.scalar_vector_multiplication(alpha, direction));
										
				}
				
				if(version == "Newton"){
					
					args_2 = MatrixOperations.add_vectors(args_1, direction);
					
				}
											
				if(Math.abs(targetFunction(args_1, further_args) - targetFunction(args_2, further_args)) < convergence_criterion){
					
					System.out.println(version + " Optimization converged after " + i + " iterations");
					
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
			
		}
	
			
		// test client
	    public static void main(String[] args) {
	    		
	    	NewtonMethod optim = new NewtonMethod(100000);
	    	
	    	double [] start_value = {-3000.05, -2000.1};
	        //double [] further_args = null;
	        
	        optim.set_version("ModNewton");
	        
	        //double [] solution = do_Newton_Optimization(start_value, TargetFunction::target_function, further_args, 100000);
	          
	    	double [] further_args = {1.0, 100.0};
	    	optim.do_Newton_Optimization(start_value, TargetFunction::target_function_with_further_args, further_args);
	                
	        MatrixOperations.print_vector(optimal_candidate);
	        	
	    }
		
}
