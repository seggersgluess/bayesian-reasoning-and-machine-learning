package Optimization;
import java.util.function.BiFunction;

import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class NewtonMethod extends UnconstrainedOptimizer{

	static String version = "Newton";
	
	// constructor
	public NewtonMethod(BiFunction<double [], double [], Double> g, int max_iterations) {
		
		super(g, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, double [] further_args, int max_iterations) {
		
		super(g, further_args, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		super(g, grad, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations) {
		
		super(g, grad, hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations) {
		
		super(g, further_args, grad, hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations) {
		
		super(g, grad, further_args_grad, hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, int max_iterations) {
		
		super(g, further_args, grad, further_args_grad, hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations) {
		
		super(g, grad, hessian, further_args_hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations) {
		
		super(g, grad, further_args_grad, hessian, further_args_hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations) {
		
		super(g, further_args, grad, hessian, further_args_hessian, max_iterations);
	
	}
	
	
	public NewtonMethod(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, BiFunction <double [], double [], double [][]> hessian, double [] further_args_hessian, int max_iterations) {
		
		super(g, further_args, grad, further_args_grad, hessian, further_args_hessian, max_iterations);
	
	}
	

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
		public void do_Newton_Optimization(double [] start_value){

			int n_args = start_value.length;
			
			optimal_candidate = new double [n_args];
			
			double alpha = 1.0;
			
			double [] grad_1 = new double [n_args];
			double [][] hessian_matrix = new double [n_args][n_args];
			double [][] inv_hessian = new double [n_args][n_args];
			
			double [] direction = new double [n_args];
			
			double [] args_1 =  new double [n_args];
			double [] args_2 = new double [n_args];
			
			double [] further_args_4_gs;
			
			args_1 = start_value;
			
			for(int i = 0; i < max_number_of_iterations; i++){
					
				if(grad == null){
					
					grad_1      = NumDeriv.gradient(f, args_1, further_args_f);
					
				}else{
					
					grad_1 = gradient(args_1);
					
				}
				
				if(hessian == null){
					
					hessian_matrix   = NumDeriv.hessian(f, args_1, further_args_f);
					
				}else{
					
					hessian_matrix = hessian(args_1);
					
				}
					
				inv_hessian = MatrixOperations.inverse(hessian_matrix);
				
				for(int j=0; j<inv_hessian.length; j++){
					for(int k=0; k<inv_hessian.length; k++){
						if(Double.isInfinite(inv_hessian[j][k]) == true){
							inv_hessian[j][k] = Double.MAX_VALUE;
						}
					}
				}
				
				direction = MatrixOperations.scalar_vector_multiplication(-1.0, MatrixOperations.multiplyMatrixWithVec(inv_hessian, grad_1));
				
				if(version == "ModNewton"){
					
					GoldenSection.lower_bound = 0.0;
					GoldenSection.upper_bound = 5.0;
					
					further_args_4_gs = MatrixOperations.combine_vectors(args_1, direction);
					
				    if(further_args_f != null){
				    	
				    	further_args_4_gs = MatrixOperations.combine_vectors(further_args_4_gs, further_args_f);
				    	
				    }
					
					alpha = GoldenSection.do_Golden_Section_Optimization(alpha, NewtonMethod::get_golden_section_opt_res, further_args_4_gs);
							
					args_2    = MatrixOperations.add_vectors(args_1, MatrixOperations.scalar_vector_multiplication(alpha, direction));
										
				}
				
				if(version == "Newton"){
					
					args_2 = MatrixOperations.add_vectors(args_1, direction);
					
				}
									
				if(Math.abs(targetFunction(args_1, further_args_f) - targetFunction(args_2, further_args_f)) < convergence_criterion){
					
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
			optimal_value     = targetFunction(optimal_candidate, further_args_f);
			
		}
	
			
		// test client
	    public static void main(String[] args) {
	    		
	    	double [] further_args = {1.0, 100.0};
	    	
	    	NewtonMethod optim = new NewtonMethod(TargetFunction::target_function_with_further_args, further_args, 100000);
	    	
	    	double [] start_value = {-3000.05, -2000.1};
	        //double [] further_args = null;
	        
	        optim.set_version("ModNewton");
	        
	        //double [] solution = do_Newton_Optimization(start_value, TargetFunction::target_function, further_args, 100000);
	          
	    	optim.do_Newton_Optimization(start_value);
	                
	        MatrixOperations.print_vector(optimal_candidate);
	        	
	    }
		
}
