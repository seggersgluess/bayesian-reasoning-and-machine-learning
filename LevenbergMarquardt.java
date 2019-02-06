import java.util.function.BiFunction;

public class LevenbergMarquardt extends UnconstrainedOptimizer{

	
	// constructor
	public LevenbergMarquardt(int max_iterations) {
		
		super(max_iterations);
	
	}
	

	// main optimization routine of the Levenberg-Marquardt optimization routine
	public void do_Levenberg_Marquardt_Optimization(double [] start_value, BiFunction<double [], double [], Double> g, double [] further_args){
		
		f = g;
		
		int n_args = start_value.length;
		
		double lambda =1.0;
		double [] grad = new double [n_args];
		double [][] hessian = new double [n_args][n_args];
		double [][] steepestDescMatrix = new double [n_args][n_args];
		double [][] inv_matrix = new double [n_args][n_args];
		double [] direction = new double [n_args];
		double [] args_1 =  new double [n_args];
		double [] args_2 = new double [n_args];
		
		double [][] identity = MatrixOperations.identity(n_args);
		
		args_1 = start_value;
		
		for(int i = 0; i < max_number_of_iterations; i++){
				
			grad      = NumDeriv.gradient(f, args_1, further_args);
			hessian   = NumDeriv.hessian(f, args_1, further_args);
			
			steepestDescMatrix = MatrixOperations.scalar_multiplication(lambda, identity);
			inv_matrix = MatrixOperations.add(hessian, steepestDescMatrix);
			
			inv_matrix = MatrixOperations.inverse(inv_matrix);
						
			direction = MatrixOperations.scalar_vector_multiplication(-1.0, MatrixOperations.multiplyMatrixWithVec(inv_matrix, grad));
					
			args_2 = MatrixOperations.add_vectors(args_1, direction);
				
			if(targetFunction(args_1, further_args) > targetFunction(args_2, further_args)){
				
				lambda = lambda/2.0;
				
			}else{
				
				lambda = lambda*2.0;
				
			}
						
			if(Math.abs(targetFunction(args_1, further_args) - targetFunction(args_2, further_args)) < convergence_criterion){
				
				System.out.println("Levenberg-Marquardt Optimization converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}else{
				
				number_of_iterations = i;
				args_1 = args_2;
				
			}
					
		}
		
		optimal_candidate = args_2;
		optimal_value     = targetFunction(args_2, further_args);
		
	}
	
	
	// test client
    public static void main(String[] args) {
    		    
    	LevenbergMarquardt optim = new LevenbergMarquardt(100000);
    	
    	double [] start_value = {3.05, 20.0};
        double [] further_args = null;
        
        optim.do_Levenberg_Marquardt_Optimization(start_value, TargetFunction::target_function, further_args);
         
        MatrixOperations.print_vector(optimal_candidate);
        	
    }
	
}
