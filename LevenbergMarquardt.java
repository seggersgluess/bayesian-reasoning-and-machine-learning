import java.util.function.BiFunction;

public class LevenbergMarquardt extends UnconstrainedOptimizer{

	
	// constructor
	public LevenbergMarquardt(BiFunction<double [], double [], Double> g, int max_iterations) {
		
		super(g, max_iterations);
	
	}
	
	
	public LevenbergMarquardt(BiFunction<double [], double [], Double> g, double [] further_args, int max_iterations) {
		
		super(g, further_args, max_iterations);
	
	}
	
	
	public LevenbergMarquardt(BiFunction<double [], double [], Double> g, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		super(g, grad, max_iterations);
	
	}
	
	
	public LevenbergMarquardt(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, int max_iterations) {
		
		super(g, further_args, grad, max_iterations);
	
	}

	
	public LevenbergMarquardt(BiFunction<double [], double [], Double> g, double [] further_args, BiFunction <double [], double [], double []> grad, double [] further_args_grad, int max_iterations) {
		
		super(g, further_args, grad, further_args_grad, max_iterations);
	
	}
	

	// main optimization routine of the Levenberg-Marquardt optimization routine
	public void do_Levenberg_Marquardt_Optimization(double [] start_value){
				
		int n_args = start_value.length;
		
		double lambda                  =1.0;
		double [] grad_1               = new double [n_args];
		double [][] hessian_matrix     = new double [n_args][n_args];
		double [][] steepestDescMatrix = new double [n_args][n_args];
		double [][] inv_matrix         = new double [n_args][n_args];
		double [] direction            = new double [n_args];
		double [] args_1               =  new double [n_args];
		double [] args_2               = new double [n_args];
		
		double [][] identity = MatrixOperations.identity(n_args);
		
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
			
			steepestDescMatrix = MatrixOperations.scalar_multiplication(lambda, identity);
			inv_matrix = MatrixOperations.add(hessian_matrix, steepestDescMatrix);
			
			inv_matrix = MatrixOperations.inverse(inv_matrix);
						
			direction = MatrixOperations.scalar_vector_multiplication(-1.0, MatrixOperations.multiplyMatrixWithVec(inv_matrix, grad_1));
					
			args_2 = MatrixOperations.add_vectors(args_1, direction);
				
			if(targetFunction(args_1, further_args_f) > targetFunction(args_2, further_args_f)){
				
				lambda = lambda/2.0;
				
			}else{
				
				lambda = lambda*2.0;
				
			}
					
			if(Math.abs(targetFunction(args_1, further_args_f) - targetFunction(args_2, further_args_f)) < convergence_criterion){
				
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
		optimal_value     = targetFunction(args_2, further_args_f);
		
	}
	
	
	// test client
    public static void main(String[] args) {
    	
    	//double [] further_args = null;
    	
    	LevenbergMarquardt optim = new LevenbergMarquardt(TargetFunction::target_function, 100000);
    	
    	double [] start_value = {3.05, 20.0};
        
        optim.do_Levenberg_Marquardt_Optimization(start_value);
         
        MatrixOperations.print_vector(optimal_candidate);
        	
    }
	
}
