import java.util.function.BiFunction;

public class GoldenSection {

	static double tau = 2.0-(1.0 + Math.sqrt(5.0))/2.0;
	static double lower_bound;
	static double upper_bound;
	
	static double convergence_criterion = 1e-45;
	static double n_iterations = 1000;
	
	
	// computes target function
	public static double targetFunction(BiFunction<Double, double [], Double> f, double opt_arg, double [] further_args){
		
		return f.apply(opt_arg, further_args);
		
	}
	
	
	//main optimization routine of the golden section optimization routine
	public static double do_Golden_Section_Optimization(double x, BiFunction<Double, double [], Double> f, double [] further_args){
			
		double alpha_1 = lower_bound*(1.0-tau) + upper_bound*tau;
		double alpha_2 = lower_bound*tau + upper_bound*(1.0-tau);
		
		for(int i = 0; i < n_iterations; i++){
			
			if(targetFunction(f, alpha_1, further_args) > targetFunction(f, alpha_2, further_args)){
				
				lower_bound = alpha_1;
				alpha_1 = alpha_2;
				alpha_2 = lower_bound*tau + upper_bound*(1.0-tau);
				
			}else{
				
				upper_bound = alpha_2;
				alpha_2= alpha_1;
				alpha_1 = lower_bound*(1.0-tau) + upper_bound*tau;
				
			}
			
			if(Math.abs(targetFunction(f, alpha_1, further_args) - targetFunction(f, alpha_2, further_args)) < convergence_criterion){
				
				//System.out.println("Golden Section converged after " + i + " iterations");
				
				break;
				
			}
			
		}
		
		return alpha_1;
		
	}
	
	
	// test client
    public static void main(String[] args) {
    		
        //double start_value = 1.0;
        //double [] further_args = null;
        
        //double solution = do_Golden_Section_Optimization(start_value, TargetFunction.target_function_4_golden_section, further_args);
                	
        //System.out.println(solution);
       
    }
	
}
