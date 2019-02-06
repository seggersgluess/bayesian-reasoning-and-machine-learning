import java.util.function.BiFunction;

public class SimulatedAnnealing extends UnconstrainedOptimizer{

	// constructor
	public SimulatedAnnealing(int max_iterations) {
		
		super(max_iterations);
	
	}
	
	static double evaluations = 100000;
	static double epsilon = 20.0;
	static double [] upper_bounds;
	static double [] lower_bounds;
	

	// main optimization routine of the simulated annealing optimization routine
	public void do_Simulated_Annealing_Optimization(double [] upper_values, double [] lower_values, BiFunction<double [], double [], Double> g, double [] further_args){

		if(upper_values.length != lower_values.length){
			
			throw new RuntimeException("Length of upper and lower values unequal.");
			
		}
		
		check_supplied_boundaries(upper_values, lower_values);
		
		upper_bounds = upper_values;
		lower_bounds = lower_values;
		
		f = g;
		
		int n_args = upper_bounds.length;
		
		double [] args_1 = new double [n_args];
		double [] args_2 = new double [n_args];
		
		double [] perturbations = new double [n_args];
		
		double energy_old = 0.0;
		double energy_new = 0.0;
		
		double [] best_candidate;
		double best_value;
		
		int flag = 0;
		
		args_1 = get_random_values_between_bounds();
		
		energy_old = targetFunction(args_1, further_args);
		
		best_candidate = args_1;
		best_value     = energy_old;
		
		for(int i = 0; i < max_number_of_iterations; i++){
			
			perturbations = get_perturbations(args_1);
			
			args_2 = MatrixOperations.add_vectors(args_1, perturbations);
			
			args_2 = correct_boundary_violation(args_2);
			
			energy_new = targetFunction(args_2, further_args);
			
			if(energy_new < energy_old){
				
				args_1     = args_2;
				energy_old = energy_new;
							
			}else{
				
				double rel_delta_energy = (-1.0)*(energy_new - energy_old)/energy_old;
				
				if(Math.exp(rel_delta_energy) > Math.random()){
					
					args_1     = args_2;
					energy_old = energy_new;
					
				}
				
			}
			
			if(energy_old < best_value){
				
				best_candidate = args_1;
				best_value = energy_old;
							
				flag = 0;
				
			}else{
				
				flag = flag + 1;
				
			}
											
			if(flag == evaluations){
				
				System.out.println("Simulated annealing converged after " + i + " iterations");
							
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}
			
			number_of_iterations = i;
			
		}
				
		
		optimal_candidate    = best_candidate;
		optimal_value        = best_value;
				
	}
	
	
	// checks supplied boundaries
	public static void check_supplied_boundaries(double [] upper_bounds, double [] lower_bounds){
		
		int n_args = upper_bounds.length;
		
		for(int i = 0; i < n_args; i++){
			
			if(upper_bounds[i] <= lower_bounds[i]){
				
				throw new RuntimeException("Upper bounds lower than lower bounds.");
				
			}
			
		}
		
	}
	
	
	// generates random values between upper and lower bound
	public static double [] get_random_values_between_bounds(){
		
		int n_args = upper_bounds.length;
		
		double [] rand_values = new double [n_args];
			
		for(int i = 0; i < n_args; i++){
			
			rand_values[i] = lower_bounds[i] + (upper_bounds[i] - lower_bounds[i])*Math.random();
			
		}
		
		return rand_values;
		
	}
	
	
	// generates (random) perturbations of x
	public static double [] get_perturbations(double [] x){
		
		int n_args = x.length;
		
		double [] perturbations = new double [n_args];
			
		for(int i = 0; i < n_args; i++){
			
			perturbations[i] = epsilon*x[i]*Math.random();
							
		}
		
		return perturbations;
		
	}
	
	
	// returns false if boundaries are violated / else true.
	public static double [] correct_boundary_violation(double [] x){
		
		int n_args = x.length;
		
		for(int i = 0; i < n_args; i++){
			
			if(x[i] < lower_bounds[i] || x[i] > upper_bounds[i]){
				
				x[i] = lower_bounds[i] + (upper_bounds[i] - lower_bounds[i])*Math.random();;
				
			}
			
		}
		
		return x;
		
	}
	
	
	// test client
    public static void main(String[] args) {
    		
    	SimulatedAnnealing optim = new SimulatedAnnealing(10000000);
    	
        double [] upper_values = {5.12, 5.12};
        double [] lower_values = {-5.12, -5.12};
        double [] further_args = null;
        
        optim.do_Simulated_Annealing_Optimization(upper_values, lower_values, TargetFunction::rastrigin_function, further_args);
        //double [] solution = do_Simulated_Annealing_Optimization(upper_values, lower_values, TargetFunction::target_function, further_args, 10000000); 
        
        MatrixOperations.print_vector(optimal_candidate);
                	
    }
	
	
}
