import java.util.function.BiFunction;

public class DifferentialEvolution extends UnconstrainedOptimizer{

	
	// constructor
	public DifferentialEvolution(int max_iterations) {
		
		super(max_iterations);
	
	}
	
	// number of objective function evaluations until convergence
	static double evaluations = 100000;

	static double [] upper_bounds;
	static double [] lower_bounds;
		
	static int NP;
	static double CR = 0.9;
	static double F = 0.5;
		
	
	// main optimization routine of the differential evolution optimization routine
	public void do_Differential_Evolution_Optimization(double [] upper_values, double [] lower_values, BiFunction<double [], double [], Double> g, double [] further_args){
	
		if(upper_values.length != lower_values.length){
			
			throw new RuntimeException("Length of upper and lower values unequal.");
			
		}
		
		check_supplied_boundaries(upper_values, lower_values);
		
		upper_bounds = upper_values;
		lower_bounds = lower_values;
		
		f = g;
		
		int n_args = upper_bounds.length;
		
		NP = 10*n_args;
		
		double [][] prev_generation = get_initial_generation();
		double [][] new_generation = new double [n_args][NP];
		double f_1;
		double f_2;
				
		int flag = 0;
		
		for(int i = 0; i < max_number_of_iterations; i++){
			
			new_generation = mutation(prev_generation);		
			new_generation = correct_boundary_violation(new_generation);
			new_generation = crossover(prev_generation, new_generation);
			
			flag = flag + 1;
			
			for(int j = 0; j < NP; j++){
				
				f_1 = targetFunction(MatrixOperations.get_column_from_matrix(prev_generation, j),further_args);			
				f_2 = targetFunction(MatrixOperations.get_column_from_matrix(new_generation, j),further_args);
				
				if(f_2 < f_1){
					
					for(int k = 0; k < n_args; k++){
						
						prev_generation[k][j] = new_generation[k][j];
						
					}
					
					optimal_value = f_2;
					
					flag = 0;
					
				}else{
					
					optimal_value = f_1;
				}
				
				optimal_candidate = MatrixOperations.get_column_from_matrix(prev_generation, j);
				
			}
			
			if(flag == evaluations){
				
				System.out.println("Differential evolution converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}
					
			number_of_iterations = i;
			
		}
				
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
	
	
	// returns corrected generation x where all members are inside the feasible set
	public static double [][] correct_boundary_violation(double [][] x){
		
		int n_args = x.length;
		
		for(int i = 0; i < NP; i++){
			
			for(int j = 0; j < n_args; j++){
				
				if(x[j][i] < lower_bounds[j] || x[j][i] > upper_bounds[j]){
					
					x[j][i] = lower_bounds[j] + (upper_bounds[j] - lower_bounds[j])*Math.random();;
					
				}
				
			}
			
		}
		
		return x;
		
	}
	
	
	// generates random initial generation of NP population members from upper and lower bounds
	public static double [][] get_initial_generation(){
		
		int n_args = upper_bounds.length;
		
		double [][] init_generation = new double [n_args][NP];
			
		for(int i = 0; i < NP; i++){
			
			for(int j = 0; j < n_args; j++){
				
				init_generation[j][i] = lower_bounds[j] + (upper_bounds[j] - lower_bounds[j])*Math.random();
				
			}
				
		}
		
		return init_generation;
		
	}
	

	// generates population of mutants
	public static double [][] mutation(double [][] prev_generation){
		
		int n_args = prev_generation.length;
		
		double [][] mutants = new double [n_args][NP];
		
		int [] rand_idxs = new int [3];
		
		for(int i = 0; i < NP; i++){
						
			rand_idxs = Utilities.getRandomIntNumbers(0,(NP-1),3);
				
			while(Utilities.get_idx(rand_idxs, i)[0] == -1.0){
					
				rand_idxs = Utilities.getRandomIntNumbers(0,(NP-1),3);;
					
			}
						
			for(int j = 0; j < n_args; j++){
				
				mutants[j][i] = prev_generation[j][rand_idxs[0]] + F*(prev_generation[j][rand_idxs[1]] - prev_generation[j][rand_idxs[2]]);
				
			}
					
		}
		
		return mutants;
		
	}
	
	
	// makes crossover between old population and mutated population
	public static double [][] crossover(double [][] prev_generation, double [][] mutated_generation){

		int n_args = prev_generation.length;
		
		double [][] cross_over_generation = new double [n_args][NP];
		
		int rand_idx;
		double rand_number;
		
		for(int i = 0; i < NP; i++){
			
			rand_idx = Utilities.getRandomIntNumbers(0,(NP-1),1)[0];
			
			for(int j = 0; j < n_args; j++){
				
				rand_number = Math.random();
				
				if(rand_number <= CR || j == rand_idx){
					
					cross_over_generation[j][i] = mutated_generation[j][i];
					
				}else{
					
					cross_over_generation[j][i] = prev_generation[j][i];
					
				}
				
			}
			
		}
		
		return cross_over_generation;
		
	}	
	
	
	// test client
    public static void main(String[] args) {
    	
    	DifferentialEvolution optim = new DifferentialEvolution(10000000);
    	
        double [] upper_values = {5.12, 5.12};
    	double [] lower_values = {-5.12, -5.12};
    	double [] further_args = null;
        
    	optim.do_Differential_Evolution_Optimization(upper_values, lower_values, TargetFunction::rastrigin_function, further_args);
  
    	MatrixOperations.print_vector(optimal_candidate);
        System.out.println(number_of_iterations);
    	
    }
	
}
