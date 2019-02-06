import java.util.function.BiFunction;

public class NelderMead extends UnconstrainedOptimizer{

	
	// constructor
	public NelderMead(int max_iterations) {
		
		super(max_iterations);
	
	}
	
	//Meta parameters of Nelder-Mead optimization routine:
	static double alpha = 1.0;
	static double gamma = 2.0;
	static double rho   = 0.5;
	static double sigma = 0.5;

	
	// returns potential candidates indices of sorted target function values
	public static int [] get_sorted_candidates_idxs(double [][] candidates, BiFunction<double [], double [], Double> g, double [] further_args){
		
		f = g;
		
		int n_candidates = candidates[0].length;
		double [] values =  new double [n_candidates];
				
		for(int i=0; i<n_candidates; i++){
			
			double [] d = MatrixOperations.get_column_from_matrix(candidates, i);
			
			values[i] = targetFunction(d, further_args);
			
		}
				
		int [] sorted_idxs = Utilities.get_idxs_for_sorted_vec(values);
		
		return sorted_idxs;
		
	}
	
	
	// calculates centroid (mean) from possible candidates x
	public static double [] get_candidates_centroid(double [][] x){
		
		int n_rows = x.length;
		
		double [] centroid = new double [n_rows];
		
		for(int i=0; i<n_rows; i++){
			
			double [] elements = MatrixOperations.get_row_from_matrix(x, i);
			
			centroid[i] = GeneralMath.mean(elements);
			
		}
		
		return centroid;
			
	}
	
	
	// returns the reflection calculated from the centroid and the worst candidate
	public static double [] get_reflection_candidate(double [] centroid, double [] worst_candidate, double alpha){
		
		double [] diff_term = MatrixOperations.add_vectors(centroid,  MatrixOperations.scalar_vector_multiplication(-1.0, worst_candidate));
		
		diff_term = MatrixOperations.scalar_vector_multiplication(alpha, diff_term);
		
		double [] reflection = MatrixOperations.add_vectors(centroid, diff_term);
		
		return reflection;
		
	}
	
	
	// returns the expansion calculated from the centroid and the reflection candidate
	public static double [] get_expansion_candidate(double [] centroid, double [] reflection, double gamma){
		
		double [] diff_term = MatrixOperations.add_vectors(reflection,  MatrixOperations.scalar_vector_multiplication(-1.0, centroid));
		
		diff_term = MatrixOperations.scalar_vector_multiplication(gamma, diff_term);
		
		double [] expansion = MatrixOperations.add_vectors(centroid, diff_term);
		
		return expansion;
		
	}
	
	
	// returns the contraction calculated from the centroid and the worst candidate
	public static double [] get_contraction_candidate(double [] centroid, double [] worst_candidate, double rho){
		
		double [] diff_term = MatrixOperations.add_vectors(worst_candidate,  MatrixOperations.scalar_vector_multiplication(-1.0, centroid));
		
		diff_term = MatrixOperations.scalar_vector_multiplication(rho, diff_term);
		
		double [] contraction = MatrixOperations.add_vectors(centroid, diff_term);
		
		return contraction;
		
	}
	
	
	// returns the shrinkage calculated from the candidates and the best candidate
	public static double [][] do_shrinkage(double [][] all_candidates, double [] best_candidate, double sigma){
		
		int n_candidates = all_candidates[0].length;
		int n_rows = all_candidates.length;
		
		double [][] new_candidates = MatrixOperations.matrix(n_rows, n_candidates);
		
		for(int i=0; i<n_candidates; i++){
			
			double [] candidate = MatrixOperations.get_column_from_matrix(all_candidates, i);
			
			double [] diff_term = MatrixOperations.add_vectors(candidate,  MatrixOperations.scalar_vector_multiplication(-1.0, best_candidate));
			
			diff_term = MatrixOperations.scalar_vector_multiplication(sigma, diff_term);
			
			candidate = MatrixOperations.add_vectors(best_candidate, diff_term);
			
			new_candidates = MatrixOperations.set_column_to_matrix(new_candidates, candidate, i);
			
		}
		
		return new_candidates;
		
	}
	
	
	//main optimization routine of the Nelder-Mead optimization routine
	public void do_Nelder_Mead_Optimization(double [][] start_candidates, BiFunction<double [], double [], Double> f, double [] further_args){
		
		System.out.println("Start Nelder Mead Optimization");
		
		int n_candidates = start_candidates[0].length;
		double [][] candidates = start_candidates;
		
		int [] sorted_candidates_idxs = get_sorted_candidates_idxs(candidates, f, further_args);
		
		int worst_candidate_idx     = sorted_candidates_idxs[(n_candidates-1)];
		int sec_worst_candidate_idx = sorted_candidates_idxs[(n_candidates-2)];
		int best_candidate_idx      = sorted_candidates_idxs[0];
		
		double [] worst_candidate     = MatrixOperations.get_column_from_matrix(candidates, worst_candidate_idx);
		double [] sec_worst_candidate = MatrixOperations.get_column_from_matrix(candidates, sec_worst_candidate_idx);
		double [] best_candidate      = MatrixOperations.get_column_from_matrix(candidates, best_candidate_idx);
		
		for(int i=0; i<max_number_of_iterations; i++){
						
			sorted_candidates_idxs = get_sorted_candidates_idxs(candidates, f, further_args);
	
			double [][] sorted_candidates = MatrixOperations.resort_matrix_columns(candidates, sorted_candidates_idxs);
			double [][] best_candidates = MatrixOperations.get_sub_matrix_between_column_idxs(sorted_candidates, 0, (n_candidates-2));
			
			double [] centroid = get_candidates_centroid(best_candidates);
			
			worst_candidate_idx     = sorted_candidates_idxs[(n_candidates-1)];
			sec_worst_candidate_idx = sorted_candidates_idxs[(n_candidates-2)];
			best_candidate_idx      = sorted_candidates_idxs[0];
			
			worst_candidate     = MatrixOperations.get_column_from_matrix(candidates, worst_candidate_idx);
			sec_worst_candidate = MatrixOperations.get_column_from_matrix(candidates, sec_worst_candidate_idx);
			best_candidate      = MatrixOperations.get_column_from_matrix(candidates, best_candidate_idx);
			
			double worst_value     = targetFunction(worst_candidate, further_args);
			double sec_worst_value = targetFunction(sec_worst_candidate, further_args);
			double best_value      = targetFunction(best_candidate, further_args);
			
			//reflection step
			double [] reflection = get_reflection_candidate(centroid, worst_candidate, alpha);
			
			double reflection_value = targetFunction(reflection, further_args);
			
			if(reflection_value >= best_value && reflection_value < sec_worst_value){
				
				candidates = MatrixOperations.set_column_to_matrix(candidates, reflection, worst_candidate_idx);
				
			}else{
				
				//expansion step
				if(reflection_value < best_value){
										
					double [] expansion = get_expansion_candidate(centroid, reflection, gamma);
					
					double expansion_value = targetFunction(expansion, further_args);
					
					if(expansion_value < reflection_value){
						
						candidates = MatrixOperations.set_column_to_matrix(candidates, expansion, worst_candidate_idx);
						
					}else{
						
						candidates = MatrixOperations.set_column_to_matrix(candidates, reflection, worst_candidate_idx);
						
					}
								
				}//end expansion step
				
				//contraction step
				if(reflection_value >= sec_worst_value){
								
					double [] contraction = get_contraction_candidate(centroid, worst_candidate, rho);
					
					double contraction_value = targetFunction(contraction, further_args);
					
					if(contraction_value < worst_value){
						
						candidates = MatrixOperations.set_column_to_matrix(candidates, contraction, worst_candidate_idx);
						
					}else{
						
						//shrinkage
						candidates = MatrixOperations.resort_matrix_columns(candidates, sorted_candidates_idxs);
						
						double [][] all_candidates = MatrixOperations.get_sub_matrix_between_column_idxs(candidates, 1, (n_candidates-1));
						candidates = do_shrinkage(all_candidates, best_candidate, sigma);
						candidates = MatrixOperations.cbind(candidates, best_candidate);
						
					}				
					
				}
							
			}
			
			if(Math.abs(best_value - worst_value) < convergence_criterion){
				
				System.out.println("Nelder Mead Optimization converged after " + i + " iterations");
				
				number_of_iterations = i;
				convergence          = true;
				
				break;
				
			}
			
			number_of_iterations = i;
			
		}
		
		optimal_candidate = best_candidate;
		optimal_value     = targetFunction(optimal_candidate, further_args);
			
	}
	
	
	// test client
    public static void main(String[] args){
    		
    	NelderMead optim = new NelderMead(100000);
    	
        double [][] start_candidates = {{10, 21.0, 81.0, 10.8, 4.8, -2.0}, {100.2, 21.05, 7.0, 10.8, 202.0, -11.5}};
        double [] further_args = null; //{1.0, 100.0};	
        
    	optim.do_Nelder_Mead_Optimization(start_candidates, TargetFunction::rastrigin_function, further_args);
        
    	MatrixOperations.print_vector(optimal_candidate);
               
    }
	
}
