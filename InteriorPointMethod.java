
import java.util.function.BiFunction;

public class InteriorPointMethod {

	//Implements the Mehrotra predictor-corrector algorithm
	//As a Primal-Dual Interior Point Method
	
	static double convergence_criterion = 1e-08;
	
	static double [][] opt_candidate;
	
	static int number_of_constraints = 0;
	
	// Vector of coefficients in objective function f(x) = c^T x
	static double [][] c;
	
	// Coefficient matrix A of constraint Ax - b = 0
	static double [][] A;
	
	// Vector of constraint Ax - b = 0
	static double [][] b;
	
	// Equality and inequality (function) values
	static double [][] constraint_values;
		
	// Lambda
	static double [][] lambda;
	
	// s
	static double [][] s;
	
	// control parameter eta
	static double eta = 0.9;
			
	// Objective function
	static BiFunction <double [], double [], Double> f;
	
	
	// computes target function
	public static double targetFunction(double [] opt_arg, double [] further_args){
		
		return f.apply(opt_arg, further_args);
		
	}
	
				
	// main optimization routine of interior point optimization routine
	public static void do_Interior_Point_Optimization(double [] start_value, BiFunction<double [], double [], Double> objective_function, double [] further_args, double [][] c_vec, double [][] A_matrix, double [][] b_vec, int n_iterations){
		
		f = objective_function;
				
		A = A_matrix;
		b = b_vec;
		c = c_vec;
		
		number_of_constraints = A.length;
			
		int n_args         = start_value.length;
		
		if(n_args != A[0].length){
			
			throw new RuntimeException("Invalid dimensions.");
			
		}
		
		if(number_of_constraints != b.length){
			
			throw new RuntimeException("Invalid dimensions.");
			
		}
		
		opt_candidate = new double [n_args][1];
		
		//ToDos: How to determine lambda and s?
		//See Paper!!!
		lambda = new double [number_of_constraints][1];
		s      = MatrixOperations.unit_vector(n_args);
		
		double [][] B               = new double [2*n_args+number_of_constraints][(2*n_args+number_of_constraints)];
		double [][] delta_vec       = new double [(2*n_args+number_of_constraints)][1];
		double [][] correction_vec  = new double [(2*n_args+number_of_constraints)][1];
		
		double [][] r_c             = new double [n_args][1];
		double [][] r_b             = new double [number_of_constraints][1];
		double [][] r_xs            = new double [n_args][1];
		
		double [][] X               = new double [n_args][n_args];
		double [][] S               = new double [n_args][n_args];
		
		double [][] delta_X         = new double [n_args][n_args];
		double [][] delta_S         = new double [n_args][n_args]; 
		double [][] delta_xs        = new double [n_args][n_args]; 
		
		double [][] e               = MatrixOperations.unit_vector(n_args);
		double [][] identity_matrix = MatrixOperations.identity(n_args);
		
		double [][] x_1        = new double [n_args][1];		
		double [][] x_2        = new double [n_args][1];
		
		double [][] s_1        = new double [n_args][1];
		double [][] s_2        = new double [n_args][1];
		
		double [][] lambda_1   = new double [number_of_constraints][1];
		double [][] lambda_2   = new double [number_of_constraints][1];
		
		double [][] delta_x      = new double [n_args][1];
		double [][] delta_lambda = new double [number_of_constraints][1];
		double [][] delta_s      = new double [n_args][1];
		
		double alpha_p;
		double alpha_d;
		
		double [] pars;
				
		for(int i = 0; i < n_iterations; i++){
			
			X = MatrixOperations.diagonal(x_1);
			S = MatrixOperations.diagonal(s_1);
			
			B = MatrixOperations.set_sub_matrix_to_matrix(B, MatrixOperations.transpose(A), 0, n_args-1, n_args-1, n_args+number_of_constraints-2);
			B = MatrixOperations.set_sub_matrix_to_matrix(B, identity_matrix, 0, n_args-1, n_args+number_of_constraints-1, 2*n_args+number_of_constraints-1);
			B = MatrixOperations.set_sub_matrix_to_matrix(B, A, n_args-1, n_args+number_of_constraints-1, 0, n_args+number_of_constraints-1);
			B = MatrixOperations.set_sub_matrix_to_matrix(B, S, n_args+number_of_constraints-1, 2*n_args+number_of_constraints-1, 0, n_args-1);
			B = MatrixOperations.set_sub_matrix_to_matrix(B, X, n_args+number_of_constraints-1, 2*n_args+number_of_constraints-1, n_args+number_of_constraints-1, 2*n_args+number_of_constraints-1);
			
			r_c  = MatrixOperations.add(MatrixOperations.add(MatrixOperations.multiplication(MatrixOperations.transpose(A), lambda_1),s_1), MatrixOperations.scalar_multiplication(-1.0, c));
			r_b  = MatrixOperations.add(MatrixOperations.multiplication(A, x_1), MatrixOperations.scalar_multiplication(-1.0, b));
			r_xs = MatrixOperations.multiplication(MatrixOperations.multiplication(X, S),e);
				
			r_c  = MatrixOperations.scalar_multiplication(-1.0, r_c);
			r_b  = MatrixOperations.scalar_multiplication(-1.0, r_b);
			r_xs = MatrixOperations.scalar_multiplication(-1.0, r_xs);
						
			correction_vec = MatrixOperations.set_sub_matrix_to_matrix(correction_vec, r_c, 0, n_args-1, 0, 0);
			correction_vec = MatrixOperations.set_sub_matrix_to_matrix(correction_vec, r_b, n_args, n_args+number_of_constraints-1, 0, 0);
			correction_vec = MatrixOperations.set_sub_matrix_to_matrix(correction_vec, r_xs, n_args+number_of_constraints, 2*n_args+number_of_constraints-1, 0, 0);
			
			delta_vec = MatrixOperations.multiplication(MatrixOperations.inverse(B), correction_vec);
			
			delta_x = MatrixOperations.get_double_sub_vec(delta_vec, 0, n_args-1);
			delta_s = MatrixOperations.get_double_sub_vec(delta_vec, n_args+number_of_constraints-1, 2*n_args+number_of_constraints-1);
			
			alpha_p = get_step_length(x_1, delta_x);
			alpha_d = get_step_length(s_1, delta_s);
			
			pars    = get_duality_and_centering_parameters(x_1, delta_x, alpha_p, s_1, delta_s, alpha_d);
			
			delta_X = MatrixOperations.diagonal(delta_x);
			delta_S = MatrixOperations.diagonal(delta_s);
			
			delta_xs = MatrixOperations.multiplication(MatrixOperations.multiplication(delta_X, delta_S),e);
			delta_xs = MatrixOperations.scalar_multiplication(-1.0, delta_xs);
			
			r_xs     = MatrixOperations.add(MatrixOperations.add(r_xs, delta_xs),MatrixOperations.scalar_multiplication(pars[0]*pars[1], e));
			
			correction_vec = MatrixOperations.set_sub_matrix_to_matrix(correction_vec, r_xs, n_args+number_of_constraints, 2*n_args+number_of_constraints-1, 0, 0);
			
			delta_vec = MatrixOperations.multiplication(MatrixOperations.inverse(B), correction_vec);
			
			alpha_p = get_step_length(x_1, delta_x);
			alpha_d = get_step_length(s_1, delta_s);
			
			alpha_p = Math.min(1.0, alpha_p);
			alpha_d = Math.min(1.0, alpha_d);
			
			delta_lambda = MatrixOperations.get_double_sub_vec(delta_vec, n_args, n_args+number_of_constraints-1);
			
			x_2      = MatrixOperations.add(x_1, MatrixOperations.scalar_multiplication(alpha_p, delta_x)); 
			lambda_2 = MatrixOperations.add(lambda_1, MatrixOperations.scalar_multiplication(alpha_d, delta_lambda)); 
			s_2      = MatrixOperations.add(s_1, MatrixOperations.scalar_multiplication(alpha_d, delta_lambda)); 
			
			double relative_duality_gap = evaluate_relative_duality_gap(x_2, lambda_2);
			
			if(relative_duality_gap <= convergence_criterion){
				
				System.out.println("Interior Point Optimization converged after " + i + " iterations");
				
				break;		
				
			}else{
				
				x_1      = x_2;
				lambda_1 = lambda_2;
				s_1      = s_2;
				
			}
			
		}
		
		opt_candidate = x_2;
        lambda        = lambda_2;
        s             = s_2;
		
	}
	
	
	// returns step length alpha = min[1, min(-arg/delta_arg)]
	public static double get_step_length(double [][] arg, double [][] delta_arg){
		
		int    n_args = arg.length;
		double c      = 1.0;
		double alpha  = 1.0;
		
		for(int i = 0; i < n_args; i++){
			
			if(delta_arg[i][1] < 0.0){
				
				c = -arg[i][1]/delta_arg[i][1];
				
			}
			
			if(c < 1.0){
				
				alpha = c;
				
			}
			
		}
		
		return alpha;
		
	}
	
	
	// returns duality measure my and centering parameter sigma
	public static double [] get_duality_and_centering_parameters(double [][] x, double [][] delta_x, double alpha_p, double [][] s, double [][] delta_s, double alpha_d){
		
		int n_args = x.length;
		
		double [] parameters = new double [2];
		
		double [][] term_1;
		double [][] term_2;
		
		term_1 = MatrixOperations.add(x, MatrixOperations.scalar_multiplication(alpha_p, delta_x));
		term_2 = MatrixOperations.add(s, MatrixOperations.scalar_multiplication(alpha_d, delta_s));
		
		parameters[0] = MatrixOperations.multiplication(MatrixOperations.transpose(term_1), term_2)[0][0]/n_args;
		
		parameters[1] = Math.pow(parameters[0]/MatrixOperations.multiplication(MatrixOperations.transpose(x), s)[0][0]/n_args, 3.0);
		
		return parameters;
		
	}
	
	
	// returns relative duality gap (c^T x - b^T lambda)/(1 + |b^T lambda|)
	public static double evaluate_relative_duality_gap(double [][] x, double [][] lambda){
		
		double[][] denominator = MatrixOperations.multiplication(MatrixOperations.transpose(b), lambda);
		double[][] nominator   = MatrixOperations.add(MatrixOperations.multiplication(MatrixOperations.transpose(c), x), MatrixOperations.scalar_multiplication(-1.0, denominator));
		
		double rel_duality_gap = nominator[0][0]/Math.abs(denominator[0][0]);
		
		return rel_duality_gap;
		
	}
	
	
	// test client
	public static void main(String[] args){
		
		
	}
	
}
