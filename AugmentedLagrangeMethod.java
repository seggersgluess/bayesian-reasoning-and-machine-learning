import java.util.ArrayList;
import java.util.function.BiFunction;

public class AugmentedLagrangeMethod {

	static double convergence_criterion = 1e-15;
	
	static String unconstrained_optimizer = "L_BFGS";
	
	static double [] opt_candidate;
	
	// Growth factor for penalty parameter at every parameter update
	static double penalty_growth_factor = 1.0;
	
	static int number_of_equality_constraints = 0;
	static int number_of_inequality_constraints = 0;
	
	// Equality and inequality constraints
	static ArrayList <BiFunction <double [], double [], Double>> equality_constraints;
	static ArrayList <BiFunction <double [], double [], Double>> inequality_constraints;
	
	// Equality and inequality (function) values
	static double [][] equality_values;
	static double [][] inequality_values;
	
	// Further arguments of equality and inequality constraints
	static GenList further_args_eq_constraints;
	static GenList further_args_ineq_constraints;
	
	// Lagrangian multipliers for equality and inequality constraints
	static double [][] equality_multiplier;
	static double [][] inequality_multiplier;
	
	// Penalty parameter c
	static double penalty = 1.0;
	
	// Objective function
	static BiFunction <double [], double [], Double> f;
	
	// Augmented Lagrangian function
	static BiFunction <double [], double [], Double> augmentedLagrangian;
	
	// computes target function
	public static double targetFunction(double [] opt_arg, double [] further_args){
		
		return f.apply(opt_arg, further_args);
		
	}
	
	
	// returns vector of constraint values at point opt_args
	public static double [][] evaluate_constraints(ArrayList <BiFunction<double [], double [], Double>> constraints, double [] opt_args, GenList list_of_further_args){
		
		int n_constraints = constraints.size();
		
		if(list_of_further_args != null){
			
			if(n_constraints != list_of_further_args.filledSlots){
				
				throw new RuntimeException("Mismatch between number of constraints and its arguments.");
				
			}
			
		}
			
		double [] constraint_values = new double [n_constraints];
		
		for(int i = 0; i < n_constraints; i++){
			
			double [] further_args = null;
			
			if(list_of_further_args != null){
				
				further_args = list_of_further_args.get_double_array(i);
				
			}
						
			constraint_values[i] = constraints.get(i).apply(opt_args, further_args);
			
		}
		
		return MatrixOperations.convArrayToVec(constraint_values);
		
	}
	
	
	// sets ArrayList of equality constraints (functions)
	public static void set_equality_constraints(ArrayList <BiFunction<double [], double [], Double>> constraints){
		
		equality_constraints = constraints;
		
		number_of_equality_constraints = constraints.size();
		
	}
	
	
	// sets ArrayList of inequality constraints (functions)
	public static void set_inequality_constraints(ArrayList <BiFunction<double [], double [], Double>> constraints){
		
		inequality_constraints = constraints;
		
		number_of_inequality_constraints = constraints.size();
		
	}
	
	
	// sets the objective function
	public static void set_objectiv_function(BiFunction<double [], double [], Double> objective_function){
		
		f = objective_function;
		
	}
	
	
	// main optimization routine of augmented Lagrangian optimization routine
	public static double [] do_Augmented_Lagrange_Optimization(double [] start_value, BiFunction<double [], double [], Double> objective_function, double [] further_args, ArrayList <BiFunction<double [], double [], Double>> equality_constraints, ArrayList <BiFunction<double [], double [], Double>> inequality_constraints, GenList further_args_eq_const, GenList further_args_ineq_const, int n_iterations){

		String [] allowed_unconst_optimizer = get_unconstrained_optimizer_4_ALM();
		
		int [] valid_unconst_optimizer = Utilities.get_idx(allowed_unconst_optimizer, unconstrained_optimizer);
		
		if(valid_unconst_optimizer[0] == -1.0){
			
			throw new RuntimeException(unconstrained_optimizer + " is not a valid optimizer for ALM.");
			
		}
		
		set_objectiv_function(objective_function);
		
		if(equality_constraints != null){
			
			set_equality_constraints(equality_constraints);
			
		}
		
		if(inequality_constraints != null){
			
			set_inequality_constraints(inequality_constraints);
			
		}
		
		if(further_args_eq_const != null){
			
			further_args_eq_constraints = further_args_eq_const;
			
		}
		
		if(further_args_ineq_const != null){
			
			further_args_ineq_constraints = further_args_ineq_const;
			
		}
		
		int n_args         = start_value.length;
		
		opt_candidate = new double [n_args];
		
		double [] args_1        =  new double [n_args];		
		double [] args_2        = new double [n_args];
				
		if(number_of_equality_constraints != 0 && number_of_inequality_constraints != 0){
			
			augmentedLagrangian = AugmentedLagrangeMethod::augmented_Lagrangian_with_equality_and_inequality_constraints;
			
			equality_multiplier = new double [number_of_equality_constraints][1];
			inequality_multiplier = new double [number_of_inequality_constraints][1];
			
		}
		
		if(number_of_equality_constraints != 0 && number_of_inequality_constraints == 0){
			
			augmentedLagrangian = AugmentedLagrangeMethod::augmented_Lagrangian_with_equality_constraints;
			
			equality_multiplier = new double [number_of_equality_constraints][1];
			
		}
		
		if(number_of_equality_constraints == 0 && number_of_inequality_constraints != 0){
			
			augmentedLagrangian = AugmentedLagrangeMethod::augmented_Lagrangian_with_inequality_constraints;
			
			inequality_multiplier = new double [number_of_inequality_constraints][1];
			
		}
		
		args_1 = start_value;
		
		double augmented_prev = augmentedLagrangian.apply(args_1, further_args);
		
		for(int i = 0; i < n_iterations; i++){
						
			if(i > 0){
					
				args_1 = args_2;
				
				augmented_prev = augmentedLagrangian.apply(args_1, further_args);
				
			}
		
			if(unconstrained_optimizer == "Newton"){
				
				NewtonMethod unconstOptimizer = new NewtonMethod(1);
				
				unconstOptimizer.set_version("ModNewton");
				unconstOptimizer.do_Newton_Optimization(args_1, augmentedLagrangian, further_args);
				
			}
			
			if(unconstrained_optimizer == "ConjugateGradient"){
				
				ConjugateGradient unconstOptimizer = new ConjugateGradient(1);
				
				unconstOptimizer.do_Conjugate_Gradient_Optimization(args_1, augmentedLagrangian, further_args);
				
				
			}
			
			if(unconstrained_optimizer == "L_BFGS"){
				
				L_BFGS  unconstOptimizer = new L_BFGS(1);
				
			    unconstOptimizer.do_LBFGS_Optimization(args_1, augmentedLagrangian, further_args);
				
			}
			
			args_2 = UnconstrainedOptimizer.optimal_candidate;
			
			if(number_of_equality_constraints != 0 && number_of_inequality_constraints != 0){
				
				equality_values   = evaluate_constraints(equality_constraints, args_2, further_args_eq_constraints);				
				inequality_values = evaluate_constraints(inequality_constraints, args_2, further_args_ineq_constraints);
				
				update_equality_multiplier();
				update_inequality_multiplier();
								
			}
			
			if(number_of_equality_constraints != 0 && number_of_inequality_constraints == 0){
				
				equality_values   = evaluate_constraints(equality_constraints, args_2, further_args_eq_constraints);				
								
				update_equality_multiplier();				
				
			}
						
			if(number_of_equality_constraints == 0 && number_of_inequality_constraints != 0){
				
				inequality_values = evaluate_constraints(inequality_constraints, args_2, further_args_ineq_constraints);
				
				update_inequality_multiplier();
				
			}
			
			update_penalty_parameter();
				
			if(Math.abs(augmented_prev - augmentedLagrangian.apply(args_2, further_args)) < convergence_criterion){
				
				System.out.println("Augmented Lagrange Optimization converged after " + i + " iterations");
								
				break;
				
			}
			
			//System.out.println(targetFunction(args_2, further_args));
			
		}
		
		System.out.println(targetFunction(args_2, further_args));
		System.out.println("");
		
		//MatrixOperations.print_matrix(inequality_multiplier);
		
		//MatrixOperations.print_matrix(equality_multiplier);
		System.out.println("");
		//MatrixOperations.print_matrix(equality_values);
		
		System.out.println("");
		
		opt_candidate = args_2;
		
		return args_2;
	
	}
	
	
	// returns value of augmented Lagrangian with equality and inequality constraints
	public static double augmented_Lagrangian_with_equality_and_inequality_constraints(double [] opt_args, double [] further_args){
		
		double objectiv_value         = targetFunction(opt_args, further_args);

		equality_values   = evaluate_constraints(equality_constraints, opt_args, further_args_eq_constraints);
		
		inequality_values = evaluate_constraints(inequality_constraints, opt_args, further_args_ineq_constraints);
		
		double term_1 = objectiv_value;
		double term_2 = MatrixOperations.multiplication(MatrixOperations.transpose(equality_multiplier), equality_values)[0][0];
		double term_3 = penalty/2.0* Math.pow(MatrixOperations.euclidian(equality_values),2.0);
		
		double term_4 = 0.0;
		
		for(int i = 0; i < number_of_equality_constraints; i++){
			
			term_4 += Math.pow(Math.max(0.0, inequality_multiplier[i][0] + penalty*inequality_values[i][0]),2.0) - Math.pow(inequality_multiplier[i][0], 2.0);
			
		}
		
		term_4 = term_4/(2.0*penalty);
			
		double augmented_lagrangian = term_1 + term_2 + term_3 + term_4; 
				
		return augmented_lagrangian;	
		
	}
	
	
	// returns value of augmented Lagrangian with equality constraints
	public static double augmented_Lagrangian_with_equality_constraints(double [] opt_args, double [] further_args){
		
		double objectiv_value         = targetFunction(opt_args, further_args);

		equality_values   = evaluate_constraints(equality_constraints, opt_args, further_args_eq_constraints);
		
		double term_1 = objectiv_value;
		double term_2 = MatrixOperations.multiplication(MatrixOperations.transpose(equality_multiplier), equality_values)[0][0];
		double term_3 = penalty/2.0* Math.pow(MatrixOperations.euclidian(equality_values),2.0);
		
		double augmented_lagrangian = term_1 + term_2 + term_3; 
		
		return augmented_lagrangian;	
		
	}
	
	
	// returns value of augmented Lagrangian with inequality constraints
	public static double augmented_Lagrangian_with_inequality_constraints(double [] opt_args, double [] further_args){
		
		double objectiv_value         = targetFunction(opt_args, further_args);

		inequality_values   = evaluate_constraints(inequality_constraints, opt_args, further_args_ineq_constraints);
		
		double term_1 = objectiv_value;
		double term_2 = MatrixOperations.multiplication(MatrixOperations.transpose(inequality_multiplier), inequality_values)[0][0];
		double term_3 = penalty/2.0* Math.pow(MatrixOperations.euclidian(inequality_values),2.0);
		
		double augmented_lagrangian = term_1 + term_2 + term_3; 
		
		return augmented_lagrangian;
		
	}
	
	
	// sets updated Lagrangian multiplier for equality constraints
	public static void update_equality_multiplier(){
		
		equality_multiplier = MatrixOperations.add(equality_multiplier, MatrixOperations.scalar_multiplication(penalty, equality_values));
		
	}
	
	
	// sets updated Lagrangian multiplier for inequality constraints
	public static void update_inequality_multiplier(){
			
		for(int i = 0; i < number_of_inequality_constraints; i++){
			
			inequality_multiplier[i][0] = Math.max(0.0,  inequality_multiplier[i][0] + penalty*inequality_values[i][0]);
			
		}
		
	}
	
	
	// sets updated penalty parameter
	public static void update_penalty_parameter(){
		
		penalty =  penalty_growth_factor*penalty;
		
	}
	
		
	// returns possible unconstrained optimizer for augmented Lagrangian method (ALM)
	public static String [] get_unconstrained_optimizer_4_ALM(){
		
		String [] unconst_optimizer = {"L_BFGS", "Newton", "ConjugateGradient"};
		
		return unconst_optimizer;
		
	}
	
	
	// Unit test no. 1: Nonlinear objective function with 2 nonlinear inequality constraints
	//    No further arguments for objective functions and constraints
	//    Optimal solution: [0.75, 4.56]
	public static void test_1(){
		
		int n_iterations = 100000;
		
		double [] start_values = {0.0, 1.0};
		double [] further_args = null;
		
		ArrayList <BiFunction<double [], double [], Double>> eq_constraints = null;
		GenList further_eq_args = null;
		
		ArrayList <BiFunction<double [], double [], Double>> ineq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(2);
        ineq_constraints.add(TargetFunction::inequality_constraint_1);
        ineq_constraints.add(TargetFunction::inequality_constraint_2);
        
        GenList further_ineq_args = null;
 
		double [] solution = do_Augmented_Lagrange_Optimization(start_values, TargetFunction::constrained_target_function_1, further_args, eq_constraints, ineq_constraints, further_eq_args, further_ineq_args, n_iterations);
		
		MatrixOperations.print_vector(solution);
		
	}
	
	
	// Unit test no. 2: Nonlinear objective function with 2 nonlinear inequality constraints
	//    Further arguments for both constraints and objective function
    //    Optimal solution: [0.75, 4.56]
	public static void test_2(){
		
		int n_iterations = 100000;
		
		double [] start_values = {0.0, 1.0};
		double [] further_args = {1.0, 2.0, 5.0, 2.0};
		
		ArrayList <BiFunction<double [], double [], Double>> eq_constraints = null;
		GenList further_eq_args = null;
		
		ArrayList <BiFunction<double [], double [], Double>> ineq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(2);
        ineq_constraints.add(TargetFunction::inequality_constraint_with_further_args_1);
        ineq_constraints.add(TargetFunction::inequality_constraint_with_further_args_2);
        
        GenList further_ineq_args = new GenList(2);
 
        double [] arg1 = {2.0, 4.0};
		double [] arg2 = {-1.0, 2.0, 2.0, 3.0};
		further_ineq_args.add(arg1);
		further_ineq_args.add(arg2);
        
		double [] solution = do_Augmented_Lagrange_Optimization(start_values, TargetFunction::constrained_target_function_with_further_args_1, further_args, eq_constraints, ineq_constraints, further_eq_args, further_ineq_args, n_iterations);
		
		MatrixOperations.print_vector(solution);
		
	}
	
	
	// Unit test no. 3: Nonlinear objective function with 1 nonlinear equality and 1 nonlinear inequality constraints
	//    Further arguments for inequality constraints and objective function
	//    Optimal solution: [0.75, 4.56]
	public static void test_3(){
		
		int n_iterations = 100000;
		
		double [] start_values = {2.0, -10.0};
		double [] further_args = {1.0, 2.0, 5.0, 2.0};
		
		ArrayList <BiFunction<double [], double [], Double>> eq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(1);
		eq_constraints.add(TargetFunction::equality_constraint_1);
		
		GenList further_eq_args = null;
		
		ArrayList <BiFunction<double [], double [], Double>> ineq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(1);
        ineq_constraints.add(TargetFunction::inequality_constraint_with_further_args_1);
        
        GenList further_ineq_args = new GenList(1);
 
        double [] arg1 = {2.0, 4.0};
		further_ineq_args.add(arg1);
      
		double [] solution = do_Augmented_Lagrange_Optimization(start_values, TargetFunction::constrained_target_function_with_further_args_1, further_args, eq_constraints, ineq_constraints, further_eq_args, further_ineq_args, n_iterations);
		
		MatrixOperations.print_vector(solution);
		
	}
	
	
	// Unit test no. 4: Nonlinear objective function with 1 nonlinear equality and 1 nonlinear inequality constraints
	//    Further arguments for equality and inequality constraints and objective function
	//    Optimal solution: [0.75, 4.56]
	public static void test_4(){
		
		int n_iterations = 100000;
		
		double [] start_values = {2.0, -10.0};
		double [] further_args = {1.0, 2.0, 5.0, 2.0};
		
		ArrayList <BiFunction<double [], double [], Double>> eq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(1);
		eq_constraints.add(TargetFunction::equality_constraint_with_further_args_1);
		
		GenList further_eq_args = new GenList(1);
		
		double [] arg1 = {2.0, 2.0, 3.0};
		further_eq_args.add(arg1);
		
		ArrayList <BiFunction<double [], double [], Double>> ineq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(1);
        ineq_constraints.add(TargetFunction::inequality_constraint_with_further_args_1);
        
        GenList further_ineq_args = new GenList(1);
 
        double [] arg2 = {2.0, 4.0};
		further_ineq_args.add(arg2);
      
		double [] solution = do_Augmented_Lagrange_Optimization(start_values, TargetFunction::constrained_target_function_with_further_args_1, further_args, eq_constraints, ineq_constraints, further_eq_args, further_ineq_args, n_iterations);
		
		MatrixOperations.print_vector(solution);
		
	}

	
	// Unit test no. 5: Nonlinear objective function with 2 nonlinear equality constraints
	public static void test_5(){
		
		int n_iterations = 100000;
		
		double [] start_values = {-2.0, -105.0};
		double [] further_args = {1.0, 2.0, 5.0, 2.0};
		
		ArrayList <BiFunction<double [], double [], Double>> eq_constraints = new ArrayList <BiFunction<double [], double [], Double>>(1);
		eq_constraints.add(TargetFunction::equality_constraint_2);
		eq_constraints.add(TargetFunction::equality_constraint_with_further_args_3);
		
		GenList further_eq_args = new GenList(2);
		
		double [] arg1 = null;
		double [] arg2 = {2.0, 14};
		further_eq_args.add(arg1);
		further_eq_args.add(arg2);
		
		ArrayList <BiFunction<double [], double [], Double>> ineq_constraints = null;
        GenList further_ineq_args = null;
 
		double [] solution = do_Augmented_Lagrange_Optimization(start_values, TargetFunction::constrained_target_function_with_further_args_1, further_args, eq_constraints, ineq_constraints, further_eq_args, further_ineq_args, n_iterations);
		
		MatrixOperations.print_vector(solution);
		
	}
	

	// test client
	public static void main(String[] args){
				
		//test_1();
			
		//unconstrained_optimizer = "ConjugateGradient";
		//test_2();
		
		//unconstrained_optimizer = "ConjugateGradient";
		
	    //test_3();
		
		//test_4();
		
	    //test_5();
	    
	}
	
}
