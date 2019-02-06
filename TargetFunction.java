
public class TargetFunction {

	static double [] gs_prev_arg = {};
	static double [] gs_direction  = {};
	
	
	public static double target_function(double [] x, double [] further_args){
		
		double a = 1.0;
		double b = 100.0;

		double target_function = Math.pow((a - x[0]),2)  + b*Math.pow(x[1]-Math.pow(x[0],2),2);
		
		return target_function;
		
	}
	
	
	public static double target_function_with_further_args(double [] x, double [] further_args){
		
		double a = further_args[0];
		double b = further_args[1];
		
		double target_function = Math.pow((a - x[0]),2)  + b*Math.pow(x[1]-Math.pow(x[0],2),2);
		
		return target_function;
		
	}
	
	
	public static double rastrigin_function(double [] x, double [] further_args){
		
		double target_function = 20.0 + Math.pow((x[0]),2)  - 10.0*Math.cos(2.0*Math.PI*x[0]) + Math.pow(x[1], 2) - 10.0*Math.cos(2.0*Math.PI*x[1]);
		
		return target_function;
		
	}
	
	
	// test target function constrained optimization problem 1 (with 2 inequality constraints)
	public static double constrained_target_function_1(double [] x, double [] further_args){
		
		double target_function = Math.pow((x[0]-1.0), 2.0) + Math.pow((x[1]-5.0), 2.0);
		
		return target_function;
		
	}
	
	
	// test equality constraint no. 1
	public static double equality_constraint_1(double [] x, double [] further_args){
		
		double constraint = -1.0*Math.pow((x[0]-2.0), 2.0) + x[1] - 3.0;
		
		return constraint;
		
	}
	
	
	// test equality constraint no. 1 with further arguments
	public static double equality_constraint_with_further_args_1(double [] x, double [] further_args){
		
		double par_1 = further_args[0]; //2.0
		double par_2 = further_args[1]; //2.0
		double par_3 = further_args[2]; //3.0
		
		double constraint = -1.0*Math.pow((x[0]-par_1), par_2) + x[1] - par_3;
		
		return constraint;
		
	}
	
	
	// test inequality constraint no. 2
	public static double equality_constraint_2(double [] x, double [] further_args){
		
		double constraint = x[1] - 10.0;
		
		return constraint;
		
	}
	
	
	// test inequality constraint no. 3
	public static double equality_constraint_with_further_args_3(double [] x, double [] further_args){
		
		double par_1 = further_args[0]; // 2.0
		double par_2 = further_args[1]; // 14.0
		
		double constraint = Math.pow(x[0], par_1) + x[1] - par_2;
		
		return constraint;
		
	}
	
	
	// test inequality constraint no. 1
	public static double inequality_constraint_1(double [] x, double [] further_args){
		
		double constraint = -1.0*Math.pow((x[0]), 2.0) + x[1] - 4.0;
		
		return constraint;
		
	}
	
	
	// test inequality constraint no. 2
	public static double inequality_constraint_2(double [] x, double [] further_args){
		
		double constraint = -1.0*Math.pow((x[0]-2.0), 2.0) + x[1] - 3.0;
		
		return constraint;
		
	}
	

	// test target function constrained optimization problem 1 with further arguments (with 2 inequality constraints)
	public static double constrained_target_function_with_further_args_1(double [] x, double [] further_args){
		
		double par_1 = further_args[0]; // 1.0
		double par_2 = further_args[1]; // 2.0
		double par_3 = further_args[2]; // 5.0
		double par_4 = further_args[3]; // 2.0
		
		double target_function = Math.pow((x[0]-par_1), par_2) + Math.pow((x[1]-par_3), par_4);
		
		return target_function;
		
	}
	
	
	// test inequality constraint no. 1 with further arguments
	public static double inequality_constraint_with_further_args_1(double [] x, double [] further_args){
		
		double par_1 = further_args[0]; //2.0
		double par_2 = further_args[1]; //4.0
		
		double constraint = -1.0*Math.pow((x[0]), par_1) + x[1] - par_2;
		
		return constraint;
		
	}
	
	
	// test inequality constraint no. 2 with further arguments
	public static double inequality_constraint_with_further_args_2(double [] x, double [] further_args){
		
		double par_1 = further_args[0]; //-1.0
		double par_2 = further_args[1]; // 2.0
		double par_3 = further_args[2]; // 2.0
		double par_4 = further_args[3]; // 3.0
		
		double constraint = par_1*Math.pow((x[0]-par_2), par_3) + x[1] - par_4;
		
		return constraint;
		
	}
	
	
	public static double target_function_4_golden_section(double alpha){
		
		int n_args = gs_prev_arg.length;
		
		double [] args = new double [n_args];
		double [] further_args = null;
		
		args    = MatrixOperations.add_vectors(gs_prev_arg, MatrixOperations.scalar_vector_multiplication(alpha, gs_direction));
		
		return target_function(args, further_args);
		 
	}
	
}
