import java.util.function.BiFunction;

public class NumDeriv {

	static double epsilon = 1e-04;
	
	
	// computes target function
	public static double targetFunction(BiFunction<double [], double [], Double> f, double [] x, double [] further_args){
		
		return f.apply(x, further_args);
		
	}
	
	
	// calculates the gradient of the target function
	public static double [] gradient(BiFunction<double [], double [], Double> f, double [] x, double [] further_args){
		
		int n_args = x.length;
		
		double [] gradient = new double [n_args];
		
		double [] x_1 = new double [n_args];
		double [] x_2 = new double [n_args];
		
		for(int i=0; i<n_args; i++){
					
			for(int j=0; j<n_args; j++){
				
				x_1[j] = x[j];
				x_2[j] = x[j];
				
			}	
				
			x_1[i] = x[i] + epsilon;
			x_2[i] = x[i] - epsilon;
			
			double value_1 = targetFunction(f, x_1, further_args);
			double value_2 = targetFunction(f, x_2, further_args);
			
			gradient[i] = (value_1 - value_2)/(2.0*epsilon);
			
		}
		
		return gradient;
		
	}
	
	
	// calculates the hessian of the target function
	public static double [][] hessian(BiFunction<double [], double [], Double> f, double [] x, double [] further_args){
			
		int n_args = x.length;
		
		double [][] hessian = MatrixOperations.matrix(n_args, n_args);
		
		for(int i = 0; i < n_args; i++){
			
			double [] x_1 = new double [n_args];
			double [] x_2 = new double [n_args];
			double [] x_3 = new double [n_args];
			double [] x_4 = new double [n_args];
			
			for(int j = 0; j < n_args; j++){
				
				x_1[j] = x[j];
				x_2[j] = x[j];
				x_3[j] = x[j];
				x_4[j] = x[j];
				
			}
			
			x_1[i] = x[i] + 2.0*epsilon;
			x_2[i] = x[i] + epsilon;
			x_3[i] = x[i] - epsilon;
			x_4[i] = x[i] - 2.0*epsilon;
			
			double value_1 = -1.0*targetFunction(f, x_1, further_args);
			double value_2 = 16.0*targetFunction(f, x_2, further_args);
			double value_3 = -30.0*targetFunction(f, x, further_args);
			double value_4 = 16.0*targetFunction(f, x_3, further_args);
			double value_5 = -1.0*targetFunction(f, x_4, further_args);
			
			hessian[i][i] = (value_1 + value_2 + value_3 + value_4 + value_5)/(12.0*Math.pow(epsilon, 2.0));
			
			//Cross derivatives
			for(int j = 0; j < n_args; j++){
				
				if(j != i){
					
					double [] x_cross_1 = new double [n_args];
					double [] x_cross_2 = new double [n_args];
					double [] x_cross_3 = new double [n_args];
					double [] x_cross_4 = new double [n_args];
					
					for(int k = 0; k < n_args; k++){
						
						x_cross_1[k] = x[k];
						x_cross_2[k] = x[k];
						x_cross_3[k] = x[k];
						x_cross_4[k] = x[k];
																	
					}
					
					x_cross_1[i] = x[i] + epsilon;
					x_cross_1[j] = x[j] + epsilon;
					
					x_cross_2[i] = x[i] + epsilon;
					x_cross_2[j] = x[j] - epsilon;
					
					x_cross_3[i] = x[i] - epsilon;
					x_cross_3[j] = x[j] + epsilon;
					
					x_cross_4[i] = x[i] - epsilon;
					x_cross_4[j] = x[j] - epsilon;
					
					double cross_value_1 = targetFunction(f, x_cross_1, further_args);
					double cross_value_2 = -1.0*targetFunction(f, x_cross_2, further_args);
					double cross_value_3 = -1.0*targetFunction(f, x_cross_3, further_args);
					double cross_value_4 = targetFunction(f, x_cross_4, further_args);
					
					hessian[i][j] = (cross_value_1 + cross_value_2 + cross_value_3 + cross_value_4)/(4.0*Math.pow(epsilon, 2.0));
					
				}
				
			}
					
		}
				
		return hessian;
		
	}
	
		
    public static void main(String[] args) {
    	
    	double [] a = new double [2];
    	a[0] = 2.0;
    	a[1] = 3.0;
    	
    	double [] further_args = null;
    	
    	//System.out.print(TargetFunction.target_function(a));
    	
    	System.out.print(gradient(TargetFunction::target_function, a, further_args)[0]);
    	//System.out.print("    ");
    	//System.out.print(gradient(a)[1]);
    	//System.out.print("    ");
    	//System.out.print(gradient(a)[2]);
    	
    	//MatrixOperations.print_matrix(hessian(a));
    	
    }
	
}
