package ComponentModels;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class ComponentModels {

    //--- parent class for component models ---	
	
	double [][] X;
	int n_variables;
	int n_observations;
	
	double [][] X_scaled;
	
	boolean center = true;
	boolean scale  = false;
	
	double [][] center_pars;
	double [][] scale_pars;
	
	//Estimated X (implied by reduced base z)
	double [][] rotated_X;

	
	public ComponentModels(double [][] X, boolean scale) {
		
		int n_vars = X[0].length;
		
		if(n_vars > X.length) {
			throw new RuntimeException("Not enough observations supplied. (Number of features > number of observations)!");
		}
		
		this.X = X;
		this.n_observations = X.length;
		this.n_variables = n_vars;
		this.scale = scale;
		
		scale_and_center_data();		
	}
	
	
	public ComponentModels(double [][] X, boolean center, boolean scale) {
		
		int n_vars = X[0].length;
		
		if(n_vars > X.length) {
			throw new RuntimeException("Not enough observations supplied. (Number of features > number of observations)!");
		}
				
		this.X = X;
		this.n_observations = X.length;
		this.n_variables = n_vars;
		this.center = center;
		this.scale = scale;
		
		scale_and_center_data();		
	}
		
	
	public void scale_and_center_data() {
		
		if(center == false && scale == false) {			
			X_scaled = X;			
		}
		
		if(center == true && scale == true) {
			center_pars = new double [1][n_variables];
			scale_pars  = new double [1][n_variables];
			X_scaled = new double [n_observations][n_variables];
			for(int i=0; i<n_variables; i++) {
				double [][] x_col = MatrixOperations.get_column_vec_from_matrix(X, i);
				center_pars[0][i] = GeneralMath.mean(x_col);
				scale_pars[0][i]  = GeneralMath.sd(x_col);
				for(int j=0; j<n_observations; j++) {
					X_scaled[j][i] = (X[j][i]-center_pars[0][i])/scale_pars[0][i];
				}		
			}			
		}
		
		if(center == true && scale == false) {
			center_pars = new double [1][n_variables];
			X_scaled = new double [n_observations][n_variables];
			for(int i=0; i<n_variables; i++) {
				double [][] x_col = MatrixOperations.get_column_vec_from_matrix(X, i);
				center_pars[0][i] = GeneralMath.mean(x_col);
				for(int j=0; j<n_observations; j++) {
					X_scaled[j][i] = X[j][i]-center_pars[0][i];
				}		
			}
		}
		
		if(center == false && scale == true) {
			scale_pars = new double [1][n_variables];
			X_scaled = new double [n_observations][n_variables];
			for(int i=0; i<n_variables; i++) {
				double [][] x_col = MatrixOperations.get_column_vec_from_matrix(X, i);
				scale_pars[0][i] = GeneralMath.sd(x_col);
				for(int j=0; j<n_observations; j++) {
					X_scaled[j][i] = X[j][i]/scale_pars[0][i];
				}		
			}			
		}	
	}
	
	
	public double [][] get_scaled_input_data() {
		if(center == false && scale == false) {
			System.out.println("Data neither centered nor scaled. No scaled input available.");
			return null;
		}else {
			return X_scaled;
		}
	}
	
	
	public double [][] get_input_data() {
		if(X == null) {
			System.out.println("No input data supplied yet.");
			return null;
		}else {
			return X;
		}
	}
	
	
	public double [][] get_center_pars() {
		
		if(center == false) {
			System.out.println("Data not centered.");
			return null;
		}else {
			return center_pars;
		}
	}
	
	
	public double [][] get_scale_pars() {
		
		if(scale == false) {
			System.out.println("Data not scaled.");
			return null;
		}else {
			return scale_pars;
		}
	}
	
	
	public double [][] get_rotated_input() {
		if(rotated_X == null) {
			System.out.println("No PCA done yet.");
			return null;
		}else {
			return rotated_X;
		}
	}
	
}
