package StateSpaceModels;

import Mathematics.MatrixOperations;

public class StateSpaceModel {

	static double [][] measurement_variables;
	static double [][] state_variables;
	
	static int n_measurement_variables;
	static int n_measurements;
	
	static int n_state_variables;
	
	static double [][] measurement_constant;
	static double [][] measurement_matrix;
	static double [][] measurement_covariance;
	
	static double [][] transition_constant;
	static double [][] transition_matrix;
	static double [][] transition_covariance;
	
	static boolean useMeasurementConstant;
	static boolean useTransitionConstant;
	
	/**
     * Constructor for State Space model class
     * @param obs_variables double vector containing data series of measurements
     * @param m_const double vector for setting measurement constants (set null if not used)
     * @param m_matrix double matrix for setting measurement matrix
     * @param m_cov double matrix for setting measurement covariance
     * @param t_const double vector for setting transtition constants (set null if not used)
     * @param t_matrix double matrix for setting transition matrix
     * @param t_cov double matrix for setting space variable´s covariance
     */
	public StateSpaceModel(double [][] obs_variables, double [][] m_const, double [][] m_matrix, double [][] m_cov, double [][] t_const, double [][] t_matrix, double [][] t_cov){
		
		measurement_variables = obs_variables;
		
		if(m_const != null){
			useMeasurementConstant = true;
		}
		
		if(t_const != null){
			useTransitionConstant = true;
		}
				
		n_measurement_variables = measurement_variables[0].length;		
		n_measurements = measurement_variables.length;
		
		n_state_variables = t_matrix.length;
		
		if(useMeasurementConstant == true){
			set_measurement_constant(m_const);
		}
		
		if(useTransitionConstant == true){
			set_transition_constant(t_const);
		}
		
		set_measurement_matrix(m_matrix);
		set_transition_matrix(t_matrix);
		
		set_measurement_covariance(m_cov);
		set_transition_covariance(t_cov);
	}
	
	
	public static void set_measurement_constant(double [][] m_constant){
		measurement_constant = m_constant;		
	}
	
	
	public static void set_measurement_matrix(double [][] m_matrix){				
		measurement_matrix = m_matrix;	
		
		if(useMeasurementConstant == true && useTransitionConstant == false){
			measurement_matrix = MatrixOperations.cbind(measurement_matrix, measurement_constant);			
		}
		
		if(useMeasurementConstant == true  && useTransitionConstant == true){			
			int nRows = measurement_matrix.length;
			int nCols = measurement_matrix[0].length+1;
			
			m_matrix = new double [nRows+n_state_variables][nCols];
			
			for(int i=0; i<nRows; i++){
				m_matrix[i][0] = measurement_constant[i][0];
			}
			
			for(int i=0; i<nRows; i++){
				for(int j=0; j<(nCols-1); j++){
					m_matrix[i][j+1] = measurement_matrix[i][j];
				}
			}
			
			measurement_matrix = m_matrix;
		}
		
	}
	
	
	public static void set_transition_constant(double [][] t_constant){
		transition_constant = t_constant;
	}
	
	
	public static void set_transition_matrix(double [][] t_matrix){
		transition_matrix = t_matrix;
		
		if(useTransitionConstant == false && useMeasurementConstant == true){
			transition_matrix = new double [n_state_variables+1][n_state_variables+1];
			transition_matrix[0][0] = 1.0;
			for(int i=0; i<n_state_variables; i++){
				for(int j=0; j<n_state_variables; j++){
					transition_matrix[i+1][j+1] = t_matrix[i][j];
				}
			}
		}
		
		if(useTransitionConstant == true && useMeasurementConstant == false){							
			transition_matrix = new double [2*n_state_variables][2*n_state_variables];
			int idx = n_state_variables-1;
			for(int i=0; i<n_state_variables; i++){
				for(int j=0; j<n_state_variables; j++){
					transition_matrix[i][j] = t_matrix[i][j];
				}
				transition_matrix[idx][idx] = 1.0;
				idx++;
			}							
		}
		
		if(useTransitionConstant == true && useMeasurementConstant == true){
			transition_matrix = new double [2*n_state_variables+1][2*n_state_variables+1];
			transition_matrix[0][0] = 1.0;
			int idx = n_state_variables;
			for(int i=0; i<n_state_variables; i++){
				for(int j=0; j<n_state_variables; j++){
					transition_matrix[i+1][j+1] = t_matrix[i][j];
				}
				transition_matrix[idx][idx] = 1.0;
				idx++;
			}	
		}
		
	}
	
	
	public static void set_measurement_covariance(double [][] m_cov_matrix){		
		measurement_covariance = m_cov_matrix;
	}
	
	
	public static void set_transition_covariance(double [][] t_cov_matrix){
		transition_covariance = t_cov_matrix;
		
		if(useTransitionConstant == false && useMeasurementConstant == true){
			double [][] t_cov = new double [n_state_variables+1][n_state_variables+1];
			
			for(int i=0; i<n_state_variables; i++){
				for(int j=0; j<n_state_variables; j++){
					t_cov[i+1][j+1] = transition_covariance[i][j];
				}
			}		
		}
		
		if(useTransitionConstant == true && useMeasurementConstant == false){
			double [][] t_cov = new double [2*n_state_variables][2*n_state_variables];
			
			for(int i=0; i<n_state_variables; i++){
				for(int j=0; j<n_state_variables; j++){
					t_cov[i][j] = transition_covariance[i][j];
				}
			}		
		}
		
		if(useTransitionConstant == true && useMeasurementConstant == true){
			double [][] t_cov = new double [2*n_state_variables+1][2*n_state_variables+1];
			
			for(int i=0; i<n_state_variables; i++){
				for(int j=0; j<n_state_variables; j++){
					t_cov[i+1][j+1] = transition_covariance[i][j];
				}
			}		
		}
		
	}
	
	
	public static double [][] get_measurement_matrix(){
			
		double [][] m_matrix = new double [n_measurement_variables][n_state_variables];
		
		if(useMeasurementConstant == true){
			m_matrix = MatrixOperations.get_sub_matrix_between_column_idxs(measurement_matrix, 1, n_state_variables);
		}else{
			m_matrix = measurement_matrix;
		}
		
		return m_matrix;
	}
	
	
	public static double [][] get_transition_matrix(){
		
		double [][] t_matrix = new double [n_state_variables][n_state_variables];
		
		if(useMeasurementConstant == true){
			int idx1 = 0;
			int idx2 = n_state_variables-1;
			t_matrix = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(transition_matrix, idx1, idx2, idx1, idx2);
		}else{
			t_matrix = transition_matrix;
		}
		
		return t_matrix;
	}
	
	
	public static double [][] get_measurement_covariance(){
		return measurement_covariance;
	}
	
	
	public static double [][] get_transition_covariance(){
		
		double [][] t_cov = new double [n_state_variables][n_state_variables];
		
		if(useMeasurementConstant == true){
			int idx1 = 0;
			int idx2 = n_state_variables-1;
			t_cov = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(transition_covariance, idx1, idx2, idx1, idx2);
		}else{
			t_cov = transition_covariance;
		}
		
		return t_cov;
	}
	
	
	public static double [][] get_state_variables(){
		
		double [][] state_vars = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(state_variables, 0, n_measurements-1, 0, n_state_variables-1);
		
		return state_vars;
		
	}
	
}
