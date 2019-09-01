package StateSpaceModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.MatrixOperations;

public class KalmanFilter extends StateSpaceModel{
	
	static ArrayList<List<Double>> pred_state_vars;
	static ArrayList<List<Double>> pred_covariance;
	
	static ArrayList<List<Double>> update_state_vars;
	static ArrayList<List<Double>> update_covariance;
	
	static ArrayList<List<Double>> kalman_gain;
	static ArrayList<List<Double>> residuals;
	static ArrayList<List<Double>> res_covariance;
	
	static List<Double> prior_mean;
	static List<Double> prior_covariance;
	
	static double log_likelihood;
	
	public KalmanFilter(double [][] obs_variables, double [][] m_const, double [][] m_matrix, double [][] m_cov, double [][] t_const, double [][] t_matrix, double [][] t_cov) {
		super(obs_variables, m_const, m_matrix, m_cov, t_const, t_matrix, t_cov);

		pred_state_vars = new ArrayList<List<Double>>(n_measurements);
		pred_covariance = new ArrayList<List<Double>>(n_measurements);
		
		update_state_vars = new ArrayList<List<Double>>(n_measurements);
		update_covariance = new ArrayList<List<Double>>(n_measurements);
		
		kalman_gain = new ArrayList<List<Double>>(n_measurements);
		residuals = new ArrayList<List<Double>>(n_measurements);
		res_covariance = new ArrayList<List<Double>>(n_measurements);
		
		first_initialization_of_kalman_filter();
		
	}

	
	public static void first_initialization_of_kalman_filter(){
			
		prior_mean = new ArrayList<Double>();
		prior_covariance = new ArrayList<Double>();
		
		if(useMeasurementConstant == true){
			prior_mean.add(1.0);
		}
			
		for(int i=0; i<n_state_variables; i++){
			prior_mean.add(0.0);
		}
		
		if(useTransitionConstant == true){
			for(int i=0; i<n_state_variables; i++){
				prior_mean.add(transition_constant[i][0]);
			}
		}

		int n=n_state_variables;
		
		if(useMeasurementConstant == true){
			n++;
		}
		
		if(useTransitionConstant == true){
			n += n_state_variables;
		}
		
		double [][] vec_prior_covariance = new double [n][n];
		
		int idx=0; 
		if(useMeasurementConstant == true){
			idx++;
		}
		
		for(int i=0; i<n_state_variables; i++){
			vec_prior_covariance[idx+i][idx+i] = 1e+02;			
		}
		
		prior_covariance = MatrixOperations.vecAsList(vec_prior_covariance);
		
	}
		
	
	public static void run_kalman_filter(){
		
		int n = prior_mean.size();
		
		List<Double> pred_s = new ArrayList<Double>(n);
		
		for(int i=0; i<n; i++){
			double value = 0.0;
			for(int j=0; j<n; j++){
				value += transition_matrix[i][j]*prior_mean.get(j);
			}
			pred_s.add(value);
		}
		
		pred_state_vars.add(pred_s);
		
		double [][] pred_cov_s = MatrixOperations.get_matrix_from_vec(prior_covariance, n, n);
		pred_cov_s = MatrixOperations.multiplication(MatrixOperations.multiplication(transition_matrix, pred_cov_s),MatrixOperations.transpose(transition_matrix));
		pred_cov_s = MatrixOperations.add(pred_cov_s, transition_covariance);
		
		pred_covariance.add(MatrixOperations.vecAsList(pred_cov_s));
		
        log_likelihood = 0.0;
		
		for(int t=0; t<n_measurements; t++){
			
			double [][] pred_cov = MatrixOperations.get_matrix_from_vec(pred_covariance.get(t), n, n);
			
			List<Double> res = new ArrayList<Double>(n_measurement_variables);
			double [][] res_vec = new double [n_measurement_variables][1];
			
			for(int i=0; i<n_measurement_variables; i++){
				for(int j=0; j<n; j++){
					res_vec[i][0] += -measurement_matrix[i][j]*pred_state_vars.get(t).get(j);
				}
				res_vec[i][0] += measurement_variables[t][i];
				res.add(res_vec[i][0]);
			}
			
			residuals.add(res);
			
			double [][] matrix = new double [n][n_measurement_variables];
			
			//Calculate S_k!
			//matrix = P^-_k x H^T_k -> is also used in Kalman gain K_k! (more efficient).
			for(int i=0; i<n_measurement_variables; i++){
				for(int j=0; j<n; j++){
					for(int k=0; k<n; k++){
						matrix[j][i] += pred_cov[j][k]*measurement_matrix[i][k];
					}					
				}
			}

		    double [][] S_k = new double [n_measurement_variables][n_measurement_variables];
			
		    for(int i=0; i<n_measurement_variables; i++){
		    	for(int j=0; j<n_measurement_variables; j++){
		    		for(int k=0; k<n; k++){
		    			S_k[j][i] += measurement_matrix[j][k]*matrix[k][i];
		    		}
		    		S_k[j][i] += measurement_covariance[j][i];
		    	}
		    }
		    
		    res_covariance.add(MatrixOperations.vecAsList(S_k));
		    
			//Calculate K_k
			double [][] S_k_inv = MatrixOperations.inverse(S_k);
		    double [][] K_k = MatrixOperations.multiplication(matrix, S_k_inv);
			
		    kalman_gain.add(MatrixOperations.vecAsList(K_k));
			
			//Calculate m_k
		    List<Double> update_s = new ArrayList<Double>(n);
		    matrix = new double [n][n_measurement_variables];
		    
			for(int i=0; i<n; i++){
				double update = 0.0;
				for(int j=0; j<n_measurement_variables; j++){					
					update += K_k[i][j]*res.get(j);
				}
				update += pred_state_vars.get(t).get(i);
				update_s.add(update);
				
				for(int j=0; j<n_measurement_variables; j++){					
					for(int k=0; k<n_measurement_variables; k++){					
						matrix[i][j] += -K_k[i][k]*S_k[k][j];
					}
				}
					
			}
			
			update_state_vars.add(update_s);
			
			//Calculate P_k & predicted s
			double [][] P_k = new double [n][n];
			pred_s = new ArrayList<Double>(n);
		    for(int i=0; i<n; i++){
		    	double value = 0.0;
		    	for(int j=0; j<n; j++){
		    		value += transition_matrix[i][j]*update_s.get(j);
		    		for(int k=0; k<n_measurement_variables; k++){
		    			P_k[j][i] += matrix[j][k]*K_k[i][k];
		    		}
		    		P_k[j][i] += pred_cov[j][i];
		    	}
		    	pred_s.add(value);
		    }
			
		    update_covariance.add(MatrixOperations.vecAsList(P_k));
		    pred_state_vars.add(pred_s);
		    
		    //Calculate pred covariance P_pred_k
		    matrix = MatrixOperations.multiplication(transition_matrix, P_k);
		    pred_cov_s = new double [n][n];
		    
		    for(int i=0; i<n; i++){
		    	for(int j=0; j<n; j++){
		    		for(int k=0; k<n; k++){
		    			pred_cov_s[j][i] += matrix[j][k]*transition_matrix[i][k];
		    		}
		    		pred_cov_s[j][i] += transition_covariance[j][i];
		    	}
		    }
		    
		    pred_covariance.add(MatrixOperations.vecAsList(pred_cov_s));
		    
			double S_k_det = MatrixOperations.determinant(S_k);
			
			log_likelihood += MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(res_vec), S_k_inv),res_vec)[0][0];
			log_likelihood += Math.log(S_k_det);
		    
		}
		
		log_likelihood *= -0.5;
		log_likelihood += -n_measurements/2.0*n_measurement_variables*Math.log(2.0*Math.PI);
		
	}
	
	
	public static void run_kalman_filter_ProtoType(){
		
		int n = prior_mean.size();
		
		List<Double> pred_s = new ArrayList<Double>(n);
		
		//Prediction with t=0
		for(int i=0; i<n; i++){
			double value = 0.0;
			for(int j=0; j<n; j++){
				value += transition_matrix[i][j]*prior_mean.get(j);
			}
			pred_s.add(value);
		}
		
		pred_state_vars.add(pred_s);
		
		double [][] pred_cov_s = MatrixOperations.get_matrix_from_vec(prior_covariance, n, n);
		pred_cov_s = MatrixOperations.multiplication(MatrixOperations.multiplication(transition_matrix, pred_cov_s),MatrixOperations.transpose(transition_matrix));
		pred_cov_s = MatrixOperations.add(pred_cov_s, transition_covariance);
		
		pred_covariance.add(MatrixOperations.vecAsList(pred_cov_s));
		
        log_likelihood = 0.0;
		
		for(int t=0; t<n_measurements; t++){
			
			double [][] m_k_pred = MatrixOperations.get_matrix_from_vec(pred_state_vars.get(t), n, 1);
			double [][] P_k_pred = MatrixOperations.get_matrix_from_vec(pred_covariance.get(t), n, n);
			double [][] y_k = MatrixOperations.get_sub_matrix_between_row_idxs(measurement_variables, t, t);
			
			//Calculate residuals v_k
			double [][] v_k = MatrixOperations.multiplication(measurement_matrix, m_k_pred);
			v_k = MatrixOperations.substract(y_k, v_k);
						
			residuals.add(MatrixOperations.vecAsList(v_k));
			
			double [][] matrix = new double [n][n_measurements];
			
			//Calculate covariance of residualsS_k!
			matrix = MatrixOperations.multiplication(P_k_pred, MatrixOperations.transpose(measurement_matrix));
			double [][] S_k = MatrixOperations.add(MatrixOperations.multiplication(measurement_matrix, matrix), measurement_covariance);
			res_covariance.add(MatrixOperations.vecAsList(S_k));
			
			//Calculate Kalman gain K_k
			double [][] S_k_inv = MatrixOperations.inverse(S_k);
			double [][] K_k = MatrixOperations.multiplication(matrix, S_k_inv);			
			kalman_gain.add(MatrixOperations.vecAsList(K_k));
			
			//Calculate update state variables m_k
			double [][] m_k = MatrixOperations.multiplication(K_k, v_k);
			m_k = MatrixOperations.add(m_k_pred, m_k);
			update_state_vars.add(MatrixOperations.vecAsList(m_k));
			
			//Calculate update covariance P_k
			double [][] P_k = MatrixOperations.multiplication(MatrixOperations.multiplication(K_k, S_k),MatrixOperations.transpose(K_k));
			P_k = MatrixOperations.substract(P_k_pred, P_k);
			update_covariance.add(MatrixOperations.vecAsList(P_k));
			
			//Calculate pred state_variables m_pred_k
			m_k_pred = MatrixOperations.multiplication(transition_matrix,m_k);
			pred_state_vars.add(MatrixOperations.vecAsList(m_k_pred));
			
			//Calculate pred covariance P_pred_k
			P_k_pred = MatrixOperations.multiplication(MatrixOperations.multiplication(transition_matrix, P_k),MatrixOperations.transpose(transition_matrix));
			P_k_pred = MatrixOperations.add(P_k_pred, transition_covariance);
			pred_covariance.add(MatrixOperations.vecAsList(P_k_pred));
			
			double S_k_det = MatrixOperations.determinant(S_k);
			
			log_likelihood += MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(v_k), S_k_inv),v_k)[0][0];
			log_likelihood += Math.log(S_k_det);
			
		}
		
		log_likelihood *= -0.5;
		log_likelihood += -n_measurements/2.0*n_measurement_variables*Math.log(2.0*Math.PI);
		
	}
		
	
	public static void set_kf_prior_mean(double [][] prior_mean_vec){
		
		prior_mean = new ArrayList<Double>();
		if(useMeasurementConstant == true){
			prior_mean.add(1.0);
		}
			
		for(int i=0; i<n_state_variables; i++){
			prior_mean.add(prior_mean_vec[i][0]);
		}
		
		if(useTransitionConstant == true){
			for(int i=0; i<n_state_variables; i++){
				prior_mean.add(transition_constant[i][0]);
			}
		}
	}
	
	
	public static void set_kf_prior_covariance(double [][] prior_cov_matrix){
		
		prior_covariance = new ArrayList<Double>();
		
		int n=n_state_variables;
		
		if(useMeasurementConstant == true){
			n++;
		}
		
		if(useTransitionConstant == true){
			n += n_state_variables;
		}
		
		double [][] vec_prior_covariance = new double [n][n];
		
		int idx=0; 
		if(useMeasurementConstant == true){
			idx++;
		}
		
		for(int i=0; i<n_state_variables; i++){
			for(int j=0; j<n_state_variables; j++){
				vec_prior_covariance[idx+i][idx+j] = prior_cov_matrix[i][j];
			}
		}
		
		prior_covariance = MatrixOperations.vecAsList(vec_prior_covariance);
	}
		
	
	public static double get_kf_log_likelihood(){
		return log_likelihood;
	}
	
	
	public static double [][] get_kf_residual_vec(int step){
		
		double [][] res_vec = MatrixOperations.get_matrix_from_vec(residuals.get(step), n_measurement_variables, 1);
		
		return res_vec;
		
	}
	
	
	public static double [][] get_kf_update_covariance(int step){
		
		int n = (int)Math.sqrt(update_covariance.get(step).size());
		
		double [][] updated_cov = MatrixOperations.get_matrix_from_vec(update_covariance.get(step), n, n);
		
		return updated_cov;
		
	}
	
	
	public static double [][] get_kf_predicted_covariance(int step){
		
		int n = (int)Math.sqrt(update_covariance.get(step).size());
		
		double [][] pred_cov = MatrixOperations.get_matrix_from_vec(pred_covariance.get(step), n, n);
		
		return pred_cov;
		
	}
	
	
}
