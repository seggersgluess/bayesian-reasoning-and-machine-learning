package NaiveBayesClassifier;

import java.util.HashMap;

import Mathematics.MatrixOperations;

public class GaussianNB extends NBC{

	static double [][] my;
	static double [][] sigma;
	static double [][] pi;
	static double [][] n_c;
	
	
	public void fit() {
		
		my    = new double [n_classes][n_explaining_variables];
		sigma = new double [n_classes][n_explaining_variables];
		pi    = new double [n_classes][1];
		n_c   = new double [n_classes][1];
		
		for(int c=0; c<n_classes; c++) {
			for(int i=0; i<n_observations; i++) {
				if(mod_explained_variable[i][0] == c) {
					n_c[c][0]++;
					for(int j=0; j<n_explaining_variables; j++) {					
						my[c][j] += explaining_variables[i][j];
					}				
				}
			}
		}
		
		for(int c=0; c<n_classes; c++) {
			for(int j=0; j<n_explaining_variables; j++) {
				my[c][j]/=n_c[c][0];
			}
		}
			
		for(int c=0; c<n_classes; c++) {
			for(int i=0; i<n_observations; i++) {
				if(mod_explained_variable[i][0] == c) {
					for(int j=0; j<n_explaining_variables; j++) {
						sigma[c][j] += Math.pow(explaining_variables[i][j]-my[c][j],2.0);
						//sigma[c][j] /= n_c[c][0];
					}				
				}
			}
			
			for(int j=0; j<n_explaining_variables; j++) {
				sigma[c][j]/=(n_c[c][0]-1.0);
			}
			
			pi[c][0] = n_c[c][0]/n_observations;
		}
		
	}
	
	
	public HashMap<String, double [][]> predict(double [][] X) {
		
		if(X[0].length != n_explaining_variables) {
			throw new RuntimeException("Number of features for prediction not compatible with trained model.");
		}
		
		int n = X.length;
		
		double [] l_ic = new double [n_classes];
		double [][] p_ic = new double [n][n_classes];
		double [][] y_hat = new double [n][1];
		
		for(int i=0; i<n; i++) {
			for(int c=0; c<n_classes; c++) {
				l_ic[c] = Math.log(pi[c][0]);
				for(int j=0; j<n_explaining_variables; j++) {
					double mean = my[c][j];
				    double variance = sigma[c][j];
					double gaussian_prob = get_gaussianPDF(X[i][j], mean, variance);
					l_ic[c] += Math.log(gaussian_prob);
				}
			}
			
			double logSum = logSumExp(l_ic);
			for(int c=0; c<n_classes; c++) {
				p_ic[i][c] = Math.exp(l_ic[c]-logSum);
			}			
			
			double [] p_i = MatrixOperations.get_row_from_matrix(p_ic, i);
			double maxProb = Utilities.Utilities.getMax(p_i);
			int [] maxClass = Utilities.Utilities.get_idx(p_i, maxProb);			
			y_hat[i][0] = classes[maxClass[0]];
		}
		
		HashMap<String, double [][]> predRes = new HashMap<String, double [][]>(2);
		
		predRes.put("Prediction", y_hat);
		predRes.put("Probability",p_ic);
		
		return predRes;		
	}
	
	
	public static double get_gaussianPDF(double x, double mean, double variance){	
		return 1.0/(Math.sqrt(2.0*Math.PI*variance))*Math.exp(-Math.pow(x-mean, 2.0)/(2.0*variance));			
	}
	

	public String [] get_par_names() {
		String [] par_names = {"my",
							   "sigma",
				               "pi"
		};
		
		return par_names;
	}
	
	
	public HashMap<String, Object> get_est_pars() {
		
		HashMap<String, Object> est_pars = new HashMap<String, Object>();
				
		est_pars.put("my", my);
		est_pars.put("sigma", sigma);
		est_pars.put("pi", pi);
		
		return est_pars;
	}
	
}
