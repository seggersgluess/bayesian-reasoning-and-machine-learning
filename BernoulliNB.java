package NaiveBayesClassifier;

import java.util.HashMap;

import Mathematics.MatrixOperations;

public class BernoulliNB extends NBC{

	double [][] pi;
	double [][] my;
	
	
	public void fit() {
		
		pi = new double [n_classes][1];
		my = new double [n_classes][n_explaining_variables];
		
		double [][] n_c = new double [n_classes][1];
		double [][] n_jc = new double [n_classes][n_explaining_variables];
		
		for(int i=0; i<n_observations; i++) {
			for(int c=0; c<n_classes; c++) {
				if(mod_explained_variable[i][0] == c) {
					n_c[c][0]++;		
					for(int j=0; j<n_explaining_variables; j++) {
						if(explaining_variables[i][j] == 1.0) {
							n_jc[c][j]++;
						}
					}
				}
				
			}
		}
		
		for(int c=0; c<n_classes; c++) {
			pi[c][0] = n_c[c][0]/n_observations;
			for(int j=0; j<n_explaining_variables; j++) {
				my[c][j] = n_jc[c][j]/n_c[c][0];
				if(my[c][j] == 0.0) {
					my[c][j] = 1.0/Double.MAX_VALUE;
				}
			}
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
					if(X[i][j] == 1.0) {
						l_ic[c] += Math.log(my[c][j]);
					}else {
						l_ic[c] += Math.log(1.0-my[c][j]);
					}
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
		
	
	public HashMap<String, Object> get_est_pars() {
		
		HashMap<String, Object> est_pars = new HashMap<String, Object>();
				
		est_pars.put("my", my);
		est_pars.put("pi", pi);
		
		return est_pars;
	}
	
	
	public String [] get_par_names() {
		String [] par_names = {"my",
				               "pi"
		};
		
		return par_names;
	}
	
	
}
