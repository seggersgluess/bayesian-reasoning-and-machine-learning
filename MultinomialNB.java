package NaiveBayesClassifier;

import java.util.HashMap;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class MultinomialNB extends NBC{

	double [][] pi;
	double [][] my;
			
	
	public void fit() {
				
		pi = new double [n_classes][1];
		my = new double [n_classes][n_explaining_variables];
		
		double [][] n_c = new double [n_classes][1];
 		
		for(int c=0; c<n_classes; c++) {	
			double [] my_jc = new double [n_explaining_variables];
			double summedFeatures = 0.0;
						
			for(int i=0; i<n_observations; i++) {
				if(mod_explained_variable[i][0] == c) {
					for(int j=0; j<n_explaining_variables; j++) {
						my_jc[j] += explaining_variables[i][j];											
					}
					n_c[c][0]++;
				}
				
			}
			
			for(int j=0; j<n_explaining_variables; j++) {
				summedFeatures += my_jc[j];
			}
					
			for(int j=0; j<n_explaining_variables; j++) {										
				my_jc[j] += alpha;
				my_jc[j] /= (summedFeatures+alpha*n_explaining_variables);	
				my[c][j] = my_jc[j];
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
			double summedFeat = GeneralMath.sum(MatrixOperations.get_row_vec_from_matrix(X, i));
			double n_fact = GeneralMath.factorial(summedFeat);
			for(int c=0; c<n_classes; c++) {
				l_ic[c] = Math.log(pi[c][0]);
				for(int j=0; j<n_explaining_variables; j++) {	
					double prob = 1.0/GeneralMath.factorial(X[i][j]);
					prob *= Math.pow(my[c][j],X[i][j]);					
					l_ic[c] += Math.log(prob);
				}					
				l_ic[c] += Math.log(n_fact);						
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
	
}
