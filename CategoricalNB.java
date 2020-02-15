package NaiveBayesClassifier;

import java.util.ArrayList;
import java.util.HashMap;

import Mathematics.MatrixOperations;

public class CategoricalNB extends NBC{

	double [][] pi;
	ArrayList<ArrayList<double []>> my;
	
	ArrayList<double[]> categories;
	int [] n_categories;
		
	
	public void fit() {
		
		handleCategories();
		
		pi = new double [n_classes][1];
		my = new ArrayList<ArrayList<double []>>(n_classes);
		
		double [][] n_c = new double [n_classes][1];
		
		for(int c=0; c<n_classes; c++) {
			ArrayList<double []> my_jc = new ArrayList<double []>(n_explaining_variables);		
			for(int j=0; j<n_explaining_variables; j++) {
				my_jc.add(new double [n_categories[j]]);
				for(int i=0; i<n_observations; i++) {
					if(mod_explained_variable[i][0] == c) {
						n_c[c][0]++;
						for(int k=0; k<n_categories[j]; k++) {
							if(explaining_variables[i][j] == categories.get(j)[k]) {
								my_jc.get(j)[k]++;
							}							
						}						
					}
				}
			}
			
			for(int j=0; j<n_explaining_variables; j++) {							
				for(int k=0; k<n_categories[j]; k++) {
					my_jc.get(j)[k] += alpha;
					my_jc.get(j)[k] /= (n_c[c][0]+alpha*n_categories[j]);
				}
			}
			my.add(my_jc);
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
					for(int k=0; k<n_categories[j]; k++) {
						if(X[i][j]==categories.get(j)[k]) {
							l_ic[c] += Math.log(my.get(c).get(j)[k]);
						}
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
	
	
	public void handleCategories() {
		
		double [] firstFeat = MatrixOperations.get_column_from_matrix(explaining_variables, 0);
		double [] uniqueFeat = Utilities.Utilities.get_unique_elements(firstFeat);
		
		n_categories = new int [n_explaining_variables];	
		categories = new ArrayList<double[]>(n_explaining_variables);
			
		for(int j=0; j<n_explaining_variables; j++) {
			
			firstFeat = MatrixOperations.get_column_from_matrix(explaining_variables, 0);
			uniqueFeat = Utilities.Utilities.get_unique_elements(firstFeat);
			
			n_categories[j] = uniqueFeat.length;
			categories.add(new double [n_categories[j]]);	
			
			for(int k=0; k<n_categories[j]; k++) {
				categories.get(j)[k] = uniqueFeat[k];
			}			
		}
	}
	
	
	public HashMap<String, Object> get_est_pars() {
		
		HashMap<String, Object> est_pars = new HashMap<String, Object>();
		
		est_pars.put("my", my);
		est_pars.put("pi", pi);
		
		return est_pars;
	}
	
	

	
}
