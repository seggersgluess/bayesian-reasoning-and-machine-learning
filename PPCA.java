package ComponentModels;

import java.util.ArrayList;
import java.util.HashMap;

import Distributions.NormalDistribution;
import Mathematics.MatrixOperations;

public class PPCA extends FactorAnalysis{
	
	public PPCA(double[][] X, int n_factors, boolean scale) {
		super(X, n_factors, scale);
	}

	
	public PPCA(double [][] X, int n_factors, boolean center, boolean scale) {		
		super(X, n_factors, center, scale);		
	}
	
	
	public PPCA(double [][] X, int n_factors, int n_analyzers, boolean center, boolean scale) {		
		super(X, n_factors, n_analyzers, center, scale);		
	}
	
	
	public void do_PPCA() {		
		em_mixtureFactorAnalyzers4PPCA();
		calc_rotated_X();
	}
	
	
	public void em_mixtureFactorAnalyzers4PPCA() {
		
		HashMap<String, ArrayList<double [][]>> pars = new HashMap<String, ArrayList<double [][]>>();
		if(ext_init_values == null)	{
			pars = get_input_pars_4_PPCA_em();
		}else {
			pars = ext_init_values;
		}
		
		double [][] x_i  = new double [n_variables][1];
		
		double [][] my_k = new double [n_variables][1];
		double [][] W_k  = new double [n_variables][n_factors];
		
		pi_k = new double [n_analyzers][1];
		r_ik = new double [n_observations][n_analyzers];
			
		double [][] diag_L = MatrixOperations.identity(n_factors);	
		
		for(int j=0; j<iterations; j++) {	
			
			factors = new ArrayList<double [][]>(n_analyzers);
			factor_covariances = new ArrayList<ArrayList<double [][]>>(n_analyzers);

			ArrayList<double [][]> myList = new ArrayList<double [][]>(n_analyzers);
			ArrayList<double [][]> W      = new ArrayList<double [][]>(n_analyzers);
			ArrayList<double [][]> piList  = new ArrayList<double [][]>(n_analyzers);
			ArrayList<double [][]> PsiList = new ArrayList<double [][]>(1);
				
			ArrayList<ArrayList<double [][]>> bList  = new ArrayList<ArrayList<double [][]>>();
			ArrayList<ArrayList<double [][]>> C_ik_List  = new ArrayList<ArrayList<double [][]>>();			

			ArrayList<double [][]> W_tilde = new ArrayList<double [][]>(n_analyzers);
			
			Psi = pars.get("Psi").get(0);
			double [][] Psi_inv = new double [n_variables][n_variables];
			for(int i=0; i<n_variables; i++) {
				Psi_inv[i][i] = 1.0/Psi[i][i];
			}
						
			double newLogLikelihood = Math.log(2.0*Math.PI)*n_observations*n_variables/2.0;
			newLogLikelihood -= n_observations*MatrixOperations.determinant(Psi)/2.0;
			
			double [] r_i_sum = new double [n_observations];
			
			for(int k=0; k<n_analyzers; k++) {
						
				my_k = pars.get("my").get(k);
				W_k  = pars.get("W").get(k);			
				pi_k[k][0] = pars.get("pi").get(k)[0][0];
					
				double [][] W_k_trans = MatrixOperations.transpose(W_k);
				double [][] Sigma = MatrixOperations.add(MatrixOperations.multiplication(W_k, W_k_trans),Psi);
							
				NormalDistribution mvGaussian = new NormalDistribution(my_k, Sigma);
				
				for(int i=0; i<n_observations; i++) {					
					x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);					
					r_ik[i][k] = pi_k[k][0]*mvGaussian.get_multivariateNormalPDF(x_i);
					r_i_sum[i] += r_ik[i][k];
				}
			}
					
			for(int k=0; k<n_analyzers; k++) {
				double r_sum = 0.0;
				for(int i=0; i<n_observations; i++) {	
					if(r_i_sum[i] == 0.0) {
						r_i_sum[i] = Double.MIN_VALUE;
					}
					r_ik[i][k] /= r_i_sum[i];
					r_sum += r_ik[i][k];
				}
				piList.add(new double [1][1]);
				piList.get(k)[0][0] = r_sum/n_observations;
			}
						
			for(int k=0; k<n_analyzers; k++ ) {
					
				ArrayList<double [][]> bList4Analyzer = new ArrayList<double [][]>();
				ArrayList<double [][]> CList4Analyzer = new ArrayList<double [][]>();
				ArrayList<double [][]> factor_cov_k_List4Analyzer = new ArrayList<double [][]>();
				
				my_k = pars.get("my").get(k);
				W_k  = pars.get("W").get(k);
				
				double [][] m_k  = new double [n_observations][n_factors];
				double [][] matrixTerm1 = new double [n_variables][(n_factors+1)];
				double [][] matrixTerm2 = new double [(n_factors+1)][(n_factors+1)];
				
				double [][] W_k_trans = MatrixOperations.transpose(W_k);		
				double [][] matrixProd = MatrixOperations.multiplication(W_k, W_k_trans);
				matrixProd = MatrixOperations.add(Psi, matrixProd);
				matrixProd = MatrixOperations.inverse(matrixProd);
				double [][] beta = MatrixOperations.multiplication(W_k_trans, matrixProd);
				double [][] const_4_factor_cov_k = MatrixOperations.substract(diag_L, MatrixOperations.multiplication(beta, W_k));

				for(int i=0; i<n_observations; i++) {				
					
					x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
					
					double [][] b_ik = new double [(n_factors+1)][1];
					b_ik[n_factors][0] = 1.0;
					
					double [][] C_ik = new double [(n_factors+1)][(n_factors+1)];
					C_ik[n_factors][n_factors] = 1.0;
					
					double [][] mean_adj_x_i = MatrixOperations.substract(x_i, my_k);
					double [][] m_ik = MatrixOperations.multiplication(beta, mean_adj_x_i);
					double [][] factor_cov_k = MatrixOperations.multiplication(m_ik, MatrixOperations.transpose(m_ik));
					factor_cov_k = MatrixOperations.add(const_4_factor_cov_k, factor_cov_k);
					
					C_ik = MatrixOperations.set_sub_matrix_to_matrix(C_ik, factor_cov_k, 0, (n_factors-1), 0, (n_factors-1));
					
					for(int l=0; l<n_factors; l++) {
						m_k[i][l] = m_ik[l][0]; 
						b_ik[l][0] = m_ik[l][0];
						C_ik[l][n_factors] = m_ik[l][0];
						C_ik[n_factors][l] = m_ik[l][0];
					}
						
					bList4Analyzer.add(b_ik);
					CList4Analyzer.add(C_ik);
					factor_cov_k_List4Analyzer.add(factor_cov_k);
					
					matrixTerm1 = MatrixOperations.add(matrixTerm1,MatrixOperations.scalar_multiplication(r_ik[i][k],MatrixOperations.multiplication(x_i, MatrixOperations.transpose(b_ik))));
					matrixTerm2 = MatrixOperations.add(matrixTerm2, MatrixOperations.scalar_multiplication(r_ik[i][k], C_ik));					
				}
				
				bList.add(bList4Analyzer);
				C_ik_List.add(CList4Analyzer);
	
				factors.add(m_k);
				factor_covariances.add(factor_cov_k_List4Analyzer);			
				
				matrixTerm2 = MatrixOperations.inverse(matrixTerm2);
				W_tilde.add(MatrixOperations.multiplication(matrixTerm1, matrixTerm2));
				
				myList.add(MatrixOperations.get_column_vec_from_matrix(W_tilde.get(k), n_factors));
				W.add(MatrixOperations.get_sub_matrix_between_row_and_col_idxs(W_tilde.get(k), 0, (n_variables-1), 0, (n_factors-1)));							
			}
	
			double [][] x_i_trans = new double [1][n_variables];
			double sigma = 0.0;
			for(int k=0; k<n_analyzers; k++) {
				
				double [][] matrixProdTerm3 = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(W_tilde.get(k)),Psi_inv),W_tilde.get(k));
				double [][] W_k_matrix = W.get(k);
				double [][] W_k_matrix_trans = MatrixOperations.transpose(W_k_matrix);
				double [][] W_prod = MatrixOperations.multiplication(W_k_matrix_trans, W_k_matrix);
				double [][] my_k_vec   = myList.get(k);
	
				for(int i=0; i<n_observations; i++) {

					//Sigma
					double [][] m_ik = MatrixOperations.get_sub_matrix_between_row_idxs(bList.get(k).get(i),0,(n_factors-1));
					double [][] factor_cov_ik = MatrixOperations.get_sub_matrix_between_row_and_col_idxs(C_ik_List.get(k).get(i), 0, (n_factors-1), 0, (n_factors-1));
				
					x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
					x_i_trans = MatrixOperations.transpose(x_i);
					double [][] mean_adj_x_i = MatrixOperations.substract(x_i, my_k_vec);
					double [][] mean_adj_x_i_trans = MatrixOperations.transpose(mean_adj_x_i);
					sigma += r_ik[i][k]*MatrixOperations.multiplication(mean_adj_x_i_trans, mean_adj_x_i)[0][0];
					sigma -= 2.0*r_ik[i][k]*MatrixOperations.multiplication(mean_adj_x_i_trans,MatrixOperations.multiplication(W_k_matrix, m_ik))[0][0];
					sigma += r_ik[i][k]*MatrixOperations.trace(MatrixOperations.multiplication(factor_cov_ik, W_prod));
			
	                //Log likelihood
					double [][] prod = MatrixOperations.multiplication(W_tilde.get(k), bList.get(k).get(i));
										
					prod =  MatrixOperations.substract(x_i,prod);	
					double [][] matrixProd = MatrixOperations.multiplication(x_i_trans, Psi_inv);
					newLogLikelihood -= MatrixOperations.multiplication(matrixProd, x_i)[0][0]/2.0;
					newLogLikelihood += MatrixOperations.multiplication(MatrixOperations.multiplication(matrixProd, W_tilde.get(k)),bList.get(k).get(i))[0][0];
				    
					matrixProd = MatrixOperations.multiplication(matrixProdTerm3, C_ik_List.get(k).get(i));
				    newLogLikelihood -= MatrixOperations.trace(matrixProd)/2.0;
				    newLogLikelihood *= r_ik[i][k];
				}
			}

			sigma /= (n_observations*n_variables);
					
			for(int i=0; i<n_variables; i++) {
				Psi[i][i] = sigma;
			}
			
			PsiList.add(Psi);
			pars.clear();
			pars.put("Psi", PsiList);
			pars.put("my", myList);
			pars.put("W", W);	
			pars.put("pi", piList);
			
			n_iterations_done = j+1;
			
			if(Math.abs(newLogLikelihood-logLikelihood)<=convergence_criterion){
				convergence = true;
				logLikelihood = newLogLikelihood;				
				break;	
			}else {
				logLikelihood = newLogLikelihood;
			}
			
		}		
		rotation_matrices = pars.get("W");
		my = pars.get("my");
	}
	
	
	public HashMap<String, ArrayList<double [][]>> get_input_pars_4_PPCA_em() {
		
		HashMap<String, ArrayList<double [][]>> input_pars = new HashMap<String, ArrayList<double [][]>>();
		
		ArrayList<double [][]> my = new ArrayList<double [][]>(n_analyzers);
		ArrayList<double [][]> W  = new ArrayList<double [][]>(n_analyzers);
		ArrayList<double [][]> pi = new ArrayList<double [][]>(n_analyzers);
		ArrayList<double [][]> Psi = new ArrayList<double [][]>(1);
		
		double [][] pi_k   = new double [1][1];
		double p = 1.0/(double) n_analyzers;
		
		for(int k=0; k<n_analyzers; k++) {

			double [][] my_k = new double [n_variables][1];
			double [][] W_k  = new double [n_variables][n_factors];
				
			//TODO: Use standard normal!
			NormalDistribution normalDist = new NormalDistribution(0.0, 0.1);
			for(int i=0; i<n_variables; i++) {
				my_k[i][0] = normalDist.sample()[0][0];
				for(int j=0; j<n_factors; j++) {
					W_k[i][j] = normalDist.sample()[0][0];
				}			
			}
			
			pi_k[0][0] = p;
 			
			my.add(my_k);
			W.add(W_k);
			pi.add(pi_k);		
		}
		
		double [][] Psi_matrix = new double [n_variables][n_variables];
		
		for(int i=0; i<n_variables; i++) {
			Psi_matrix[i][i] = 1.0;			
		}
		Psi.add(Psi_matrix);
		
		input_pars.put("my", my);
		input_pars.put("W", W);
		input_pars.put("pi", pi);
		input_pars.put("Psi", Psi);
		
		return input_pars;
	}
	
}
