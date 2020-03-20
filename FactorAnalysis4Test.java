package ComponentModels;

import java.util.ArrayList;
import java.util.HashMap;

import Distributions.NormalDistribution;
import Mathematics.MatrixOperations;

public class FactorAnalysis4Test extends FactorAnalysis{

	ArrayList<HashMap<String, ArrayList<double [][]>>> EMoutput;
	ArrayList<ArrayList<double [][]>> EMfactors;
	ArrayList<ArrayList<ArrayList<double[][]>>> EMfactorCovs;
	ArrayList<Double> trace_logLikelihoods;
	
	
	//Class only for testing: em_mixtureFactorAnalyzers4Test() modified for storing data from the iterations of EM.
	public FactorAnalysis4Test(double[][] X, int n_factors, int n_analyzers, boolean center, boolean scale) {
		super(X, n_factors, n_analyzers, center, scale);
	}

	
	public void em_mixtureFactorAnalyzers4Test() {
		
		HashMap<String, ArrayList<double [][]>> pars = new HashMap<String, ArrayList<double [][]>>();
		if(ext_init_values == null)	{
			pars = get_input_pars_4_em();
		}else {
			pars = ext_init_values;
		}
		
		//--- Here modification to org. algorithm because of em output after different steps
		EMoutput = new ArrayList<HashMap<String, ArrayList<double[][]>>>();
		EMfactors = new ArrayList<ArrayList<double [][]>>();
		EMfactorCovs = new ArrayList<ArrayList<ArrayList<double[][]>>>();	
		trace_logLikelihoods = new ArrayList<Double>();
		double [][] copyPsi = new double [n_variables][n_variables];
		for(int p=0; p<n_variables; p++) {
			copyPsi[p][p] = pars.get("Psi").get(0)[p][p];
		}
		//MatrixOperations.print_matrix(copyPsi);
		//MatrixOperations.print_matrix(pars.get("W").get(0));
		ArrayList<double [][]> PsiListCopy = new ArrayList<double [][]>();
		PsiListCopy.add(copyPsi);
		HashMap<String, ArrayList<double[][]>> parsCopy = new HashMap<String, ArrayList<double[][]>>();
		parsCopy.put("Psi", PsiListCopy);
		parsCopy.put("W", pars.get("W"));
		parsCopy.put("my", pars.get("my"));
		EMoutput.add(parsCopy);
		//---------------------------------------------------------------------------------------
		
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
			double [][] Psi_inv = MatrixOperations.inverse(Psi);
			double [][] Psi_sum = new double [n_variables][n_variables];
			
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
					double [][] factor_cov_k = MatrixOperations.multiplication(m_ik,MatrixOperations.transpose(m_ik));
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
			for(int k=0; k<n_analyzers; k++) {
				double [][] matrixProdTerm3 = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(W_tilde.get(k)),Psi_inv),W_tilde.get(k));
				for(int i=0; i<n_observations; i++) {
					
					x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
					x_i_trans = MatrixOperations.transpose(x_i);
					
					double [][] prod = MatrixOperations.multiplication(W_tilde.get(k), bList.get(k).get(i));
										
					prod =  MatrixOperations.substract(x_i,prod);
					Psi_sum = MatrixOperations.add(Psi_sum, MatrixOperations.scalar_multiplication(r_ik[i][k],MatrixOperations.multiplication(prod,x_i_trans)));
						
					double [][] matrixProd = MatrixOperations.multiplication(x_i_trans, Psi_inv);
					newLogLikelihood -= MatrixOperations.multiplication(matrixProd, x_i)[0][0]/2.0;
					newLogLikelihood += MatrixOperations.multiplication(MatrixOperations.multiplication(matrixProd, W_tilde.get(k)),bList.get(k).get(i))[0][0];
				    
					matrixProd = MatrixOperations.multiplication(matrixProdTerm3, C_ik_List.get(k).get(i));
				    newLogLikelihood -= MatrixOperations.trace(matrixProd)/2.0;
				    newLogLikelihood *= r_ik[i][k];
				}
			}

			double [][] Psi_copy = new double [n_variables][n_variables];
			for(int i=0; i<n_variables; i++) {
				Psi[i][i] = Psi_sum[i][i]/n_observations;
				Psi_copy[i][i] = Psi_sum[i][i]/n_observations;
			}
			
			PsiList.add(Psi_copy);
			pars = new HashMap<String, ArrayList<double[][]>>();
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
			
			trace_logLikelihoods.add(logLikelihood);
			copyPsi = new double [n_variables][n_variables];
			for(int p=0; p<n_variables; p++) {
				copyPsi[p][p] = pars.get("Psi").get(0)[p][p];
			}
			PsiList = new ArrayList<double [][]>();
			PsiList.add(copyPsi);
			pars.put("Psi", PsiList);
			EMoutput.add(pars);		
			EMfactors.add(factors);
			EMfactorCovs.add(factor_covariances);
		}		
		
		rotation_matrices = pars.get("W");
		my = pars.get("my");
		
		EMfactors.add(factors);
		EMfactorCovs.add(factor_covariances);
		
	}		
	
}
