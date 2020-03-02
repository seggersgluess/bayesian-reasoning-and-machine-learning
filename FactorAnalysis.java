package ComponentModels;

import java.util.ArrayList;
import java.util.HashMap;

import Distributions.NormalDistribution;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class FactorAnalysis extends ComponentModels{

	int n_factors = 1; 
	int n_analyzers = 1;

	ArrayList<double [][]> factors;
	ArrayList<double [][]> factor_covariances;
	
	ArrayList<double [][]> rotation_matrices;
	ArrayList<double [][]> my;	
	
	double [][] Psi;
	
	double [][] pi_k;
	double [][] r_ik;
	
	int iterations = 500;
	boolean convergence = false;
	double convergence_criterion = 1e-04;
	
	double logLikelihood = Double.MIN_VALUE;
	
	
	public FactorAnalysis(double [][] X, int n_factors,  boolean scale) {
		
		super(X, scale);
		
		int n_vars = X[0].length;
		
		if(n_factors<1 || n_factors>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		this.n_factors = n_factors;

	}
	
	
	public FactorAnalysis(double [][] X, int n_factors, boolean center, boolean scale) {
		
		super(X, center, scale);
		
		int n_vars = X[0].length;

		if(n_factors<1 || n_factors>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		this.n_factors = n_factors;
	
	}
	
	
	public FactorAnalysis(double [][] X, int n_factors, int n_analyzers, boolean center, boolean scale) {
		
		super(X, center, scale);
		
		int n_vars = X[0].length;

		if(n_factors<1 || n_factors>n_vars) {
			throw new RuntimeException("Invalid number of principal components supplied. Only values between 1 and " + n_vars + " allowed.");
		}
		
		if(n_analyzers < 1) {
			throw new RuntimeException("Number of factor analyzers has to be larger than 1.");
		}
		
		this.n_factors = n_factors;
	    this.n_analyzers = n_analyzers;
		
	}
	
	
	public void em_mixtureFactorAnalyzers() {
		
		factors = new ArrayList<double [][]>(n_analyzers);
		factor_covariances = new ArrayList<double [][]>(n_analyzers);
		
		HashMap<String, ArrayList<double [][]>> pars = get_input_pars_4_em();
				
		double [][] x_i  = new double [n_variables][1];
		
		double [][] my_k = new double [n_variables][1];
		double [][] W_k  = new double [n_variables][n_factors];
		
		pi_k = new double [n_factors][1];
		r_ik = new double [n_observations][n_analyzers];
		double []   r_i_sum = new double [n_observations];
		
		double [][] diag_L = MatrixOperations.identity(n_factors);
		double [][] m_k  = new double [n_observations][n_factors];	
		
		double [][] b_ik = new double [(n_factors+1)][1];
		b_ik[n_factors][0] = 1.0;
		
		double [][] C_ik = new double [(n_factors+1)][(n_factors+1)];
		C_ik[n_factors][n_factors] = 1.0;
		
		for(int j=0; j<iterations; j++) {	
			
			factors.clear();
			factor_covariances.clear();

			ArrayList<double [][]> myList = new ArrayList<double [][]>(n_analyzers);
			ArrayList<double [][]> W      = new ArrayList<double [][]>(n_analyzers);
			ArrayList<double [][]> piList  = new ArrayList<double [][]>(n_analyzers);
			ArrayList<double [][]> PsiList = new ArrayList<double [][]>(1);
				
			ArrayList<ArrayList<double [][]>> bList  = new ArrayList<ArrayList<double [][]>>();
			ArrayList<ArrayList<double [][]>> C_ik_List  = new ArrayList<ArrayList<double [][]>>();
			
			bList.add(new ArrayList<double [][]>());
			ArrayList<double [][]> W_tilde = new ArrayList<double [][]>(n_analyzers);
			
			Psi = pars.get("Psi").get(0);
			double [][] Psi_inv = MatrixOperations.inverse(Psi);
			double [][] Psi_sum = new double [n_variables][n_variables];
			
			double newLogLikelihood = Math.log(2.0*Math.PI)*n_observations*n_variables/2.0;
			newLogLikelihood -= n_observations*MatrixOperations.determinant(Psi)/2.0;
			
			for(int k=0; k<n_analyzers; k++) {
	
				my_k = pars.get("my").get(k);
				W_k  = pars.get("W").get(k);
				pi_k[k][0] = pars.get("pi").get(k)[0][0];
					
				double [][] W_k_trans = MatrixOperations.transpose(W_k);
				double [][] Sigma = MatrixOperations.add(MatrixOperations.multiplication(W_k, W_k_trans),Psi);
				
				NormalDistribution mvGaussian = new NormalDistribution(my_k, Sigma);
				
				for(int i=0; i<n_observations; i++) {
					if(k==0) {
						x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
					}
					r_ik[i][k] = pi_k[k][0]*mvGaussian.get_multivariateNormalPDF(x_i);
					r_i_sum[i] += r_ik[i][k];
				}
			}
					
			double r_sum = 0.0;
			for(int k=0; k<n_analyzers; k++) {
				for(int i=0; i<n_observations; i++) {							
					r_ik[i][k] /= r_i_sum[i];
					r_sum += r_ik[i][k];
				}
				piList.add(new double [1][1]);
				piList.get(k)[0][0] = r_sum/n_observations;
			}
			
			for(int k=0; k<n_analyzers; k++ ) {
								
				my_k = pars.get("my").get(k);
				W_k  = pars.get("W").get(k);
				
				double [][] matrixTerm1 = new double [n_variables][(n_factors+1)];
				double [][] matrixTerm2 = new double [(n_factors+1)][(n_factors+1)];
				
				double [][] W_k_trans = MatrixOperations.transpose(W_k);	
				double [][] matrixProd =  MatrixOperations.multiplication(W_k_trans, Psi_inv);
				double [][] factor_cov_k = MatrixOperations.add(diag_L,MatrixOperations.multiplication(matrixProd,W_k));
				factor_cov_k = MatrixOperations.inverse(factor_cov_k);
				
				C_ik = MatrixOperations.set_sub_matrix_to_matrix(C_ik, factor_cov_k, 0, (n_factors-1), 0, (n_factors-1));
				
				for(int i=0; i<n_observations; i++) {				
					if(k==0) {
						x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
					}
					
					double [][] m_ik = MatrixOperations.multiplication(factor_cov_k, MatrixOperations.multiplication(matrixProd,MatrixOperations.substract(x_i, my_k)));				
					
					for(int l=0; l<n_factors; l++) {
						m_k[i][l] = m_ik[l][0]; 
						b_ik[l][0] = m_ik[l][0];
						C_ik[l][n_factors] = m_ik[l][0];
						C_ik[n_factors][l] = m_ik[l][0];
					}
					
					factors.add(m_k);
					factor_covariances.add(factor_cov_k);
					bList.get(i).add(b_ik);
					C_ik_List.get(i).add(C_ik);
					
					matrixTerm1 = MatrixOperations.add(matrixTerm1,MatrixOperations.scalar_multiplication(r_ik[i][k],MatrixOperations.multiplication(x_i, MatrixOperations.transpose(b_ik))));
					matrixTerm2 = MatrixOperations.add(matrixTerm2, MatrixOperations.scalar_multiplication(r_ik[i][k], C_ik));					
				}

				matrixTerm2 = MatrixOperations.inverse(matrixTerm2);
				W_tilde.add(MatrixOperations.multiplication(matrixTerm1, matrixTerm2));
				
				myList.add(MatrixOperations.get_column_vec_from_matrix(W_tilde.get(k), n_factors));
				W.add(MatrixOperations.get_sub_matrix_between_row_and_col_idxs(W_tilde.get(k), 0, (n_variables-1), 0, (n_factors-1)));
					
			}
						
			double [][] x_i_trans = new double [1][n_variables];
			for(int k=0; k<n_analyzers; k++) {
				double [][] matrixProdTerm3 = MatrixOperations.multiplication(MatrixOperations.multiplication(MatrixOperations.transpose(W_tilde.get(k)),Psi_inv),W_tilde.get(k));
				for(int i=0; i<n_observations; i++) {
					if(k==0) {
						x_i = MatrixOperations.get_row_vec_from_matrix(X_scaled, i);
						x_i_trans = MatrixOperations.transpose(x_i);
					}
					Psi_sum = MatrixOperations.add(Psi_sum, MatrixOperations.scalar_multiplication(r_ik[i][k],MatrixOperations.multiplication(MatrixOperations.substract(x_i,MatrixOperations.multiplication(W_tilde.get(k), bList.get(i).get(k))),x_i_trans)));
					
					double [][] matrixProd = MatrixOperations.multiplication(x_i_trans, Psi_inv);
					newLogLikelihood -= MatrixOperations.multiplication(matrixProd, x_i)[0][0]/2.0;
					newLogLikelihood += MatrixOperations.multiplication(MatrixOperations.multiplication(matrixProd, W_tilde.get(k)),bList.get(i).get(k))[0][0];
				    
					matrixProd = MatrixOperations.multiplication(matrixProdTerm3, C_ik_List.get(i).get(k));
				    newLogLikelihood -= MatrixOperations.trace(matrixProd)/2.0;
				    newLogLikelihood *= r_ik[i][k];
				}
			}

			for(int i=0; i<n_variables; i++) {
				Psi[i][i] = Psi_sum[i][i]/n_observations;
			}
			
			PsiList.add(Psi);
			//TODO: Check clearing!
			pars.clear();
			pars.put("Psi", PsiList);
			pars.put("my", myList);
			pars.put("W", W);	
			pars.put("pi", piList);
			
			if((newLogLikelihood-logLikelihood)<=convergence_criterion){
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
	
	
	private HashMap<String, ArrayList<double [][]>> get_input_pars_4_em() {
		
		HashMap<String, ArrayList<double [][]>> input_pars = new HashMap<String, ArrayList<double [][]>>();
		
		ArrayList<double [][]> my = new ArrayList<double [][]>(n_analyzers);
		ArrayList<double [][]> W  = new ArrayList<double [][]>(n_analyzers);
		ArrayList<double [][]> pi = new ArrayList<double [][]>(n_analyzers);
		ArrayList<double [][]> Psi = new ArrayList<double [][]>(1);
		
		double [][] pi_k   = new double [1][1];
		double p = 1.0/(double) n_analyzers;
		
		for(int k=0; k<n_analyzers; k++) {

			//TODO: How to set initial mu_i, W_i?
			double [][] my_k = new double [n_variables][1];
			double [][] W_k  = new double [n_variables][n_factors];
				
			pi_k[0][0] = p;
 			
			my.add(my_k);
			W.add(W_k);
			pi.add(pi_k);		
		}
		
		double [][] Psi_matrix = new double [n_variables][n_variables];
		
		for(int i=0; i<n_variables; i++) {
			Psi_matrix[i][i] = GeneralMath.variance(MatrixOperations.get_column_from_matrix(X_scaled, i));			
		}

		input_pars.put("my", my);
		input_pars.put("W", W);
		input_pars.put("pi", pi);
		input_pars.put("Psi", Psi);
		
		return input_pars;
	}
	
	
	public void calc_rotated_X() {
		
		if(factors == null) {
			throw new RuntimeException("Factor Analysis not done yet. No factors found.");
		}
		
		rotated_X = new double [n_observations][n_variables];
		
		for(int i=0; i<n_observations; i++) {
			for(int k=0; k<n_analyzers; k++) {
				double [][] factor = MatrixOperations.get_row_vec_from_matrix(factors.get(k), i);
				double [][] x_ik = MatrixOperations.scalar_multiplication(r_ik[i][k],MatrixOperations.add(my.get(k), MatrixOperations.multiplication(rotation_matrices.get(k), factor)));
				for(int j=0; j<n_variables; j++) {
					rotated_X[i][j] = x_ik[j][0];
				}
			}
		}	
	}
	
	//k=0,1,...,n_analyzers-1
	public double [][] get_roation_matrix(int k) {
		if(rotation_matrices == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated rotation matrix found.");
		}
		if(k<1 || k>=n_analyzers) {
			throw new RuntimeException("Invalid no. of factor analyzer supplied.");
		}
		return rotation_matrices.get(k);
	}
	
	
	public double [][] get_my(int k) {
		if(my == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated my vector found.");
		}
		if(k<1 || k>=n_analyzers) {
			throw new RuntimeException("Invalid no. of factor analyzer supplied.");
		}
		return my.get(k);
	}
	
	
	public double [][] get_factors(int k) {
		if(factors == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated factors found.");
		}
		if(k<1 || k>=n_analyzers) {
			throw new RuntimeException("Invalid no. of factor analyzer supplied.");
		}
		return factors.get(k);
	}
	
	
	public double [][] get_factor_covariance(int k) {
		if(factor_covariances == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated factor covariance found.");
		}
		if(k<1 || k>=n_analyzers) {
			throw new RuntimeException("Invalid no. of factor analyzer supplied.");
		}
		return factor_covariances.get(k);
	}
	
	
	public double [][] get_psi_matrix() {
		if(Psi == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated Psi matrix found.");
		}
		return Psi;
	}
	

	public double [][] get_pi_k(){
		if(pi_k == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated probability vector pi_k found.");
		}
		return pi_k;
	}
	
	
	public double [][] get_r_ik(){
		if(r_ik == null) {
			throw new RuntimeException("Factor Analysis not done yet. No calculated probability vector r_ik found.");
		}
		return r_ik;
	}
	
	
	public int get_number_of_factor_analyzers() {
		return n_analyzers;
	}
	
	
	public int get_number_of_factors() {
		return n_factors;
	}
	
	
	public boolean isConverged() {
		return convergence;
	}
	
	
	public void set_iterations4EM(int iterations) {
		if(iterations < 1) {
			throw new RuntimeException("Invalid iterations for EM supplied.");
		}
		this.iterations = iterations;
	}
	
	
	public void set_convergence_criterion4EM(double criterion) {
		this.convergence_criterion = criterion;
	}
		
}
