package MixtureModels;

import java.util.ArrayList;
import java.util.List;

import DataManagement.InputDataManager;
import Distributions.NormalDistribution;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class MixtureOfGaussian {

	static InputDataManager inputData;
	
	static double [][] observed_variables;
	
	static int n_observations;
	static int n_variables;

	static int n_clusters;
	
	//Lists with number of slots eq. n_clusters
	static ArrayList<List<Double>> my;
	static ArrayList<List<Double>> Sigma;  //Vectorized sigma matrix
	static ArrayList<List<Double>> probs;
	
	//responsibilities (i=observation index, k=cluster index)
	static ArrayList<List<Double>> r_ik;
	
	static boolean useMAPEstimation;
	
	static int max_iterations;
	static double convergence_criterion;
	
	//EM stats
	static boolean convergence_reached = false;
	static int n_iterations;
	
	static double log_likelihood;
	
	//Hyperparameters for MAP estimation
	static List<Double> alpha;
	static double kappa_0;
	static double v_0;
	static double [][] m_0;
	static double [][] S_0;
	
	
	public static ArrayList<List<Double>> get_gaussian_mean_from_par_vec(double [][] par_vec){
		
		ArrayList<List<Double>> means = new ArrayList<List<Double>>(n_clusters);
		
		int idx=0;
		
		for(int c=0; c<n_clusters; c++){
			List<Double> mean = new ArrayList<Double>(n_variables);
			for(int i=0; i<n_variables; i++){
				mean.add(par_vec[idx][0]);
				idx++;
			}
			means.add(mean);
		}
					
		return means;
				
	}
	
	
	public static ArrayList<List<Double>> get_gaussian_sigma_elements_from_par_vec(double [][] par_vec){
				
		int n_elements = (int)Math.pow(n_variables, 2.0);
		
		ArrayList<List<Double>> sigma_vecs = new ArrayList<List<Double>>(n_clusters);
		
		int idx = n_variables*n_clusters;
		
		for(int c=0; c<n_clusters; c++){
			List<Double> sigma_vec = new ArrayList<Double>(n_elements);
			for(int i=0; i<n_elements; i++){
				sigma_vec.add(par_vec[idx][0]);
				idx++;
			}
			sigma_vecs.add(sigma_vec);
		}
		
		return sigma_vecs;
			
	}
	
	
	public static ArrayList<List<Double>> get_gaussian_densities(double [][] par_vec){
		
		ArrayList<List<Double>> means = get_gaussian_mean_from_par_vec(par_vec);
		ArrayList<List<Double>> sigmas = get_gaussian_sigma_elements_from_par_vec(par_vec);
		
		ArrayList<List<Double>> densities = new ArrayList<List<Double>>(n_clusters);
		
		for(int c=0; c<n_clusters; c++){
			
			List<Double> density = new ArrayList<Double>(n_observations);
			
			double [][] mean  = new double [n_variables][0]; 
			double [][] sigma = new double [n_variables][n_variables];
			
			int idx=0;
			
			for(int i=0; i<n_variables; i++){
				
				mean[i][0] = means.get(c).get(i);
				
				for(int j=0; j<n_variables; j++){					
					sigma[j][i]=sigmas.get(c).get(idx);					
				    idx++;										
				}
				
			}
			
			NormalDistribution normDist = new NormalDistribution(mean, sigma);
			
			for(int i=0; i<n_observations; i++){
				
				double [][] x = MatrixOperations.get_row_vec_from_matrix(observed_variables, i);
				
				density.add(normDist.get_multivariateNormalPDF(x));
				
			}
	
		}
			
		return densities;
		
	}
	
	
	public static ArrayList<List<Double>> get_gaussian_mixture_expectation(List<Double> cluster_probs, double [][] par_vec){
		
		ArrayList<List<Double>> responsibilities = new ArrayList<List<Double>>(n_clusters);
		
		ArrayList<List<Double>> gaussian_densities = get_gaussian_densities(par_vec);
        List<Double> densitySums = new ArrayList<Double>(n_observations);
		
        for(int i=0; i<n_observations; i++){       	
        	for(int c=0; c<n_clusters; c++){           	
            	densitySums.add(cluster_probs.get(c)*gaussian_densities.get(c).get(i));            	
            }        	
        }
  
		for(int c=0; c<n_clusters; c++){
			
			List<Double> responsibility = new ArrayList<Double>(n_observations);			
			
			for(int i=0; i<n_observations; i++){				
				responsibility.add(cluster_probs.get(c)*gaussian_densities.get(c).get(i)/densitySums.get(i));				
			}
			
			responsibilities.add(responsibility);
			
		}
		
		return responsibilities;
		
	}
	
	
	public static ArrayList<ArrayList<List<Double>>> get_gaussian_mixture_maximization_results(ArrayList<List<Double>> responsibilities, double [][] par_vec){
		
		//1st slot probs, 2nd slot mean vectors, 3rd slot sigma matrices (vectorized)
		ArrayList<ArrayList<List<Double>>> max_results = new ArrayList<ArrayList<List<Double>>>(3);
		ArrayList<List<Double>> means  = new ArrayList<List<Double>>(n_clusters);
		ArrayList<List<Double>> sigmas = new ArrayList<List<Double>>(n_clusters);
		ArrayList<List<Double>> probs  = new ArrayList<List<Double>>(n_clusters);
		
		List<Double> r_k = new ArrayList<Double>(n_clusters);
		
		for(int c=0; c<n_clusters; c++){			
			r_k.add(GeneralMath.sumDblList(responsibilities.get(c)));	
			
			List<Double> prob = new ArrayList<Double>(1);
			prob.add(r_k.get(c)/n_observations);
			probs.add(prob);
					
		}
		
		max_results.add(probs);
		
		for(int c=0; c<n_clusters; c++){
			
			List<Double> mean_vec = new ArrayList<Double>(n_variables);
			List<Double> sigma_vec = new ArrayList<Double>((int)Math.pow(n_variables,2.0));
			
			double [][] x_prod  = new double [n_variables][n_variables];
			double [][] my      = new double [n_variables][1];
			double [][] my_prod = new double [n_variables][n_variables];
			double [][] sigma_matrix = new double [n_variables][n_variables];
			
			for(int i=0; i<n_observations; i++){
				
				double [][] x  = MatrixOperations.get_row_vec_from_matrix(observed_variables, i);	
				
				if(i==0){					
					for(int j=0; j<n_variables; j++){
						my_prod[j][0] = x[j][0]/r_k.get(c);
						mean_vec.add(my_prod[j][0]);
					}					
				}else{					
					for(int j=0; j<n_variables; j++){
						my_prod[j][0] = mean_vec.get(j)+x[j][0]/r_k.get(c);
						mean_vec.add(j,my_prod[j][0]);
					}					
				}
				
				double r = responsibilities.get(c).get(i)/r_k.get(c);
				
				x_prod = MatrixOperations.add(x_prod,MatrixOperations.scalar_multiplication(r,MatrixOperations.multiplication(x, MatrixOperations.transpose(x))));
				
			}
			
			means.add(mean_vec);
			
			my_prod      = MatrixOperations.multiplication(my, MatrixOperations.transpose(my));			
			sigma_matrix = MatrixOperations.substract(x_prod, my_prod);
			
			for(int i=0; i<n_variables; i++){				
				for(int j=0; j<n_variables; j++){					
					sigma_vec.add(sigma_matrix[j][i]);									
				}				
			}
			
			sigmas.add(sigma_vec);
			
		}
		
		max_results.add(means);
		max_results.add(sigmas);
		
		return max_results;
		
	}
	
	
	public static ArrayList<ArrayList<List<Double>>> get_gaussian_mixture_maximization_results_4_MAP_estimation(ArrayList<List<Double>> responsibilities, double [][] par_vec){
		
		//1st slot probs, 2nd slot mean vectors, 3rd slot sigma matrices (vectorized)
		ArrayList<ArrayList<List<Double>>> max_results = new ArrayList<ArrayList<List<Double>>>(3);
		ArrayList<List<Double>> means  = new ArrayList<List<Double>>(n_clusters);
		ArrayList<List<Double>> sigmas = new ArrayList<List<Double>>(n_clusters);
		ArrayList<List<Double>> probs  = new ArrayList<List<Double>>(n_clusters);
		
		List<Double> r_k = new ArrayList<Double>(n_clusters);
		
		double alphaSum = GeneralMath.sumDblList(alpha);
		
		for(int c=0; c<n_clusters; c++){			
			r_k.add(GeneralMath.sumDblList(responsibilities.get(c)));	
			
			List<Double> prob = new ArrayList<Double>(1);
			prob.add((r_k.get(c)+alpha.get(c)-1)/(n_observations+alphaSum-n_clusters));
			probs.add(prob);
					
		}
		
		max_results.add(probs);
	
		for(int c=0; c<n_clusters; c++){
			
			List<Double> mean_vec = new ArrayList<Double>(n_variables);
			List<Double> sigma_vec = new ArrayList<Double>((int)Math.pow(n_variables,2.0));
			
			double [][] x_mean  = new double [n_variables][1];
			double [][] my_prod = new double [n_variables][n_variables];
			double [][] sigma_matrix = new double [n_variables][n_variables];
			
			for(int i=0; i<n_observations; i++){
				
				double [][] x  = MatrixOperations.get_row_vec_from_matrix(observed_variables, i);	
				
				if(i==0){					
					for(int j=0; j<n_variables; j++){
						my_prod[j][0] = x[j][0]/r_k.get(c);
						mean_vec.add(my_prod[j][0]);
					}					
				}else{					
					for(int j=0; j<n_variables; j++){
						my_prod[j][0] = mean_vec.get(j)+x[j][0]/r_k.get(c);
						mean_vec.add(j,my_prod[j][0]);
					}					
				}
				
			}
			
			for(int i=0; i<n_variables; i++){
				double my4MAP = (r_k.get(c)*mean_vec.get(i)+kappa_0*m_0[i][0])/(r_k.get(c)+kappa_0);
				mean_vec.add(i,my4MAP);
				x_mean[i][0] = my4MAP;
			}
			
			means.add(mean_vec);
			
			double [][] S_k = new double [n_variables][n_variables];
			double [][] x_diff = new double [n_variables][1];
			double [][] x_prod = new double [n_variables][n_variables];
			
			for(int i=0; i<n_observations; i++){
				
				double [][] x  = MatrixOperations.get_row_vec_from_matrix(observed_variables, i);
				double r       = responsibilities.get(c).get(i);
				
				x_diff = MatrixOperations.substract(x, x_mean);
				x_prod = MatrixOperations.scalar_multiplication(r,MatrixOperations.multiplication(x_diff, MatrixOperations.transpose(x_diff)));
				
				for(int j=0; j<n_variables; j++){
					S_k = MatrixOperations.add(S_k, x_prod);
				}
				
			}

			double [][] x_adj = MatrixOperations.substract(x_mean, m_0);
			x_adj = MatrixOperations.multiplication(x_adj, MatrixOperations.transpose(x_adj));
			x_adj = MatrixOperations.scalar_multiplication((kappa_0*r_k.get(c))/(kappa_0+r_k.get(c)), x_adj);
			
			sigma_matrix = MatrixOperations.add(MatrixOperations.add(S_0, S_k),x_adj);
			sigma_matrix = MatrixOperations.scalar_multiplication((v_0+r_k.get(c)+n_variables+2), sigma_matrix);
			
			for(int i=0; i<n_variables; i++){				
				for(int j=0; j<n_variables; j++){					
					sigma_vec.add(sigma_matrix[j][i]);									
				}				
			}
			
			sigmas.add(sigma_vec);
			
		}
		
		max_results.add(means);
		max_results.add(sigmas);
			
		return max_results;
		
	}
	
	
	public static double get_gaussian_mixture_log_likelihood(ArrayList<ArrayList<List<Double>>> maximization_results){
		
		double logLik = 0.0;
		
		for(int i=0; i<n_observations; i++){
			
			for(int c=0; c<n_clusters; c++){
				
				double [][] my    = new double [n_variables][1];
				double [][] sigma = new double [n_variables][n_variables];
				int idx = 0;
				
				for(int j=0; j<n_variables; j++){
					my[j][0] = maximization_results.get(1).get(c).get(j);
					for(int k=0; k<n_variables; k++){
						sigma[k][j] = maximization_results.get(2).get(c).get(idx);
						idx++;
					}
				}
				
				NormalDistribution normDist = new NormalDistribution(my, sigma);
				
				double [][] x  = MatrixOperations.get_row_vec_from_matrix(observed_variables, i);
				double prob    = maximization_results.get(0).get(c).get(0);
				
				logLik = prob*normDist.get_multivariateNormalPDF(x);
				
			}
			
		}
		
		return logLik;
		
	}
	
	
	public static double [][] get_par_vec_from_maximization_result(ArrayList<ArrayList<List<Double>>> maximization_results){
		
		int n = n_clusters*(n_variables + (int)Math.pow(n_variables,2.0));
		
		double [][] par_vec = new double [n][1];
		
		int idx = 0;
		
		for(int c=0; c<n_clusters; c++){
			for(int i=0; i<n_variables; i++){
				par_vec[idx][0] = maximization_results.get(1).get(c).get(i);
			}	
		}
		
		for(int c=0; c<n_clusters; c++){
			int sigmaIdx = 0;
			for(int i=0; i<n_variables; i++){
				for(int j=0; j<n_variables; j++){
					par_vec[idx][0] = maximization_results.get(1).get(c).get(sigmaIdx);
				}				
			}	
		}
		
		return par_vec;
		
	}
	
	
	public static List<Double> get_cluster_probs_from_maximization_result(ArrayList<ArrayList<List<Double>>> maximization_results){
		
		List<Double> probs = new ArrayList<Double>(n_clusters);
		
		for(int c=0; c<n_clusters; c++){
			probs.add(maximization_results.get(0).get(c).get(0));
		}
		
		return probs;
		
	}
	
	
	public static void expectation_maximization_algorithm(){
		
		//Initial parameters
		List<Double> cluster_probs = new ArrayList<Double>(n_clusters);
		double [][] par_vec = new double [(n_variables+(int)Math.pow(n_variables, 2.0))][1];
			
		ArrayList<List<Double>> resp = null;
		ArrayList<ArrayList<List<Double>>> maxRes = null;
		
		double prevlogLik = Double.MIN_VALUE;
		double newlogLik  = 0.0;
		
		for(int i=0; i<max_iterations; i++){
			
			resp = get_gaussian_mixture_expectation(cluster_probs, par_vec);
			
			if(useMAPEstimation == false){
				maxRes = get_gaussian_mixture_maximization_results(resp, par_vec);
			}else{
				maxRes = get_gaussian_mixture_maximization_results_4_MAP_estimation(resp, par_vec);
			}
					
			par_vec = get_par_vec_from_maximization_result(maxRes);
			
			cluster_probs = get_cluster_probs_from_maximization_result(maxRes);
			
			newlogLik = get_gaussian_mixture_log_likelihood(maxRes);
			
			if(Math.abs(prevlogLik-newlogLik)<=convergence_criterion){
				System.out.println("EM algorithm for Gaussian mixture model has reached convergence after " + i + " iterations.");
				convergence_reached = true;
				n_iterations = i;					
				break;
			}
			
		}
		
		if(convergence_reached == false){
			n_iterations = max_iterations;
		}
		
		log_likelihood = newlogLik;
		
		probs = maxRes.get(0);
		my    = maxRes.get(1);
		Sigma = maxRes.get(2);
		r_ik  = resp;
			
	}
	
	
	public static double [][] get_my(int clusterNumber){
		
		if(clusterNumber>(n_clusters-1)){
			throw new RuntimeException("Invalid cluster number supplied for getting the mean of the GMM");
		}
		
		double [][] my_vec = new double [n_variables][1];
		
		for(int i=0; i<n_variables; i++){
			my_vec[i][0] = my.get(clusterNumber).get(i);
		}
		
		return my_vec;
		
	}
	
	
	public static double [][] get_sigma_matrix(int clusterNumber){
		
		if(clusterNumber>(n_clusters-1)){
			throw new RuntimeException("Invalid cluster number supplied for getting the sigma matrix of the GMM");
		}
		
		double [][] sigma_matrix = new double [n_variables][n_variables];
		
		int idx=0;
		
		for(int i=0; i<n_variables; i++){
			for(int j=0; j<n_variables; j++){
				sigma_matrix[j][i] = Sigma.get(clusterNumber).get(idx);
				idx++;
			}
			
		}
		
		return sigma_matrix;
		
	}
	
	
	public static double get_cluster_prob(int clusterNumber){
		
		if(clusterNumber>(n_clusters-1)){
			throw new RuntimeException("Invalid cluster number supplied for getting the cluster probs of the GMM");
		}
		
		return probs.get(clusterNumber).get(0);
		
	}
	
	
	//getter for r_ik = p(z_i = k | x_i, ...)
	public static double [][] get_responsibilities(int clusterNumber){
		
		if(clusterNumber>(n_clusters-1)){
			throw new RuntimeException("Invalid cluster number supplied for getting the responsibilities of the GMM");
		}
		
		double [][] resp = new double [n_observations][1];
		
		for(int i=0; i<n_observations; i++){			
			resp[i][0] = r_ik.get(clusterNumber).get(i);			
		}
		
		return resp;
		
	}
	
	
	public static void use_MAP_estimation(boolean useMAP){
		
		useMAPEstimation = useMAP;
		
	}
	
	
	public static void set_default_hyperparameters_4_MAP_est(){
		
		kappa_0 = 0.0;
		v_0     = n_variables/2.0;
		m_0     = new double [n_variables][1];
		
		for(int c=0; c<n_clusters; c++){
			alpha.add(1.0);
		}
		
		double [] s       = new double [n_variables];
		
		for(int i=0; i<n_variables; i++){
			double [][] x = MatrixOperations.get_column_vec_from_matrix(observed_variables, i);
			double x_means = GeneralMath.mean(x);
			
			for(int j=0; j<n_observations; j++){
				s[i] = s[i] + Math.pow(x[j][0]-x_means,2.0);
			}
			
			s[i] = s[i]/n_observations;
			s[i] = s[i]/Math.pow(n_clusters, 1.0/n_variables);
			
		}
		
		S_0 = MatrixOperations.diagonal(s);
				
	}
	
		
	@SuppressWarnings("static-access")
	public static void read_input_data(String fileName, boolean hasRowNames, boolean hasColNames) throws Exception{
		
		inputData = new InputDataManager();
		
		inputData.fileReader(fileName, false, hasRowNames, hasColNames);
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void select_input_data(String [] rownames, String [] colnames, String ref_class){
		
		inputData.selectLoadedData(rownames, colnames);
		
		observed_variables = inputData.selectedDblFileData;		
		n_observations     = observed_variables.length;
		n_variables        = observed_variables[0].length;
		
	}
	
	
	public static void remove_input_data_manager(){	
		
		inputData = null;		
		
	}
	
}
