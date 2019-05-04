package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Distributions.NormalDistribution;
import Mathematics.MatrixOperations;

public class ForwardBackwardAlgorithm extends HMM{
	
	public void runForwardAlgorithm(){
		
		filteredProbs = new ArrayList<List<Double>>(n_observations);
	
		double [][] T = MatrixOperations.transpose(transMatrix);
		
		List<Double> probs = new ArrayList<Double>(n_states);
		double summedProbs = 0.0;
				
		//t=1
		double [][] x_t = MatrixOperations.get_row_vec_from_matrix(observed_variables, 0);
		
		for(int s=0; s<n_states; s++){
			
			double [][] my_vec = get_my(s);
			double [][] sigma_matrix = get_sigma_matrix(s);
					
			NormalDistribution normal = new NormalDistribution(my_vec, sigma_matrix);
			obs_model_probs[s][0] = normal.get_multivariateNormalPDF(x_t);
			probs.add(s,obs_model_probs[s][0]*initProbs[s][0]);
				
			summedProbs = summedProbs + probs.get(s);
			
		}
		
		//Normalize probs
		for(int s=0; s<n_states; s++){
			probs.add(s,probs.get(s)/summedProbs);
		}
		
		filteredProbs.add(probs);
		
		//t=2,3,...,T
		for(int t=1; t<n_observations; t++){
			
			x_t = MatrixOperations.get_row_vec_from_matrix(observed_variables, t);
			
			double [][] probVec = get_filtered_probs(t-1);
			double [][] term = MatrixOperations.multiplication(T, probVec);
			summedProbs = 0.0;
			
			for(int s=0; s<n_states; s++){
				
				double [][] my_vec = get_my(s);
				double [][] sigma_matrix = get_sigma_matrix(s);
				
				NormalDistribution normal = new NormalDistribution(my_vec, sigma_matrix);
				obs_model_probs[s][t] = normal.get_multivariateNormalPDF(x_t);
				probs.add(s,obs_model_probs[s][t]*term[s][0]);
					
				summedProbs = summedProbs + probs.get(s);
				
			}
			
			//Normalize probs
			for(int s=0; s<n_states; s++){
				probs.add(s,probs.get(s)/summedProbs);
			}
			
			filteredProbs.add(probs);
			
		}
		
	}
	
	
	public void runBackwardAlgorithm(){
		
		smoothedProbs = new ArrayList<List<Double>>(n_observations);
		twoSliceMarginals = new ArrayList<List<Double>>(n_observations);
		
		double [][] betas = new double [n_states][1];
		
		//t=T
		for(int s=0; s<n_states; s++){
			betas[s][0] = 1.0;
		}
		
		//smoothedProbs = filteredProbs*betas
		smoothedProbs.add(n_observations-1,filteredProbs.get(n_observations-1));
		
		//t=T-1,T-2,...,1
		for(int t=1; t<n_observations; t++){
			
			List<Double> probs        = new ArrayList<Double>(n_states);
			List<Double> twoSliceMarg = new ArrayList<Double>(n_states);
			int idx1 = n_observations-1-t;
			int idx2 = 0;
			double summedProbs = 0.0;
					
			double [][] term = new double [n_states][1];
			
			for(int s=0; s<n_states; s++){				
				term[s][0] = obs_model_probs[s][idx1+1]*betas[s][0];	
				
				//Calculate the two-slice marginal matrices
				for(int k=0; k<n_states; k++){
					double eps_t_t_plus_1 = filteredProbs.get(idx1).get(k)*obs_model_probs[s][idx1+1]*betas[s][0]*transMatrix[k][s];
					twoSliceMarg.add(idx2,eps_t_t_plus_1);
				}				
			}
			
			//Set two-slice marignals at position T,T-1,T-2,...,2
			twoSliceMarginals.add(idx1+1, twoSliceMarg);
			
			betas = MatrixOperations.multiplication(transMatrix, term);
			
			for(int s=0; s<n_states; s++){
				probs.add(s,betas[s][0]*filteredProbs.get(idx1).get(s));
				summedProbs = summedProbs+probs.get(s);
			}
			
			//Normalize probs
			for(int s=0; s<n_states; s++){
				probs.add(s,probs.get(s)/summedProbs);
			}
			
			smoothedProbs.add(idx1,probs);
			
		}
			
	}
	
	
	public double [][] get_obs_model_probs(){
		
		return obs_model_probs;
		
	}
	
	
}
