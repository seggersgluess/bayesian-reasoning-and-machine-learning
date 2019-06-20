package HiddenMarkovModels;

import java.util.ArrayList;
import java.util.List;

import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class ForwardBackwardAlgorithm extends HMM{
	
	public void runForwardAlgorithm(){
		
		int n_obs = obs_model_probs[0].length;
		
		filteredProbs = new ArrayList<List<Double>>(n_obs);
	
		double [][] T = MatrixOperations.transpose(transMatrix);
		
		List<Double> probs = new ArrayList<Double>(n_states);
		double summedProbs = 0.0;
				
		//t=1
		for(int s=0; s<n_states; s++){			
			probs.add(s,obs_model_probs[s][0]*initProbs[s][0]);				
			summedProbs += probs.get(s);		
		}
		
		//Normalize probs
		for(int s=0; s<n_states; s++){
			probs.set(s,probs.get(s)/summedProbs);
		}
		
		filteredProbs.add(probs);
		
		//t=2,3,...,T
		for(int t=1; t<n_obs; t++){

			probs = new ArrayList<Double>(n_states);
			
			double [][] probVec = get_filtered_probs(t-1);
			double [][] term = MatrixOperations.multiplication(T, probVec);
			summedProbs = 0.0;
			
			for(int s=0; s<n_states; s++){				
				probs.add(s,obs_model_probs[s][t]*term[s][0]);		
				summedProbs += probs.get(s);			
			}
							
			//Normalize probs
			for(int s=0; s<n_states; s++){
				probs.set(s,probs.get(s)/summedProbs);
			}
			
			filteredProbs.add(probs);
			
		}
		
	}
	
	
	public void runBackwardAlgorithm(){
		
		int n_obs = obs_model_probs[0].length;
		
		smoothedProbs = new ArrayList<List<Double>>(n_obs);
		twoSliceMarginals = new ArrayList<List<Double>>(n_obs);
		
		double [][] betas = new double [n_states][1];
		
		//t=T
		for(int s=0; s<n_states; s++){
			betas[s][0] = 1.0;
		}
		
		//smoothedProbs = filteredProbs*betas
		smoothedProbs.add(filteredProbs.get(n_obs-1));
		
		//t=T-1,T-2,...,1
		for(int t=1; t<n_obs; t++){
			
			List<Double> probs        = new ArrayList<Double>(n_states);
			List<Double> twoSliceMarg = new ArrayList<Double>(n_states);
			int idx1 = n_obs-1-t;
			int idx2 = idx1+1;
			double summedProbs = 0.0;
					
			double [][] term = new double [n_states][1];
			
			for(int s=0; s<n_states; s++){				
				term[s][0] = obs_model_probs[s][idx2]*betas[s][0];	
				
				//Calculate the two-slice marginal matrices
				for(int k=0; k<n_states; k++){
					double eps_t_t_plus_1 = filteredProbs.get(idx1).get(k)*obs_model_probs[s][idx2]*betas[s][0]*transMatrix[k][s];
					twoSliceMarg.add(k,eps_t_t_plus_1);
				}				
			}
			
			//Set two-slice marignals at position T,T-1,T-2,...,2
			//twoSliceMarginals.add(idx2, twoSliceMarg);
			
			//---> Achtung: Stimmen idxs der Liste für Folgeverarbeitung? Checken! 
			//              (Werte werden jetzt einfach nach rechts verschoben).
			twoSliceMarginals.add(0, twoSliceMarg);
			//
			
			betas = MatrixOperations.multiplication(transMatrix, term);
			double summedBetas = GeneralMath.sum(betas);	
			
			for(int s=0; s<n_states; s++){
				probs.add(s,betas[s][0]*filteredProbs.get(idx1).get(s));
				summedProbs = summedProbs+probs.get(s);
			}
			
			//Normalize probs
			for(int s=0; s<n_states; s++){
				probs.set(s,probs.get(s)/summedProbs);
				//Normalizing betas
				betas[s][0] = betas[s][0]/summedBetas;
			}
			
			//every element is shifted to the right
			smoothedProbs.add(0,probs);
						
		}
			
	}
		
}
