package HiddenMarkovModels;

public class Viterbi extends HMM{

	static double [][] deltas;
	static int [][] alphas;
	
	
	public void runViterbi(){
		
		stateSequence  = new int [n_observations][1];
		stateProbs = new double [n_observations][n_states];
		
		calc_alphas_and_deltas();

		double delta = deltas[n_observations-1][0];
		
		for(int s=1; s<n_states; s++){
			
			if(deltas[n_observations-1][s]>=delta){
				delta = deltas[n_observations-1][s];
				stateSequence[n_observations-1][0] = s;
				stateProbs[n_observations-1][s] = 1.0;
			}
			
		}
		
		for(int t=1; t<n_observations; t++){
			
			int idx = n_observations-1-t;
			stateSequence[idx][0] = alphas[idx+1][stateSequence[idx+1][0]];
			stateProbs[idx][stateSequence[idx][0]] = 1.0;			
		}
			
	}
	
	
	private void calc_alphas_and_deltas(){
		
		deltas = new double [n_observations][n_states];
		alphas = new int [n_observations][n_states];
		
		//Initialization t=1
		for(int s=0; s<n_states; s++){
			deltas[0][s] = initProbs[s][0]*obs_model_probs[s][0];
		}
		
		//t=2,3,...,T
		for(int t=1; t<n_observations; t++){
			
			for(int s=0; s<n_states; s++){
				
				deltas[t][s] = deltas[t-1][0]*transMatrix[0][s]*obs_model_probs[s][t];
				alphas[t][s] = 0;
				
				for(int k=1; k<n_states; k++){
					
					double delta_new = deltas[t-1][k]*transMatrix[k][s]*obs_model_probs[s][t];
					
					if(delta_new>=deltas[t][s]){
						deltas[t][s] = delta_new;
						alphas[t][s] = k;
					}
					
				}
				
			}
			
		}
		
	}
	
}
