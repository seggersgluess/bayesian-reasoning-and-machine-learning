package ObjectDetection;

import java.util.ArrayList;
import java.util.List;

import Mathematics.MatrixOperations;

public class ViolaJonesCascadeClassifier {

	ArrayList<List<Double>> images = new ArrayList<List<Double>>();
	int imageWidth;
	int imageHeight;
	int nImages;
	List<Double> labels = new ArrayList<Double>();
	
	ArrayList<ViolaJones> vjClassifiers;
	int [] layers;
	
	
	public void cascadeClassifier() {
		
		int nLayers = layers.length;
		vjClassifiers = new ArrayList<ViolaJones>(nLayers);
		
		List<Integer> posIdxs = new ArrayList<Integer>();
		List<Integer> negIdxs = new ArrayList<Integer>();
		
		for(int i=0; i<nImages; i++) {
			double label = labels.get(i);
			if(label == 1) {
				posIdxs.add(i);
			}else {
				negIdxs.add(i);
			}
		}
		
		int n_posIdxs = posIdxs.size();

		for(int l=0; l<nLayers; l++) {
			
			List<Integer> falsePositives = new ArrayList<Integer>();
			int n_negIdxs = negIdxs.size();
			
			if(n_negIdxs == 0) {
				break;
			}
			
			ViolaJones vj = new ViolaJones();
			vj.set_imageHeigth(imageHeight);
			vj.set_imageWidth(imageWidth);
			vj.set_numberOfIterations(layers[l]);
			
			for(int i=0; i<n_posIdxs; i++) {
				int idx = posIdxs.get(i);
				vj.setImageAndLabel(images.get(idx),labels.get(idx));
			}
			
			for(int i=0; i<n_negIdxs; i++) {
				int idx = negIdxs.get(i);
				vj.setImageAndLabel(images.get(idx),labels.get(idx));
			}
				
			//TODO: features of images do not change-> not calculating for every layer!
			//Has to be calculated one time!
			//Then idx selection!
			vj.set_features4Images();
			vj.calc_ViolaJonesObjectDetection();
			
			vjClassifiers.add(vj);
			
			for(int i=0; i<n_negIdxs; i++) {
				falsePositives = new ArrayList<Integer>();
				int idx = negIdxs.get(i);
				List<Double> vecImage = images.get(idx);
				double [][] image = MatrixOperations.get_matrix_from_vec(vecImage, imageHeight, imageWidth);
				double classification = classify(image);
				if(classification == 1.0) {
					falsePositives.add(idx);
				}
			}
			
			negIdxs = falsePositives;					
		}		
	}
	
	
	public double classify(double [][] newImage) {
		
		int nClassifiers = vjClassifiers.size();
				
		for(int i=0; i<nClassifiers; i++) {
			double c = vjClassifiers.get(i).classify(newImage);
			if(c==0) {
				return 0.0;
			}
		}
		
		return 1.0;
		
	}
	
	
	public void set_defaultLayers() {
		layers = new int [5];
		layers[0] = 1;
		layers[1] = 5;
		layers[2] = 10;
		layers[3] = 50;
	}
	
	
}
