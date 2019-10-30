package ObjectDetection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import AdaptiveBasisModels.AdaBoost;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class ViolaJones extends AdaBoost{

	ArrayList<List<Double>> images = new ArrayList<List<Double>>();
	int imageWidth;
	int imageHeight;
	int nImages;
	List<Double> labels = new ArrayList<Double>();
	
	//rectangle features
	ArrayList<List<Double>> featuresOfImages;
	
	//weak classifier results
	ArrayList<HashMap<String,List<String>>> stump_infos;
	
	
	public void trainViolaJonesObjectDetection() {
		
		if(stump_infos == null) {
			featuresOfImages = new ArrayList<List<Double>>();			
			for(int i=0; i<nImages; i++) {
				
				List<Double> vecImage = images.get(0);
				
				double [][] image = MatrixOperations.get_matrix_from_vec(vecImage, imageHeight, imageWidth);
				
				RectangleFeatures4ViolaJones rectFeatures = new RectangleFeatures4ViolaJones(image);
				rectFeatures.create_rectangular_features();
				featuresOfImages.add(rectFeatures.get_imageFeatures());
			}
		}
				
		n_observations = nImages;
		//classes {1;0} if left then 1, else 0
		nClasses = 2;
		classes = new double [nClasses];
		classes[0] = 1.0;
		
		doAdaBoost4ViolaJones();
		
	}
		
	
	public void doAdaBoost4ViolaJones() {
		
		alphas = new ArrayList<Double>();
		
		List<Integer> rootIdxs = new ArrayList<Integer>(n_observations);
		List<Double> weights = new ArrayList<Double>(nImages);
		
		int nPos = 0;
		int nNeg = 0;
		
		int [] idx = Utilities.get_idx(labels,1.0);
		if(idx[0] == -1) {
			throw new RuntimeException("No positive examples found.");
		}else {
			nPos = idx.length;
		}
		
		idx = Utilities.get_idx(labels,0.0);
		if(idx[0] == -1) {
			throw new RuntimeException("No negative examples found.");
		}else {
			nNeg = idx.length;
		}
		
		for(int i=0; i<nImages; i++) {
			rootIdxs.add(i);
			if(explained_variable[i][0] == 1) {
				weights.add(1.0/(2.0*nPos));
			}else {
				weights.add(1.0/(2.0*nNeg));
			}			
		}
			
		if(stump_infos == null) {
			if(featuresOfImages == null) {
				throw new RuntimeException("No rectangle features are setted for running AdaBoost4Viola.");
			}
			stump_infos = get_decision_stump(featuresOfImages, true);
		}
		
		int nFeatures = stump_infos.size();
						
		for(int i=0; i<iterations; i++) {
				
			double minError = Double.MAX_VALUE;
			int minIdx = 0;
			int idx1 = 0;
			
			double summedWeights = 0.0;
			for(int j=0; j<nImages; j++) {
				summedWeights += weights.get(j);
			}
			
			for(int j=0; j<nImages; j++) {
				double label = labels.get(j);
				double weight = weights.get(j)/summedWeights;
				for(int f=0; f<nFeatures; f++) {
					
					double error = 1.0;
					List<String> leftIdxs = stump_infos.get(f).get("LeftIdxs");
					double y_pred = classes[0];
					boolean isLeft = leftIdxs.contains(Integer.toString(j));
					
					if(isLeft == false) {
						y_pred = classes[1];
					}		
					
					if(label != y_pred) {
						error = 1.0;
					}
					
					double weightedError = weight*error;
					if(weightedError<minError) {
						minError = weightedError;
						minIdx = idx1;
					}
				}
				idx1++;
			}
			
			classifierInfos.add(stump_infos.get((int) minIdx));
								
			double beta = minError/(1.0-minError);
			double alpha = Math.log(1.0/beta);
			alphas.add(alpha);	
			
			List<Double> updatedWeights = new ArrayList<Double>(nImages);
			
			for(int j=0; j<nImages; j++) {
				double e = minError;	
				double curWeight = weights.get(j);
				updatedWeights.add(curWeight*Math.pow(beta, (1.0-e)));				
			}
			
			weights = updatedWeights;
						
		}
			
	}
	
	
	public void set_pre_calculated_features(ArrayList<List<Double>> inputImages, int inputImageWidth, int inputImageHeight, ArrayList<HashMap<String,List<String>>> inputStump_infos, List<Double> inputLabels) {
		
		images = inputImages;
		nImages = images.size();
		imageWidth = inputImageWidth;
		imageHeight = inputImageHeight;
		stump_infos = inputStump_infos;
		labels = inputLabels;
		
	}
	
	
	public void setImageAndLabel(ArrayList<List<Double>> imageSelection, int inputImageWidth, int inputImageHeight, List<Double> labelSelection) {
		
		labels = labelSelection;
		images = imageSelection;
		nImages = labels.size();
		imageWidth = inputImageWidth;
		imageHeight = inputImageHeight;
		
	}
	
	
	public void set_imageWidth(int width) {
		imageWidth = width;
	}
	
	
	public void set_imageHeigth(int height) {
		imageHeight = height;
	}
	
	
	public void setImageAndLabel(double [][] image, double label) {
		
		int width = image[0].length;
		int height = image.length;
		
		if(imageWidth == 0) {
			imageWidth = width;
		}else {
			if(imageWidth != width) {
				throw new RuntimeException("Inconsitent image width detected. Check your input images.");
			}
		}
		
		if(imageHeight == 0) {
			imageHeight = height;
		}else {
			if(imageHeight != height) {
				throw new RuntimeException("Inconsitent image height detected. Check your input images.");
			}
		}
		
		labels.add(label);
		images.add(MatrixOperations.vecAsList(image));
		
	}
	
	
	public double classify(double [][] newImage) {
		
		double classification = 0.0;
		
		RectangleFeatures4ViolaJones rectFeatures = new RectangleFeatures4ViolaJones(newImage);
		rectFeatures.create_rectangular_features();
		List<Double> featuresOfImage = rectFeatures.get_imageFeatures();
		int nFeatures = featuresOfImage.size();
		
		double [][] newImageFeatures = new double [1][nFeatures];
		double summedAlphas = 0.0;
		names_of_explaining_variables = new String [nFeatures];
		
		for(int f=0; f<nFeatures; f++) {
			newImageFeatures[0][f] = featuresOfImage.get(f);
			summedAlphas = alphas.get(f);
			names_of_explaining_variables[f] = Integer.toString(f);
		}
		
		summedAlphas *= 0.5;
		
		double alphaSummedWeakClassifiers = make_predictionWithDecisionStumpInfos(newImageFeatures);
		
		if(alphaSummedWeakClassifiers <= summedAlphas) {
			classification = 1.0;
		}
		
		return classification;
		
	}
	
	
	public void set_numberOfIterations(int n_iterations) {
		set_numberOfIterations4AdaBoost(n_iterations);
	}

	
	public ArrayList<HashMap<String,List<String>>> get_weakClassifierResults() {
		return stump_infos;
	}

}
