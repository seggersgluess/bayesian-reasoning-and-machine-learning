package ObjectDetection;

import java.awt.Color;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import AdaptiveBasisModels.AdaBoost;
import Images.Image;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

@SuppressWarnings("serial")
public class ViolaJones extends AdaBoost implements Serializable {

	ArrayList<List<Double>> images = new ArrayList<List<Double>>();
	int imageWidth;
	int imageHeight;
	int nImages;
	List<Double> labels = new ArrayList<Double>();
	
	public static ArrayList<HashMap<String,List<String>>> classifierInfos;
	
	//properties necessary for strorage of objects
	List<Double> trainedAlphas;
	ArrayList<HashMap<String,List<String>>> trainedClassifierInfos;
	int nUsedClasses;
	double [] usedClasses;
	
	//rectangle features
	ArrayList<List<Double>> featuresOfImages;
	
	boolean featurePreSelection;
	double percentage4PreSelection = 0.1;
	List<Integer> selectedIdxs;
	
	
	@SuppressWarnings("static-access")
	public void trainViolaJonesObjectDetection(boolean logStats, boolean storeTrainedObject, boolean featurePreSelect) throws IOException {
		
		long startTime = System.currentTimeMillis();
		featurePreSelection = featurePreSelect;
		
		featuresOfImages = new ArrayList<List<Double>>();	
						
		for(int i=0; i<nImages; i++) {
					
			List<Double> vecImage = images.get(i);
				
			double [][] image = MatrixOperations.get_matrix_from_vec(vecImage, imageHeight, imageWidth);
				
			RectangleFeatures4ViolaJones rectFeatures = new RectangleFeatures4ViolaJones(image);
								
			rectFeatures.create_rectangular_features();
			if(selectedIdxs == null) {
				featuresOfImages.add(rectFeatures.get_imageFeatures());	
			}else {
				int nRedFeatures = selectedIdxs.size();
				List<Double> recFeature = rectFeatures.get_imageFeatures();
				List<Double> redFeaturesOfImage = new ArrayList<Double>();
				for(int f=0; f<nRedFeatures; f++) {
					redFeaturesOfImage.add(recFeature.get(selectedIdxs.get(f)));
				}
				featuresOfImages.add(redFeaturesOfImage);				
			}				
		}
			
		if(selectedIdxs != null) {
			featurePreSelection = true;
		}
		
		if(featurePreSelection == true && selectedIdxs == null) {
			System.out.println("Do preselection of rectangle features.");
			int nRectangleFeatures = featuresOfImages.get(0).size();
			int nBestFeatures = (int) Math.round(percentage4PreSelection*nRectangleFeatures);
			FeaturePreSelector preSelector = new FeaturePreSelector(labels, featuresOfImages, nBestFeatures);
			selectedIdxs = preSelector.get_selectedBest_f_scores_Idxs();
			
			ArrayList<List<Double>> reducedFeaturesOfImages = new ArrayList<List<Double>>();
			for(int s=0; s<nImages; s++) {
				List<Double> redFeaturesOfImage = new ArrayList<Double>();
				for(int i=0; i<nBestFeatures; i++) {
					redFeaturesOfImage.add(featuresOfImages.get(s).get(selectedIdxs.get(i)));
				}
				reducedFeaturesOfImages.add(redFeaturesOfImage);
			}
			featuresOfImages = reducedFeaturesOfImages;
		}
				
		n_observations = nImages;
		//classes {1;0} if left then 1, else 0 --> 1 = face/ 0 = non-face
		nClasses = 2;
		classes = new double [nClasses];
		classes[0] = 1.0;
				
		doAdaBoost4ViolaJones(logStats, storeTrainedObject);
		
		//Set trained classification infos to Viola-Jones-Object
		setTrainedParameters2Object();
		
		long endTime = System.currentTimeMillis();
		System.out.println("Finished training Viola Jones after " + ((endTime-startTime)/1000.0) + " secs.");
	}
	
	
	public void doAdaBoost4ViolaJones(boolean logStats, boolean storeTrainedObjects) throws IOException {
		
		alphas = new ArrayList<Double>();

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
			double label = labels.get(i);
			if(label == 1) {
				weights.add(1.0/(2.0*nPos));
			}else {
				weights.add(1.0/(2.0*nNeg));
			}			
		}
			
		classifierInfos = new ArrayList<HashMap<String,List<String>>>();
		
		int counter = 1;
		String dir4Save = dirManager4ViolaJones().get("dirName4Iterations");
		
		for(int i=0; i<iterations; i++) {
				
			System.out.println("Start iteration " + (i+1));
			
		    double sumOfWeights = 0.0;
		    for(int j=0; j<nImages; j++) {
		    	sumOfWeights += weights.get(j);
		    }
		    
			HashMap<String,List<String>> minWeakClassifier = get_decision_stump(weights, sumOfWeights);
			
			classifierInfos.add(minWeakClassifier);
					
			double minError = Double.parseDouble(minWeakClassifier.get("MinError").get(0));
					
			System.out.println("Number of Errors: " + minWeakClassifier.get("NumberOfErrors").get(0));
			
			double beta = minError/(1.0-minError);
			double alpha = 0.0;
			if(beta == 0.0) {
				alpha = Double.MAX_VALUE;
			}else {
				alpha = Math.log(1.0/beta);
			}
			alphas.add(alpha);	
			
			if(logStats == true) {				
				System.out.println("Classification Statistics for iteration "+(i+1));
				HashMap<String, Double> accStats = getAccuracyRatesOfViolaJones();
			    double acc = accStats.get("totalAccuracy");
			    System.out.println("Total Accuracy Rate: " + acc);
			    
			    double fp = accStats.get("falsePositive");
			    System.out.println("False Positive Rate: " + fp);
			}
						
			if(beta == 0.0) {
				System.out.println("Found perfect weak classifier. Break after " + (i+1) + " iterations.");
				iterations = i+1;
				break;
			}
			
			List<Double> updatedWeights = new ArrayList<Double>(nImages);
						
			double e = 0.0;
			
			for(int j=0; j<nImages; j++) {				

				boolean isCorrectClass = checkCorrectClassification(minWeakClassifier, j);
				if(isCorrectClass == false) {
					e = 1.0;
				}else {
					e = 0.0;
				}				
				
				double curWeight = weights.get(j)/sumOfWeights;
				updatedWeights.add(curWeight*Math.pow(beta, (1.0-e)));				
			}
			
			weights = updatedWeights;
			
			if(storeTrainedObjects == true) {
				if(counter == 5 || i == 0) {
					setTrainedParameters2Object();
			        System.out.println("Store trained vj object after " + (i+1) + " iterations.");
					String directory = dir4Save + (i+1) + ".txt";
					saveViolaJonesObject(this,directory);
					if(counter == 5) {
						counter = 0;
					}					
				}
				
				counter++;
			}
			
		}
			
	}

	
	public void setTrainedParameters2Object() {
		trainedAlphas = alphas;
		trainedClassifierInfos = classifierInfos;
		nUsedClasses = nClasses;
		usedClasses = classes;
	}
	
	
	public HashMap<String,List<String>> get_decision_stump(List<Double> weights, double sumOfWeights){
		
		HashMap<String,List<String>> splittingInfos = new HashMap<String,List<String>>();
		
		int nSamples = featuresOfImages.size();
		int nExpFeatures = featuresOfImages.get(0).size();
		
		double minError = Double.MAX_VALUE;
		int optSplitFeature = 0;
    	double optThreshold = 0.0;
    	double optPolarity  = 0.0;
    	int numberOfErrors  = 0;
			
		for(int j=0; j<nExpFeatures; j++) {
			
			//System.out.println(j);

			for(int k=0; k<nSamples; k++) {
				
				double t = featuresOfImages.get(k).get(j);
								
				double error_pos = 0.0;
				double error_neg = 0.0;
				double error = 0.0;
				double polarity  = 1.0;
				int nPosErrors = 0;
				int nNegErrors = 0;
				int nErrors = 0;
				
	            for(int s=0; s<nSamples; s++) {	
	            	
	            	double y_pred_pos = 0.0;
	            	double y_pred_neg = 0.0;
	            	
	            	double feature = featuresOfImages.get(s).get(j);
	            	//polarity = 1.0
	            	if(feature<=t) {
	            		y_pred_pos = 1.0;
	            		
	            	}else {
	            		y_pred_pos = 0.0;
	            		
	            	}
	            	
	            	//polarity = -1.0
	            	if((-1.0*feature)<=(-1.0*t)) {
	            		y_pred_neg = 1.0;
	            	}else {
	            		y_pred_neg = 0.0;
	            	}
	            	
	            	double label = labels.get(s);
	            	
	            	if(y_pred_pos != label) {
	            		error_pos += weights.get(s)/sumOfWeights;
	            		nPosErrors++;
	            	}
	            	if(y_pred_neg != label) {
	            		error_neg += weights.get(s)/sumOfWeights;
	            		nNegErrors++;
	            	}
	            }
	           	  
	            if(error_pos < error_neg) {
	            	error = error_pos;
	            	nErrors = nPosErrors;
	            }else {
	            	error = error_neg;
	            	nErrors = nNegErrors;
	            	polarity = -1.0;
	            }
	            
	            if(error<minError) {	            	
	            	minError = error;
	            	optSplitFeature = j;
	            	optThreshold = t;
	            	optPolarity = polarity;
	                numberOfErrors = nErrors;
	            }
   	               
			}
					
		}
		
        List<String> splitFeat  = new ArrayList<String>(1);
        List<String> threshold  = new ArrayList<String>(1);
        List<String> pol        = new ArrayList<String>(1);
        List<String> splitError = new ArrayList<String>(1);
        List<String> nErr       = new ArrayList<String>(1);
        splitFeat.add(Integer.toString(optSplitFeature));			        
        threshold.add(Double.toString(optThreshold));
        pol.add(Double.toString(optPolarity));
        splitError.add(Double.toString(minError));
        nErr.add(Integer.toString(numberOfErrors));    
        
        splittingInfos.put("SplittingFeature", splitFeat);   
        splittingInfos.put("Threshold", threshold);
    	splittingInfos.put("Polarity", pol);
    	splittingInfos.put("MinError", splitError);
    	splittingInfos.put("NumberOfErrors", nErr);
		
		return splittingInfos;
		
	}	
	
	
	public boolean checkCorrectClassification(HashMap<String, List<String>> splittingInfos, int imageIdx) {
		
		boolean isCorrect = true;
		
		int splittingFeature = Integer.parseInt(splittingInfos.get("SplittingFeature").get(0));
		double t = Double.parseDouble(splittingInfos.get("Threshold").get(0));
		double p = Double.parseDouble(splittingInfos.get("Polarity").get(0));
		
		double label = labels.get(imageIdx);
		double rectFeature = featuresOfImages.get(imageIdx).get(splittingFeature);
		
		double y_pred = classes[0];
		
		if((p*rectFeature)>(p*t)) {
			y_pred = classes[1];
		}
			
		if(label != y_pred) {
			isCorrect = false;
		}
		
		return isCorrect;
		
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
		nImages++;
		
	}
	
	
	public void useTrainedData4Classification() {
		
		if(trainedAlphas == null) {
			System.out.println("Viola-Jones not trained yet. No trained parameters available.");
			return;
		}
		
		alphas = trainedAlphas;
		classifierInfos = trainedClassifierInfos;
		classes = usedClasses;
		nClasses = nUsedClasses;
		
	}
	
	
	@SuppressWarnings("static-access")
	public double classify(double [][] newImage) {
		
		double classification = 0.0;
		
		RectangleFeatures4ViolaJones rectFeatures = new RectangleFeatures4ViolaJones(newImage);
		rectFeatures.create_rectangular_features();
		List<Double> featuresOfImage = rectFeatures.get_imageFeatures();
		
		if(featurePreSelection == true) {
			List<Double> redFeaturesOfImage = new ArrayList<Double>();
			int nRedFeatures = selectedIdxs.size();
			for(int i=0; i<nRedFeatures; i++) {
				int idx = selectedIdxs.get(i);
				redFeaturesOfImage.add(featuresOfImage.get(idx));
			}
			featuresOfImage = redFeaturesOfImage;
		}
		
		int nFeatures = featuresOfImage.size();
		int nSelFeatures = alphas.size();
		
		double [][] newImageFeatures = new double [1][nFeatures];
		double summedAlphas = 0.0;
		names_of_explaining_variables = new String [nFeatures];
		
		for(int f=0; f<nFeatures; f++) {
			newImageFeatures[0][f] = featuresOfImage.get(f);			
			names_of_explaining_variables[f] = Integer.toString(f);
		}
		
		for(int f=0; f<nSelFeatures; f++) {
			summedAlphas += alphas.get(f);
		}

		summedAlphas *= 0.5;
		
		double alphaSummedWeakClassifiers = make_predictionFromViolaJonesAdaBoost(newImageFeatures);
		
		if(alphaSummedWeakClassifiers >= summedAlphas) {
			classification = 1.0;
		}
		
		return classification;
		
	}
	
	
	//Make prediction for y with new_x (has to be a 1xn_explaining_vars vector)
	public static double make_predictionFromViolaJonesAdaBoost(double [][] new_x) {
		
		if(classifierInfos == null) {
			System.out.println("AdaBoost not trained yet. Can´t make prediction from new input.");
		}
		
		double prediction = 0.0;
		int nIterations = alphas.size();
		
		for(int i=0; i<nIterations; i++) {
			
			String splitFeat = get_splittingFeature4Iteration(i);
			double threshold = get_threshold4Iteration(i);
			double polarity  = get_polarity4Iteration(i);
			double p = 1.0;

			if(polarity != 0.0) {
				p = polarity;
			}
			
			int idx = Utilities.get_idx(names_of_explaining_variables, splitFeat)[0];
			double inputValue4Feature = new_x[0][idx];
			
			inputValue4Feature *=p;
			threshold *=p;
			
			if(inputValue4Feature <= threshold) {
				prediction += classes[0]*alphas.get(i);
			}else {
				prediction += classes[1]*alphas.get(i);
			}
			
		}
		
		return prediction;
		
	}
	
	
	public static String get_splittingFeature4Iteration(int iteration) {
		String feature = classifierInfos.get(iteration).get("SplittingFeature").get(0);
		return feature;
	}
	
	
	public static double get_threshold4Iteration(int iteration) {
		double threshold = Double.parseDouble(classifierInfos.get(iteration).get("Threshold").get(0));
		return threshold;
	}
	
	
	public static double get_polarity4Iteration(int iteration) {
		double polarity = Double.parseDouble(classifierInfos.get(iteration).get("Polarity").get(0));
		return polarity;
	}
	
	
	public HashMap<String, Double> getAccuracyRatesOfViolaJones() {
		
		double truePositive = 0.0;
		double falsePositive = 0.0;
		double trueNegative = 0.0;
		double falseNegative = 0.0;
		double totalCorrect = 0.0;
		double posCount = 0.0;
		double negCount = 0.0;
		
		for(int i=0; i<nImages; i++) {
	    	double [][] image = getImageFromSample(i);
	    	double label = labels.get(i);
		    double res = classify(image);	
		    
		    if(label == 1.0) {
		    	posCount++;
		    	if(res == 1.0) {
		    		truePositive++;
		    		totalCorrect++;
		    	}else {
		    		falsePositive++;
		    	}
		    }
		    
		    if(label == 0) {
		    	negCount++;
		    	if(res == 1.0) {
		    		falseNegative++;
		    	}else {
		    		trueNegative++;
		    		totalCorrect++;
		    	}
		    }		    
	    }
			
		HashMap<String,Double> accuracyStats = new HashMap<String,Double>();
		accuracyStats.put("truePositive", truePositive/posCount);  //-> also named detection rate!
		accuracyStats.put("falsePositive", falsePositive/negCount);
		accuracyStats.put("trueNegative", trueNegative/negCount);
		accuracyStats.put("falseNegative", falseNegative/posCount);
		accuracyStats.put("totalAccuracy", totalCorrect/nImages);
		
		return accuracyStats;
		
	}
	
	
	public void resetImagesAndLabels() {
		images = new ArrayList<List<Double>>();
		labels = new ArrayList<Double>();
		nImages = 0;
	}
	
	
	public void set_numberOfIterations(int n_iterations) {
		set_numberOfIterations4AdaBoost(n_iterations);
	}

	
	public double [][] getImageFromSample(int imageNumber) {		
		if(imageNumber > (images.size()-1)) {
			System.out.println(imageNumber + " invalid image number.");
		}
		double [][] scaledImage = MatrixOperations.get_matrix_from_vec(images.get(imageNumber), imageHeight, imageWidth);
		return scaledImage;
	}
	
	
	@SuppressWarnings("static-access")
	public void showPGMImage(int imageNumber) {
		if(imageNumber > (images.size()-1)) {
			System.out.println(imageNumber + " invalid image number.");
		}
		
		Image im = new Image();
		
		double [][] scaledImage = getImageFromSample(imageNumber);
		Color [][] image4Plot = new Color [imageHeight][imageWidth];
		
		for(int i=0; i<imageWidth; i++) {
			for(int j=0; j<imageHeight; j++) {
				int rgbVal = (int) (scaledImage[j][i]*255);
				image4Plot[j][i] = new Color(rgbVal, rgbVal, rgbVal);
			}
		}
		
		im.setPlotWidth(imageWidth+10);
		im.setPlotHeight(imageHeight+10);
		im.plotImage(image4Plot);
			
	}
	
	
	public static HashMap<String,String> dirManager4ViolaJones() {
		HashMap<String,String> dirNames = new HashMap<String,String>();
		dirNames.put("dirName4Iterations", "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/vjFullSampleTraining_IterNo_");
		
		return dirNames;
	}
	
	
	public static void saveViolaJonesObject(ViolaJones vj, String fileName) throws IOException{
	    FileOutputStream fos = new FileOutputStream(fileName);
	    ObjectOutputStream oos = new ObjectOutputStream(fos);
	    oos.writeObject(vj);
	    oos.close();
	}

	
	public static ViolaJones loadViolaJonesObject(String fileName) throws IOException, ClassNotFoundException{
	   FileInputStream fin = new FileInputStream(fileName);
	   ObjectInputStream ois = new ObjectInputStream(fin);
	   ViolaJones iHandler= (ViolaJones) ois.readObject();
	   ois.close();
	   return iHandler;
	}
	
}
