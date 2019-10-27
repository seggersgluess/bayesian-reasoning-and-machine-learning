package ObjectDetection;

import java.util.ArrayList;
import java.util.List;

import AdaptiveBasisModels.AdaBoost;
import Mathematics.MatrixOperations;

public class ViolaJones extends AdaBoost{

	ArrayList<List<Double>> images = new ArrayList<List<Double>>();
	int imageWidth;
	int imageHeight;
	int nImages;
	List<Double> labels = new ArrayList<Double>();
	AdaBoost adaBoost;
	
	boolean featuresSetted = false;
	
	@SuppressWarnings("static-access")
	public void calc_ViolaJonesObjectDetection() {
		
		if(featuresSetted == false) {
			System.out.println("Features not yet setted!");
			return;
		}
		
		nImages = images.size();
		adaBoost.explained_variable = new double [nImages][1];
		
		for(int i=0; i<nImages; i++) {
			adaBoost.explained_variable[i][0] = labels.get(i);
		}
		
		adaBoost.n_observations = nImages;
		adaBoost.nClasses = 2;
		
		//classes {1;0} if left then 1, else 0
		double [] inputClasses = new double [nClasses];
		inputClasses[0] = 1.0;
		
		adaBoost.set_classes4AdaBoost(inputClasses);
		
		adaBoost.doAdaBoost4ViolaJones();
		
	}
	
	
	@SuppressWarnings("static-access")
	public void set_features4Images() {
		
		adaBoost = new AdaBoost();
		
		List<Double> vecImage = images.get(0);
		
		double [][] image = MatrixOperations.get_matrix_from_vec(vecImage, imageHeight, imageWidth);
		
		RectangleFeatures4ViolaJones rectFeatures = new RectangleFeatures4ViolaJones(image);
		rectFeatures.create_rectangular_features();
		List<Double> featuresOfImage = rectFeatures.get_imageFeatures();
		int nFeatures = featuresOfImage.size();
		
		adaBoost.explaining_variables = new double [nImages][nFeatures];
		adaBoost.n_explaining_variables = nFeatures;
		adaBoost.names_of_explaining_variables = new String [nFeatures];
		
		for(int f=0; f<nFeatures; f++) {
			adaBoost.explaining_variables[0][f] = featuresOfImage.get(f);
			adaBoost.names_of_explaining_variables[f] = Integer.toString(f);
		}
		
		for(int i=1; i<nImages; i++) {
			
			rectFeatures = new RectangleFeatures4ViolaJones(image);
			rectFeatures.create_rectangular_features();
			featuresOfImage = rectFeatures.get_imageFeatures();
			
			for(int f=0; f<nFeatures; f++) {
				adaBoost.explaining_variables[i][f] = featuresOfImage.get(f);
			}
			
		}
		
		featuresSetted = true;
		
	}
	
	
	public void setImageAndLabel(List<Double> image, double label) {
		labels.add(label);
		images.add(image);
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
	
	
	@SuppressWarnings("static-access")
	public double classify(double [][] newImage) {
		
		double classification = 0.0;
		
		RectangleFeatures4ViolaJones rectFeatures = new RectangleFeatures4ViolaJones(newImage);
		rectFeatures.create_rectangular_features();
		List<Double> featuresOfImage = rectFeatures.get_imageFeatures();
		int nFeatures = featuresOfImage.size();
		
		double [][] newImageFeatures = new double [1][nFeatures];
		double summedAlphas = 0.0;
		
		for(int f=0; f<nFeatures; f++) {
			newImageFeatures[0][f] = featuresOfImage.get(f);
			summedAlphas = alphas.get(f);
		}
		
		summedAlphas *= 0.5;
		
		double alphaSummedWeakClassifiers = adaBoost.make_predictionWithDecisionStumpInfos(newImageFeatures);
		
		if(alphaSummedWeakClassifiers <= summedAlphas) {
			classification = 1.0;
		}
		
		return classification;
		
	}
	
	
	public void set_numberOfIterations(int n_iterations) {
		set_numberOfIterations4AdaBoost(n_iterations);
	}


	
}
