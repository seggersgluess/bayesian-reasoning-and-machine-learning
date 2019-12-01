package ObjectDetection;

import java.io.IOException;
import java.util.HashMap;

import DataManagement.ParserPGM;

public class TestViolaJones {

	
	@SuppressWarnings("static-access")
	public static void calcAccuracyOfTrainedViolaJones() throws ClassNotFoundException, IOException {
		
		int nFaceImages = 20;
		int nNonFaceImages = 20;
		
		String directoryNameFaceImages = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/Training/face/"; 
		String directoryNameNonFaceImages = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/Training/non-face/";
		
		//String fileName4UploadTrainedObj = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/vjFullSampleTraining_IterNo_20.txt";
		String fileName4UploadTrainedObj = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/vjWithFullSampleTraining.txt";
		ViolaJones vj = new ViolaJones();
		vj = vj.loadViolaJonesObject(fileName4UploadTrainedObj);
		vj.resetImagesAndLabels();
		vj.useTrainedData4Classification();
		vj = set_face_and_non_face_images4Training(vj, nFaceImages, nNonFaceImages, directoryNameFaceImages, directoryNameNonFaceImages);
		
		HashMap<String, Double> accStats = vj.getAccuracyRatesOfViolaJones();
	    double acc = accStats.get("totalAccuracy");
	    System.out.println("Total Accuracy Rate: " + acc);
	    
	    double fp = accStats.get("falsePositive");
	    System.out.println("False Positive Rate: " + fp);
		
	}
	
	
	public static ViolaJones set_face_and_non_face_images4Training(ViolaJones vj, int numberOfFaceImages, int numberOfNonFaceImages, String directoryNameFaceImages, String directoryNameNonFaceImages) throws IOException {
		
		if(directoryNameFaceImages == null) {
			System.out.println("Directory to face images not set yet. Please set the path.");
			return null;
		}
		
		if(directoryNameNonFaceImages == null) {
			System.out.println("Directory to non-face images not set yet. Please set the path.");
			return null;
		}
		
		if(vj==null) {
			vj = new ViolaJones();
		}
		
		//set face images to vj-object
		long startTime = System.currentTimeMillis();
		
		for(int i=0; i<numberOfFaceImages; i++) {
			String imageNumber = String.format("%05d", (i+1));
			String fileName = directoryNameFaceImages + "face"+ imageNumber + ".pgm";
			//System.out.println("face"+ imageNumber + ".pgm");
			ParserPGM pgm = new ParserPGM(fileName, false, true);
			double [][] grayScaledImage = pgm.getScaledPGMImage();
			vj.setImageAndLabel(grayScaledImage , 1.0);
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("Loading " + numberOfFaceImages + " face data after: " + ((endTime-startTime)/1000.0) + " secs.");
			
		//set non-face images to vj-object
		startTime = System.currentTimeMillis();
		for(int i=0; i<numberOfNonFaceImages; i++) {
			String pgmName = "";
			if(i<559) {
				String imageNumber = String.format("%05d", (i+1));
				pgmName = "B1_" + imageNumber + ".pgm";
			}
			if(i>= 559 && i<899) { //899
				String imageNumber = String.format("%05d", (i-559));
				pgmName = "B5_" + imageNumber + ".pgm";
			}
			if(i>=899) {
				String imageNumber = String.format("%05d", (i+608));
				pgmName = "B20_" + imageNumber + ".pgm";
			}
				
			//System.out.println(pgmName);
			String fileName = directoryNameNonFaceImages + pgmName;
			ParserPGM pgm = new ParserPGM(fileName, false, true);
			double [][] grayScaledImage = pgm.getScaledPGMImage();			
			vj.setImageAndLabel(grayScaledImage , 0.0);
		}
		
		endTime = System.currentTimeMillis();
		System.out.println("Loading " + numberOfNonFaceImages + " non-face data after: " + ((endTime-startTime)/1000.0) + " secs.");
		
		return vj;
		
	}
	
	
	public static void main(String[] args) throws IOException, ClassNotFoundException {
    	
		calcAccuracyOfTrainedViolaJones();
	    
    }
	
}
