package ObjectDetection;

import java.io.IOException;

import DataManagement.ParserPGM;

public class TrainViolaJones {

	static ViolaJones vj;
	static String directoryNameFaceImages;
	static String directoryNameNonFaceImages;
	static int numberOfFaceImages = 20; //2429;
	static int numberOfNonFaceImages = 20; //3698;
	
	public static void set_face_and_non_face_images4Training() throws IOException {
		
		if(directoryNameFaceImages == null) {
			System.out.println("Directory to face images not set yet. Please set the path.");
			return;
		}
		
		if(directoryNameNonFaceImages == null) {
			System.out.println("Directory to non-face images not set yet. Please set the path.");
			return;
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
		
	}
	
	
	public static void setPath2FaceImages(String directory) {
		directoryNameFaceImages = directory;
	}
	
	
	public static void setPath2NonFaceImages(String directory) {
		directoryNameNonFaceImages = directory;
	}
	
	
	@SuppressWarnings("static-access")
	public static ViolaJones loadViolaJonesObject(String fileName) throws IOException {
		
		ViolaJones vj = new ViolaJones();		
		try {
			vj = vj.loadViolaJonesObject(fileName);
			return vj;
		} catch (ClassNotFoundException e) {
			System.out.println("Can´t load Viola Jones object.");
			return null;
		}
	}
	
	
	@SuppressWarnings("static-access")
	public static void doTwoStepTraining() throws IOException {
		
		TrainViolaJones trainVJ = new TrainViolaJones();
		
		String dirFaceImages = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/Training/face/";
		String dirNonFaceImages = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/Training/non-face/";
			
		//Step 1: Make feature preselection based on F-scores to reduce number of features (only 10% of the original number)
		numberOfFaceImages    = 400;
		numberOfNonFaceImages = 400;
		
		trainVJ.setPath2FaceImages(dirFaceImages);
	    trainVJ.setPath2NonFaceImages(dirNonFaceImages);
	    trainVJ.set_face_and_non_face_images4Training();
	    
	    vj.set_numberOfIterations(0);
	    vj.trainViolaJonesObjectDetection(false, false, true);
	    
	    String fileName2Save = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/vjWithImagesAndPreSelection.txt";
	    vj.saveViolaJonesObject(vj, fileName2Save);
	    
	    //Step 2: Train the VJ with preselected features on the whole sample (6000 images)
	    vj.resetImagesAndLabels();
	    numberOfFaceImages = 1500;
		numberOfNonFaceImages = 1700;
	    trainVJ.set_face_and_non_face_images4Training();
	    
	    vj.set_numberOfIterations(50);
	    vj.trainViolaJonesObjectDetection(true, true, false);
	    
	    fileName2Save = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/vjWithFullSampleTraining.txt";
	    //Delete images from trained object (for reducing size of stored object)
	    vj.resetImagesAndLabels();
	    vj.saveViolaJonesObject(vj, fileName2Save);
	    
	    System.out.println("Finished training Viola Jones with " + (numberOfFaceImages+numberOfNonFaceImages) + " images.");
	    	    
	}
	
	
	public static void main(String[] args) throws IOException {
    	
		doTwoStepTraining();
	    
    }
	
}
