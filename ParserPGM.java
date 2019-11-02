package DataManagement;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;

import Images.Image;


public class ParserPGM {

	String imageName;
	String file;
	
	double [][] image;
	int imageWidth;
	int imageHeight;
	int maxValue;
	String magicNumber; //PGM format type (P2 ASCII || P5 Binary)
	
	double [][] scaledImage;
	Color [][] rgbImage;
	
	public ParserPGM(String fileName, boolean calcGrayScaledImage, boolean calcRGBImage) throws IOException {
		
		file = fileName;
		 
		BufferedReader br = new BufferedReader(new FileReader(file));
		
		//read in header infos
		int headerLength = 0;
		String line = br.readLine();
		magicNumber = line;
		headerLength += line.length();
		
		line = br.readLine();
		int idx = line.indexOf(" ");
		imageHeight = Integer.parseInt(line.substring(0, idx));
		imageWidth = Integer.parseInt(line.substring(idx+1, line.length()));
		headerLength += line.length();
		
		line = br.readLine();
		maxValue = Integer.parseInt(line);
		headerLength += line.length();
		headerLength += 2; //two blank separators for lines
		
		br.close();
		
		int nPixels = imageWidth*imageHeight;
		
		image = new double[imageHeight][imageWidth];
		
		if(calcGrayScaledImage == true || calcRGBImage == true) {
			scaledImage = new double [imageHeight][imageWidth];
		}
		
		if(calcRGBImage == true) {
			rgbImage = new Color[imageHeight][imageWidth];
		}
		
		InputStream fileInputStream = new FileInputStream(file);
		DataInputStream dis = new DataInputStream(fileInputStream);
		
		int counter = 0;
		int c; 
				
		while(counter<headerLength) {
			c = dis.readUnsignedByte();
			counter++;
		}
		counter = 0;
		while(counter<nPixels) {
			
		   //TODO: Implemented only for P5 --> open P2!
		   c = dis.readUnsignedByte(); 
		   if(c<0) {
			   break;
		   }		   
		   if (c == '\n') { 
			   // do nothing if new line encountered			   
           }else {

        	   byte b = (byte)(c & 0xFF);
        	   if(b<0) {
        		   image[counter/imageWidth][counter % imageWidth] = (128-Math.abs(b))+128;
        	   }else {
        		   image[counter/imageWidth][counter % imageWidth] = b;
        	   }
        	   if(calcGrayScaledImage == true) {
        		   scaledImage[counter/imageWidth][counter % imageWidth] = image[counter/imageWidth][counter % imageWidth]/maxValue;
        	   }
        	   if(calcRGBImage == true) {
        		   if(calcGrayScaledImage == false) {
        			   scaledImage[counter/imageWidth][counter % imageWidth] = image[counter/imageWidth][counter % imageWidth]/maxValue;       			   
        		   }
        		   int rgbVal = (int) ((scaledImage[counter/imageWidth][counter % imageWidth])*255.0);
    			   rgbImage[counter/imageWidth][counter % imageWidth] = new Color(rgbVal, rgbVal, rgbVal);
        	   }
               counter++;
           }		   
		}
		
		fileInputStream.close();
		dis.close();
		
	}
	
	
	public double [][] getPGMImage() {
		return image;
	}
	
	
	public double [][] getScaledPGMImage() {
		
		if(scaledImage == null) {
			System.out.println("No scaled image calculated.");
			return null;
		}
		
		return scaledImage;
		
	}
	
	
	public Color [][] getPGMImageWithRGBColors() {
		
		if(rgbImage == null) {
			System.out.println("No rgb image calculated.");
			return null;
		}
		
		return rgbImage;
	}
	
	
	public String getFileName() {
		return file;
	}
	
	
	public int getImageHeight() {
		return imageHeight;
	}
	
	
	public int getImageWidth() {
		return imageWidth;
	}
	
	
	public int getMaxValue() {
		return maxValue;
	}
	
	
	public String getPGMType() {
		return magicNumber;
	}
	
	
	@SuppressWarnings("static-access")
	public void plotImage() {
		
		if(rgbImage == null) {
			System.out.println("No RGB information for image plot created.");
			return;
		}
		
		Image obj_image = new Image();
		
		obj_image.plotImage(rgbImage);
		//obj_image.plot();
	}
	
	
	@SuppressWarnings("static-access")
	public void plotImage(int width, int height) {
		
		if(rgbImage == null) {
			System.out.println("No RGB information for image plot created.");
			return;
		}
		
		Image obj_image = new Image();
		
		obj_image.setPlotWidth(width);
		obj_image.setPlotHeight(height);
		obj_image.plotImage(rgbImage);
	}
	
	
	public static void main(String[] args) throws Exception {
    	
    	String fileName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/Boosting/ViolaJones/Test/face/cmu_0226.pgm";
    	
    	ParserPGM pgm = new ParserPGM(fileName, true, true);   	
    	pgm.plotImage(40,40);
	    
    }
	
}
