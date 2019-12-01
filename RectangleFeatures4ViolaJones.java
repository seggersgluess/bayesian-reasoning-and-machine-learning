package ObjectDetection;

import java.util.ArrayList;
import java.util.List;


public class RectangleFeatures4ViolaJones {

	//nxn pixel (processing (n+1)x(n+1) pixel e.g. 24x24 pixel image as a 25x25 image)
	static double [][] image;
	static int imageWidth;
	static int imageHeight;
	
	static ArrayList<List<Double>> integratedImage;
    static List<Double> features;
	
	
	public RectangleFeatures4ViolaJones(double [][] image) {
		setImage4ViolaJones(image);
	}
	

	public static void create_rectangular_features() {
		
		calculate_integrated_images();
		
		features = new ArrayList<Double>();
		
		for(int w=1; w<imageWidth; w++) {		
			for(int h=1; h<imageHeight; h++) {	
				int x=0;				
				while((x+w)<imageWidth) {		
					int y=0;
					while((y+h)<imageHeight) {
												
						double baseRect = calculate_rectangle(x,y,w,h);
						double rightRect = 0.0;
						double rightRect2 = 0.0;
						double bottomRect = 0.0;
						double bottomRightRect = 0.0;
						double bottomRect2 = 0.0;
						double diff = 0;
						
						//Two rectangle features (horizontal)
						if((x+2*w)<imageWidth) {							
							rightRect = calculate_rectangle(x+w,y,w,h);			
							diff = baseRect - rightRect;						
							features.add(diff);
						}					
						
						//Two rectangle features (vertical)
						if((y+2*h)<imageHeight) {
							bottomRect = calculate_rectangle(x,y+h,w,h);
							diff = bottomRect-baseRect;
							features.add(diff);
						}
						
						//Three rectangle features (horizontal)
						if((x+3*w)<imageWidth) {
							rightRect2 = calculate_rectangle(x+2*w,y,w,h);
							diff = baseRect+rightRect2-rightRect;
							features.add(diff);						
						}
						
						//Three rectangle features (vertical)
						if((y+3*h)<imageHeight) {
							bottomRect2 = calculate_rectangle(x,y+2*h,w,h);
							diff = baseRect+bottomRect2-bottomRect;							
							features.add(diff);
						}
						
						//Four rectangle features
						if((x+2*w)<imageWidth && (y+2*h)<imageHeight) {
							bottomRightRect = calculate_rectangle(x+w,y+h,w,h);
							diff = baseRect+bottomRightRect-(rightRect+bottomRect);
							features.add(diff);
						}
						y++;			
					}
					x++;
				}							
			}		
		}
	}
	
	
	public static double calculate_rectangle(int x, int y, int w, int h) {
		
		if(w<0 && h<0) {
			throw new RuntimeException("Invalid width or height supplied for rectangle calculation.");
		}
		 		
		double rect1 = integratedImage.get(x).get(y);
		double rect2 = integratedImage.get(x+w).get(y);
		double rect3 = integratedImage.get(x).get(y+h);
		double rect4 = integratedImage.get(x+w).get(y+h);
			
        double rect = rect4+rect1-(rect2+rect3);
		
		return rect;
	
	}
	
	
	//cumulated pixels in (x,y) position
	public static void calculate_integrated_images() {
		
		if(image == null) {
			System.out.println("No image supplied yet for calculating integrated image.");
			return;
		}
		
		integratedImage = new ArrayList<List<Double>>();
				
		double prevCum = 0.0;
		
		for(int x=0; x<imageWidth; x++) {	
			double s = 0.0;
			List<Double> cumRow = new ArrayList<Double>(imageHeight);			
			for(int y=0; y<imageHeight; y++) {
	            s+= image[y][x];	
	            if(x!=0) {
	            	prevCum = integratedImage.get(x-1).get(y);
	            }
				cumRow.add(s+prevCum);
			}
			integratedImage.add(cumRow);
		}
		
	}
	
	
	//Attention: nxn images has to be (n+1)x(n+1) images with 0 in first row and column
	//Features are not correctly calculated if nxn frame is used.
	public static void setImage4ViolaJones(double [][] inputImage) {
		
		int nRows = inputImage.length;
		int nCols = inputImage[0].length;
		
		image = new double [(nRows+1)][(nCols+1)];
		
		//1st rows and cols are 0
		for(int i=0; i<nRows; i++) {
			for(int j=0; j<nCols; j++) {
				image[(i+1)][(j+1)] = inputImage[i][j];
			}
		}
		
		imageWidth = nCols+1;
		imageHeight = nRows+1;
	}
		
	
	public List<Double> get_imageFeatures(){
		return features;
	}
	
	
	@SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
		
		double [][] testImage = new double [3][3];

		int n = 1;
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				testImage[i][j] = n;
				n++;
			}
		}
		
		
		RectangleFeatures4ViolaJones rectFeat = new RectangleFeatures4ViolaJones(testImage);
		rectFeat.create_rectangular_features();
		
		System.out.println(rectFeat.calculate_rectangle(0, 2, 2, 1));
		
		int nFeatures = rectFeat.get_imageFeatures().size();
		System.out.println(nFeatures);
	}
		
}
