package ObjectDetection;

import java.util.ArrayList;
import java.util.List;


public class RectangleFeatures4ViolaJones {

	double [][] image;
	int imageWidth;
	int imageHeight;
	
	ArrayList<List<Double>> integratedImage;
    List<Double> features;
	
	
	public RectangleFeatures4ViolaJones(double [][] image) {
		this.image = image;
		imageWidth = image[0].length;
		imageHeight = image.length;
	}
	

	public void create_rectangular_features() {
		
		calculate_integrated_images();
		
		features = new ArrayList<Double>();
		
		for(int w=1; w<imageWidth; w++) {
			int x=0;	
			for(int h=1; h<imageHeight; h++) {				
				while((x+w)<imageWidth) {
					int y=0;					
					while((y+h)<imageHeight) {
						
						//Two rectangle features (horizontal)
						double baseRect = calculate_rectangle(x,y,w,h);
						double rightRect = 0.0;
						double rightRect2 = 0.0;
						double bottomRect = 0.0;
						double bottomRightRect = 0.0;
						double diff = 0;
						
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
						
						//Three rectangle features
						if((x+2*w)<imageWidth) {
							rightRect2 = calculate_rectangle(x+2*w,y,w,h);
							diff = baseRect+rightRect2-rightRect;
							features.add(diff);
						}
						
						//Four rectangle features
						if((x+2*w)<imageWidth && (y+2*h)<imageHeight) {
							bottomRightRect = calculate_rectangle(x+2*w,y+2*h,w,h);
							diff = baseRect+bottomRightRect-(rightRect+bottomRect);
						}
						y++;			
					}
					x++;
				}							
			}		
		}
	}
	
	
	public double calculate_rectangle(int x, int y, int w, int h) {
		
		if(w<=0 && h<=0) {
			throw new RuntimeException("Invalid witdth or height supplied for rectangle calculation.");
		}
		
		double rect1 = integratedImage.get(x).get(y);
		double rect2 = integratedImage.get(x+w).get(y);
		double rect3 = integratedImage.get(x).get(y+h);
		double rect4 = integratedImage.get(x+w).get(y+h);
		
		double rect = rect4+rect1-(rect2+rect3);
	
		return rect;
	
	}
	
	
	public void calculate_integrated_images() {
		
		if(image == null) {
			System.out.println("No image supplied yet for calculating integrated image.");
			return;
		}
		
		integratedImage = new ArrayList<List<Double>>();
		
		for(int x=0; x<imageWidth; x++) {	
			List<Double> cumRow = new ArrayList<Double>(imageHeight);
			double s = 0.0;
			for(int y=0; y<imageHeight; y++) {
				s+= image[x][y];
				cumRow.add(s);
			}
			integratedImage.add(cumRow);
		}
		
	}
	
	
	public void setImage4ViolaJones(double [][] inputImage) {
		image = inputImage;
		imageWidth = image[0].length;
		imageHeight = image.length;
	}
		
	
	public List<Double> get_imageFeatures(){
		return features;
	}
	
}
