package Images;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;


@SuppressWarnings("serial")
public class Image extends JPanel {

	public static int pref_w = 400;
	public static int pref_h = 400;
	public static Color [][] image;
	public static int imageWidth;
	public static int imageHeight;
	
	public static Graphics2D g2;
	
	public static int pixelWidth = 0;
	public static int pixelHeight = 0;
	
	@Override
	protected void paintComponent(Graphics g) {
		
		super.paintComponent(g);
        g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                
        createImage();
         
	}
	
	
	public static void plotImage(Color [][] image4plot) {
		image = image4plot;
		imageWidth = image[0].length;
		imageHeight = image.length;
		
		pixelWidth = (int) pref_w/imageWidth;
		pixelHeight = (int) pref_h/imageHeight;
		
		pref_w = pixelWidth*imageWidth;
		pref_h = pixelHeight*imageHeight;
		
		plot();
	}
		
	
	public static void createImage() {
			
		int xPos = 0;
		int yPos = 0;
		
		for(int i=0; i<imageWidth; i++) {
			yPos = 0;
			for(int j=0; j<imageHeight; j++) {
				Color pixelRGBcolor = image[j][i];
				g2.setColor(pixelRGBcolor);				
		        g2.fillRect(xPos, yPos, pixelWidth, pixelHeight);
		        yPos += pixelHeight;
			}
			xPos += pixelWidth;
		}
		
		System.out.println(xPos);
		
	}
	
	
	public static void setPlotWidth(int plotWidth) {
		pref_w = plotWidth;
	}
	
	
	public static void setPlotHeight(int plotHeight) {
		pref_h = plotHeight;
	}
	
	
	private static void createAndShowGui() {

    	Image mainPanel = new Image();

    	String frameLabel;
    	
    	frameLabel = "Graph";    		

    	JFrame frame = new JFrame(frameLabel);
    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	frame.getContentPane().add(mainPanel);
    	frame.pack();
    	frame.setLocationByPlatform(true);
    	frame.setVisible(true);
      
    }
    
    
   	public static void plot(){
 	   
   		SwingUtilities.invokeLater(new Runnable() {
   			
   			public void run() {  				
   				createAndShowGui();
   			}
   			
   		});
	   
   	}
	
   	
    @Override
    public Dimension getPreferredSize() {	   
    	return new Dimension(pref_w, pref_h);
    }
   	
   	
}
