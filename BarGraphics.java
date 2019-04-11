package Graphics;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import Utilities.Utilities;

@SuppressWarnings("serial")
public class BarGraphics extends GraphicDevice{

	public static List<Color> barColor = new ArrayList<Color>();
	
	
	@Override
	protected void paintComponent(Graphics g) {
		
	    super.paintComponent(g);
	    	   
	    g2 = (Graphics2D) g;     
	    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	    
	    set_rectangular();
	     
	    if(useDifferentColors == true){	    	
	    	barColor = new ArrayList<Color>();	
	    	
	    	for(int i=0; i<y.size(); i++){
	
					int red   = getDefaultColors()[i][0];
					int green = getDefaultColors()[i][1];
					int blue  = getDefaultColors()[i][2];
					
					barColor.add(new Color(red,green,blue));							
	    	}
	    }
	    
	    createBarPlot();
	    
	    createGroupedBarPlot();
	    
	    grid = false;
	      
	    set_title2plot();
	    set_xLabel2plot();
	    set_yLabel2plot();
	        
	}
	
	
	public void createBarPlot(){
		
		int [] idxs = get_idxs_4_plot_type("B");
		
		if(idxs[0] == -1){			
			return;
		}
		
    	int nSamples = idxs.length; 
		
	    int rectangularWidth  = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
    	set_x_axis();
	    set_y_axis();
	       
    	int counter = 0;   	
    	
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(counter>=nSamples){	    		
	    		break;	    		
	    	}
    		
    		int yBottom = (border_gap+rectangularHeight)*(r+1);
    		
    	    for(int c=0; c<numberOfPlotColumns; c++){
  				
    	    	double y_max = 0.0;
    	    	int n = 1;
    	    	int [] plotIdxs = new int[0]; 
    	    	
    	    	int plotInfoNumber = plotInfo.get(idxs[counter]);
    	    	
    	    	if(counter != 0){
    	    		
    	    		if(plotInfoNumber != plotInfo.get(idxs[counter-1])){
    	    			y_max = get_y_max(plotInfoNumber);
    	    			plotIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
    	    			n = plotIdxs.length;
    	    		}
    	    		
    	    	}else{   	    		
    	    		y_max    = get_y_max(plotInfoNumber);
    	    		plotIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
    	    		n        = plotIdxs.length;
    	    	}
    	    	 		 	
    	    	List<Double> yValues = y.get(idxs[counter]);    	    	
    	    	int sampleLength     = yValues.size();
    	    	
    	    	double stepSize = (double) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/(n*sampleLength);    	    	   	    		    		
    	    	double scale    = (double)rectangularHeight/(double)y_max;
    	    			    		
    	    	int x0 = 0;
    	    	int y0 = 0;
    	    	int x1 = 0;
    	    	int y1 = 0;
				
    	    	int nIterations = 0;
    	    	
	    	    for(int i=0; i<sampleLength; i++){
	    		
	    	    	for(int p=0; p<n; p++){
	    	    		
	    	    		x0 = (int) (border_gap*(c+1)+rectangularWidth*c + nIterations*stepSize);
	    	    		x1 = (int) stepSize;
	    	    			
	    	    		y0 = (int)((border_gap*(r+1)+rectangularHeight*(r+1))-scale*(y.get(plotIdxs[p])).get(i));						
	    	    		y1 = (int)(scale*(y.get(plotIdxs[p])).get(i));  						
	    				
	    				g2.setColor(barColor.get(plotIdxs[p]));
	    				g2.fillRect(x0, y0, x1, y1);
	    					
	    				g2.setColor(rectColor);
	    				g2.drawRect(x0, y0, x1, y1);
	    				
	    				int xLeft   = border_gap*(c+1)+rectangularWidth*c;
	        	    	int xRight  = (border_gap+rectangularWidth)*(c+1);
	     			
	        	    	g2.setColor(Color.BLACK);
	        	    	g2.drawLine(xLeft, yBottom, xRight, yBottom);
	        	    	
	        	    	nIterations++;
	        	    	
	        	    	if(i==sampleLength-1){
	        	    		counter++;
	        	    	}
	        	    			
	    	    	}
	    	    		    	   	    	    	
				}			
    	    		    	    
	    	    if(counter>=nSamples){    	    		
    	    		break;    	    		
    	    	}
	    	    
    	    }
    		
    	}
		
	}
	
	
	public void createGroupedBarPlot(){
		
		int [] idxs = get_idxs_4_plot_type("BGroup");
		
		if(idxs[0] == -1){
			return;
		}
		
    	int nSamples = idxs.length; 
    	    	
	    int rectangularWidth  = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
    	set_x_axis();
	    set_y_axis();	    
	    
    	int counter = 0;   	
    	
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(counter>=nSamples){	    		
	    		break;	    		
	    	}
    		
    		int yBottom = (border_gap+rectangularHeight)*(r+1);
    		
    	    for(int c=0; c<numberOfPlotColumns; c++){
  				
    	    	double y_max = 0.0;
    	    	int n = 1;
    	    	int [] plotIdxs = new int[0]; 
    	    	
    	    	int plotInfoNumber = plotInfo.get(idxs[counter]);
    	    	
    	    	if(counter != 0){
    	    		
    	    		if(plotInfoNumber != plotInfo.get(idxs[counter-1])){
    	    			
    	    			plotIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
    	    			n = plotIdxs.length;
    	    			
    	    			for(int i=0; i<n; i++){
    	    				y_max = y_max + Utilities.getMax(y.get(plotIdxs[i]));
    	    			}  	    			
    	    		}
    	    		
    	    	}else{   	    		
    	    		
    	    		plotIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
    	    		n        = plotIdxs.length;
    	    		
    	    		for(int i=0; i<n; i++){
	    				y_max = y_max + Utilities.getMax(y.get(plotIdxs[i]));
	    			}
    	    	}
    	    	 		 	
    	    	List<Double> yValues = y.get(idxs[plotIdxs[0]]);    	    	
    	    	int sampleLength     = yValues.size();
    	    	
    	    	double stepSize = (double) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/sampleLength;    	    	   	    		    		
    	    	double scale    = (double)rectangularHeight/(double)y_max;
    	    			    		
    	    	int x0 = 0;
    	    	int y0 = 0;
    	    	int x1 = 0;
    	    	int y1 = 0;
				
    	    	//int nIterations = 0;
    	    	
	    	    for(int i=0; i<sampleLength; i++){
	    		
	    	    	for(int p=0; p<n; p++){
	    	    		
	    	    		x0 = (int) (border_gap*(c+1)+rectangularWidth*c + i*stepSize);
	    	    		x1 = (int) stepSize;
	    	    		
	    	    		if(p==0){
	    	    			y0 = (int)((border_gap*(r+1)+rectangularHeight*(r+1))-scale*(y.get(plotIdxs[p])).get(i));								    	    		 
	    	    		}else{
	    	    			y0 = (int)(y0-scale*(y.get(plotIdxs[p])).get(i));
	    	    		}
	    	    		 	
	    	    		y1 = (int)(scale*(y.get(plotIdxs[p])).get(i));
	    				
	    				g2.setColor(barColor.get(plotIdxs[p]));
	    				g2.fillRect(x0, y0, x1, y1);
	    					
	    				g2.setColor(rectColor);
	    				g2.drawRect(x0, y0, x1, y1);
	    				
	    				int xLeft   = border_gap*(c+1)+rectangularWidth*c;
	        	    	int xRight  = (border_gap+rectangularWidth)*(c+1);
	     			
	        	    	g2.setColor(Color.BLACK);
	        	    	g2.drawLine(xLeft, yBottom, xRight, yBottom);
	        	    	
	        	    	//nIterations++;
	        	    	
	        	    	if(i==sampleLength-1){
	        	    		counter++;
	        	    	}
	        	    			
	    	    	}
	    	    		    	   	    	    	
				}			
    	    		    	    
	    	    if(counter>=nSamples){    	    		
    	    		break;    	    		
    	    	}
	    	    
    	    }
    		
    	}
		
	}

	
	public static void plotBars(String [][] xValues, double [][] yValues, boolean newPlot){
		
		setDefaultDesign();
		
		convert_input_data(xValues, yValues, "B");
		
	    int nSamples = xValues[0].length;
	    
		int nPlotInfos = plotInfo.size();
		
		for(int i=0; i<nSamples; i++){
			
			if(nPlotInfos == 0){
				plotInfo.add(0);
			}else{
				if(newPlot == true){				
					plotInfo.add(plotInfo.get(nPlotInfos-1)+1);				
				}else{			
					plotInfo.add(plotInfo.get(nPlotInfos-1));			
				}
			}
			
			nPlotInfos++;
			
			barColor.add(new Color(44, 102, 230, 180));
			
		}
		
	}
	
	
	public static void plotGroupedBars(String [][] xValues, double [][] yValues, boolean newPlot){
		
		setDefaultDesign();
		
		convert_input_data(xValues, yValues, "BGroup");
		
	    int nSamples = xValues[0].length;
	    
		int nPlotInfos = plotInfo.size();
		
		for(int i=0; i<nSamples; i++){
			
			if(nPlotInfos == 0){
				plotInfo.add(0);
			}else{
				if(newPlot == true){				
					plotInfo.add(plotInfo.get(nPlotInfos-1)+1);				
				}else{			
					plotInfo.add(plotInfo.get(nPlotInfos-1));			
				}
			}
			
			nPlotInfos++;
			
			barColor.add(new Color(44, 102, 230, 180));
			
		}
		
	}
	
	
	public void set_x_axis(){
        
		Font orgFont = g2.getFont();
		g2.setFont(fontOfXAxisUnits);
		
		int nPlots              = Utilities.getMaxFromIntList(plotInfo)+1;		 
		int [] sampleNumbers    = new int [nPlots];
		
		for(int i=0; i<nPlots; i++){			
			sampleNumbers[i] = Utilities.get_idx(plotInfo, i)[0];			
		}
				
		int counter = 0;
			
		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
	    int r = 0;
	       
		while(r<numberOfPlotRows){
			
			if(counter>=nPlots){    	    		
	    		break;    	    		
	    	}
						
			int yUp     = border_gap*(r+1)+rectangularHeight*r;
			int yBottom = (border_gap+rectangularHeight)*(r+1);
			
			int c = 0;
			
			while(c<numberOfPlotColumns){
				
				int xLeft   = border_gap*(c+1)+rectangularWidth*(c);
				int xRight  = (border_gap+rectangularWidth)*(c+1);
				
				g2.setColor(Color.BLACK);
				g2.drawLine(xLeft, yBottom, xRight, yBottom);

				int sampleLength   = y.get(sampleNumbers[counter]).size(); 
				               
				double barWidth =  (double)((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/sampleLength; ;
								
			    double width      = 0.05*border_gap;
				   
			    int pos = 0;
			    
			    while(pos<sampleLength){
				    
			    	int x0 = (int) (border_gap*(c+1)+rectangularWidth*c + pos*barWidth+barWidth/2.0);			    	
			    	int x1 = x0;
			    	int y0 = yBottom;
			    	int y1 = y0;
			    	
			    	g2.setColor(Color.BLACK);
			    	g2.drawLine(x0, y0, x1, y0+(int)width);
				       
			    	String xlabel = String.format(xStr.get(sampleNumbers[counter]).get(pos))+ "";
			    	FontMetrics metrics = g2.getFontMetrics();
			    	int labelWidth = metrics.stringWidth(xlabel);
			    	
			    	g2.setColor(textColor);
			    	g2.drawString(xlabel, x0-labelWidth/2, y0+metrics.getHeight() + 3);
				       
			    	if(grid == true){
				                 
			    		g2.setColor(grid_color);
			    		g2.drawLine(x0, y1, x0, yUp); 
				    		   	       	       	   
			    	}
				    
			    	g2.setColor(Color.BLACK);
			    	
			    	pos = pos+numberXDivisions;
			    	
			    }
				
				counter++;
				
				if(counter>=nPlots){    	    		
		    		break;    	    		
		    	}
				
				c++;
								
			}
			
			r++;
							
		}
		 
		g2.setFont(orgFont);
	    g2.setColor(Color.BLACK);
	    
	}
	
	
	public void set_y_axis(){

		Font orgFont = g2.getFont();
	    g2.setFont(fontOfYAxisUnits);
		
	    int nPlots = Utilities.getMaxFromIntList(plotInfo)+1;
	    int [] plotInfoNumbers = Utilities.intGenerator(nPlots);
	    
		int counter = 0;
			
		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
	    int r = 0;
	    
		while(r<numberOfPlotRows){
			
			if(counter>=nPlots){	    		
	    		break;	    		
	    	}
				
			int yUp     = border_gap*(r+1)+rectangularHeight*r;
			int yBottom = (border_gap+rectangularHeight)*(r+1);
			
			int c = 0;
			
			while(c<numberOfPlotColumns){
				
				double y_max = 0.0;
				
				int plotInfoNumber = plotInfoNumbers[counter];
				
				int xLeft   = border_gap*(c+1)+rectangularWidth*(c);
				int xRight  = (border_gap+rectangularWidth)*(c+1);
				
			    g2.setColor(Color.BLACK);
			    g2.drawLine(xLeft-1, yBottom, xLeft-1, yUp);
				
			    int step_size     = rectangularHeight/(numberYDivisions); 
			    double width      = 0.05*border_gap;
				   
			    double y_min = 0.0;
			    
			    if(plotType.get(counter).equals("B")==true){
			    	y_max = get_y_max(plotInfoNumber);
			    }
				
			    if(plotType.get(counter).equals("BGroup")==true){
			    	
			    	int [] plotIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
    	    		int n           = plotIdxs.length;
    	    		
    	    		for(int i=0; i<n; i++){
	    				y_max = y_max + Utilities.getMax(y.get(plotIdxs[i]));
	    			}
			    		    	
			    }
				   
			    double l = (y_max-y_min)/(numberYDivisions);
				   
			    for(int i = 0; i<numberYDivisions+1; i++) {
				       
			    	int x0 = xLeft;
			    	int x1 = x0-(int) width;
			    	int y0 = (step_size*i)+yUp;
			    	int y1 = y0;
			    	
			    	g2.setColor(Color.BLACK);
			    	g2.drawLine(x0, y0, x1, y1);
				       			    	 
			    	String ylabel = String.format(yAxisUnitsFormat, ((int) y_max-l*i))+ "";
			    	FontMetrics metrics = g2.getFontMetrics();
			        int labelWidth = metrics.stringWidth(ylabel);
			        
			        g2.setColor(textColor);
			        g2.drawString(ylabel, x0 - labelWidth - 5, y0 + (metrics.getHeight() / 2) - 3);
				       
			        if(grid == true){
				    	   
			           	if(i != numberYDivisions){
				    		   
			           		g2.setColor(grid_color);
			           		g2.drawLine(xLeft, y1, xRight, y1); 
				    		   
			           	}
				    	      
			        }
				          
			    }
				
				counter++;
				
				if(counter>=nPlots){    	    		
		    		break;    	    		
		    	}
			    	
				c++;	
			}
			
			r++;
								
		}
		  
		g2.setFont(orgFont);
		g2.setColor(Color.BLACK);
		
	}
	
	
    public static void setDefaultDesign(){
    	
	    setNumberOfXDivisions(5);
	    setNumberOfYDivisions(5);
    	
	    setRectColor(new Color(200, 200, 200, 200));
	    setGridColor(new Color(238,238,238));
	    setTextColor(new Color(128,128,128));
	    
	    setFontOfXAxisUnits("bold",12);
	    setColorOfXAxisUnits(new Color(128,128,128));
	    setFontOfYAxisUnits("bold",12);
	    setColorOfYAxisUnits(new Color(128,128,128));
    	
    }
	
	
    private static void createAndShowGui() {

    	BarGraphics mainPanel = new BarGraphics();

    	String frameLabel;
    	
    	frameLabel = "Bar Plot";    		
    	
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
	
   	
   	//Example 1: Bar plots
   	public static void parametrizationExample1(){
   		
		int maxDataPoints = 20;
	      
	 	String [][] x_values = new String [maxDataPoints][2];
	 	double [][] y_values = new double [maxDataPoints][2];
	 	      
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 	   x_values[i][0] = Integer.toString(i);
	 	   y_values[i][0] = Math.pow(i, 2.0);
	 	   x_values[i][1] = Integer.toString(i);
	 	   y_values[i][1] = i+60;
	 	}
   		
	 	String [] title    = {"Bar Plot 1", "Bar Plot 2"};
	 	String [] subTitle = {"SubTitle1", "SubTitle2"};
	 	String [] yLabel   = {"Bar Plot 1", "Bar Plot 2"};
	 	String [] xLabel   = {"x1", "x2"};
	 	
	 	//List<Color> lineColor = new ArrayList<Color>();
	 	//lineColor.add(Color.RED);
	 	//lineColor.add(Color.BLACK);
	 	
	 	setNumberOfPlotColums(2);
	 	setNumberOfPlotRows(1);
	 	
	 	setGraphWidth(900);
	 	setGraphHeight(500);
	 	
	 	plotBars(x_values,y_values, true);
	 	plotBars(x_values,y_values, false); 
	 	
 	 	setTitle(title, null, "12");
	 	setSubTitle1(subTitle, null, "10");
	 	//setSubTitle2("SubTitle2", "bold", "10");
	 	setYLabel(yLabel, null, "12");
	 	setXLabel(xLabel, null, "10"); 
	 	setNumberOfDigits4YAxis(2);
	 	setFontOfXAxisUnits("plain", 10);
	 	setFontOfYAxisUnits("plain", 10);
	 	
	 	plot();
	 	
   	}
   	
   	
   	//Example 2: Bar plots
   	public static void parametrizationExample2(){
   		
		int maxDataPoints = 20;
	      
	 	String [][] x_values1 = new String [maxDataPoints][1];
	 	double [][] y_values1 = new double [maxDataPoints][1];
	 	  
	 	String [][] x_values2 = new String [maxDataPoints][1];
	 	double [][] y_values2 = new double [maxDataPoints][1]; 
	 	
	 	String [][] x_values3 = new String [maxDataPoints][1];
	 	double [][] y_values3 = new double [maxDataPoints][1]; 
	 	
	 	String [][] x_values4 = new String [maxDataPoints][1];
	 	double [][] y_values4 = new double [maxDataPoints][1]; 
	 	
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 	   x_values1[i][0] = Integer.toString(i);
	 	   x_values2[i][0] = "Data " + Integer.toString(i+1);
	 	   x_values3[i][0] = "Data " + Integer.toString(i+1);
	 	   x_values4[i][0] = "Data " + Integer.toString(i+1);
	 	   y_values1[i][0] = Math.pow(i, 2.0);
	 	   y_values2[i][0] = i+1;
	 	   y_values3[i][0] = i+5;
	 	   y_values4[i][0] = i+16;
	 	}
	 	   
	 	String [] title    = {"Bar Plot 1", "Bar Plot 2", "Bar Plot 3"};
	 	String [] subTitle = {"SubTitle1", "SubTitle2"};
	 	String [] yLabel   = {"Bar Plot 1", "Bar Plot 2", "Bar Plot 3"};
	 	String [] xLabel   = {"x1", "x2", "x3"};
	 	
	 	//List<Color> lineColor = new ArrayList<Color>();
	 	//lineColor.add(Color.RED);
	 	//lineColor.add(Color.BLACK);
	 	
	 	setNumberOfPlotColums(3);
	 	setNumberOfPlotRows(1);
	 	
	 	setGraphWidth(1200);
	 	setGraphHeight(500);
	 	
	 	plotBars(x_values1,y_values1, true);
	 	plotBars(x_values2,y_values2, true); 
	 	plotBars(x_values3,y_values3, false); 
	 	plotBars(x_values4,y_values4, true);
	 	
 	 	setTitle(title, null, "12");
	 	setSubTitle1(subTitle, null, "10");
	 	//setSubTitle2("SubTitle2", "bold", "10");
	 	setYLabel(yLabel, null, "12");
	 	setXLabel(xLabel, null, "10");
	 	setNumberOfDigits4YAxis(0);   
	 	setFontOfXAxisUnits("plain", 11);
	 	setFontOfYAxisUnits("plain", 11);
	 	
	 	plot();
	 	
   	}


    //Example 3: Grouped bar plots
	public static void parametrizationExample3(){
		
		int maxDataPoints = 18;
      
		String [][] x_values1 = new String [maxDataPoints][1];
		double [][] y_values1 = new double [maxDataPoints][1];
 	  
		String [][] x_values2 = new String [maxDataPoints][1];
		double [][] y_values2 = new double [maxDataPoints][1]; 
 	
		String [][] x_values3 = new String [maxDataPoints][1];
		double [][] y_values3 = new double [maxDataPoints][1]; 
 	
		String [][] x_values4 = new String [maxDataPoints][1];
		double [][] y_values4 = new double [maxDataPoints][1]; 
 	
 		for(int i = 0; i < maxDataPoints ; i++) {
 		   	x_values1[i][0] = Integer.toString(i);
 	   		x_values2[i][0] = "Data " + Integer.toString(i+1);
 	   		x_values3[i][0] = "Data " + Integer.toString(i+1);
 	   		x_values4[i][0] = "Data " + Integer.toString(i+1);
 	   		y_values1[i][0] = i+6; //Math.pow(i, 2.0);
 	   		y_values2[i][0] = i+1;
 	   		y_values3[i][0] = i+5;
 	   		y_values4[i][0] = i+16;	   
 		}
 	   
 		int maxDataPoints2 = 55;
 		
 		String [][] x_values5 = new String [maxDataPoints2][1];
		double [][] y_values5 = new double [maxDataPoints2][1];
 		
 		for(int i = 0; i < maxDataPoints2 ; i++) {
 		   	x_values5[i][0] = Integer.toString(i);
 	   		y_values5[i][0] = i+16;
   
 		}
		
 		String [] title    = {"Bar Plot 1", "Bar Plot 2", "Bar Plot 3", "Bar Plot 4"};
 		String [] subTitle = {"SubTitle1", "SubTitle2"};
 		String [] yLabel   = {"Bar Plot 1", "Bar Plot 2", "Bar Plot 3", "Bar Plot 4"};
 		String [] xLabel   = {"x1", "x2", "x3"};
 	
 		//List<Color> lineColor = new ArrayList<Color>();
 		//lineColor.add(Color.RED);
 		//lineColor.add(Color.BLACK);
 	 		
 		setNumberOfPlotColums(3);
 		setNumberOfPlotRows(2);
 	
 		setGraphWidth(1200);
 		setGraphHeight(500);
 	
 		plotGroupedBars(x_values2,y_values2, true); 
 		plotGroupedBars(x_values3,y_values3, false);
 		
 		plotGroupedBars(x_values1,y_values1, true); 
 		plotGroupedBars(x_values2,y_values2, false);
 		
 		plotGroupedBars(x_values1,y_values1, true); 
 		plotGroupedBars(x_values3,y_values3, false); 		
 		
 		plotGroupedBars(x_values5,y_values5, true);
 		
 		setTitle(title, null, "12");
 		setSubTitle1(subTitle, null, "10");
 		//setSubTitle2("SubTitle2", "bold", "10");
 		setYLabel(yLabel, null, "12");
 		setXLabel(xLabel, null, "10");
 		setNumberOfDigits4YAxis(0);   
 		setFontOfXAxisUnits("plain", 11);
 		setFontOfYAxisUnits("plain", 11);
 	
 		plot();
 	
	}

	
	public static void main(String[] args) {
		
		parametrizationExample3();
				
	}
	
	
}
