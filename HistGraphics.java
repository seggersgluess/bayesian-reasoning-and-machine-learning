package Graphics;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import Distributions.NormalDistribution;
//import Distributions.NormalDistribution;
import Mathematics.GeneralMath;
import Utilities.Utilities;

@SuppressWarnings("serial")
public class HistGraphics extends GraphicDevice{

    public static int nBins;
    public static ArrayList<List<Double>> breakpoints = new ArrayList<List<Double>>();
    public static ArrayList<List<Integer>> counts  = new ArrayList<List<Integer>>();
    public static List<Double> maxValue = new ArrayList<Double>();
    public static List<Double> minValue = new ArrayList<Double>();
	public static boolean freq = true;
	
	public static Color barColor = new Color(135,206,250,180);
	public static Color lineColor = new Color(0,0,128,180); //Color.WHITE;

	public static boolean plotNormalPDF = false;
	public static List<Integer> plotIdxs4NormalPDF = new ArrayList<Integer>();
	
	
	@Override
	protected void paintComponent(Graphics g) {
		
		if(maxValue.size() != minValue.size()){
			throw new RuntimeException("Number of max. x-values equals not number of min. x-values.");
		}
		
		for(int i=0; i<maxValue.size(); i++){
			if(maxValue.get(i) <= minValue.get(i)){
				throw new RuntimeException("Max. x-value lower then min. x-value. Correct your input.");	
			}
		}
		
	    super.paintComponent(g);
	    	   
	    g2 = (Graphics2D) g;     
	    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	    
	    set_rectangular();
	        
	    createHistogram();
	    
	    if(plotNormalPDF == true){
	    	
	    	if(freq == false){
	    		create_normalPDF();   		
	    	}else{
	    		System.out.println("Cannot plot normal density for frequencies.");
	    	}
	    	
	    }
	    
	    set_title2plot();
	    set_xLabel2plot();
	    set_yLabel2plot();
	        
	}
	
	
    public void createHistogram(){
    	
    	int nSamples = y.size();
    	
	    set_breakpoints();
	    count_elements_between_breakpoints();
	    	
	    int rectangularWidth  = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
    	set_x_axis();
    	set_y_axis();
	    
    	int idx = 0;
    	
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(idx>=nSamples){	    		
	    		break;	    		
	    	}
    		
    		int yBottom = (border_gap+rectangularHeight)*(r+1);
    		
    	    for(int c=0; c<numberOfPlotColumns; c++){
  				
    	    	double stepSize = (double) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/nBins;
    	    	int maxFreq     = Utilities.getMaxFromIntList(counts.get(idx));
    	    			    		
    	    	double scale = (double)rectangularHeight/(double)maxFreq;
    	    			    		
    	    	int x0 = 0;
    	    	int y0 = 0;
    	    	int x1 = 0;
    	    	int y1 = 0;
    	    		
    	    	for(int i=0; i<nBins; i++){

    	    		x0 = (int) (border_gap*(c+1)+rectangularWidth*c + i*stepSize);
    	    		x1 = (int) stepSize;
    	    			
    	    		y0 = (int)((border_gap*(r+1)+rectangularHeight*(r+1))-scale*(counts.get(idx)).get(i));						
    	    		y1 = (int)(scale*(counts.get(idx)).get(i));
    							
    				if(useDifferentColors == true){
    					
    					int red   = getDefaultColors()[idx][0];
    					int green = getDefaultColors()[idx][1];
    					int blue  = getDefaultColors()[idx][2];
    					
    					barColor = new Color(red,green,blue);
    					
    				}
    				
    				g2.setColor(barColor);
    				g2.fillRect(x0, y0, x1, y1);
    					
    				g2.setColor(lineColor);
    				g2.drawRect(x0, y0, x1, y1);
    				
    				int xLeft   = border_gap*(c+1)+rectangularWidth*c;
        	    	int xRight  = (border_gap+rectangularWidth)*(c+1);
     			
        	    	g2.setColor(grid_color);
        	    	g2.drawLine(xLeft, yBottom, xRight, yBottom);
        	    	g2.setColor(Color.BLACK);
        	    	
    	    	}
    	    	   
    	    	idx++;
    	    	
    	    	if(idx>=nSamples){    	    		
    	    		break;    	    		
    	    	}
    	    	
    	    }
    		
    	}
   
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
    
    
    public static void freqHist(boolean freqFlag){
    	
    	//if true then frequencies else relative counts (which sum to one).
    	freq = freqFlag;
    	
    }
    
    
    public static void set_breakpoints(){
    	
    	if(nBins == 0){
    		throw new RuntimeException("Number of bins not set.");
    	}
    	
    	int nSamples = y.size();
    	
    	for(int i=0; i<nSamples; i++){
    		
        	double length = maxValue.get(i)-minValue.get(i);
        	double stepSize = length/nBins;
        	
        	List<Double> newBreaks = new ArrayList<Double>();
        		
        	for(int j=0; j<nBins+1; j++){
        		
        		newBreaks.add(minValue.get(i)+j*stepSize);
        		
        	}
    		
        	breakpoints.add(newBreaks);
        	
    	}
     	
    }
    
    
    public static void count_elements_between_breakpoints(){
    	
    	int nSamples = y.size();
    	
    	for(int i=0; i<nSamples; i++){
    		
        	List<Double> sample = y.get(i);
        	int sampleLength    = sample.size();
        	
        	List<Integer> newCounts = new ArrayList<Integer>(nBins);
 
        	for(int j=0; j<nBins; j++){
        		
        		int counter = 0;
        		
        		for(int k=0; k<sampleLength; k++){
        			
        			double curVal = sample.get(k);
        			
        			if(curVal>=(breakpoints.get(i)).get(j) && curVal<(breakpoints.get(i)).get(j+1)){
        				
        				counter++;
        				
        			}

        		}
        		
        		newCounts.add(counter);
        		
        	}
    		
    		counts.add(newCounts);
    		
    	}
    	  		
    }
    
    
    public static void set_numberOfBins(int numberOfBins){
    	
    	if(numberOfBins<= 0){
    		throw new RuntimeException("Invalid number of bins supplied.");
    	}
    	
    	nBins = numberOfBins;
    	
    }
    

    public static void plotHistogram(double [][] values, boolean newPlot){
  	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	    
    	maxValue.add(Utilities.getMax(y.get(y.size()-1)));
    	minValue.add(Utilities.getMin(y.get(y.size()-1)));
    	
    	int nPlotInfos = plotInfo.size();
    	
    	if(nPlotInfos == 0){
			plotInfo.add(0);
		}else{
			if(newPlot == true){				
				plotInfo.add(plotInfo.get(nPlotInfos-1)+1);				
			}else{			
				plotInfo.add(plotInfo.get(nPlotInfos-1));			
			}
		}
    	
    }
	
    
    public static void plotHistogram(double [][] values, boolean newPlot, boolean plotNormalPDF){
   	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	    
    	maxValue.add(Utilities.getMax(y.get(y.size()-1)));
    	minValue.add(Utilities.getMin(y.get(y.size()-1)));
    	
    	int nPlotInfos = plotInfo.size();
    	
    	if(nPlotInfos == 0){
			plotInfo.add(0);
		}else{
			if(newPlot == true){				
				plotInfo.add(plotInfo.get(nPlotInfos-1)+1);				
			}else{			
				plotInfo.add(plotInfo.get(nPlotInfos-1));			
			}
		}
    	
    	if(plotNormalPDF == true){
    		
    		plotIdxs4NormalPDF.add(plotInfo.get(plotInfo.size()-1));
    		
    		plotNormalPDF(true);
    		
    	}
    	
    }
    
    
    public static void plotHistogram(double [][] values, double max_value, double min_value, boolean newPlot){
   	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	   
    	maxValue.add(max_value);
    	minValue.add(min_value);
    	
    	int nPlotInfos = plotInfo.size();
    	
    	if(nPlotInfos == 0){
			plotInfo.add(0);
		}else{
			if(newPlot == true){				
				plotInfo.add(plotInfo.get(nPlotInfos-1)+1);				
			}else{			
				plotInfo.add(plotInfo.get(nPlotInfos-1));			
			}
		}
    	  	
    }
    
    
    public static void plotHistogram(double [][] values, double max_value, double min_value, boolean newPlot, boolean plotNormalPDF){
    	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	   
    	maxValue.add(max_value);
    	minValue.add(min_value);
    	
    	int nPlotInfos = plotInfo.size();
    	
    	if(nPlotInfos == 0){
			plotInfo.add(0);
		}else{
			if(newPlot == true){				
				plotInfo.add(plotInfo.get(nPlotInfos-1)+1);				
			}else{			
				plotInfo.add(plotInfo.get(nPlotInfos-1));			
			}
		}
    	
    	if(plotNormalPDF == true){
    		
    		plotIdxs4NormalPDF.add(plotInfo.get(plotInfo.size()-1));
    		
    		plotNormalPDF(true);
    		
    	}
    	
    }
    
    
    private static void createAndShowGui() {

    	HistGraphics mainPanel = new HistGraphics();

    	String frameLabel;
    	
    	frameLabel = "Histogram";    		
    	
    	JFrame frame = new JFrame(frameLabel);
    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	frame.getContentPane().add(mainPanel);
    	frame.pack();
    	frame.setLocationByPlatform(true);
    	frame.setVisible(true);
      
    }
    
    
    public static void noLinesAroundBars(){
    	
    	lineColor = barColor;
    	
    }
    
    
   	public static void plot(){
 	   
   		SwingUtilities.invokeLater(new Runnable() {
   			
   			public void run() {
   				
   				createAndShowGui();
	       	
   			}
   			
   		});
	   
   	}
   	
   	
   	public static void set_min_x_value(double [] min){
   	
   		if(min.length != y.size()){
   			throw new RuntimeException("Number of supplied min. x-values not equal to number of samples.");
   		}
   		
   		minValue = new ArrayList<Double>();
   		
   		for(int i=0; i<min.length; i++){
   			
   			minValue.add(min[i]);
   			
   		}
   			
   	}
   	
   	
   	public static void set_max_x_value(double [] max){
   		
   		if(max.length != y.size()){
   			throw new RuntimeException("Number of supplied max. x-values not equal to number of samples.");
   		}
   		
   		maxValue = new ArrayList<Double>();
   		
   		for(int i=0; i<max.length; i++){
   			
   			maxValue.add(max[i]);
   			
   		}
   		
   	}
   	
	
   	public void set_x_axis(){
   		
   		//g2.setColor(Color.BLACK); 
	    //g2.drawLine(border_gap, getHeight()-border_gap, getWidth()-border_gap, getHeight()-border_gap);

   		int nSamples = y.size();
   		
	    Font orgFont = g2.getFont();
	    g2.setFont(fontOfXAxisUnits);
	    
	    double step_size      = (double) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/numberXDivisions;
	    int rectangularWidth  = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
	    int idx = 0;
	    
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(idx>=nSamples){	    		
	    		break;	    		
	    	}
	    
    	    for(int c=0; c<numberOfPlotColumns; c++){
    	    	
    		    double width      = 0.0;//0.01*(getHeight()-2*border_gap);
    			double xLength    = (maxValue.get(idx)-minValue.get(idx))/numberXDivisions; 
    		    			
    		    int y0 = border_gap*(r+1)+rectangularHeight*(r+1);
    	    	int y1 = y0 - (int) width;
    		    
    		    for(int i = 0; i<numberXDivisions+1; i++) {
    			       
    		    	g2.setColor(Color.BLACK);
    		    	int x0 = (int) (step_size*(i)+border_gap*(c+1)+rectangularWidth*c);
    		    	int x1 = x0;
    		    	
    		    	//g2.drawLine(x0, y0, x1, y1);
    			    
    		    	if(grid == true){
    	                
    		    		Stroke oldStroke = g2.getStroke();
    		        	
    		        	g2.setStroke(new BasicStroke(4));	
    		    		g2.setColor(grid_color);
    		    		g2.drawLine(x0, y1, x1, (border_gap*(r+1)+rectangularHeight*r)); 
    		    		
    		    		if(i != numberXDivisions){
    		    			
    		    			g2.setStroke(new BasicStroke(2));
    			    		
    			    		int x0_half = (int) (x0+0.5*step_size);
    			    		int x1_half = x0_half;
    			    		g2.drawLine(x0_half, y1, x1_half, (border_gap*(r+1)+rectangularHeight*r)); 
    		    			
    		    		}
    		    		
    		    		g2.setStroke(oldStroke);
    			    	
    		    	}
    		    	
    		    	g2.setColor(colorOfXAxisUnits );
    		    	String xlabel = String.format(xAxisUnitsFormat, minValue.get(idx)+xLength*(i))+ "";
    			    FontMetrics metrics = g2.getFontMetrics();
    			    int labelWidth = metrics.stringWidth(xlabel);
    			    g2.drawString(xlabel, x0 - labelWidth/2, y0 + metrics.getHeight() + 3);  	    	
    		    	
    		    }
    	    	
    		    idx++;
    		    
    		    if(idx>=nSamples){    	    		
    	    		break;    	    		
    	    	}
    		    
    	    }
    		
    		
    	}	
    		
	    g2.setFont(orgFont);
	    g2.setColor(Color.BLACK);
	    
   	}
   	
   	
   	public void set_y_axis(){
   		
	    //g2.setColor(Color.BLACK);
	    //g2.drawLine(border_gap, getHeight()-border_gap, border_gap, border_gap);
		
   		int nSamples = y.size();
   		
    	Font orgFont = g2.getFont();
	    g2.setFont(fontOfYAxisUnits);
   		
	    int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
	    int idx = 0;
	    
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(idx>=nSamples){	    		
	    		break;	    		
	    	}
    		
    		for(int c=0; c<numberOfPlotColumns; c++){
    	    	
    		    int sampleLength = GeneralMath.sum(counts.get(idx));
    		    double maxFreq   = Utilities.getMaxFromIntList(counts.get(idx));	    		    
    			double scale     = (double) ((rectangularHeight)/maxFreq);
    		    int yMaxPos      = (int) ((double)(border_gap*(r+1)+rectangularHeight*(r+1))-scale*maxFreq);
    		    int yMinPos      = (int) ((double)yMaxPos-scale*maxFreq);
    			
    		    double step_size  = (yMaxPos-yMinPos)/(numberYDivisions); 
    		    double width      = 0.01*(rectangularWidth);
    			 
    		    double y_min = 0.0;
    		    double y_max = maxFreq;
    			   
    		    double l = (y_max-y_min)/(numberYDivisions);
    	        int label;
    		    String ylabel;
    	        
    		    for(int i = 0; i<numberYDivisions+1; i++){
    			       
    		    	g2.setColor(Color.BLACK);
    		    	int x0 = border_gap*(c+1)+rectangularWidth*c;
    		    	int x1 = x0+(int) width;
    		    	
    		    	if(i==numberYDivisions){
    		    		x1 = x0; //Scaling problem (only integers allowed): 0.00 not at origin of x/y axis!
    		    	}
    		    	
    		    	int y0 = (int) ((step_size*i)+yMaxPos);
    		    	int y1 = y0;
    		    	//g2.drawLine(x0, y0, x1, y1);
    			    
    		        if(grid == true){
    			    	
    		        	Stroke oldStroke = g2.getStroke();
    		        	
    		        	g2.setStroke(new BasicStroke(4));
    		           	g2.setColor(grid_color);
    		           	g2.drawLine(x1, y1, (border_gap+rectangularWidth)*(c+1), y1); 
    			    	
    		    		if(i != numberYDivisions){
    		    			
    		    			g2.setStroke(new BasicStroke(2));
    			    		
    			    		int y0_half = (int) (y0+0.5*step_size);
    			    		int y1_half = y0_half;
    			    		g2.drawLine(x1, y1_half, (border_gap+rectangularWidth)*(c+1), y1_half); 
    		    			
    		    		}
    		           	
    		           	g2.setStroke(oldStroke);
    		           	
    		        }
    		    	 
    		        g2.setColor(colorOfYAxisUnits );
    		        
    		    	if(freq == true){
    		    		label = (int)(y_max-l*i);
    		    		ylabel = String.format((String.valueOf(label))+ "");
    		    	}else{
    		    		double dblLabel = (y_max-l*i)/sampleLength;
    		    		ylabel = String.format(yAxisUnitsFormat, dblLabel)+ "";
    		    	}
    	 	    	
    		    	 
    		    	FontMetrics metrics = g2.getFontMetrics();
    		        int labelWidth = metrics.stringWidth(ylabel);
    		        
    		        g2.drawString(ylabel, x0 - labelWidth - 5, y0 + (metrics.getHeight() / 2) - 3);
    		         
    		    }
    	   		
    		    if(yMaxPos>(border_gap*(r+1)+rectangularHeight*r)){
    	            	    	
    		    	int i = 0;
    		    	
    		    	while(yMaxPos-(step_size*i)>(border_gap*(r+1)+rectangularHeight*r)){
    	                
    		    		g2.setColor(Color.BLACK);
    		 	    	int x0 = border_gap*(c+1)+rectangularWidth*c;
    		 	    	int x1 = x0+(int) width;
    		 	    	int y0 = (int) (yMaxPos-(step_size*i));
    		 	    	int y1 = y0;
    		 	    	//g2.drawLine(x0, y0, x1, y1);
    		    		 
    			        if(grid == true){
    				    	
    			        	Stroke oldStroke = g2.getStroke();
    			        	
    			        	g2.setStroke(new BasicStroke(4));
    			           	g2.setColor(grid_color);
    			           	g2.drawLine(x1, y1, (border_gap+rectangularWidth)*(c+1), y1); 
    			           	
    			    		if(i != numberYDivisions){
    			    			
    			    			g2.setStroke(new BasicStroke(2));
    				    		
    				    		int y0_half = (int) (y0+0.5*step_size);
    				    		int y1_half = y0_half;
    				    		g2.drawLine(x1, y1_half, (border_gap+rectangularWidth)*(c+1), y1_half); 
    			    			
    			    		}
    			           		           	
    			           	g2.setStroke(oldStroke);
    			           	
    			        }
    		 	    	
    		 	    	if(freq == true){
    			    		label = (int)(y_max+l*i);
    			    		ylabel = String.format((String.valueOf(label))+ "");
    			    	}else{
    			    		double dblLabel = (y_max+l*i)/sampleLength;
    			    		ylabel = String.format(yAxisUnitsFormat, dblLabel)+ "";
    			    	}
    		 	    	
    		 	    	g2.setColor(colorOfYAxisUnits );
    		 	    	FontMetrics metrics = g2.getFontMetrics();
    		 	        int labelWidth = metrics.stringWidth(ylabel);
    		 	        g2.drawString(ylabel, x0 - labelWidth - 5, y0 + (metrics.getHeight() / 2) - 3);
    		 	    	
    		 	        i++;
    		 	        
    		 	    }
    		    	
    		    }
    	    	
        		idx++;
        		
    		    if(idx>=nSamples){    	    		
    	    		break;    	    		
    	    	}
    		    
    	    }
    		
    	}	
    			
	    g2.setColor(Color.BLACK);
	    g2.setFont(orgFont);
	    
   	}
   	
   	
   	@SuppressWarnings("static-access")
	public void create_normalPDF(){
   		
   		if(freq == false){
   			
   			List<Integer> plotInfoOrg = plotInfo;
   			List<String>  plotTypeOrg = plotType;
   			ArrayList<List<Double>> xOrg = x;
   			ArrayList<List<Double>> yOrg = y;
   			
   		    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
   			
   	   		double [][] my    = new double [1][1];
   	   		double [][] sigma = new double [1][1];
   	   		
   	   		GenGraphics linePlot = new GenGraphics();

	   		int graphWidth  = getGraphWidth();
	   		int graphHeight = getGraphHeight();
   	   		
	   		linePlot.setGraphWidth(graphWidth);
   	   		linePlot.setGraphHeight(graphHeight);
	   		
   	   		int nSamples = x.size();
   	   		int nPlotsWithNormalPDF = plotIdxs4NormalPDF.size();
   	   		int counter = 0;
   	   		 
   	   		List<Integer> rowIdxs = new ArrayList<Integer>();
   	   		List<Integer> colIdxs = new ArrayList<Integer>();
   	   		
   	   		int idx = 0;
   	   		
   	   	    for(int i=0; i<numberOfPlotRows; i++){
   	   	    	
   	   	    	for(int j=0; j<numberOfPlotColumns; j++){
   	   	    		
   	   	    		if(counter == plotIdxs4NormalPDF.get(idx)){
   	   	    			
   	   	    			rowIdxs.add(i);
   	   	    			colIdxs.add(j);
   	   	    			idx++;
   	   	    			
   	   	    			if(idx>=nPlotsWithNormalPDF){
   	   	    				break;
   	   	    			}
   	   	    			
   	   	    		}
   	   	    		
   	   	    		if(idx>=nPlotsWithNormalPDF){
  	    				break;
  	    			}
   	   	    		
   	   	    		counter++;
   	   	    		
   	   	    	}
   	   	    	
   	   	    }
   	   				  	   		
   	   		for(int i=0; i<nPlotsWithNormalPDF; i++){
	   	   			
   	    	   	my[0][0]    = GeneralMath.mean(y.get(plotIdxs4NormalPDF.get(i)));
   	    	   	sigma[0][0] = Math.sqrt(GeneralMath.variance(y.get(plotIdxs4NormalPDF.get(i))));
   	    	   		
   	    	   	NormalDistribution normalDist = new NormalDistribution(my, sigma);
   	    	   		
   	    	   	int n_breakpoints           = breakpoints.get(i).size();
   	    	   	double [][] normalDensities = new double [n_breakpoints][1];
   	    	   	double [][] x               = new double [n_breakpoints][1];
   	    	   	
   	    	   	int sampleLength = y.get(i).size();
   	    	   			
   	    	   	for(int j=0; j<n_breakpoints; j++){
   	    	   			
   	    	   		double breakpoint = breakpoints.get(plotIdxs4NormalPDF.get(i)).get(j);
   	    	   			
   	    	   		normalDensities[j][0] = normalDist.get_univariateNormalPDF(breakpoint)*sampleLength;
   	    	   		x[j][0]               = breakpoint;
   	    	   			
   	    	   	}   	    	   		
   	    	   		
   	    	   	//prepare scaling of line plot for PDF.
   	    		int maxFreq = Utilities.getMaxFromIntList(counts.get(plotIdxs4NormalPDF.get(i)));				
   	    			
   	    		double scale   = (double)rectangularHeight/(double)maxFreq;          
   	    		double yMaxPos = (double)((border_gap+rectangularHeight)*(rowIdxs.get(i)+1))-scale*maxFreq;  
	  	       		  	
   	    	   	linePlot.plotLines(x,normalDensities, true, Color.RED);
   	    	   	
   	    	   	List<Point> graphPoints = new ArrayList<Point>();
   	      
			    graphPoints = linePlot.calc_scaled_points_4_plot_cell(rowIdxs.get(i), colIdxs.get(i), (nSamples+i), (nSamples+i), (double) maxFreq, yMaxPos);
		    
			    linePlot.graphPoints4Samples.add(graphPoints);
   	    	   	   	    	   	
   	   		}		
   	   		 	
   	   		int [] idx2remove = get_idxs_4_plot_type("H");	
   	   		remove_input_data(idx2remove);
   	   		
   	   		linePlot.createLinePlot();
   	   		
   	   		plotInfo = plotInfoOrg;
   	   		plotType = plotTypeOrg;
   	   		x        = xOrg;
   	   		y        = yOrg;
   	   		
   	   		idx2remove = get_idxs_4_plot_type("L");	
	   		remove_input_data(idx2remove);
   	   		
   		}

   	}
   	
   	
   	public static void plotNormalPDF(boolean doNormalPDF){
   		
   		plotNormalPDF = doNormalPDF;
   			
   	}
   	
   	
	public static void main(String[] args) {
		
		int maxDataPoints = 1000;
		
	 	double [][] x = new double [maxDataPoints][1];
	 	double [][] y = new double [maxDataPoints][1];
	 		
	 	Random r = new Random();
	 	   
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x[i][0] = i;
	 		y[i][0] = 10.0*r.nextGaussian();
	 	}
		
	 	double [][] x1 = new double [maxDataPoints][1];
	 	double [][] y1 = new double [maxDataPoints][1];
	 	
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x1[i][0] = i;
	 		y1[i][0] = 100.0+8.0*r.nextGaussian();
	 	}
	 	
	 	
	 	//double [] min = {70.0, -50.0, -15.0, -50.0, -50.0, -50.0};
	 	//double [] max = {130.0, 50.0, 25.0, 50.0, 50.0, 50.0};
	 	
	 	String [] titles = {"Hist1", "Hist2", "Hist3", "Hist4"};
	 	String [] subTitles = {"Sub1", "Sub2", "Sub3", "Sub4"};
	 	String [] yLabels = {"yLabel1", "yLabel2", "yLabel3", "yLabel4"};
	 	String [] xLabels = {"xLabel1", "xLabel2", "xLabel3", "xLabel4"};
	 	
	 	plotHistogram(y1,true,true);
	 	plotHistogram(y,true);
	 	plotHistogram(y,true,true);
	 	plotHistogram(y,true, true);
	 	plotHistogram(y1,true);
	 	plotHistogram(y,true, true);
	 	setNumberOfPlotColums(3);
	 	setNumberOfPlotRows(2);
	 	set_numberOfBins(70);
	 	setTitle(titles, null, null);
	 	setSubTitle1(subTitles, null, null);
	 	setYLabel(yLabels, null, null);
	 	setXLabel(xLabels, null, null);
	 	setGraphWidth(1000);
	 	setGraphHeight(600);
	 	setNumberOfDigits4XAxis(0);
	 	setNumberOfDigits4YAxis(2);
	 	setFontOfXAxisUnits("bold",11);
	 	setFontOfYAxisUnits("bold",11);
	 	//set_max_x_value(max);
	 	//set_min_x_value(min);
	 	freqHist(false);	 	
	 	//noLinesAroundBars();
	 	plot();
	 	
	}
	
}
