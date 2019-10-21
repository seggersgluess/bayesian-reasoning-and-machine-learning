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
    public static ArrayList<List<Double>> counts     = new ArrayList<List<Double>>();
    public static List<Double> maxValue = new ArrayList<Double>();
    public static List<Double> minValue = new ArrayList<Double>();
	public static boolean freq = true;
	
	public static List<Color> barColor = new ArrayList<Color>();
	public static List<Color> lineColor = new ArrayList<Color>() ; //Color.WHITE;

	public static boolean plotNormalPDF = false;
	public static List<Integer> plotIdxs4NormalPDF = new ArrayList<Integer>();
	public static ArrayList<ArrayList<List<Double>>> normalPDFs = new ArrayList<ArrayList<List<Double>>>();
	public static int pdfLineWidth = 2;
	
	
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
	    		createNormalPDF();   		
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
    	
		if(useDifferentColors == true){			
			barColor = new ArrayList<Color>(nSamples);			
			for(int i=0; i<nSamples; i++) {
				int red   = getDefaultColors()[i][0];
				int green = getDefaultColors()[i][1];
				int blue  = getDefaultColors()[i][2];
				
				barColor.add(new Color(red,green,blue,210));
			}					
		}
    	
	    set_breakpoints();
	    count_elements_between_breakpoints();
	    	
	    int rectangularWidth  = (int) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((pref_h-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
    	set_x_axis();
    	set_y_axis();
	    
    	int idx       = 0;
    	int plotIdx   = -1;
    	int nPlotIdxs = 0;
    	
    	ArrayList<List<Double>> listOfCounts = new ArrayList<List<Double>>();
    	List<Integer> pdfIdxs = null;
    	
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(idx>=nSamples){	    		
	    		break;	    		
	    	}
    		
    		int yBottom = (border_gap+rectangularHeight)*(r+1);
    		
    	    for(int c=0; c<numberOfPlotColumns; c++){
  				
    	    	if(idx>=nSamples){    	    		
    	    		break;    	    		
    	    	}
    	    	
    	    	int curPlotIdx = plotInfo.get(idx);
    	    	
    	    	if(plotIdx != curPlotIdx){
    	    		plotIdx = curPlotIdx; 
    	    		
    	    		int [] plotIdxs = Utilities.get_idx(plotInfo, plotIdx);
    	    		nPlotIdxs  = plotIdxs.length;
    	    		listOfCounts = new ArrayList<List<Double>>(nPlotIdxs); 
    	    		
    	    		for(int p=0; p<nPlotIdxs; p++){
    	    			listOfCounts.add(counts.get(idx+p));
    	    		}
    	    		
    	    		//Calculate pdf for plot selection
    	    		//int [] normalPlotIdx = Utilities.get_idx(plotIdxs4NormalPDF, idx);
        	    	if(plotNormalPDF==true){
        	    		int nPDFPlots = plotIdxs4NormalPDF.size();
        	    		pdfIdxs = new ArrayList<Integer>();
        	    		for(int i=0; i<nPDFPlots; i++){
        	    			int pdfPlotIdx = plotInfo.get(plotIdxs4NormalPDF.get(i));
        	    			if(pdfPlotIdx==plotIdx){    
        	    				pdfIdxs.add(i);
        	    				ArrayList<List<Double>> normalPDF = calcNormalPDF(plotIdxs4NormalPDF.get(i));
                	    		normalPDFs.add(normalPDF);
        	    			}    	    			
        	    		}
        	    	}
    	    		    	    		
    	    	}
    	    	   	    	
    	    	double stepSize = (double) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/nBins;
    	    	double maxFreq  = 0.0;
    	    	
    	    	//Determine maximum regarding pdfs
    	    	int [] normalPlotIdx = Utilities.get_idx(plotIdxs4NormalPDF, idx);
    	    	if(normalPlotIdx[0] != -1){
    	    		ArrayList<List<Double>> countsAndPDf = listOfCounts;
    	    		int startIdx = normalPDFs.size()-pdfIdxs.size();
    	    		for(int i=0; i<pdfIdxs.size(); i++){
    	    			countsAndPDf.add(normalPDFs.get(startIdx+i).get(1));
    	    		}
    	    		
    	    		maxFreq = Utilities.getMaxFromDblList(countsAndPDf);
    	    		
    	    	}else{
    	    		maxFreq = Utilities.getMaxFromDblList(listOfCounts);
    	    	}
    	    	   		
    	    	double scale = (double)rectangularHeight/maxFreq;
    	    			    		
    	    	int x0 = 0;
    	    	int y0 = 0;
    	    	int x1 = 0;
    	    	int y1 = 0;
    	    	
    	    	for(int p=0; p<nPlotIdxs; p++){
    	    		  	    		
        	    	for(int i=0; i<nBins; i++){

        	    		x0 = (int) (leftPosition+border_gap*(c+1)+rectangularWidth*c + i*stepSize);
        	    		x1 = (int) stepSize;
        	    			
        	    		y0 = (int)(upperPosition+(border_gap*(r+1)+rectangularHeight*(r+1))-scale*(listOfCounts.get(p)).get(i));						
        	    		y1 = (int)(scale*(listOfCounts.get(p)).get(i));
        								
        				g2.setColor(barColor.get(idx));
        				g2.fillRect(x0, y0, x1, y1);
        					
        				g2.setColor(lineColor.get(idx));
        				g2.drawRect(x0, y0, x1, y1);
        				
        				int xLeft   = leftPosition+border_gap*(c+1)+rectangularWidth*c;
            	    	int xRight  = leftPosition+(border_gap+rectangularWidth)*(c+1);
         			
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
   
	}
	
    
    public static Color getDefaultBarColor() {
    	
    	Color defaultBarColor = new Color(135,206,250,180);
    	
    	return defaultBarColor;
    	
    }
    
    
    public static Color getDefaultBarLineColor() {
    	
    	Color defaultLineColor = new Color(0,0,128,180);
    	
    	return defaultLineColor;
    	
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
    
    
    public static void setPDFLineWidth(int width){
    	pdfLineWidth = width;
    }
    
    
    public static void freqHist(boolean freqFlag){
    	
    	//if true then frequencies else relative counts (which sum to one).
    	freq = freqFlag;
    	
    }
    
    
    public static void set_breakpoints(){
    	
    	if(nBins == 0){
    		throw new RuntimeException("Number of bins not set.");
    	}
    	
    	int nSamples  = y.size();
    	
    	int plotIdx   = -1;
    	int nPlotIdxs = 0;
    	ArrayList<List<Double>> listOfObs = new ArrayList<List<Double>>();
    	
    	boolean boundsSetted = false;
    	
    	if(maxValue.size()!=0){
    		boundsSetted = true;
    	}
    	
		double max = 0.0;
		double min = 0.0;
    	
    	for(int i=0; i<nSamples; i++){
    		       		
    		if(boundsSetted == false){
    			
    	    	int curPlotIdx = plotInfo.get(i);
    	    	
    	    	if(plotIdx != curPlotIdx){
    	    		plotIdx = curPlotIdx; 
    	    		
    	    		int [] plotIdxs = Utilities.get_idx(plotInfo, plotIdx);
    	    		nPlotIdxs  = plotIdxs.length;
    	    		listOfObs = new ArrayList<List<Double>>(nPlotIdxs); 
    	    		
    	    		for(int p=0; p<nPlotIdxs; p++){    	    			
    	    			listOfObs.add(y.get(i+p));
    	    		}
    	    			
        	    	max = Utilities.getMaxFromDblList(listOfObs);
        	    	min = Utilities.getMinFromDblList(listOfObs);
        	    	
        	    	for(int p=0; p<nPlotIdxs; p++){
        	    		maxValue.add(max);
            	    	minValue.add(min);
    	    		}
    	    		
    	    	}
    			   		
    		}else{
    			
    			max = maxValue.get(i);
    			min = minValue.get(i);
    			
    		}
    		
        	double length = max-min;
        	double stepSize = length/nBins;
        	
        	List<Double> newBreaks = new ArrayList<Double>();
        		
        	for(int j=0; j<nBins+1; j++){
        		
        		newBreaks.add(min+j*stepSize);
        		
        	}
    		
        	breakpoints.add(newBreaks);
        	
    	}
     	
    }
    
    
    public static void count_elements_between_breakpoints(){
    	
    	int nSamples = y.size();
    	
    	for(int i=0; i<nSamples; i++){
    		
        	List<Double> sample = y.get(i);
        	int sampleLength    = sample.size();
        	
        	List<Double> newCounts = new ArrayList<Double>(nBins);
 
        	for(int j=0; j<nBins; j++){
        		
        		double counter = 0;
        		
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
	    
    	//maxValue.add(Utilities.getMax(y.get(y.size()-1)));
    	//minValue.add(Utilities.getMin(y.get(y.size()-1)));
    	
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
    	
    	barColor.add(getDefaultBarColor());
    	lineColor.add(getDefaultBarLineColor());
    	
    }
	
    
    public static void plotHistogram(double [][] values, boolean newPlot, Color barCol){
   	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	    
    	//maxValue.add(Utilities.getMax(y.get(y.size()-1)));
    	//minValue.add(Utilities.getMin(y.get(y.size()-1)));
    	
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
    	
    	barColor.add(barCol);
    	lineColor.add(getDefaultBarLineColor());
    	
    }
    
    
    public static void plotHistogram(double [][] values, boolean newPlot, boolean plotNormalPDF){
   	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	    
    	//maxValue.add(Utilities.getMax(y.get(y.size()-1)));
    	//minValue.add(Utilities.getMin(y.get(y.size()-1)));
    	
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
    		
    		plotIdxs4NormalPDF.add(plotInfo.size()-1);
    		
    		plotNormalPDF(true);
    		
    	}
    	
    	barColor.add(getDefaultBarColor());
    	lineColor.add(getDefaultBarLineColor());
    	
    }
     
    
    public static void plotHistogram(double [][] values, boolean newPlot, boolean plotNormalPDF, Color barCol){
    	   
    	setDefaultDesign();
    	
    	double [][] x_values = new double [values.length][1];
    	
    	convert_input_data(x_values, values, "H");
	    
    	//maxValue.add(Utilities.getMax(y.get(y.size()-1)));
    	//minValue.add(Utilities.getMin(y.get(y.size()-1)));
    	
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
    		
    		plotIdxs4NormalPDF.add(plotInfo.size()-1);
    		
    		plotNormalPDF(true);
    		
    	}
    	
    	barColor.add(barCol);
    	lineColor.add(getDefaultBarLineColor());
    	
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
    	
    	int nColors = barColor.size();
    	
    	for(int i=0; i<nColors; i++) {
    		lineColor.add(barColor.get(i));
    	}
    	
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
	    
	    double step_size      = (double) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns)/numberXDivisions;
	    int rectangularWidth  = (int) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((pref_h-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
	    int idx       = 0;
	    int plotIdx   = -1;
	    int nPlotIdxs = 0;
	    
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(idx>=nSamples){	    		
	    		break;	    		
	    	}
	    	
    	    for(int c=0; c<numberOfPlotColumns; c++){
    	    		
    		    int curPlotIdx = plotInfo.get(idx);
    	    	
    	    	if(plotIdx != curPlotIdx){
    	    		plotIdx = curPlotIdx; 
    	    		
    	    		int [] plotIdxs = Utilities.get_idx(plotInfo, plotIdx);
    	    		nPlotIdxs = nPlotIdxs+plotIdxs.length;	
    	    		
    	    	}
    	    	
    		    double width      = 0.0;//0.01*(getHeight()-2*border_gap);
    			double xLength    = (maxValue.get(idx)-minValue.get(idx))/numberXDivisions; 
    		    			
    		    int y0 = upperPosition+border_gap*(r+1)+rectangularHeight*(r+1);
    	    	int y1 = y0 - (int) width;
    		    
    		    for(int i = 0; i<numberXDivisions+1; i++) {
    			       
    		    	g2.setColor(Color.BLACK);
    		    	int x0 = (int) (leftPosition+step_size*(i)+border_gap*(c+1)+rectangularWidth*c);
    		    	int x1 = x0;
    		    	
    		    	//g2.drawLine(x0, y0, x1, y1);
    			    
    		    	if(grid == true){
    	                
    		    		Stroke oldStroke = g2.getStroke();
    		        	
    		        	g2.setStroke(new BasicStroke(4));	
    		    		g2.setColor(grid_color);
    		    		g2.drawLine(x0, y1, x1, (upperPosition+border_gap*(r+1)+rectangularHeight*r)); 
    		    		
    		    		if(i != numberXDivisions){
    		    			
    		    			g2.setStroke(new BasicStroke(2));
    			    		
    			    		int x0_half = (int) (x0+0.5*step_size);
    			    		int x1_half = x0_half;
    			    		g2.drawLine(x0_half, y1, x1_half, (upperPosition+border_gap*(r+1)+rectangularHeight*r)); 
    		    			
    		    		}
    		    		
    		    		g2.setStroke(oldStroke);
    			    	
    		    	}
    		    	
    		    	g2.setColor(colorOfXAxisUnits);
    		    	String xlabel = String.format(xAxisUnitsFormat, minValue.get(idx)+xLength*(i))+ "";
    			    FontMetrics metrics = g2.getFontMetrics();
    			    int labelWidth = metrics.stringWidth(xlabel);
    			    g2.drawString(xlabel, x0 - labelWidth/2, y0 + metrics.getHeight() + 3);  	    	
    		    	
    		    }
    	    	
    		    idx = nPlotIdxs;
    		    
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
   		
	    int rectangularWidth = (int) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((pref_h-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
	    
	    int idx       = 0;
	    int plotIdx   = -1;
	    int nPlotIdxs = 0;
	    
    	for(int r=0; r<numberOfPlotRows; r++){
    		
    		if(idx>=nSamples){	    		
	    		break;	    		
	    	}
    	    		
    		for(int c=0; c<numberOfPlotColumns; c++){
    	    	
    		    int curPlotIdx = plotInfo.get(idx);
    	    	
    	    	if(plotIdx != curPlotIdx){
    	    		plotIdx = curPlotIdx; 
    	    		
    	    		int [] plotIdxs = Utilities.get_idx(plotInfo, plotIdx);
    	    		nPlotIdxs = nPlotIdxs+plotIdxs.length;	
    	    		
    	    	}
    			
    		    int sampleLength = (int)GeneralMath.sumDblList(counts.get(idx));
    		    double maxFreq = 0.0;
    		    
    		    int [] normalPlotIdx = Utilities.get_idx(plotIdxs4NormalPDF, idx);
    	    	if(normalPlotIdx[0] != -1){
    	    		ArrayList<List<Double>> normalPDF = calcNormalPDF(idx);
    	    		normalPDFs.add(normalPDF);
    	    		ArrayList<List<Double>> countsAndPDf = new ArrayList<List<Double>>(2);
    	    		countsAndPDf.add(counts.get(idx));
    	    		countsAndPDf.add(normalPDF.get(1));
    	    		maxFreq = Utilities.getMaxFromDblList(countsAndPDf);
    	    		normalPDFs = new ArrayList<ArrayList<List<Double>>>();
    	    	}else{
    	    		maxFreq     = Utilities.getMaxFromDblList(counts.get(idx));
    	    	}
    		        	    		    
    			double scale     = (double) ((rectangularHeight)/maxFreq);
    		    int yMaxPos      = (int) ((double)(upperPosition+border_gap*(r+1)+rectangularHeight*(r+1))-scale*maxFreq);
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
    		    	int x0 = leftPosition+border_gap*(c+1)+rectangularWidth*c;
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
    		           	g2.drawLine(x1, y1, (leftPosition+border_gap+rectangularWidth)*(c+1), y1); 
    			    	
    		    		if(i != numberYDivisions){
    		    			
    		    			g2.setStroke(new BasicStroke(2));
    			    		
    			    		int y0_half = (int) (y0+0.5*step_size);
    			    		int y1_half = y0_half;
    			    		g2.drawLine(x1, y1_half, (leftPosition+border_gap+rectangularWidth)*(c+1), y1_half); 
    		    			
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
    	   		
    		    if(yMaxPos>(upperPosition+border_gap*(r+1)+rectangularHeight*r)){
    	            	    	
    		    	int i = 0;
    		    	
    		    	while(yMaxPos-(step_size*i)>(upperPosition+border_gap*(r+1)+rectangularHeight*r)){
    	                
    		    		g2.setColor(Color.BLACK);
    		 	    	int x0 = leftPosition+border_gap*(c+1)+rectangularWidth*c;
    		 	    	int x1 = x0+(int) width;
    		 	    	int y0 = (int) (yMaxPos-(step_size*i));
    		 	    	int y1 = y0;
    		 	    	//g2.drawLine(x0, y0, x1, y1);
    		    		 
    			        if(grid == true){
    				    	
    			        	Stroke oldStroke = g2.getStroke();
    			        	
    			        	g2.setStroke(new BasicStroke(4));
    			           	g2.setColor(grid_color);
    			           	g2.drawLine(x1, y1, (leftPosition+border_gap+rectangularWidth)*(c+1), y1); 
    			           	
    			    		if(i != numberYDivisions){
    			    			
    			    			g2.setStroke(new BasicStroke(2));
    				    		
    				    		int y0_half = (int) (y0+0.5*step_size);
    				    		int y1_half = y0_half;
    				    		g2.drawLine(x1, y1_half, (leftPosition+border_gap+rectangularWidth)*(c+1), y1_half); 
    			    			
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
 	        
    		 	    }
    		    	
    		    }
    	    	
    		    idx = nPlotIdxs;
        		
    		    if(idx>=nSamples){    	    		
    	    		break;    	    		
    	    	}
    		    
    	    }
    		
    	}	
    			
	    g2.setColor(Color.BLACK);
	    g2.setFont(orgFont);
	    
   	}
   	
   	
   	@SuppressWarnings("static-access")
	public void createNormalPDF(){
   		
   		if(freq == false){
   			
   			List<Integer> plotInfoOrg = plotInfo;
   			List<String>  plotTypeOrg = plotType;
   			ArrayList<List<Double>> xOrg = x;
   			ArrayList<List<Double>> yOrg = y;
   			
   		    int rectangularHeight = (int)((pref_h-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
   			
   	   		GenGraphics linePlot = new GenGraphics();

	   		int graphWidth  = getGraphWidth();
	   		int graphHeight = getGraphHeight();
   	   		
	   		linePlot.setGraphWidth(graphWidth);
   	   		linePlot.setGraphHeight(graphHeight);
	   		
   	   		int nSamples = x.size();
   	   		int nPlotsWithNormalPDF = plotIdxs4NormalPDF.size();
   	   		int nPlots = numberOfPlotRows*numberOfPlotColumns; 
   	   		int maxPlotInfo = Utilities.getMaxFromIntList(plotInfo)+1; 
   	   		
   			for(int i=0; i<nPlotsWithNormalPDF; i++){
	   	   					
   	   			int normalPDFidx = plotIdxs4NormalPDF.get(i);
   	   			
   	   			if(plotInfo.get(normalPDFidx)>nPlots-1){
   	   				break;
   	   			}
   	   			
   	   			ArrayList<List<Double>> pdf = normalPDFs.get(i);
   	   				
   	    	   	int n_breakpoints           = breakpoints.get(normalPDFidx).size();
   	    	   	double [][] normalDensities = new double [n_breakpoints][1];
   	    	   	double [][] x               = new double [n_breakpoints][1];
 	    	   		
   	    	   	for(int j=0; j<n_breakpoints; j++){  	    	   				    	   			 	    	   		
   	    	   		x[j][0]               = pdf.get(0).get(j);
   	    	   		normalDensities[j][0] = pdf.get(1).get(j);	   	    	   		
   	    	   	}   	    	   		

   	    		linePlot.plotLines(x, normalDensities, true, Color.RED);
  	    	   	
   	    		//Determine position of PDF
   	    		int pdfPlotIdx = plotInfo.get(plotIdxs4NormalPDF.get(i));
   	   			List<Integer> plotPos = getMultiPlotPositions(pdfPlotIdx);
   	    		int rowIdx = plotPos.get(0);
   	    		int colIdx = plotPos.get(1);
   	   			
   	    		int [] graphIdxs = Utilities.get_idx(plotInfo, pdfPlotIdx);
   	    		ArrayList<List<Double>> countsAndPDFs = new ArrayList<List<Double>>();
   	    		for(int j=0; j<graphIdxs.length; j++){
   	    			countsAndPDFs.add(counts.get(graphIdxs[j]));
   	    		}
   	    		
   	    		for(int j=0; j<nPlotsWithNormalPDF; j++){
   	    			if(plotInfo.get(plotIdxs4NormalPDF.get(i))==pdfPlotIdx){
   	    				countsAndPDFs.add(normalPDFs.get(i).get(1));
   	    			}
   	    		}
   	    		
   	    		//prepare scaling of line plot for PDF.
   	    		double maxFreq    = Utilities.getMaxFromDblList(countsAndPDFs);
   	    		double maxDensity = Utilities.getMax(normalDensities);
   	    		
   	    	   	List<Point> graphPoints = new ArrayList<Point>();
   	    		   	    	   	
   	    		if(maxFreq>maxDensity){
   	    			double scale   = (double)rectangularHeight/maxFreq;          
   	   	    		double yMaxPos = (double)((border_gap+rectangularHeight)*(rowIdx+1))-scale*maxFreq;
   	   	    		
   	   	    		graphPoints = linePlot.calc_scaled_points_4_plot_cell(rowIdx, colIdx, maxPlotInfo+i, nSamples+i, maxFreq, yMaxPos);
   	    		
   	    		}else{
   	    			
   	    			graphPoints = linePlot.calc_scaled_points_4_plot_cell(rowIdx, colIdx, maxPlotInfo+i, nSamples+i);
   				    
   	    		}
   	    		
			    linePlot.graphPoints4Samples.add(graphPoints);
   	    	   	   	    	   	
   	   		}		
   	   		 	
   	   		int [] idx2remove = get_idxs_4_plot_type("H");	
   	   		remove_input_data(idx2remove);
   	   		
   	   		Stroke orgStroke = linePlot.g2.getStroke();
   	   		
   	   		linePlot.graph_stroke = new BasicStroke(pdfLineWidth);
   	   		linePlot.createLinePlot();
   	   		linePlot.graph_stroke = orgStroke;
   	   		
   	   		plotInfo = plotInfoOrg;
   	   		plotType = plotTypeOrg;
   	   		x        = xOrg;
   	   		y        = yOrg;
   	   		
   	   		idx2remove = get_idxs_4_plot_type("L");	
	   		remove_input_data(idx2remove);
   	   		
   		}

   	}
   	
   	
   	public static List<Integer> getMultiPlotPositions(int plotIdx){
   		
   		int r=0;
   		int c=0;
   		int counter = 0;
   	    List<Integer> multiPlotPos = new ArrayList<Integer>(2);
   		
   		while(r<numberOfPlotRows){
	   	    		
   			c=0;
   			
	   	    while(c<numberOfPlotColumns){
	   	    	
	   	    	if(plotIdx==counter){
	   	    		multiPlotPos.add(r);
	   	    		multiPlotPos.add(c);
	   	    	}
	   	    	
	   	    	counter++;
	   	    	c++;
	   	    		
	   	    }
   		
	   	    r++;
	   	    	
   		}
	   	 
   		return multiPlotPos;
   		
   	}
   	
   	
   	@SuppressWarnings("static-access")
	public static ArrayList<List<Double>> calcNormalPDF(int sampleIdx){
   		
	   	double [][] my    = new double [1][1];
	   	double [][] sigma = new double [1][1];
   		
   	   	my[0][0]    = GeneralMath.mean(y.get(sampleIdx));
   	   	sigma[0][0] = Math.sqrt(GeneralMath.variance(y.get(sampleIdx)));
   	   		
   	   	NormalDistribution normalDist = new NormalDistribution(my, sigma);
   	   		
   	   	int n_breakpoints             = breakpoints.get(sampleIdx).size();
   	    
   	    ArrayList<List<Double>> res   = new ArrayList<List<Double>>(2);
   	   	List<Double> normalDensities  = new ArrayList<Double>(n_breakpoints);
   	   	List<Double> scaledDensities  = new ArrayList<Double>(n_breakpoints);
   	    List<Double> x                = new ArrayList<Double>(n_breakpoints);
   	   	
   	   	int sampleLength = y.get(sampleIdx).size();
   	   		
   	   	for(int j=0; j<n_breakpoints; j++){
   	   			
   	   		double breakpoint = breakpoints.get(sampleIdx).get(j);
   	   			
   	   		x.add(breakpoint);
   	   		normalDensities.add(normalDist.get_univariateNormalPDF(breakpoint));
   	   		   			
   	   	} 
   		
   	   	//Scale the normal pdfs such that sum of pdfs equals 1.0
   	   	double pdfSum = GeneralMath.sumDblList(normalDensities);
   	   	 	
   	   	for(int j=0; j<n_breakpoints; j++){	   	
   	   		double scaledDensity = (normalDensities.get(j)/pdfSum)*sampleLength;
   	   		scaledDensities.add(scaledDensity);	   		   			
   	   	} 
   	   	
   	   	res.add(x);
   	   	res.add(scaledDensities);
   	   	
   	   	return res;
   	   	
   	}
   	
   	
   	public static void plotNormalPDF(boolean doNormalPDF){
   		
   		plotNormalPDF = doNormalPDF;
   			
   	}
   	
   	
   	public static void histExample1() {
   		
   		int maxDataPoints = 1000;
		
	 	double [][] x = new double [maxDataPoints][1];
	 	double [][] y = new double [maxDataPoints][1];
	 		
	 	Random r = new Random();
	 	   
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x[i][0] = i;
	 		y[i][0] = 15.0+1.0*r.nextGaussian();
	 	}
		
	 	double [][] x1 = new double [maxDataPoints][1];
	 	double [][] y1 = new double [maxDataPoints][1];
	 	
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x1[i][0] = i;
	 		y1[i][0] = 10.0+2.0*r.nextGaussian(); //Math.exp(0.9*r.nextGaussian());
	 	}
	 	
	 	double [][] x2 = new double [maxDataPoints][1];
	 	double [][] y2 = new double [maxDataPoints][1];
	 	
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x2[i][0] = i;
	 		y2[i][0] = Math.exp(0.5*r.nextGaussian());
	 	}
	 	
	 		 	
	 	//double [] min = {-3.0, -5.0, -15.0, -5.0, -5.0, -10.0};
	 	//double [] max = {5.0, 5.0, 25.0, 1.0, 5.0, 10.0};
	 	
	 	String [] titles = {"Hist1", "Hist2", "Hist3", "Hist4"};
	 	String [] subTitles = {"Sub1", "Sub2", "Sub3", "Sub4"};
	 	String [] yLabels = {"yLabel1", "yLabel2", "yLabel3", "yLabel4"};
	 	String [] xLabels = {"xLabel1", "xLabel2", "xLabel3", "xLabel4"};
	 	
	 	plotHistogram(y,true,true);
	 	plotHistogram(y1,false,true);
	 	plotHistogram(y2,true,true);
	 	plotHistogram(y,true,false);
	 	plotHistogram(y1,false,true);
	 	plotHistogram(y2,true,true);
	 	plotHistogram(y,true,true);
	 	//plotHistogram(y1,false,true);
	 	setNumberOfPlotColums(3);
	 	setNumberOfPlotRows(1);
	 	set_numberOfBins(50);
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
	 	setPDFLineWidth(2);
	 	//noLinesAroundBars();
	 	plot();
   		
   	}
   	
   	
	public static void main(String[] args) {
		
		histExample1();
	 	
	}
	
}
