package Graphics;
import java.awt.FontMetrics;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.List;
import javax.swing.*;

import Utilities.Utilities;

@SuppressWarnings("serial")
public class GraphicDevice extends JPanel{
	
	public static Color backroundColor  = new Color(200, 200, 200, 200);
	public static Color rectColor       = Color.WHITE; 
	
    public static boolean grid = true;
    public static  Color grid_color = new Color(200, 200, 200, 200);
	
    public static int pref_w = 400;
    public static int pref_h = 400;
    protected static int border_gap = (int) (pref_h*0.15);
    
    protected static ArrayList<List<Double>> x    = new ArrayList<List<Double>>();
    protected static ArrayList<List<String>> xStr = new ArrayList<List<String>>();
    protected static ArrayList<List<Double>> y    = new ArrayList<List<Double>>();
    protected static List<String> plotType        = new ArrayList<String>();  
    
    public static List<Integer> plotInfo       = new ArrayList<Integer>();
    
    public static int numberOfPlotColumns = 1;
    public static int numberOfPlotRows   = 1;
    
    static Graphics2D g2;
    
    static Color textColor = Color.BLACK;
    
    static Font fontOfXAxisUnits = new Font(null, Font.PLAIN, 12);
    static Color colorOfXAxisUnits = Color.BLACK;
    
    static Font fontOfYAxisUnits = new Font(null, Font.PLAIN, 12);
    static Color colorOfYAxisUnits = Color.BLACK;
    
    static int numberOfDigitsXAxisUnits = 2;
    static String xAxisUnitsFormat      = "%.2f%n";
    static int numberOfDigitsYAxisUnits = 2;  
    static String yAxisUnitsFormat      = "%.2f%n";
     
    static List<String> title = new ArrayList<String>();
    static Font titleFont     = new Font(null, Font.BOLD, 12); 
    
    static List<String> subTitle1 = new ArrayList<String>();
    static String subTitle2       = null;
    static Font subTitle1Font     = new Font(null, Font.PLAIN, 10);
    static Font subTitle2Font     = new Font(null, Font.PLAIN, 10);
    
    static List<String> xLabel = new ArrayList<String>();   
    static String xSubLabel    = null;
    static Font xLabelFont     = new Font(null, Font.BOLD, 10);
    static Font xSubLabelFont  = new Font(null, Font.PLAIN, 10);
    
    static List<String> yLabel = new ArrayList<String>();
    static Font yLabelFont     = new Font(null, Font.BOLD, 10);
    
    public static int numberYDivisions = 5;
    public static int numberXDivisions = 10;
    
    public static boolean useDifferentColors = true;
    
    
    public void set_rectangular(){
    	
    	g2.setColor(rectColor);
    	
    	int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
    	int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
    	
    	int counter = 0;
    	
    	int nPlots = Utilities.getMaxFromIntList(plotInfo)+1;
    	
    	for(int i=0; i<numberOfPlotRows; i++){
    		
    		if(counter>=nPlots){
    			break;
    		}
    		
    		int y0 = border_gap*(i+1)+rectangularHeight*i;
    		
    		for(int j=0; j<numberOfPlotColumns; j++){
    				
    			int x0 = border_gap*(j+1)+rectangularWidth*j;
    			
        		g2.fillRect(x0, y0, rectangularWidth, rectangularHeight);
    		
        		counter++;
        		
        		if(counter>=nPlots){
        			break;
        		}
        		
    		}
    			
    	}
    	
    }
    
    
    public static int [][] getDefaultColors(){
    	
    	int [][] defaultColors = {	{230,159,0},
    								{86,180,233},
    								{0,158,115},
    								{240,228,66},
    								{0,114,178},
    								{213,94,0},
    								{204,121,167},	
    												};
    	
    	return defaultColors;
    	
    }
    
    
    public static void setRectColor(Color color){
    	
    	rectColor = color;
    	
    }
    
    
    public static void setNumberOfPlotColums(int nCols){
    	
    	numberOfPlotColumns = nCols;
    	
    }
    
    
    public static void setNumberOfPlotRows(int nRows){
    	
    	numberOfPlotRows = nRows;
    	
    }
   
    
    public static void setGridColor(Color color){
    	
    	grid_color = color;
    	
    }
    
    
    public static void setTextColor(Color color){
    	
    	textColor = color;
    	
    }
    
    
    public void set_title2plot(){
    	
    	if(title != null){
    		
        	int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
        	int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
    		
        	int nTitles        = title.size();
        	int nPlots         = Utilities.getMaxFromIntList(plotInfo)+1;
    		String title4plot  = "";
        	
    		int scale = 0;
    		
    		if(subTitle1 != null){
    			
    			scale = -1;
    			
    			if(subTitle2 != null){
    				
    				scale = 2;
    				
    			}
    			
    			scale = scale + 1;
    			
    		}
    		
    		g2.setColor(textColor); 
    		FontMetrics metrics = g2.getFontMetrics();
    		
    		int counter = 0;
    		
    		for(int r=0; r<numberOfPlotRows; r++){
    				
	    		if(counter >= nPlots){
	    			break;
	    		}
    			
	    		int yPos = (int)((border_gap*(r)+rectangularHeight*r)+border_gap/(scale+9.0/5.5));	    		
	    		
    			for(int c=0; c<numberOfPlotColumns; c++){
    				
    				if(counter<nTitles){
    					title4plot = title.get(counter);
    				}else{
    					title4plot = "";
    				}
    				
    	    		int xPos = border_gap*(c+1)+rectangularWidth*c + rectangularWidth/2;
    	    			    		
    	    		int titleWidth = metrics.stringWidth(title4plot)/2;
    	    		 
    	    		Font orgFont = g2.getFont();
    	    		    		
    	    		g2.setFont(titleFont);
    	    		g2.drawString(title4plot, xPos-titleWidth, yPos + (metrics.getHeight()/2));    		
    	    		g2.setFont(orgFont);
    	    		
    	    		if(subTitle1 != null){
    	    			
    	    			int nSubTitles = subTitle1.size();
    	    			
    	    			String subTitle4plot = "";
    	    			
    	    			if(counter<nSubTitles){
        					subTitle4plot = subTitle1.get(counter);
        				}else{
        					subTitle4plot = "";
        				}
    	    			
    	    			int subTitle1Width = metrics.stringWidth(subTitle4plot)/2;
    	        		
    	    			g2.setFont(subTitle1Font);
    	        		g2.drawString(subTitle4plot, xPos - subTitle1Width, yPos + (int)(metrics.getHeight()*6.0/5.0));
    	        		g2.setFont(orgFont);
    	        		
    	    		}
    	    		
    	    		if(subTitle2 != null){
    	    			
    	    			int subTitle1Width = metrics.stringWidth(subTitle2)/2;
    	        		
    	    			g2.setFont(subTitle2Font);
    	        		g2.drawString(subTitle2, xPos - subTitle1Width, yPos + (metrics.getHeight()*5/2));
    	        		g2.setFont(orgFont);
    	        		
    	    		}
    				
    	    		counter++;
    	    		
    	    		if(counter >= nPlots){
    	    			break;
    	    		}
    	    			    		
    			}    			
    			
    		}
    		
    	}
    	
    	g2.setColor(Color.BLACK);
    	
    }
    
    
    public void set_xLabel2plot(){
    	
    	if(xLabel != null){
    		
    		int scale = 0;
    		
    		if(xSubLabel != null){
    			
    			scale = 1/2;
    			
    		}
    		
    		g2.setColor(textColor);
    		
    		int nLabels        = xLabel.size();
    		int nPlots         = Utilities.getMaxFromIntList(plotInfo)+1;
    		
    		String label       = "";
    		
    		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
        	int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
    		
    		int counter = 0;
    		
    		for(int r=0; r<numberOfPlotRows; r++){
				
	    		if(counter >= nPlots){
	    			break;
	    		}
	    		
	    		for(int c=0; c<numberOfPlotColumns; c++){
    				
	    			if(counter < nLabels){
	    				label = xLabel.get(counter);
	    			}else{
	    				label = "";
	    			}
	    		
	        		FontMetrics metrics = g2.getFontMetrics();
	        		
	        		int xPos = border_gap*(c+1)+rectangularWidth*c + (rectangularWidth)/2;
	        		int yPos = ((border_gap+rectangularHeight)*(r+1)+border_gap/(scale+3));
	        		
	        		int labelWidth = metrics.stringWidth(label)/2;
	        		
	        		Font orgFont = g2.getFont();
	        		
	        		g2.setFont(xLabelFont);
	        		g2.drawString(label, xPos-labelWidth, yPos+(metrics.getHeight()/2));    		
	        		g2.setFont(orgFont);
	        		
	        		if(xSubLabel != null){
	        			
	        			int subLabelWidth = metrics.stringWidth(xSubLabel)/2;
	            		
	        			g2.setFont(xSubLabelFont);
	            		g2.drawString(xSubLabel, xPos-subLabelWidth, yPos + (metrics.getHeight()*3/2));
	            		g2.setFont(orgFont);
	            		
	        		}
	    			
	    			counter++;
		    		
	    			if(counter >= nPlots){
	    				break;
	    			}
	    			
	    		}
	    		
	    		
    		}
    		
    		g2.setColor(Color.BLACK);
    		
    	}
    	
    }
    
    
    public void set_yLabel2plot(){
    	
    	if(yLabel != null){
    		
    		g2.setColor(textColor);  
    		
    		Font orgFont = g2.getFont();
    		
    		g2.setFont(fontOfYAxisUnits);   		
    		FontMetrics metrics4Label = g2.getFontMetrics();
    		int labelUnitWidth = metrics4Label.stringWidth(String.format(yAxisUnitsFormat, 100.0));
    		int labelHeight    = metrics4Label.stringWidth(String.format("A"))*6/5;
    		int labelWidth     = metrics4Label.stringWidth(yLabel.get(0))/2;
    		
    		g2.setFont(yLabelFont);
    		
    		int nLabels        = yLabel.size();
    		int nPlots         = Utilities.getMaxFromIntList(plotInfo)+1;
    		String label       = "";
    		
    		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
        	int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
    		
    		int counter = 0;
    		 	
    		for(int r=0; r<numberOfPlotRows; r++){
    				
	    		if(counter >= nPlots){
	    			break;
	    		}
    		
	    		int yPos = (border_gap*(r+1) + rectangularHeight*r)+rectangularHeight/2;
	    		
	    		for(int c=0; c<numberOfPlotColumns; c++){
    				
	    			if(counter < nLabels){
	    				label = yLabel.get(counter);
	    			}else{
	    				label = "";
	    			}
	    			
	    			int xPos = border_gap*(c+1)+rectangularWidth*c-labelUnitWidth;//border_gap/2;
	        		        		    				
	        		AffineTransform affineTransform = new AffineTransform();
	        		affineTransform.rotate(Math.toRadians(-90), 0, 0);

	        		Font rotatedFont = yLabelFont.deriveFont(affineTransform);
	        		
	        		g2.setFont(rotatedFont);
	        		g2.drawString(label, xPos-labelHeight, yPos + (labelWidth/2));
	        			
	    			counter++;
	    		
	    			if(counter >= nPlots){
	    				break;
	    			}
	    		
	    		}
			
    		}
	    		
    		//g2.dispose();
    		g2.setFont(orgFont);
    		g2.setColor(Color.BLACK); 
    		
    	}
    	
    }
    
    
    @Override
    public Dimension getPreferredSize() {
	   
    	return new Dimension(pref_w, pref_h);
      
    }
    
      
    public static void convert_input_data(double [][] x_values, double [][] y_values, String pType){
	   
    	int [] validPlotType = Utilities.get_idx(get_plot_types(), pType);
	   
    	if(validPlotType[0] == -1){
    		throw new RuntimeException(pType + " is not a valid plot type.");
    	}
	   
    	for(int i=0; i<x_values[0].length; i++){
		   
    		List<Double> x_list = new ArrayList<Double>(x_values.length);
		   
    		for(int j=0; j<x_values.length; j++){
			   
    			x_list.add(x_values[j][i]);
			   
    		}
		   
    		x.add(x_list);
    		plotType.add(pType);
		   
    	}
	    
    	for(int i=0; i<y_values[0].length; i++){
		   
    		List<Double> y_list = new ArrayList<Double>(x_values.length);
		   
    		for(int j=0; j<y_values.length; j++){
			   
    			y_list.add(y_values[j][i]);
			   
    		}
		   
    		y.add(y_list);
		   
    	}
	   
    }
    
    
    public static void convert_input_data(String [][] x_values, double[][] yValues, String pType){
   	   
    	int [] validPlotType = Utilities.get_idx(get_plot_types(), pType);
	   
    	if(validPlotType[0] == -1){
    		throw new RuntimeException(pType + " is not a valid plot type.");
    	}
	   
    	for(int i=0; i<x_values[0].length; i++){
		   
    		List<String> x_list = new ArrayList<String>(x_values.length);
		   
    		for(int j=0; j<x_values.length; j++){
			   
    			x_list.add(x_values[j][i]);
			   
    		}
		   
    		xStr.add(x_list);
    		plotType.add(pType);
		   
    	}
	    
    	for(int i=0; i<yValues[0].length; i++){
		   
    		List<Double> y_list = new ArrayList<Double>(x_values.length);
		   
    		for(int j=0; j<yValues.length; j++){
			   
    			y_list.add(yValues[j][i]);
			   
    		}
		   
    		y.add(y_list);
		   
    	}
	   
    }
    
    
    public static String [] get_plot_types(){
	 	   
    	String [] pTypes = {"L", "P", "H", "B", "BGroup"}; 
	   
    	return pTypes;
	   
    }
    
    
	public static int [] get_idxs_4_plot_type(String pType){
	 	   
	    int [] plotTypeIdx = Utilities.get_idx(plotType, pType);
		   
		return plotTypeIdx;
		   
	}
    
	
	public static int [] get_idxs_4_plot_info(int infoNumber){
	 	   
	    int [] plotInfoIdx = Utilities.get_idx(plotInfo, infoNumber);
		   
		return plotInfoIdx;
		   
	}
	
	
	public static void remove_input_data(int [] idxs2remove){
		
		int n_idxs = idxs2remove.length;
		
		if(n_idxs != 0){
			
			int n_samples = x.size();
			
		    ArrayList<List<Double>> x_new = new ArrayList<List<Double>>();
		    ArrayList<List<Double>> y_new = new ArrayList<List<Double>>();;
		    List<String> plotType_new     = new ArrayList<String>();
			List<Integer> plotInfo_new    = new ArrayList<Integer>();
		    
		    for(int i=0; i<n_samples; i++){
		    	
		    	int [] idx2remove = Utilities.get_idx(idxs2remove, i);
		    	
		    	if(idx2remove[0] == -1){
		    		
		    		x_new.add(x.get(i));
		    		y_new.add(y.get(i));
		    		plotType_new.add(plotType.get(i));
		    		plotInfo_new.add(plotInfo.get(i));
		    		
		    	}
		    	
		    }
			
		    x = x_new;
		    y = y_new;
		    plotType = plotType_new;
		    plotInfo = plotInfo_new;
				
		}
		
	}
	
	
	public static int [] get_row_and_col_idx_4_plotPosition(int plotIdx){
		
		int [] row_col_idxs = new int [2];
		row_col_idxs[0] = -1;
		row_col_idxs[1] = -1;
		
	   	int counter = 0;
	   	int doBreak = 0;	
	   	
	   	for(int i=0; i<numberOfPlotRows; i++){
	   	    	
	   	    for(int j=0; j<numberOfPlotColumns; j++){
	   	    		
	   	    	if(counter == plotIdx){
	   	    			
	   	    		row_col_idxs[0] = i;
	   	    		row_col_idxs[1] = j;
	   	    			
	   	    		doBreak = 1;	   	    		
	   	    		break;
	   	    			
	   	    	}	
	   	    	
	   	    	counter++;
	   	    	
	   	    }
	   	    
   	    	if(doBreak==1){
    			break;
    		}
   	    			   	    
	   	}
		
	   	return row_col_idxs;
	   	    
	}
	
    
    public static double get_x_min(int plotInfoNumber){
        	   
	   	int [] infoIdxs = get_idxs_4_plot_info(plotInfoNumber);

	   	if(infoIdxs[0] == -1){
	   		throw new RuntimeException("No data found for plot no. " + plotInfoNumber);
	   	}
	   	
       int nSamples = infoIdxs.length;
	   
	   	double minValue = (x.get(infoIdxs[0])).get(0); 

	   	for(int i=0; i<nSamples;i++){ 
		    
	   		double curMin = Utilities.getMin(x.get(infoIdxs[i]));
	   		
	   		if(minValue>curMin){
	   			
	   			minValue = curMin;
	   			
	   		}
			  		      
	   	} 
		    
	   	return minValue;
          
   }
    
   
   public static double get_y_min(int plotInfoNumber){
	   
	   	int [] infoIdxs = get_idxs_4_plot_info(plotInfoNumber);

	   	if(infoIdxs[0] == -1){
	   		throw new RuntimeException("No data found for plot no. " + plotInfoNumber);
	   	}
	   	
       int nSamples = infoIdxs.length;
	  	   
	   	double minValue = (y.get(infoIdxs[0])).get(0); 

	   	for(int i=0; i<nSamples;i++){ 
		    
	   		double curMin = Utilities.getMin(y.get(infoIdxs[i]));
	   		
	   		if(minValue>curMin){
	   			
	   			minValue = curMin;
	   			
	   		}
			  		      
	   	} 
	   	
	   	return minValue;
         
   }
   
   
   public static double get_x_max(int plotInfoNumber){
	   
	   	int [] infoIdxs = get_idxs_4_plot_info(plotInfoNumber);

	   	if(infoIdxs[0] == -1){
	   		throw new RuntimeException("No data found for plot no. " + plotInfoNumber);
	   	}
	   	
        int nSamples = infoIdxs.length;
	   	
	   	double maxValue = (x.get(infoIdxs[0])).get(0); 

	   	for(int i=0; i<nSamples;i++){ 
		    
	   		double curMax = Utilities.getMax(x.get(infoIdxs[i]));
	   		
	   		if(maxValue<curMax){
	   			
	   			maxValue = curMax;
	   			
	   		}
			  		      
	   	} 
		    
	   	return maxValue;
         
   }
   
   
   public static double get_y_max(int plotInfoNumber){
	   
	   	int [] infoIdxs = get_idxs_4_plot_info(plotInfoNumber);

	   	if(infoIdxs[0] == -1){
	   		throw new RuntimeException("No data found for plot no. " + plotInfoNumber);
	   	}
	   	
       int nSamples = infoIdxs.length;

	   	double maxValue = (y.get(infoIdxs[0])).get(0); 

	   	for(int i=0; i<nSamples;i++){ 
		    
	   		double curMax = Utilities.getMax(y.get(infoIdxs[i]));
	   		
	   		if(maxValue<curMax){
	   			
	   			maxValue = curMax;
	   			
	   		}
	   			  		      
	   	} 
		    
	   	return maxValue;
        
   	}
   	
   
    public static int getGraphWidth(){
    	
    	return pref_w;
    	
    }
   
    
    public static int getGraphHeight(){
    	
    	return pref_h;
    	
    }
   
    
   	public static void setGraphWidth(int width){
   		
   		pref_w = width;
     	
   	}
   	
    
   	public static void setGraphHeight(int height){
   		
   		pref_h = height;
   		
   	}
    
   	
    public static void setTitle(String [] title4plot, String font, String size){
    	
    	for(int i=0; i<title4plot.length; i++){
    		title.add(title4plot[i]);
    	}
    	
    	int fontSize;
    	
    	if(size == null){    		
    		fontSize = 12;    		
    	}else{  
    		
    		fontSize = Integer.parseInt(size);
        	titleFont = new Font (null, Font.BOLD, fontSize);	
    	}
    	
    	if(font != null){
    		
    		if(font == "bold"){
        		titleFont = new Font (null, Font.BOLD, fontSize);
        	}
        	
        	if(font == "plain"){
        		titleFont = new Font (null, Font.PLAIN, fontSize);
        	}   		
    	}
    			
    }
    
    
    public static void setSubTitle1(String [] subTitle4plot1, String font, String size){
    	
    	for(int i=0; i<subTitle4plot1.length; i++){
    		subTitle1.add(subTitle4plot1[i]);
    	}
    	
    	int fontSize;
    	
    	if(size == null){    		
    		fontSize = 12;    		
    	}else{  
    		
    		fontSize = Integer.parseInt(size);
        	subTitle1Font = new Font (null, Font.PLAIN, fontSize);	
    	}
    	
    	if(font != null){
    		
    		if(font == "bold"){
    			subTitle1Font = new Font (null, Font.BOLD, fontSize);
        	}
        	
        	if(font == "plain"){
        		subTitle1Font = new Font (null, Font.PLAIN, fontSize);
        	}   		
    	}
    	
    }
    
    
    public static void setSubTitle2(String subTitle4plot2, String font, String size){
    	
    	subTitle2 = subTitle4plot2;
    	  
    	int fontSize;
    	
    	if(size == null){    		
    		fontSize = 12;    		
    	}else{  
    		
    		fontSize = Integer.parseInt(size);
        	subTitle2Font = new Font (null, Font.PLAIN, fontSize);	
    	}
    	
    	if(font != null){
    		
    		if(font == "bold"){
    			subTitle2Font = new Font (null, Font.BOLD, fontSize);
        	}
        	
        	if(font == "plain"){
        		subTitle2Font = new Font (null, Font.PLAIN, fontSize);
        	}   		
    	}
    	
    }
    
    
    public static void setXLabel(String [] x_label, String font, String size){
    	
    	for(int i=0; i<x_label.length; i++){
    		xLabel.add(x_label[i]);
    	}  	
    	
    	int fontSize;
    	
    	if(size == null){    		
    		fontSize = 12;    		
    	}else{  
    		
    		fontSize = Integer.parseInt(size);
        	xLabelFont = new Font (null, Font.PLAIN, fontSize);	
    	}
    	
    	if(font != null){
    		
    		if(font == "bold"){
    			xLabelFont = new Font (null, Font.BOLD, fontSize);
        	}
        	
        	if(font == "plain"){
        		xLabelFont = new Font (null, Font.PLAIN, fontSize);
        	}   		
    	}
    	
    }
    
    
    public static void setXSubLabel(String x_subLabel, String font, String size){
    	
    	xSubLabel = x_subLabel;
    	
    	int fontSize;
    	
    	if(size == null){    		
    		fontSize = 12;    		
    	}else{  
    		
    		fontSize = Integer.parseInt(size);
    		xSubLabelFont = new Font (null, Font.PLAIN, fontSize);	
    	}
    	
    	if(font != null){
    		
    		if(font == "bold"){
    			xSubLabelFont = new Font (null, Font.BOLD, fontSize);
        	}
        	
        	if(font == "plain"){
        		xSubLabelFont = new Font (null, Font.PLAIN, fontSize);
        	}   		
    	}
    	
    }
    
    
    public static void setYLabel(String [] y_label, String font, String size){
    	
    	for(int i=0; i<y_label.length; i++){
    		yLabel.add(y_label[i]);
    	}
    	
    	int fontSize;
    	
    	if(size == null){    		
    		fontSize = 12;    		
    	}else{  
    		
    		fontSize = Integer.parseInt(size);
    		yLabelFont = new Font (null, Font.PLAIN, fontSize);	
    	}
    	
    	if(font != null){
    		
    		if(font == "bold"){
    			yLabelFont = new Font (null, Font.BOLD, fontSize);
        	}
        	
        	if(font == "plain"){
        		yLabelFont = new Font (null, Font.PLAIN, fontSize);
        	}   		
    	}
    	
    }

    
    public static void setFontOfXAxisUnits(String font, int fontSize){
    	
    	if(font == "bold"){  		
    		fontOfXAxisUnits = new Font(null, Font.BOLD, fontSize);   		
    	}
    	
    	if(font == "plain"){  		
    		fontOfXAxisUnits = new Font(null, Font.PLAIN, fontSize);   		
    	}
    	
    }
    
    
    public static void setFontOfYAxisUnits(String font, int fontSize){
    	
    	if(font == "bold"){  		
    		fontOfYAxisUnits = new Font(null, Font.BOLD, fontSize);   		
    	}
    	
    	if(font == "plain"){  		
    		fontOfYAxisUnits = new Font(null, Font.PLAIN, fontSize);   		
    	}
    	
    }

    
    public static void  setNumberOfDigits4XAxis(int nDigits){
    	
    	numberOfDigitsXAxisUnits = nDigits;
    	xAxisUnitsFormat = "%."+nDigits+"f%n";
    	
    }
        
    
    public static void setColorOfXAxisUnits(Color color){
    	
    	colorOfXAxisUnits = color;
    	
    }
    
    
    public static void setColorOfYAxisUnits(Color color){
    	
    	colorOfYAxisUnits = color;
    	
    }
    
    
    public static void  setNumberOfDigits4YAxis(int nDigits){
    	
    	numberOfDigitsYAxisUnits = nDigits;
    	yAxisUnitsFormat = "%."+nDigits+"f%n";
    	
    }
    
    
    public static void setNumberOfXDivisions(int n_xDivisions){
    	
    	numberXDivisions = n_xDivisions;
    	
    }
    
    
    public static void setNumberOfYDivisions(int n_yDivisions){
    	
    	numberYDivisions = n_yDivisions;
    	
    }
    

	public static void main(String[] args) {
		g2.setColor(Color.BLACK);
	}
    
}