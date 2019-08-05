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
import java.awt.geom.GeneralPath;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import Utilities.Utilities;

//Class for generating line and point plots
@SuppressWarnings("serial")
public class GenGraphics extends GraphicDevice{

    public static List<Color> graph_line_color = new ArrayList<Color>();
    public static List<Color> graph_point_color = new ArrayList<Color>();
    public static Stroke graph_stroke = new BasicStroke(3f);
    public static int graph_point_width = 6;
   	
    public static List<List<Point>> graphPoints4Samples = new ArrayList<List<Point>>();
    
    static boolean drawErrorBars           = false;
    static List<Integer> errorBarPlotInfos = new ArrayList<Integer>();
    static Color color4ErrorBars = Color.RED;
    
    static boolean drawShadedAreas                 = false;
    static List<Integer> shadedAreasPlotInfos      = new ArrayList<Integer>();
    static ArrayList<List<Point>> ShadedAreaPoints = new ArrayList<List<Point>>();
    static List<Color> shadedAreaColor             =  new ArrayList<Color>();

    public static List<Double> extYScaleFactor = new ArrayList<Double>();
    public static List<Double> extYMaxPos = new ArrayList<Double>();
    
	@Override
	protected void paintComponent(Graphics g) {
		
	    super.paintComponent(g);
	      
	    g2 = (Graphics2D) g;     
	    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	    
	    calc_and_set_scaled_points();
	        
	    //Set white rectangular
	    set_rectangular();
	    
	    if(drawShadedAreas == true){
	    	createShadedArea();
	    }
	    
	    set_x_axis();
	    set_y_axis();
	                
	    createLinePlot();
	        
	    createPointPlot();
	        
	    if(drawErrorBars == true){	        	
	    	draw_error_bars();	        	
	    }
	        
	    set_title2plot();
	    set_xLabel2plot();
	    set_yLabel2plot();
	     
	}
	
	
	//2d array yValues = [upperValues, lowerValues]
	public static void plotShadedArea(double [][] xValues, double [][] yValues){
		
		int nValues = xValues.length;
		
		if(nValues != yValues.length){
			throw new RuntimeException("Unequal number of supplied x and y values for shaded area.");
		}
		
		drawShadedAreas = true;
		
		if(yValues[0][0] < yValues[0][1]){
			
			for(int i=0; i<nValues; i++){
				
				double yVal = yValues[i][0];
				yValues[i][0] = yValues[i][1];
				yValues[i][1] = yVal;
				
			}
			
		}
				
		plotLines(xValues, yValues, false, new Color(211, 211, 211));
		
		int plotInfoNumber = plotInfo.get(plotInfo.size()-1);
		int sampleNumber   = x.size()-1;
		
		int [] posIdxs = get_row_and_col_idx_4_plotPosition(plotInfoNumber);
		
		if(posIdxs[0] != -1){
			
			for(int i=0; i<2; i++){
				List<Point> shadedAreaPoints = calc_scaled_points_4_plot_cell(posIdxs[0], posIdxs[1], plotInfoNumber, sampleNumber-i);
				ShadedAreaPoints.add(shadedAreaPoints);
				shadedAreasPlotInfos.add(sampleNumber-i);
			}
			
			shadedAreaColor.add(new Color(211, 211, 211));
			
		}
			
	}
	
	
	public static void plotShadedArea(double [][] xValues, double [][] yValues, Color areaColor){
		
		int nValues = xValues.length;
		
		if(nValues != yValues.length){
			throw new RuntimeException("Unequal number of supplied x and y values for shaded area.");
		}
		
		drawShadedAreas = true;
		
		if(yValues[0][0]<yValues[0][1]){
			
			for(int i=0; i<nValues; i++){
				
				double yVal = yValues[i][0];
				yValues[i][0] = yValues[i][1];
				yValues[i][1] = yVal;
				
			}
			
		}
				
		plotLines(xValues, yValues, false, areaColor);
		
		int plotInfoNumber = plotInfo.get(plotInfo.size()-1);
		int sampleNumber   = x.size()-1;
		
		int [] posIdxs = get_row_and_col_idx_4_plotPosition(plotInfoNumber);
		
		if(posIdxs[0] != -1){
			
			for(int i=0; i<2; i++){
				List<Point> shadedAreaPoints = calc_scaled_points_4_plot_cell(posIdxs[0], posIdxs[1], plotInfoNumber, sampleNumber-i);
				ShadedAreaPoints.add(shadedAreaPoints);
				shadedAreasPlotInfos.add(sampleNumber-i);
			}
			
			shadedAreaColor.add(areaColor);
			
		}
			
	}
	
	
	
	public static void createShadedArea(){
		
	    int nPlots = ShadedAreaPoints.size()/2;
            	
	    if(nPlots != 0){
            
	    	int idx = 0;
	    	
	    	for(int i = 0; i<nPlots; i++){	 	        
	    		 
	 	    	GeneralPath path = new GeneralPath();

	 	    	List<Point> upperBounds = ShadedAreaPoints.get(idx);
	 	    	List<Point> lowerBounds = ShadedAreaPoints.get(idx+1);
	 	    	
	 	    	int nCoordinates = upperBounds.size();
	 	    	
	 	    	Point p = upperBounds.get(0);	 	    	
	 	    	path.moveTo(p.getX(), p.getY());
	 	    		
	 	    	for(int j=1; j<nCoordinates; j++){
	 	    		
	 	    		p = upperBounds.get(j);		 	    		
	 	    		path.lineTo(p.getX(), p.getY());
		 	    		
	 	    	}

	 	    	for(int j=0; j<nCoordinates; j++){
	 	    		
	 	    		int backwardIdx = nCoordinates-j-1;
	 	    		p = lowerBounds.get(backwardIdx);		 	    		
	 	    		path.lineTo(p.getX(), p.getY());
	 	    	}
	 	    	
	 	    	g2.setColor(shadedAreaColor.get(i));
	 		    g2.fill(path);
	 		    path.reset();
                
	 		    idx = idx+2;
	 		    
	 	    }
	    	
	    }
	    		
	}
	
	
	public static void createLinePlot(){
		
		int [] linePlotIdxs = get_idxs_4_plot_type("L");
        
	    //Stroke oldStroke = g2.getStroke();
	    if(linePlotIdxs[0] != -1){
	         	
	    	for(int i = 0; i<linePlotIdxs.length; i++){
	          	
	    		boolean isAreaLine = false;
	    		
	    		if(shadedAreasPlotInfos.size() != 0){
	    			int [] shadedAreaLine = Utilities.get_idx(shadedAreasPlotInfos, linePlotIdxs[i]);
	    			if(shadedAreaLine[0] != -1){
	    				isAreaLine = true;
	    			}
	    		}
	    		
	    		if(isAreaLine == false){
	    			
		    		g2.setColor(graph_line_color.get(i));
		            g2.setStroke(graph_stroke);
		        	           
		            List<Point> graphPoints = graphPoints4Samples.get(linePlotIdxs[i]);
		        	int sampleLength        = graphPoints.size();
		            	
		        	for(int j=0; j<sampleLength-1; j++){
		        		 
		            	int x1 = graphPoints.get(j).x;
		            	int y1 = graphPoints.get(j).y;
		            	int x2 = graphPoints.get(j+1).x;
		            	int y2 = graphPoints.get(j+1).y;
		            	g2.drawLine(x1, y1, x2, y2); 
		        		 
		            }
	    			
	    		}	    	
	        	   
	        }
	        	
	    }
		
	}
	
	
	public static void createPointPlot(){
		
	    int [] pointPlotIdxs = get_idxs_4_plot_type("P");
	    //g2.setStroke(oldStroke);   
	        
	    if(pointPlotIdxs[0] != -1){        	
	        	
	        for(int i=0; i<pointPlotIdxs.length; i++){
	          	
	        	g2.setColor(graph_point_color.get(i));
	        	
	            List<Point> graphPoints = graphPoints4Samples.get(pointPlotIdxs[i]);
	            int sampleLength        = graphPoints.size();
	            	
	            for(int j=0; j<sampleLength; j++) {
	                        
	            	int x = graphPoints.get(j).x - graph_point_width / 2;
	            	int y = graphPoints.get(j).y - graph_point_width / 2;;
	            	int ovalW = graph_point_width;
	            	int ovalH = graph_point_width;
	            	g2.fillOval(x, y, ovalW, ovalH);
	        	         
	            }
	        	  
	        }
	        	
	    }
		
	}
	 
	
	public static void draw_error_bars(){
	    		    		
	    int nPlots = Utilities.getMaxFromIntList(plotInfo)+1;
	    	
	    for(int i = 0; i<nPlots; i++){
	        
	    	int [] plotIdxs = get_idxs_4_plot_info(i);
	    	
	    	int [] errorBarIdxs = Utilities.get_idx(errorBarPlotInfos, i);
	    	
	    	if(errorBarIdxs[0] != -1){
	    		if(plotIdxs.length >= 2){
		    			
		    		if(plotType.get(plotIdxs[0]) != plotType.get(plotIdxs[1])){
		    			
		    			int pointIdx = 0;
		    			if(plotType.get(plotIdxs[0]) != "P"){
		    				pointIdx = 1;
		    			}
		    			
			  		    g2.setStroke(graph_stroke);
			  		    				  		    
			  		    List<Point> pointPlotsSample = graphPoints4Samples.get(plotIdxs[pointIdx]);
			  		    Color pointColor = graph_point_color.get(pointIdx);
			  		    
			  		    List<Point> GraphPoints1 = graphPoints4Samples.get(plotIdxs[0]);
			  		    List<Point> GraphPoints2 = graphPoints4Samples.get(plotIdxs[1]);
			  		        	
			  		    if(GraphPoints1.size() != GraphPoints2.size()){
			  		        throw new RuntimeException("Supplied data for error bar plot have unequal length."); 
			  		    }
			  		        	
			  		    for(int j=0; j<GraphPoints1.size(); j++){
			  		    		 
			  		    	g2.setColor(color4ErrorBars);
			  		        int x1 = GraphPoints1.get(j).x;
			  		        int y1 = GraphPoints1.get(j).y;
			  		        int x2 = GraphPoints2.get(j).x;
			  		        int y2 = GraphPoints2.get(j).y;
			  		        g2.drawLine(x1, y1, x2, y2); 
			  		    		
			  		        //
			  		        g2.setColor(pointColor);
			  		        int x = pointPlotsSample.get(j).x - graph_point_width / 2;
			  		        int y = pointPlotsSample.get(j).y - graph_point_width / 2;;
			  		        int ovalW = graph_point_width;
			  		        int ovalH = graph_point_width;
			  		        g2.fillOval(x, y, ovalW, ovalH);
			  		        
			  		    } 
		    			
		    		}    		
		    		
		    	}
	    		
	    	}
	    		     
	    }
	    		       	   	
	}
	
	
	public void calc_and_set_scaled_points(){
	    
		int nSamples = x.size();
		int counter  = 0;
		
		List<String> plotTypes = new ArrayList<String>(numberOfPlotRows*numberOfPlotColumns);
		
		int r = 0;
				
		while(r<numberOfPlotRows){
			
			if(counter>=nSamples){	    		
	    		break;	    		
	    	}
			
			int plotInfoNumber = plotInfo.get(counter);
			
			int c = 0;
			
			while(c<numberOfPlotColumns){
					
				plotTypes.add(plotType.get(counter));
		    
			    List<Point> graphPoints = new ArrayList<Point>();
  
			    graphPoints = calc_scaled_points_4_plot_cell(r, c, plotInfoNumber, counter);
		    
			    graphPoints4Samples.add(graphPoints);
			      			
				counter++;
				
				if(counter>=nSamples){	    		
		    		break;	    		
		    	}
				
				plotInfoNumber = plotInfo.get(counter);
				
				if(counter != 0){
					
					if(plotInfoNumber != plotInfo.get(counter-1)){
						c++;
					}
					
				}
				
			}
			
			if(counter != 0){
				
				if(plotInfoNumber != plotInfo.get(counter-1)){
					r++;
				}
				
			}
			
		}
		  
		plotType = plotTypes;
		
	}
	
	
	public static List<Point> calc_scaled_points_4_plot_cell(int plotRow, int plotColumn, int plotInfoNumber, int sampleNumber){
		
		int rectangularWidth = (int) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((pref_h-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
		double x_min = get_x_min(plotInfoNumber);
		double x_max = get_x_max(plotInfoNumber);
		double y_min = get_y_min(plotInfoNumber);
		double y_max = get_y_max(plotInfoNumber);
		
	    int sample_length = x.get(sampleNumber).size();
	      
	    double xScale = ((double) rectangularWidth)/(x_max-x_min); //((double) pref_w-2*border_gap)/(x_max-x_min);
	    
	    double yScaleFactor = y_max-y_min;
	        
	    double yMaxPos = border_gap*(plotRow+1)+rectangularHeight*plotRow;	    
	    
        double yScale = (double) (((border_gap+rectangularHeight)*(plotRow+1))-yMaxPos)/yScaleFactor; //(double) ((pref_h-border_gap) - yMaxPos)/yScaleFactor;
	    
	    List<Point> graphPoints = new ArrayList<Point>();
	    	  
	    for (int i=0; i<sample_length; i++) {
	                 
	        double xVal = x.get(sampleNumber).get(i);
	        double yVal = y.get(sampleNumber).get(i);
	        	 
	        int x1 = (int) ((xVal-x_min)*xScale + border_gap*(plotColumn+1)+rectangularWidth*plotColumn);  //(int) ((xVal-x_min)*xScale + border_gap);       			
            int y1 = (int) ((border_gap+rectangularHeight)*(plotRow+1)+(y_min-yVal)*yScale); //((pref_h-border_gap) + (y_min-yVal)*yScale);  (y_min-yVal)*yScale)
	        		        	
	        graphPoints.add(new Point(x1, y1));
	             
	    }
		
	    return graphPoints;
	    
	}
	
	
	public static List<Point> calc_scaled_points_4_plot_cell(int plotRow, int plotColumn, int plotInfoNumber, int sampleNumber, double extYScaleFactor, double extYMaxPos){
		
		int rectangularWidth = (int) ((pref_w-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((pref_h-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
		double x_min = get_x_min(plotInfoNumber);
		double x_max = get_x_max(plotInfoNumber);
		double y_min = get_y_min(plotInfoNumber);
		//double y_max = get_y_max(plotInfoNumber);
		
	    int sample_length = x.get(sampleNumber).size();
	      
	    double xScale = ((double) rectangularWidth)/(x_max-x_min); //((double) pref_w-2*border_gap)/(x_max-x_min);
	    
	    double yScaleFactor = extYScaleFactor;	    	
	    double yMaxPos = extYMaxPos; //border_gap;	    	
 
        double yScale = (double) (((border_gap+rectangularHeight)*(plotRow+1))-yMaxPos)/yScaleFactor; //(double) ((pref_h-border_gap) - yMaxPos)/yScaleFactor;
	    
	    List<Point> graphPoints = new ArrayList<Point>();
	    	  
	    for (int i=0; i<sample_length; i++) {
	                 
	        double xVal = x.get(sampleNumber).get(i);
	        double yVal = y.get(sampleNumber).get(i);
	        	 
	        int x1 = (int) ((xVal-x_min)*xScale + border_gap*(plotColumn+1)+rectangularWidth*plotColumn);  //(int) ((xVal-x_min)*xScale + border_gap);       			
            int y1 = (int) ((border_gap+rectangularHeight)*(plotRow+1)+(y_min-yVal)*yScale); //((pref_h-border_gap) + (y_min-yVal)*yScale);  (y_min-yVal)*yScale)
	        		        	
	        graphPoints.add(new Point(x1, y1));
	             
	    }
		
	    return graphPoints;
	    
	}
	
	
	public void set_yScaleFactor(double scaleFactor){
		
		extYScaleFactor.add(scaleFactor);
		
	}
	
	
	public void set_yMaxPos(double y_maxPos){
		
		extYMaxPos.add(y_maxPos);
		
	}
	 
	
	public static void setDefaultDesign(){
	    	
		setNumberOfXDivisions(10);
		setNumberOfYDivisions(5);
	    	
		setTextColor(new Color(128,128,128));
		    
		setFontOfXAxisUnits("bold",12);
		setColorOfXAxisUnits(new Color(128,128,128));
		setFontOfYAxisUnits("bold",12);
		setColorOfYAxisUnits(new Color(128,128,128));
	    	
	}
	
	
	public void set_point_width(int point_width){
		graph_point_width = point_width;
	}
	
	
	public void set_x_axis(){
        
		Font orgFont = g2.getFont();
		g2.setFont(fontOfXAxisUnits);
		
		int nSamples = x.size();
		
		int counter = 0;
			
		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
	    int plotInfoNumber = -1;	    	    
	    double x_min = 0.0;
	    double x_max = 0.0;
	    int nPlotInfoIdxs = 0;
	    int r = 0;
	    
		while(r<numberOfPlotRows){
			
			if(counter>=nSamples){    	    		
	    		break;    	    		
	    	}
					
			int yUp     = border_gap*(r+1)+rectangularHeight*r;
			int yBottom = (border_gap+rectangularHeight)*(r+1);
			
			int c = 0;
			
			while(c<numberOfPlotColumns){
				
				int curPlotInfoNumber = plotInfo.get(counter);
				if(curPlotInfoNumber != plotInfoNumber){
					plotInfoNumber = curPlotInfoNumber;
					
					int [] plotInfoIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
					nPlotInfoIdxs = plotInfoIdxs.length;
					ArrayList<List<Double>> xCollection = new ArrayList<List<Double>>(nPlotInfoIdxs);
					
					for(int i=0; i<nPlotInfoIdxs; i++){
						xCollection.add(x.get(plotInfoIdxs[i]));
					}
					
					x_min = Utilities.getMinFromDblList(xCollection);
					x_max = Utilities.getMaxFromDblList(xCollection);
					
				}
				
				int xLeft   = border_gap*(c+1)+rectangularWidth*(c);
				int xRight  = (border_gap+rectangularWidth)*(c+1);
				
				g2.setColor(Color.BLACK);
				g2.drawLine(xLeft, yBottom, xRight, yBottom);

			    int step_size  = (rectangularWidth)/(numberXDivisions);
			    double width   = 0.05*border_gap;
				   
			    double l = (x_max-x_min)/(numberXDivisions);
				    
			    for(int i = 0; i<numberXDivisions+1; i++) {
				       
			    	int x0 = step_size*(i)+xLeft;
			    	int x1 = x0;
			    	int y0 = yBottom;
			    	int y1 = y0;
			    	
			    	g2.setColor(Color.BLACK);
			    	g2.drawLine(x0, y0, x1, y0+(int)width);
				       
			    	String xlabel = String.format(xAxisUnitsFormat, x_min+l*(i))+ "";
			    	FontMetrics metrics = g2.getFontMetrics();
			    	int labelWidth = metrics.stringWidth(xlabel);
			    	
			    	g2.setColor(textColor);
			    	g2.drawString(xlabel, x0-labelWidth/2, y0+metrics.getHeight() + 3);
				       
			    	if(grid == true){				                 
			    		g2.setColor(grid_color);
			    		g2.drawLine(x0, y1, x0, yUp); 				    		   	       	       	   
			    	}
			    	
			    	g2.setColor(Color.BLACK);
			    	
			    }
				
				counter=counter+nPlotInfoIdxs;
				
				if(counter>=nSamples){    	    		
		    		break;    	    		
		    	}
				
				if(counter != 0){					
					if(plotInfo.get(counter)!=plotInfo.get(counter-1)){
						c++;
					}					
				}
				
			}
			
			if(counter != 0){	
				
				if(counter>=nSamples){    	    		
		    		break;    	    		
		    	}
				
				if(plotInfo.get(counter)!=plotInfo.get(counter-1)){
					r++;
				}				
			}
			
		}
		 
		g2.setFont(orgFont);
	    g2.setColor(Color.BLACK);
	    
	}
	
	   
	public void set_y_axis(){

		Font orgFont = g2.getFont();
	    g2.setFont(fontOfYAxisUnits);
		
		int nSamples = x.size();
		
		int counter = 0;
			
		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
	    int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
	    int plotInfoNumber = -1;	    
	    double y_min = 0.0;
	    double y_max = 0.0;
	    int nPlotInfoIdxs = 0;
	    int r = 0;
	    
		while(r<numberOfPlotRows){
			
			if(counter>=nSamples){	    		
	    		break;	    		
	    	}
					
			int yUp     = border_gap*(r+1)+rectangularHeight*r;
			int yBottom = (border_gap+rectangularHeight)*(r+1);
			
			int c = 0;
			
			while(c<numberOfPlotColumns){
				
				int xLeft   = border_gap*(c+1)+rectangularWidth*(c);
				int xRight  = (border_gap+rectangularWidth)*(c+1);
				
			    g2.setColor(Color.BLACK);
			    g2.drawLine(xLeft, yBottom, xLeft, yUp);
				
			    int step_size     = rectangularHeight/(numberYDivisions); 
			    double width      = 0.05*border_gap;
				
				int curPlotInfoNumber = plotInfo.get(counter);
				
				if(curPlotInfoNumber != plotInfoNumber){
					plotInfoNumber = curPlotInfoNumber;
					int [] plotInfoIdxs = Utilities.get_idx(plotInfo, plotInfoNumber);
					nPlotInfoIdxs = plotInfoIdxs.length;
					
					y_min = get_y_min(plotInfoNumber);
					y_max = get_y_max(plotInfoNumber);
					
				}
			    			    				   
			    double l = (y_max-y_min)/(numberYDivisions);
				   
			    for(int i = 0; i<numberYDivisions+1; i++) {
				       
			    	int x0 = xLeft;
			    	int x1 = x0-(int) width;
			    	int y0 = (step_size*i)+yUp;
			    	int y1 = y0;
			    	
			    	g2.setColor(Color.BLACK);
			    	g2.drawLine(x0, y0, x1, y1);
				       			    	 
			    	String ylabel = String.format(yAxisUnitsFormat, (y_max-l*i))+ "";
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
				
				counter = counter+nPlotInfoIdxs;
				
				if(counter>=nSamples){    	    		
		    		break;    	    		
		    	}
			    				
				if(counter != 0){					
					if(plotInfo.get(counter)!=plotInfo.get(counter-1)){
						c++;
					}					
				}
				
			}
			
			if(counter != 0){
				
				if(counter>=nSamples){    	    		
		    		break;    	    		
		    	}
				
				if(plotInfo.get(counter)!=plotInfo.get(counter-1)){
					r++;
				}
				
			}
			
		}
		  
		g2.setFont(orgFont);
		g2.setColor(Color.BLACK);
		
	}
	

	public static void setLineColor(List<Color> lineColor){
		
		if(graph_line_color.size() != 0){
			for(int i=0; i<graph_line_color.size(); i++){
				
				if(i<lineColor.size()){
					graph_line_color.set(i, lineColor.get(i));
				}	
			}
		}
			
	}
	
	
	public static void setColor4ErrorBars(Color errorBarColor){
		color4ErrorBars = errorBarColor;
	}
	
	   
	public static void plotLines(double [][] x_values, double [][] y_values, boolean newPlot){
	 	   
		convert_input_data(x_values, y_values, "L");
		
	    int nSamples = x_values[0].length;
	    
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
			
			graph_line_color.add(new Color(44, 102, 230, 180));
			
		}
		
		setDefaultDesign();	
		
	}
	   
	
	public static void plotLines(double [][] x_values, double [][] y_values, boolean newPlot, Color lineColor){
	 	   
		convert_input_data(x_values, y_values, "L");
		
	    int nSamples = x_values[0].length;
	    
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
			
			graph_line_color.add(lineColor);
			
		}
		
		setDefaultDesign();	
		
	}
	
	   	
	public static void plotPoints(double [][] x_values, double [][] y_values, boolean newPlot){
		   
	    convert_input_data(x_values, y_values, "P");
		
	    int nSamples = x_values[0].length;
	    
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
			
			graph_point_color.add(new Color(100, 100, 100, 180));
			
		}
		
		setDefaultDesign();
		
	}
	  
	
	public static void plotPoints(double [][] x_values, double [][] y_values, boolean newPlot, Color pointColor){
		   
	    convert_input_data(x_values, y_values, "P");
		
	    int nSamples = x_values[0].length;
	    
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
			
			graph_point_color.add(pointColor);
			
		}
				
		setDefaultDesign();
		
	}
	
	    
    public static void drawErrorBars(boolean errorBars){
	    	
	    drawErrorBars = errorBars;
	    
	    errorBarPlotInfos.add(plotInfo.get(plotInfo.size()-1));
	    
	}
	
    
    private static void createAndShowGui() {

    	GenGraphics mainPanel = new GenGraphics();

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
    
   	
   	public static void parametrizationExample1(){
   		
   		Random r = new Random();
		
		int maxDataPoints = 60;
	     	        	
	 	double [][] x2 = new double [maxDataPoints][1];
	 	double [][] y2 = new double [maxDataPoints][1];
	 	 
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x2[i][0] = i;
	 		y2[i][0] = i-20;
	 	}
	 	
	 	double [][] xBand = new double [maxDataPoints][2];
	 	double [][] yBand = new double [maxDataPoints][2];
	 	 
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		xBand[i][0] = i;
	 		xBand[i][1] = i;
	 		yBand[i][0] = i-10;
	 		yBand[i][1] = i-30;
	 	}
	 	
	 	
	 	double [][] x3 = new double [maxDataPoints][1];
	 	double [][] y3 = new double [maxDataPoints][1];
	 		   
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		x3[i][0] = i;
	 		y3[i][0] = i-20+20.0*r.nextGaussian();
	 	}
	 	  
	 	
	 	double [][] xBand2 = new double [maxDataPoints][2];
	 	double [][] yBand2 = new double [maxDataPoints][2];
	 	 
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 		xBand2[i][0] = i;
	 		xBand2[i][1] = i;
	 		yBand2[i][0] = y3[i][0]+1;
	 		yBand2[i][1] = Utilities.getMin(y3,0);
	 	}
	 	
	 	setNumberOfPlotColums(3);
	 	setNumberOfPlotRows(2);
	 	
	 	setGraphWidth(1000);
	 	setGraphHeight(600);

	 	plotLines(x3,y3,true);
	 	plotShadedArea(xBand2, yBand2, Color.GRAY);
	 	
	 	plotLines(x2,y2,true);
	 	plotPoints(x3,y3,false, Color.RED);
	 	plotShadedArea(xBand, yBand);
	 	
	    plotLines(x2,y2,true);
	    plotPoints(x3,y3,false);
	    drawErrorBars(true);
	    
	 	plotLines(x2,y2,true);
	 	plotShadedArea(xBand, yBand, Color.GRAY);
	 	
	 	plotLines(x3,y3,true);
	 	plotShadedArea(xBand2, yBand2, Color.GRAY);
	 	
	 	plotLines(x2,y2,true,Color.BLACK);	
	    plotPoints(x3,y3,false,Color.BLACK);
	 	drawErrorBars(true);
	 	
	 	
	 	String [] title    = {"Title1", "Title2", "Title3", "Title4", "Title5", "Title6"};
	 	String [] subTitle = {"SubTitle1", "SubTitle2", "SubTitle3", "SubTitle4", "SubTitle5", "SubTitle6"};
	 	String [] yLabel   = {"y", "y2", "y3", "y4", "y5", "y6"};
	 	String [] xLabel   = {"x1", "x2", "x3", "x4", "x5", "x6"};
	 	
 	 	setTitle(title, null, "12");
	 	setSubTitle1(subTitle, null, "10");
	 	//setSubTitle2("SubTitle2", "bold", "10");
	 	//setXLabel("xLabel", "bold", null);
	 	//setXSubLabel("xSubLabel", null, null);
	 	setYLabel(yLabel, null, "10");
	 	setXLabel(xLabel, null, "10");
	 	setNumberOfDigits4XAxis(0);   
	 	setNumberOfDigits4YAxis(2);
	 	setFontOfXAxisUnits("plain", 10);
	 	setFontOfYAxisUnits("plain", 10);
	 	
	 	plot();
   		
   	}
   	
   	
   	public static void parametrizationExample2(){
   		
		int maxDataPoints = 61;
	      
	 	double [][] x_values = new double [maxDataPoints][2];
	 	double [][] y_values = new double [maxDataPoints][2];
	 	      
	 	for(int i = 0; i < maxDataPoints ; i++) {
	 	   x_values[i][0] = i;
	 	   y_values[i][0] = i;
	 	   x_values[i][1] = i;
	 	   y_values[i][1] = i-60;
	 	}
	 	      		
	 	String [] title    = {"Title1", "Title2"};
	 	String [] subTitle = {"SubTitle1", "SubTitle2"};
	 	String [] yLabel   = {"y", "y2"};
	 	String [] xLabel   = {"x1", "x2"};
	 	
	 	List<Color> lineColor = new ArrayList<Color>();
	 	lineColor.add(Color.RED);
	 	lineColor.add(Color.BLACK);
	 	
	 	setNumberOfPlotColums(2);
	 	setNumberOfPlotRows(1);
	 	
	 	setGraphWidth(900);
	 	setGraphHeight(500);
	 	
	 	plotLines(x_values,y_values, true);
	 	
	 	setLineColor(lineColor);	 	
 	 	setTitle(title, null, "12");
	 	setSubTitle1(subTitle, null, "10");
	 	setSubTitle2("SubTitle2", "bold", "10");
	 	setYLabel(yLabel, null, "10");
	 	setXLabel(xLabel, null, "10");
	 	setNumberOfDigits4XAxis(0);   
	 	setNumberOfDigits4YAxis(2);
	 	setFontOfXAxisUnits("plain", 10);
	 	setFontOfYAxisUnits("plain", 10);
	 	
	 	plot();
	 	
   	}
   	
   	
   	public static void parametrizationExample3(){
   		
		int maxDataPoints1 = 60;
	      
	 	double [][] x1_values = new double [maxDataPoints1][1];
	 	double [][] y1_values = new double [maxDataPoints1][1];
	 	      
	 	for(int i=0; i<maxDataPoints1; i++) {
	 		x1_values[i][0] = i+10;
	 	    y1_values[i][0] = i;
	 	}
	 	  
	 	int maxDataPoints2 = 60;
	 	
	 	double [][] x2_values = new double [maxDataPoints2][1];
	 	double [][] y2_values = new double [maxDataPoints2][1];
	 	
	 	for(int i=0; i<maxDataPoints2; i++){
	 		x2_values[i][0] = i+1;
	 	    y2_values[i][0] = i-20;
	 	}
	 	
	 	
	 	int maxDataPoints3 = 160;
	 	
	 	double [][] x3_values = new double [maxDataPoints3][1];
	 	double [][] y3_values = new double [maxDataPoints3][1];
	 	double [][] x3_band   = new double [maxDataPoints3][2];
	 	double [][] y3_band   = new double [maxDataPoints3][2];
	 	
	 	for(int i=0; i<maxDataPoints3; i++){
	 		x3_values[i][0] = -(i+10);
	 		x3_band[i][0]   = -(i+10);
	 		x3_band[i][1]   = -(i+10);
	 	    y3_values[i][0] = Math.sqrt(i);
	 	    y3_band[i][0]   = y3_values[i][0]-5.0;
	 	    y3_band[i][1]   = y3_values[i][0]+5.0;	 		
	 	}
	 	
	 	String [] title    = {"Title1", "Title2"};
	 	String [] subTitle = {"SubTitle1", "SubTitle2"};
	 	String [] yLabel   = {"y", "y2"};
	 	String [] xLabel   = {"x1", "x2"};
	 	
	 	List<Color> lineColor = new ArrayList<Color>();
	 	lineColor.add(Color.RED);
	 	lineColor.add(Color.BLACK);
	 	
	 	setNumberOfPlotColums(2);
	 	setNumberOfPlotRows(2);
	 	
	 	setGraphWidth(900);
	 	setGraphHeight(500);
	 	
	 	plotLines(x1_values,y1_values, true);
	 	plotLines(x2_values,y2_values,false);
	 	plotLines(x3_values,y3_values,false);
	 	plotShadedArea(x3_band, y3_band, Color.GRAY);
	 	plotLines(x1_values,y1_values,true);
	 	plotLines(x2_values,y2_values,false);
	 	plotPoints(x1_values,y1_values,true);
	 	plotLines(x2_values,y2_values,true);
	 	plotPoints(x2_values,y2_values,false);
	 	plotLines(x2_values,y2_values,false,Color.BLACK);
	 	
	 	setLineColor(lineColor);	 	
 	 	setTitle(title, null, "12");
	 	setSubTitle1(subTitle, null, "10");
	 	//setSubTitle2("SubTitle2", "bold", "10");
	 	setYLabel(yLabel, null, "10");
	 	setXLabel(xLabel, null, "10");
	 	setNumberOfDigits4XAxis(0);   
	 	setNumberOfDigits4YAxis(2);
	 	setFontOfXAxisUnits("plain", 10);
	 	setFontOfYAxisUnits("plain", 10);
	 	
	 	plot();
	 	
   	}
   	
   	
   	public static void parametrizationExample4(){
   		
   		Random r = new Random();
	
		int maxDataPoints1 = 60;
	      
	 	double [][] x1_values = new double [maxDataPoints1][1];
	 	double [][] y1_values = new double [maxDataPoints1][1];
	 	      
	 	for(int i=0; i<maxDataPoints1; i++) {
	 		x1_values[i][0] = i+1;	 		
	 	    y1_values[i][0] = i+5.0*r.nextGaussian();
	 	}
	 	 
	 	double [][] x2_values = new double [maxDataPoints1][1];
	 	double [][] y2_values = new double [maxDataPoints1][1];
	 	      
	 	for(int i=0; i<maxDataPoints1; i++) {
	 		x2_values[i][0] = i+1;	 		
	 	    y2_values[i][0] = i+15.0*r.nextGaussian();
	 	}
	 	
	 	double [] y_max = new double [2];
	 	double [] y_min = new double [2];
	 	y_max[0] = Utilities.getMax(y1_values);
	 	y_max[1] = Utilities.getMax(y2_values);
	 	y_min[0] = Utilities.getMin(y1_values);
	 	y_min[1] = Utilities.getMin(y2_values);
	 	
	 	
	 	double [][] xrecession1 = new double [4][2];
	 	double [][] yrecession1 = new double [4][2];
	 	
	 	double [][] xrecession2 = new double [2][2];
	 	double [][] yrecession2 = new double [2][2];
	 	
	 	double [][] xrecession3 = new double [8][2];
	 	double [][] yrecession3 = new double [8][2];
	 	
	 	double [][] xrecession4 = new double [3][2];
	 	double [][] yrecession4 = new double [3][2];
	 	
	 	double [][] xrecession5 = new double [5][2];
	 	double [][] yrecession5 = new double [5][2];
	 	
	 	//-----------------------
	 	//Number of Plots and width/height has to be defined before plotLines/plotPoints is used! 
	 	setNumberOfPlotColums(2);
	 	setNumberOfPlotRows(1);	 	
	 	setGraphWidth(700);
	 	setGraphHeight(500);
	 	//-----------------------
	 	
	 	for(int p=0; p<2; p++){
	 			 	
		 	for(int i=0; i<4; i++){	 		
		 		xrecession1[i][0] = i+4;
		 		xrecession1[i][1] = i+4;
		 		yrecession1[i][0] = y_min[p];
		 		yrecession1[i][1] = y_max[p];	 		
		 	}
		 		
		 	for(int i=0; i<2; i++){	 		
		 		xrecession2[i][0] = i+11;
		 		xrecession2[i][1] = i+11;
		 		yrecession2[i][0] = y_min[p];
		 		yrecession2[i][1] = y_max[p];	 		
		 	}
		 	
		 	for(int i=0; i<8; i++){	 		
		 		xrecession3[i][0] = i+21;
		 		xrecession3[i][1] = i+21;
		 		yrecession3[i][0] = y_min[p];
		 		yrecession3[i][1] = y_max[p];	 		
		 	}
		 	
		 	for(int i=0; i<3; i++){	 		
		 		xrecession4[i][0] = i+41;
		 		xrecession4[i][1] = i+41;
		 		yrecession4[i][0] = y_min[p];
		 		yrecession4[i][1] = y_max[p];	 		
		 	}
		 		
		 	for(int i=0; i<5; i++){	 		
		 		xrecession5[i][0] = i+56;
		 		xrecession5[i][1] = i+56;
		 		yrecession5[i][0] = y_min[p];
		 		yrecession5[i][1] = y_max[p];	 		
		 	}
	 		
		 	if(p==0){
		 		plotLines(x1_values,y1_values, true, Color.BLACK);
		 	}
		 	if(p==1){
		 		plotLines(x2_values,y2_values, true, Color.RED);
		 	}
		 	
		 	plotShadedArea(xrecession1, yrecession1, Color.GRAY);
		 	plotShadedArea(xrecession2, yrecession2, Color.GRAY);
		 	plotShadedArea(xrecession3, yrecession3, Color.GRAY);
		 	plotShadedArea(xrecession4, yrecession4, Color.GRAY);
		 	plotShadedArea(xrecession5, yrecession5, Color.GRAY);
		 	
	 	}
	 	
	 	String [] title    = {"Real GDP", "Inflation"};
	 	String [] subTitle = {"(with Recessions)", "SubTitle2"};
	 	String [] yLabel   = {"Real GDP (in %)", "Inflation"};
	 	String [] xLabel   = {"Dates", "Dates"};
	 		 	
 	 	setTitle(title, null, "14");
	 	setSubTitle1(subTitle, null, "12");
	 	//setSubTitle2("SubTitle2", "bold", "10");
	 	setYLabel(yLabel, "bold", "12");
	 	setXLabel(xLabel, "bold", "12");
	 	setNumberOfDigits4XAxis(0);   
	 	setNumberOfDigits4YAxis(1);
	 	setFontOfXAxisUnits("bold", 12);
	 	setFontOfYAxisUnits("bold", 12);
	 	
	 	plot();
	 	
   	}
   	
   		
	public static void main(String[] args) {
		
		parametrizationExample1();
		//parametrizationExample2();
		//parametrizationExample3();
		//parametrizationExample4();
		
	}
		
}
