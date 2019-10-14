package Graphics;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import AdaptiveBasisModels.CART;
import Regression.LinearRegression;
import Utilities.Utilities;

@SuppressWarnings("serial")
public class DecisionTreeGraphics extends GraphicDevice{

	public static Stroke graph_stroke = new BasicStroke(3f);
    public static int graph_point_width = 10;
	public static Color pointColor = Color.BLUE;
	public static Color leafColor  = Color.green;
	public static Color lineColor = Color.BLACK;
	
	public static boolean showLeafs = true;
	
	public static Font infoFont  = new Font(null, Font.PLAIN, 9);
	public static int decimalPlaces4InfoBox = 2;
	
	public static boolean plotLayerGrid = false;
    
	public static boolean rotateTree = false;
	
	public static CART cart_obj;
	
	HashMap<String, HashMap<String, Integer>> internal_knot_positions;

	
	@Override
	protected void paintComponent(Graphics g) {
		
	    super.paintComponent(g);
	    
	    if(cart_obj == null) {
	    	throw new RuntimeException("No CART object supplied for plotting the tree.");
	    }
	    
	    if(isTreeSetted() == false) {
	    	throw new RuntimeException("No tree found in supplied CART object.");
	    }
	    
	    g2 = (Graphics2D) g;     
	    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	    
	    set_rectangular();
	   	    
	    calc_and_set_tree_knots();
    
	    create_and_set_info_box(2,2);
	    
	    set_title2plot();
	    	
	}
	
	
	@SuppressWarnings("static-access")
	public void calc_and_set_tree_knots() {
		
		g2.setColor(pointColor);
		g2.setStroke(graph_stroke);
		
		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
		int rectangularHeight = (int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
		HashMap <String, HashMap<String, HashMap<String,List<String>>>> tree = cart_obj.get_tree();
		
		HashMap<String, List<String>> leafs4Layers = null;	
		if(showLeafs == true) {
			leafs4Layers = cart_obj.get_leafs_4_layers();
		}
		
		int treeDepth = tree.size();				
		get_internal_knot_positions();
		
		int maxNumberOfKnots = get_maxNumberOfKnotsInLayer();
		maxNumberOfKnots = 2*maxNumberOfKnots+1;
		
		int distBetweenLayers = 0;
		int distBetweenKnots  = 0;
		
		if(rotateTree == false) {
			distBetweenLayers = (int) rectangularHeight/(treeDepth-1);
		}
		if(rotateTree == true) {
			distBetweenLayers = (int) rectangularWidth/(treeDepth-1);
		}

		ArrayList<List<Point>> points= new ArrayList<List<Point>>();
			
		for(int l=0; l<treeDepth; l++) {
			
			String layerLabel = "Layer"+(l+1);
			int nKnots = tree.get(layerLabel).size();											
			int x = 0;
			int y = 0;
			int startY4LayerLine = 0;
			int endY4LayerLine   = 0;
			int startX4LayerLine = 0;
			int endX4LayerLine   = 0;
			
			if(rotateTree == false) {
				distBetweenKnots = (int) rectangularWidth/(maxNumberOfKnots+1);
				x = border_gap;
				y = border_gap+distBetweenLayers*l-graph_point_width/2;
				startX4LayerLine = x;
				endX4LayerLine   = x+rectangularWidth;
				startY4LayerLine = y+graph_point_width/2;
				endY4LayerLine   = startY4LayerLine;
			}
			
			if(rotateTree == true) {
				distBetweenKnots = (int) rectangularHeight/(maxNumberOfKnots+1);
				y = border_gap;
				x = border_gap+distBetweenLayers*l-graph_point_width/2;
				startX4LayerLine = x+graph_point_width/2;
				endX4LayerLine   = startX4LayerLine;
				startY4LayerLine = y;
				endY4LayerLine   = y+rectangularHeight; 
			}
			
			HashMap<String,Integer> int_knot_pos_4_layer = internal_knot_positions.get(layerLabel);
			List<Point> layerKnotPoints = new ArrayList<Point>();
						
			if(plotLayerGrid == true) {
				if(l!=0 && l!=(treeDepth-1)) {
					Stroke dashed = new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{9}, 0);				
					g2.setStroke(dashed);
					g2.setColor(backroundColor);
					g2.drawLine(startX4LayerLine, startY4LayerLine, endX4LayerLine, endY4LayerLine);
					g2.setStroke(graph_stroke);
				}
			}
								
			for(int k=0; k<nKnots; k++) {		
				
				String knotLabel = "Knot"+(k+1);
				int internalPos = int_knot_pos_4_layer.get(knotLabel);
				
				if(rotateTree == false) {
					x = border_gap+distBetweenKnots*internalPos-graph_point_width/2;
				}
				if(rotateTree == true) {
					y = border_gap+distBetweenKnots*internalPos-graph_point_width/2;
				}
				
				layerKnotPoints.add(new Point(x,y));
				
				int parentKnotNumber = Integer.parseInt(tree.get(layerLabel).get(knotLabel).get("ParentKnot").get(0).substring(4));
							
				if(l!=0) {
						
					int x1 = points.get((l-1)).get(parentKnotNumber-1).x+graph_point_width/2;
					int y1 = points.get((l-1)).get(parentKnotNumber-1).y+graph_point_width/2;
					int x2 = layerKnotPoints.get(k).x+graph_point_width/2;	
					int y2 = layerKnotPoints.get(k).y+graph_point_width/2;
					g2.setColor(lineColor);
					g2.drawLine(x1, y1, x2, y2);
					g2.setStroke(graph_stroke);
					g2.setColor(pointColor);
					g2.fillOval(x1-graph_point_width/2, y1-graph_point_width/2, graph_point_width, graph_point_width);
					g2.setColor(pointColor);
					g2.fillOval(x2-graph_point_width/2, y2-graph_point_width/2, graph_point_width, graph_point_width);
					
					//Check if knot is leaf
					if(showLeafs == true) {
						String prevLayer = "Layer"+l;
						List<String> leafs = leafs4Layers.get(prevLayer);
						int nLeafs = leafs.size();
						if(nLeafs != 0) {
							int n_prevKnots = tree.get(prevLayer).size();
							for(int i=0; i<n_prevKnots; i++) {
								String prevKnot = "Knot"+(i+1);
								int [] leafIdx = Utilities.get_idx(leafs, prevKnot);
								if(leafIdx[0] != -1) {
									//If knot is leaf then change color to leaf color!
									g2.setColor(leafColor);									
									x1 = points.get((l-1)).get(i).x+graph_point_width/2;
									y1 = points.get((l-1)).get(i).y+graph_point_width/2;
									g2.fillOval(x1-graph_point_width/2, y1-graph_point_width/2, graph_point_width, graph_point_width);
								}
							}
							g2.setColor(pointColor);
						}
					}	
				}	
			}
			
			points.add(layerKnotPoints);
					
		}
		

		if(showLeafs == true) {
			//Leafs in last layer (at the end of tree)
			g2.setColor(leafColor);			
			String layerLabel = "Layer"+treeDepth;
			int nKnots = tree.get(layerLabel).size();
			for(int i=0; i<nKnots; i++) {
				int x1 = points.get((treeDepth-1)).get(i).x+graph_point_width/2;
				int y1 = points.get((treeDepth-1)).get(i).y+graph_point_width/2;
				g2.fillOval(x1-graph_point_width/2, y1-graph_point_width/2, graph_point_width, graph_point_width);	
			}	
			g2.setColor(pointColor);
		}
		
	}
	
	
	@SuppressWarnings("static-access")
	public void create_and_set_info_box(int layerNumber, int knotNumber) {
		
		g2.setStroke(graph_stroke);
		
		int treeDepth = cart_obj.getTreeDepth();
		
		if(layerNumber>treeDepth) {
			throw new RuntimeException("Invalid layer number supplied for info box in tree.");
		}
		
		String layerLabel = "Layer"+layerNumber;
		
		int nKnots = internal_knot_positions.get(layerLabel).size();
		
		if(knotNumber>nKnots) {
			throw new RuntimeException("Invalid knot number supplied for info box in tree.");
		}
			
		String knotLabel  = "Knot"+knotNumber;
		int internalPos = internal_knot_positions.get(layerLabel).get(knotLabel);
		
		int rectangularWidth = (int) ((getWidth()-(numberOfPlotColumns+1)*border_gap)/numberOfPlotColumns);
		int rectangularHeight =(int) ((getHeight()-(numberOfPlotRows+1)*border_gap)/numberOfPlotRows);
		
		//int dummyNumberOfKnots = get_dummy_knots_4_layers(treeDepth)[(layerNumber-1)][0];
		
		int maxNumberOfKnots = get_maxNumberOfKnotsInLayer();
		maxNumberOfKnots = 2*maxNumberOfKnots+1;
		
		int distBetweenKnots = 0;
		int distBetweenLayers = 0;
		int xKnotPos = 0;
		int yKnotPos = 0;
		
		if(rotateTree == false) {
			distBetweenLayers = (int) ((double) rectangularHeight/((double) treeDepth-1));
			distBetweenKnots = (int) rectangularWidth/(maxNumberOfKnots+1);
			xKnotPos = border_gap+distBetweenKnots*internalPos;
			yKnotPos = border_gap+distBetweenLayers*(layerNumber-1);
		}
		if(rotateTree == true) {
			distBetweenLayers = (int) rectangularWidth/(treeDepth-1);
			distBetweenKnots = (int) rectangularHeight/(maxNumberOfKnots+1);
			xKnotPos = border_gap+distBetweenLayers*(layerNumber-1);
			yKnotPos = border_gap+distBetweenKnots*internalPos;
		}
		
		String [] infoStrings = get_info_str_4_box(layerNumber,knotNumber);
		
		HashMap<String, Integer> box_height_and_width = get_box_height_and_width();
		int boxWidth = box_height_and_width.get("Width");
		int boxHeight = box_height_and_width.get("Height");
		
		int vertDistFromKnot = (int) graph_point_width/2;
		int horDistFromKnot  = (int) graph_point_width/2;
		int xArrowWidth = (int) graph_point_width/2;
		xArrowWidth -= 2;
		int yArrowWidth = xArrowWidth; 
		
		//Info Box:		
		int x = xKnotPos+horDistFromKnot;
		int y = yKnotPos+vertDistFromKnot;
		int xPoly = x;
		int yPoly = y;
		
		//Check Distance to plot frame (leftside/ downside):
		int leftSideDist = getWidth()-(x+boxWidth);
		int downSideDist = getHeight()-(y+boxHeight);
			
		if(leftSideDist<0) {
			x -= boxWidth+2*horDistFromKnot;
			xPoly = xKnotPos-horDistFromKnot;
			xArrowWidth *=-1;
		}
		if(downSideDist<0) {
			y -=boxHeight+2*vertDistFromKnot;
			yPoly = yKnotPos-vertDistFromKnot;
			yArrowWidth *= -1;
		}
		
		g2.setColor(Color.WHITE);
		g2.fillRect(x, y, boxWidth,boxHeight);
		g2.setColor(lineColor);
		g2.drawRect(x, y, boxWidth,boxHeight);
		//Arrow to knot:
		Polygon p = new Polygon();
		p.addPoint(xPoly, yPoly);
		p.addPoint(xPoly, yPoly+yArrowWidth);
		p.addPoint(xKnotPos,yKnotPos);
		p.addPoint(xPoly+xArrowWidth, yPoly);
		
		g2.setColor(Color.WHITE);
		g2.fillPolygon(p);
		g2.setColor(lineColor);
		g2.drawPolygon(p);
		g2.setColor(Color.WHITE);
		g2.drawLine(xPoly, yPoly, xPoly+xArrowWidth, yPoly);
		g2.drawLine(xPoly, yPoly, xPoly, yPoly+yArrowWidth);
		g2.setColor(lineColor);
		g2.drawLine(xPoly+xArrowWidth, yPoly, xKnotPos, yKnotPos);
		g2.drawLine(xPoly, yPoly+yArrowWidth, xKnotPos, yKnotPos);
		
		Point boxCoordinates = new Point(x, y);
		set_info_str_2_box(boxCoordinates, infoStrings);
		
		//TODO: Scale distance between Boxes (here 5)
		//Flag for BarPlot!
		//int x4DistributionPlot = x+boxWidth+5;
		//int y4DistributionPlot = y;
		//Point distPlotCoordinates = new Point(x4DistributionPlot, y4DistributionPlot);
		//createAndSetBarPlotOfClassDist(distPlotCoordinates, layerNumber, knotNumber);
		
	}
	
	
	public HashMap<String, Integer> get_box_height_and_width(){
		
		HashMap<String, Integer> box_height_and_width = new HashMap<String, Integer>();
		
		String [] infoStrings = get_info_str_4_box(2,1);
		int nInfoKeys = infoStrings.length;
		
		Font orgFont = g2.getFont();
		g2.setFont(infoFont);		
    	FontMetrics metrics = g2.getFontMetrics();
		
    	int maxStrWidth = 0;
    	for(int s=0; s<nInfoKeys; s++) {
    		int strWidth = metrics.stringWidth(infoStrings[s]);
    		if(strWidth>maxStrWidth) {
    			maxStrWidth = strWidth;
    		}
    	}
    	
    	maxStrWidth += 6;
    	
    	g2.setFont(orgFont);
    	
    	box_height_and_width.put("Width", maxStrWidth);
    	box_height_and_width.put("Height", (nInfoKeys+1)*metrics.getHeight());
		
    	return box_height_and_width;
    	
	}
	
	
	@SuppressWarnings("static-access")
	public void createAndSetBarPlotOfClassDist(Point plotCoordinates, int layerNumber, int knotNumber){
				
		HashMap<String, Integer> box_height_and_width = get_box_height_and_width();
		int boxWidth = box_height_and_width.get("Width");
		int boxHeight = box_height_and_width.get("Height");
	    int frameWidth = 18;
		
		int [][] classDist = cart_obj.get_knot_class_distribution(layerNumber, knotNumber);
	    String [] classLabels = {"A", "B", "C"};
	    //--------------
	    //TODO: Delete!
	    //int n=10;
	    //int [][] classDist = new int [n][1];
	    //String [] classLabels = new String [n];
	    //for(int i=0; i<n; i++) {
	    //	classDist[i][0]=i+11;
	    //	classLabels[i] = "A";
	    //}
	  //--------------
	    
		int nClasses = classDist.length;
		int maxNumber = Utilities.getMax(classDist);
		
		//TODO:maxNumber < 4 (see also yAxis)		
		if(maxNumber % 4 != 0) {
			while(maxNumber % 4 != 0) {
				maxNumber++;
			}
			
		}
		
   	    double barWidth = (double) (boxWidth-frameWidth)/nClasses;	
   	    double scale    = (double) (boxHeight-frameWidth)/(double) maxNumber;
	    
   	    int x = plotCoordinates.x;
   	    int y = plotCoordinates.y;
   	    
   	    //Generate Box for Distribution plot:
   	    g2.setColor(Color.WHITE);
		g2.fillRect(x, y, boxWidth, boxHeight);
		g2.setColor(lineColor);
		g2.drawRect(x, y, boxWidth, boxHeight);
		
		//TODO: Generic procedure for bar colors!
		Color [] colors = new Color [3];
		colors[0] = Color.RED;
		colors[1] = Color.GREEN;
		colors[2] = Color.YELLOW;
		
		Stroke stroke = new BasicStroke(1f);
		
		Font baseFont = g2.getFont();
		Font labelFont = new Font(null, Font.BOLD, 7);	
		g2.setFont(labelFont);
    	FontMetrics metrics = g2.getFontMetrics();
    	int strHeight = metrics.getHeight();
			
		for(int c=0; c<nClasses; c++) {
			//Draw bar for class:
			g2.setColor(colors[c]);
			int barLeftXPos  = (int) ((x+frameWidth*3/4)+c*barWidth);
			int barHeight = (int) scale*classDist[c][0];
			int barUpperYPos = y+boxHeight-frameWidth*3/4-barHeight;			
			g2.fillRect(barLeftXPos, barUpperYPos, (int) barWidth, barHeight);
			g2.setColor(Color.WHITE);
			g2.setStroke(stroke);
			g2.drawRect(barLeftXPos, barUpperYPos, (int) barWidth, barHeight);
			//Set class label:
			int xLabelPos = (int) ((x+frameWidth/2)+c*barWidth+barWidth*1/2);
			int yLabelPos = barUpperYPos+barHeight+strHeight*9/10;
			g2.setFont(labelFont);
			g2.setColor(lineColor);
			g2.drawString(classLabels[c],xLabelPos,yLabelPos); 
		}
		
		g2.setColor(lineColor);
		g2.setFont(baseFont);
		
		//Set x-Axis:
		int lowerYPos = y+boxHeight-frameWidth*7/10;
		int xPos = (int)(x+frameWidth*7/10);
		int x2Pos = (int) ((x+frameWidth*3/4)+(nClasses-1)*barWidth)+(int) barWidth;
		
		g2.drawLine(xPos, lowerYPos, x2Pos, lowerYPos);
		g2.setColor(Color.WHITE);
		g2.drawLine(xPos, lowerYPos-1, x2Pos, lowerYPos-1);
		g2.setStroke(graph_stroke);
		
		//Set y-Axis:	
		int scaledMaxNumber = (int) scale*maxNumber;
		int upperYPos = y+boxHeight-frameWidth*3/4-scaledMaxNumber;
		
		set_y_axis4DistPlot(xPos,lowerYPos,upperYPos,maxNumber);
		
	}
	
	
	public void set_y_axis4DistPlot(int x, int lowerYpos, int upperYpos, int maxNumber) {
		
		//TODO: Cases: maxNumber<4!
		int nSteps  = 4;
		int axisLength = lowerYpos-upperYpos;
		int stepLength = (int) axisLength/nSteps;
		
		//Is max number even?
		if(maxNumber % 4 != 0) {
			throw new RuntimeException("Only max. number that is multiple of 4 allowed for scaling the y axis.");
		}
		
		Stroke stroke = new BasicStroke(1f);
		g2.setStroke(stroke);
		
		Font baseFont = g2.getFont();
		Font labelFont = new Font(null, Font.BOLD, 7);	
		g2.setFont(labelFont);
		
		FontMetrics metrics = g2.getFontMetrics();
    	int strWidth = metrics.stringWidth(Integer.toString(maxNumber));
		
		g2.setColor(lineColor);
		g2.drawLine(x, lowerYpos, x, upperYpos);
		
		int stepLength4Labels = maxNumber/nSteps;
		String yLabel = "";
		int xPos4Labels = x-strWidth-2;
		for(int i=0; i<(nSteps+1); i++) {
			yLabel = Integer.toString(stepLength4Labels*i);
			int y = (int) lowerYpos-stepLength*i;
			g2.drawString(yLabel,xPos4Labels,y); 
		}
		
		g2.setFont(baseFont);
		g2.setStroke(graph_stroke);
		
	}
	
	
	@SuppressWarnings("static-access")
	public String [] get_info_str_4_box(int layerNumber, int knotNumber) {
		
		HashMap<String, List<String>> box_infos = get_knot_infos_4_info_box(layerNumber, knotNumber);
		
		Object [] infoKeys = box_infos.keySet().toArray();
		boolean isCategoricalTree = cart_obj.isCategoricalTree();
		int nKeys = infoKeys.length;
		int nKey4Box = nKeys;
		if(isCategoricalTree == false) {
			nKey4Box --;
		}
		
		String [] keys2Round = {"Cost", "Cost Reduction", "Threshold"};
		
		String [] infoStrings = new String [nKey4Box];
		int idx = 0;
		
		for(int i=0; i<nKeys; i++) {
			
			String key = (String) infoKeys[i];
			if(key.contentEquals("nElementsOfClasses")==true) {				
				if(isCategoricalTree == true) {
					int nClasses = box_infos.get(key).size();
					infoStrings[idx] = "Class Dist.: [";
					for(int c=0; c<(nClasses-1); c++) {
						infoStrings[idx] += box_infos.get(key).get(c) + ", ";
					}
					infoStrings[idx] += box_infos.get(key).get((nClasses-1))+ "]";
					idx++;
				}				
			}else {
				
				int [] roundIdx = Utilities.get_idx(keys2Round, key);
				if(roundIdx[0] != -1) {
					String formatPattern = getDecimalPlacesFormat();
					DecimalFormat df = new DecimalFormat(formatPattern);
					df.setRoundingMode(RoundingMode.CEILING);
					double keyValue = Double.parseDouble(box_infos.get(key).get(0));
					String str = df.format(keyValue);
					if(key.contentEquals("Threshold")==true) {
						if(knotNumber %2 != 0) {
							infoStrings[idx] = key + " <= " + str;
						}else {
							infoStrings[idx] = key + " > " + str;
						}
					}else {
						infoStrings[idx] = key + ": " + str;
					}					
				}else {
					infoStrings[idx] = key + ": " + box_infos.get(key).get(0);
				}	
				idx++;
			}	
			
		}
			
		return infoStrings;
		
	}
	
	
	public void set_info_str_2_box(Point boxCoordinates, String [] infoStrings) {
		
		int x = boxCoordinates.x;
		int y = boxCoordinates.y;
		int distFromBoxBound = 3;
		
		Font orgFont = g2.getFont();
		g2.setFont(infoFont);
		
		String str = "";
    	FontMetrics metrics = g2.getFontMetrics();
    	int strHeight = metrics.getHeight();
		
        int nKeys = infoStrings.length;
        for(int i=0; i<nKeys; i++) {   
        	y += strHeight;
        	str = infoStrings[i];
        	g2.drawString(str,x+distFromBoxBound,y);    	
        }
        
    	g2.setFont(orgFont);
    	
	}
	
	
	@SuppressWarnings("static-access")
	public static HashMap<String, List<String>> get_knot_infos_4_info_box(int layerNumber, int knotNumber){
		
		HashMap<String, List<String>> knot_infos = new HashMap<String, List<String>>();
		
		//Costs
		List<String> info = new ArrayList<String>();
		info.add(Double.toString(cart_obj.get_knot_cost(layerNumber, knotNumber)));
		knot_infos.put("Cost", info);
		
		//Number of elements
		info = new ArrayList<String>();
		info.add(Integer.toString(cart_obj.get_knot_number_of_elements(layerNumber, knotNumber)));
		knot_infos.put("Sample Size", info);	
		
		//Distribution of classes
		info = new ArrayList<String>();
		int [][] classDist = cart_obj.get_knot_class_distribution(layerNumber, knotNumber);
		int nClasses = classDist.length;
		
		for(int c=0; c<nClasses; c++) {
			info.add(Integer.toString(classDist[c][0]));
		}
		
		knot_infos.put("nElementsOfClasses", info);	
		
		//SplittingFeature
		info = new ArrayList<String>();
		info.add(cart_obj.get_knot_splitFeature(layerNumber, knotNumber));
		knot_infos.put("Feature", info);	
		
		//Threshold
		info = new ArrayList<String>();
		info.add(Double.toString(cart_obj.get_knot_threshold(layerNumber, knotNumber)));
		knot_infos.put("Threshold", info);
		
		//CostReduction
		info = new ArrayList<String>();
		info.add(Double.toString(cart_obj.get_knot_costReduction(layerNumber, knotNumber)));
		knot_infos.put("Cost Reduction", info);
		
		return knot_infos;
		
	}
	
	
	@SuppressWarnings("static-access")
	public boolean isTreeSetted() {
		
		boolean treeCheck = true;
		HashMap <String, HashMap<String, HashMap<String,List<String>>>> tree4Check = cart_obj.get_tree();
		
		if(tree4Check==null) {
			treeCheck = false;
		}
		
		return treeCheck;
		
	}
	
	
	@SuppressWarnings("static-access")
	public int get_maxNumberOfKnotsInLayer() {
		
		HashMap <String, HashMap<String, HashMap<String,List<String>>>> tree = cart_obj.get_tree();
		int nLayers = tree.size();
		int maxKnots = 0;
		for(int l=0; l<nLayers; l++) {
			String layerLabel = "Layer"+(l+1);
			int nKnots = tree.get(layerLabel).size();
			if(nKnots>maxKnots) {
				maxKnots = nKnots;
			}
		}
				
		return maxKnots;	
	}
	
	
	public void get_internal_knot_positions(){
		
		int maxNumberOfKnotsInLayer = get_maxNumberOfKnotsInLayer();
			
		HashMap<String, HashMap<String, Integer>> int_knot_positions = new HashMap<String, HashMap<String, Integer>>();
		
		int rootPos = maxNumberOfKnotsInLayer+1;
		
		HashMap<String, Integer> int_pos = new HashMap<String, Integer>();
		int_pos.put("Knot1",rootPos);
		int_knot_positions.put("Layer1",int_pos);
		
		@SuppressWarnings("static-access")
		HashMap <String, HashMap<String, HashMap<String,List<String>>>> tree = cart_obj.get_tree();
		int treeDepth = tree.size();
		
		for(int l=1; l<treeDepth; l++) {
			int_pos = new HashMap<String, Integer>();
			String layerLabel = "Layer" + (l+1);
			HashMap<String, HashMap<String,List<String>>> layer = tree.get(layerLabel);
			int nKnots = layer.size();
			String prevLayerLabel = "Layer"+l;
			
			HashMap<String, Integer> prev_int_pos = int_knot_positions.get(prevLayerLabel);
			List<Integer> leftTreeIdxs = new ArrayList<Integer>();
			List<Integer> rightTreeIdxs = new ArrayList<Integer>();
			
			if(l!=1) {
				for(int k=0; k<nKnots; k++) {
					int knotNumber = k+1;
					String knotLabel = "Knot"+knotNumber;
					String parentKnot = tree.get(layerLabel).get(knotLabel).get("ParentKnot").get(0);
					int prev_pos = prev_int_pos.get(parentKnot);
					if(prev_pos<rootPos) {
						leftTreeIdxs.add(k);
					}else {
						rightTreeIdxs.add(k);
					}
				}
			}else {
				leftTreeIdxs.add(0);
				rightTreeIdxs.add(1);
			}
			
			int nLeftTreeKnots = leftTreeIdxs.size();
			int nRightTreeKnots = rightTreeIdxs.size();
			int knot_int_pos = 0;
			
			for(int k=0; k<nLeftTreeKnots; k++) {
				int idx = nLeftTreeKnots-k-1;
				int knotNumber = leftTreeIdxs.get(idx)+1;
				String knotLabel = "Knot"+knotNumber;
				knot_int_pos = rootPos-k-1;			
				int_pos.put(knotLabel, knot_int_pos);
			}
			
			for(int k=0; k<nRightTreeKnots; k++) {
				int knotNumber = rightTreeIdxs.get(k)+1;
				String knotLabel = "Knot"+knotNumber;
				knot_int_pos = rootPos+k+1;
				int_pos.put(knotLabel, knot_int_pos);
			}
			
			int_knot_positions.put(layerLabel, int_pos);
		}
		
		internal_knot_positions = int_knot_positions;
		
	}
	
	
	public void setCARTObject(CART cart_object) {
		cart_obj = cart_object;
	}
	
	
	public void plotDecisionTree(){
		  		
		plotInfo.add(0);		
		plot();
		
	}
	
	
	public void rotateTree() {
		rotateTree = true;
	}
	
	
	public void setDecimalPlaces4InfoBox(int decPlaces) {
		decimalPlaces4InfoBox = decPlaces;
	}
	
	
	public String getDecimalPlacesFormat() {
		
		String pattern = "#.";
		for(int i=0; i<decimalPlaces4InfoBox; i++) {
			pattern += "#";
		}
		
		return pattern;
		
	}
	
	
	public void showLeafs(boolean show) {
		showLeafs = show;
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotFittedExplainedVariable() {
		
		LinearRegression obj_linearReg = cart_obj.get_linearRegObject();
		
		if(obj_linearReg == null) {
			System.out.println("No results of linear regression supplied for plot.");
			return;
		}
		
		double [][] explained_var        = cart_obj.get_explained_variables();
		double [][] fitted_explained_var = obj_linearReg.get_fitted_values();
		
		int n=explained_var.length;
		
		double [][] x = new double [n][1];
		for(int i=0; i<n; i++) {
			x[i][0]=i+1;
		}
		
		GenGraphics regPlot = new GenGraphics();
		
		regPlot.setGraphWidth(pref_w);
		regPlot.setGraphHeight(pref_h);
		
	    regPlot.plotLines(x,fitted_explained_var,true);
	    regPlot.plotPoints(x, explained_var,false);
	    
	    String [] title = {"Fitted vs. Observed"};
	    String var_name = cart_obj.get_name_of_explained_variable();
	    String [] yLabel = {var_name};
	    String [] subTitle = {var_name};
	    String [] xLabel = {"No. of Observation"};
	    regPlot.setTitle(title, null, "12");
	    regPlot.setSubTitle1(subTitle, null, "11");
	    regPlot.setYLabel(yLabel, null, "11");
	    regPlot.setXLabel(xLabel, null, "10");
	    regPlot.setNumberOfDigits4YAxis(2);
	    regPlot.setNumberOfDigits4XAxis(0);   
	    regPlot.setFontOfXAxisUnits("plain", 11);
	    regPlot.setFontOfYAxisUnits("plain", 11);
 	
	    regPlot.plot();
	    
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotHistOfResiduals() {
		
		LinearRegression obj_linearReg = cart_obj.get_linearRegObject();
		
		if(obj_linearReg == null) {
			System.out.println("No results of linear regression supplied for plot.");
			return;
		}
		
		double [][] residuals = obj_linearReg.get_residuals();
					 	
	 	String [] titles = {"Error Distribution"};
	 	String var_name = cart_obj.get_name_of_explained_variable();	 	
	 	String [] subTitles = {var_name};	 	
	 	String [] yLabels = {"Empirical Distribution"};
	 	String [] xLabels = {"Error"};
	 	
	 	HistGraphics histPlot = new HistGraphics();
	 	
	 	histPlot.plotHistogram(residuals,true,true);
	 	histPlot.useDifferentColors = false;
	 	histPlot.barColor = new Color(44, 102, 230, 180);
	 	histPlot.setNumberOfPlotColums(1);
	 	histPlot.setNumberOfPlotRows(1);
	 	histPlot.set_numberOfBins(50);
	 	histPlot.setTitle(titles, null, null);
	 	histPlot.setSubTitle1(subTitles, null, null);
	 	histPlot.setYLabel(yLabels, null, null);
	 	histPlot.setXLabel(xLabels, null, null);
	 	histPlot.setGraphWidth(pref_w);
	 	histPlot.setGraphHeight(pref_h);
	 	histPlot.setNumberOfDigits4XAxis(0);
	 	histPlot.setNumberOfDigits4YAxis(2);
	 	histPlot.setFontOfXAxisUnits("bold",11);
	 	histPlot.setFontOfYAxisUnits("bold",11);
	 	//set_max_x_value(max);
	 	//set_min_x_value(min);
	 	histPlot.freqHist(false);	
	 	histPlot.setPDFLineWidth(2);
	 	//noLinesAroundBars();
	 	histPlot.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotLeafWeights() {
		
		LinearRegression obj_linearReg = cart_obj.get_linearRegObject();
		
		if(obj_linearReg == null) {
			System.out.println("No results of linear regression supplied for plot.");
			return;
		}
		
		double [][] weights = obj_linearReg.get_est_parameters();
		
		int n=weights.length;
		
		double [][] x = new double [n][1];
		for(int i=0; i<n; i++) {
			x[i][0]=i+1;
		}
		
		GenGraphics weightPlot = new GenGraphics();
		
		weightPlot.setGraphWidth(pref_w);
		weightPlot.setGraphHeight(pref_h);
		
		weightPlot.plotLines(x,weights,true);
	    
	    String [] title = {"Estimated Leaf Weights"};
	    String [] yLabel = {"Leaf Weights"};
	    String [] xLabel = {"No. of Leaf"};
	    weightPlot.setTitle(title, null, "12");
	    weightPlot.setYLabel(yLabel, null, "11");
	    weightPlot.setXLabel(xLabel, null, "10");
	    weightPlot.setNumberOfDigits4YAxis(2);
	    weightPlot.setNumberOfDigits4XAxis(0);   
	    weightPlot.setFontOfXAxisUnits("plain", 11);
	    weightPlot.setFontOfYAxisUnits("plain", 11);
 	
	    weightPlot.plot();
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotDistOfSplittingFeature(){
		
		BarGraphics barGraph = new BarGraphics();		
		
		String [] explaining_vars = cart_obj.get_names_of_explaining_variables();
		int n_explaining_vars = explaining_vars.length;
		
		HashMap<String,Integer> dist = cart_obj.get_distOfSplittingFeatures();
		
		if(dist.size() != n_explaining_vars) {
			throw new RuntimeException("Mismatch between number of explaining variables and distribution of splitting features.");
		}
		
		String [][] x_values = new String [n_explaining_vars][1];
		double [][] y_values = new double [n_explaining_vars][1];
		
		for(int i=0; i<n_explaining_vars; i++) {
			x_values[i][0] = explaining_vars[i];
			y_values[i][0] = dist.get(explaining_vars[i]);
		}
				
 		String [] title    = {"Splitting Features"};
 		//String [] subTitle = {"SubTitle1"};
 		String [] yLabel   = {"Number of Splits"};
 		String [] xLabel   = {"Feature"};
 	 	 		
 		barGraph.setNumberOfPlotColums(1);
 		barGraph.setNumberOfPlotRows(1);
 	
 		barGraph.setGraphWidth(pref_w);
 		barGraph.setGraphHeight(pref_h);
 	 		
 		barGraph.plotGroupedBars(x_values,y_values, true);
 		
 		barGraph.setTitle(title, null, "12");
 		if(subTitle1.size() != 0 ) {
 			String [] subTitle = {subTitle1.get(0)};
 			barGraph.setSubTitle1(subTitle, null, "10");
 		}
 			
 		barGraph.setYLabel(yLabel, null, "12");
 		barGraph.setXLabel(xLabel, null, "10");
 		barGraph.setNumberOfXDivisions(1);
 		barGraph.setNumberOfDigits4YAxis(0);		   
 		barGraph.setFontOfXAxisUnits("plain", 11);
 		barGraph.setFontOfYAxisUnits("plain", 11);
 	
 		barGraph.plot();
 	
	}
	
	
    private static void createAndShowGui() {

    	DecisionTreeGraphics mainPanel = new DecisionTreeGraphics();

    	String frameLabel;
    	
    	frameLabel = "Decision Tree Plot";    		
    	
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
	
}
