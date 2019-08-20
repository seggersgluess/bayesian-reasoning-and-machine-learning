package Graphics;

import java.awt.Color;

import Distributions.NormalDistribution;

public class QQ_Plot {

	public static int width = 1000;
	public static int heigth = 600;
	
	@SuppressWarnings("static-access")
	public static void create_qq_plot(double [][] series_4_qq_plot){
		
		int n = series_4_qq_plot.length;
		double [] data = new double [n];
		
		for(int i=0; i<n; i++){
			data[i] = series_4_qq_plot[i][0];
		}
		
		double [] sortedData = Utilities.Utilities.get_sorted_vec(data);
		
    	double [][] mean = new double [1][1];
    	double [][] sd = new double [1][1];   	
    	
    	mean[0][0] = 0.0;
    	sd[0][0] = 1.0;
    	
    	NormalDistribution normalDist = new NormalDistribution(mean, sd);
		
    	double quantile = 0.0;
		double [][] normalQuantiles = new double [n][1];
		
		for(int i=0; i<n; i++){
			series_4_qq_plot[i][0] = sortedData[i];
			quantile = (i+0.5)/n;
			normalQuantiles[i][0] = normalDist.get_univariateNormalQuantile(quantile);
		}
		
        GenGraphics obj_graph = new GenGraphics();
	       
        obj_graph.setNumberOfPlotColums(1);
        obj_graph.setNumberOfPlotRows(1);
	 	
        obj_graph.setGraphWidth(width);
        obj_graph.setGraphHeight(heigth);

	 	obj_graph.plotPoints(normalQuantiles, series_4_qq_plot, true, Color.BLUE);	
	 	obj_graph.plotLines(normalQuantiles, normalQuantiles, false, Color.RED);
	 	 	
	 	String [] title  = {"QQ-Plot"};
	 	String [] subTitle = {"Standard Normal"};
	 	String [] yLabel = {"Sample Quantiles"};
	 	String [] xLabel = {"Standard Normal Quantiles"};
	 		
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setSubTitle1(subTitle, null, "11");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setXLabel(xLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(2);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("bold", 10);
	 	obj_graph.setFontOfYAxisUnits("bold", 10);
	 	
	 	obj_graph.plot();
		
	}
	
	
	public static void set_width(int used_width){
		width = used_width;
	}
	
	
	public static void set_heigth(int used_heigth){
		heigth = used_heigth;
	}
	
}
