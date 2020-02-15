package Kernels;

import java.awt.Color;
import java.util.HashMap;
import java.util.List;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class TestKernelDensities {

	public static void test1KernelDensities() {
		
		int n = 100;
		
		NormalDistribution objND1 = new NormalDistribution(0.0, 1.0);
		double [][] sample1 = objND1.sample((int)0.5*n);
	
		NormalDistribution  objND2 = new NormalDistribution(5.0, 1.0);
		double [][] sample2 = objND2.sample((int)0.5*n);

		//Bimodal distribution
		double [][] sample = MatrixOperations.rbind(sample1, sample2);
		
		HashMap<String, List<Double>> sortRes = Utilities.get_sorted_elements_and_idxs_of_double_vector(sample);
		
		double [][] sortedSample = new double [n][0];
		for(int i=0; i<n; i++) {
			sortedSample[i][0] = sortRes.get("SortedValues").get(i);
		}
		
		sample = sortedSample;
		
		double [][] truePDF1 = objND1.get_univariateNormalPDF(sample);
		double [][] truePDF2 = objND2.get_univariateNormalPDF(sample);
		double [][] truePDF  = MatrixOperations.add(truePDF1, truePDF2);
		
		KernelDensities kd = new KernelDensities(sample, 0.5, "gaussian", "euclidian");
		kd.calcKernelDensity();
		double [][] fittedPDF = kd.getFittedKernelDensities();
				
		plotDensities4Test1(sample, truePDF, fittedPDF);
			
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotDensities4Test1(double [][] sample, double [][] truePDF, double [][] fittedPDF){
   			 	
	 	GenGraphics graph = new GenGraphics();
	 	graph.setNumberOfPlotColums(2);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(1000);
	 	graph.setGraphHeight(600);

	 	graph.plotLines(sample,truePDF,true);
	 	graph.plotShadedArea(sample, truePDF, Color.GRAY);
	 	
	 	graph.plotLines(sample,fittedPDF,true);
	 	graph.plotPoints(sample,fittedPDF,false, Color.RED);
	 	graph.plotShadedArea(sample, fittedPDF);
	 	
	 	String [] title    = {"True PDF", "Fitted PDF"};
	 	String [] yLabel   = {"True PDF", "Fitted PDF"};
	 	String [] xLabel   = {"Value", "Value"};
	 	
	 	graph.setTitle(title, null, "12");
	 	graph.setYLabel(yLabel, null, "10");
	 	graph.setXLabel(xLabel, null, "10");
	 	graph.setNumberOfDigits4XAxis(2);   
	 	graph.setNumberOfDigits4YAxis(2);
	 	graph.setFontOfXAxisUnits("plain", 10);
	 	graph.setFontOfYAxisUnits("plain", 10);
	 	
	 	graph.plot();
   		
   	}
	
}
