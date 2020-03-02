package Kernels;

import java.awt.Color;
import java.util.ArrayList;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.GeneralMath;

public class TestKernelRegression {

	//Test 1: 
	public static void test1KernelRegression() {
		
		int n = 300;
		ArrayList<double [][]> sincf = get_sinc_function_4_test(n);
		
		double [][] x = sincf.get(0);
		double [][] y = sincf.get(1);
		
		double [][] y_noisy = new double [n][1];
		
		for(int i=0; i<n; i++) {
			NormalDistribution nd = new NormalDistribution(y[i][0],0.01);
			y_noisy[i][0] = nd.sample(1)[0][0]; 
		}
		
		KernelRegression kr = new KernelRegression(y, x);
		String [] kernels = kr.getValidKernels();
		
		int nKernels = kernels.length;
		
		ArrayList<double [][]> y_fitted = new ArrayList<double [][]>();
		
		for(int i=0; i<nKernels; i++) {
			kr = new KernelRegression(y, x, 0.5, kernels[i]);
			kr.inSamplePredictionKernelRegression();
			y_fitted.add(kr.get_fitted_values());
		}
			
		plot4Test(x, y, y_noisy, y_fitted);
		
	}
	
	
	public static ArrayList<double [][]> get_sinc_function_4_test(int nSteps) {
		
		double l = -4.0;
		double u = 4.0;
		
		double step = (u-l)/(nSteps-1);
		
		double [][] steps = new double [nSteps][1];
		double [][] sincf = new double [nSteps][1];
		
		for(int i=0; i<nSteps; i++) {
			steps[i][0] = l+step*i;
			sincf[i][0] = GeneralMath.sincf(steps[i][0]);
		}
		
		ArrayList<double [][]> sinf4Test = new ArrayList<double [][]>();
		sinf4Test.add(steps);
		sinf4Test.add(sincf);
				
		return sinf4Test;
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot4Test(double [][] x, double [][] y, double [][] y_noisy, ArrayList<double [][]> y_fitted){
		 	
	 	GenGraphics graph = new GenGraphics();
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(500);
	 	graph.setGraphHeight(400);
	 		 	
	 	String [] title    = new String [1];
	 	String [] xLabel   = new String [1];
	 	String [] yLabel   = new String [1];
	 	
	 	graph.plotPoints(x,y_noisy, true, Color.RED);
	 	
	 	ArrayList<Color> defColors = graph.getDefaultColorsAsList();
	 	
	 	int nKernels = y_fitted.size();
	 	
	 	for(int i=0; i<nKernels; i++) {
	 		graph.plotLines(x,y_fitted.get(i), false, defColors.get(i));
	 	}
	 			
	 	graph.plotLines(x,y,false, Color.BLACK);
		
		title[0] = "Kernel Regression";
		yLabel[0]   = "Observed & fitted values";
		xLabel[0]   = "x values";
	 	
	 	graph.setTitle(title, null, "12");
	 	graph.setYLabel(yLabel, null, "10");
	 	graph.setXLabel(xLabel, null, "10");
	 	graph.setNumberOfDigits4XAxis(1);   
	 	graph.setNumberOfDigits4YAxis(2);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	
	 	graph.plot();
   		
   	}
	
	
	public static void main(String[] args) {
		test1KernelRegression();
	}
	
}
