package Kernels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;
import Utilities.Utilities;

public class TestKernelDensities {

	//Test 1: Usage of different kernels and bandwidths
	public static void test1KernelDensities() {
		
		double [] bandwidth = {0.25, 0.50, 0.75};
		
		long startTime = System.currentTimeMillis();
		
		HashMap<String, double [][]> testDist = getMultimodalTestDist();
		double [][] sample = testDist.get("sample");
		double [][] truePDF = testDist.get("probs");
		
		HashMap<Double, ArrayList<double [][]>> estRes = new HashMap<Double, ArrayList<double [][]>>();
		int nBws = bandwidth.length;
		
		KernelDensities kd = new KernelDensities(sample, bandwidth[0]);
		
		String [] kernels = kd.getValidKernels();
		int nKernels = kernels.length;
		
		for(int b=0; b<nBws; b++) {
			ArrayList<double [][]> fittedPDFs = new ArrayList<double [][]>(nKernels);			
			for(int i=0; i<nKernels; i++) {
				kd = new KernelDensities(sample, bandwidth[b], kernels[i], "euclidian");
				kd.calcKernelDensity();
				fittedPDFs.add(kd.getFittedKernelDensities());
			}
			estRes.put(bandwidth[b], fittedPDFs);
		}
		
		
		plotDensities4Test(sample, truePDF, estRes);
			
		long endTime = System.currentTimeMillis();
		System.out.println("Test 1 for kernel density estimation (KDE) finished after: " + ((endTime-startTime)/1000.0) + " secs.");			
	}
	
	
	//Test 2: Usage of different distance metrics (Euclidian, squared Euclidian, weighted Minkowski)
	public static void test2KernelDensities() {
		
		double bw = 0.5;
		
		HashMap<String, double [][]> testDist = getMultimodalTestDist();
		double [][] sample = testDist.get("sample");
		double [][] truePDF = testDist.get("probs");
		
		HashMap<Double, ArrayList<double [][]>> estRes = new HashMap<Double, ArrayList<double [][]>>();
		ArrayList<double [][]> fittedPDFs = new ArrayList<double [][]>(2);
		
		double [][] w = new double [1][1];
		w[0][0] = 1.1;
		
		KernelDensities kd = new KernelDensities(sample, bw, "gaussian", "euclidian");
		kd.calcKernelDensity();
		fittedPDFs.add(kd.getFittedKernelDensities());
		
		kd = new KernelDensities(sample, bw, "gaussian", "sqeuclidian");
		kd.calcKernelDensity();
		fittedPDFs.add(kd.getFittedKernelDensities());
		
		kd = new KernelDensities(sample, bw, "gaussian", "wminkowski");
		kd.setAdditionalDistancePars(2.0, w);
		kd.calcKernelDensity();
		fittedPDFs.add(kd.getFittedKernelDensities());
		
		estRes.put(bw, fittedPDFs);
		
		plotDensities4Test(sample, truePDF, estRes);
		
	}
	
	
	//Test 3: Usage of multivariate input data for kernel density estimation
	public static void test3KernelDensities() {
		
		double bw = 0.25;
		int dim = 5;
        int nRand = 200;
		
		double [][] xLabels = new double [nRand][1];
				
		for(int i=0; i<nRand; i++) {
			
			xLabels[i][0] = (i+1);				
		}
		
		double [][] my = new double [dim][1];
		
		for(int i=0; i<dim; i++) {
			my[i][0] = i;
		}
		
		double [][] Sigma = MatrixOperations.identity(dim);
		
		NormalDistribution mnd = new NormalDistribution(my,Sigma);
		double [][] x = mnd.sample(nRand);
		
		double [][] mvPDF = new double [nRand][1];
		for(int i=0; i<nRand; i++) {
			double [][] x_sel = MatrixOperations.get_column_vec_from_matrix(x, i);
			mvPDF[i][0] = mnd.get_multivariateNormalPDF(x_sel);
		}
		
		HashMap<Double, ArrayList<double [][]>> estRes = new HashMap<Double, ArrayList<double [][]>>();
		ArrayList<double [][]> fittedPDFs = new ArrayList<double [][]>(2);

		x = MatrixOperations.transpose(x);
				
		KernelDensities kd = new KernelDensities(x, bw, "gaussian");
		kd.calcKernelDensity();
		fittedPDFs.add(kd.getFittedKernelDensities());
		
		estRes.put(bw, fittedPDFs);
		
		plotDensities4Test(xLabels, mvPDF, estRes);   	
		
	}
	
	
	public static HashMap<String, double [][]> getMultimodalTestDist() {
		
		int n = 100;
		double w1 = 0.3;
		double w2 = (1.0-w1);
		int n1 = (int) (w1*n);
		int n2 = (int) (w2*n);
		
		NormalDistribution objND1 = new NormalDistribution(0.0, 1.0);
		double [][] sample1 = objND1.sample(n1);
	
		NormalDistribution  objND2 = new NormalDistribution(5.0, 1.0);
		double [][] sample2 = objND2.sample(n2);

		//Bimodal distribution
		double [][] sample = MatrixOperations.rbind(sample1, sample2);
		
		HashMap<String, List<Double>> sortRes = Utilities.get_sorted_elements_and_idxs_of_double_vector(sample);
		
		double [][] sortedSample = new double [n][1];
				
 		for(int i=0; i<n; i++) {
			sortedSample[i][0] = sortRes.get("SortedValues").get(i);
		}
		
		sample = sortedSample;
		
		double [][] truePDF1 = objND1.get_univariateNormalPDF(sample);
		double [][] truePDF2 = objND2.get_univariateNormalPDF(sample);
		double [][] truePDF  = MatrixOperations.add(MatrixOperations.scalar_multiplication(w1,truePDF1), MatrixOperations.scalar_multiplication(w2,truePDF2));
		
		HashMap<String, double [][]> testDist = new HashMap<String, double [][]>();
		testDist.put("sample", sample);
		testDist.put("probs", truePDF);
		
		return testDist;
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotDensities4Test(double [][] sample, double [][] truePDF, HashMap<Double, ArrayList<double [][]>> fittedPDFs){
   			 	
	 	int nBws = fittedPDFs.size();
	 	Object [] keys = fittedPDFs.keySet().toArray();
	 	
	 	int nFittedPDFs = fittedPDFs.get(keys[0]).size();
	 		
	 	GenGraphics graph = new GenGraphics();
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setNumberOfPlotColums(nBws);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(650);
	 	graph.setGraphHeight(400);
	 	
	 	ArrayList<Color> defColors = graph.getDefaultColorsAsList();
	 	
	 	String [] title    = new String [nBws];
	 	String [] subTitle = new String [nBws];
	 	String [] xLabel   = new String [nBws];
	 	String [] yLabel   = new String [nBws];
	 	
	 	for(int b=0; b<nBws; b++) {
	 		boolean newPlot = true;
	 		for(int i=0; i<nFittedPDFs; i++) {
	 			if(i!=0) {
	 				newPlot = false;
	 			}
		 		graph.plotLines(sample,fittedPDFs.get(keys[b]).get(i),newPlot, defColors.get(i));
		 	}
	 		graph.plotLines(sample,truePDF,false, Color.BLACK);
		 	graph.plotPoints(sample,truePDF,false, Color.RED);
		 	
		 	title[b] = "Kernel Density Estimation";
		 	subTitle[b] = "Bandwidth: " + keys[b];
		 	yLabel[b]   = "Observed & estimated PDF";
		 	xLabel[b]   = "Random Numbers";
	 	}
	 	
	 	graph.setTitle(title, null, "12");
	 	graph.setSubTitle1(subTitle, null, "10");
	 	graph.setYLabel(yLabel, null, "10");
	 	graph.setXLabel(xLabel, null, "10");
	 	graph.setNumberOfDigits4XAxis(2);   
	 	graph.setNumberOfDigits4YAxis(2);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	
	 	graph.plot();
   		
   	}
	
	
	public static void main(String[] args) {
		test1KernelDensities();
		//test2KernelDensities();
		//test3KernelDensities();
	}
	
}
