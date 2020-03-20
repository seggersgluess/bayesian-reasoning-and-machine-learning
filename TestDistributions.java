package Distributions;

import java.awt.Color;
import java.util.ArrayList;

import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class TestDistributions {

	public static void test1BinomialDist() {
		
		double [] n = new double [2];
		n[0] = 30.0;
		n[1] = 60.0;
		
		double [] p = new double [3];
		p[0] = 0.2;
		p[1] = 0.4;
		p[2] = 0.6;
		
		double [][] x = new double [61][1];
		ArrayList<double [][]> probs = new ArrayList<double [][]>();
		
		for(int i=0; i<n.length; i++) {
			int l = (int)n[i]+1;
			double [][] x_case = new double [l][1];
			for(int j=0; j<l; j++) {
				x_case[j][0] = j;
			}
			for(int j=0; j<p.length; j++) {
				x = new double [61][1];
				double [][] probs_long = new double [61][1];
				Binomial bn = new Binomial(x_case,n[i],p[j]);
				double [][] bn_probs = bn.pmf();
				for(int k=0; k<l; k++) {
					x[k][0] = x_case[k][0];
					probs_long[k][0] = bn_probs[k][0];
				}
				probs.add(probs_long);
			}
		}
		
		plotBinomial(x ,probs);		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotBinomial(double [][] x, ArrayList<double [][]> probs){
   	 
		int n_dist = probs.size();
		
	 	GenGraphics obj_graph = new GenGraphics();
	 	
	 	obj_graph.setNumberOfPlotColums(1);
	 	obj_graph.setNumberOfPlotRows(1);
	 	obj_graph.setGraphWidth(500);
	 	obj_graph.setGraphHeight(300);
	 	
	 	ArrayList<Color> defColors = obj_graph.getDefaultColorsAsList();
	 	
	 	String [] title = {"Binomial Distribution"};
	 	
	 	for(int i=0; i<n_dist; i++) {
	 		if(i==0) {
	 			obj_graph.plotLines(x, probs.get(i), true, defColors.get(i));
	 		}else {
	 			obj_graph.plotLines(x, probs.get(i), false, defColors.get(i));
	 		}
	 		//obj_graph.plotPoints(x, probs.get(i), false, defColors.get(i));
	 	}
	 	
	 	String [] yLabel   = {"PMF"};
	 	String [] xLabel   = {"Number of Success"};
	 	
	 	obj_graph.setBackroundColor(Color.WHITE);
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setXLabel(xLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("plain", 10);
	 	obj_graph.setFontOfYAxisUnits("plain", 10);
	 	
	 	obj_graph.plot();  		
   	}
	
	
	public static void test2BinomialDist() {
		int n_x = 6;
		double [][] x = new double [n_x][1];
		for(int i=0; i<n_x; i++) {
			x[i][0] = i;
		}
		double n = 5.0;
		double p = 0.2;
		
		Binomial bin = new Binomial(x,n,p);
		
		MatrixOperations.print_matrix(bin.pmf());		
	}
	
	
	public static void test1PoissonDist() {
		
		int nGoals = 16;
		
		double [][] x = new double [nGoals][1];

		for(int i=0; i<nGoals; i++) {
			x[i][0] = i;
		}
		
		ArrayList<double [][]> probs = new ArrayList<double [][]>();
		
		double lambda1954 = 5.38;
		double lambda2018 = 2.64;
		
		Poisson poisson = new Poisson(x,lambda1954);		
		probs.add(poisson.pmf());
	
		poisson = new Poisson(x,lambda2018);		
		probs.add(poisson.pmf());
	
		plotPoisson(x ,probs);		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotPoisson(double [][] x, ArrayList<double [][]> probs){
   	 
		int n_dist = probs.size();
		
	 	GenGraphics obj_graph = new GenGraphics();
	 	
	 	obj_graph.setNumberOfPlotColums(1);
	 	obj_graph.setNumberOfPlotRows(1);
	 	obj_graph.setGraphWidth(500);
	 	obj_graph.setGraphHeight(300);
	 	
	 	ArrayList<Color> defColors = obj_graph.getDefaultColorsAsList();
	 	
	 	String [] title = {"Poisson Distribution"};
	 	String [] subTitle = {"       Soccer World Cups 1954 vs. 2008"};
	 	
	 	for(int i=0; i<n_dist; i++) {
	 		if(i==0) {
	 			obj_graph.plotLines(x, probs.get(i), true, defColors.get(i));
	 		}else {
	 			obj_graph.plotLines(x, probs.get(i), false, defColors.get(i));
	 		}
	 		//obj_graph.plotPoints(x, probs.get(i), false, defColors.get(i));
	 	}
	 	
	 	String [] yLabel   = {"PMF"};
	 	String [] xLabel   = {"Number of Goals in Game"};
	 	
	 	obj_graph.setBackroundColor(Color.WHITE);
	 	obj_graph.setTitle(title, null, "12");
	 	obj_graph.setSubTitle1(subTitle, null, "10");
	 	obj_graph.setYLabel(yLabel, null, "10");
	 	obj_graph.setXLabel(xLabel, null, "10");
	 	obj_graph.setNumberOfDigits4XAxis(0);   
	 	obj_graph.setNumberOfDigits4YAxis(2);
	 	obj_graph.setFontOfXAxisUnits("plain", 10);
	 	obj_graph.setFontOfYAxisUnits("plain", 10);
	 	
	 	obj_graph.plot();
   		
   	}
	
	
	public static void test1MVGaussian() {
		
		int nRandNumbers = 1000;
		
		ArrayList<double [][]> randNumbers = new ArrayList<double [][]>();
		ArrayList<double [][]> ellipses = new ArrayList<double [][]>();
		
		double [][] mean = new double [2][1];
		mean [0][0] = 12.0;
		mean [1][0] = 8.0;
		
		double [][] Sigma = new double [2][2];
		Sigma[0][0] = 2.0;
		Sigma[1][0] = -1.2;
		Sigma[0][1] = Sigma[1][0];
		Sigma[1][1] = 6.0;
		
		NormalDistribution nDist = new NormalDistribution(mean, Sigma);
		randNumbers.add(MatrixOperations.transpose(nDist.sample(nRandNumbers)));		
		ellipses.add(nDist.get_2d_confidence_ellipse(nRandNumbers, mean, Sigma, 0.9));
		
		mean [0][0] = 3.0;
		mean [1][0] = 10.0;
		
		Sigma[0][0] = 6.0;
		Sigma[1][0] = 1.2;
		Sigma[0][1] = Sigma[1][0];
		Sigma[1][1] = 6.0;
		
		nDist = new NormalDistribution(mean, Sigma);
		randNumbers.add(MatrixOperations.transpose(nDist.sample(nRandNumbers)));	
		ellipses.add(nDist.get_2d_confidence_ellipse(nRandNumbers, mean, Sigma, 0.9));
		
		mean [0][0] = 0.0;
		mean [1][0] = 0.0;
		Sigma = MatrixOperations.identity(2);
		nDist = new NormalDistribution(mean, Sigma);
		randNumbers.add(MatrixOperations.transpose(nDist.sample(nRandNumbers)));		
		ellipses.add(nDist.get_2d_confidence_ellipse(nRandNumbers, mean, Sigma, 0.9));
		
		plotMVGaussian(randNumbers, ellipses);
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotMVGaussian(ArrayList<double [][]> randNumbers, ArrayList<double [][]> ellipses) {
		
		GenGraphics graph = new GenGraphics();
		
		int nSamples = randNumbers.size();
      		
	 	String [] title    = {"Multivariate Normal Distribution"};
	 	String [] yLabel   = {"x2"};
	 	String [] xLabel   = {"x1"};
	 	
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(500);
	 	graph.setGraphHeight(400);
	 	
	 	ArrayList<Color> defCols = graph.getDefaultColorsAsList();
	 	
	 	boolean newPlot = true;
	 	for(int i=0; i<nSamples; i++) {
			double [][] x1 = MatrixOperations.get_column_vec_from_matrix(randNumbers.get(i), 0);
			double [][] x2 = MatrixOperations.get_column_vec_from_matrix(randNumbers.get(i), 1);
	 		graph.plotPoints(x1,x2, newPlot, defCols.get(i));
	 		newPlot = false;	 		
	 	}
	 	
	 	for(int i=0; i<nSamples; i++) {
			double [][] x1 = MatrixOperations.get_column_vec_from_matrix(ellipses.get(i), 0);
			double [][] x2 = MatrixOperations.get_column_vec_from_matrix(ellipses.get(i), 1);
	 		graph.plotLines(x1,x2, newPlot, Color.RED);
	 	}
	 	
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(title, null, "12");
	 	graph.setYLabel(yLabel, null, "10");
	 	graph.setXLabel(xLabel, null, "10");
	 	graph.setNumberOfDigits4XAxis(1);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 10);
	 	graph.setFontOfYAxisUnits("plain", 10);
	 	graph.set_point_width(4);
	 	
	 	graph.plot();
	}
	

	public static void main(String[] args) {
		//test1BinomialDist();
		//test2BinomialDist();
		//test1PoissonDist();
		test1MVGaussian();
	}
	
}
