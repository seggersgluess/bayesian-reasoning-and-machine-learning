package ComponentModels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;

import DataManagement.InputDataManager;
import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.GeneralMath;
import Mathematics.MatrixOperations;

public class TestFactorAnalysis {

	
	public static void test1() {
		
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/DecisionTrees/Classification/IrisData.txt";
		InputDataManager inputData = new InputDataManager();	
		try {
			inputData.fileReader(dirName, false, true, true);
		} catch (Exception e) {
			System.out.println("Can´t find data set for upload.");
			e.printStackTrace();
		}
		
		String [][] input = inputData.strFileData;
		
		int n_rows = input.length;
		int n_cols = input[0].length-1;
		double [][] data = new double [n_rows][n_cols];
		for(int i=0; i<n_cols; i++) {
			for(int j=0; j<n_rows; j++) {
				data[j][i] =  Double.parseDouble(input[j][i]);
			}			
		}
		
		long startTime = System.currentTimeMillis();
		FactorAnalysis fa = new FactorAnalysis(data, 2, false, false); 
		fa.do_factorAnalysis();
		long endTime = System.currentTimeMillis();
		System.out.println("FA done after " + ((endTime-startTime)/1000.0) + " secs.");
				
		System.out.println("");
		System.out.println("Rotation Matrix:");
		MatrixOperations.print_matrix(fa.get_rotation_matrix(0));
		System.out.println("");
		System.out.println("Factors:");
		MatrixOperations.print_matrix(fa.get_factors(0));
		
		//MatrixOperations.print_matrix(fa.get_my(0));
		System.out.println("Psi matrix:");
		MatrixOperations.print_matrix(fa.get_psi_matrix());
			
		double [][] prod = MatrixOperations.multiplication(fa.get_rotation_matrix(0), MatrixOperations.transpose(fa.get_rotation_matrix(0)));
		prod = MatrixOperations.add(fa.get_psi_matrix(), prod);
		MatrixOperations.print_matrix(prod);

	}
	
	
	public static void test2() {
		
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/DecisionTrees/Classification/IrisData.txt";
		InputDataManager inputData = new InputDataManager();	
		try {
			inputData.fileReader(dirName, false, true, true);
		} catch (Exception e) {
			System.out.println("Can´t find data set for upload.");
			e.printStackTrace();
		}
		
		String [][] input = inputData.strFileData;
		
		int n_rows = input.length;
		int n_cols = input[0].length-1;
		double [][] data = new double [n_rows][n_cols];
		for(int i=0; i<n_cols; i++) {
			for(int j=0; j<n_rows; j++) {
				data[j][i] =  Double.parseDouble(input[j][i]);
			}			
		}
		
		long startTime = System.currentTimeMillis();
		FactorAnalysis fa = new FactorAnalysis(data, 3, 1, true, true); 
		fa.do_factorAnalysis();
		long endTime = System.currentTimeMillis();
		System.out.println("FA done after " + ((endTime-startTime)/1000.0) + " secs.");
				
		System.out.println("");
		System.out.println("Rotation Matrix:");
		MatrixOperations.print_matrix(fa.get_rotation_matrix(0));
		System.out.println("");
		System.out.println("Factors:");
		MatrixOperations.print_matrix(fa.get_factors(0));
		
		//MatrixOperations.print_matrix(fa.get_my(0));
		System.out.println("Psi matrix:");
		MatrixOperations.print_matrix(fa.get_psi_matrix());
			
		double [][] prod = MatrixOperations.multiplication(fa.get_rotation_matrix(0), MatrixOperations.transpose(fa.get_rotation_matrix(0)));
		prod = MatrixOperations.add(fa.get_psi_matrix(), prod);
		MatrixOperations.print_matrix(prod);

	}
	
	
	public static void test3() {
		
		int nRandNumbers = 500;
		
		ArrayList<double [][]> randNumbers = new ArrayList<double [][]>();
		ArrayList<double [][]> ellipses = new ArrayList<double [][]>();
		
		double [][] mean = new double [2][1];
		mean [0][0] = 0.0;
		mean [1][0] = 0.0;
		
		double [][] Sigma = new double [2][2];
		Sigma[0][0] = 8.0;
		Sigma[1][0] = -2.2;
		Sigma[0][1] = Sigma[1][0];
		Sigma[1][1] = 4.0;
		
		NormalDistribution nDist = new NormalDistribution(mean, Sigma);
		double [][] org_randNumbers = MatrixOperations.transpose(nDist.sample(nRandNumbers));		
		double [][] org_ellipse = nDist.get_2d_confidence_ellipse(200, mean, Sigma, 0.9);
			
		long startTime = System.currentTimeMillis();
		FactorAnalysis4Test faTest = new FactorAnalysis4Test(org_randNumbers, 1, 1, false, false); 
		faTest.em_mixtureFactorAnalyzers4Test();
		long endTime = System.currentTimeMillis();
		System.out.println("FA done after " + ((endTime-startTime)/1000.0) + " secs.");
		
		int n = faTest.EMoutput.size();
		
		ArrayList<double [][]> W = new ArrayList<double [][]>();
		ArrayList<double [][]> EM_ellipses = new ArrayList<double [][]>();
		ArrayList<double [][]> reconstructions = new ArrayList<double [][]>();
		
		
		int [] iterations = {0, 2, (int) (n/4), n-1};
		int n_sel_iterations = iterations.length;
		
		for(int i=0; i<n_sel_iterations; i++) {
			
			randNumbers.add(org_randNumbers);
			ellipses.add(org_ellipse);
			
			int iter = iterations[i];
			W.add(faTest.EMoutput.get(iter).get("W").get(0));
			ArrayList<ArrayList<double[][]>> factors = faTest.EMfactors;
			reconstructions.add(MatrixOperations.multiplication(factors.get(iter).get(0),MatrixOperations.transpose(faTest.EMoutput.get(iter).get("W").get(0))));
			
			double [][] cov = MatrixOperations.multiplication(faTest.EMoutput.get(iter).get("W").get(0), MatrixOperations.transpose(faTest.EMoutput.get(iter).get("W").get(0)));
			cov = MatrixOperations.add(faTest.EMoutput.get(iter).get("Psi").get(0), cov);
			MatrixOperations.print_matrix(cov);
			
			EM_ellipses.add(nDist.get_2d_confidence_ellipse(200, faTest.EMoutput.get(iter).get("my").get(0), cov, 0.9));
		}
			
		plot_test3(randNumbers, ellipses, W, EM_ellipses, reconstructions, iterations);
		//plot_logLikelihoodOfEM(faTest.trace_logLikelihoods);
		//plot_variance_contributions(faTest.EMoutput, faTest.trace_logLikelihoods);
	}	

	
	public static void test4() {
		
		int n_analyzers    = 4;
		int nRandNumInCluster = 300;
		
		ArrayList<double [][]> randNumbers = new ArrayList<double [][]>();
		ArrayList<double [][]> ellipses = new ArrayList<double [][]>();
		ArrayList<double [][]> em_ellipses = new ArrayList<double [][]>();
		
	    double [] mean_x = {10, 20, 40, 20}; 
	    double [] mean_y = {20, 40, 20, 10};
		
	    int n_cluster = mean_x.length;
	    int nRandNumbers = nRandNumInCluster*n_cluster;
	    double [][] overall_randNumbers = new double [nRandNumbers][2];

		double [][] Sigma_org = new double [2][2];
		Sigma_org[0][0] = 8.0;
		Sigma_org[1][0] = -2.2;
		Sigma_org[0][1] = Sigma_org[1][0];
		Sigma_org[1][1] = 4.0;
	    
		int counter = 0;
		
		for(int i=0; i<n_cluster; i++) {
			double [][] mean = new double [2][1];
			mean [0][0] = mean_x[i];
			mean [1][0] = mean_y[i];
			
			double [][] Sigma = Sigma_org;
			Sigma[1][0] *= -1.0;
			Sigma[0][1] = Sigma[1][0];
			
			NormalDistribution nDist = new NormalDistribution(mean, Sigma);
			double [][] org_randNumbers = MatrixOperations.transpose(nDist.sample(nRandNumInCluster));		
			double [][] org_ellipse = nDist.get_2d_confidence_ellipse(200, mean, Sigma, 0.9);
			randNumbers.add(org_randNumbers);
			ellipses.add(org_ellipse);
			for(int j=0; j<nRandNumInCluster; j++) {
				for(int k=0; k<2; k++) {
					overall_randNumbers[counter][k] = org_randNumbers[j][k];
				}
				counter++;
			}		
		}
	    		
		long startTime = System.currentTimeMillis();
		FactorAnalysis4Test faTest = new FactorAnalysis4Test(overall_randNumbers, 1, n_analyzers, false, false); 
		faTest.em_mixtureFactorAnalyzers4Test();
		long endTime = System.currentTimeMillis();
		System.out.println("FA done after " + ((endTime-startTime)/1000.0) + " secs.");
			
		System.out.println(faTest.get_numberOfIterationsDone());
		//FactorAnalysis4Test faTest2 = new FactorAnalysis4Test(overall_randNumbers, 1, n_analyzers, false, false); 
		
		NormalDistribution nDist = new NormalDistribution(0.0, 0.01);
		ArrayList<double [][]> W_list = new ArrayList<double [][]>();
		for(int i=0; i<n_analyzers; i++) {
			double [][] W = new double [2][1];
			for(int j=0; j<2; j++) {
				W[j][0] = faTest.rotation_matrices.get(i)[j][0]+nDist.sample()[0][0];
			}
			W_list.add(W);
		}
		 
		faTest.set_external_init_pars(faTest.my, W_list, faTest.pi_k, faTest.Psi);
		faTest.em_mixtureFactorAnalyzers4Test();			
		System.out.println(faTest.get_numberOfIterationsDone());
		
		double [][] psi = faTest.get_psi_matrix();
		
		for(int i=0; i<n_analyzers; i++) {
			double [][] W = faTest.get_rotation_matrix(i);			
			double [][] my = faTest.get_my(i);
			double [][] cov = MatrixOperations.multiplication(W, MatrixOperations.transpose(W));
			cov = MatrixOperations.add(psi, cov);
			em_ellipses.add(nDist.get_2d_confidence_ellipse(200, my, cov, 0.9));
		}
			
		//plot_multiModalDist(overall_randNumbers); 
		
		//TODO: -> Show rotated W-axis over every cluster
		//      -> Show how clustering takes place over the EM steps
		//      -> Compare clustered (n_analyzers>1) with unclustered FA (n_analyzers = 1)
		plot_test4(randNumbers, ellipses, em_ellipses, faTest.rotation_matrices, faTest.my);
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_test3(ArrayList<double [][]> randNumbers, ArrayList<double [][]> ellipses, ArrayList<double [][]> W, ArrayList<double [][]> EM_ellipses, ArrayList<double [][]> reconstructions, int [] iterations) {
		
		GenGraphics graph = new GenGraphics();
		
		int nSamples = randNumbers.size();
      		
		String [] titles    = new String [nSamples];
	 	String [] subTitles = new String [nSamples];
	 	String [] yLabels   = new String [nSamples];;
	 	String [] xLabels   = new String [nSamples];;
	 	
	 	graph.setNumberOfPlotColums(2);
	 	graph.setNumberOfPlotRows(2);
	 	
	 	graph.setGraphWidth(600);
	 	graph.setGraphHeight(600);
	 	
	 	ArrayList<Color> defCols = graph.getDefaultColorsAsList();
	 	
	 	for(int i=0; i<nSamples; i++) {
	 		
		 	titles[i]    = "Factor Analysis";
		 	subTitles[i] = "    EM Iteration No." + iterations[i];
	 		xLabels[i]   = "x1";
	 		yLabels[i]   = "x2";
		 	
	 		boolean newPlot = true;
	 		
			int nRandNumbers = randNumbers.get(i).length;
	 		double [][] x1 = MatrixOperations.get_column_vec_from_matrix(randNumbers.get(i), 0);
			
			double x1_max = Utilities.Utilities.getMax(x1);
			double x1_min = Utilities.Utilities.getMin(x1);
			
			double [][] w1 = new double [2*nRandNumbers][1];
			double [][] w2 = new double [2*nRandNumbers][1];
			
			if(x1_max>=0.0 && x1_min>=0.0) {
				double s = x1_max/W.get(i)[0][0];
				double step1 = W.get(i)[0][0]*s/nRandNumbers;
				double step2 = W.get(i)[1][0]*s/nRandNumbers;
				int n= nRandNumbers;
				for(int j=0; j<n; j++) {
					w1[j][0] = +step1*(j+1);
					w2[j][0] = +step2*(j+1);
				}
			}
			
			if(x1_max<=0.0 && x1_min<=0.0) {		
				double s = -x1_min/W.get(i)[0][0];
				double step1 = W.get(i)[0][0]*s/nRandNumbers;
				double step2 = W.get(i)[1][0]*s/nRandNumbers;
				int n= nRandNumbers;
				for(int j=0; j<n; j++) {
					w1[j][0] = -step1*(j+1);
					w2[j][0] = -step2*(j+1);
				}
			}
			
			if(x1_max>=0.0 && x1_min<=0.0) {
				double x4Scaling = x1_max;
				if(Math.abs(x1_min)>x1_max) {
					x4Scaling = -x1_min;
				}
				double s = x4Scaling/W.get(i)[0][0];
				double step1 = W.get(i)[0][0]*s/nRandNumbers;
				double step2 = W.get(i)[1][0]*s/nRandNumbers;
				int n= nRandNumbers;
				for(int j=0; j<n; j++) {
					w1[j][0] = -step1*(n-j);
					w2[j][0] = -step2*(n-j);
					w1[n+j][0] = step1*(j+1);
					w2[n+j][0] = step2*(j+1);
				}
			}
			
			double [][] x2 = MatrixOperations.get_column_vec_from_matrix(randNumbers.get(i), 1);
	 		graph.plotPoints(x1,x2, newPlot, defCols.get(0));
	 		newPlot = false;	 
	 		graph.plotLines(w1,w2, newPlot, new Color(44, 102, 230, 180));
	 		
	 		//Plot reconstruction
	 		x1 = MatrixOperations.get_column_vec_from_matrix(reconstructions.get(i), 0);
			x2 = MatrixOperations.get_column_vec_from_matrix(reconstructions.get(i), 1);
	 		graph.plotPoints(x1,x2, newPlot, Color.RED);
	 		
	 		//Plot org. confidence ellipse
			x1 = MatrixOperations.get_column_vec_from_matrix(ellipses.get(i), 0);
			x2 = MatrixOperations.get_column_vec_from_matrix(ellipses.get(i), 1);
	 		graph.plotLines(x1,x2, newPlot, Color.GREEN);
	 		
	 		//Plot fitted EM confidence ellipse
			x1 = MatrixOperations.get_column_vec_from_matrix(EM_ellipses.get(i), 0);
			x2 = MatrixOperations.get_column_vec_from_matrix(EM_ellipses.get(i), 1);
	 		graph.plotLines(x1,x2, newPlot, Color.RED);
	 		
	 	}
	 		
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "11");
	 	graph.setSubTitle1(subTitles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(0);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
	 		 	
	}
	
	
	//Test Probabilistic PCA (PPCA)
	public static void test5() {
		
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/DecisionTrees/Classification/IrisData.txt";
		InputDataManager inputData = new InputDataManager();	
		try {
			inputData.fileReader(dirName, false, true, true);
		} catch (Exception e) {
			System.out.println("Can´t find data set for upload.");
			e.printStackTrace();
		}
		
		String [][] input = inputData.strFileData;
		
		int n_rows = input.length;
		int n_cols = input[0].length-1;
		double [][] data = new double [n_rows][n_cols];
		for(int i=0; i<n_cols; i++) {
			for(int j=0; j<n_rows; j++) {
				data[j][i] =  Double.parseDouble(input[j][i]);
			}			
		}
	
		long startTime = System.currentTimeMillis();
		PPCA ppca = new PPCA(data, 4, 1, false, false); 
		ppca.do_PPCA();
		long endTime = System.currentTimeMillis();
		System.out.println("PPCA done after " + ((endTime-startTime)/1000.0) + " secs.");
				
		System.out.println("");
		System.out.println("Rotation Matrix:");
		MatrixOperations.print_matrix(ppca.get_rotation_matrix(0));
		System.out.println("");
		//System.out.println("Factors:");
		//MatrixOperations.print_matrix(ppca.get_factors(0));
		
		//MatrixOperations.print_matrix(fa.get_my(0));
		System.out.println("Psi matrix:");
		MatrixOperations.print_matrix(ppca.get_psi_matrix());
		
		System.out.println("");
		System.out.println("Est. covariance (W^T W + Psi):");
		double [][] prod = MatrixOperations.multiplication(ppca.get_rotation_matrix(0), MatrixOperations.transpose(ppca.get_rotation_matrix(0)));
		prod = MatrixOperations.add(ppca.get_psi_matrix(), prod);
		MatrixOperations.print_matrix(prod);
		
		System.out.println("");
		System.out.println("Org. covariance:");
		MatrixOperations.print_matrix(GeneralMath.cov(data));		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_variance_contributions(ArrayList<HashMap<String, ArrayList<double [][]>>> pars, ArrayList<Double> logLiks) {
		
		GenGraphics graph = new GenGraphics();
	
	 	graph.setNumberOfPlotColums(3);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(600);
	 	graph.setGraphHeight(400);
	 	
	 	int n_variables = pars.get(0).get("W").get(0).length;
	 	int n_iterations = pars.size();
	 	
	 	String [] titles    = new String [n_variables+1];
	 	String [] subTitles = new String [n_variables+1];
	 	String [] yLabels   = new String [n_variables+1];
	 	String [] xLabels   = new String [n_variables+1];
	 		
	 	for(int l=0; l<n_variables; l++) {
	 		titles[l] = "  Variance Decomposition";
	 		subTitles[l] = "Var(x"+(l+1)+")";
	 		yLabels[l] = "x"+(l+1);
	 		xLabels[l] = "           EM Iteration No.";
		 	
	 		double [][] x = new double [n_iterations-1][1];
		 	double [][] factorContribution = new double [n_iterations-1][1];
		 	double [][] psiContribution = new double [n_iterations-1][1];
		 	double [][] totalVariance   = new double [n_iterations-1][1];		 	
		 	
		 	for(int i=1; i<n_iterations; i++) {
		 		int idx = i-1;
		 		x[idx][0] = i;
		 		factorContribution[idx][0] = Math.pow(pars.get(i).get("W").get(0)[l][0],2.0);
		 		psiContribution[idx][0] = pars.get(i).get("Psi").get(0)[l][l];
		 		totalVariance[idx][0] = factorContribution[idx][0]+psiContribution[idx][0];		 		
		 	}
		 	
	 		graph.plotLines(x,totalVariance, true, Color.GRAY);
	 		graph.plotLines(x,psiContribution, false, Color.GREEN);
	 		graph.plotLines(x,factorContribution, false, Color.RED);
	 	}

	 	n_iterations -= 1;
	 	double [][] logLik  = new double [n_iterations][1];
	 	double [][] x = new double [n_iterations][1];
	 	for(int i=0; i<n_iterations; i++) {
	 		logLik[i][0] = logLiks.get(i);
	 		x[i][0] = (i+1);
	 	}
	 	
 		titles[n_variables] = " Expected log Likelihood";
 		subTitles[n_variables] = "";
 		yLabels[n_variables] = "Q";
 		xLabels[n_variables] = "EM Iteration No.";
	 	graph.plotLines(x,logLik, true, new Color(44, 102, 230, 180));
	 	
 		
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "10");
	 	graph.setSubTitle1(subTitles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);

	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
		
	}
	
		
	@SuppressWarnings("static-access")
	public static void plot_logLikelihoodOfEM(ArrayList<Double> trace_logLik) {
		
		GenGraphics graph = new GenGraphics();
		
		String [] titles    = {"Trace EM Log Likelihood"};
	 	String [] yLabels   = {"Log Likelihood"};
	 	String [] xLabels   = {"Iteration No."};
	 	
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(400);
	 	graph.setGraphHeight(400);
	 	
	 	int n_iterations = trace_logLik.size();
	 	double [][] x = new double [n_iterations][1];
	 	double [][] logLik = new double [n_iterations][1];
	 	for(int i=0; i<n_iterations; i++) {
	 		x[i][0] = i+1;
	 		logLik[i][0] = trace_logLik.get(i);
	 	}
 		graph.plotLines(x,logLik, true, Color.RED);
	 		
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(0);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);

	 	graph.set_line_widht(2);
	 	
	 	graph.plot();
	 	
	}

	
	@SuppressWarnings("static-access")
	public static void plot_multiModalDist(double [][] multModalDist) {
		
		GenGraphics graph = new GenGraphics();
		
		String [] titles    = {"Multimodal Gaussian"};
	 	String [] yLabels   = {"x2"};
	 	String [] xLabels   = {"x1."};
	 	
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(400);
	 	graph.setGraphHeight(400);
	 	
	 	double [][] x1 = MatrixOperations.get_column_vec_from_matrix(multModalDist, 0);
	 	double [][] x2 = MatrixOperations.get_column_vec_from_matrix(multModalDist, 1);
 		graph.plotPoints(x1, x2, true, Color.RED);
	 		
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(0);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);

	 	graph.set_line_widht(2);
	 	
	 	graph.plot();
	 	
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_test4(ArrayList<double [][]> randNumbers, ArrayList<double [][]> org_ellipses, ArrayList<double [][]> em_ellipses, ArrayList<double [][]> W, ArrayList<double [][]> my) {
		
		int n_samples = randNumbers.size();
		
		GenGraphics graph = new GenGraphics();
		
		String [] title    = {"Mixture of Factor Analyzers"};
		String [] subTitle  = {"Multimodal Gaussian"};
	 	String [] yLabels   = {"x2"};
	 	String [] xLabels   = {"x1."};
	 	
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(500);
	 	graph.setGraphHeight(300);
	 	
	 	ArrayList<Color> defCols = graph.getDefaultColorsAsList();
	 	
	 	boolean newPlot = true;
	 	
	 	for(int i=0; i<n_samples; i++) {
	 		double [][] x1 = MatrixOperations.get_column_vec_from_matrix(randNumbers.get(i), 0);
		 	double [][] x2 = MatrixOperations.get_column_vec_from_matrix(randNumbers.get(i), 1);
	 		graph.plotPoints(x1, x2, newPlot, defCols.get(i));
	 		
	 		newPlot = false;
	 		
	 		x1 = MatrixOperations.get_column_vec_from_matrix(org_ellipses.get(i), 0);
		 	x2 = MatrixOperations.get_column_vec_from_matrix(org_ellipses.get(i), 1);
	 		graph.plotLines(x1, x2, newPlot, defCols.get(i)); 		
	 	}
	 	
	 	int n_analyzers = em_ellipses.size();
	 	
	 	for(int i=0; i<n_analyzers; i++) {	 	 			 		
	 		double [][] x1 = MatrixOperations.get_column_vec_from_matrix(em_ellipses.get(i), 0);
	 		double [][] x2 = MatrixOperations.get_column_vec_from_matrix(em_ellipses.get(i), 1);
	 		graph.plotLines(x1, x2, newPlot, Color.RED); 		
	 	}
	 	 	
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(title, null, "10");
	 	graph.setSubTitle1(subTitle, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(0);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);

	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
	 	
	}
	
	
	public static void main(String[] args) throws Exception {		
		//test1();
		//test2();
		//test3();
		//test4();
		test5();
	}
	
}
