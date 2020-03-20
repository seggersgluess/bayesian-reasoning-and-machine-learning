package ComponentModels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class TestICA {

	
	public static void test1() {
		
		//Plot: obs, ica, pca		
		String plotData = "pca";
		
		ArrayList<double [][]> org_signals = get_orginial_signals();
		HashMap<String, ArrayList<double [][]>> obs_signalList = get_observed_signals(org_signals);
		ArrayList<double [][]> obs_signals = obs_signalList.get("list");
		double [][] obs_signal_X = obs_signalList.get("array").get(0);
		
		ICA ica = new ICA(obs_signal_X, 4, true, true);
		ica.do_ICA();
		double [][] icaReconstruction = ica.calc_uncentered_and_unwhitened_rotated_X();
		double [][] icaFactors = ica.get_factors();
		
		int n_variables = icaReconstruction[0].length;
		ArrayList<double [][]> icaReconsList = new ArrayList<double [][]>();
		ArrayList<double [][]> icaFactorList = new ArrayList<double [][]>();
		for(int i=0; i<n_variables; i++) {
			icaReconsList.add(MatrixOperations.get_column_vec_from_matrix(icaReconstruction,i));
			icaFactorList.add(MatrixOperations.get_column_vec_from_matrix(icaFactors,i));
		}
	
		PCA pca = new PCA(obs_signal_X, 4, false, false);
		pca.do_PCA();
		double [][] pcaReconstruction = pca.get_rotated_input();
		double [][] pcaFactors = pca.get_factors();

		n_variables = pcaReconstruction[0].length;
		ArrayList<double [][]> pcaReconsList = new ArrayList<double [][]>();
		ArrayList<double [][]> pcaFactorList = new ArrayList<double [][]>();
		for(int i=0; i<n_variables; i++) {
			pcaReconsList.add(MatrixOperations.get_column_vec_from_matrix(pcaReconstruction,i));
			pcaFactorList.add(MatrixOperations.get_column_vec_from_matrix(pcaFactors,i));
		}
	
		if(plotData.contentEquals("ica")) {
			plot_ica(icaFactorList,icaReconsList);
		}
		if(plotData.contentEquals("pca")) {
			plot_pca(pcaFactorList,pcaReconsList);
		}
		if(plotData.contentEquals("obs")) {
			plotSignals(org_signals, obs_signals, null);
		}
	}
	
	
	@SuppressWarnings("static-access")
	public static void test2() {
		
		int n=100;
		int max = 10;
		int min = -10;
		double diff = (max-min)/((double)(n));
		
		ArrayList<double [][]> org_signals = get_orginial_signals();
		HashMap<String, ArrayList<double [][]>> obs_signalList = get_observed_signals(org_signals);
		double [][] obs_signal_X = obs_signalList.get("array").get(0);
		
		double [][] x = new double [n][1];
		double [][] dist = new double [n][1];
		
		ICA ica = new ICA(obs_signal_X, 4, true, true);
		for(int i=0; i<n; i++) {
			x[i][0] = min+diff*i;
			dist[i][0] = ica.g_superGaussian(x[i][0], 0);
		}
		
		plotDistribution4ICA(x, dist, "Sub Gaussian");
		
	}
	

	public static HashMap<String, ArrayList<double [][]>> get_observed_signals(ArrayList<double [][]> org_signals) {
		
		HashMap<String, ArrayList<double [][]>> obs_signalList = new HashMap<String, ArrayList<double [][]>>();
		
		int n_vars = org_signals.size();
		int sampleLength = org_signals.get(0).length;
		
		double [][] mixing_matrix = new double [4][4];
		for(int i=0; i<n_vars; i++) {
			for(int j=0; j<n_vars; j++) {
				if(i==j) {
					mixing_matrix[i][j] = 1.5;
				}else {			
					if(i==2) {
						//Weighting of the white noise term
						mixing_matrix[i][j] = 1.2;
					}else {
						mixing_matrix[i][j] = 0.8;
					}														
				}				
			}
		}
		
		double [][] obs_X = new double [sampleLength][n_vars];
		for(int i=0; i<n_vars; i++) {
			for(int j=0; j<sampleLength; j++) {
				obs_X[j][i] = org_signals.get(i)[j][0];
			}
		}
		
		obs_X = MatrixOperations.multiplication(obs_X, mixing_matrix);
		
		ArrayList<double [][]> obs_signals = new ArrayList<double [][]>();
		for(int i=0; i<n_vars; i++) {
			obs_signals.add(MatrixOperations.get_column_vec_from_matrix(obs_X, i));
		}
		
		ArrayList<double [][]> arrayList = new ArrayList<double [][]>();
		arrayList.add(obs_X);
		
		obs_signalList.put("list", obs_signals);
		obs_signalList.put("array", arrayList);		
		
		return obs_signalList;
	}
	
	
	public static ArrayList<double [][]> get_orginial_signals() {
		
		int signalLength = 501;
		
		ArrayList<double [][]> org_signals = new ArrayList<double [][]>();
		
		double [][] sawtoothWave = sawtoothWave(signalLength, 25, 4.0, -2.0);		
		org_signals.add(sawtoothWave);
		
		double [][] tangWave = tangensWave(signalLength, 20, 4.0, 0.0);
		org_signals.add(tangWave);
		
		NormalDistribution nDist = new NormalDistribution(0.0, 2.0);
		double [][] whiteNoise = nDist.sample(signalLength);
		org_signals.add(whiteNoise);
		
		double [][] sinusWave = sinusWave(signalLength, 40, 1.0, 0.0);
		org_signals.add(sinusWave);
		
		return org_signals;
	}
	
	
	public static double [][] sawtoothWave(int length, int period, double amplitude, double yShift) {
		double [][] sawtoothWave = new double [length][1];
		double stepFrac = amplitude/((double) period);
		int counter = 0;
		for(int i=0; i<length; i++) {					
			sawtoothWave[i][0] = (counter*stepFrac)+yShift;
			if(counter == period) {
				counter = 0;
			}
			counter++;
		}
		return sawtoothWave;
	}

	
	public static double [][] tangensWave(int length, int period, double amplitude, double yShift) {
		double min = -0.45*Math.PI;
		double max = 0.45*Math.PI;
		double s   = Math.tan(max);
		double diff = max-min;
		double step = diff/((double) period);
		double [][] tangensWave = new double [length][1];
		int counter = 0;
		for(int i=0; i<length; i++) {
			double x = min+counter*step;
			tangensWave[i][0] = Math.tan(x)/s*amplitude+yShift;		
			if(counter == period) {
				counter = -1;
			}
			counter++;
		}
		return tangensWave;
	}
	
	
	public static double [][] sinusWave(int length, int period, double amplitude, double yShift) {
		double min = 0.0;
		double max = 2.0*Math.PI;
		double step = max/((double) period);
		double [][] sinusWave = new double [length][1];
		int counter = 0;
		for(int i=0; i<length; i++) {
			double x = min+counter*step;
			sinusWave[i][0] = Math.sin(x)*amplitude+yShift;		
			if(counter == period) {
				counter = 0;
			}
			counter++;
		}
		return sinusWave;
	}
	
	
	public static void plot_ica(ArrayList<double [][]> factors, ArrayList<double [][]> reconstruction) {
		
		ArrayList<String []> titleList = new ArrayList<String []>();
		String [] titles = {"ICA", "ICA"};
		String [] subTitles = {"Extracted Factors", "Reconstructed Observations"};
		
		titleList.add(titles);
		titleList.add(subTitles);
		
		plotSignals(factors, reconstruction, titleList);	
	}
	
	
	public static void plot_pca(ArrayList<double [][]> factors, ArrayList<double [][]> reconstruction) {
		
		ArrayList<String []> titleList = new ArrayList<String []>();
		String [] titles = {"PCA", "PCA"};
		String [] subTitles = {"Extracted Factors", "Reconstructed Observations"};
		
		titleList.add(titles);
		titleList.add(subTitles);
		
		plotSignals(factors, reconstruction, titleList);
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotSignals(ArrayList<double [][]> orgSignals, ArrayList<double [][]> obsSignals, ArrayList<String []> titleList) {
		
		GenGraphics graph = new GenGraphics();
		
		int nSignals = orgSignals.size();
		int nSamples = 2*nSignals;
      		
	 	String [] yLabels   = new String [nSamples];
	 	String [] xLabels   = new String [nSamples];
	 	
	 	graph.setNumberOfPlotColums(2);
	 	graph.setNumberOfPlotRows(nSignals);
	 	
	 	graph.setGraphWidth(700);
	 	graph.setGraphHeight(500);
	 	List<Color> defColors = graph.getDefaultColorsAsList();
	 	
	 	int n = orgSignals.get(0).length;
	 	double [][] x = graph.get_default_x_axis_labels(n);
	 	int counter = 0;
	 	int idx = 0;
	 	for(int i=0; i<nSamples; i++) {
	 		
	 		xLabels[i]   = "";
	 		yLabels[i]   = "";
		 	 
	 		double [][] y = new double [1][1];
	 		
	 		if(counter == 0) {
	 			y = orgSignals.get(idx);
	 		}
	 		if(counter == 1) {
	 			y = obsSignals.get(idx);
	 		}
	 				 		
	 		graph.plotLines(x, y, true, defColors.get(idx));
	 		
	 		if(counter == 1) {
	 			idx++;
	 			counter = -1;
	 		}
	 		
	 		counter++;
	 	}
	 	
	 	String [] titles = new String [2];
	 	String [] subTitles = new String [2];
	 	
	 	if(titleList == null) {
	 		titles[0] = "Source Signals";
		 	titles[1] = "Observed Signals";
	 	}else {
	 		if(titleList.get(0) == null) {
		 		titles[0] = "Source Signals";
			 	titles[1] = "Observed Signals";
	 		}else {
		 		titles[0] = titleList.get(0)[0];
			 	titles[1] = titleList.get(0)[1];
	 		}
	 	}
	 	if(titleList == null) {
	 		subTitles[0] = "";
	 		subTitles[1] = "";	
	 	}else {
	 		if(titleList.get(1) == null) {
		 		subTitles[0] = "";
		 		subTitles[1] = "";	
	 		}else {
		 		subTitles[0] = titleList.get(1)[0];
		 		subTitles[1] = titleList.get(1)[1];	
	 		}
	 	}
	 	
	 	
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "11");
	 	graph.setSubTitle1(subTitles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
	 		 	
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotDistribution4ICA(double [][] x, double [][] distValues, String distribution) {
		
		GenGraphics graph = new GenGraphics();		
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1); 	
	 	graph.setGraphWidth(700);
	 	graph.setGraphHeight(500);
	 	
	 	graph.plotLines(x, distValues, true, Color.BLUE);
	 		
	 	String titles [] = new String [1];
	 	String yLabels [] = new String [1];
	 	String xLabels [] = new String [1];
	 	titles[0] = distribution;
	 	yLabels[0] = distribution;
	 	xLabels[0] = "";
	 			
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "11");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
	 		 	
	}
	
	
	public static void main(String[] args) throws Exception {		
		//test1();
		test2();
	}
	
}
