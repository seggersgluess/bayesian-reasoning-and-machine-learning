package Clustering;

import java.awt.Color;
import java.util.ArrayList;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class TestOPTICS {


	public static void test_OPTICS() {
		
		int nRandNumInCluster = 250;
		
		ArrayList<double [][]> randNumbers = new ArrayList<double [][]>();

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
			randNumbers.add(org_randNumbers);
			for(int j=0; j<nRandNumInCluster; j++) {
				for(int k=0; k<2; k++) {
					overall_randNumbers[counter][k] = org_randNumbers[j][k];
				}
				counter++;
			}		
		}
	    		
		long startTime = System.currentTimeMillis();
		Optics optics = new Optics(overall_randNumbers, 2.5, 10);
		//optics.set_metric("chebyshev");
		optics.do_OPTICS();
		long endTime = System.currentTimeMillis();
		System.out.println("Clustering with OPTICS done after " + ((endTime-startTime)/1000.0) + " secs.");
		n_cluster = optics.get_number_of_clusters();
		ArrayList<double [][]> clusterData = new ArrayList<double [][]>(); 
		for(int i=0; i<n_cluster; i++) {
			clusterData.add(optics.get_X_of_cluster(i));
		}
		ArrayList<double [][]> noise = new ArrayList<double [][]>();
		noise.add(optics.get_X_of_noise());
		
		//optics.reachability_plot();
		
		plot_clusterData(overall_randNumbers, clusterData, noise); 
		
	}
	
	
	@SuppressWarnings("static-access")
	public static void plot_clusterData(double [][] orgData, ArrayList<double [][]> clusterData, ArrayList<double [][]> noise) {
		
		GenGraphics graph = new GenGraphics();
		
		String [] titles    = new String [2];
	 	String [] yLabels   = new String [2];
	 	String [] xLabels   = new String [2];
	 	
	 	titles[0] = "Original Data";
	 	titles[1] = "      OPTICS Clustered Data";
	 	yLabels[0] = "x2";
	 	yLabels[1] = yLabels[0];
	 	xLabels[0] = "x1";
	 	xLabels[1] = xLabels[0];
	 	
	 	graph.setNumberOfPlotColums(2);
	 	graph.setNumberOfPlotRows(1);
	 	
	 	graph.setGraphWidth(600);
	 	graph.setGraphHeight(600);
	 	
	 	double [][] x1 = MatrixOperations.get_column_vec_from_matrix(orgData, 0);
	 	double [][] x2 = MatrixOperations.get_column_vec_from_matrix(orgData, 1);
 		graph.plotPoints(x1, x2, true, Color.RED);
	 	
 		
 		ArrayList<Color> defColors = graph.getDefaultColorsAsList();
 		int n_clusters = clusterData.size();
 		boolean newPlot = true;
 		for(int i=0; i<n_clusters; i++) {
 		 	x1 = MatrixOperations.get_column_vec_from_matrix(clusterData.get(i), 0);
 		 	x2 = MatrixOperations.get_column_vec_from_matrix(clusterData.get(i), 1);
 	 		graph.plotPoints(x1, x2, newPlot, defColors.get(i));
 	 		newPlot = false;
 		}
 		
 		if(noise != null && noise.get(0) != null) {
 		 	x1 = MatrixOperations.get_column_vec_from_matrix(noise.get(0), 0);
 		 	x2 = MatrixOperations.get_column_vec_from_matrix(noise.get(0), 1);
 	 		graph.plotPoints(x1, x2, newPlot, Color.gray);
 		}
 			 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(0);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.set_point_width(4);	 	
	 	graph.plot();
	 	
	}
	
	
	public static void main(String[] args) throws Exception {	
		test_OPTICS();
	}
	
}



