package Clustering;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import Graphics.GenGraphics;
import Mathematics.DistanceMetrics;
import Mathematics.MatrixOperations;

public class Optics extends Clustering{

	ArrayList<HashMap<String, double [][]>> orderedList;
	
	
	public Optics(double [][] X, double epsilon, int minPoints) {		
		super(X,epsilon,minPoints);		
		convertInputIntoNeededStructure();	
	}
	
	
	public void do_OPTICS() {
		expandClusterOrder();
		extract_OPTICS_clustering();
	}
	
	
	public void convertInputIntoNeededStructure() {
		
		inputStructure = new ArrayList<HashMap<String, double [][]>>();
		
		for(int i=0; i<n_observations; i++) {
			HashMap<String, double [][]> struct = new HashMap<String, double [][]>();
			double [][] x = MatrixOperations.get_row_vec_from_matrix(X, i);
			double [][] idx = new double [1][1];
			double [][] p = new double [1][1];
			idx[0][0] = i;
			p[0][0] = 0.0;
			struct.put("features", x);
			struct.put("idx", idx);
			struct.put("coreDist", null);
			struct.put("reachDist", null);
			struct.put("processed", p);
			struct.put("clusterID", null);
			inputStructure.add(struct);
		}	
	}
	
	
	public void expandClusterOrder() {
		
		orderedList = new ArrayList<HashMap<String, double [][]>>();
		
		for(int i=0; i<n_observations; i++) {
			HashMap<String, double [][]> point = inputStructure.get(i);
			if(point.get("processed")[0][0] == 0.0) {
				ArrayList<HashMap<String, double [][]>> N = get_neighborhood(point);
				markAsProcessed((int)(point.get("idx")[0][0]));			
			    double [][] cd = coreDistance(N, point);
			    point.put("coreDist", cd);
			    orderedList.add(point);
			    if(cd != null) {
			    	ArrayList<HashMap<String, double [][]>> seeds = new ArrayList<HashMap<String, double [][]>>();
			    	seeds = updateSeed(N, point, seeds);
			    	int nSeedPts = seeds.size();
			    	while(nSeedPts>0) {
				    	HashMap<String, double [][]> q_point = seeds.get(0); 
				    	seeds.remove(0);
				    	ArrayList<HashMap<String, double [][]>> N_q = get_neighborhood(q_point);
						markAsProcessed((int)(q_point.get("idx")[0][0]));						
						cd = coreDistance(N_q, q_point);						
						q_point.put("coreDist", cd);
						orderedList.add(q_point);
						if(cd != null) {
							seeds = updateSeed(N_q, q_point, seeds);
						}
						nSeedPts = seeds.size();
			    	}
			    }
			}
		}
		
	}
	
	
	public ArrayList<HashMap<String, double [][]>> updateSeed(ArrayList<HashMap<String, double [][]>> N, HashMap<String, double [][]> point, ArrayList<HashMap<String, double [][]>> seeds) {
		
		double [][] cd = coreDistance(N, point);
		if(cd != null) {
			int n = N.size();
			for(int i=0; i<n; i++) {
				HashMap<String, double [][]> selPoint = N.get(i);
				if(selPoint.get("processed")[0][0] == 0.0) {					
					double [][] reachDist = reachabilityDistance(cd[0][0], point, selPoint);
					if(selPoint.get("reachDist") == null) {
						selPoint=setReachDist(selPoint, reachDist);
						seeds = insert2Seeds(selPoint, seeds);
					}else {
						if(selPoint.get("reachDist")[0][0]<reachDist[0][0]) {
							selPoint=setReachDist(selPoint, reachDist);
							seeds = resortSeeds(selPoint, seeds);
						}
					}
				}
			}
		}

		return seeds;
	}
	
	
	public HashMap<String, double [][]> setReachDist(HashMap<String, double [][]> point, double [][] reachDist) {
		point.put("reachDist", reachDist);
		inputStructure.get((int)(point.get("idx")[0][0])).put("reachDist", reachDist);
		return point;
	}
	
	
	public ArrayList<HashMap<String, double [][]>> insert2Seeds(HashMap<String, double [][]> point, ArrayList<HashMap<String, double [][]>> seeds) {
		
		int n = seeds.size();
		
		if(n==0) {
			seeds.add(point);
		}else {
			ArrayList<HashMap<String, double [][]>> newSeeds = new ArrayList<HashMap<String, double [][]>>();
			boolean inserted = false;
			for(int i=0; i<n; i++) {
				double reachDist = seeds.get(i).get("reachDist")[0][0];
				if(inserted == false) {
					if(point.get("reachDist")[0][0]<=reachDist) {
						newSeeds.add(point);
						inserted = true;
					}
				}
				newSeeds.add(seeds.get(i));					
			}	
			if(inserted == false) {
				newSeeds.add(point);
				inserted = true;
			}			
			seeds = newSeeds;
		}		
		return seeds;
	}
	

	public ArrayList<HashMap<String, double [][]>> resortSeeds(HashMap<String, double [][]> point, ArrayList<HashMap<String, double [][]>> seeds) {
		
		int n = seeds.size();
		
		if(n==0) {
			seeds.add(point);
		}else {
			ArrayList<HashMap<String, double [][]>> newSeeds = new ArrayList<HashMap<String, double [][]>>();
			boolean inserted = false;
			for(int i=0; i<n; i++) {
				double [][] idx = point.get("idx");
				double reachDist = seeds.get(i).get("reachDist")[0][0];
				if(inserted == false) {
					if(point.get("reachDist")[0][0]<=reachDist) {
						newSeeds.add(point);
						inserted = true;
					}
				}
				if(idx[0][0] != seeds.get(i).get("idx")[0][0]) {
					newSeeds.add(seeds.get(i));
				}									
			}
			if(inserted == false) {
				newSeeds.add(point);
			}
			seeds = newSeeds;
		}		
		return seeds;
	}
	
	
	public void markAsProcessed(int idx) {
		double [][] p = new double [1][1];
		p[0][0] = 1.0;
		inputStructure.get(idx).put("processed", p);
	}
		

	public double [][] coreDistance(ArrayList<HashMap<String, double [][]>> N, HashMap<String, double [][]> point) {
		
		double [][] coreDistance = new double [1][1];
		int n= N.size();
		if(n<minPts) {
			coreDistance = null;
		}else {
			double [][] x = point.get("features");
			ArrayList<Double> distList = new ArrayList<Double>(n);
			for(int i=0; i<n; i++) {
				double [][] y = N.get(i).get("features");
				distList.add(new DistanceMetrics(x,y).calcDistanceMetric(metric));
			}
			Collections.sort(distList);
			coreDistance[0][0] = distList.get(minPts-1);
		}
		return coreDistance;
	}
	
	
	public double [][] reachabilityDistance(double coreDist, HashMap<String, double [][]> point1, HashMap<String, double [][]> point2) {
		
		double [][] reachDist = new double [1][1];
		reachDist[0][0] = coreDist;
		double dist = new DistanceMetrics(point1.get("features"),point2.get("features")).calcDistanceMetric(metric);		
		if(dist>reachDist[0][0]) {
			reachDist[0][0] = dist;
		}
		return reachDist;	
	}
	
	
	//Identify/ extract clusters
	public void extract_OPTICS_clustering() {
		if(orderedList == null)	{
			if(orderedList.size() == 0) {
				throw new RuntimeException("OPTICS not run yet. Cannot create a reachability plot.");
			}
		}
		
		clusterMemberIdxs = new ArrayList<ArrayList<Integer>>();
		noiseIdxs = new ArrayList<Integer>();
		
		double e4Clustering = e;
		int n = orderedList.size();
		double [][] clusterID = new double [1][1];
		double c = 0.0;
		
		for(int i=0; i<n; i++) {	
			HashMap<String, double [][]> point = orderedList.get(i);
			if(point.get("reachDist") == null || point.get("reachDist")[0][0]>e4Clustering) {
				if(point.get("coreDist") == null || point.get("coreDist")[0][0]>e4Clustering) {
					clusterID[0][0] = -1.0;
					orderedList.get(i).put("clusterID", clusterID);
					noiseIdxs.add((int)(point.get("idx")[0][0]));
				}else {
					ArrayList<Integer> clusterIdxList = new ArrayList<Integer>();
					clusterID[0][0] = c;
					orderedList.get(i).put("clusterID", clusterID);
					clusterIdxList.add((int)(point.get("idx")[0][0]));
					clusterMemberIdxs.add(clusterIdxList);
					c++;
				}
			}else {
				orderedList.get(i).put("clusterID", clusterID);
				clusterMemberIdxs.get((int)(clusterID[0][0])).add((int)(point.get("idx")[0][0]));
			}			
		}
		n_clusters = (int) c;
	}
		
	
	@SuppressWarnings("static-access")
	public void reachability_plot() {
		if(orderedList == null)	{
			if(orderedList.size() == 0) {
				throw new RuntimeException("OPTICS not run yet. Cannot create a reachability plot.");
			}
		}
		
		double [][] reachDist = new double [n_observations][1];
		double [][] epsDist   = new double [n_observations][1];
		for(int i=0; i<n_observations; i++) {
			double [][] reachDistOfPoint = orderedList.get(i).get("reachDist");
			if(reachDistOfPoint == null) {
				reachDist[i][0] = 5.0;
			}else {
				reachDist[i][0] = reachDistOfPoint[0][0];
			}	
			epsDist[i][0] = e;
		}
				
	 	String yLabels [] = {"Reachability Distance"};
	 	String xLabels [] = {"Ordered Index"};
	 	String title [] = {"Reachability Distance"};
	 	
		GenGraphics graph = new GenGraphics();		
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1); 	
	 	graph.setGraphWidth(900);
	 	graph.setGraphHeight(600);
	 	
	 	double [][] x = graph.get_default_x_axis_labels(n_observations);

	 	graph.plotLines(x, reachDist, true, Color.RED);
	 	graph.plotLines(x, epsDist, false, Color.BLACK);

	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(title, null, "9");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.setNumberOfXDivisions(4);
	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
			
	}
	
	
	public int get_number_of_clusters() {
		if(orderedList == null)	{
			if(orderedList.size() == 0) {
				throw new RuntimeException("OPTICS not run yet. Cannot create a reachability plot.");
			}
		}
		return n_clusters;	
	}

	
}
