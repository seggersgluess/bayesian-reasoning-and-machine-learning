package Clustering;

import java.util.ArrayList;
import java.util.HashMap;

import Mathematics.MatrixOperations;

public class DBScan extends Clustering{

	public DBScan(double [][] X, double epsilon, int minPoints) {
		super(X,epsilon, minPoints);		
		convertInputIntoNeededStructure();
	}
	
	
	public void do_DBScan() {
		do_DBScan_clustering();
		extract_DBScan_clustering();
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
			struct.put("processed", p);
			struct.put("clusterID", null);
			inputStructure.add(struct);
		}	
	}
	
	
	public void do_DBScan_clustering() {
		double c = 0.0;
		for(int i=0; i<n_observations; i++) {
			HashMap<String, double [][]> point = inputStructure.get(i);
			if(point.get("processed")[0][0]==0.0) {
				ArrayList<HashMap<String, double [][]>> N = get_neighborhood(point);
				markAsProcessed((int)(point.get("idx")[0][0]));	
				if(N.size()<minPts) {
					point = setClusterIdx(point, -1.0);					
				}else {
					expandCluster(point, c, N);
					c++;
				}
			}		
		}	
		n_clusters = (int) c;
	}
	
	
	public void expandCluster(HashMap<String, double [][]> point, double clusterID, ArrayList<HashMap<String, double [][]>> N) {
		
		setClusterIdx(point, clusterID);
		for(int i=0; i<N.size(); i++) {
			HashMap<String, double [][]> q_Point = N.get(i);
			if(q_Point.get("processed")[0][0] == 0.0) {
				markAsProcessed((int)(q_Point.get("idx")[0][0]));
				ArrayList<HashMap<String, double [][]>> N_q = get_neighborhood(q_Point);
				if(N_q.size()>=minPts) {
					N=join_neighborhoods(N,N_q);
				}
			}
			if(q_Point.get("clusterID")==null) {
				q_Point = setClusterIdx(q_Point, clusterID);	
			}
		}
		
	}
	
	
	public ArrayList<HashMap<String, double [][]>> join_neighborhoods(ArrayList<HashMap<String, double [][]>> N, ArrayList<HashMap<String, double [][]>> N_q) {
		
		int n = N.size();
		int n_q = N_q.size();
		
		ArrayList<HashMap<String, double [][]>> joint_N = N;
		for(int i=0; i<n_q; i++) {
			double idx2 = N_q.get(i).get("idx")[0][0];
			boolean add = true;
			for(int j=0; j<n; j++) {
				double idx1 = N.get(j).get("idx")[0][0];
				if(idx1 == idx2) {
					add = false;
				}
			}
			if(add == true) {
				joint_N.add(N_q.get(i));
			}
		}
		return joint_N;
	}
	
	
	public void markAsProcessed(int idx) {
		double [][] p = new double [1][1];
		p[0][0] = 1.0;
		inputStructure.get(idx).put("processed", p);
	}
	
	
	public HashMap<String, double [][]> setClusterIdx(HashMap<String, double [][]> point, double clusterID) {
		double [][] id = new double [1][1];;
		id[0][0] = clusterID;
		point.put("clusterID",id);
		inputStructure.get((int)(point.get("idx")[0][0])).put("clusterID", id);
		return point;
	}
	
	
	public void extract_DBScan_clustering() {
		
		clusterMemberIdxs = new ArrayList<ArrayList<Integer>>();
		for(int i=0; i<n_clusters; i++) {
			ArrayList<Integer> list = new ArrayList<Integer>();
			clusterMemberIdxs.add(list);
		}
		noiseIdxs = new ArrayList<Integer>();
		
		for(int i=0; i<n_observations; i++) {
			int clusterID = (int) (inputStructure.get(i).get("clusterID")[0][0]);
			int pointIdx = (int)(inputStructure.get(i).get("idx")[0][0]);
			if(clusterID == -1.0) {
				noiseIdxs.add(pointIdx);
			}else {
				clusterMemberIdxs.get(clusterID).add(pointIdx);
			}
		}		
	}
	
	
}
