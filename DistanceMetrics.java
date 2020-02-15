package Mathematics;

import java.util.HashMap;

public class DistanceMetrics {

	public double [][] x;
	public double [][] y;
	
	HashMap<String, double [][]> additionalInputPars = new HashMap<String, double [][]>();
	
	
	public DistanceMetrics(double [][] x, double [][] y) {
		this.x = x;
		this.y = y;
	}
	
	
	public double euclidian() {
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
		
		double d = 0.0;
		
		for(int i=0; i<n; i++) {
			d += (x[i][0]-y[i][0])*(x[i][0]-y[i][0]);
		}
		
		d = Math.sqrt(d);
		
		return d;
	}
	
	
	public double seuclidian() {
		
		if(additionalInputPars == null) {
			throw new RuntimeException("No additional input parameters for metric set.");
		}
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
		
		double [][] s = additionalInputPars.get("s");
		
		if(n != s.length) {
			throw new RuntimeException("Incorrect number of elements in supplied variances.");
		}
		
		double d = 0.0;
		
		for(int i=0; i<n; i++) {
			d += (x[i][0]-y[i][0])*(x[i][0]-y[i][0])/s[i][0];
		}

		d = Math.sqrt(d);
		
		return d;
	}
	
	
	public double manhatten() {
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
		
		double d = 0.0;
		
		for(int i=0; i<n; i++) {
			d += Math.abs((x[i][0]-y[i][0]));
		}

		return d;
	}
	
	
	public double chebyshev() {
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
		
		double d = 0.0;
		
		for(int i=0; i<n; i++) {
			double absDiff = Math.abs((x[i][0]-y[i][0]));
			if(absDiff > d) {
				d = absDiff;
			}
		}

		return d;
	}
	
	
	public double minkowski() {
		
		if(additionalInputPars == null) {
			throw new RuntimeException("No additional input parameters for metric set.");
		}
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
		
		double p = additionalInputPars.get("p")[0][0];
		
		double d = 0.0;
		
		for(int i=0; i<n; i++) {
			d += Math.pow(Math.abs((x[i][0]-y[i][0])), p);
		}

		d = Math.pow(d, 1.0/p);
		
		return d;
	}
	
	
	public double wminkowski() {
		
		if(additionalInputPars == null) {
			throw new RuntimeException("No additional input parameters for metric set.");
		}
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
		
		double p      = additionalInputPars.get("p")[0][0];
		double [][] w = additionalInputPars.get("w");
		
		if(n != w.length) {
			throw new RuntimeException("Invalid number of elements in supplied vector of weights.");
		}
		
		double d = 0.0;
		
		for(int i=0; i<n; i++) {
			d += Math.pow(Math.abs(w[i][0]*(x[i][0]-y[i][0])), p);
		}

		d = Math.pow(d, 1.0/p);
		
		return d;
	}
	
	
	public double mahalanobis() {
		
		if(additionalInputPars == null) {
			throw new RuntimeException("No additional input parameters for metric set.");
		}
		
		int n = x.length;
		
		if(n != y.length) {
			throw new RuntimeException("Unequal length of supplied vectors x and y.");
		}
			
		double [][] V = additionalInputPars.get("V");
		
		if(n != V.length || n != V[0].length) {
			throw new RuntimeException("Invalid number of elements in supplied vector of weights.");
		}
		
		double d = 0.0;
		
		double [][] V_inv = MatrixOperations.inverse(V);
		double [][] diff  = MatrixOperations.substract(x, y);
		double [][] diff_t = MatrixOperations.transpose(diff);
		
		d = MatrixOperations.multiplication(MatrixOperations.multiplication(diff_t, V_inv),diff)[0][0];
		
		d = Math.sqrt(d);
		
		return d;
	}
	
	
	public void setInput4Minkowski(double p) {
		
		double [][] p_mod = new double [1][1];
		p_mod[0][0] = p;
		
		additionalInputPars.put("p", p_mod);
	}
	
	
	public void setInput4WeightedMinkowski(double p, double [][] w) {
		
		double [][] p_mod = new double [1][1];
		p_mod[0][0] = p;
		
		additionalInputPars.put("p", p_mod);
		additionalInputPars.put("w", p_mod);
	}
	
	
	public void setInput4Mahalanobis(double [][] V) {
		additionalInputPars.put("V", V);
	}
	
	
	public void setInput4StandardizedEuclidian(double [][] s) {
		additionalInputPars.put("s", s);
	}

	
	public double calcDistanceMetric(String metric) {
		
		metric = metric.toLowerCase();
		
		String [] validMetrics = getListOfDistanceMetrics();
		int [] validIdx = Utilities.Utilities.get_idx(validMetrics, metric);
		
		if(validIdx[0] == -1) {
			throw new RuntimeException(metric + " is not a valid distance metric.");
		}
				
		double d = 0.0;
		
		if(metric.contentEquals("euclidian")) {
			d = euclidian(); 
		}
		
		if(metric.contentEquals("manhatten")) {
			d = manhatten(); 
		}
		
		if(metric.contentEquals("chebyshev")) {
			d = chebyshev(); 
		}
		
		if(metric.contentEquals("seuclidean")) {
			d = seuclidian(); 
		}
		
		if(metric.contentEquals("mahalanobis")) {
			d = mahalanobis(); 
		}
		
		if(metric.contentEquals("minkowski")) {
			d = minkowski(); 
		}
		
		if(metric.contentEquals("wminkowski")) {
			d = wminkowski(); 
		}
		
		return d;
	}

	
	public double calcDistanceMetric(String metric, double p) {
		
		setInput4Minkowski(p);
		double d = calcDistanceMetric(metric);
		
		return d;
	}
	
	
	public double calcDistanceMetric(String metric, double p, double [][] w) {
		
		setInput4WeightedMinkowski(p,w);
		double d = calcDistanceMetric(metric);
		
		return d;
	}
	
	
	public double calcDistanceMetric(String metric, double [][] w) {
		
		String [] validMetrics = getListOfDistanceMetrics();
		int [] validIdx = Utilities.Utilities.get_idx(validMetrics, metric);
		if(validIdx[0] == -1) {
			throw new RuntimeException(metric + " is not a valid distance metric.");
		}
		
		if(metric.contentEquals("seuclidian")) {
			setInput4StandardizedEuclidian(w);
		}
		
		if(metric.contentEquals("mahalanobis")) {
			setInput4Mahalanobis(w);
		}
		
		double d = calcDistanceMetric(metric);
		
		return d;
	}
	
	
	public void setVectors(double [][] x, double [][] y) {
		this.x = x;
		this.y = y;
	}
	
	
	public String [] getListOfDistanceMetrics() {
		
		String [] metrics = new String [7];
		
		metrics[0] = "euclidian";
		metrics[1] = "manhatten";
		metrics[2] = "chebyshev";
		metrics[3] = "minkowski";
		metrics[4] = "wminkowski";
		metrics[5] = "seuclidean";
		metrics[6] = "mahalanobis";
		
		return metrics;
	}
	
}
