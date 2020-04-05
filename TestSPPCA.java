package ComponentModels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;

import DataManagement.InputDataManager;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;


public class TestSPPCA {

	
	//Test for classification (with Iris data)
	public static void test1_SPPCA() {
		
		HashMap<String, double [][]> iris_data = get_iris_data_4_classification_task();
		
		SPPCA sppca = new SPPCA(iris_data.get("X"), iris_data.get("y"), 4, false);
		HashMap<String, double [][]> init_pars = get_dummy_init_pars();
		sppca.set_external_init_pars(init_pars);
		sppca.do_SPPCA();
					
		ArrayList<double[][]> trained_data = new ArrayList<double[][]>();
		ArrayList<double[][]> org_data     = new ArrayList<double[][]>();
		
		double [][] X_trained = sppca.predict(iris_data.get("X")).get("X");//sppca.get_rotated_input();	
		double [][] y_trained = sppca.predict(iris_data.get("X")).get("Y");//sppca.get_rotated_Y();
				
		trained_data.add(y_trained);
		org_data.add(iris_data.get("y"));		
		int n = X_trained[0].length;
		for(int i=0; i<n; i++) {
			trained_data.add(MatrixOperations.get_column_vec_from_matrix(X_trained, i));
			org_data.add(MatrixOperations.get_column_vec_from_matrix(iris_data.get("X"), i));
		}
		
		String [] titles = {"Species", 
				            "Sepal Length",
				            "Sepal Width",
				            "Petal Length",
				            "Petal Width"};
		
		HashMap<String, Integer> config = new HashMap<String, Integer>();
		config.put("n_cols",1);
		config.put("n_rows",n+1);
		
		plotSPPCATrainedData(org_data, trained_data, titles, config);
		System.out.println(sppca.get_sigma_y());
		System.out.println(sppca.get_sigma_x());
	}
	

	//Test for regression (with Boston housing data)
	public static void test2_SPPCA() {
		
		String target = "MEDV";
		
		HashMap<String, Object [][]> boston_data = get_boston_housing_data_4_regression_task();
		
		String [][] featNames = (String[][]) boston_data.get("var_names");
		int n_vars = featNames.length-1;
		int n_obs  = boston_data.get("X").length;
		int idx = Utilities.Utilities.get_idx(featNames, target)[0];
		
		double [][] X = new double [n_obs][n_vars];
		double [][] y = new double [n_obs][1];
		String [] titles = new String [(n_vars+1)];
		titles[0] = target;		
		for(int i=0; i<n_obs; i++) {
			int c = 0;
			for(int j=0; j<n_vars; j++) {
				if(j!=idx) {
					X[i][c] = (double) boston_data.get("X")[i][j];
					titles[(c+1)] = featNames[j][0];
					c++;
				}
			}
			y[i][0] = (double) boston_data.get("X")[i][idx];
		}
		
		SPPCA sppca = new SPPCA(X, y, 13, false);
		sppca.do_SPPCA();
		
		ArrayList<double[][]> trained_data = new ArrayList<double[][]>();
		ArrayList<double[][]> org_data     = new ArrayList<double[][]>();
		
		double [][] X_trained = sppca.get_rotated_input();	
		double [][] y_trained = sppca.get_rotated_Y();
		
		trained_data.add(y_trained);
		org_data.add(y);
		
		int n = X_trained[0].length;
		for(int i=0; i<n; i++) {
			trained_data.add(MatrixOperations.get_column_vec_from_matrix(X_trained, i));
			org_data.add(MatrixOperations.get_column_vec_from_matrix(X, i));
		}
			
		HashMap<String, Integer> config = new HashMap<String, Integer>();
		config.put("n_cols",7);
		config.put("n_rows",2);
		
		plotSPPCATrainedData(org_data, trained_data, titles, config);
	}
	
	
	//Test for classification (with Iris data)
	public static void test1_SSPPCA() {
		
		HashMap<String, double [][]> iris_data = get_iris_data_4_classification_task();
		
		double labeledFrac = 0.8;
		int n_labeled = (int) (50*labeledFrac);
		int n_unlabeled = (int) (50.0-n_labeled);
		
		int n_features = iris_data.get("X")[0].length;
		
		int n_labeled_data = 3*n_labeled;
		int n_unlabeled_data = 3*n_unlabeled;
		
		double [][] labeled_X = new double [n_labeled_data][n_features];
		double [][] unlabeled_X = new double [n_unlabeled_data][n_features];
		double [][] y = new double [n_labeled_data][1];
		double [][] unlabeled_y = new double [n_unlabeled_data][1];
		
		int c = 0;
		int s = 0;
		//Labeled data X_1 & y
		for(int i=0; i<n_labeled_data; i++) {
			for(int j=0; j<n_features; j++) {
				labeled_X[i][j] = iris_data.get("X")[i+s][j];
				y[i][0] = iris_data.get("y")[i+s][0];
			}
			c++;
			if(c==n_labeled) {
				s+=n_unlabeled;
				c=0;
			}
		}
		
		c = 0;
		s = n_labeled;
		//Unlabeled data X_2
		for(int i=0; i<n_unlabeled_data; i++) {
			for(int j=0; j<n_features; j++) {				
				unlabeled_X[i][j] = iris_data.get("X")[i+s][j];
			}
			unlabeled_y[i][0] = iris_data.get("y")[i+s][0];
			c++;
			if(c==n_unlabeled) {
				s+=n_labeled;
				c=0;
			}
		}
				
		SSPPCA ssppca = new SSPPCA(labeled_X, y, unlabeled_X, 4, false);
		//HashMap<String, double [][]> init_pars = get_dummy_init_pars();
		//ssppca.set_external_init_pars(init_pars);
		ssppca.do_SSPPCA();
	
		ArrayList<double[][]> trained_data_X1 = new ArrayList<double[][]>();
		ArrayList<double[][]> trained_data_X2 = new ArrayList<double[][]>();
		ArrayList<double[][]> labeled_org_data = new ArrayList<double[][]>();
		ArrayList<double[][]> unlabeled_org_data = new ArrayList<double[][]>();
		
		double [][] X_1_trained = ssppca.get_rotated_X_1_and_X_2().get("X_1");//ssppca.predict(labeled_X).get("X");
		double [][] X_2_trained = ssppca.get_rotated_X_1_and_X_2().get("X_2");	
		double [][] y_trained = ssppca.get_rotated_Y();//ssppca.predict(labeled_X).get("Y");
		double [][] y_unlabeled_predicted = ssppca.predict(unlabeled_X).get("Y");
		
		trained_data_X1.add(y_trained);
		trained_data_X2.add(y_unlabeled_predicted);
		labeled_org_data.add(y);
		unlabeled_org_data.add(unlabeled_y);
		
		int n = X_1_trained[0].length;
		for(int i=0; i<n; i++) {
			trained_data_X1.add(MatrixOperations.get_column_vec_from_matrix(X_1_trained, i));
			trained_data_X2.add(MatrixOperations.get_column_vec_from_matrix(X_2_trained, i));
			labeled_org_data.add(MatrixOperations.get_column_vec_from_matrix(labeled_X,i));
			unlabeled_org_data.add(MatrixOperations.get_column_vec_from_matrix(unlabeled_X,i));
		}
		
		String [] titleNames = {"Species", 
				            	"Sepal Length",
				            	"Sepal Width",
				            	"Petal Length",
				            	"Petal Width"};
		
		n+=1;
		String [] titles = new String [2*n];
		c=0;
		for(int i=0; i<2*n; i++) {
			titles[i] = titleNames[c];
			i++;
			titles[i] = titleNames[c];
			c++;
		}
		
		HashMap<String, Integer> config = new HashMap<String, Integer>();
		config.put("n_cols",2);
		config.put("n_rows",n);
		
		plotSSPPCATrainedData(labeled_org_data, unlabeled_org_data, trained_data_X1, trained_data_X2, titles, config);
			
		System.out.println(ssppca.get_sigma_y());
		System.out.println(ssppca.get_sigma_x());
		
	}
	
	
	@SuppressWarnings("static-access")
	public static HashMap<String, double [][]> get_iris_data_4_classification_task() {
		TestPC obj_pcaTest = new TestPC();
		double [][] irisData = obj_pcaTest.get_iris_data();
		int n = irisData.length;
		double [][] y = new double [n][1];
		
		for(int i=0; i<50; i++) {
			y[i][0] = 0.0;
		}
		for(int i=50; i<100; i++) {
			y[i][0] = 1.0;
		}
		for(int i=100; i<150; i++) {
			y[i][0] = 2.0;
		}
		
		HashMap<String, double [][]> data = new HashMap<String, double [][]>();
		data.put("X", irisData);
		data.put("y", y);
		return data;
	}
	
	
	public static HashMap<String, Object[][]> get_boston_housing_data_4_regression_task() {
		
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/DecisionTrees/BostonHousePriceData.txt";
		InputDataManager inputData = new InputDataManager();	
		try {
			inputData.fileReader(dirName, false, true, true);
		} catch (Exception e) {
			System.out.println("Can´t find data set for upload.");
			e.printStackTrace();
		}
		
		String [][] input = inputData.strFileData;
		
		String [] colnames = inputData.colnames;
		String [][] labels = Utilities.Utilities.convertStringVec(colnames);
		
		int n_rows = input.length;
		int n_cols = input[0].length;
		Double [][] data = new Double [n_rows][n_cols];
		for(int i=0; i<n_cols; i++) {
			for(int j=0; j<n_rows; j++) {
				data[j][i] =  Double.parseDouble(input[j][i]);
			}			
		}
		
		HashMap<String, Object[][]> boston_data = new HashMap<String, Object[][]>();
		boston_data.put("X",data);
		boston_data.put("var_names", labels);
		return boston_data;
	}
	

	@SuppressWarnings("static-access")
	public static void plotSPPCATrainedData(ArrayList<double[][]> org_data, ArrayList<double[][]> trained_data, String [] titles, HashMap<String, Integer> config) {
		
		int n = org_data.size();
		
	 	String yLabels [] = new String [n];
	 	String xLabels [] = new String [n];

		GenGraphics graph = new GenGraphics();		
	 	graph.setNumberOfPlotColums(config.get("n_cols"));
	 	graph.setNumberOfPlotRows(config.get("n_rows")); 	
	 	graph.setGraphWidth(900);
	 	graph.setGraphHeight(600);
	 		 	
	 	double [][] x = graph.get_default_x_axis_labels(org_data.get(0).length);
	 	for(int i=0; i<n; i++) {
		 	yLabels[i] = "";
		 	xLabels[i] = "";
	 		graph.plotLines(x, trained_data.get(i), true, Color.RED);
	 		graph.plotPoints(x, org_data.get(i), false, Color.BLUE);
	 	}
	 		
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "9");
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
	
	
	@SuppressWarnings("static-access")
	public static void plotSSPPCATrainedData(ArrayList<double[][]> labeled_org_data, ArrayList<double[][]> unlabeled_org_data, ArrayList<double[][]> trained_data_X1, ArrayList<double[][]> trained_data_X2, String [] titles, HashMap<String, Integer> config) {
		
		int n = labeled_org_data.size();
		
	 	String yLabels [] = new String [n];
	 	String xLabels [] = new String [n];

		GenGraphics graph = new GenGraphics();		
	 	graph.setNumberOfPlotColums(config.get("n_cols"));
	 	graph.setNumberOfPlotRows(config.get("n_rows")); 	
	 	graph.setGraphWidth(900);
	 	graph.setGraphHeight(600);
	 		 	
	 	double [][] x1 = graph.get_default_x_axis_labels(labeled_org_data.get(0).length);
	 	double [][] x2 = graph.get_default_x_axis_labels(unlabeled_org_data.get(0).length);
	 	for(int i=0; i<n; i++) {
		 	yLabels[i] = "";
		 	xLabels[i] = "";
	 		graph.plotLines(x1, trained_data_X1.get(i), true, Color.RED);
	 		graph.plotPoints(x1, labeled_org_data.get(i), false, Color.BLUE);
	 		graph.plotLines(x2, trained_data_X2.get(i), true, Color.RED);
	 		graph.plotPoints(x2, unlabeled_org_data.get(i), false, Color.BLUE);
	 	}
	 		
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "9");
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
	
	
	public static HashMap<String, double [][]> get_dummy_init_pars() {
		
		double [] w_y = {-1.0549515894240349E-5, 5.291941949226044E-6, -4.084462024939727E-6, 1.0665868571878978E-5};

		double [] w_x = {1.6152233921566986E-5, 9.397460082849387E-6, -1.7293922811048556E-5, -7.055389791567451E-6, 
						-1.5812953841767713E-5, -1.0806995209241715E-5, -1.2753525267814404E-5, -1.7940281501595482E-5,
						4.019157361372715E-6, -1.597281455303131E-5, -7.859564434533353E-6, 1.5542827205446673E-5, 
						-1.7033894050265344E-5, 1.5528886335297018E-5, -7.808588025635451E-6, 9.132192068773571E-6, 
						};

		double [][] W_y = MatrixOperations.transpose(MatrixOperations.convArrayToVec(w_y));
		double [][] W_x = MatrixOperations.transpose(MatrixOperations.convVecToMatrix(w_x, 4, 4));
		
		double [][] sigma_x = new double [1][1];
		double [][] sigma_y = new double [1][1];
		sigma_x[0][0] = 1.1955648915318308E-10; 	
		sigma_y[0][0] = 0.0463850852614263;
		
		double [][] n_factors = new double [1][1];
		n_factors[0][0] = 4;
		
		HashMap<String, double [][]> init_pars = new HashMap<String, double [][]>();
		init_pars.put("n_factors", n_factors);
		init_pars.put("sigma_x", sigma_x);
		init_pars.put("sigma_y", sigma_y);
		init_pars.put("W_x", W_x);
		init_pars.put("W_y", W_y);
		
		return init_pars;
	}
	
	
	public static void main(String[] args) throws Exception {	
		//test1_SPPCA();
		//test2_SPPCA();
		test1_SSPPCA();
	}
	
}
