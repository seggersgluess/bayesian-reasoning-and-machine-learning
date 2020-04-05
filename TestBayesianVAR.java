package TimeSeriesAnalysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import DataManagement.InputDataManager;
import Mathematics.MatrixOperations;

public class TestBayesianVAR {

	static String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Data/MacroData/Data_txt";
	static String [][] data_matrix;
	static String [][] data_labels;
	static int n_data_series = 0;
	
	static HashMap<String, String [][]> prepared_data;
	
	
	public static void load_macro_data() {
		
		ArrayList<String> seriesNames = getSeriesNamesInFolder();
			
		setDefaultDataMatrix(seriesNames);
		
		n_data_series = seriesNames.size();
		System.out.println("--- Found " + n_data_series + " macroeconomic time series ---");
		Utilities.Utilities.print_string_array(data_labels);
		
		int n_series = seriesNames.size();
		
		for(int i=0; i<n_series; i++) {
			
			String seriesName = seriesNames.get(i);
			String fileName = dirName +"/"+seriesName;
			InputDataManager inputData = new InputDataManager();		
			try {
				inputData.fileReader(fileName, false, false, false);
			} catch (Exception e) {
				e.printStackTrace();
			}
	    		
			String [][] data = inputData.strFileData;
			data = filterUltimoDatesFromSeries(data);
			//Utilities.Utilities.print_string_array(data);
			fill_data_matrix(data, data_labels[i][0]);
		}		
	}
	
	
	public static String [][] upload_information_for_log_trans() {
		String dirName = "C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Data/MacroData/LogRateControl.txt";
		InputDataManager inputData = new InputDataManager();		
		try {
			inputData.fileReader(dirName, false, false, false);
		} catch (Exception e) {
			e.printStackTrace();
		}
    		
		String [][] infos = inputData.strFileData;
		return infos;
	}
	
	
	public static void prepare_loaded_data_matrix_4_analysis() {
		String [][] infos4log = upload_information_for_log_trans();
		if(infos4log.length < n_data_series) {
			throw new RuntimeException("Only " + infos4log.length + " infos for log found for " + n_data_series + " uploaded series.");
		}
		if(data_matrix == null) {
			throw new RuntimeException("No macroeconomic data uploaded yet. Cannot found data matrix for preparation.");
		}
		int n = data_matrix.length;
		int startIdx = 0;
		int endIdx = n-1;
		for(int i=0; i<n_data_series;i++) {
			boolean newSeries = true;
			for(int j=0; j<n; j++) {
				if(data_matrix[j][i+1]!=null) {
					if(newSeries == true) {
						if(j>startIdx) {;
							startIdx = j;
							newSeries = false;
							break;
						}
					}
				}
			}
			for(int j=startIdx; j<n; j++) {
				if(data_matrix[j][i+1] == null) {
					if((j-1)<endIdx) {
						endIdx = j-1;
						break;
					}
				}
			}
		}
		
		n = endIdx-startIdx+1;
		String [][] prep_data_matrix = new String [n][n_data_series];
		String [][] used_dates = new String [n][1];
				
		for(int i=0; i<n_data_series; i++) {
			String seriesName = data_labels[i][0];
			int [] idxs = Utilities.Utilities.get_idx(infos4log,seriesName);
			if(idxs[0] == -1) {
				throw new RuntimeException("No log infos found for data series " + seriesName);
			}
			int log = Integer.parseInt(infos4log[idxs[0]][1]);
			for(int j=0; j<n; j++) {
				
				String x = data_matrix[startIdx+j][i+1];
				if(x.contentEquals(".")) {
					System.out.println("Found missing data for " + data_labels[i][0] + " at data " + data_matrix[startIdx+j][0]);
					double x1 = Double.parseDouble(data_matrix[startIdx+j-1][i+1]);
					double x2 = Double.parseDouble(data_matrix[startIdx+j+1][i+1]);
					double x_interPol = x1+(x2-x1)/2.0;	
					x = Double.toString(x_interPol);
					System.out.println("Set linearly interpolated value at t:" + x + "(t-1: "+ x1 +" & t+1: " + x2 + ")");					
				}
				
				if(log==1) {
					prep_data_matrix[j][i] = Double.toString(Math.log(Double.parseDouble(x))); 
				}else {
					prep_data_matrix[j][i] = x;
				}	
				used_dates[j][0] = data_matrix[startIdx+j][0];
			}	
		}
		
		prepared_data = new HashMap<String, String [][]>();
		prepared_data.put("X", prep_data_matrix);
		prepared_data.put("dates", used_dates);
		prepared_data.put("logInfos", infos4log);
	}
	
	
	public static String [][] filterUltimoDatesFromSeries(String [][] dataSeries) {
		
		int n = dataSeries.length-1;
		String [] years = new String [n];
		String [] month = new String [n];
		String [] days  = new String [n];
 		
		ArrayList<Integer> ultimoIdxs = new ArrayList<Integer>();
		
		for(int i=0; i<n; i++) {
			years[i] = dataSeries[i+1][0].substring(0,4);
			month[i] = dataSeries[i+1][0].substring(5,7);
			days[i]  = dataSeries[i+1][0].substring(8,10);
		}
		
		String [] uniqueYears = Utilities.Utilities.get_unique_elements(years);
		int nYears = uniqueYears.length;
		
		for(int i=0;i<nYears; i++) {
			int [] yearIdx = Utilities.Utilities.get_idx(years, uniqueYears[i]);
			int nSelYearIdxs = yearIdx.length;
			String [] selMonthInYear = new String [nSelYearIdxs];
			for(int j=0; j<nSelYearIdxs; j++) {
				selMonthInYear[j] = month[yearIdx[j]];
			}
			String [] uniqueMonth = Utilities.Utilities.get_unique_elements(selMonthInYear);
			int nMonth = uniqueMonth.length;
			for(int j=0; j<nMonth; j++) {
				int [] monthIdx = Utilities.Utilities.get_idx(selMonthInYear, uniqueMonth[j]);
				int nSelMonthIdxs = monthIdx.length;
				ultimoIdxs.add(yearIdx[monthIdx[nSelMonthIdxs-1]]);
			}
		}
		
		n = ultimoIdxs.size();
		String [][] ultimoData = new String [n][2];
		for(int i=0; i<n; i++) {
			int idx = ultimoIdxs.get(i);
			ultimoData[i][0] = dataSeries[idx+1][0];
			ultimoData[i][1] = dataSeries[idx+1][1];
		}
		
		return ultimoData;
	}
	
	
	public static ArrayList<String> getSeriesNamesInFolder() {
		
		ArrayList<String> fileNames = new ArrayList<String>();
		File folder = new File(dirName);
		File[] listOfFiles = folder.listFiles();
		int nFiles = listOfFiles.length;
		for(int i=0; i<nFiles; i++) {
			if(listOfFiles[i].isFile()) {
				fileNames.add(listOfFiles[i].getName());
			}
		}
				
		return fileNames;
	}
	
	
	public static void setDefaultDataMatrix(ArrayList<String> seriesNames) {
		
		//Jan 1900 - Feb 2020
		int n_defaultMonth = 121*12;
		
		int n_series = seriesNames.size();
		
		data_matrix = new String [n_defaultMonth+1][n_series+1];
		data_matrix[0][0] = "Dates";
		
		int year = 1900;
		int monthCount = 0;
		
		for(int i=0; i<n_defaultMonth; i++) {
			monthCount++;
			int month = monthCount;
			String strMonth = (month < 10 ? "0" : "") + month;
			String strDate = year+"-"+strMonth;
			data_matrix[i+1][0] = strDate;
			if(monthCount==12) {
				monthCount = 0;
				year++;
			}
		}
		
		data_labels = new String [n_series][1];
		for(int i=0; i<n_series; i++) {
			String seriesName = seriesNames.get(i);
			seriesName = seriesName.substring(0, (seriesName.length()-4));
			data_labels[i][0] = seriesName;
			data_matrix[0][i+1] = seriesName;
		}
	}
	
	
	public static void fill_data_matrix(String [][] series, String name) {
		
		if(data_matrix == null) {
			throw new RuntimeException("No (default) data matrix found. Set default data matrix.");
		}
		
		int n_def = data_matrix.length-1;
		int n = series.length;
		int startIdx = 0;
		int colIdx = -1;
		
		int n_variables = data_matrix[0].length;
		
		for(int i=0; i<n_variables; i++) {
			if(data_matrix[0][i].contentEquals(name)==true) {
				colIdx = i;
			}
		}
		if(colIdx == -1) {
			throw new RuntimeException("Series label " + name + " not found in default data matrix.");
		}
		
		for(int i=0; i<n_def; i++) {
			String firstDate = series[0][0].substring(0,7);
			if(data_matrix[i+1][0].contentEquals(firstDate)==true) {
				startIdx = i+1;
				break;
			}
		}
		
		for(int i=0; i<n; i++) {
			String firstDate = series[i][0].substring(0,7);
			if(data_matrix[startIdx+i][0].contentEquals(firstDate)==true) {
				data_matrix[startIdx+i][colIdx] = series[i][1];
			}else {
				data_matrix[startIdx+i][colIdx] = "NA";
				System.out.println(series[i][0] + " invalid Date.");
			}
		}	
	}
	
	
	public static HashMap<String, String [][]> get_prepared_data_4_analysis() {
		if(prepared_data == null) {
			throw new RuntimeException("No macroecnomic data prepared yet.");
		}
		return prepared_data;
	}
	
	
	public static double [][] get_prep_data_matrix_X_for_analysis(){
		if(prepared_data == null) {
			throw new RuntimeException("No macroecnomic data prepared yet.");
		}
		int n_rows = prepared_data.get("X").length;
		int n_cols = prepared_data.get("X")[0].length;
		double [][] X = new double [n_rows][n_cols];
		for(int i=0; i<n_rows; i++) {
			for(int j=0; j<n_cols; j++) {
				X[i][j] = Double.parseDouble(prepared_data.get("X")[i][j]);
			}
		}
		return X;
	}
	

	public static String [][] get_macro_series_labels() {
		if(data_labels == null) {
			throw new RuntimeException("No macroecnomic data uploaded yet.");
		}
		return data_labels;
	}
	
	
	public static void set_dirName2MacroData(String path) {
		dirName = path;
	}
	
	
	public static void testBayesianVAR() {
		set_dirName2MacroData("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Data/MacroData/Data_txt/");
		load_macro_data();
		prepare_loaded_data_matrix_4_analysis();
		//Utilities.Utilities.print_string_array(prepared_data.get("dates"));
		double [][] X = get_prep_data_matrix_X_for_analysis();
		
		BayesianVAR bVAR = new BayesianVAR(X, 3, 0.03);
		bVAR.set_deltas(MatrixOperations.unit_vector(X[0].length));
		bVAR.set_variable_names(data_labels);
		bVAR.estimate_bayesianVAR();
		bVAR.plot_time_series_and_fitted_values();
	}
	
	
	public static void main(String[] args) {
		testBayesianVAR();		
	}

}
