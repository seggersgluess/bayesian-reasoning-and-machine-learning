
import java.io.*; 

public class InputDataManager { 
	
	static String file;
	static String [] colnames;
	static String [] rownames;
	static int numberOfColumns;
	static int numberOfRows;
	static String seperator;
	static String [][] strFileData;
	static double [][] dblFileData;
	
	//if data selection is done
	static String [] selected_colnames;
	static String [] selected_rownames;
	static String [][] selectedStrFileData;
	static double [][] selectedDblFileData;
	
	//reads input data as string array
    public static void fileReader(String fileName, boolean dblData, boolean hasRowNames, boolean hasColNames) throws Exception{
    	
 	    file = fileName;  
	  	
		numberOfLines();
		
		BufferedReader reader = getBufferedReader();
			
		String st;
		String [] strArray = null;
		double [] dblArray = null;
		
  		int rowCorrIdx = 0;
  		int colCorrIdx = 0;
		  		
	  	if((st = reader.readLine()) != null){
	  		
	  		st = Utilities.trim(st);
	  		determine_and_set_seperator(st);
	  		numberOfColumns(st);
	  			
	  	}
	  	
	  	reader.close();
	  	
	  	reader = getBufferedReader();
	  	
	  	if(hasRowNames == true){
	  		
	  		colCorrIdx = 1;
	
	  	}
	  	
	  	if(hasColNames == true){
	  		
	  		rowCorrIdx = 1;
	  		
	  	}
	  	
	  	if(dblData == true){
	  		
	  		dblFileData = new double [numberOfRows-rowCorrIdx][numberOfColumns-colCorrIdx];
	  		
	  	}else{
	  		
	  		strFileData = new String [numberOfRows-rowCorrIdx][numberOfColumns-colCorrIdx];
	  		
	  	}
	  	
	  	colnames = new String [numberOfColumns-colCorrIdx];
	  	rownames = new String [numberOfRows-rowCorrIdx];
	
	  	int lineCounter = 0;
	  	
	  	while((st = reader.readLine()) != null){
	  				
	  		st = Utilities.trim(st);
	  		
	  		if(st.length() == 0){
	  			
	  			break;
	  			
	  		}
	  		
	  		if(hasColNames == true && lineCounter == 0){
	  			
	  			strArray = convertLineIntoStringArray(st);
	  			
	  			if(hasRowNames == true){
	  				
	  				for(int i=0; i<numberOfColumns-1; i++){
	  					
	  					colnames[i] = strArray[i+1];
	  					
	  				}
	  				
	  			}else{
	  				
	  				colnames = strArray;
	  				
	  			}
	  				
	  		}else{
	  			
	  			if(dblData == true){
	  				
	  				dblArray = convertLineIntoDoubleArray(st);
	  				
	  			}else{
	  				
	  				strArray = convertLineIntoStringArray(st);
	  				
	  			}
	  				
	  			if(hasRowNames == true){
  					
  					rownames[lineCounter-rowCorrIdx] = strArray[0];
  					
  				}
	  			
	  			for(int i=0; i<numberOfColumns-colCorrIdx; i++){
		  			
	  				if(dblData == true){
	  					
	  					if(i==0){
	  						
	  						int n_missing_values = dblArray.length - colCorrIdx - (numberOfColumns-colCorrIdx);
	  						
	  						if(n_missing_values != 0){	  							
	  							throw new RuntimeException(n_missing_values + " in line " + (lineCounter-rowCorrIdx) + ". Check your input data.");		  							
	  						}
	  						
	  					}
	  					
	  					dblFileData[lineCounter-rowCorrIdx][i] = dblArray[i+colCorrIdx];
	  					
	  				}else{
	  					
	  					if(i==0){
	  						
	  						int n_missing_values = strArray.length - colCorrIdx - (numberOfColumns-colCorrIdx);
	  						
	  						if(n_missing_values != 0){	  							
	  							throw new RuntimeException(n_missing_values + " in line " + (lineCounter-rowCorrIdx) + ". Check your input data.");		  							
	  						}
	  						
	  					}
	  					
	  					strFileData[lineCounter-rowCorrIdx][i] = strArray[i+colCorrIdx];
	  					
	  				}
	  					
		  		}
	  					
	  		}

	  		lineCounter++;
	  		
	  	} 
	  	
	  	reader.close();
	  	
	}

    
    public static void show_loaded_data(){
    	
    	if(strFileData != null){
    		
    		Utilities.print_string_array(strFileData);
    		
    	}else{
    		
    		MatrixOperations.print_matrix(dblFileData);
    		
    	}
    		
    }
    
    
    public static void show_selected_loaded_data(){
    	
    	if(selectedStrFileData != null){
    		
    		Utilities.print_string_array(selectedStrFileData);
    		
    	}else{
    		
    		MatrixOperations.print_matrix(selectedDblFileData);
    		
    	}
    		
    }
    
    
    public static void show_colnames(){
    	
    	if(colnames != null){
    		
    		Utilities.print_string_array(colnames);
    		
    	}else{    		
    		throw new RuntimeException("No column names supplied in loaded data.");    		
    	}
    		
    }
    
    
    public static BufferedReader getBufferedReader(){
    	
		BufferedReader reader = null;
		
		try {
			reader = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
    	
		return reader;
		
    }
    
    
    public static String [] get_allowed_seperators(){
    	
    	String [] allowed_seperators = {",", ";", " "};
    	
    	return allowed_seperators;
    	
    }
    
    
    public static void determine_and_set_seperator(String strLine){
    	
    	String [] allowed_seps = get_allowed_seperators();
    	
        int n_allowed_seps = allowed_seps.length;
    	int posIdx;
        int match = 0;
    	
        for(int i=0; i<n_allowed_seps; i++){
        	
        	posIdx = strLine.indexOf(allowed_seps[i], 0);
        	
        	if(posIdx != -1){
        		
        		seperator = allowed_seps[i];
        		match = 1;
        		
        		break;
        		
        	}
        	
        }
        
    	if(match == 0){    		
    		throw new RuntimeException("Can´t determine seperator of columns.");    		
    	}
    	
    }
    
    
    public static void numberOfLines() throws Exception{
    	
    	String st;
    	
		BufferedReader reader = getBufferedReader();

    	int numberOfLines = 0;

	  	while((st = reader.readLine()) != null){
		  	
	  		if(Utilities.trim(st).length() == 0){
	  			
	  			break;
	  			
	  		}
	  			
	  		numberOfLines++;
	  						
	  	} 
    	
	  	reader.close();
	  	
	  	numberOfRows = numberOfLines;
	  	
    }
    
    
    public static void numberOfColumns(String strLine){
    	
    	int n_cols = 0;
    	int posIdx1 = 0;
    	int posIdx2 = strLine.indexOf(seperator, 0);;
    	
    	while(posIdx2 != -1){
    		
    		n_cols  = n_cols  + 1;
		
    		posIdx1 = posIdx2+1;
    		posIdx2 = strLine.indexOf(seperator, posIdx1);
    			
    	}
    	
    	if(posIdx1 != strLine.length()){
    		
    		n_cols = n_cols+1;
    		
    	}
    	
    	numberOfColumns = n_cols ;
    	
    }
    
    
    //returns double values from string line
    public static double [] convertLineIntoDoubleArray(String strLine){
    	    	
    	int sepCount = 0;
    	  		
    	int posIdx1;
    	int posIdx2;
    	
    	double dbl = 0;
    	double [] dblArray = new double [numberOfColumns];
    	
    	posIdx1 = 0;
    	posIdx2 = strLine.indexOf(seperator, 0);
    
    	while(posIdx2 != -1){
    		
    		try{
    			dbl = Double.valueOf(strLine.substring(posIdx1, posIdx2));
    		}catch(Exception e){			
    			throw new RuntimeException("'" + strLine.substring(posIdx1, posIdx2) + "' is not a number. Please check your input data.");			
    		}
    		
    		if(sepCount == 0){
    			
    			dblArray[0] = dbl;
    			
    		}else{
    			
    			dblArray[sepCount] =  dbl;
    			
    		}
    		
    		posIdx1 = posIdx2+1;
    		posIdx2 = strLine.indexOf(seperator, posIdx1);
    		
    		sepCount     = sepCount + 1;
    		
    	}
    	
    	if(sepCount != 0){
    		
    		try{
    			dbl = Double.valueOf(strLine.substring(posIdx1, strLine.length()));
    		}catch(Exception e){			
    			throw new RuntimeException("'" + strLine.substring(posIdx1, posIdx2) + "' is not a number. Please check your input data.");			
    		}
    		
    		dblArray[sepCount]  = dbl;
    		
    	}
    	
    	if(sepCount == 0){
    		
    		dblArray[0] = Double.valueOf(strLine);
    		
    	}
    	
    	return dblArray;
    	
    }
    
    
    //returns double values from string line
    public static String [] convertLineIntoStringArray(String strLine){
    	    	
    	int sepCount = 0;
    	  		
    	int posIdx1;
    	int posIdx2;
    	
    	String str;
    	String [] strArray = new String [numberOfColumns];
    	
    	posIdx1 = 0;
    	posIdx2 = strLine.indexOf(seperator, 0);
    
    	while(posIdx2 != -1){
    		  
    		str = strLine.substring(posIdx1, posIdx2);
    		
    		if(sepCount == 0){
    			
    			strArray[0] = str;
    			
    		}else{
    			
    			strArray[sepCount] =  str;
    			
    		}
    		
    		posIdx1 = posIdx2+1;
    		posIdx2 = strLine.indexOf(seperator, posIdx1);
    		
    		sepCount     = sepCount + 1;
    		
    	}
    	
    	if(sepCount != 0){
    		
    		str = strLine.substring(posIdx1, strLine.length());
    		strArray[sepCount]  = str;
    		
    	}
    	
    	if(sepCount == 0){
    		
    		strArray[0] = strLine;
    		
    	}
    	
    	return strArray;
    	
    }
    
    
    public static void selectLoadedData(String [] rowNames, String [] colNames){
    	
    	int n_rownames = rownames.length;
    	int n_colnames = colnames.length;
    	
    	int[] row_idxs = Utilities.intGenerator(n_rownames);
    	int[] col_idxs = Utilities.intGenerator(n_colnames);
    	
    	if(rowNames != null){
    		
    		if(rownames == null){
    			throw new RuntimeException("No row names supplied in loaded data. Can´t do any row selection.");    			
    		}
    		
    		row_idxs = new int [rowNames.length];
    		
    		for(int i=0; i<rowNames.length; i++){
    			
    			row_idxs[i] = Utilities.get_idx(rownames, rowNames[i])[0];
    			
    			if(row_idxs[i] == -1){
    				throw new RuntimeException(rowNames[i] + " not found as row name of loaded data.");
    			}
    			
    		}

    		n_rownames = rowNames.length;
    		
    		selected_rownames = rowNames;
    		
    	}	
    		
    	if(colNames != null){
    		
    		if(colnames == null){
    			throw new RuntimeException("No column names supplied in loaded data. Can´t do any column selection.");    			
    		}
    		
    		col_idxs = new int [colNames.length];
    		
    		for(int i=0; i<colNames.length; i++){
    			
    			col_idxs[i] = Utilities.get_idx(colnames, colNames[i])[0];
    			
    			if(col_idxs[i] == -1){
    				throw new RuntimeException(colNames[i] + " not found as column name of loaded data.");
    			}
    			
    		}
    		
    		n_colnames = colNames.length;
    		
    		selected_colnames = colNames;
    		
    	}
    	
    	if(dblFileData != null){
    		
    		selectedDblFileData = new double [n_rownames][n_colnames];
    		
    		for(int i=0; i<n_rownames; i++){
    			
    			for(int j=0; j<n_colnames; j++){
    				
					selectedDblFileData[i][j] = dblFileData[row_idxs[i]][col_idxs[j]];
    				
    			}
    			
    		}
    		
    	}else{
    		
    		selectedStrFileData = new String [n_rownames][n_colnames];
    		
    		for(int i=0; i<n_rownames; i++){
    			
    			for(int j=0; j<n_colnames; j++){
    				
					selectedStrFileData[i][j] = strFileData[row_idxs[i]][col_idxs[j]];
    				
    			}
    			
    		}
    		
    	}
    	
    }
    
    
    //test client
    public static void main(String[] args) throws Exception {
    	
	    fileReader("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Tests/LogisticRegression/Data.txt", false, true, true);
	    
	    show_loaded_data();
	    
	    String [] col = {"id", "prog"};
	    
	    selectLoadedData(null, col);
	    
	    show_selected_loaded_data();
	    
	    //Utilities.print_string_array(colnames);
	    //Utilities.print_string_array(rownames);
	    
    }
    
}
	
