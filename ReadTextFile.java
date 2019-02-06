
import java.io.*; 

public class ReadTextFile { 
	
	//returns double array from loaded text file
    public static double [][] readfile(String fileName) throws Exception{
    	
		BufferedReader br = null;
		
		try {
			br = new BufferedReader(new FileReader(fileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} 	  		    
	  	
	  	String st; 
	    int lineCount = 0;
	  	
	  	double [][] fileData = new double [1][1];
	  	
	  	while((st = br.readLine()) != null){
	  			  		
	  		if(lineCount == 0){
	  			
	  			fileData = MatrixOperations.transpose(convertLineIntoDoubleArray(st));
	  			
	  		}else{
	  			
	  			fileData = MatrixOperations.cbind(MatrixOperations.transpose(convertLineIntoDoubleArray(st)),fileData);
	  			
	  		}
	  		
	  		lineCount = lineCount + 1;
	  		
	  	} 
	  	
	  	fileData = MatrixOperations.transpose(fileData);
	  	
	  	MatrixOperations.print_matrix(fileData);
	    
		return fileData;
    
	}

    
    //returns double values from string line
    public static double [][] convertLineIntoDoubleArray(String strLine){
    	
    	int sepCount = 0;
    	  		
    	int posIdx1;
    	int posIdx2;
    	double [][] dbl      = new double [1][1];
    	double [][] dblArray = new double [1][1];
    	
    	posIdx1 = 0;
    	posIdx2 = strLine.indexOf(",", 0);
    	
    	while(posIdx2 != -1){
    		  
    		dbl[0][0] = Double.valueOf(strLine.substring(posIdx1, posIdx2));
    		
    		if(sepCount == 0){
    			
    			dblArray[0][0] = dbl[0][0];
    			
    		}else{
    			
    			dblArray = MatrixOperations.cbind(dblArray, dbl);
    			
    		}
    		
    		posIdx1 = posIdx2+1;
    		posIdx2 = strLine.indexOf(",", posIdx1);
    		
    		sepCount     = sepCount + 1;
    		
    	}
    	
    	if(sepCount != 0){
    		
    		dbl[0][0] = Double.valueOf(strLine.substring(posIdx1, strLine.length()-1));
    		dblArray  = MatrixOperations.cbind(dblArray, dbl);
    		dblArray  = MatrixOperations.reverse(dblArray);
    		
    	}
    	
    	if(sepCount == 0){
    		
    		dblArray[0][0] = Double.valueOf(strLine);
    		
    	}
    	
    	return dblArray;
    	
    }
    
    
    //test client
    public static void main(String[] args) throws Exception {
    	
	    double [][] a = readfile("C:/Users/sven_/Documents/Bayesian_Reasoning_and_ML/Test.txt");

	    System.out.println(a.length);
	    
    	
    }
    
}
	


  

  

