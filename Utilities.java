package Utilities;
import java.util.ArrayList;
import java.util.List;

import Mathematics.MatrixOperations;

public class Utilities {

	// returns sorted vector
	public static double [] get_sorted_vec(double [] x){
		
		int n_elements = x.length;
		int idx = 0;
		
		double [] sorted_vec = new double [n_elements];		
		double [] unique_vec = get_unique_elements(x);
		
		for(int i=0; i<unique_vec.length; i++){
			
			double min      = getMin(x);
			int [] min_idxs = get_idx(x, min);
			
			for(int j=0; j<min_idxs.length; j++){
				
				sorted_vec[idx] = min;
				idx = idx + 1;
			}
						
			x = MatrixOperations.get_double_comp_sub_vec_4_indices(x, min_idxs);
			
		}
		
		return sorted_vec;		
		
	}
	
	
	// returns indices of sorted elements in unsorted vector
	public static int [] get_idxs_for_sorted_vec(double [] x){
		
		double [] sorted_vec = get_sorted_vec(x);
		
		int [] idxs = new int [sorted_vec.length];
		
		double [] unique_vec = get_unique_elements(sorted_vec);
		
		int idx = 0;
		
		for(int i=0; i<unique_vec.length; i++){
			
			int [] search_idx = get_idx(x, unique_vec[i]);
			
            for(int j=0; j<search_idx.length; j++){            	
            	idxs[idx] = search_idx[j];            	
            	idx++;            	
            }
			
		}
		
		return idxs;
		
	}
	
	
	// returns unique elements of a supplied integer list
	public static List<Integer> get_unique_elements(List<Integer> x){
		
		int n = x.size();
		
		List<Integer> unique_elements = new ArrayList<Integer>();
		unique_elements.add(x.get(0));
		
		for(int i=1; i<n; i++){
			
			int [] idx = get_idx(unique_elements,x.get(i));
			
			if(idx[0] == -1){
				unique_elements.add(x.get(i));
			}
			
		}
		
		return unique_elements;
		
	}
	
	
	// returns unique elements of a supplied vector
	public static double [] get_unique_elements(double [] x){
		
		int n_elements = x.length;
		int idx = 0;
		
		double [] unique_vec = new double [n_elements];		
				
		for(int i=0; i<n_elements; i++){
			
			int [] idxs = get_idx(x, x[0]);
				
			unique_vec[idx] = x[0];
			idx = idx + 1;
						
			if(x.length == 1){
				
				unique_vec = MatrixOperations.get_double_sub_vec(unique_vec, 0, idx-1);
				break;
				
			}
			
			x = MatrixOperations.get_double_comp_sub_vec_4_indices(x, idxs);
			
			if(x.length == 0){
				
				unique_vec = MatrixOperations.get_double_sub_vec(unique_vec, 0, idx-1);
				break;
				
			}
			
		}
		
		return unique_vec;		
		
	}
	
	
	// returns unique elements of a supplied array
	public static String [] get_unique_elements(String [] x){
		
		int n_elements = x.length;
		int idx = 0;
		
		String [] unique_vec = new String [n_elements];		
				
		unique_vec[0] = x[0];
		
		for(int i=1; i<n_elements; i++){
			
			int match = 0;
			
			for(int j=0; j<idx; j++){
				
				if(unique_vec[j].contentEquals(x[i]) == true){
					
					match = 1;
					break;
					
				}
				
			}
			
			if(match == 0){
				
				unique_vec[idx] = x[i];
				idx++;
					
			}
			
		}
	
		String [] unique_vec_2 = new String[idx];
		
		for(int i=0; i<idx; i++){
			
			unique_vec_2[i] = unique_vec[i];
			
		}
		
		return unique_vec_2;		
		
	}
	
	
	// returns unique elements of a supplied n x 1 column vector
	public static String [] get_unique_elements(String [][] x){
		
		int n_elements = x.length;
		int idx = 0;
		
		String [] unique_vec = new String [n_elements];		
				
		unique_vec[0] = x[0][0];
		
		for(int i=1; i<n_elements; i++){
			
			int match = 0;
			
			for(int j=0; j<idx; j++){
				
				if(unique_vec[j].contentEquals(x[i][0]) == true){
					
					match = 1;
					break;
					
				}
				
			}
			
			if(match == 0){
				
				unique_vec[idx] = x[i][0];
				idx++;
					
			}
			
		}
	
		String [] unique_vec_2 = new String[idx];
		
		for(int i=0; i<idx; i++){
			
			unique_vec_2[i] = unique_vec[i];
			
		}
		
		return unique_vec_2;		
		
	}
	
	
	// returns maximum of a supplied vector
	public static double getMax(double[] x){ 
		    
		double maxValue = x[0]; 
		int n_elements = x.length;   
		
		for(int i=1;i < n_elements;i++){ 
		    	
		  if(x[i] > maxValue){ 
		    	  
		     maxValue = x[i]; 
		         
		  } 
		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of a supplied vector
	public static int getMax(int[] x){ 
		    
		int maxValue = x[0]; 
		int n_elements = x.length;   
		
		for(int i=1;i < n_elements;i++){ 
		    	
		  if(x[i] > maxValue){ 
		    	  
		     maxValue = x[i]; 
		         
		  } 
		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of a supplied vector
	public static double getMax(double[][] x){ 
		    
		double maxValue = x[0][0]; 
		int nRows = x.length; 
		int nCols = x[0].length;
		
		for(int i=0; i<nRows;i++){ 
		    
			for(int j=0; j<nCols; j++){
				
				if(x[i][j] > maxValue){ 
			    	  
					maxValue = x[i][j]; 
				         
				}
				
			}
			  		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of a column in supplied matrix X
	public static double getMax(double[][] X, int colIdx){ 
		    
		if(colIdx > X[0].length){
			throw new RuntimeException("Invalid column index supplied.");
		}
		
		double maxValue = X[0][colIdx]; 
		int n_elements = X.length;   
		
		for(int i=1 ;i < n_elements;i++){ 
		    	
		  if(X[i][colIdx] > maxValue){ 
		    	  
		     maxValue = X[i][colIdx]; 
		         
		  } 
		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of a column in supplied list x
	public static double getMax(List<Double> x){ 
		    
		int n_elements = x.size();
		
		double maxValue = x.get(0); 
		   
		for(int i=1 ;i < n_elements;i++){ 
		  
		  double curVal = x.get(i);
		  if(curVal > maxValue){ 
		    	  
		     maxValue = curVal; 
		         
		  } 
		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of a column in supplied list x
	public static int getMaxFromIntList(List<Integer> x){ 
		    
		int n_elements = x.size();
		
		int maxValue = x.get(0); 
		   
		for(int i=1 ;i < n_elements;i++){ 
		  
		  int curVal = x.get(i);
		  if(curVal > maxValue){ 
		    	  
		     maxValue = curVal; 
		         
		  } 
		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of integer array of arrays
	public static int getMaxFromIntList(ArrayList<List<Integer>> x){ 
		    
		int n_elements = x.size();		
		int maxValue   = x.get(0).get(0); 
		   
		for(int i=0; i<n_elements; i++){ 
		  
			for(int j=0; j<x.get(i).size(); j++){
				
				int curVal = x.get(i).get(j);
				
				  if(curVal > maxValue){ 				    	  
				     maxValue = curVal; 			         
				  }
				
			}
			   
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of a column in supplied list x
	public static double getMaxFromDblList(List<Double> x){ 
		    
		int n_elements = x.size();
		
		double maxValue = x.get(0); 
		   
		for(int i=1; i<n_elements; i++){ 
		  
		  double curVal = x.get(i);
		  if(curVal > maxValue){ 
		    	  
		     maxValue = curVal; 
		         
		  } 
		      
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns maximum of integer array of arrays
	public static double getMaxFromDblList(ArrayList<List<Double>> x){ 
		    
		int n_elements = x.size();		
		double maxValue   = x.get(0).get(0); 
		   
		for(int i=0; i<n_elements; i++){ 
		  
			for(int j=0; j<x.get(i).size(); j++){
				
				double curVal = x.get(i).get(j);
				
				  if(curVal > maxValue){ 				    	  
				     maxValue = curVal; 			         
				  }
				
			}
			   
		} 
		    
		return maxValue; 
		    
	}
	
	
	// returns minimum of a supplied vector
	public static double getMin(double[] x){ 
	    
		double minValue = x[0]; 
		int n_elements = x.length;   
		
		for(int i=1;i < n_elements;i++){ 
		    	
		  if(x[i] < minValue){ 
		    	  
		     minValue = x[i]; 
		         
		  } 
		      
		} 
		    
		return minValue; 
		    
	}
	
	
	// returns minimum of a supplied vector
	public static double getMin(int[] x){ 
	    
		double minValue = x[0]; 
		int n_elements = x.length;   
		
		for(int i=1;i < n_elements;i++){ 
		    	
		  if(x[i] < minValue){ 
		    	  
		     minValue = x[i]; 
		         
		  } 
		      
		} 
		    
		return minValue; 
		    
	}
	
	
	// returns minimum of a supplied vector/ matrix
	public static double getMin(double[][] X){ 
		    
		double minValue = X[0][0]; 
		int nRows = X.length; 
		int nCols = X[0].length;
		
		for(int i=0; i<nRows;i++){ 
		    
			for(int j=0; j<nCols; j++){
				
				if(X[i][j] < minValue){ 
			    	  
					minValue = X[i][j]; 
				         
				}
				
			}
			  		      
		} 
		    
		return minValue; 
		    
	}
	
	
	// returns minimum of a column in supplied matrix X
	public static double getMin(double[][] X, int colIdx){ 
		    
		if(colIdx > X[0].length){
			throw new RuntimeException("Invalid column index supplied.");
		}
		
		double minValue = X[0][colIdx]; 
		int n_elements = X.length;   
		
		for(int i=1 ;i < n_elements;i++){ 
		    	
		  if(X[i][colIdx] < minValue){ 
		    	  
		     minValue = X[i][colIdx]; 
		         
		  } 
		      
		} 
		    
		return minValue; 
		    
	}
	
	
	// returns minimum of a column in supplied list x
	public static double getMin(List<Double> x){ 
		    
		int n_elements = x.size();
		
		double minValue = x.get(0); 
		   
		for(int i=1 ;i < n_elements;i++){ 
		  
		  double curVal = x.get(i);
		  if(curVal < minValue){ 
		    	  
		     minValue = curVal; 
		         
		  } 
		      
		} 
		    
		return minValue; 
		    
	}
	
	
	// returns maximum of integer array of arrays
	public static double getMinFromDblList(ArrayList<List<Double>> x){ 
		    
		int n_elements  = x.size();		
		double minValue = x.get(0).get(0); 
		   
		for(int i=0; i<n_elements; i++){ 
		  
			for(int j=0; j<x.get(i).size(); j++){
				
				double curVal = x.get(i).get(j);
				
				  if(curVal < minValue){ 				    	  
				     minValue = curVal; 			         
				  }
				
			}
			   
		} 
		    
		return minValue; 
		    
	}
	
	
	// returns position index of an element in supplied vector
	public static int [] get_idx(double [] x, double search_element){
		
		int [] search_idxs = new int [x.length];
		
		int idx = 0;		
				
		for(int i=0; i<x.length; i++){
			
			if(x[i] == search_element){
				
				search_idxs[idx] = i;
				idx = idx+1;
				
			}
			
		}
		
		if(idx == 0){
			
			int [] idxs = new int [1];
			idxs[0] = -1;
			return idxs;
			
		}else{
			
			int [] idxs = MatrixOperations.get_int_sub_vec(search_idxs, 0, (idx-1));
			return idxs;
			
		}
	   	
	}
	
	
	// returns position index of an element in supplied vector
	public static int [] get_idx(String [] x, String search_element){
		
		int [] search_idxs = new int [x.length];
		
		int idx = 0;		
				
		for(int i=0; i<x.length; i++){						
			if(x[i].contentEquals(search_element) == true){				
				search_idxs[idx] = i;
				idx = idx+1;				
			}			
		}
		
		if(idx == 0){			
			int [] idxs = new int [1];
			idxs[0] = -1;
			return idxs;			
		}else{			
			int [] idxs = MatrixOperations.get_int_sub_vec(search_idxs, 0, (idx-1));
			return idxs;			
		}
	   	
	}
	
	
	// returns position index of an element in supplied vector
	public static int [] get_idx(List<String> x, String search_element){
		
		int [] search_idxs = new int [x.size()];
		
		int idx = 0;		
				
		for(int i=0; i<x.size(); i++){
						
			if(x.get(i).contentEquals(search_element) == true){
				
				search_idxs[idx] = i;
				idx = idx+1;
				
			}
			
		}
		
		if(idx == 0){
			
			int [] idxs = new int [1];
			idxs[0] = -1;
			return idxs;
			
		}else{
			
			int [] idxs = MatrixOperations.get_int_sub_vec(search_idxs, 0, (idx-1));
			return idxs;
			
		}
	   	
	}
	
	
	// returns position index of an element in supplied vector
	public static int [] get_idx(List<Integer> x, int search_element){
		
		int [] search_idxs = new int [x.size()];
		
		int idx = 0;		
				
		for(int i=0; i<x.size(); i++){
						
			if(x.get(i) == search_element){
				
				search_idxs[idx] = i;
				idx = idx+1;
				
			}
			
		}
		
		if(idx == 0){
			
			int [] idxs = new int [1];
			idxs[0] = -1;
			return idxs;
			
		}else{
			
			int [] idxs = MatrixOperations.get_int_sub_vec(search_idxs, 0, (idx-1));
			return idxs;
			
		}
	   	
	}
	
	
	// returns position index of an element in supplied vector
	public static int [] get_idx(int [] x, int search_element){
		
		int [] search_idxs = new int [x.length];
		
		int idx = 0;		
				
		for(int i=0; i<x.length; i++){
			
			if(x[i] == search_element){
				
				search_idxs[idx] = i;
				idx = idx+1;
				
			}
			
		}
		
		if(idx == 0){
			
			int [] idxs = new int [1];
			idxs[0] = -1;
			return idxs;
			
		}else{
			
			int [] idxs = MatrixOperations.get_int_sub_vec(search_idxs, 0, (idx-1));
			return idxs;
			
		}
	   	
	}
	
	
	// converts double array into integer array
	public static int [] convert_double_to_int_array(double [] x){
		
		int [] int_array = new int [x.length];
		
		for(int i=0; i<x.length; i++){
			
			int_array[i] = (int) x[i];
			
		}
				
		return int_array;		
		
	}
	
	
	// converts integer array into double array
	public static double [] convert_int_to_double_array(int [] x){
		
		double [] int_array = new double [x.length];
		
		for(int i=0; i<x.length; i++){
			
			int_array[i] = (double) x[i];
			
		}
				
		return int_array;		
		
	}
	
	
	//returns array x in reverse order
	public static double [] reverse_array(double [] x){
		
		int n = x.length;
		double [] rev_x = new double [n];
		
		for(int i=0; i<n; i++){
			
			rev_x[i] = x[n-i-1];
			
		}
		
		return rev_x;
		
	}
	
	
	//returns n integer number between 0 and n-1
	public static int [] intGenerator(int n){
		
		int[] a = new int[n];
		
	    for (int i = 0; i < n; ++i) {
	    	
	        a[i] = i;
	        
	    }
	    
	    return a;
		
	}
	
	
	//draws random int numbers over supplied interval
	public static int [] getRandomIntNumbers(int lb, int ub, int n){
		
		int n_ints = ub - lb + 1;
		
		double int_length = 1.0/n_ints;
		
		int [] random_int_numbers = new int [n];
		
		for(int i = 0; i < n; i++){
			
			double rand_number = Math.random();
			int idx = 0;
			
			for(int j = 0; j < n_ints; j++){
				
				idx = idx + 1;
				
				if(rand_number <= idx*int_length){
					
					random_int_numbers[i] = lb + idx - 1;
					break;
				}
				
			}
					
		}
		
		return random_int_numbers;
		
	}
 	
	
	public static List<Integer> get_idx_of_char_in_str(String str, String searchSymbol){
        		
		List<Integer> pos = new ArrayList<Integer>();
		
		if(str.length() == 1){
			if(str == searchSymbol){
				pos.add(0);
			}
		}else{
			for (int i=1; i<str.length()+1; i++){
		    	
		    	System.out.println(str.substring(i-1,i));
		    	
		    	if (str.substring(i-1,i).contentEquals(searchSymbol)){
		    		pos.add(i-1);
		        }
		    	
		    }
		}
		
	    
	    	
	    return pos;
	    
	}
	
	
	public static String trim(String str){
		    	  		
    	String trim_str = "";
    	
       	int posIdx1;
    	int posIdx2;
    	int sepCount = 0;
    	
    	posIdx1 = 0;
    	posIdx2 = str.indexOf(" ", 0);
    
    	while(posIdx2 != -1){
    		
    		if(posIdx1 != posIdx2 && str.substring(posIdx1, posIdx2) != " "){
    			
    			if(sepCount == 0){
    				
    				trim_str = str.substring(posIdx1, posIdx2);
    				
    			}else{
    				
    				trim_str = trim_str + " " + str.substring(posIdx1, posIdx2);
    			 				
    			}
    			
    			sepCount++;
    			
    		}
    		    		
    		posIdx1 = posIdx2+1;
    		posIdx2 = str.indexOf(" ", posIdx1);
    			
    	}
		
    	if(sepCount != 0){
    		
    		if(posIdx1 != str.length() && str.substring(posIdx1, str.length()) != " "){
    			
    			trim_str = trim_str + " " + str.substring(posIdx1, str.length());
    			
    		}
		
    	}
    	
    	if(sepCount == 0){
    		
    		trim_str = str;
    		
    	}
    	
    	trim_str = trim_str.trim();
    	
    	return trim_str;
    	
	}
	
	
	// prints supplied string array A
	public static void print_string_array(String [][] A){
		
		for (int i = 0; i < A.length; i++) {
			
		    for (int j = 0; j < A[0].length; j++) {
		    	
		        System.out.print(A[i][j] + " ");
		        
		    }
		    
		    System.out.println();
		}
		
	}
	
	
	// prints supplied string array A
	public static void print_string_array(String [] A){
		
		for (int i = 0; i < A.length; i++) {
			  	
			System.out.print(A[i] + " "); 
		    
		}
		
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	//double [] a = {-2.0, -1.0, 8.0, -7.0};
    	//double [] sorted_vec = get_sorted_vec(a);
    	//double [] sorted_vec = get_unique_elements(a);
    	//int [] idxs =  get_idxs_for_sorted_vec(a);
    	
    	//double [] x = MatrixOperations.get_double_sub_vec_4_indices(a, idxs);
    	
    	//double [] dbl_idxs = convert_int_to_double_array(idxs);
    	
    	//System.out.print(sorted_vec.length);
    	
    	//String [] b= {"Hallo", "Bla", "Hallo"};
    	//int [] idxs = get_idx(b, "H");
    	
    	//int [] a = getRandomIntNumbers(1,5,10);
    	
       	//MatrixOperations.print_vector(reverse_array(a));
    	//System.out.print(sorted_vec[1]);
    	
    	//String str = " No     id female    ses  schtyp     prog read write math science socst       honors awards cid ";
    	
    	//String [] a = {"a", "a", "b", "a", "c", "d", "c"};
    	
    	//String [] b = get_unique_elements(a);
    	  	
    	@SuppressWarnings("unused")
		List<Integer> a = get_idx_of_char_in_str("2,3,", ",");
    	
    	
    }
	
}
