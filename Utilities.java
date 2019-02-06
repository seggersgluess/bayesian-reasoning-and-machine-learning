
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
            	
            	idx = idx + 1;
            	
            }
			
		}
		
		return idxs;
		
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
	
	
	// draws random int numbers over supplied interval
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
 	
	
    // test client
    public static void main(String[] args) {
    	
    	
    	//double [] a = {-1.0, -1.0, 8.0, -7.0};
    	//double [] sorted_vec = get_sorted_vec(a);
    	//double [] sorted_vec = get_unique_elements(a);
    	//int [] idxs =  get_idxs_for_sorted_vec(a);
    	
    	//double [] x = MatrixOperations.get_double_sub_vec_4_indices(a, idxs);
    	
    	//double [] dbl_idxs = convert_int_to_double_array(idxs);
    	
    	//System.out.print(sorted_vec.length);
    	
    	//String [] b= {"Hallo", "Bla", "Hallo"};
    	//int [] idxs = get_idx(b, "H");
    	
    	//int [] a = getRandomIntNumbers(1,5,10);
    	
       	//MatrixOperations.print_vector(convert_int_to_double_array(idxs));
    	//System.out.print(sorted_vec[1]);
    	
    }
	
}
