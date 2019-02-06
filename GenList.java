
import java.util.ArrayList;
import java.util.function.BiFunction;

public class GenList {

	int length;
	int filledSlots = 0;
	
	@SuppressWarnings("rawtypes")
	ArrayList[] List;
	
	//Constructor
	public GenList(int length){
		
		this.length = length;
		List = new ArrayList[length];
		
	}
	
	
	public void add(double [] arg){
		
		ArrayList<Double> new_element;
		
		if(arg != null){
			
			int n_args = arg.length;
			
			new_element = new ArrayList<Double>(n_args);
			
			for(int i=0; i<n_args; i++){
				
				new_element.add(arg[i]);
				
			}
			
		}else{
			
			new_element = null;
			
		}
		
		List[filledSlots] = new_element;
		filledSlots = filledSlots + 1;
		
	}
	
	
	public double [] get_double_array(int slot_idx){
		
		if(slot_idx > filledSlots-1){
			
			throw new RuntimeException("List has only " + filledSlots + " elements.");
			
		}
		
		double [] array; 
		
		if(List[slot_idx] != null){
			
			int n_args = List[slot_idx].size();
			
			array = new double [n_args];
			
			for(int i = 0; i < n_args; i++){
				
				array[i] = (double) List[slot_idx].get(i);
				
			}
			
		}else{
			
			array = null;
			
		}
		
		return array;
		
	}
	
	
	public void add(int [] arg){
		
		int n_args = arg.length;
		
		ArrayList<Integer> new_element = new ArrayList<Integer>(n_args);
		
		for(int i=0; i<n_args; i++){
			
			new_element.add(arg[i]);
			
		}
		
		List[filledSlots] = new_element;
		filledSlots = filledSlots + 1;
		
	}
	
	
	public void add(String [] arg){
		
		int n_args = arg.length;
		
		ArrayList<String> new_element = new ArrayList<String>(n_args);
		
		for(int i=0; i<n_args; i++){
			
			new_element.add(arg[i]);
			
		}
		
		List[filledSlots] = new_element;
		filledSlots = filledSlots + 1;
		
	}
		
	
	public int length(){
		
		return length;
		
	}
	
	
	public void main(String[] args){
		
    	//GenList genericList = new GenList(2);
    	//double [] arg1 = {0.0, 1.0};
    	//String [] arg2 = {"Hallo", "Blabla", "120349"};
        //GenList.add(arg1);
        //GenList.add(arg2);
    	
        //double [] a;
        //a = genericList.List[0];
        
        //System.out.println(genericList.List[0].get(1));
        //System.out.println(genericList.get_double_array(0));
        //System.out.println(genericList.List[1].size());
        //System.out.println(genericList.filledSlots);
		
        ArrayList <BiFunction<double [], double [], Double>> constraints = new ArrayList <BiFunction<double [], double [], Double>>(1);
        
        constraints.add(TargetFunction::target_function);
        
        double [] arg_1 = {30.05, -2000.1};
        double [] arg_2 = null;
        
        System.out.println(constraints.get(0).apply(arg_1, arg_2));
        
	}
	
}
