
public class Cholesky {

	public static double [][] decompose(double [][] A){
		
		int n = A.length;
		
		if(n != A[0].length){
			
			throw new RuntimeException("Suplied matrix is not symmetric.");
			
		}
		
		double [][] L = new double [n][n];
		
		int nCols = 0;
		double sum = 0.0;
		
		for(int i=0; i<n; i++){
		
			nCols = nCols+1;
		    	
			for(int j=0; j<nCols; j++){
				
				sum   = 0.0;
				
				if(i == j){
					
					for(int k=0; k<j; k++){
						
						sum = sum + Math.pow(L[i][k],2.0);
						
					}
					
					L[i][j] = Math.sqrt(A[i][j] - sum);
					
				}
				
				if(i > j){
					
					for(int k=0; k<j; k++){
						
						sum = sum + L[i][k]*L[j][k];
						
					}
					
					L[i][j] = 1.0/L[j][j]*(A[i][j]-sum);
					
				}
				
			}
			
		}
		
		return L;
		
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	double [][] B = {	{4.0, 12.0, -16.0},
    						{12.0, 37.0, -43.0},
    						{-16.0, -43.0, 98.0}
    					};
    	
    	MatrixOperations.print_matrix(B);
    	
    	MatrixOperations.print_matrix(decompose(B));
    	
    }
	
}
