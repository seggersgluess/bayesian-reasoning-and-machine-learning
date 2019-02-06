
public class MatrixOperations {

	
	// returns matrix A of zeros
	public static double [][] matrix(int n_row, int n_cols){
		
		double [][] A = new double [n_row][n_cols];
		
		return A;
		
	}
	
	
	// returns matrix A from a 2D data frame
	public static double [][] fillMatrix(double[][] data){
		
		int n_rows = data.length;
		int n_cols = data[0].length;
		
		double [][] A = MatrixOperations.matrix(n_rows, n_cols);
		
		for(int i = 0; i < n_rows; i++){
			
			for(int j = 0; j < n_cols; j++){
				
				A[i][j] = data[i][j];
				
			}
			
		}
		
		return A;
		
	}
	
	
	// sets a supplied vector c as column into matrix A
	public static double [][] set_column_to_matrix(double[][] A, double [] c, int col_number){
		
		int n_rows = A.length;
		
		for(int i = 0; i < n_rows; i++){
				
			A[i][col_number] = c[i];
				
		}
		
		return A;
		
	}
	
	
	// sets matrix A[start_row, end_row][start_col, end_col] = B
	public static double [][] set_sub_matrix_to_matrix(double [][] A, double [][] B, int start_row, int end_row, int start_col, int end_col){
		
		int n_rows = B.length;
		int n_cols = B[0].length;
		
		if(n_rows != (end_row-start_row+1) || n_cols != (end_col-start_col+1)){
			
			throw new RuntimeException("Mismatch between supplied indices and matrix dimensionality.");
			
		}
		
		for(int i = 0; i < n_rows; i++){
			
			for(int j = 0; j < n_cols; j++){
				
				A[start_row+i][start_col+j] = B[i][j];
				
			}
			
		}
		
		return A;
		
	}
	
	// returns columns of a matrix between start and end indices
	public static double [][] get_sub_matrix_between_column_idxs(double[][] A, int start_col_idx, int end_col_idx){
		
		int n_rows = A.length;
		int n_cols = end_col_idx - start_col_idx + 1;
		
		double [][] sub_matrix = matrix(n_rows, n_cols);
		
		int col_idx = start_col_idx;
		
		for(int i = 0; i < n_cols; i++){
			
			for(int j = 0; j<n_rows; j++){
				
				sub_matrix[j][i] = A[j][col_idx];	
				
			}
				
			col_idx = col_idx + 1;
			
		}
		
		return sub_matrix;
		
	}
	
	
	//sorts matrix columns
	public static double [][] resort_matrix_columns(double[][] A, int [] col_idxs){
		
		int n_rows = A.length;
		int n_cols = col_idxs.length;
		
		double [][] resorted_matrix = matrix(n_rows, n_cols);
				
		for(int i = 0; i < n_rows; i++){
			
			for(int j = 0; j<n_cols; j++){
				
				resorted_matrix[i][j] = A[i][col_idxs[j]];	
				
			}
				
		}
		
		return resorted_matrix;
		
	}
	
	
	// create and return diagonal matrix D with elements of vector d
	public static double [][] diagonal(double [] d){
		
		int n = d.length;
		
		double [][] D = MatrixOperations.matrix(n, n);
		
		for(int i = 0; i < n; i++){
			
			D[i][i] = d[i];
			
		}
				
		return D;
		
	}
	
	
	// create and return diagonal matrix D with elements of vector d
	public static double [][] diagonal(double [][] d){
		
		int n = d.length;
		
		double [][] D = MatrixOperations.matrix(n, n);
		
		for(int i = 0; i < n; i++){
			
			D[i][i] = d[i][1];
			
		}
				
		return D;
		
	}
	
	
	// creates and returns identity matrix I with dimension n x n
	public static double [][] identity(int n){
		
		double [][] I = MatrixOperations.matrix(n, n);
		
		for(int i = 0; i < n; i++){
			
			I[i][i] = 1;
			
		}
				
		return I;
	
	}
	
	
	// creates and returns the transpose of a matrix A
	public static double [][] transpose(double[][] A){
		
		int n_rows = A.length;
		int n_cols = A[0].length;
		
		double [][] A_T = MatrixOperations.matrix(n_cols, n_rows);
		
		for(int i = 0; i < n_rows; i++){
			
			for(int j = 0; j < n_cols; j++){
				
				A_T[j][i] = A[i][j];
				
			}
			
		}
		
		return A_T;
		
	}
	
	
	// returns n x 1 unit vector e
	public static double [][] unit_vector(int n){
		
		double [][] unit_vector = new double [n][1];
		
		for(int i = 0; i < n; i++){
			
			unit_vector[i][0] = 1.0;
			
		}
		
		return unit_vector;
		
	}
	
	
	// creates and returns inverse B^-1 of quadratic matrix B by applying Gaussian-elimination (Gauss-Jordan-algorithm)
	public static double [][] inverse(double [][] B){
				
		int n_rows = B.length;
		int n_cols = B[0].length;
		
		if( n_rows != n_cols){
			
			throw new RuntimeException("Cannot calculate inverse of matrix. No quadratic Matrix supplied.");
			
		}
		
		double [][] A = matrix(n_rows, n_cols);
		A = fillMatrix(B);
		
		double [][] I = MatrixOperations.identity(n_cols);
		double [][] A_inverse = MatrixOperations.matrix(n_rows, n_cols);
		
		//Case of all elements in matrix unequal 0
				
		double [] col = MatrixOperations.get_column_from_matrix(A, 0);
		int [] zero_idxs = Utilities.get_idx(col, 0.0);
					
		if(zero_idxs[0] != -1){
				
			int non_zero_idx = 0;
			int n_zero_idxs = zero_idxs.length;	
			
			int [] non_zero_idxs = new int [n_rows-n_zero_idxs];				
			int [] resorted_idxs = new int [n_rows];
				
			for(int i = 0; i < n_rows; i++){
					
				int check = 0;
					
				for(int j = 0; j < n_zero_idxs; j++){
																		
					if(zero_idxs[j] == i){
							
						check = 1;
						break;
							
					}
															
				}
					
				if(check == 0){
						
					non_zero_idxs[non_zero_idx] = i;
					non_zero_idx = non_zero_idx + 1;
						
				}
					
			}
				
			double [][] A_new = MatrixOperations.matrix(n_rows, n_cols);
			double [][] I_new = MatrixOperations.matrix(n_rows, n_cols);
				
			for(int i = 0; i < n_rows; i++){
					
				int zero_idx_counter = 0;
					
				if(i < non_zero_idxs.length){
						
					resorted_idxs[i] = non_zero_idxs[i];
						
				}else{
						
					resorted_idxs[i] = zero_idxs[zero_idx_counter];
						
					zero_idx_counter = zero_idx_counter + 1;
						
				}
					
				for(int j = 0; j < n_rows; j++){
						
					double [] A_row = MatrixOperations.get_row_from_matrix(A, resorted_idxs[i]);
						
					A_new[i][j] = A_row[j];
						
					double [] I_row = MatrixOperations.get_row_from_matrix(I,resorted_idxs[i]);
						
					I_new[i][j] = I_row[j];
						
				}
										
			}
							
			A = A_new;
			I = I_new;
				
		}
							
		//Two steps: 1. Create lower triangular, 2. Create upper triangular
		
		//Step 1:
		for(int i = 0; i < n_cols-1; i++){
			
			int n_row_counter = n_rows - i - 1;
			
			for(int j = 0; j < n_row_counter; j++){
				
				int row_idx_1 = n_rows-j-1;
				
				if(A[row_idx_1][i] != 0.0){
					
					int row_idx_2 = row_idx_1-1;
						
					double element_1 = A[row_idx_1][i];
					double element_2 = A[row_idx_2][i];
					
					for(int k = 0; k < n_cols; k++){
											
						A[row_idx_1][k] = A[row_idx_1][k]*element_2;
						A[row_idx_2][k] = A[row_idx_2][k]*element_1;
							
						double a_1 = A[row_idx_1][i];
						double a_2 = A[row_idx_2][i];
																					
						if((a_1 >= 0 && a_2 >= 0) || (a_1 <= 0 && a_1 <= 0)){
							
							A[row_idx_1][k] = A[row_idx_1][k] - A[row_idx_2][k];
							
						}else{
							
							A[row_idx_1][k] = A[row_idx_1][k] + A[row_idx_2][k];
							
						}
							
						I[row_idx_1][k] = I[row_idx_1][k]*element_2;
						I[row_idx_2][k] = I[row_idx_2][k]*element_1;
						
						if((a_1 >= 0 && a_2 >= 0) || (a_1 <= 0 && a_1 <= 0)){
							
							I[row_idx_1][k] = I[row_idx_1][k] - I[row_idx_2][k];
							
						}else{
							
							I[row_idx_1][k] = I[row_idx_1][k] + I[row_idx_2][k];
							
						}
										
					}
					
				}
												
			}
			
		}
		
		//Step 2:
		for(int i = 0; i < n_cols-1; i++){
			
			int n_row_counter = n_rows - i - 1;
			
			int col_idx = n_cols-i-1;
			
			for(int j = 0; j < n_row_counter; j++){
				
				int row_idx_1 = j;
				
				if(Math.abs(A[row_idx_1][col_idx]) != 0.0){
					
					int row_idx_2 = row_idx_1;	
					
					for(int k = 1; k < n_rows; k++){
							
						row_idx_2 = row_idx_2 + 1;
												
						if(Math.abs(A[row_idx_2][col_idx]) != 0.0){
							
							break;
								
						}
											
					}
													
					double element_1 = A[row_idx_1][col_idx];
					double element_2 = A[row_idx_2][col_idx];
									
					for(int k = 0; k < n_cols; k++){
											
						A[row_idx_1][k] = A[row_idx_1][k]*element_2;
						A[row_idx_2][k] = A[row_idx_2][k]*element_1;
						
						double a_1 = 0.0;
						double a_2 = 0.0;
						
						if(k == 0){
							
							a_1 = A[row_idx_1][col_idx]*element_2;
							a_2 = A[row_idx_2][col_idx]*element_1;
							
						}
											
						if((a_1 >= 0 && a_2 >= 0) || ((a_1 <= 0 && a_2 <= 0))){
							
							A[row_idx_1][k] = A[row_idx_1][k] - A[row_idx_2][k];
							
						}else{
							
							A[row_idx_1][k] = A[row_idx_1][k] + A[row_idx_2][k];
							
						}
						
						I[row_idx_1][k] = I[row_idx_1][k]*element_2;
						I[row_idx_2][k] = I[row_idx_2][k]*element_1;
						
						if((a_1 >= 0 && a_2 >= 0) || ((a_1 <= 0 && a_2 <= 0))){
							
							I[row_idx_1][k] = I[row_idx_1][k] - I[row_idx_2][k];
							
						}else{
							
							I[row_idx_1][k] = I[row_idx_1][k] + I[row_idx_2][k];
							
						}
						
					}
					
				}
																
			}
			
		}
		
		for(int i = 0; i < n_rows; i++){
			
			for(int j=0; j < n_cols; j++){
				
				I[i][j] = I[i][j]/A[i][i];
				
			}
			
		}
		
		A_inverse = I;
			
		return A_inverse;
		
	}
	
	
	// returns C = A + B
	public static double [][] add(double [][] A, double [][] B){
		
		int n_rows_a = A.length;
		int n_rows_b = B.length;
		int n_cols_a = A[0].length;
		int n_cols_b = B[0].length;
		
		if( n_rows_a != n_rows_b || n_cols_a != n_cols_b){
			
			throw new RuntimeException("Illegal matrix dimensions.");
			
		}
		
		double [][] C = MatrixOperations.matrix(n_rows_a, n_cols_a);
		
		for(int i = 0; i < n_rows_a; i++){
			
			for(int j = 0; j < n_cols_a; j++){
				
				C[i][j] = A[i][j] + B[i][j];
				
			}
			
		}
				
		return C;
		
	}
	
	
	// returns C = A + B
	public static double [] add_vectors(double [] a, double [] b){
		
		int n_rows_a = a.length;
		int n_rows_b = b.length;
		
		if( n_rows_a != n_rows_b){
			
			throw new RuntimeException("Illegal matrix dimensions.");
			
		}
		
		double [] c = new double [n_rows_a];
		
		for(int i = 0; i < n_rows_a; i++){
						
			c[i] = a[i] + b[i];
				
		}
				
		return c;
		
	}
	
	
	// returns b = s*a
	public static double [] scalar_vector_multiplication(double s, double[] a){
		
		int n_rows = a.length;
		
		double [] b = new double [n_rows];
		
		for(int i = 0; i < n_rows; i++){
							
			b[i] = s*a[i];
							
		}
		
		return b;
		
	}
	
	
	// returns B = s*A
	public static double [][] scalar_multiplication(double s, double[][] A){
		
		int n_rows = A.length;
		int n_cols = A[0].length;
		
		double [][] B = MatrixOperations.matrix(n_rows, n_cols);
		
		for(int i = 0; i < n_rows; i++){
			
			for(int j = 0; j < n_cols; j++){
				
				B[i][j] = s*A[i][j];
				
			}
			
		}
		
		return B;
		
	}
	
	
	// returns B = A*B
	public static double [][] multiplication(double[][] A, double[][] B){
		
		int n_rows_a = A.length;
		int n_rows_b = B.length;
		int n_cols_a = A[0].length;
		int n_cols_b = B[0].length;
		
		if( n_cols_a != n_rows_b){
			
			throw new RuntimeException("Illegal matrix dimensions.");
			
		}
		
		double [][] C = MatrixOperations.matrix(n_rows_a, n_cols_b);
		
		for(int i = 0; i < n_rows_a; i++){
			
			for(int j = 0; j < n_cols_b; j++){
											
				for(int k = 0; k < n_cols_a; k++){
						
					C[i][j] += A[i][k]*B[k][j];
						
				}
				
			}
			
		}
		
		return C;
		
	}

	
	// returns n x 1 vector b
	public static double [][] convArrayToVec(double [] a){
		
		int n = a.length;
		double [][] b = new double [n][1];
		
		for(int i = 0; i < n; i++){
			
			b[i][0] = a[i];
			
		}
		
		return b;
		
	}
	
	
	// returns c = A*b
	public static double [] multiplyMatrixWithVec(double [][] A, double [] b){
		
		int n_rows_A = A.length;
		int n_cols_A = A[0].length;
		int n_rows_b = b.length;
		
		if( n_cols_A != n_rows_b){
			
			throw new RuntimeException("No valid dimensions.");
			
		}
		
		double [] c = new double [n_rows_A];
		
		for(int i = 0; i < n_rows_A; i++){
			
			for(int j = 0; j < n_cols_A; j++){
				
				c[i] += A[i][j]*b[j];
				
			}
			
		}
		
		return c;
		
	}
	
	
	// binds column x vector to matrix A
	public static double [][] cbind(double [][] A, double [] x){
		
		int n_rows = x.length;
		int n_cols = A[0].length + 1;
		
		double [][] merged_matrix = matrix(n_rows, n_cols);
		
		for(int i=0; i<n_cols; i++){
			
			for(int j=0; j<n_rows; j++){
				
				if(i == 0){
					
					merged_matrix[j][i] = x[j];
					
				}else{
					
					merged_matrix[j][i] = A[j][(i-1)];
					
				}
					
			}
				
		}
		
		return merged_matrix;
		
	}
	
	
	// binds matrix X vector to matrix A
	public static double [][] cbind(double [][] A, double [][] X){
		
		int n_rows = X.length;
		int n_cols = A[0].length + X[0].length;
		
		double [][] merged_matrix = matrix(n_rows, n_cols);
		
		for(int i=0; i<n_cols; i++){
			
			for(int j=0; j<n_rows; j++){
				
				if(i < X[0].length){
					
					merged_matrix[j][i] = X[j][i];
					
				}else{
					
					merged_matrix[j][i] = A[j][(i-X[0].length)];
					
				}
					
			}
				
		}
		
		return merged_matrix;
		
	}
	
	
	// returns diagonal from matrix A
	public static double [][] get_diagonal_from_matrix(double [][] A){
		
		int n_cols = A[0].length;
		int n_rows = A.length;
			
		if(n_rows != n_cols){
			
			throw new RuntimeException("No quadratic matrix supplied.");
			
		}
		
		double [][] diagonal = new double [n_rows][1];
		
		for(int i=0; i<n_rows; i++){
			
			diagonal[i][0] = A[i][i];
			
		}
		
		return diagonal;
		
	}
	
	
	// returns specific column of a matrix
	public static double [] get_column_from_matrix(double [][] A, int col_number){
		
		int n_cols = A[0].length;
		int n_rows = A.length;
		double [] column_vec = new double [n_rows];
		
		if( col_number > n_cols){
			
			throw new RuntimeException("Mismatch between column number and number of matrix columns.");
			
		}
		
		for(int i=0; i<n_rows; i++){
			
			column_vec[i] = A[i][col_number];
			
		}
		
		return column_vec;
		
	}
	
	
	// returns specific row of a matrix
	public static double [] get_row_from_matrix(double [][] A, int row_number){
		
		int n_cols = A[0].length;
		int n_rows = A.length;
		double [] row_vec = new double [n_cols];
		
		if( row_number > n_rows){
			
			throw new RuntimeException("Mismatch between row number and number of matrix rows.");
			
		}
		
		for(int i=0; i<n_cols; i++){
			
			row_vec[i] = A[row_number][i];
			
		}
		
		return row_vec;
		
	}
	
	
	// returns a sub vector between start and end index of a supplied vector
	public static double [] get_double_sub_vec(double [] x, int start_idx, int end_idx){
		
		if(x.length < start_idx){
			
			throw new RuntimeException("Start index larger than vector length.");
			
		}
		
		if(x.length < end_idx){
			
			throw new RuntimeException("End index larger than vector length.");
			
		}
		
		int n_sub_elements = end_idx - start_idx + 1;
		
		double [] sub_vec = new double [n_sub_elements];
		
		int idx = 0	;	
				
		for(int i=start_idx; i<(end_idx+1); i++){
			
			sub_vec[idx] = x[i];
			
			idx  =  idx + 1;
					
		}
		
		return sub_vec;
		
	}
	
	
	// returns a sub vector between start and end index of a supplied vector
	public static double [][] get_double_sub_vec(double [][] x, int start_idx, int end_idx){
		
		if(x.length < start_idx){
			
			throw new RuntimeException("Start index larger than vector length.");
			
		}
		
		if(x.length < end_idx){
			
			throw new RuntimeException("End index larger than vector length.");
			
		}
		
		int n_sub_elements = end_idx - start_idx + 1;
		
		double [][] sub_vec = new double [n_sub_elements][1];
		
		int idx = 0	;	
				
		for(int i = start_idx; i < (end_idx+1); i++){
			
			sub_vec[idx][0] = x[i][0];
			
			idx  =  idx + 1;
					
		}
		
		return sub_vec;
		
	}
	
	
	// returns a vector with specific elements of a supplied vector
	public static double [] get_double_sub_vec_4_indices(double [] x, int [] idxs){
		
		int n_idxs = idxs.length;
		double [] sub_vec = new double [n_idxs];
		
		for(int i=0; i<n_idxs; i++){
			
			sub_vec[i] = x[idxs[i]];
			
		}
				
		return sub_vec;
		
	}
	
	
	// returns (complementary) sub vector with elements not in index set idxs
	public static double [] get_double_comp_sub_vec_4_indices(double [] x, int [] idxs){
		
		int n_comps = x.length-idxs.length;
		double [] comp_sub_vec = new double [n_comps];
		int idx = 0;
		
		double [] dbl_idxs = Utilities.convert_int_to_double_array(idxs);
		
		for(int i=0; i<x.length; i++){
			
			int [] valid_idxs = Utilities.get_idx(dbl_idxs, (double) i);
				
			if(valid_idxs[0] == -1){
					
				comp_sub_vec[idx] = x[i];
				idx = idx + 1;
					
			}
				 			
		}
				
		return comp_sub_vec;
		
	}
	
	
	// returns a integer sub vector between start and end index of a supplied vector
	public static int [] get_int_sub_vec(int [] x, int start_idx, int end_idx){
		
		if(x.length < start_idx){
			
			throw new RuntimeException("Start index larger than vector length.");
			
		}
		
		if(x.length < end_idx){
			
			throw new RuntimeException("End index larger than vector length.");
			
		}
		
		int n_sub_elements = end_idx - start_idx + 1;
		
		int [] sub_vec = new int [n_sub_elements];
		
		int idx = 0	;	
				
		for(int i=start_idx; i<(end_idx+1); i++){
			
			sub_vec[idx] = x[i];
			
			idx  =  idx + 1;
					
		}
		
		return sub_vec;
		
	}
	
	
	// reverse ordering of column vector x
	public static double [][] reverse(double [][] x){
		
		int n_cols = x[0].length;
		
		double [][] reverse_x = new double [1][n_cols];
		
		for(int i=0; i<n_cols; i++){
			
			reverse_x[0][i] = x[0][n_cols-(i+1)];
			
		}
				
		return reverse_x;
		
	}
	
	
	// prints supplied double matrix A
	public static void print_matrix(double [][] A){
		
		for (int i = 0; i < A.length; i++) {
			
		    for (int j = 0; j < A[i].length; j++) {
		    	
		        System.out.print(A[i][j] + " ");
		        
		    }
		    
		    System.out.println();
		}
		
	}
	
	
	// prints supplied string matrix A
	public static void print_matrix(String [][] A){
		
		for (int i = 0; i < A.length; i++) {
			
		    for (int j = 0; j < A[i].length; j++) {
		    	
		        System.out.print(A[i][j] + " ");
		        
		    }
		    
		    System.out.println();
		}
		
	}
	

	// prints supplied vector x
	public static void print_vector(double [] x){
		
		for (int i = 0; i < x.length; i++) {
					    	
		    System.out.print(x[i] + " ");
		         		    
		}
		
		System.out.println();
		
	}
	
	
	// combines two vectors a and b in vector c
	public static double [] combine_vectors(double [] a, double [] b){
		
		int n_a = a.length;
		int n_b = b.length;
		
		int n_elements = n_a + n_b;
		int b_idx = 0;
		
		double [] c = new double [n_elements];
		
		for(int i = 0; i < n_a; i++){
			
			c[i] = a[i];
			
		}
				
		for(int i = n_a; i < n_elements; i++){
			
			c[i] = b[b_idx];
			
			b_idx = b_idx + 1;
			
		}
		
		return c;
		
	}
	
	
	// returns the Euclidian norm of vector a
	public static double euclidian(double [] a){
		
		int n = a.length;
		double sum = 0.0;
		
		for(int i = 0; i < n; i++){
			
			sum += a[i]*a[i];
			
		}
		
		sum = Math.sqrt(sum);
		
		return sum;
		
	}
	
	
	// returns the Euclidian norm of vector a
	public static double euclidian(double [][] a){
		
		int n = a.length;
		
		if(a[0].length != 1){
			
			throw new RuntimeException("No column vector supplied for Euclidian norm.");
			
		}
		
		double sum = 0.0;
		
		for(int i = 0; i < n; i++){
			
			sum += a[i][0]*a[i][0];
			
		}
		
		sum = Math.sqrt(sum);
		
		return sum;
		
	}
	
	
    // test client
    public static void main(String[] args) {
    	
    	//double [][] A = MatrixOperations.matrix(4, 4);
    	//A[0][0] = 0;
    	//A[0][1] = 0;
    	//A[0][2] = -2;
    	//A[0][3] = 12;
    	//A[1][0] = 4;
    	//A[1][1] = 0;
    	//A[1][2] = 16;
    	//A[1][3] = -15;
    	//A[2][0] = 7;
    	//A[2][1] = 8;
    	//A[2][2] = -9;
    	//A[2][3] = 8;
    	//A[3][0] = 10;
    	//A[3][1] = 81;
    	//A[3][2] = -9;
    	//A[3][3] = 0;
        
        //double [][] A = {{1,1},{1,1}};
        
        //double [][] A_inverse = inverse(A);
        
    	//print_matrix(A_inverse);

    	//double [] a = {18,3};
    	//double [] b = {2,4};
    	
    	//double [][] a_new = convArrayToVec(a);
    	//double [][] b_new = convArrayToVec(b);
    	
    	
    	double [][] A = new double [4][4];
    	double [][] B = {{1,2},{3,4}};
    	
    	print_matrix(set_sub_matrix_to_matrix(A, B, 2, 3, 2, 3));
    	
    	//print_vector(multiplyMatrixWithVec(A,b));
    	
    	//print_matrix(multiplication(a_new,transpose(b_new)));
    	//System.out.println(Math.pow(2.0,3.0));
    	
    	//print_matrix(a_new);
    	
    	//print_vector(get_row_from_matrix(C,1));
    	
    }
	
	
}
