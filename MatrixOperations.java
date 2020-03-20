package Mathematics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Singular;
import org.jblas.Solve;

import Utilities.Utilities;

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
	
	
	// returns columns of a matrix between start and end indices
	public static double [][] get_sub_matrix_4_column_idxs(double[][] A, int [] colIdxs){
		
		int n_rows = A.length;
		int n_cols = colIdxs.length;
		
		double [][] sub_matrix = matrix(n_rows, n_cols);
		
		for(int i = 0; i<n_cols; i++){		
			for(int j = 0; j<n_rows; j++){				
				sub_matrix[j][i] = A[j][colIdxs[i]];					
			}			
		}
		
		return sub_matrix;	
	}
	
	
	// returns columns of a matrix between start and end indices
	public static double [][] get_sub_matrix_4_row_idxs(double[][] A, int [] rowIdxs){
		
		int n_cols = A[0].length;
		int n_rows = rowIdxs.length;
		
		double [][] sub_matrix = matrix(n_rows, n_cols);
		
		for(int i = 0; i<n_cols; i++){		
			for(int j = 0; j<n_rows; j++){				
				sub_matrix[j][i] = A[rowIdxs[j]][i];					
			}			
		}
		
		return sub_matrix;	
	}
	
	
	// returns rows of a matrix between start and end indices
	public static double [][] get_sub_matrix_between_row_idxs(double[][] A, int start_row_idx, int end_row_idx){
		
		int n_rows = end_row_idx - start_row_idx + 1;
		int n_cols = A[0].length;
		
		double [][] sub_matrix = matrix(n_rows, n_cols);
		
		int row_idx = start_row_idx;
		
		for(int i = 0; i < n_rows; i++){			
			for(int j = 0; j<n_cols; j++){				
				sub_matrix[i][j] = A[row_idx][j];					
			}
			
			row_idx = row_idx + 1;
			
		}
		
		return sub_matrix;	
	}
	
	
	// returns rows of a matrix between start and end indices
	public static double [][] get_sub_matrix_between_row_and_col_idxs(double[][] A, int start_row_idx, int end_row_idx, int start_col_idx, int end_col_idx){
		
		int n_rows = end_row_idx - start_row_idx + 1;
		int n_cols = end_col_idx - start_col_idx + 1;
		
		double [][] sub_matrix = matrix(n_rows, n_cols);
		
		int row_idx = start_row_idx;
		
		for(int i = 0; i < n_rows; i++){	
			
			int col_idx = start_col_idx;
			
			for(int j = 0; j<n_cols; j++){				
				sub_matrix[i][j] = A[row_idx][col_idx];	
				col_idx++;
			}
			
			row_idx++;
			
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
			
			D[i][i] = d[i][0];
			
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
	
	
	/**
	 * Method for generating nx1 unit (column) vector
	 * @param n number of rows
	 * @return nx1 unit vector
	 */
	public static double [][] unit_vector(int n){
		
		double [][] unit_vector = new double [n][1];
		
		for(int i = 0; i < n; i++){
			
			unit_vector[i][0] = 1.0;
			
		}
		
		return unit_vector;		
	}
	
	
	/**
	 * Method for calculating the determinant of a given matrix
	 *
	 * @param matrix matrix for which determinant should be calculated
	 *
	 * @return determinant of given matrix
	 */
	public static double determinant (double[][] matrix) {
		
		double temporary[][];
		double result = 0;

		if (matrix.length == 1) {
			result = matrix[0][0];
			return (result);
		}

		if (matrix.length == 2) {
			result = ((matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]));
			return (result);
		}

		for (int i = 0; i < matrix[0].length; i++) {
			temporary = new double[matrix.length - 1][matrix[0].length-1];

			for (int j = 1; j < matrix.length; j++) {
				for (int k = 0; k < matrix[0].length; k++) {
					if (k < i) {
						temporary[j-1][k] = matrix[j][k];
					} else if (k > i) {
						temporary[j-1][k-1] = matrix[j][k];
					}
				}
			}

			result += matrix[0][i]*Math.pow (-1, (double) i)*determinant (temporary);
		}
		
		return (result);
	}
	
	
	@SuppressWarnings("static-access")
	public static double [][] inverse_fast(double [][] X) {
		
		DoubleMatrix X_new = new DoubleMatrix(X);		
		DoubleMatrix X_new_inv = Solve.solve(X_new, X_new.eye(X_new.rows));
		double [][] X_inv = X_new_inv.toArray2();
		
		return X_inv;		
	}
	
	
	public static double [][] inverse(double [][] X){
		
		int n_rows = X.length;
		int n_cols = X[0].length;
		
		if( n_rows != n_cols){			
			throw new RuntimeException("Cannot calculate inverse of matrix. No quadratic Matrix supplied.");			
		}
		
		double [][] A = matrix(n_rows, n_cols);
		A = fillMatrix(X);
		
		if(is_diagonal(A) == true){
			
			for(int i=0; i<n_rows; i++){				
				A[i][i] = 1.0/A[i][i];				
			}
		
			return A;
			
		}	
			
		double d = 1.0;
		
		for(int p=0; p<n_rows; p++){
			
			double pivot = A[p][p];
			
			if(pivot == 0.0){			
				throw new RuntimeException("Inverse matrix cannot calculated.");						
			}
			
			d = d*pivot;
			
			for(int i=0; i<n_rows; i++){				
				A[i][p] =  -A[i][p]/pivot;										
			}
			
			for(int i=0; i<n_rows; i++){
		
				if(i != p){					
					for(int j=0; j<n_cols; j++){						
						if(j != p){							
							A[i][j] =  A[i][j] + A[p][j]*A[i][p]; 							
						}								
					}						
				}				
			}
			
			for(int j=0; j<n_cols; j++){						
				A[p][j] = A[p][j]/pivot;												
			}
			
			A[p][p] = 1.0/pivot;
			
		}
		
		return A;		
	}
	
	
	// checks if matrix A is diagonal
	public static boolean is_diagonal(double [][] A){
		
		int nCols = A[0].length;
		
		int zeroCount = 0;
		
		boolean isDiagonal = true;
		
		for(int i=0; i<nCols; i++){
			
			double [] col = MatrixOperations.get_column_from_matrix(A, i);
			
			if(i==0){
				
				int[] idx = Utilities.get_idx(MatrixOperations.get_double_sub_vec(col, 1, nCols-1), 0.0);
				
				if(idx[0] != -1){
					
					zeroCount = idx.length;
					
				}else{
					
					zeroCount = 0;
					
				}
					
			}else{
				
				if(i != nCols-1){
					
					int[] idx = Utilities.get_idx(MatrixOperations.get_double_sub_vec(col, 0, i-1), 0.0);
					
					if(idx[0] != -1){
						
						zeroCount = idx.length;
						
					}else{
						
						zeroCount = 0;
						
					}
					
					idx = Utilities.get_idx(MatrixOperations.get_double_sub_vec(col, i+1, nCols-1), 0.0);
					
					if(idx[0] != -1){
						
						zeroCount = idx.length + zeroCount;
						
					}else{
						
						zeroCount = 0;
						
					}
					
				}else{
					
					int[] idx = Utilities.get_idx(MatrixOperations.get_double_sub_vec(col, 0, nCols-2), 0.0);
					
					if(idx[0] != -1){
						
						zeroCount = idx.length;
						
					}else{
						
						zeroCount = 0;
						
					}
					
				}
				
			}
			
			if(zeroCount != nCols-1){
				
				isDiagonal = false;
				break;
				
			}
				
		}
		
		return isDiagonal;		
	}
	

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
	
	
	// returns C = A - B
	public static double [][] substract(double [][] A, double [][] B){
		
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
				
				C[i][j] = A[i][j] - B[i][j];
				
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
	
	
	// returns B = A*B
	public static double [][] multiplication_fast(double[][] A, double[][] B){
		
		int n_rows_b = B.length;
		int n_cols_a = A[0].length;
		
		if( n_cols_a != n_rows_b){			
			throw new RuntimeException("Illegal matrix dimensions.");			
		}
		
		DoubleMatrix A_dbl = new DoubleMatrix(A);
		DoubleMatrix B_dbl = new DoubleMatrix(B);
		
		DoubleMatrix C_dbl = A_dbl.mmul(B_dbl);
		
		double [][] C = C_dbl.toArray2();
		
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
	
	
	// returns n array b
	public static double [] convVecToArray(double [][] a){
		
		int n = a.length;
		double [] b = new double [n];
		
		for(int i = 0; i < n; i++){			
			b[i] = a[i][0];			
		}
		
		return b;		
	}
	
	
	public static double [][] convNumberToArray(double x) {
		double [][] array = new double [1][1];
		array[0][0] = x;
		
		return array;
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
	
	
	public static double [][] rbind(double [][] A, double [][] X){
		
		int n_x_rows = X.length;
		int n_rows = A.length + n_x_rows;
		int n_cols = A[0].length;
		
		if(n_cols != X[0].length) {
			throw new RuntimeException("Supplied matrices have unequal number of columns.");
		}
		
		double [][] merged_matrix = matrix(n_rows, n_cols);
		
		for(int i=0; i<n_cols; i++){		
			for(int j=0; j<n_rows; j++){			
				if(j<n_x_rows){				
					merged_matrix[j][i] = X[j][i];				
				}else{				
					merged_matrix[j][i] = A[(j-n_x_rows)][i];				
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
		
		if( col_number>n_cols){			
			throw new RuntimeException("Mismatch between column number and number of matrix columns.");			
		}
		
		for(int i=0; i<n_rows; i++){		
			column_vec[i] = A[i][col_number];			
		}
		
		return column_vec;		
	}
	
	
	// returns specific column of a matrix
	public static double [][] get_column_vec_from_matrix(double [][] A, int col_number){
		
		int n_cols = A[0].length;
		int n_rows = A.length;
		double [][] column_vec = new double [n_rows][1];
		
		if( col_number>n_cols){			
			throw new RuntimeException("Mismatch between column number and number of matrix columns.");			
		}
		
		for(int i=0; i<n_rows; i++){			
			column_vec[i][0] = A[i][col_number];			
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
	
	
	// returns specific row of a matrix
	public static double [][] get_row_vec_from_matrix(double [][] A, int row_number){
		
		int n_cols = A[0].length;
		int n_rows = A.length;
		double [][] row_vec = new double [n_cols][1];
		
		if( row_number > n_rows){			
			throw new RuntimeException("Mismatch between row number and number of matrix rows.");			
		}
		
		for(int i=0; i<n_cols; i++){
			
			row_vec[i][0] = A[row_number][i];
			
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
	
	
	// returns (complementary) sub vector as list with elements not in index set idxs
	public static List<Double> get_double_comp_sub_vec_4_indices(List<Double> x, int [] idxs){
		
		int n = x.size();
		int n_comps = n-idxs.length;
		List<Double> comp_sub_vec = new ArrayList<Double>(n_comps);
		int idx = 0;
		
		double [] dbl_idxs = Utilities.convert_int_to_double_array(idxs);
		
		for(int i=0; i<n; i++){
			
			int [] valid_idxs = Utilities.get_idx(dbl_idxs, (double) i);
				
			if(valid_idxs[0] == -1){					
				comp_sub_vec.add(x.get(i));
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
			idx++;					
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
			b_idx++;			
		}
		
		return c;	
	}
	
	
	// combines two vectors a and b in vector c
	public static int [] combine_vectors(int [] a, int [] b){
		
		int n_a = a.length;
		int n_b = b.length;
		
		int n_elements = n_a + n_b;
		int b_idx = 0;
		
		int [] c = new int [n_elements];
		
		for(int i = 0; i < n_a; i++){			
			c[i] = a[i];			
		}
				
		for(int i = n_a; i < n_elements; i++){			
			c[i] = b[b_idx];			
			b_idx = b_idx + 1;
			
		}
		
		return c;		
	}
	
	
	// returns trace of matrix (sum of diagonal elements)
	public static double trace(double [][] A){
		
		if(A.length != A[0].length){
			throw new RuntimeException("Supplied matrix is not a square matrix.");
		}
		
		double trace = 0.0;
		int n = A.length;
		
		for(int i=0; i<n; i++){
			trace = trace + A[i][i];			
		}
		
		return trace;
		
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
	
	
	//Element wise multiplication of two nx1 vectors a, b
	public static double [][] multiplyVecByElement_old(double [][] a, double [][] b){
		
		if(a.length != b.length){
			throw new RuntimeException("Supplied vectors are of unequal length.");
		}
		
		if(a[0].length != 1 || b[0].length != 1){
			throw new RuntimeException("No column vectors supplied.");
			
		}
		
		int nRows = a.length;
		
		double [][] c = new double [nRows][1];
		
		for(int i=0; i<nRows; i++){
			c[i][0] = a[i][0]*b[i][0];
		}
		
		return c;
		
	}
	
	
	public static double [][] multiplyVecByElement(double [][] a, double [][] b){
		
		if(a.length != b.length){
			throw new RuntimeException("Supplied vectors are of unequal length.");
		}
		
		if(a[0].length != 1 || b[0].length != 1){
			throw new RuntimeException("No column vectors supplied.");
			
		}
		
		DoubleMatrix a_dbl = new DoubleMatrix(a);
		DoubleMatrix b_dbl = new DoubleMatrix(b);
		
		double [][] c = a_dbl.mul(b_dbl).toArray2();
		
		return c;
		
	}
	
	
	//Element wise division of two nx1 vectors a, b
	public static double [][] divideVecByElement(double [][] a, double [][] b){
		
		if(a.length != b.length){
			throw new RuntimeException("Supplied vectors are of unequal length.");
		}
		
		if(a[0].length != 1 || b[0].length != 1){
			throw new RuntimeException("No column vectors supplied.");
			
		}
		
		int nRows = a.length;
		
		double [][] c = new double [nRows][1];
		
		for(int i=0; i<nRows; i++){
			c[i][0] = a[i][0]/b[i][0];
		}
		
		return c;		
	}
	
	
	// vec-operator for vectorizing matrix A
	public static double [] vec(double [][] A){
		
		int nRows = A.length;
		int nCols = A[0].length;
		
		int idx   = 0;
		
		double [] vec = new double [(nRows*nCols)];
		
		for(int i=0; i<nCols; i++){
			
			for(int j=0; j<nRows; j++){				
				vec[idx] = A[j][i];				
				idx++;				
			}			
		}
		
		return vec;		
	}
	
	
	// vec-operator for vectorizing matrix A
	public static List<Double> vecAsList(double [][] A){
		
		int nRows = A.length;
		int nCols = A[0].length;
		
		List<Double> vec = new ArrayList<Double>(nRows*nCols);
		
		for(int i=0; i<nCols; i++){			
			for(int j=0; j<nRows; j++){				
				vec.add(A[j][i]);								
			}			
		}
		
		return vec;	
	}
	
	// vec-operator for vectorizing matrix A
	public static double [][] vecAs2dArray(double [][] A){
		
		int nRows = A.length;
		int nCols = A[0].length;
		
		int idx   = 0;
		
		double [][] vec = new double [(nRows*nCols)][1];
		
		for(int i=0; i<nCols; i++){			
			for(int j=0; j<nRows; j++){				
				vec[idx][0] = A[j][i];				
				idx++;				
			}			
		}
		
		return vec;		
	}
	
		
	// returns matrix A from column vector x
	public static double [][] get_matrix_from_vec(double [] x, int n_rows, int n_cols){
		
		double [][] A = new double [n_rows][n_cols];
		
		int idx = 0;
		
		for(int i=0; i<n_cols; i++){			
			for(int j=0; j<n_rows; j++){				
				A[j][i] = x[idx];				
				idx++;						
			}			
		}
		
		return A;		
	}
	
	
	// returns matrix A from column vector x
	public static double [][] get_matrix_from_vec(double [][] x, int n_rows, int n_cols){
		
		double [][] A = new double [n_rows][n_cols];
		
		int idx = 0;
		
		for(int i=0; i<n_cols; i++){			
			for(int j=0; j<n_rows; j++){				
				A[j][i] = x[idx][0];				
				idx++;						
			}			
		}
		
		return A;		
	}
	
	
	// returns matrix A from column vector x
	public static double [][] get_matrix_from_vec(List<Double> x, int n_rows, int n_cols){
		
		double [][] A = new double [n_rows][n_cols];
		
		int idx = 0;
		
		for(int i=0; i<n_cols; i++){			
			for(int j=0; j<n_rows; j++){				
				A[j][i] = x.get(idx);				
				idx++;						
			}			
		}
		
		return A;		
	}
	
	
	// returns Kronecker product of matrices A and B
	public static double [][] kronecker(double [][] A, double [][] B){
		
		int n_rows_A = A.length;
	    int n_cols_A = A[0].length;
		int n_rows_B = B.length;
		int n_cols_B = B[0].length;
	    
		double [][] C = new double [n_rows_A*n_rows_B][n_cols_A*n_cols_B];
		
		for(int i=0; i<n_rows_A;i++){			
			for(int j=0; j<n_cols_A; j++){				
				for(int k=0; k<n_rows_B; k++){					
					for(int l=0; l<n_cols_B; l++){						
						C[(i*n_rows_B)+k][(j*n_cols_B)+l] = A[i][j]*B[k][l];						
					}									
				}					
			}				
		}
		
		return C;
	}
	
	
	 /**
     * Method for calculating eigendecomposition of a symmetric matrix X.
     * (JBLAS library is used)
     * @param X symmetric double matrix
     * @return HashMap of double arrays with keys eigenvalues and eigenvectors. 
     * Key eigenvalues returns diagonal matrix with eigenvalues on its diagonal.
     * Key eigenvectors returns matrix with eigenvectors in its columns)
     */
	@SuppressWarnings("static-access")
	public static HashMap<String, double [][]> get_eigen_dec_4_symmetric_matrix(double [][] X) {
		
		HashMap<String, double [][]> eigenDecRes = new HashMap<String, double [][]>();
		
		Eigen eigen = new Eigen();
		DoubleMatrix X_new = new DoubleMatrix(X);	
		DoubleMatrix [] eigenDec = eigen.symmetricEigenvectors(X_new);
		
		eigenDecRes.put("eigenvectors", eigenDec[0].toArray2());
		eigenDecRes.put("eigenvalues", eigenDec[1].toArray2());
		
		return eigenDecRes;
	}
	
	
	 /**
     * Method for calculating square root X^1/2 of a square matrix X by
     * eigendecomposition (JBLAS library is used)
     * @param X square double matrix
     * @return double [][] X^1/2. 
     */
	public static double [][] get_square_root_of_matrix(double [][] X) {
		
		int n_vars = X.length;
		
		if(n_vars != X[0].length) {
			throw new RuntimeException("No square matrix supplied. Can´t compute square root of matrix.");
		}
		
		HashMap<String, double [][]> eigenDec = get_eigen_dec_4_symmetric_matrix(X);
		
		double [][] e_val = eigenDec.get("eigenvalues");
		double [][] e_vec = eigenDec.get("eigenvectors");
		
		for(int i=0; i<n_vars; i++) {
			e_val[i][i] = Math.sqrt(e_val[i][i]);
		}
		
		double [][] sq_root_X = MatrixOperations.multiplication(MatrixOperations.multiplication(e_vec, e_val),MatrixOperations.transpose(e_vec));
		
		return sq_root_X;		
	}
	
	
	 /**
     * Method for calculating inverse square root X^(-1/2) of a square matrix X by
     * eigendecomposition (JBLAS library is used)
     * @param X square double matrix
     * @return double [][] X^(-1/2). 
     */
	public static double [][] get_inv_square_root_of_matrix(double [][] X) {
		
		int n_vars = X.length;
		
		if(n_vars != X[0].length) {
			throw new RuntimeException("No square matrix supplied. Can´t compute square root of matrix.");
		}
		
		HashMap<String, double [][]> eigenDec = get_eigen_dec_4_symmetric_matrix(X);
		
		double [][] e_val = eigenDec.get("eigenvalues");
		double [][] e_vec = eigenDec.get("eigenvectors");
		
		for(int i=0; i<n_vars; i++) {
			e_val[i][i] = 1.0/Math.sqrt(e_val[i][i]);
		}
		
		double [][] sq_root_X = MatrixOperations.multiplication(MatrixOperations.multiplication(e_vec, e_val),MatrixOperations.transpose(e_vec));
		
		return sq_root_X;		
	}
	
	
    /**
     * Method for calculating singular value decomposition of matrix X.
     * (JBLAS library is used)
     * @param X double matrix
     * @return HashMap of double arrays with keys U,S,V for the respective matrices U,S,V
     */
	@SuppressWarnings("static-access")
	public static HashMap<String, double [][]> get_fullSVD(double [][] X) {
		
		HashMap<String, double [][]> res_SVD = new HashMap<String, double [][]>();
		
		Singular singular = new Singular();		
		DoubleMatrix X_new = new DoubleMatrix(X);		
		DoubleMatrix[] singularValueDec = singular.fullSVD(X_new);
		
		res_SVD.put("U", singularValueDec[0].toArray2());
		res_SVD.put("S", singularValueDec[1].toArray2());
		res_SVD.put("V", singularValueDec[2].toArray2());
		
		return res_SVD;
	}
	
	
	public static double [][] get_matrix_with_equal_elments(double x, int n_rows, int n_cols) {
		
		double [][] m = new double [n_rows][n_cols];
		for(int i=0; i<n_rows; i++) {
			for(int j=0; j<n_cols; j++) {
				m[i][j] = x;
			}
			
		}
		return m;
	}
	
	
    // test client
    public static void main(String[] args) {
    	    	
    	double [][] A = {{1,4,3},{4,5,6},{3,6,10}};
    	
    	print_matrix(A);
    	
    	print_matrix(get_eigen_dec_4_symmetric_matrix(A).get("eigenvectors"));
    	
    	print_matrix(get_eigen_dec_4_symmetric_matrix(A).get("eigenvalues"));
    	
    }
	
}
