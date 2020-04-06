package Mathematics;

import org.jblas.Decompose;
import org.jblas.Decompose.LUDecomposition;
import org.jblas.DoubleMatrix;

public class LU {

	double [][] L;
	double [][] U;
	double [][] P;
	
	public LU(double [][] X) {
		
		DoubleMatrix X_new = new DoubleMatrix(X);		
		LUDecomposition<DoubleMatrix> lu_dec = Decompose.lu(X_new);
		L = lu_dec.l.toArray2();
		U = lu_dec.u.toArray2();
		P = lu_dec.p.toArray2();
	}
	
	
	public double [][] get_permutation_matrix() {
		if(P == null) {
			throw new RuntimeException("No LU-decomposition done yet. Found new permutation matrix.");
		}
		return P;
	}

	
	public double [][] get_lower_triangular() {
		if(L == null) {
			throw new RuntimeException("No LU-decomposition done yet. Found new lower triangular matrix.");
		}
		return L;
	}
	
	
	public double [][] get_upper_triangular() {
		if(U == null) {
			throw new RuntimeException("No LU-decomposition done yet. Found new upper triangular matrix.");
		}
		return U;
	}
	
}
