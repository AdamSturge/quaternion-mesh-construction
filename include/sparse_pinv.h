#include <Eigen\Sparse>
#include <Eigen\SVD>

// Computes a pseudo-inverse of the provided matrix using QR factorization
// Inputs:    
//      X  matrix
// Outputs:
//      X_pinv  pseudo-inverse of X
void sparse_pinv(
	const Eigen::SparseMatrix<double> X,
	Eigen::SparseMatrix<double>& X_pinv);