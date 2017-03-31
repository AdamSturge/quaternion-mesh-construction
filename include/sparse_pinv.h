#include <Eigen\Sparse>
#include <Eigen\SVD>

void sparse_pinv(
	const Eigen::SparseMatrix<double> X,
	Eigen::SparseMatrix<double>& X_pinv);