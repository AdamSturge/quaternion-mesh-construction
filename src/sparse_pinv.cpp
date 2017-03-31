#include <sparse_pinv.h>


void sparse_pinv(
	Eigen::SparseMatrix<double> X,
	Eigen::SparseMatrix<double>& X_pinv)
{
	X.makeCompressed();

	Eigen::SparseMatrix<double> I(X.rows(), X.cols());
	I.setIdentity();

	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int>> QR(X);
	X_pinv = QR.solve(I);
}