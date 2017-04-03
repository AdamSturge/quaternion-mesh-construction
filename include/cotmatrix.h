#include <Eigen\Core>
#include <Eigen\SparseCore>

void cotmatrix(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	Eigen::SparseMatrix<double>& L);