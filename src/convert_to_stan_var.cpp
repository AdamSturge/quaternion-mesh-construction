#include <convert_to_stan_var.h>

void convert_to_stan_var(
	const Eigen::MatrixXd X,
	Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>& X_var)
{
	for(int i = 0; i < X.rows(); ++i)
	{
		for(int j = 0; j < X.cols(); ++j)
		{
			X_var(i, j) = stan::math::var(X(i, j));
		}
	}
}

void convert_to_stan_var_sparse(
	const Eigen::SparseMatrix<double> X,
	Eigen::SparseMatrix<stan::math::var>& X_var)
{
	std::vector<Eigen::Triplet<stan::math::var>> triplets;
	for(int k = 0; k < X.outerSize(); ++k)
	{
		for(Eigen::SparseMatrix<double>::InnerIterator it(X,k); it; ++it)
		{
			triplets.push_back(Eigen::Triplet<stan::math::var>(it.row(), it.col(), stan::math::var(it.value())));
		}
	}
	X_var.setFromTriplets(triplets.begin(), triplets.end());
}