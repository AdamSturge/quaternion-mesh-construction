#include <stan/math.hpp>
#include <Eigen/Core>

void convert_to_stan_var(
	 const Eigen::MatrixXd X, 
	 Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>& X_var);

void convert_to_stan_var_sparse(
	const Eigen::SparseMatrix<double> X,
	Eigen::SparseMatrix<stan::math::var>& X_var);