#include <stan/math.hpp>
#include <Eigen/Core>

// Converts a matrix of doubles to a matrix of stan::math::var
// Inputs:
//   X  NxM matrix of doubles
// Outputs:
//   X_var  NxM matrix of stan::math:var whose values correspond to the values in X
void convert_to_stan_var(
	 const Eigen::MatrixXd X, 
	 Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>& X_var);

 // Converts a sparse matrix of doubles to a matrix of stan::math::var
// Inputs:
//   X  NxM sparse matrix of doubles
// Outputs:
//   X_var  sparse NxM matrix of stan::math:var whose values correspond to the values in X
void convert_to_stan_var_sparse(
	const Eigen::SparseMatrix<double> X,
	Eigen::SparseMatrix<stan::math::var>& X_var);