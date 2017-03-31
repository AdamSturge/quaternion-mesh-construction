#include <Eigen/Core>

void solve_for_lambda(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::MatrixXd Y,
	const Eigen::MatrixXd N_Y,
	Eigen::VectorXd& lam);