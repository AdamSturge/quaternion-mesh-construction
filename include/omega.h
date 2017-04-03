#include <Eigen\Core>

void omega(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::VectorXd lambda,
	Eigen::VectorXd& omega
);