#include <Eigen/Core>

void new_vertex_positions(const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::VectorXd lambda,
	Eigen::MatrixXd& U);
