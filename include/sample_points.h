#include <Eigen\Core>
#include "random_points_on_mesh.h"

void sample_points(
	const int k,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & Y,
	const Eigen::MatrixXd & N_Y,
	Eigen::MatrixXd & X,
	Eigen::MatrixXd & P,
	Eigen::MatrixXd & N,
	std::vector<int>& FL);