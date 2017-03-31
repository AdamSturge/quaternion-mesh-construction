#include <Eigen/Core>
#include <vector>

void closest_points_cloud_to_mesh(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::MatrixXd Y,
  Eigen::MatrixXd& X,
	std::vector<int>& FL);