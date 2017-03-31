#include <Eigen\Core>
#include <vector>

void closest_points_cloud_to_cloud(
	const Eigen::MatrixXd X,
	const Eigen::MatrixXd Y,
	const Eigen::MatrixXd N_Y,
	Eigen::MatrixXd& P, 
	Eigen::MatrixXd& N);