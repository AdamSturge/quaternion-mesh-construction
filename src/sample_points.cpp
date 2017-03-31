#include <sample_points.h>
#include <closest_points_cloud_to_cloud.h>
#include <random_points_on_mesh.h>

void sample_points(
	const int k,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & Y,
	const Eigen::MatrixXd & N_Y,
	Eigen::MatrixXd & X,
	Eigen::MatrixXd & P,
	Eigen::MatrixXd & N,
	std::vector<int>& FL)
{
	random_points_on_mesh(k, V, F, X, FL);

	closest_points_cloud_to_cloud(X, Y, N_Y, P, N);
}