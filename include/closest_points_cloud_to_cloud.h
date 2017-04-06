#include <Eigen\Core>
#include <vector>

// Finds the closest points between two point clouds as well 
// as the associated normals for the second point cloud
// Inputs:
//   X  Kx3 point cloud
//   Y  Mx3 point bloud
//   N_Y  normals to points in Y
// Outputs:
//   P Kx3 matrix of points in Y that are closest to the corresponding points in X
//   N Kx3 matrix of normals to points in Y that are closest to the corresponding points in X
void closest_points_cloud_to_cloud(
	const Eigen::MatrixXd X,
	const Eigen::MatrixXd Y,
	const Eigen::MatrixXd N_Y,
	Eigen::MatrixXd& P, 
	Eigen::MatrixXd& N);