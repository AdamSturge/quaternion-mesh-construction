#include <Eigen/Core>
#include <vector>

// Finds the closest points between a mesh and a point cloud
// Inputs:
//   V  vertex matrix of mesh
//   F  face matrix of mesh
//   Y  Mx3 point bloud
// Outputs:
//   X  Mx3 matrix of points on the mesh that are closest to the corresponding point in Y
//   FL list of faces for each point in X. [3,4] would mean the first sample came from face 3 for example
void closest_points_cloud_to_mesh(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::MatrixXd Y,
    Eigen::MatrixXd& X,
	std::vector<int>& FL);