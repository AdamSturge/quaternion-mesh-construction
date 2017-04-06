#include <Eigen\Core>
#include "random_points_on_mesh.h"

// Randomly samples points from a mesh and computes their closest point in a point cloud
// Inputs:    
//      k  number of samples to draw
//      V  vertex matrix for mesh
//      F  face matrix for mesh
//      Y  point cloud
//      N_Y  point cloud normals
// Outputs:
//      X  kx3 matrix of points sampled from the mesh
//      P  kx3 matrix of closest points in Y to the points in X
//      N  kx3 matrix of normals to closest points in Y to the points in X
//      FL list of faces for each point in X. [3,4] would mean the first sample came from face 3 for example
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