#include <Eigen/Core>
#include <vector>
// RANDOM_POINTS_ON_MESH Randomly sample a mesh (V,F) n times.
//
// Inputs:
//   n  number of samples
//   V  #V by dim list of mesh vertex positions
//   F  #F by 3 list of mesh triangle indices
// Outputs:
//   X  n by 3 list of random points on (V,F)
//
void random_points_on_mesh(
  const int n,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & X,
	std::vector<int>& FL);

