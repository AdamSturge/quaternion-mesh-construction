#include <Eigen/Core>

//  Normalized a mesh so that its vertices lie centered at the origin
// Inputs:
//      V   Nx3 original mesh
// Outputs:
//      U  Nx3 normalized mesh
void normalize_solution(
	const Eigen::MatrixXd V,
	Eigen::MatrixXd& U);