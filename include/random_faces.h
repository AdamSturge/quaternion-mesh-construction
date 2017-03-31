#pragma once
#include <Eigen/Core>
#include <vector>

/**
Selects random faces from a tiangle mesh
Inputs :
n  number of faces to select
V  verices of mesh
F  faces of mesh
Outputs :
LF List of indices into F
**/

void random_faces(
	const int n,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	std::vector<int> & FL);