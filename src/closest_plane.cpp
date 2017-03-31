#include <closest_plane.h>
#include <Eigen/SVD>
#include <iostream>

void closest_plane_normal(
	Eigen::MatrixXd X,
	Eigen::Vector3d& n)
{
	
	Eigen::VectorXd centroid = X.rowwise().sum()/X.cols();

	X = X.colwise() - centroid;
	Eigen::JacobiSVD<Eigen::MatrixXd> sol(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
	sol.computeU();
	Eigen::MatrixXd U = sol.matrixU();
	//std::cout << U << std::endl << std::endl;
	//std::cout << V << std::endl << std::endl;
	n = U.col(2);
	n.normalize();

}