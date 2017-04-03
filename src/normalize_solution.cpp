#include <normalize_solution.h>
#include <Eigen/Core>

void normalize_solution(
	const Eigen::MatrixXd V,
	Eigen::MatrixXd& U)
{
	U = Eigen::MatrixXd(V.rows(), V.cols());
	
	// center vertices around the origin
	U = V - V.mean()*Eigen::MatrixXd::Ones(V.rows(),V.cols());

	// rescale vertices to be at most norm 1
	double r = V.rowwise().norm().maxCoeff();
	U = U / r;

}