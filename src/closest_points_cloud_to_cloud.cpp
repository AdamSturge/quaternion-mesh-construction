#include <closest_points_cloud_to_cloud.h>

void closest_points_cloud_to_cloud(
	const Eigen::MatrixXd X,
	const Eigen::MatrixXd Y,
	const Eigen::MatrixXd N_Y,
	Eigen::MatrixXd& P,
	Eigen::MatrixXd& N)
{
	int n_X = X.rows();
	int n_Y = Y.rows();

	P = Eigen::MatrixXd(n_X, 3);
	N = Eigen::MatrixXd(n_X, 3);

	Eigen::Vector3d x_point;
	Eigen::Vector3d y_point;
	Eigen::Vector3d y_min_point;
	Eigen::Vector3d n_min_point;
	double min;
	double dist;
	for(int i = 0; i < n_X; ++i)
	{
		x_point = X.row(i);
		y_min_point = Eigen::Vector3d();
		n_min_point = Eigen::Vector3d();
		min = INT32_MAX;
		dist = 0.0;
		for(int j = 0; j < n_Y; ++j)
		{
			y_point = Y.row(j);
			dist = (x_point - y_point).squaredNorm();
			if(dist < min)
			{
				y_min_point = y_point;
				n_min_point = N_Y.row(j);
				min = dist;
			}
		}
		P.row(i) = y_min_point;
		N.row(i) = n_min_point;
	}
}