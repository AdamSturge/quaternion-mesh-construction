#include <closest_points_cloud_to_mesh.h>
#include <point_triangle_distance.h>

void closest_points_cloud_to_mesh(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::MatrixXd Y, //points
	Eigen::MatrixXd& X,
	std::vector<int>& FL)
{
	for (int i = 0; i < Y.rows(); ++i)
	{
		Eigen::RowVector3d y = Y.row(i);
		double dMin = INT_MAX;
		Eigen::RowVector3d xMin;
		int jMin = 0;
		for (int j = 0; j < F.rows(); ++j)
		{
			double d = 0.0f;
			Eigen::RowVector3d x;

			Eigen::VectorXi vIndices = F.row(j);
			Eigen::RowVector3d a = V.row(vIndices(0));
			Eigen::RowVector3d b = V.row(vIndices(1));
			Eigen::RowVector3d c = V.row(vIndices(2));

			point_triangle_distance(y, a, b, c, d, x);
			if (d < dMin) {
				dMin = d;
				xMin = x;
				jMin = j;
			}
			if (dMin == 0) {
				break; // break out of inner loop if x is already on p
			}
		}
		X.row(i) = xMin;
		FL.push_back(jMin);
	}
}