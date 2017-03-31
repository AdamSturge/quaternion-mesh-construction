#include <solve_for_lambda.h>
#include <closest_points_cloud_to_mesh.h>
#include <unordered_map>
#include <closest_plane.h>
#include <quaternion.h>
#include <iostream>
#include <sample_points.h>

void solve_for_lambda(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::MatrixXd Y,
	const Eigen::MatrixXd N_Y,
	Eigen::VectorXd& lam)
{
	int sample_number = Y.rows();
	Eigen::MatrixXd X(sample_number,3);
	Eigen::MatrixXd P(sample_number, 3);
	Eigen::MatrixXd N(sample_number, 3);
	std::vector<int> FL;
	sample_points(sample_number, V, F, Y, N_Y, X, P, N, FL);
	//closest_points_mesh_to_cloud(V, F, Y, X, FL);

	std::unordered_map<int, std::vector<int>> face_counts;
	for(int i = 0; i < FL.size(); ++i)
	{
		if(face_counts.count(FL[i]) == 0)
		{
			face_counts[FL[i]] = std::vector<int>();
		}
		face_counts[FL[i]].push_back(i);
	}

	Eigen::MatrixXd Mesh_Points;
	Eigen::MatrixXd Cloud_Points;
	for(auto it : face_counts)
	{
		if(it.second.size() >= 2)
		{
			Mesh_Points = Eigen::MatrixXd(it.second.size(), 3); 
			Cloud_Points = Eigen::MatrixXd(it.second.size(), 3);
			Eigen::Vector3d n;
			for(int i = 0; i < it.second.size(); ++i)
			{
				Mesh_Points.row(i) = X.row(it.second[i]);
				Cloud_Points.row(i) = P.row(it.second[i]);
			}
			Eigen::MatrixXd Data = Eigen::MatrixXd(Mesh_Points.rows() + Cloud_Points.rows(), 3);
			Data << Mesh_Points,Cloud_Points;

			closest_plane_normal(Data.transpose(), n);

			int count = Mesh_Points.rows();
			for(int i = 0; i < Mesh_Points.rows(); ++i)
			{
				Eigen::Vector3d u = Mesh_Points.row(i);
				Eigen::Vector3d w = Cloud_Points.row(i);
				double u_norm = u.norm();
				double w_norm = w.norm();
				u = u / u_norm;
				w = w / w_norm;
				double angle = acos(u.dot(w)) / 2.0;
				int sign = u.cross(w).dot(n) > 0 ? 1 : -1;
				angle = sign*angle;
				Quaternion q(cos(angle), sin(angle)*n(0), sin(angle)*n(1), sin(angle)*n(2));
				q = sqrt((w_norm / u_norm))*q;

				int I[3] = { F(it.first,0), F(it.first,1), F(it.first,2) };
				for(int j = 0; j < 3; ++j)
				{
					lam(4 * I[j] + 0) += (1.0 / 3.0) * (1.0 / (double)count) * q.re();
					lam(4 * I[j] + 1) += (1.0 / 3.0) * (1.0 / (double)count) * q.im()(0);
					lam(4 * I[j] + 2) += (1.0 / 3.0) * (1.0 / (double)count) * q.im()(1);
					lam(4 * I[j] + 3) += (1.0 / 3.0) * (1.0 / (double)count) * q.im()(2);
				}
			}
		}
	}
}