#include <solve_for_lambda.h>
#include <closest_points_cloud_to_mesh.h>
#include <unordered_map>
#include <closest_plane.h>
#include <quaternion.h>
#include <iostream>
#include <sample_points.h>
#include <vector>
#include <numeric>
#include <functional>

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
	//closest_points_could_to_mesh(V, F, Y, X, FL);

	std::unordered_map<int, std::vector<int>> face_counts;
	for(int i = 0; i < FL.size(); ++i)
	{
		if(face_counts.count(FL[i]) == 0)
		{
			face_counts[FL[i]] = std::vector<int>();
		}
		face_counts[FL[i]].push_back(i);
	}

	std::unordered_map<int, std::vector<double>> q_magnitude;
	std::unordered_map<int, std::vector<double>> q_angle;
	for(int i = 0; i < V.rows(); ++i)
	{
		q_magnitude[i] = std::vector<double>();
		q_angle[i] = std::vector<double>();
		q_magnitude[i].reserve(25);
		q_angle[i].reserve(25);
	}


	Eigen::MatrixXd Mesh_Points;
	Eigen::MatrixXd Cloud_Points;
	Eigen::MatrixXd N_LS(F.rows(),3);
	N_LS.setZero();
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
			N_LS.row(it.first) = n.transpose();

			int count = Mesh_Points.rows();
			for(int i = 0; i < Mesh_Points.rows(); ++i)
			{
				Eigen::Vector3d u = Mesh_Points.row(i);
				Eigen::Vector3d w = Cloud_Points.row(i);
				double u_norm = u.norm();
				double w_norm = w.norm();
				u = u / u_norm;
				w = w / w_norm;
				double angle = acos(u.dot(w));
				int sign = u.cross(w).dot(n) > 0 ? 1 : -1;
				angle = sign*angle;

				double mag = sqrt((w_norm / u_norm));
				//Quaternion q(cos(angle), sin(angle)*n(0), sin(angle)*n(1), sin(angle)*n(2));
				//q = mag*q;

				int I[3] = { F(it.first,0), F(it.first,1), F(it.first,2) };
				for(int j = 0; j < 3; ++j)
				{
					q_magnitude[I[j]].push_back(mag);
					//std::cout << I[j] << "->" << angle << std::endl;
					q_angle[I[j]].push_back(angle);
					//lam(4 * I[j] + 0) += (1.0 / 3.0) * (1.0 / (double)count) * q.re();
					//lam(4 * I[j] + 1) += (1.0 / 3.0) * (1.0 / (double)count) * q.im()(0);
					//lam(4 * I[j] + 2) += (1.0 / 3.0) * (1.0 / (double)count) * q.im()(1);
					//lam(4 * I[j] + 3) += (1.0 / 3.0) * (1.0 / (double)count) * q.im()(2);
				}
			}
		}
	}

	for(int i = 0; i < F.rows(); ++i)
	{
		int I[3] = { F(i,0), F(i,1), F(i,2) };
		Eigen::Vector3d n = N_LS.row(i);
		if (n.isApprox(Eigen::Vector3d::Zero(),0.1))
		{
			continue;
		}
	
		for(int j = 0; j < 3; ++j)
		{
			if(q_magnitude.count(I[j]) == 1 && q_angle.count(I[j]) == 1)
			{
				double mean_mag = std::accumulate(q_magnitude[I[j]].begin(), q_magnitude[I[j]].end(), 0.0) / (double)q_magnitude[I[j]].size();
				double mean_angle = std::accumulate(q_angle[I[j]].begin(), q_angle[I[j]].end(), 0.0) / (double)q_angle[I[j]].size();
				
				Eigen::Vector3d angle_vec(1.0, 1.0, 1.0);
				angle_vec = sin(mean_angle / 2.0)*angle_vec;
				Eigen::Vector3d prod = angle_vec.cwiseProduct(n);

				lam(4 * I[j] + 0) =  mean_mag * cos(mean_angle / 2.0);
				lam(4 * I[j] + 1) =  mean_mag * prod(0);
				lam(4 * I[j] + 2) =  mean_mag * prod(1);
				lam(4 * I[j] + 3) =  mean_mag * prod(2);

				q_magnitude.erase(I[j]);
				q_angle.erase(I[j]);
			
			}
			
		}
		
	}
}