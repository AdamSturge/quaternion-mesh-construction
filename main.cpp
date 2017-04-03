#define GLOG_NO_ABBREVIATED_SEVERITIES // Has a clash with windows.h. Not sure where that's being included. This fixes it?
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#define NO_STRICT
#include <igl/list_to_matrix.h>
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <eigen_fileio.h>
#include <igl/cotmatrix.h>
#include <gradient_ascent.h>
#include <igl/barycentric_coordinates.h>
#include <closest_points_cloud_to_mesh.h>
#include <sample_points.h>
#include <solve_for_lambda.h>
#include <quaternion_matrix.h>
#include <convert_to_stan_var.h>
#include <normalize_solution.h>
#include <thread>
#include <mutex>
#include <solve_for_lambda.h>
#include <new_vertex_positions.h>

int main(int argc, char *argv[])
{
	// predefined colors
	const Eigen::RowVector3d white(1.0, 1.0, 1.0); 
  const Eigen::RowVector3d red(1, 0, 0);
	const Eigen::RowVector3d green(0, 1, 0);
	
  //Load in points + normals from .pwn file
	Eigen::MatrixXd target_points;
	Eigen::MatrixXd target_normals;
  {
    Eigen::MatrixXd D;
    std::vector<std::vector<double> > vD;
    std::string line;
    std::fstream in;
    //in.open(argc>1?argv[1]:"../shared/data/sphere-noisy.pwn");
		//in.open(argc>1 ? argv[1] : "../shared/data/hand.pwn");
		in.open(argc>1 ? argv[1] : "C:\\UoT-Masters\\Geometry-Processing\\geometry-processing-project\\shared\\data\\head.pwn");
		
    while(in)
    {
      std::getline(in, line);
      std::vector<double> row;
      std::stringstream stream_line(line);
      double value;
      while(stream_line >> value) row.push_back(value);
      if(!row.empty()) vD.push_back(row);
    }
    if(!igl::list_to_matrix(vD,D)) return EXIT_FAILURE;
    target_points = D.leftCols(3);
		target_normals = D.rightCols(3);
  }
	int n_TP = target_points.rows();

  // Load seed mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh("../shared/data/sphere.obj", V, F);
	//igl::read_triangle_mesh("../shared/data/24rqomjz3oqo-face/face.obj", V, F);
	//igl::read_triangle_mesh("C:\\Users\\adams\\Desktop\\deformed_mesh_head_overnight.obj", V, F);

	normalize_solution(V, V);

	Eigen::MatrixXd U = V; // matrix used to store new vertex locations

  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::viewer::Viewer viewer;
  std::cout<<R"(
  )";

	// load seed mesh into viewer
	viewer.data.set_mesh(V, F);
	viewer.core.point_size = 10;

	// load target points in viewer
	Eigen::MatrixXd points = target_points;
	viewer.data.set_points(points, white);				

	int num_samples = n_TP; // NOTE: THIS HAS TO BE NUMBER OF TARGET POINTS IF YOU ARE ALERTERNATING SAMPLING TECHNIQUES

	Eigen::MatrixXd C(n_TP + 2 * num_samples, 3);
	C = Eigen::MatrixXd(n_TP + 2 * num_samples, 3);
	C.topRows(2 * num_samples).rowwise() = red; // X color
	C.topRows(num_samples).rowwise() = green; // P color
	C.bottomRows(n_TP).rowwise() = white; // target point color

  #pragma region cotmatrix_computation
	Eigen::SparseMatrix<double> L(V.rows(), V.cols());
	igl::cotmatrix(V, F, L);
	L = -L;

	Eigen::SparseMatrix<double> L_big(4 * L.rows(), 4 * L.cols());
	real_quaternion_sparse_matrix_to_real_sparse_matrix(L, L_big);
	Eigen::SparseMatrix<stan::math::var> L_big_var(4 * L.rows(), 4 * L.cols());
	convert_to_stan_var_sparse(L_big, L_big_var);
	
  Eigen::SimplicialLLT<Eigen::SparseMatrix<stan::math::var>> L_sol;
	L_sol.compute(L_big_var.transpose()*L_big_var);
#pragma endregion

  #pragma region sampling_functions
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> P_vec(4 * num_samples, 1);
	Eigen::SparseMatrix<double> X_b_big = Eigen::SparseMatrix<double>(4 * num_samples, 4 * V.rows());
	Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> N_var(num_samples, 3);
	//std::mutex mutex;
	const auto compute_grad_decent_matricies = [&](Eigen::MatrixXd P, Eigen::MatrixXd X_w, Eigen::MatrixXd N, std::vector<int> FL)
	{
		// Vectorize P
		Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> P_vec_tmp(4 * num_samples, 1);
		for (int i = 0; i < num_samples; ++i)
		{
			P_vec_tmp(4 * i) = stan::math::var(0.0);
			for (int j = 1; j < 4; ++j)
			{
				P_vec_tmp(4 * i + j) = stan::math::var(P(i, j - 1));
			}
		}
		// Create barycentric coordinate matrix X_b from X_w
		Eigen::MatrixXd A = Eigen::MatrixXd(num_samples, 3);
		Eigen::MatrixXd B = Eigen::MatrixXd(num_samples, 3);
		Eigen::MatrixXd D = Eigen::MatrixXd(num_samples, 3);
		for (int i = 0; i < num_samples; ++i)
		{
			int I[3] = { F(FL[i],0), F(FL[i],1), F(FL[i],2) };
			A.row(i) = U.row(I[0]);
			B.row(i) = U.row(I[1]);
			D.row(i) = U.row(I[2]);
		}
		Eigen::MatrixXd X_b = Eigen::MatrixXd(num_samples, 3);
		igl::barycentric_coordinates(X_w, A, B, D, X_b);

		Eigen::SparseMatrix<double> X_b_big_tmp(4 * num_samples, 4 * U.rows());
		std::vector<Eigen::Triplet<double>> triplets;
		triplets.reserve(4 * 3 * U.rows());
		for (int i = 0; i < num_samples; ++i)
		{
			int I[3] = { F(FL[i],0), F(FL[i],1), F(FL[i],2) };

			for (int j = 0; j < 4; ++j)
			{
				triplets.push_back(Eigen::Triplet<double>(4 * i + j, 4 * I[0] + j, X_b(i, 0)));
				triplets.push_back(Eigen::Triplet<double>(4 * i + j, 4 * I[1] + j, X_b(i, 1)));
				triplets.push_back(Eigen::Triplet<double>(4 * i + j, 4 * I[2] + j, X_b(i, 2)));
			}
		}
		X_b_big_tmp.setFromTriplets(triplets.begin(), triplets.end());

		Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> N_var_tmp(num_samples, 3);
		convert_to_stan_var(N, N_var_tmp);

		//mutex.lock();
		P_vec = P_vec_tmp;
		X_b_big = X_b_big_tmp;
		N_var = N_var_tmp;
		//mutex.unlock();
	
	};
	
  //bool sampling_running = false;
	const auto mesh_to_cloud_sample = [&]()
	{
		//sampling_running = true;
		std::cout << "sampling: mesh to cloud" << std::endl;
		Eigen::MatrixXd X_w = Eigen::MatrixXd(num_samples, 3);
		Eigen::MatrixXd P = Eigen::MatrixXd(num_samples, 3);
		Eigen::MatrixXd N = Eigen::MatrixXd(num_samples, 3);
		std::vector<int> FL;

		sample_points(num_samples, U, F, target_points, target_normals, X_w, P, N, FL);
		compute_grad_decent_matricies(P, X_w, N, FL);
		
		//sampling_running = false;
		std::cout << "done sampling: mesh to cloud" << std::endl;
	};

	const auto cloud_to_mesh_sample = [&]()
	{
		//sampling_running = true;
		std::cout << "sampling: cloud to mesh" << std::endl;
		Eigen::MatrixXd X_w = Eigen::MatrixXd(num_samples, 3);
		Eigen::MatrixXd P = Eigen::MatrixXd(num_samples, 3);
		Eigen::MatrixXd N = Eigen::MatrixXd(num_samples, 3);
		std::vector<int> FL;

		P = target_points;
		N = target_normals;
		closest_points_cloud_to_mesh(U, F, P, X_w, FL);
		
		compute_grad_decent_matricies(P, X_w, N, FL);
		//sampling_running = false;
		std::cout << "done sampling: cloud to mesh" << std::endl;
	};
  #pragma endregion

	cloud_to_mesh_sample();

	#pragma region variables_for_grad_decent
	Eigen::VectorXd lam(4 * V.rows());
	double step_size = -1.0;
	double step_decay = 0.9;
	bool points_visible = true;
	bool exit = true;
	Eigen::VectorXd u(4 * V.rows()); //vectorized new vertex positions (as quaternions)
	int sampling_trigger = 1;
	bool cloud_to_mesh = false;
	std::vector<double> fx_vals;
	double fx_min_old = INT_MAX;
	double fx_min_new = 0;
	
	lam.setConstant(0.00);
	for (int i = 0; i < 4 * V.rows(); ++i)
	{
		if (i % 4 == 0)
		{
			lam(i) = 1.0;
		}
	}
 #pragma endregion  

	viewer.callback_pre_draw = [&](igl::viewer::Viewer &)->bool
	{
		if (viewer.core.is_animating)
		{
			if (exit == false)
			{
				sampling_trigger++;
				if((sampling_trigger % 10) == 0)
				{
					exit = true;
					fx_min_new = *std::min(fx_vals.begin(),fx_vals.end());
					if(fx_min_new - fx_min_old > 0 || abs(fx_min_new - fx_min_old) < 0.0001)
					{
						cloud_to_mesh = !cloud_to_mesh;
						fx_min_old = INT_MAX;
					}else
					{
						fx_min_old = fx_min_new;
					}
					
					fx_vals.clear();
					if(cloud_to_mesh)
					{
						cloud_to_mesh_sample();
						/*std::cout << "starting thread" << std::endl;
						std::thread t(cloud_to_mesh_sample);
						t.detach();*/
					}
					else 
					{
						mesh_to_cloud_sample();
						/*std::cout << "starting thread" << std::endl;
						std::thread t(mesh_to_cloud_sample);
						t.detach();*/
					}
				}
				if ((sampling_trigger % 100) == 0)
				{
					std::cout << "reducing step size from " << step_size << " to " << step_decay*step_size << std::endl;
					step_size = step_decay*step_size;
					sampling_trigger = 0;
				}

				Energy E(V, F, L_big_var, L_sol, X_b_big, P_vec, N_var);
				double fx;
				Eigen::Matrix<double, Eigen::Dynamic, 1> grad_fx;
				stan::math::gradient(E, lam, fx, grad_fx);
				fx_vals.push_back(fx);

				grad_fx.normalize();
				
				lam = lam + step_size*grad_fx;
				std::cout << "energy: " << fx << std::endl;
				E.get_new_vertices(u);

				Eigen::VectorXd u_block(4);
				for (int i = 0; i < U.rows(); ++i)
				{
					u_block = u.block(4 * i, 0, 4, 1);
					U.row(i) = u_block.bottomRows(3).transpose();
				}
				
				viewer.data.set_mesh(U, F);
			}
		}
		return false;
	};
	
	viewer.callback_key_pressed = [&](igl::viewer::Viewer&, unsigned int key,int)
  {
    switch(key)
    {
		case 'C':
		case 'c':
			cloud_to_mesh_sample();
			return true;
		case 'H':
		case 'h': 
			exit = !exit; 
			return true;
		case 'I':
    case 'i':
			mesh_to_cloud_sample();
			return true;
		case 'P':
		case 'p':
			if(points_visible)
			{
				Eigen::MatrixXd empty;
				viewer.data.set_points(empty, white);
			}
			else
			{
				Eigen::MatrixXd C2 = C.bottomRows(points.rows());
				viewer.data.set_points(points, C2);
			}
			points_visible = !points_visible;
			return true;
		case 'V':
		case 'v':
			std::cout << "saving" << std::endl;
			viewer.save_mesh_to_file("C:\\Users\\adams\\Desktop\\deformed_mesh.obj");
			std::cout << "saving complete" << std::endl;
			return true;
    }
    return false;
  };
	viewer.core.is_animating = true;
  viewer.launch();
	

  return EXIT_SUCCESS;
}

