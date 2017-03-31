#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stan\math.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

struct normal_ll {
	const Eigen::MatrixXd V_;
	const Eigen::MatrixXi F_;
	boost::shared_ptr<Eigen::SparseMatrix<double>> L_pinv_ptr_;
	boost::shared_ptr<Eigen::SparseMatrix<double>> X_ptr_;
	boost::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, 1>> P_ptr_;

	normal_ll(
		const Eigen::MatrixXd V, 
		const Eigen::MatrixXi F, 
		boost::shared_ptr<Eigen::SparseMatrix<double>> L_pinv_ptr,
		boost::shared_ptr<Eigen::SparseMatrix<double>> X_ptr,
		boost::shared_ptr<Eigen::Matrix<double,Eigen::Dynamic,1>> P_ptr) : V_(V), F_(F)
	{
		L_pinv_ptr_ = L_pinv_ptr;
		X_ptr_ = X_ptr;
		P_ptr_ = P_ptr;
	}

	template <typename T>
	void quaternion_block(const T a, const T b, const T c, const T d, Eigen::Matrix<T, 4, 4>& Q, Eigen::Matrix<T, 4, 4>& Q_bar) const
	{
		Q << a, -b, -c, -d,
				 b,  a, -d,  c,
			   c,  d,  a, -b,
				 d, -c,  b,  a;

		Q_bar <<  a,  b,  c,  d,
						 -b,  a,  d, -c,
						 -c, -d,  a,  b,
						 -d,  c, -b,  a;
	}

	template <typename T>
	T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& lambda) const {

		int n_V = V_.rows();
		int n_F = F_.rows();
		Eigen::Matrix<T,Eigen::Dynamic,1> omega(4 * n_V, 1);
		omega.setConstant(T(0.0));

		for (int i = 0; i < n_F; ++i) 
		{
			// get indices of the vertices of this face
			int v[3] = { F_(i,0), F_(i,1), F_(i,2) };

			for(int j = 0; j < 3; ++j)
			{
				Eigen::RowVector3d f0 = V_.row(v[(j + 0) % 3]);
				Eigen::RowVector3d f1 = V_.row(v[(j + 1) % 3]);
				Eigen::RowVector3d f2 = V_.row(v[(j + 2) % 3]);

				// determine orientation of this edge
				int a = v[(j + 1) % 3];
				int b = v[(j + 2) % 3];
				if (a > b)
				{
					std::swap(a, b);
				}

				Eigen::Matrix<T, 4, 4> l_a;
				Eigen::Matrix<T, 4, 4> l_a_bar;
				quaternion_block(lambda(4 * a + 0), lambda(4 * a + 1), lambda(4 * a + 2), lambda(4 * a + 3), l_a, l_a_bar);
				Eigen::Matrix<T, 4, 4> l_b;
				Eigen::Matrix<T, 4, 4> l_b_bar;
				quaternion_block(lambda(4 * b + 0), lambda(4 * b + 1), lambda(4 * b + 2), lambda(4 * b + 3), l_b, l_b_bar);

				Eigen::Matrix<double, 4, 4> e_mat;
				Eigen::Matrix<double, 4, 4> e_mat_bar;
				const Eigen::RowVector3d e = V_.row(b) - V_.row(a);
				quaternion_block(0.0, e(0), e(1), e(2), e_mat, e_mat_bar);

				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> e_tilde =
					(1.0 / 3.0) * stan::math::multiply(stan::math::multiply(l_a_bar, e_mat), l_a) +
					(1.0 / 6.0) * stan::math::multiply(stan::math::multiply(l_a_bar, e_mat), l_b) +
					(1.0 / 6.0) * stan::math::multiply(stan::math::multiply(l_b_bar, e_mat), l_a) +
					(1.0 / 3.0) * stan::math::multiply(stan::math::multiply(l_b_bar, e_mat), l_b) ;

				Eigen::Vector3d u1 = (f1 - f0);
				Eigen::Vector3d u2 = (f2 - f0);
				double cotAlpha = u1.dot(u2) / u1.cross(u2).norm();

				e_tilde = (1.0/2.0)*cotAlpha*e_tilde;

				omega(4 * a + 0) -= e_tilde(0,0);
				omega(4 * a + 1) -= e_tilde(1,0);
				omega(4 * a + 2) -= e_tilde(2,0);
				omega(4 * a + 3) -= e_tilde(3,0);

				omega(4 * b + 0) += e_tilde(0,0);
				omega(4 * b + 1) += e_tilde(1,0);
				omega(4 * b + 2) += e_tilde(2,0);
				omega(4 * b + 3) += e_tilde(3,0);
			}
		}

		omega = omega - Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(omega.rows(), omega.cols(), omega.mean());

		// HUGE RAM SPIKE INCOMING
		Eigen::Matrix<T, Eigen::Dynamic, 1> f_tilde = (*L_pinv_ptr_)*omega;
		Eigen::Matrix<T, Eigen::Dynamic, 1> x_tilde = (*X_ptr_)*f_tilde;

		return stan::math::dot_product(x_tilde,*P_ptr_);
		//return 1;
	}
};

void autograd_test(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::MatrixXd Y,
  Eigen::VectorXd& lam);
