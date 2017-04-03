#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stan\math.hpp>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

struct Energy {
	const Eigen::MatrixXd& V_;
	const Eigen::MatrixXi& F_;
	const Eigen::SparseMatrix<stan::math::var>& L_;
	const Eigen::SimplicialLLT<Eigen::SparseMatrix<stan::math::var>>& L_sol_;
	const Eigen::SparseMatrix<double>& Baryl_coords_;
	const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& P_;
	const Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>& N_;
	mutable Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> f_tilde_;

	Energy(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<stan::math::var>& L,
		const Eigen::SimplicialLLT<Eigen::SparseMatrix<stan::math::var>>& L_sol,
		const Eigen::SparseMatrix<double>& BC,
		const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& P,
		const Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>& N) : V_(V), F_(F), L_(L), L_sol_(L_sol), Baryl_coords_(BC), P_(P), N_(N)
	{
		f_tilde_ = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>();
	}

	void get_new_vertices(Eigen::VectorXd& u)
	{
		int n_f = f_tilde_.rows();
		u = Eigen::VectorXd(n_f);
		for(int i = 0; i < n_f; ++i)
		{
			u(i) = f_tilde_(i).val();
		}
	}

	template <typename T>
	void quaternion_block(const T a, const T b, const T c, const T d, Eigen::Matrix<T, 4, 4>& Q, Eigen::Matrix<T, 4, 4>& Q_bar) const
	{
		Q << a, -b, -c, -d,
			   b,  a, -d,  c,
			   c,  d,  a, -b,
			   d, -c,  b,  a;

		Q_bar << a,  b,  c,  d,
			      -b,  a,  d, -c,
			      -c, -d,  a,  b,
			      -d,  c, -b,  a;
	}

	stan::math::var operator()(const Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>& lambda) const {

		int n_V = V_.rows();
		int n_F = F_.rows();
		Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> omega(4 * n_V, 1);
		omega.setConstant(stan::math::var(0.0));
		for (int i = 0; i < n_F; ++i)
		{
			// get indices of the vertices of this face
			int v[3] = { F_(i,0), F_(i,1), F_(i,2) };

			for (int j = 0; j < 3; ++j)
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

				Eigen::Matrix<stan::math::var, 4, 4> l_a;
				Eigen::Matrix<stan::math::var, 4, 4> l_a_bar;
				quaternion_block(lambda(4 * a + 0), lambda(4 * a + 1), lambda(4 * a + 2), lambda(4 * a + 3), l_a, l_a_bar);
				Eigen::Matrix<stan::math::var, 4, 4> l_b;
				Eigen::Matrix<stan::math::var, 4, 4> l_b_bar;
				quaternion_block(lambda(4 * b + 0), lambda(4 * b + 1), lambda(4 * b + 2), lambda(4 * b + 3), l_b, l_b_bar);

				Eigen::Matrix<double, 4, 4> e_mat;
				Eigen::Matrix<double, 4, 4> e_mat_bar;
				const Eigen::RowVector3d e = V_.row(b) - V_.row(a);
				quaternion_block(0.0, e(0), e(1), e(2), e_mat, e_mat_bar);

				Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> e_tilde =
					(1.0 / 3.0) * stan::math::multiply(stan::math::multiply(l_a_bar, e_mat), l_a) +
					(1.0 / 6.0) * stan::math::multiply(stan::math::multiply(l_a_bar, e_mat), l_b) +
					(1.0 / 6.0) * stan::math::multiply(stan::math::multiply(l_b_bar, e_mat), l_a) +
					(1.0 / 3.0) * stan::math::multiply(stan::math::multiply(l_b_bar, e_mat), l_b);

				Eigen::Vector3d u1 = (f1 - f0);
				Eigen::Vector3d u2 = (f2 - f0);
				double cotAlpha = u1.dot(u2) / u1.cross(u2).norm();

				e_tilde = (1.0 / 2.0)*cotAlpha*e_tilde;

				omega(4 * a + 0) -= e_tilde(0, 0);
				omega(4 * a + 1) -= e_tilde(1, 0);
				omega(4 * a + 2) -= e_tilde(2, 0);
				omega(4 * a + 3) -= e_tilde(3, 0);

				omega(4 * b + 0) += e_tilde(0, 0);
				omega(4 * b + 1) += e_tilde(1, 0);
				omega(4 * b + 2) += e_tilde(2, 0);
				omega(4 * b + 3) += e_tilde(3, 0);
			}
		}
		omega = omega - Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>::Constant(omega.rows(), omega.cols(), omega.mean());

		Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> f_tilde = L_sol_.solve(L_.transpose()*omega);
		f_tilde_ = f_tilde;
		Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x_tilde = Baryl_coords_*f_tilde;

		Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> diff = P_ - x_tilde;
		int num_samples = 0.25*Baryl_coords_.rows();
		Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> plane_dist(num_samples, 1);
		
		for(int i = 0; i < num_samples; ++i)
		{
			Eigen::Matrix<stan::math::var, 3, 1> u = diff.block(4 * i + 1, 0, 3, 1);
			Eigen::Matrix<stan::math::var, 1, 3> v = N_.row(i);
			plane_dist(i) = v*u;
		}
		return (stan::math::dot_product(plane_dist, plane_dist))/ num_samples;
	}
};

