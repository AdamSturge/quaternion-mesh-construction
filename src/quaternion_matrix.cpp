#include <Eigen\Dense>
#include <quaternion_matrix.h>
#include <vector>

void mesh_vertices_to_quaternions(const Eigen::MatrixXd V, VectorXq& U) 
{
	int n = V.rows();
	for (int i = 0; i < n; ++i) 
	{
		U(i) = Quaternion(0.0, V.row(i));
	}
}

void quaternions_to_mesh_vertices(const VectorXq U, Eigen::MatrixXd& V)
{
	int n = U.rows();
	for (int i = 0; i < n; ++i) 
	{
		V.row(i) = U(i, 0).im();
	}
}

void quaternion_matrix_to_real_matrix(const MatrixXq QM, Eigen::MatrixXd& RM) 
{
	double Q[4][4];

	int n_rows = QM.rows();
	int n_cols = QM.rows();
	RM = Eigen::MatrixXd(4 * n_rows * n_cols, 4 * n_rows * n_cols);
	Eigen::MatrixXd B;
	for (int i = 0; i < n_rows; ++i) 
	{
		for (int j = 0; j < n_cols; ++j) 
		{
			QM(i, j).toMatrix(Q);
			insert_quaternion_block(Q, i, j, RM);
		}
	}
}

void sparse_quaternion_matrix_to_real_matrix(const SparseMatrixXq QM, Eigen::SparseMatrix<double>& R)
{
	R = Eigen::SparseMatrix<double>(4 * QM.rows(), 4 * QM.cols());
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(16 * QM.nonZeros());
	for (int k = 0; k < QM.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<Quaternion>::InnerIterator it(QM, k); it; ++it)
		{
			int i = it.row();
			int j = it.col();
			Quaternion q = it.value();
			double Q[4][4];
			q.toMatrix(Q);
			insert_quaternion_block(Q, i, j, triplets);
		}
	}
	R.setFromTriplets(triplets.begin(), triplets.end());

}

void real_quaternion_matrix_to_real_sparse_matrix(const Eigen::MatrixXd QM, Eigen::SparseMatrix<double>& R)
{
	int n = QM.rows();
	int m = QM.cols();
	R = Eigen::SparseMatrix<double>(4 * n, 4 * m);
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(16 * n * m);
	for(int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j) 
		{
			for(int k = 0; k < 4; ++k)
			{
				triplets.push_back(Eigen::Triplet<double>(4 * i + k, 4 * j + k, QM(i, j)));
			}
			
		}
	
	}
	R.setFromTriplets(triplets.begin(), triplets.end());
}

void real_quaternion_sparse_matrix_to_real_sparse_matrix(const Eigen::SparseMatrix<double> QM, Eigen::SparseMatrix<double>& R)
{
	int n = QM.rows();
	int m = QM.cols();
	R = Eigen::SparseMatrix<double>(4 * n, 4 * m);
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(16 * n * m);
	for(int k = 0; k < QM.outerSize(); ++k)
	{
		for(Eigen::SparseMatrix<double>::InnerIterator it(QM,k) ; it; ++it)
		{
			for(int l = 0; l < 4; ++l)
			{
				triplets.push_back(Eigen::Triplet<double>(4 * it.row() + l, 4 * it.col() + l, it.value()));
			}
		}
	}
	R.setFromTriplets(triplets.begin(), triplets.end());
}

void insert_quaternion_block(double Q[4][4], int i, int j, Eigen::MatrixXd& RM)
{
	for (int l = 0; l < 4; ++l)
	{
		for (int k = 0; k < 4; ++k)
		{
			RM(4 * i + l, 4 * j + k) = Q[l][k];
		}
	}
}

void insert_quaternion_block(double Q[4][4], int i , int j, std::vector<Eigen::Triplet<double>>& triplets)
{
	for (int l = 0; l < 4; ++l)
	{
		for (int k = 0; k < 4; ++k)
		{
			if (Q[l][k] != 0.0) 
			{
				triplets.push_back(Eigen::Triplet<double>(4 * i + l, 4 * j + k, Q[l][k]));
			}	
		}
	}
}
