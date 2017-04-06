#include <cotmatrix.h>
#include <igl\cotmatrix.h>
#include <quaternion.h>
#include <quaternion_matrix.h>

void cotmatrix(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	Eigen::SparseMatrix<double>& L) 
{
	int n_V = V.rows();
	L = Eigen::SparseMatrix<double>(4 * n_V, 4 * n_V);

	Eigen::SparseMatrix<double> L_small(n_V,n_V);
	igl::cotmatrix(V, F, L_small);
	L_small = -L_small;

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(4 * 4 * n_V);
	double val;
	double Q[4][4];
	Quaternion q;
	int i;
	int j;
	for (int k = 0; k < L_small.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(L_small, k); it; ++it)
		{
			val = it.value();
			i = it.row();
			j = it.col();
			q = Quaternion(val);
			q.toMatrix(Q);
			insert_quaternion_block(Q, i, j, triplets);
		}
	}
	L.setFromTriplets(triplets.begin(), triplets.end());
}