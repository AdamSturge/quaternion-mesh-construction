#include <new_vertex_positions.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <cotmatrix.h>
#include <omega.h>
#include <normalize_solution.h>

void new_vertex_positions(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::VectorXd lam,
	Eigen::MatrixXd& U)
{
	int n = V.rows();
	
	Eigen::SparseMatrix<double> L;
	cotmatrix(V, F, L);

	Eigen::VectorXd ome(4*n);
	omega(V, F, lam, ome);

	std::cout << "Begin solving for vertex positions" << std::endl;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> solver(L.transpose()*L);
	//solver.setMaxIterations(2000);
	//solver.setTolerance(1e-7);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(L.transpose()*L);
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver(L);
	Eigen::VectorXd u_vec = solver.solve(L.transpose()*ome);

	std::cout << u_vec.topRows(12) << std::endl;

	/*std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error() << std::endl;
	std::cout << "Done solving for vertex positions" << std::endl;
*/
	if (solver.info() != Eigen::Success)// || solver.error() > 0.1)
	{
		std::cout << "SOLVING FOR NEW VERTEX POSITIONS FAILED" << std::endl;
		U = V;
	}
	else 
	{
		for (int i = 0; i < n; ++i)
		{
			U.row(i) = u_vec.block(4 * i + 1, 0, 3, 1).transpose();
		}
	}


}