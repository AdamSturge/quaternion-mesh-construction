#include <random_faces.h>
#include <vector>
#include "igl/doublearea.h"
#include "igl/cumsum.h"
#include <random>

void random_faces(
	const int n,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	std::vector<int> & FL) 
{
	// Define standard uniform R.V sampler
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0, 1);

	// Compute cumsum
	Eigen::VectorXd dblA;
	Eigen::VectorXd C;
	igl::doublearea(V, F, dblA);
	igl::cumsum(dblA, 1, C);

	C = C / C.maxCoeff();

	FL.reserve(n);
	for (int i = 0; i < n; ++i) 
	{
		// Generate 1 sample from standard uniform R.V.
		double gamma = dist(gen);

		// Randomly select triangle
		// Feels like Eigen should have a nicer way to do this
		for (int j = 0; j < C.rows(); ++j) {
			if (C(j) > gamma) {
				FL.push_back(j);
				break;
			}
		}
	}
}