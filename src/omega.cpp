#include <omega.h>
#include <algorithm>
#include <quaternion.h>
#include <iostream>

void omega(
	const Eigen::MatrixXd V,
	const Eigen::MatrixXi F,
	const Eigen::VectorXd lam,
	Eigen::VectorXd& om
) 
{
	int n_F = F.rows();
	Eigen::RowVector3d f0;
	Eigen::RowVector3d f1;
	Eigen::RowVector3d f2;
	Quaternion lambda_a;
	Quaternion lambda_b;
	Quaternion e;
	Quaternion e_tilde;
	Eigen::Vector3d u1;
	Eigen::Vector3d u2;
	double d_prod;
	double c_prod_norm;
	double cotAlpha;
	Quaternion result;
	for (int i = 0; i < n_F; ++i)
	{
		// get indices of the vertices of this face
		int v[3] = { F(i,0), F(i,1), F(i,2) };
		
		for(int j = 0; j < 3; j++)
		{
			f0 = V.row(v[(j + 0) % 3]);
			f1 = V.row(v[(j + 1) % 3]);
			f2 = V.row(v[(j + 2) % 3]);

			// determine orientation of this edge
			int a = v[(j + 1) % 3];
			int b = v[(j + 2) % 3];
			if (a > b)
			{
				std::swap(a, b);
			}

			lambda_a = Quaternion(lam(4 * a + 0), Eigen::Vector3d(lam(4 * a + 1), lam(4 * a + 2), lam(4 * a + 3)));
			lambda_b = Quaternion(lam(4 * b + 0), Eigen::Vector3d(lam(4 * b + 1), lam(4 * b + 2), lam(4 * b + 3)));
			e = Quaternion(0.0, V.row(b)) - Quaternion(0.0, V.row(a));
			e_tilde = 
				(1.0 / 3.0) * (~lambda_a) * e * lambda_a +
				(1.0 / 6.0) * (~lambda_a) * e * lambda_b +
				(1.0 / 6.0) * (~lambda_b) * e * lambda_a +
				(1.0 / 3.0) * (~lambda_b) * e * lambda_b ;

			u1 = (f1 - f0);
			u2 = (f2 - f0);
			d_prod = u1.dot(u2);
			c_prod_norm = u1.cross(u2).norm();
			cotAlpha = d_prod / c_prod_norm;

			result = cotAlpha * e_tilde / 2.0;

			om(4 * a + 0) -= result.re();
			om(4 * a + 1) -= result.im()(0);
			om(4 * a + 2) -= result.im()(1);
			om(4 * a + 3) -= result.im()(2);

			om(4 * b + 0) += result.re();
			om(4 * b + 1) += result.im()(0);
			om(4 * b + 2) += result.im()(1);
			om(4 * b + 3) += result.im()(2);
		}

		// DO I NEED TO REMOVE MEAN?
		om = om - Eigen::VectorXd::Constant(om.rows(),om.cols(),om.mean());
		
	}
}