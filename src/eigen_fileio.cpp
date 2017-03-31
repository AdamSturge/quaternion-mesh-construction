#include <eigen_fileio.h>
//namespace Eigen {
//	
//	template<class Matrix>
//	void write_binary(const std::string filename, const Matrix& matrix) {
//		std::ofstream out(filename, ios::out | ios::binary | ios::trunc);
//		typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
//		out.write((char*)(&rows), sizeof(typename Matrix::Index));
//		out.write((char*)(&cols), sizeof(typename Matrix::Index));
//		out.write((char*)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
//		out.close();
//	}
//
//	template<class Matrix>
//	void read_binary(const std::string filename, Matrix& matrix) {
//		std::ifstream in(filename, ios::in | std::ios::binary);
//		typename Matrix::Index rows = 0, cols = 0;
//		in.read((char*)(&rows), sizeof(typename Matrix::Index));
//		in.read((char*)(&cols), sizeof(typename Matrix::Index));
//		matrix.resize(rows, cols);
//		in.read((char *)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
//		in.close();
//	}
//}