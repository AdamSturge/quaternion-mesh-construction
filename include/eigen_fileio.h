#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>
#include <iostream>
#include <fstream>

//namespace Eigen {
//	template<class Matrix>
//	void write_binary(const std::string filename, const Matrix& matrix);
//	template<class Matrix>
//	void read_binary(const std::string filename, Matrix& matrix);
//}

namespace Eigen {
	
	template<class T>
	void write_binary_dense(const std::string filename, const Eigen::Matrix<T,-1,-1>& matrix) {
		std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		typename Eigen::Matrix<T, -1, -1>::Index rows = matrix.rows(), cols = matrix.cols();
		out.write((char*)(&rows), sizeof(typename Eigen::Matrix<T, -1, -1>::Index));
		out.write((char*)(&cols), sizeof(typename Eigen::Matrix<T, -1, -1>::Index));
		out.write((char*)matrix.data(), rows*cols * sizeof(typename Eigen::Matrix<T, -1, -1>::Scalar));
		out.close();
	}

	template<class T>
	void read_binary_dense(const std::string filename, Eigen::Matrix<T,-1,-1>& matrix) {
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		typename Eigen::Matrix<T,-1,-1>::Index rows = 0, cols = 0;
		in.read((char*)(&rows), sizeof(typename Eigen::Matrix<T, -1, -1>::Index));
		in.read((char*)(&cols), sizeof(typename Eigen::Matrix<T,-1,-1>::Index));
		matrix.resize(rows, cols);
		in.read((char *)matrix.data(), rows*cols * sizeof(typename Eigen::Matrix<T,-1,-1>::Scalar));
		in.close();
	}

	template <typename T>
	void write_binary_sparse(const std::string filename, Eigen::SparseMatrix<T>& X) {
		std::vector<Eigen::Triplet<T>> res;
		int sz = X.nonZeros();
		X.makeCompressed();

		std::fstream writeFile;
		writeFile.open(filename, std::ios::binary | std::ios::out);

		if (writeFile.is_open())
		{
			Eigen::SparseMatrix<T>::Index rows, cols, nnzs, outS, innS;
			rows = X.rows();
			cols = X.cols();
			nnzs = X.nonZeros();
			outS = X.outerSize();
			innS = X.innerSize();

			writeFile.write((const char *)&(rows), sizeof(Eigen::SparseMatrix<T>::Index));
			writeFile.write((const char *)&(cols), sizeof(Eigen::SparseMatrix<T>::Index));
			writeFile.write((const char *)&(nnzs), sizeof(Eigen::SparseMatrix<T>::Index));
			writeFile.write((const char *)&(outS), sizeof(Eigen::SparseMatrix<T>::Index));
			writeFile.write((const char *)&(innS), sizeof(Eigen::SparseMatrix<T>::Index));

			writeFile.write((const char *)(X.valuePtr()), sizeof(T) * X.nonZeros());
			writeFile.write((const char *)(X.outerIndexPtr()), sizeof(Eigen::SparseMatrix<T>::Index) * X.outerSize());
			writeFile.write((const char *)(X.innerIndexPtr()), sizeof(Eigen::SparseMatrix<T>::Index) * X.nonZeros());

			writeFile.close();
		}
	}

	template <typename T>
	void read_binary_sparse(const std::string filename, Eigen::SparseMatrix<T>& X) {
		std::fstream readFile;
		readFile.open(filename, std::ios::binary | std::ios::in);
		if (readFile.is_open())
		{
			Eigen::SparseMatrix<T>::Index rows, cols, nnz, inSz, outSz;
			readFile.read((char*)&rows, sizeof(Eigen::SparseMatrix<T>::Index));
			readFile.read((char*)&cols, sizeof(Eigen::SparseMatrix<T>::Index));
			readFile.read((char*)&nnz, sizeof(Eigen::SparseMatrix<T>::Index));
			readFile.read((char*)&inSz, sizeof(Eigen::SparseMatrix<T>::Index));
			readFile.read((char*)&outSz, sizeof(Eigen::SparseMatrix<T>::Index));

			X.resize(rows, cols);
			X.makeCompressed();
			X.resizeNonZeros(nnz);

			readFile.read((char*)(X.valuePtr()), sizeof(T) * nnz);
			readFile.read((char*)(X.outerIndexPtr()), sizeof(Eigen::SparseMatrix<T>::Index) * outSz);
			readFile.read((char*)(X.innerIndexPtr()), sizeof(Eigen::SparseMatrix<T>::Index) * nnz);

			X.finalize();
			readFile.close();

		} // file is open
	}

}