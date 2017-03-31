#include <quaternion.h>
#include <Eigen/SparseCore>

typedef Eigen::Matrix<Quaternion, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign> MatrixXq;
typedef Eigen::SparseMatrix<Quaternion> SparseMatrixXq;
typedef Eigen::Matrix<Quaternion, Eigen::Dynamic, 1, Eigen::AutoAlign, Eigen::Dynamic, 1> VectorXq;
typedef Eigen::Matrix<Quaternion, 1, Eigen::Dynamic, Eigen::AutoAlign, 1, Eigen::Dynamic> RowVectorXq;

void mesh_vertices_to_quaternions(const Eigen::MatrixXd V, VectorXq& U);

void quaternions_to_mesh_vertices(const VectorXq U, Eigen::MatrixXd& V);

void sparse_quaternion_matrix_to_real_matrix(const SparseMatrixXq Q, Eigen::SparseMatrix<double>& R);

void real_quaternion_matrix_to_real_sparse_matrix(const Eigen::MatrixXd QM, Eigen::SparseMatrix<double>& R);

void real_quaternion_sparse_matrix_to_real_sparse_matrix(const Eigen::SparseMatrix<double> QM, Eigen::SparseMatrix<double>& R);

void quaternion_matrix_to_real_matrix(const MatrixXq Q, Eigen::MatrixXd& R);

void insert_quaternion_block(double Q[4][4], int i, int j, Eigen::MatrixXd& RM);

void insert_quaternion_block(double Q[4][4], int i, int j, std::vector<Eigen::Triplet<double>>& triplets);
