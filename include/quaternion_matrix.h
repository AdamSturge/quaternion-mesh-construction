#include <quaternion.h>
#include <Eigen/SparseCore>

typedef Eigen::Matrix<Quaternion, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign> MatrixXq;
typedef Eigen::SparseMatrix<Quaternion> SparseMatrixXq;
typedef Eigen::Matrix<Quaternion, Eigen::Dynamic, 1, Eigen::AutoAlign, Eigen::Dynamic, 1> VectorXq;
typedef Eigen::Matrix<Quaternion, 1, Eigen::Dynamic, Eigen::AutoAlign, 1, Eigen::Dynamic> RowVectorXq;

// Converts a vertex matrix into a quaternion vector
// Inputs:
//      V   Nx3 original mesh
// Outputs:
//      U  4*Nx1 quaternionic vector representation
void mesh_vertices_to_quaternions(const Eigen::MatrixXd V, VectorXq& U);

// Converts a quaternion vector to a vertex matrix
// Inputs:    
//      U  4*Nx1 quaternionic vector representation
// Outputs:
//      V   Nx3 original mesh
void quaternions_to_mesh_vertices(const VectorXq U, Eigen::MatrixXd& V);

// Converts a sparse quaternion matrix to a sparse real matrix
// Inputs:    
//      Q  NxM quaternion matrix
// Outputs:
//      R   4Nx4M real matrix
void sparse_quaternion_matrix_to_real_matrix(const SparseMatrixXq Q, Eigen::SparseMatrix<double>& R);

// Converts a dense quaternion matrix to a sparse real matrix
// Inputs:    
//      Q  NxM quaternion matrix
// Outputs:
//      R   4Nx4M real matrix
void real_quaternion_matrix_to_real_sparse_matrix(const Eigen::MatrixXd QM, Eigen::SparseMatrix<double>& R);

// Converts a sparse real matrix to its quaternion representation
// Inputs:    
//      Q  NxM quaternion matrix
// Outputs:
//      R   4Nx4M real matrix
void real_quaternion_sparse_matrix_to_real_sparse_matrix(const Eigen::SparseMatrix<double> QM, Eigen::SparseMatrix<double>& R);

// Converts a dense quaternion matrix to a dense real matrix
// Inputs:    
//      Q  NxM quaternion matrix
// Outputs:
//      R   4Nx4M real matrix
void quaternion_matrix_to_real_matrix(const MatrixXq Q, Eigen::MatrixXd& R);

// Inserts a 4x4 real block representation of a quaternion into a dense matrix
// Inputs:    
//      Q	4x4 real matrix representation of quaternion
//			i	row index of top left corner of matrix block to be inserted
//		  j	column index of top left corner of matrix block to be inserted
// Outputs:
//      RM same as original matrix RM but with a 4x4 block inserted
void insert_quaternion_block(double Q[4][4], int i, int j, Eigen::MatrixXd& RM);

// Inserts a 4x4 real block representation of a quaternion into a sparse matrix
// Inputs:    
//      Q	4x4 real matrix representation of quaternion
//			i	row index of top left corner of matrix block to be inserted
//		  j	column index of top left corner of matrix block to be inserted
// Outputs:
//      triplets same as original triplets but with entries added for the block matrix
void insert_quaternion_block(double Q[4][4], int i, int j, std::vector<Eigen::Triplet<double>>& triplets);
