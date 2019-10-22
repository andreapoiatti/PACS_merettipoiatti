#ifndef FDAPDE_H_
#define FDAPDE_H_

// Libraries
// R
#ifdef R_VERSION_
        #define  R_NO_REMAP
        #include <R.h>
        #include <Rdefines.h>
        #include <Rinternals.h>
#endif

// C/C++
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <vector>
// #include <iomanip>

// EIGEN
//Take the code from the linked RcppEigen
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#define  EIGEN_MPL2_ONLY

// For debugging purposes
// #define  EIGEN_MPL2_ONLY
// #include "Eigen/Eigen/Dense"
// #include "Eigen/Eigen/Sparse"
// #include <Eigen/StdVector>

//----------------------------------------------------------------------------//

// Typedefs
// C++
typedef double  Real;
typedef int     UInt;

// EIGEN
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>     MatrixXr;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>                  VectorXr;
typedef Eigen::Matrix<UInt, Eigen::Dynamic, 1>                  VectorXi;
typedef Eigen::SparseMatrix<Real>                               SpMat;
typedef Eigen::SparseVector<Real>                               SpVec;
typedef Eigen::Triplet<Real>                                    coeff;

#endif /* FDAPDE_H_ */
