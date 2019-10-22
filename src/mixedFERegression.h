#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

// Headers
#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include <memory>

// Classes
//! A LinearSystem class: A class for the linear system construction and resolution.
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegressionBase
{
	protected:
		// *** DATA ***
		// Techinical coefficients
		static constexpr Real dirichlet_penalization 	= 10e18;
		static constexpr Real pruning_coeff 		= 2.2204e-013;

		// Mesh related data
		const MeshHandler<ORDER, mydim, ndim> & mesh_;
		const InputHandler & 			regressionData_;
		std::vector<coeff> 			tripletsData_;

		// Matrices for linar system and estimates
		SpMat 		A_;				// System matrix, with psi^T*psi in north-west block
		SpMat 		R1_;				// North-east block of system matrix A_
		SpMat 		R0_;				// South-east block of system matrix A_
		SpMat 		psi_;				// Psi(p): num_nodes x [n] matrix
		MatrixXr 	U_;				// Psi^t*W padded with zeros

		// Other utility matrices
		MatrixXr 	H_;				// stores the value of W*(W^t*W)^{-1}*W^t [Hat matrix, Col(W) projector]
		MatrixXr 	Q_;				// stores the value of Identity-H_ [Orth(Col(W)) projector]

		// Sparse matrices for generic setting and lambda
		SpMat 		DMat_;				// Data matrix: top-left block of coeffmatrix_
		SpMat 		AMat_;				// Data matrix: diagonal of coeffmatrix_
		SpMat 		MMat_;				// Data matrix: down-right block of coeffmatrix_
		Real		lambda_;			// The value of the best penalization according to some optimization criterion

		SpMat 		_coeffmatrix;        		// A Eigen::SpMat: Stores the system right hand side.
		VectorXr 	_b;                    		// A Eigen::VectorXr: Stores the system right hand side.
		VectorXr	_solution; 			// A vector of solutions: Stores the system solution for the given lambda

		// Factorized matrices and controllers
		Eigen::SparseLU<SpMat> 		Adec_; 		// Stores the factorization of A_
		Eigen::PartialPivLU<MatrixXr> 	Gdec_;		// Stores factorization of G =  C + [V * A^-1 * U]
		Eigen::PartialPivLU<MatrixXr> 	WTWinv_;	// Stores the factorization of W^t * W, NOT of the inverse

		bool 			isWTWfactorized_; 	// Checks if the factorization WTWinv_ is already computed, false in constructor
		bool			isRegression_;		// Checks if the model has covariates or not
		bool			hasLocations_;		// Checks if the model has locations or uses nodes

		// *** METHODS ***
		// Boundary conditions
		void 	addDirichletBC();

		// Setters and builders
		// Note that builders are slightly more general allowing to pass them data different
		// from the ones that are already stord inside the class
		void 	setPsi();									// Computes psi_
		void 	setH();										// Computes H_
		void 	setQ();										// Computes Q_

		void 	buildA(const SpMat & Psi,  const SpMat & R1,  const SpMat & R0);		// std builder of matrix A_
		void 	getDataMatrix(SpMat & DMat);							// std builder of DMat_
		void 	getDataMatrixByIndices(SpMat & DMat);						// std builder of DMat_ without using phi
		void 	buildCoeffMatrix(const SpMat & DMat,  const SpMat & AMat,  const SpMat & MMat);	// std builder of coeffmatrix_
		void 	getRightHandData(VectorXr & rightHandData);					// std builder of b_

		// Utilities for faster computation
		MatrixXr LeftMultiplybyQ(const MatrixXr & u);

		// Factorizer
		void 	system_factorize();

		// Solver
		template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived> &);

	public:
		// Constructors
		//! A Constructor.
		MixedFERegressionBase(const MeshHandler<ORDER, mydim, ndim> & mesh, const InputHandler & regressionData):
			mesh_(mesh), regressionData_(regressionData), isRegression_(false), isWTWfactorized_(false), hasLocations_(false) {};

		void 	set_lambda (Real lambda) {lambda_ = lambda;}

		// Main function
		template<typename A>
		void apply(EOExpr<A> oper);

		void apply_2();

		// Getters
		//! A inline member that returns a VectorXr, returns the whole solution_.
		inline VectorXr			getSolution(void)		const {return _solution;} 	//returns f and g
		inline const InputHandler *	getData(void)			      {return &regressionData_;}
		inline SpMat &			getPsi_(void)			      {return psi_;}
		inline SpMat &			getR0_(void)			      {return R0_;}
		inline SpMat &			getR1_(void)			      {return R1_;}
		inline SpMat &			get_DMat_(void)			      {return DMat_;}
		inline MatrixXr &		getQ_(void)			      {return Q_;}
		inline MatrixXr	&		getH_(void)			      {return H_;}
		inline Real 			getlambda_(void)		      {return lambda_;}
		inline bool			checkisRegression_(void)	const {return isRegression_;}	// Checks if the model has covariates or not
		inline bool			checkhasLocations_(void)	const {return hasLocations_;}	// Checks if the model has locations or uses nodes
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>
{
	public:
		MixedFERegression(const MeshHandler<ORDER, ndim, mydim> & mesh, const InputHandler & regressionData):
			MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, regressionData) {};

		void apply()
		{
			std::cout << "Option not implemented! \n";
		}
};

#include "mixedFERegression_imp.h"

#endif
