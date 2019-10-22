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
	        VectorXr 	z_hat_; 			// [z_hat] computed in Computedegreesoffreedomexact for the gcv

		// Other utility matrices
		MatrixXr 	H_;				// Hat matrix, Col(W) projector:. W*(W^t*W)^{-1}*W^t
		MatrixXr 	Q_;				// Orth(Col(W)) projector = Identity-H_
		MatrixXr 	R_; 				// R1^t * R0^{-1} * R1  already stored
		MatrixXr 	SS_; 				//_stores Psi^t*Q*Psi+lambda*R
		MatrixXr 	V_; 				//_stores the values of SS^{-1}*Psi^t*Q
                MatrixXr        S_;                             //_matrix S of Stu-Hunter Sangalli
		// [[old deprecated matrices]]
		// SpMat 	Psi_;
		// SpMat 	NWblock_;
		// MatrixXr 	P_;

		// Sparse matrices for fast solution
		SpMat 			DMat_;			// Data matrix: top-left block of coeffmatrix_
		SpMat 			AMat_;			// Data matrix: diagonal of coeffmatrix_
		SpMat 			MMat_;			// Data matrix: down-right block of coeffmatrix_

		SpMat 			_coeffmatrix;        	//! A Eigen::SpMat: Stores the system right hand side.
		VectorXr 		_b;                     //! A Eigen::VectorXr: Stores the system right hand side.
		std::vector<VectorXr> 	_solution; 		//! A vector of solutions: Stores the system solution for various lambda [[!!! TODO !!!]]
		std::vector<Real> 	_dof;			// Stores tr([???]) for various lambda [[!!! TODO !!!]]

		// Factorized matrices and controllers
		Eigen::SparseLU<SpMat> 		Adec_; 		// Stores the factorization of A_
		Eigen::PartialPivLU<MatrixXr> 	Gdec_;		// Stores factorization of G =  C + [V * A^-1 * U]
		Eigen::PartialPivLU<MatrixXr> 	WTWinv_;	// Stores the factorization of W^t * W, NOT of the inverse

		bool 			isWTWfactorized_; 	// Checks if the factorization WTWinv_ is already computed, false in constructor
		bool 			isRcomputed_;		// Checks if matrix R is computed

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

		// DOF utilities
		void 	computeDegreesOfFreedom(UInt output_index, Real lambda);
		void 	computeDegreesOfFreedomExact(UInt output_index, Real lambda);
		void 	computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);

		// GCV utilities [POIATTI]
	        Real 	computeGCV(UInt output_index); 			//_usarla in apply, per testare se va, facendo stampare qualcosa
		Real 	computeGCV_derivative(UInt output_index); 	//_usarla in apply, per testare se va, facendo stampare qualcosa

		// Factorizer
		void 	system_factorize();

		// Solver
		template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived> &);

		// [[old deprecated uilities]]
		//void computeBasisEvaluations();
		//void computeProjOnCovMatrix();

		// [[old deprecated parts]]
		//! A normal member taking two arguments: Dirichlet Boundary condition
		/*!
		 * This member applies Dirichlet boundary conditions on the linear system with the penalization method.
		  \param bcindex is a const reference to vector<int> : the global indexes of the nodes to which the boundary condition has to be applied.
		  \param bcvalues is a const reference to vector<double> : the values of the boundary conditions relative to bcindex.
		/*
		// void applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues);
		// void computeDataMatrix(SpMat& DMat);
		// void computeDataMatrixByIndices(SpMat& DMat);
		// void computeRightHandData(VectorXr& rightHandData);
		// void computeDegreesOfFreedom(UInt output_index);

		//! A template for the system resolution: SpLu, SpQR, SpCholesky,SpConjGrad
		// template<typename P>
		// void solve(UInt output_index);
		*/

	public:
		// Constructors
		//! A Constructor.
		MixedFERegressionBase(const MeshHandler<ORDER, mydim, ndim> & mesh, const InputHandler & regressionData):
			mesh_(mesh), regressionData_(regressionData), isRcomputed_(false), isWTWfactorized_(false) {};

		// Main function
		template<typename A>
		void apply(EOExpr<A> oper);

		// Getters
		//! A inline member that returns a VectorXr, returns the whole solution_.
		inline std::vector<VectorXr> const & 	getSolution()	const {return _solution;}
		inline std::vector<Real> const &	getDOF() 	const {return _dof;} //_serve solo se si caclola la GCV in esterno a C++, noi possiamo non esportare per risparmiare
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
