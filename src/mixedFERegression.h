#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "optimizationData.h"
#include "solver.h"
#include "integratePsi.h"
#include "auxiliary_optimizer.h"
#include <memory>

/*! A base class for the smooth regression.
*/
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegressionBase
{
	protected:
		// Fundamental data
		const MeshHandler<ORDER, mydim, ndim> & mesh_;
		const InputHandler & regressionData_;
		const OptimizationData & optimizationData_;

		//  system matrix = 	|psi^T * A *psi | lambda R1^T  |   +  |psi^T * A * (-H) * psi |  O |   =  matrixNoCov + matrixOnlyCov
		//	                |     R1        | R0	       |      |         O             |  O |

		// Data matrices
		SpMat 		matrixNoCov_;	//! System matrix with psi^T*psi [ptw data] or psi^T*A*psi [areal data] in north-west block  (full system matrix in case of no covariates)
		SpMat		R1_;		//! North-east block of system matrix matrixNoCov_
		SpMat 		R0_;		//! South-east block of system matrix matrixNoCov_
		SpMat 		psi_;		//! Psi matrix of the model
		SpMat 		psi_t_; 	//! Transpose o Psi matrix of the model
		MatrixXr 	H_; 		//! The hat matrix of the regression
		MatrixXr	Q_; 		//! Identity - H, projects onto the orthogonal subspace
		VectorXr 	A_; 		//! A_.asDiagonal() = diag(|D_1|,...,|D_N|) areal matrix, = identity nnodesxnnodes if pointwise data
		MatrixXr 	U_;		//! psi^T * W [prw] or psi^T * A * W [areal] padded with zeros, needed for Woodbury decomposition
		MatrixXr 	V_;   		//! W^T*psi, if pointwise data is U^T, needed for Woodbury decomposition

		// Factorizaions
		Eigen::SparseLU<SpMat> matrixNoCovdec_; // Stores the factorization of matrixNoCov_
		Eigen::PartialPivLU<MatrixXr> Gdec_;	// Stores factorization of G =  C + [V * matrixNoCov^-1 * U]
		Eigen::PartialPivLU<MatrixXr> WTW_;	// Stores the factorization of W^T * W
		bool isWTWfactorized_ = false;

		// Terms for linear system
		VectorXr 	_rightHandSide; 	//!A Eigen::VectorXr: Stores the system right hand side.
		VectorXr 	_forcingTerm;
		VectorXr 	_solution; 		//!A Eigen::VectorXr: Stores the system solution.

		// Additional conditions
		bool isSpaceVarying = false; // used to distinguish whether to use the forcing term u in apply() or not

		// Setters
		//! A member function computing the Psi matrix
		void setPsi(void);
		//! A member function which builds the Q matrix
		void setQ();
		//! A member function which builds the H matrix
		void setH();
		//! A member function which builds the A vector containing the areas of the regions in case of areal data
		void setA();

		// Utility
		//! A function that given a vector u, performs Q*u efficiently
		MatrixXr LeftMultiplybyQ(const MatrixXr& u);

		// Builders
		//! A function which adds Dirichlet boundary conditions before solving the system ( Remark: BC for areal data are not implemented!)
		void addDirichletBC();
		//! A member function computing the no-covariates version of the system matrix
		void buildMatrixNoCov(const SpMat & Psi, const SpMat & R1, const SpMat & R0);
		//! A member function returning the system right hand data
		void getRightHandData(VectorXr& rightHandData);

		// Factorizer
	    	//! A function to factorize the system, using Woodbury decomposition when there are covariates
		void system_factorize();

		// Solver
		//! A function which solves the factorized system
		template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived> &);

	public:
		//!A Constructor.
		MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim> & mesh, const InputHandler & regressionData, const OptimizationData & optimizationData):
		 mesh_(mesh), regressionData_(regressionData), optimizationData_(optimizationData){};

		// [[TODO]]
		//! The function solving the system, used by the children classes. Saves the result in _solution
		/*!
		    \param oper an operator, which is the Stiffness operator in case of Laplacian regularization
		    \param u the forcing term, will be used only in case of anysotropic nonstationary regression
		*/
		template<typename A>
		void preapply(EOExpr<A> oper,const ForcingTerm & u);
		void apply(Real lambda);

		// Getters
		//! A inline member that returns a VectorXr, returns the whole solution_.
		inline const InputHandler *	getData(void)			      {return &regressionData_;}
		inline const SpMat *		getPsi_(void)			      {return &psi_;} //[[& PER OTTIMIZZARE, CONST???]]
		inline const SpMat *		getPsi_t(void)			      {return &psi_t_;} //[[& PER OTTIMIZZARE, CONST???]]
		inline const SpMat *		getR0_(void)			      {return &R0_;}
		inline const SpMat *		getR1_(void)			      {return &R1_;}
		//[[??]]inline const SpMat *		get_DMat_(void)			      {return &DMat_;}
		inline const MatrixXr *		getQ_(void)			      {return &Q_;}
		inline const MatrixXr *		getH_(void)			      {return &H_;}
		inline bool			checkisRegression_(void)	const {return (this->regressionData_.getCovariates()->cols()!=0 && this->regressionData_.getCovariates()->rows()!=0);}	// Checks if the model has covariates or not
		inline bool			check_is_loc_by_n(void)		const {return regressionData_.isLocationsByNodes();}	// Checks if the model has locations or uses nodes

		//! A inline member that returns a VectorXr, returns the whole solution_.
		inline VectorXr const & getSolution() const{return _solution;};
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>
{
	public:
		MixedFERegression(const MeshHandler<ORDER, ndim, mydim> & mesh, const InputHandler & regressionData, const OptimizationData & optimizationData):
			MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, regressionData, optimizationData){};

		void preapply()
		{
			std::cout << "Option not implemented! \n";
		}
};

#include "mixedFERegression_imp.h"

#endif
