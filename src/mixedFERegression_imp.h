#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>
#include <algorithm>

#include "R_ext/Print.h"

//#include <libseq/mpi.h>
#include "../inst/include/dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

//----------------------------------------------------------------------------//
// Setters
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setPsi(void)
{
	// Psi is a nlocations x nnodes  matrix, first fetch the dimensions
	UInt nnodes 	= mesh_.num_nodes();				// #cols
	UInt nlocations = regressionData_.getNumberofObservations();	// #rows
	psi_.resize(nlocations, nnodes);				// resize the empty matrix

	// [[POINTWISE DATA case]]
	// Optimized strategies according to the presence of locations
	if (regressionData_.isLocationsByNodes())
	{
		// THEORETICAL REMARK:
		// If isLocationsByNodes is active it entails that nodes are used as locations
		// However, maybe, NOT IN ALL nodes evaluations are performed: some might be utility nodes not
		// associated with values:
		// e.g. think about a triangle with 6 nodes, probably only the main 3 might be associated
		// with a value with the others not, but this is only one possible choice:
		// The efficent fact is that anyway Phi_j(p_i) will be either 1 if i is the Lagrange
		// node of j or 0 if it is not. Note that since Phi are as many as the nodes, but only
		// some of the nodes are locations some columns of Phi = [Phi_j(p_i)]_ij might be of
		// all zeros: i.e. that Phi_j is linked to a node which is not a location since it's value is NA.
		// In such a case Phi is NOT a sqare matrix and there are more culumns than rows
		// If all nodes == locations then Phi == square Identity, but this is just a particular case
		// and, even though isLocationsByNodes is true, this might not be veriefied

		std::vector<coeff> tripletAll;
		const std::vector<UInt> * k = regressionData_.getObservationsIndices();
		UInt k_size = k->size();

		tripletAll.reserve(k_size);
		for(UInt i=0; i<k_size; ++i)
		{
			// Add a value 1 for each valid index in row i
			// and column k[i] (under psi_k[i], the associated node)
			tripletAll.push_back(coeff(i,(*k)[i],1.0));
		}

		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	}
	// Still pointwse case but with autonomous locations
	else if(regressionData_.getNumberOfRegions() == 0)
	{
		// THEORETICAL REMARK:
		// If isLocationsByNodes is false locations are unspecified points
		// so we have to evaluate them directly

		constexpr UInt 			Nodes = (mydim==2) ? (3*ORDER) : (6*ORDER-2);
		Element<Nodes, mydim, ndim> 	tri_activated;
		Eigen::Matrix<Real,Nodes,1> 	coefficients;
		Real 				evaluator;
		const std::vector<Point> *	pl = regressionData_.getLocations();

		for(UInt i=0; i<nlocations; i++)
		{
			// Search the element containing the point
			tri_activated = mesh_.findLocationNaive((*pl)[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				// If not found
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}
			else
			{
				// Found, it's action might be felt a priori by all the psi of the elemnt, one for each node
				for(UInt node=0; node<Nodes; ++node)
				{
					// Define vector of all zeros but "node" component (necessary for function evaluate_point)
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base

					// Evaluate psi_node in the node and put it in the right place in the matrix according to global index
					evaluator = evaluate_point<Nodes,mydim,ndim>(tri_activated, (*pl)[i], coefficients);
					// Insert the value in the column of the global indexing of the evaluated node
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}
	}
	else //areal data
	{
		constexpr UInt 		Nodes = (mydim==2) ? (3*ORDER) : (6*ORDER-2);
		Real * 			tab = (Real*) malloc(sizeof(Real)*nnodes); // Psi_i temporary storage
		const MatrixXi * 	pi  = regressionData_.getIncidenceMatrix();

		for(UInt i=0; i<nlocations; i++) //nlocations = number of regions
		{
			for(UInt k=0; k<nnodes; k++)
				tab[k]=0;

			for(UInt j=0; j<mesh_.num_elements(); j++)
			{
				if((*pi)(i,j) == 1) //element j is in region i
				{
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j);
					for(UInt k=0; k<Nodes; k++)
					{
						tab[tri[k].getId()] += integratePsi(tri,k); // integral over tri of psi_k
					}
				}
			}

			for(int k=0; k<nnodes; k++)
				if(tab[k] != 0)
					psi_.insert(i,k) = tab[k]/A_(i); //divide by |D_i|
		}
		free(tab);
	}
	psi_.makeCompressed(); // Compress sparse matrix
	psi_t_ = SpMat(psi_.transpose());
	psi_t_.makeCompressed(); // Compress sparse matrix
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setH(void)
{
	// Debug text
	// std::cout << "Computing Projection Matrix" << std::endl;
	UInt nlocations = regressionData_.getNumberofObservations();
	const MatrixXr * Wp(this->regressionData_.getCovariates());
	bool ilbn = regressionData_.isLocationsByNodes();

	if(ilbn)
	{
		const std::vector<UInt> * k = regressionData_.getObservationsIndices();

		// Some rows might be discarded [[we probably have data for every node not only the points ???]]
		UInt n_Wcols = Wp->cols();
		MatrixXr * Wp_reduced  = new MatrixXr;
		Wp_reduced->resize(regressionData_.getNumberofObservations(), n_Wcols);

		for (UInt i=0; i<nlocations; ++i)
		{
		       UInt index_i = (*k)[i];
		       for (auto j=0; j<n_Wcols; ++j)
		       {
			       (*Wp_reduced)(i,j) = (*Wp)(index_i,j);
		       }
		}
		Wp = Wp_reduced;
	}

	MatrixXr Wt(Wp->transpose());		// Compute W^t
	H_ = (*Wp)*(Wt*(*Wp)).ldlt().solve(Wt);	// using cholesky LDLT decomposition for computing hat matrix

	if(ilbn)
		delete Wp;
}

template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setQ(void)
{
	// Debug text
 	// std::cout << "Computing Orthogonal Space Projection Matrix" << std::endl;

	// Remember Q = I-H
	Q_.resize(H_.rows(), H_.cols());	// Resizing dimensions as H
	Q_ = -H_;
	for (UInt i=0; i<H_.rows(); ++i)
	{
		Q_(i,i) += 1;			// Summing the identity by rows (or columns)
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setA(void)
{
	UInt nRegions = regressionData_.getNumberOfRegions();
	const MatrixXi * pi  = regressionData_.getIncidenceMatrix();
	A_.resize(nRegions,1);

	for (UInt i=0; i<nRegions; i++)
	{
		A_(i) = 0;
		for (UInt j=0; j<pi->cols(); j++)
		{
			if ((*pi)(i,j) == 1)
			{
				A_(i) += mesh_.elementMeasure(j);
			}
		}
	}
}

//----------------------------------------------------------------------------//
// Utilities [[NOT VERY OPTMIZED, SENSE??, we have Q and P...]]

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::LeftMultiplybyQ(const MatrixXr & u)
{
	const MatrixXr * Wp = regressionData_.getCovariates();
	if (Wp->rows() == 0)
	{
		// Q is the projecton on Col(W) if W == 0 => Q = Identity
		return u;
	}
	else
	{
		MatrixXr Wt = Wp->transpose();
		// Check factorization, if not present factorize the matrix W^t*W
		if (isWTWfactorized_ == false )
		{
			WTW_.compute(Wt*(*Wp));
			isWTWfactorized_ = true;
		}
		// Compute H (or I-Q) the projection on Col(W) and multiply it times u
		return u - (*Wp)*WTW_.solve(Wt*u);  //note (I-P)*u == Q*u
	}

}

//----------------------------------------------------------------------------//
// Builders

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::addDirichletBC()
{
	UInt id1,id3;

	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt> & bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real> & bc_values = regressionData_.getDirichletValues();

	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for(UInt i=0; i<nbc_indices; i++)
	{
		id1 = bc_indices[i];
		id3 = id1+nnodes;

		matrixNoCov_.coeffRef(id1,id1) = pen;
		matrixNoCov_.coeffRef(id3,id3) = pen;

		_rightHandSide(id1) = bc_values[i]*pen;
		_rightHandSide(id3) = 0;
	}

	matrixNoCov_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildMatrixNoCov(const SpMat & Psi,  const SpMat & R1,  const SpMat & R0) {

	UInt nnodes = mesh_.num_nodes();

	SpMat DMat;
	if(regressionData_.getNumberOfRegions()==0) // pointwise data
		DMat = Psi.transpose()*Psi;
	else                                        // areal data: need to add the diag(|D_1|,...,|D_N|)
		DMat = Psi.transpose()*A_.asDiagonal()*Psi;

	    // Vector to be filled with the triplets used to build _coeffmatrix (reserved with the right dimension)
	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());

	// Parsing all matrices, reading the values to be put inside _coeffmatrix, coordinates according to the rules
	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(), it.value()));
		}
	for (int k=0; k<R0.outerSize(); ++k)
		for (SpMat::InnerIterator it(R0,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes, it.value()));
		}
	for (int k=0; k<R1.outerSize(); ++k)
		for (SpMat::InnerIterator it(R1,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.col(), it.row()+nnodes, it.value()));
		}
	for (int k=0; k<R1.outerSize(); ++k)
	  	for (SpMat::InnerIterator it(R1,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
		}

	// Define, resize, fill and compress _coeffmatrix
	matrixNoCov_.setZero();
	matrixNoCov_.resize(2*nnodes,2*nnodes);
	matrixNoCov_.setFromTriplets(tripletAll.begin(), tripletAll.end());
	matrixNoCov_.makeCompressed();

	// Debugging purpose
	// std::cout << "Coefficients' Matrix Set Correctly" << std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::getRightHandData(VectorXr & rightHandData)
{
	UInt 			nnodes = mesh_.num_nodes();
	UInt 			nlocations = regressionData_.getNumberofObservations();
	UInt 			nregions = regressionData_.getNumberOfRegions();
	const MatrixXr * 	Wp = regressionData_.getCovariates();
	const VectorXr * 	op = regressionData_.getObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if (Wp->rows() == 0) //no covariate
	{
		if (regressionData_.isLocationsByNodes())
		{
			const std::vector<UInt> * k = regressionData_.getObservationsIndices();

			for (auto i=0; i<nlocations; ++i)
			{
				UInt index_i = (*k)[i];
				rightHandData(index_i) = (*op)[i];
			}
		}
		else if (nregions == 0) //pointwise data
		{
			rightHandData = psi_t_*(*op);
		}
		else //areal data
		{
			rightHandData = psi_t_*A_.asDiagonal()*(*op);
		}
	}
	else if (nregions == 0) //with covariates, pointwise data
	{
		rightHandData = psi_t_*LeftMultiplybyQ(*op);
	}
	else //with covariates, areal data
	{
		rightHandData = psi_t_*A_.asDiagonal()*LeftMultiplybyQ(*op);
	}
}

//----------------------------------------------------------------------------//
// Factorizer

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::system_factorize(void) {

	UInt 			nnodes = mesh_.num_nodes();
	const MatrixXr * 	Wp = regressionData_.getCovariates();

	// First phase: Factorization of matrixNoCov
	matrixNoCovdec_.compute(matrixNoCov_);

	// If Q!=Identity
	if (Wp->rows() != 0)
	{
		// Second phase: factorization of matrix  G =  C + [V * matrixNoCov^-1 * U]= C + D

		// Definition of matrix U = [ psi^T * A * W | 0 ]^T and V= [ W^T*psi| 0]
		U_ = MatrixXr::Zero(2*nnodes, Wp->cols());
		V_ = MatrixXr::Zero(Wp->cols(),2*nnodes);
		MatrixXr Wt(Wp->transpose());

		if(regressionData_.getNumberOfRegions() == 0) // pointwise data
		{
			MatrixXr utility = Wt*psi_;
			V_.leftCols(nnodes) = utility;
		 	U_.topRows(nnodes) = utility.transpose();
		}
		else	//areal data
		{
			V_.leftCols(nnodes) = Wt*psi_;
		 	U_.topRows(nnodes) = psi_t_*A_.asDiagonal()*(*Wp);
		}

		MatrixXr D = V_*matrixNoCovdec_.solve(U_);

		// G = C + D
		MatrixXr G = -Wt*(*Wp) + D;
		Gdec_.compute(G);

	}
}

//----------------------------------------------------------------------------//
// Solver

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> & b)
{
	const MatrixXr * Wp = regressionData_.getCovariates();

	// Resolution of the system matrixNoCov * x1 = b
	MatrixXr x1 = matrixNoCovdec_.solve(b);

	if (Wp->rows() != 0)
	{
		// Resolution of G * x2 = V * x1
		MatrixXr x2 = Gdec_.solve(V_*x1);

		// Resolution of the system matrixNoCov * x3 = U * x2
		x1 -= matrixNoCovdec_.solve(U_*x2);
	}

	return x1;
}

//----------------------------------------------------------------------------//
// Main functions

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::preapply(EOExpr<A> oper, const ForcingTerm & u)
{
	const MatrixXr * Wp = regressionData_.getCovariates();
	UInt nnodes = mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	setA();		// Set matrix A
	setPsi();	// Set matrix Psi

	// If there are covariates in the model set H and Q
	if(Wp->rows() != 0)
	{
		setH();
		setQ();
	}

	//Assemble matrices R0-R1
	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);

	// Define Right hand data
	VectorXr rightHandData;
	getRightHandData(rightHandData);
	this->_rightHandSide = VectorXr::Zero(2*nnodes);
	this->_rightHandSide.topRows(nnodes) = rightHandData;

	VectorXr _forcingTerm;
	if(this->isSpaceVarying)
	{
		Assembler::forcingTerm(mesh_, fe, u, _forcingTerm);
	}

	// Debugging purpose
	//std::cout << "Preliminary -- apply() -- phase completed" << std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(Real lambda)
{
	//[[PART TO COMPUTE LAMBDA IS OUT NOW]]
	SpMat R1_lambda = (-lambda)*R1_;
	SpMat R0_lambda = (-lambda)*R0_;

	this->buildMatrixNoCov(psi_, R1_lambda, R0_lambda);

	UInt nnodes = mesh_.num_nodes();

	// Define forcing term
	if(this->isSpaceVarying)
	{
	    	_rightHandSide.bottomRows(nnodes)=lambda*_forcingTerm;
	}

	//Applying boundary conditions if necessary
	if(regressionData_.getDirichletIndices().size() != 0)  // if areal data NO BOUNDARY CONDITIONS
	{
		addDirichletBC();
	}

	// Factorize the system
	system_factorize();

	// Solve the system
    	_solution = this->template system_solve(this->_rightHandSide);

	//Debugging purpose
	//Rprintf("\nsolution\n");
	//for (int i = 0; i < _solution.size(); i++)
	//	Rprintf("%f ", _solution[i]);
	//Rprintf("\n");
}

//----------------------------------------------------------------------------//

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionData, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim> & mesh, const RegressionData & regressionData, const OptimizationData & optimizationData):
		MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>(mesh, regressionData, optimizationData){};

	void preapply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::preapply(stiff, ForcingTerm(std::vector<Real>(1)));
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataElliptic, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim> & mesh, const RegressionDataElliptic & regressionData, const OptimizationData & optimizationData):
		MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>(mesh, regressionData, optimizationData){};

	void preapply()
	{
		if(mydim!=2 || ndim !=2)
		{
			#ifdef R_VERSION_
				Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
			#else
				std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
			#endif
		}
		else
		{
			typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

			const Real & c = this->regressionData_.getC();
			const Eigen::Matrix<Real,2,2> & K = this->regressionData_.getK();
			const Eigen::Matrix<Real,2,1> & b = this->regressionData_.getBeta();

			MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>::preapply(c*mass+stiff[K]+dot(b,grad), ForcingTerm(std::vector<Real>(1)));
		}
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim> & mesh, const RegressionDataEllipticSpaceVarying & regressionData, const OptimizationData & optimizationData):
		MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, regressionData, optimizationData){};

	void preapply()
	{
		if(mydim!=2 || ndim !=2)
		{
			#ifdef R_VERSION_
				Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
			#else
				std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
			#endif
		}
		else
		{
			typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

			const Reaction & c = this->regressionData_.getC();
			const Diffusivity & K = this->regressionData_.getK();
			const Advection & b = this->regressionData_.getBeta();
			const ForcingTerm & u= this->regressionData_.getU();

			this->isSpaceVarying = TRUE;

			MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>::preapply(c*mass+stiff[K]+dot(b,grad), u);
		}
	}
};

#endif
