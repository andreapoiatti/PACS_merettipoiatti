#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

// [TO DEBUG] Probably fixed in new versions of code
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::addDirichletBC()
{
	UInt id1;
	UInt id3;
	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt> & bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real> & bc_values  = regressionData_.getDirichletValues();

	UInt nbc_indices = bc_indices.size();

	Real pen = 10e20;

	for (auto i=0; i<nbc_indices; i++)
	{
		id1 = bc_indices[i];
		id3 = id1+nnodes;

		A_.coeffRef(id1, id1) = pen;
		A_.coeffRef(id3, id3) = pen;


		_b(id1) += bc_values[i]*pen;
		_b(id3)  = 0;
	 }

	A_.makeCompressed();
}

//----------------------------------------------------------------------------//
// Setters

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setPsi()
{
	// Psi is a nnodes x nlocations matrix, first fetch the dimensions and resize the empty matrix
	UInt nnodes 	= mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();  // Total number of valid observations (non NAs)
	psi_.resize(nlocations, nnodes);

	// Different management strategies according to the presence of locations
	if (!hasLocations_)
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

		// Empty vector of triplets is resized as the number of valid (not NA) observations
		std::vector<coeff> tripletAll;
		auto k = regressionData_.getObservationsIndices();  // contains the number of the nodes linked with the valid observations
		tripletAll.reserve(k.size());

		for (int i = 0; i < k.size(); ++i)
		{
			// Add a value 1 for each valid index in row i
			// and column k[i] (under psi_k[i], the associated node)
			tripletAll.push_back(coeff(i, k[i], 1.0));
		}

		psi_.setFromTriplets(tripletAll.begin(), tripletAll.end());
		psi_.makeCompressed();  // Compress sparse matrix
	}
	else
	{
		// THEORETICAL REMARK:
		// If isLocationsByNodes is false locations are unspecified points
		// so we have to evaluate them directly

		constexpr UInt 			Nodes = (mydim==2)? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> 	tri_activated;
		Eigen::Matrix<Real, Nodes, 1> 	coefficients;
		Real 				evaluator;

		for (UInt i = 0; i < nlocations; i++)
		{
			// Search the element containing the point
			tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			if (tri_activated.getId() == Identifier::NVAL)
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
				for (UInt node = 0; node < Nodes; ++node)
				{
					// Define vector of all zeros but "node" component (necessary for function evaluate_point)
					coefficients = Eigen::Matrix<Real, Nodes, 1>::Zero();
					coefficients(node) = 1; //Activates only current base

					// Evaluate psi_node in the node and put it in the right place in the matrix according to global index
					evaluator = evaluate_point<Nodes, mydim, ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
					// Insert the value in the column of the global indexing of the evaluated node
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}

		psi_.makeCompressed(); // Compress sparse matrix
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setH()
{

	// Debug text
 	// std::cout << "Computing Projection Matrix" << std::endl;
 	UInt nlocations = regressionData_.getNumberofObservations();

	// Debug text
 	// regressionData_.printCovariates(std::cout);
 	MatrixXr W(this->regressionData_.getCovariates());

	// Debug text
 	//std::cout << "W "<< W <<std::endl;

 	if (!hasLocations_)
 	{
		// Some rows might be discarded [[we probably have data for every node not only the points ???]]
 		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
 		for (auto i = 0; i < nlocations; ++i)
 		{
 			auto index_i = regressionData_.getObservationsIndices()[i];
 			for (auto j = 0; j < W.cols(); ++j)
 			{
 				W_reduced(i,j) = W(index_i,j);
 			}
 		}
 		W = W_reduced;
 	}

 	MatrixXr WTW(W.transpose()*W);		// Compute W^t*W
 	H_ = W*WTW.ldlt().solve(W.transpose()); // Using cholesky LDLT decomposition for computing hat matrix
}

template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::setQ()
{
	// Debug text
 	// std::cout << "Computing Orthogonal Space Projection Matrix" << std::endl;

	// Remember Q = I-H
 	Q_.resize(H_.rows(), H_.cols());  // Resizing dimensions as H
 	Q_ = -H_;
 	for (int i = 0; i < H_.rows(); ++i)
 	{
 		Q_(i,i) += 1;  // Summing the identity by rows (or columns)
 	}
}

//----------------------------------------------------------------------------//
// Builders

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::buildA(const SpMat & Psi,  const SpMat & R1,  const SpMat & R0)
{
	// Most general A formulation [*]
	// A == [ Psi^t*Q*Psi	, -lambda*R1^t	] ==: [alpha, beta]
	//	[-lambda*R1 	, -lambda*R0	]     [beta, gamma]

	// IMPORTANT REMARK:
	// the parameters passed to the function have an ambiguous meaning:
	// R0 is actually gamma
	// R1 is actually beta

	UInt nnodes = mesh_.num_nodes();

	// Sparse matrix of the top-left block alpha [Q is NOT present]
	SpMat DMat = Psi.transpose()*Psi;

	// Vector to be filled with the triplets used to build A (reserved with the right dimension)
	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());

	// Parsing all matrices, reading the values to be put inside A, coordinates according to [*]
	for (int k = 0; k < DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), 	    it.col(),	     it.value()));
		}
	for (int k = 0; k < R0.outerSize(); ++k)
		for (SpMat::InnerIterator it(R0, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes, it.value()));
		}
	for (int k = 0; k < R1.outerSize(); ++k)
		for (SpMat::InnerIterator it(R1, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.col(), 	    it.row()+nnodes, it.value()));
		}
	for (int k = 0; k < R1.outerSize(); ++k)
	  	for (SpMat::InnerIterator it(R1, k); it; ++it)
	  	{
		  	tripletAll.push_back(coeff(it.row()+nnodes, it.col(), 	     it.value()));
	  	}

	// Define, resize, fill and compress A
	A_.setZero();
	A_.resize(2*nnodes,2*nnodes);
	A_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	A_.makeCompressed();

	// Debug message
	// std::cout << "Coefficients' Matrix Set Correctly" << std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::getDataMatrix(SpMat & DMat)
{
	// Top-left block definition
	Rprintf("Hello\n");

	UInt nnodes = mesh_.num_nodes();
	DMat.resize(nnodes, nnodes);


	// Check if Q == Identity

	if (!isRegression_)
		DMat = psi_.transpose()*psi_;
	else
	{
		Rprintf("Hello\n");
		DMat = (SpMat(psi_.transpose())*Q_*psi_).sparseView();
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::getDataMatrixByIndices(SpMat & DMat)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	DMat.resize(nnodes, nnodes);

	if (!isRegression_)
	{
		DMat.reserve(1);
		for (auto i = 0; i < nlocations; ++i)
		{
			auto index = regressionData_.getObservationsIndices()[i];
			DMat.insert(index, index) = 1;
		}
	}
	else
	{
		//May be inefficient
		for (auto i = 0; i < nlocations; ++i)
		{
			auto index_i = regressionData_.getObservationsIndices()[i];
			for (auto j = 0; j < nlocations; ++j)
			{
				auto index_j = regressionData_.getObservationsIndices()[j];
				DMat.insert(index_i, index_j) = Q_(i, j);
			}
		}
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat & DMat,  const SpMat & AMat,  const SpMat & MMat)
{
 	// I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
 	// _coeffmatrix.setFromTriplets(tripletA.begin(), tripletA.end());

	// Builder of the coefficients matrix
	//	[ DMat | AMat^t ]
	//	[ AMat | MMat   ]

 	UInt nnodes = mesh_.num_nodes();

	// Vector to be filled with the triplets used to build _coeffmatrix (reserved with the right dimension)
 	std::vector<coeff> tripletAll;
 	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	// Parsing all matrices, reading the values to be put inside _coeffmatrix, coordinates according to [*]
 	for (int k = 0; k < DMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(DMat, k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.row(), 		it.col(),		it.value()));
 	  }
 	for (int k = 0; k < MMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(MMat, k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.row()+nnodes, 	it.col()+nnodes,	it.value()));
 	  }
 	for (int k = 0; k < AMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(AMat, k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.col(), 		it.row()+nnodes,	it.value()));
 	  }
 	for (int k = 0; k < AMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(AMat, k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.row()+nnodes, 	it.col(),		it.value()));
 	  }

	// Define, resize, fill and compress _coeffmatrix
 	_coeffmatrix.setZero();
 	_coeffmatrix.resize(2*nnodes, 2*nnodes);
 	_coeffmatrix.setFromTriplets(tripletAll.begin(), tripletAll.end());
 	_coeffmatrix.makeCompressed();
 	//std::cout << "Coefficients' Matrix Set Correctly" << std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER,mydim,ndim>::getRightHandData(VectorXr & rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if (!isRegression_)
	{
		if (!hasLocations_)
		{
			for (auto i = 0; i < nlocations; ++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData = psi_.transpose()*regressionData_.getObservations();
		}
	}
	else
	{
		rightHandData = psi_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
}

//----------------------------------------------------------------------------//
// Utilities

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
MatrixXr MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::LeftMultiplybyQ(const MatrixXr & u)
{
	if (!isRegression_)
	{
		// Q is the projecton on Col(W) if W == 0 => Q = Identity
		return u;
	}
	else
	{
		// [[ Non abbimo H già in memoria, perhè non usarla al posto di ricalcolarla ogni volta???]]
		MatrixXr W(this->regressionData_.getCovariates());
		// Check factorization, if not present factorize the matrix W^t*W
		if (isWTWfactorized_ == false )
		{
			WTWinv_.compute(W.transpose()*W);
			isWTWfactorized_ = true;
		}
		// Compute H (or I-Q) the projection on Col(W) and multiply it times u
		MatrixXr Pu = W*WTWinv_.solve(W.transpose()*u); // P*u = W*(W^t*W)^{-1}*W^t*u
		return (u - Pu); // note (I-P)*u == Q*u
	}
}
//----------------------------------------------------------------------------//
// Factorizer

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::system_factorize()
{
	UInt nnodes = mesh_.num_nodes();

	// First phase: Factorization of matrix A
	Adec_.compute(A_);

	// If Q!=Identity
	if (isRegression_)
	{
		// Second phase: factorization of matrix  G =  C + [V * A^-1 * U]

		// Definition of matrices
		// U = [ psi * W | 0 ]^t
		MatrixXr W(this->regressionData_.getCovariates());
		U_ = MatrixXr::Zero(2*nnodes, W.cols());
		U_.topRows(nnodes) = psi_.transpose()*W;

		// D = U^T * A^-1 * U
		MatrixXr D = U_.transpose()*Adec_.solve(U_);

		// G = C + D
		MatrixXr G = -W.transpose()*W + D;
		Gdec_.compute(G);  // Factorize matix G
	}
}

//----------------------------------------------------------------------------//
// Solver

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> & b)
{
	// Resolution of the system A * x1 = b
	MatrixXr x1 = Adec_.solve(b);

	if (isRegression_)
	{
		// Resolution of G * x2 = U^T * x1
		MatrixXr x2 = Gdec_.solve(U_.transpose()*x1);
		// Resolution of the system A * x3 = U * x2
		x1 -= Adec_.solve(U_*x2);
	}

	return x1;
}

//----------------------------------------------------------------------------//
// Main function
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper)
{
	UInt nnodes = mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	if (regressionData_.getCovariates().rows() == 0)
	{
		isRegression_ = false;
		Rprintf("No Regression\n");
	}
	else
	{
		isRegression_ = true;
		Rprintf("Regression\n");
	}

	if (regressionData_.isLocationsByNodes() == true)
	{
		hasLocations_ = false;
		Rprintf("Locations by nodes\n");
	}
	else
	{
		hasLocations_ = true;
		Rprintf("Locations are not nodes\n");
	}

	// Set matrix Psi
	setPsi();

	// If there are covariates in the model set H and Q
	if (isRegression_)
 	{
 		setH();
 		setQ();
 	}

	// Set the data matrix according to the nature of Psi
	if (hasLocations_)
	{
		Rprintf("Hello\n");
		getDataMatrix(DMat_);
	}
	else
	{
		Rprintf("Hello man\n");
		getDataMatrixByIndices(DMat_);
	}

	//Assemble matrices R0-R1
	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);

	// Define Right hand data
	VectorXr rightHandData;
	getRightHandData(rightHandData);

	// Define vector b
	this->_b = VectorXr::Zero(2*nnodes);
	this->_b.topRows(nnodes) = rightHandData;

	// Fix matrices from a certain lambda0
	if (regressionData_.getLambda().size()!=0)
		lambda_ = regressionData_.getLambda()[0];
	else
		lambda_ = 0;
	Rprintf("Done :)");
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply_2()
{
	//[[PART TO COMPUTE LAMBDA IS OUT NOW]]

	// Solver
	SpMat R1_lambda = (-lambda_) * R1_;
	SpMat R0_lambda = (-lambda_) * R0_;

	// Build A and coeffmatrix [[perchè due versioni???]]
	this->buildA(		psi_,  R1_lambda, R0_lambda);
	this->buildCoeffMatrix(	DMat_, R1_lambda, R0_lambda);

	// Debug text
	// std::cout << coeffmatrix_ << std::endl;

	// Applying boundary conditions if necessary
	if( regressionData_.getDirichletIndices().size() != 0)
		addDirichletBC();

	// Factorize the system
	system_factorize();

	// Solve the system
	_solution = this->template system_solve(this->_b);
}

//----------------------------------------------------------------------------//

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionData, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionData& regressionData):MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff;
		Stiff EStiff;
		ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::apply(stiff);
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataElliptic, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataElliptic& regressionData):MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = this->regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& beta = this->regressionData_.getBeta();

	    MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = this->regressionData_.getC();
		const Diffusivity& K = this->regressionData_.getK();
		const Advection& beta = this->regressionData_.getBeta();

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
	}
};

#endif
