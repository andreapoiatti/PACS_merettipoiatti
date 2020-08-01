#ifndef __MIXEDFEREGRESSION_IMP_H__
#define __MIXEDFEREGRESSION_IMP_H__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

//----------------------------------------------------------------------------//
// Dirichlet BC

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::addDirichletBC() //adds boundary conditions to all
{
	UInt id1,id3;

	UInt nnodes = N_*M_;

	const std::vector<UInt> * bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real> * bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices->size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=(*bc_indices)[i];
			id3=id1+nnodes;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;


			_rightHandSide(id1)=(*bc_values)[i]*pen;
			_rightHandSide(id3)=0;
	 }

	matrixNoCov_.makeCompressed();
}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::addDirichletBC_matrix() //adds boundary conditions to all
{
	UInt id1,id3;

	UInt nnodes = N_*M_;

	const std::vector<UInt> * bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real> * bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices->size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=(*bc_indices)[i];
			id3=id1+nnodes;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;
	 }

	matrixNoCov_.makeCompressed();
}




//Add NA
template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::addNA()
{

	const std::vector<UInt> * observations_na= regressionData_.getObservationsNA();

	for(UInt id:*observations_na)
	{
		for(UInt j=0; j<psi_.cols(); ++j)
		{
			if(psi_.coeff(id,j)!=0)
				psi_.coeffRef(id, j) = 0;
		}
	}
	psi_.pruned();
	psi_.makeCompressed();
}


//----------------------------------------------------------------------------//
// Setters

template<typename InputHandler>
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler>::setPsi(const MeshHandler<ORDER, mydim, ndim> & mesh_)
{
	// Psi is a nlocations x nnodes  matrix, first fetch the dimensions
	UInt nnodes = mesh_.num_nodes();	// Define the number of nodes
	// Set the number of locations depending on presence or not of temporal data
	UInt nlocations = regressionData_.isSpaceTime() ? regressionData_.getNumberofSpaceObservations() : regressionData_.getNumberofObservations();
	psi_.resize(nlocations, nnodes);	// Resize the matrix

	// Optimized strategies according to the presence of locations
	if(regressionData_.isLocationsByNodes() & !regressionData_.isLocationsByBarycenter()) // .Pointwise data -- no barycenters
	{
		std::vector<coeff> tripletAll;
		if(!regressionData_.isSpaceTime()) // Just spatial case
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
			const std::vector<UInt> * k = regressionData_.getObservationsIndices();
			UInt k_size = k->size();
			tripletAll.reserve(k_size);

			for (UInt i=0; i< k_size; ++i)
			{
				// Add a value 1 for each valid index in row i
				// and column k[i] (under psi_k[i], the associated node)
				tripletAll.push_back(coeff(i,(*k)[i],1.0));
			}
		}
		else
		{
			tripletAll.reserve(nlocations);
			for (UInt i=0; i<nlocations; ++i)
			{
				tripletAll.push_back(coeff(i,i,1.0));
			}
		}
		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	}
	else if (regressionData_.isLocationsByBarycenter() && (regressionData_.getNumberOfRegions() == 0)) //pointwise data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.getElement(regressionData_.getElementId(i));

			if(tri_activated.getId() == Identifier::NVAL)
			{
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			}else //tri_activated.getId() found
			{
				for(UInt node = 0; node < Nodes ; ++node)
				{
					evaluator = regressionData_.getBarycenter(i,node);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop
	}
	else if ((!regressionData_.isLocationsByBarycenter()) && (regressionData_.getNumberOfRegions() == 0))
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,Nodes,1> coefficients;

		Real evaluator;
		this->barycenters_.resize(nlocations, Nodes);
		this->element_ids_.resize(nlocations);
		for(UInt i=0; i<nlocations;i++)
		{
			if (regressionData_.getSearch() == 1) { //use Naive search
				tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			} else if (regressionData_.getSearch() == 2) { //use Tree search (default)
				tri_activated = mesh_.findLocationTree(regressionData_.getLocations()[i]);
			}

			if(tri_activated.getId() == Identifier::NVAL)
			{
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			}else //tri_activated.getId() found
			{
				element_ids_(i)=tri_activated.getId();
				for(UInt node = 0; node < Nodes ; ++node)
				{
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<Nodes,mydim,ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
					barycenters_(i,node)=tri_activated.getBaryCoordinates(regressionData_.getLocations()[i])[node];
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		} //end of for loop
	}
	else //areal data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;

		Real *tab; //Psi_i
		tab = (Real*) malloc(sizeof(Real)*nnodes);
		for(UInt i=0; i<nlocations;i++) //nlocations = number of regions
		{
			for (UInt k=0; k<nnodes; k++) {tab[k]=0;}
			for (UInt j=0; j<mesh_.num_elements(); j++)
			{
				if ((*(regressionData_.getIncidenceMatrix()))(i,j) == 1) //element j is in region i
				{
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j);
					for (UInt k=0; k<Nodes; k++)
					{
						tab[tri[k].getId()] += integratePsi(tri,k); // integral over tri of psi_k
					}
				}
			}
			for (int k=0; k<nnodes; k++)
			{
				if (tab[k] != 0)
				{
					psi_.insert(i,k) = tab[k]/A_(i); //divide by |D_i|
				}
			}
		}
		free(tab);
	}

	psi_.makeCompressed();

	// Additional storage of the transposed
	psi_t_ = SpMat(psi_.transpose());
	psi_t_.makeCompressed(); // Compress sparse matrix

}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setH(void)
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

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setQ(void)
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



template<typename InputHandler>
template<UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler>::setA(const MeshHandler<ORDER, mydim, ndim> & mesh_)
{
	UInt nRegions = regressionData_.getNumberOfRegions();	// Number of regions for areal partition
	// If the problem is temporal, m stores the number of temporal nodes, else is defaulted as 1
	UInt m = regressionData_.isSpaceTime() ? regressionData_.getNumberofTimeObservations():1;
	if(!this->regressionData_.isArealDataAvg())
	{ //areal data for FPIRLS
		A_ = VectorXr::Ones(m*nRegions);	// vector of pure ones
	}
	else
	{
		A_ = VectorXr::Zero(m*nRegions);	// neutral vector to be filled
		for(UInt i=0; i<nRegions; i++)		// fill the vector
		{
			const MatrixXi * imp = regressionData_.getIncidenceMatrix();
			for(UInt j=0; j<imp->cols(); j++)
			{
				if((*imp)(i,j) == 1)	// Valid input for region
				{
					A_(i) += mesh_.elementMeasure(j); // Add area
				}
			}
			for(UInt k=1; k<m; k++) // if m=1 we avoid the step
			{
				A_(i+k*nRegions) = A_(i); // Replicate the vector m times
			}
		}
	}
}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::setDMat(void)
{
        //Additional storage of the transpose, which changes in case of spacetime
	psi_t_ = SpMat(psi_.transpose());
	psi_t_.makeCompressed(); // Compress sparse matrix

	if(regressionData_.getWeightsMatrix()->size() == 0) // no weights
		DMat_ = psi_;
	else
		DMat_ = regressionData_.getWeightsMatrix()->asDiagonal()*psi_;


	if(regressionData_.getNumberOfRegions() == 0) // pointwise data
		DMat_ = psi_t_*DMat_;
	else                                        // areal data: need to add the diag(|D_1|,...,|D_N|)
		DMat_ = psi_t_*A_.asDiagonal()*DMat_;
}




//----------------------------------------------------------------------------//
// Utilities [[NOT VERY OPTMIZED, SENSE??, we have Q and P...]]

template<typename InputHandler>
MatrixXr MixedFERegressionBase<InputHandler>::LeftMultiplybyQ(const MatrixXr& u)
{
	const VectorXr* P = this->regressionData_.getWeightsMatrix();

	if (regressionData_.getCovariates()->rows() == 0){
		if(P->size()==0)
			return u;
		else
			return P->asDiagonal()*u;
	}
	else{
		MatrixXr W(*(this->regressionData_.getCovariates()));
		if (isWTWfactorized_ == false ){
			if(P->size()==0){
				WTW_.compute(W.transpose()*W);
			}else{
				WTW_.compute(W.transpose()*P->asDiagonal()*W);
			}
			isWTWfactorized_=true;
		}

		MatrixXr Pu;
		if(P->size()==0){
			Pu = W*WTW_.solve(W.transpose()*u);
		}else{
			Pu= W*WTW_.solve(W.transpose()*P->asDiagonal()*u);
		}
		if(P->size()==0)
			return u-Pu;
		else
			return P->asDiagonal()*(u - Pu);
	}

}

//----------------------------------------------------------------------------//
// Builders


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

        const VectorXr* obsp=regressionData_.getObservations();


	if (regressionData_.getCovariates()->rows() == 0) //no covariate
	{
		if (regressionData_.isLocationsByNodes() && !regressionData_.isSpaceTime())
		{
				VectorXr tmp = LeftMultiplybyQ(*obsp);
				for (auto i=0; i<nlocations;++i)
				{
					auto index_i = (*(regressionData_.getObservationsIndices()))[i];
					rightHandData(index_i) = tmp(i);
				}
		}
		else if (regressionData_.isLocationsByNodes() && regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
		{
			for (auto i=0; i<regressionData_.getObservationsIndices()->size();++i)
			{
				auto index_i = (*(regressionData_.getObservationsIndices()))[i];
				rightHandData(index_i) = (*obsp)[index_i];
			}
		}
		else if (regressionData_.getNumberOfRegions() == 0) //pointwise data
		{
			rightHandData=psi_.transpose()*LeftMultiplybyQ(*obsp);
		}
		else //areal data
		{
			rightHandData=psi_.transpose()*A_.asDiagonal()*LeftMultiplybyQ(*obsp);
		}
	}
	else if (regressionData_.getNumberOfRegions() == 0) //with covariates, pointwise data
	{
		rightHandData=psi_.transpose()*LeftMultiplybyQ(*obsp);
	}
	else //with covariates, areal data
	{
		rightHandData=psi_.transpose()*A_.asDiagonal()*LeftMultiplybyQ(*obsp);
	}
}


template<typename InputHandler>
template<typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedFERegressionBase<InputHandler>::buildSpaceTimeMatrices()
{
	SpMat IM(M_,M_);
	SpMat phi;

	if(regressionData_.getFlagParabolic()) // Parabolic case
	{
		MixedFDRegression <InputHandler> FiniteDifference(mesh_time_,regressionData_);
		FiniteDifference.setDerOperator();
		SpMat L = FiniteDifference.getDerOpL(); // Matrix of finite differences
		IM.setIdentity();
		LR0k_ = kroneckerProduct(L,R0_);
		phi = IM;
		//! right hand side correction for the initial condition:
		rhs_ic_correction_ = (1/(mesh_time_[1]-mesh_time_[0]))*(R0_*(*(regressionData_.getInitialValues())));
	}
	else	// Separable case
	{
		MixedSplineRegression <InputHandler, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE> Spline(mesh_time_,regressionData_);
		SpMat IN(N_,N_);
		Spline.setPhi();
		Spline.setTimeMass();
		Spline.smoothSecondDerivative();
		if(regressionData_.getFlagMass()) // Mass penalization
		{
			IM = Spline.getTimeMass();
			IN = R0_;
		}
		else	// Identity penalization
		{
			IM.setIdentity();
			IN.setIdentity();
		}
		phi = Spline.getPhi();
		SpMat Pt = Spline.getPt();
		Ptk_ = kroneckerProduct(Pt,IN);
	}
	// Make the Kronecker product to tensorize the system
	SpMat psi_temp =  psi_;
	SpMat R1_temp = R1_;
	SpMat R0_temp = R0_;
	psi_.resize(N_*M_,N_*M_);
	psi_ = kroneckerProduct(phi,psi_temp);
	addNA();
	R1_.resize(N_*M_,N_*M_);
	R1_ = kroneckerProduct(IM,R1_temp);
	R1_.makeCompressed();
	R0_.resize(N_*M_,N_*M_);
	R0_ = kroneckerProduct(IM,R0_temp);
	R0_.makeCompressed();

	//! right hand side correction for the forcing term:

	if(this->isSpaceVarying)
	{
		VectorXr forcingTerm = rhs_ft_correction_;
		rhs_ft_correction_.resize(M_*N_);
		for(UInt i=0; i<N_; i++)
		{
			for(UInt j=0; j<M_; j++)
			{
				rhs_ft_correction_(i+j*N_) = forcingTerm(i);
			}
		}
	}
}



template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildMatrixNoCov(const SpMat& R1,  const SpMat& R0)
{

	UInt nnodes = N_*M_;



	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat_.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());


	for (int k=0; k<DMat_.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat_,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for (int k=0; k<R0.outerSize(); ++k)
		for (SpMat::InnerIterator it(R0,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	matrixNoCov_.setZero();
	matrixNoCov_.resize(2*nnodes,2*nnodes);
	matrixNoCov_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	matrixNoCov_.makeCompressed();
}


//----------------------------------------------------------------------------//
// Main functions



template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::system_factorize() {

	UInt nnodes = N_*M_;

	const VectorXr* P = regressionData_.getWeightsMatrix(); // matrix of weights

	// First phase: Factorization of matrixNoCov
	matrixNoCovdec_->compute(matrixNoCov_);


	if (regressionData_.getCovariates()->rows() != 0) {
		// Second phase: factorization of matrix  G =  C + [V * matrixNoCov^-1 * U]= C + D

		// Definition of matrix U = [ psi^T * A * W | 0 ]^T and V= [ W^T*psi| 0]

		MatrixXr W(*(this->regressionData_.getCovariates()));

		U_ = MatrixXr::Zero(2*nnodes, W.cols());
		V_ = MatrixXr::Zero(W.cols(),2*nnodes);

		if(P->size()==0){
			V_.leftCols(nnodes) = W.transpose()*psi_;
		}else{
			V_.leftCols(nnodes) = W.transpose()*P->asDiagonal()*psi_;
		}
		// build "right side" of U_
		if(P->size() == 0)
			U_.topRows(nnodes) = W;
		else
			U_.topRows(nnodes) = P->asDiagonal()*W;

		// build "left side" of U_
		if(regressionData_.getNumberOfRegions()==0){ // pointwise data
		  U_.topRows(nnodes) = psi_.transpose()*U_.topRows(nnodes);
		}
		else{                                          //areal data
		  U_.topRows(nnodes) = psi_.transpose()*A_.asDiagonal()*U_.topRows(nnodes);
    	}
		MatrixXr D = V_*matrixNoCovdec_->solve(U_);
		// G = C + D
		MatrixXr G;
		if(P->size()==0){
			G = -W.transpose()*W + D;
		}else{
			G = -W.transpose()*P->asDiagonal()*W + D;
		}
		Gdec_->compute(G);

	}
}

template<typename InputHandler>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler>::system_solve(const Eigen::MatrixBase<Derived> &b)
{

	// Resolution of the system matrixNoCov * x1 = b
	MatrixXr x1 = matrixNoCovdec_->solve(b);
	if (regressionData_.getCovariates()->rows() != 0) {
		// Resolution of G * x2 = V * x1
		MatrixXr x2 = Gdec_->solve(V_*x1);
		// Resolution of the system matrixNoCov * x3 = U * x2
		x1 -= matrixNoCovdec_->solve(U_*x2);
	}
	return x1;
}



template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	int GCVmethod = regressionData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
	}
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeGeneralizedCrossValidation(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	VectorXr dataHat;
	const VectorXr* z = regressionData_.getObservations();
	if(regressionData_.getCovariates()->rows()==0) //Data estimated from the model
		dataHat = psi_*_solution(output_indexS,output_indexT).topRows(psi_.cols());
	else
		dataHat = *z - LeftMultiplybyQ(*z) + LeftMultiplybyQ(psi_*_solution(output_indexS,output_indexT).topRows(psi_.cols()));
	UInt n = dataHat.rows();
	if(regressionData_.isSpaceTime())
		{
			const std::vector<UInt> * observations_na= regressionData_.getObservationsNA();
			for(UInt id:*observations_na)
			{
				dataHat[id]=0;
			}
			n-=observations_na->size();
		}
//! GCV computation
	_GCV(output_indexS,output_indexT) = (n / ((n - regressionData_.getTuneParam()*_dof(output_indexS, output_indexT)) * (n - regressionData_.getTuneParam()*_dof(output_indexS, output_indexT)))) * (*z-dataHat).dot(*z-dataHat);
	if (_GCV(output_indexS,output_indexT) < _bestGCV)
	{
		bestLambdaS_ = output_indexS;
		bestLambdaT_ = output_indexT;
		_bestGCV = _GCV(output_indexS,output_indexT);
	}
}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	std::string file_name;
	UInt nnodes = N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees=0;


	MatrixXr X1;
	if (regressionData_.getNumberOfRegions() == 0){ //pointwise data
		X1 = psi_.transpose() * LeftMultiplybyQ(psi_);
	}else{ //areal data
		X1 = psi_.transpose() * A_.asDiagonal() * LeftMultiplybyQ(psi_);
	}


	if (isRcomputed_ == false)
	{
		isRcomputed_ = true;
		//take R0 from the final matrix since it has already applied the dirichlet boundary conditions
		SpMat R0 = matrixNoCov_.bottomRightCorner(nnodes,nnodes)/lambdaS;

		R0dec_.compute(R0);
		if(!regressionData_.isSpaceTime() || !regressionData_.getFlagParabolic())
		{
			MatrixXr X2 = R0dec_.solve(R1_);
			R_ = R1_.transpose() * X2;
		}
	}

	MatrixXr P;
	MatrixXr X3=X1;

	//define the penalization matrix: note that for separable smoothin should be P=lambdaS*Psk+lambdaT*Ptk
	// but the second term has been added to X1 for dirichlet boundary conditions
	if (regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
	{
		SpMat X2 = R1_+lambdaT*LR0k_;
		P = lambdaS*X2.transpose()*R0dec_.solve(X2);
	}
	else
	{
		P = lambdaS*R_;
	}


	if(regressionData_.isSpaceTime() && !regressionData_.getFlagParabolic())
		X3 += lambdaT*Ptk_;

	//impose dirichlet boundary conditions if needed
	if(regressionData_.getDirichletIndices()->size()!=0)
	{
		const std::vector<UInt> * bc_indices = regressionData_.getDirichletIndices();
		UInt nbc_indices = bc_indices->size();

		Real pen=10e20;
		for(UInt i=0; i<nbc_indices; i++)
		{
			UInt id = (*bc_indices)[i];
			X3(id,id)=pen;
		}
	}

	X3 -= P;
	Eigen::PartialPivLU<MatrixXr> Dsolver(X3);

	const auto k = regressionData_.getObservationsIndices();

	if(!regressionData_.isSpaceTime() && regressionData_.isLocationsByNodes()) {
		if(regressionData_.getCovariates()->rows() != 0)
			degrees += regressionData_.getCovariates()->cols();

		// Setup rhs B
		MatrixXr B;
		B = MatrixXr::Zero(nnodes,nlocations);
		// B = I(:,k) * Q
		for (auto i=0; i<nlocations;++i) {
			VectorXr ei = VectorXr::Zero(nlocations);
			ei(i) = 1;
			VectorXr Qi = LeftMultiplybyQ(ei);
			for (int j=0; j<nlocations; ++j) {
				B((*k)[i], j) = Qi(j);
			}
		}
		// Solve the system TX = B
		MatrixXr X;
		X = Dsolver.solve(B);
		// Compute trace(X(k,:))
		for (int i = 0; i < k->size(); ++i) {
			degrees += X((*k)[i], i);
		}
	}

	if (regressionData_.isSpaceTime() || !regressionData_.isLocationsByNodes())
	{
		MatrixXr X;
		X = Dsolver.solve(X1);

		if (regressionData_.getCovariates()->rows() != 0) {
			degrees += regressionData_.getCovariates()->cols();
		}
		for (int i = 0; i<nnodes; ++i) {
			degrees += X(i,i);
		}
	}

	_dof(output_indexS,output_indexT) = degrees;
}

template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{

	UInt nnodes = N_*M_;
	UInt nlocations = regressionData_.getNumberofObservations();

	// std::random_device rd;
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = regressionData_.getNrealizations();
	MatrixXr u(nlocations, nrealizations);
	for (int j=0; j<nrealizations; ++j) {
		for (int i=0; i<nlocations; ++i) {
			if (distribution(generator)) {
				u(i,j) = 1.0;
			}
			else {
				u(i,j) = -1.0;
			}
		}
	}

	// Define the first right hand side : | I  0 |^T * psi^T * A * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	if (regressionData_.getNumberOfRegions() == 0){
		b.topRows(nnodes) = psi_.transpose() * LeftMultiplybyQ(u);
	}else{
		b.topRows(nnodes) = psi_.transpose() * A_.asDiagonal() * LeftMultiplybyQ(u);
	}

	// Resolution of the system
	MatrixXr x = system_solve(b);

	MatrixXr uTpsi = u.transpose()*psi_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (regressionData_.getCovariates()->rows() != 0) {
		q = regressionData_.getCovariates()->cols();
	}
	// For any realization we compute the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {

		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof(output_indexS,output_indexT) = mean;
}




template<typename InputHandler>
template<UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, typename A>
void MixedFERegressionBase<InputHandler>::preapply(EOExpr<A> oper, const ForcingTerm & u, const MeshHandler<ORDER, mydim, ndim> & mesh_)
{
	const MatrixXr * Wp = regressionData_.getCovariates();

	UInt nnodes = N_*M_;	// total number of spatio-temporal nodes
	FiniteElement<IntegratorSpace, ORDER, mydim, ndim> fe;

	// Set Areal data if present and no already done
	if(regressionData_.getNumberOfRegions()>0 && !isAComputed)
	{
		this->template setA<ORDER, mydim, ndim>(mesh_);
		isAComputed = true;
	}

	// Set psi matrix if not already done
	if(!isPsiComputed){
		this->template setPsi<ORDER, mydim, ndim>(mesh_);
		isPsiComputed = true;
	}

	// If there are covariates in the model set H and Q
	if(Wp->rows() != 0)
	{
		setH();
		setQ();
		std::cout<<"H and Q set"<<std::endl;
	}

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	if(!isR1Computed){
		Assembler::operKernel(oper, mesh_, fe, R1_);
		isR1Computed = true;
	}
	if(!isR0Computed){
		Assembler::operKernel(mass, mesh_, fe, R0_);
		isR0Computed = true;
	}

	if(this->isSpaceVarying)
	{
        //u=ForcingTerm(VectorXr::Zero(nnodes));
		Assembler::forcingTerm(mesh_, fe, u, rhs_ft_correction_);

	}

	if (regressionData_.isSpaceTime())
	{
		this->template buildSpaceTimeMatrices<IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>();
	}

	setDMat();	// Set matrix DMat for all cases

	// Debugging purpose
	Rprintf("Preliminary problem matrices building phase completed\n");
        }





//----------------------------------------------------------------------------//
// Composed operations

//to be general, it takes in input the value of lambda
template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildSystemMatrix(Real lambda)
{
        this->R1_lambda = (-lambda)*(R1_);
        this->R0_lambda = (-lambda)*(R0_);

        this->buildMatrixNoCov(this->R1_lambda, this->R0_lambda);
}


template<typename InputHandler>
void MixedFERegressionBase<InputHandler>::buildSystemMatrix(Real lambdaS, Real lambdaT)
{
	this-> R0_lambda = (-lambdaS)*R0_; // build the SouthEast block of the matrix
	this->R1_lambda = (-lambdaS)*R1_;
	if(regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
		R1_lambda -= lambdaS*(lambdaT*LR0k_); // build the SouthWest block of the matrix (also the NorthEast block transposed)

	//SpMat NWblock=DMat_;

	// build right side of NWblock
	//if(regressionData_.getWeightsMatrix()->size() == 0) // no weights
	//	NWblock = psi_;
	//else
	//	NWblock = regressionData_.getWeightsMatrix()->asDiagonal()*psi_;

	// build left side of NWblock
	//if(regressionData_.getNumberOfRegions()==0) // pointwise data
	    //NWblock=psi_.transpose()*NWblock;
	//else                                        // areal data: need to add the diag(|D_1|,...,|D_N|)
	  //  NWblock=psi_.transpose()*A_.asDiagonal()*NWblock;

	if(regressionData_.isSpaceTime() && !regressionData_.getFlagParabolic())
		DMat_+=lambdaT*Ptk_;


	this->buildMatrixNoCov(R1_lambda, R0_lambda);

	//Build again the original DMat_ for the next lambdaT
	if(regressionData_.isSpaceTime() && !regressionData_.getFlagParabolic())
		DMat_-=lambdaT*Ptk_;
}





//----------------------------------------------------------------------------//
// Public solvers
template<typename InputHandler>
MatrixXr MixedFERegressionBase<InputHandler>::apply_to_b(const MatrixXr & b)
{
        if(lambda_ != last_lambda)
        {
                this->buildSystemMatrix(lambda_);
                this->systemFactorize();
        }

	//Applying boundary conditions if necessary
	//DA sistemare
	if(lambda_ != last_lambda && regressionData_.getDirichletIndices()->size() > 0)  // if areal data NO BOUNDARY CONDITIONS
	{
		this->addDirichletBC_matrix();
	}

        last_lambda = lambda_;

        return this->template system_solve(b);
}


template<typename InputHandler>
MatrixXv  MixedFERegressionBase<InputHandler>::apply(void)

	{
if (isGAMData||regressionData_.isSpaceTime()||lambda_!=last_lambda)
{
	UInt nnodes = N_*M_;
	VectorXr rightHandData;
	getRightHandData(rightHandData); //updated
	this->_rightHandSide = VectorXr::Zero(2*nnodes);
	this->_rightHandSide.topRows(nnodes)=rightHandData;

        UInt sizeLambdaS=regressionData_.getLambdaS()->size();
	UInt sizeLambdaT=regressionData_.getLambdaT()->size();

	this->_solution.resize(sizeLambdaS,sizeLambdaT);
	this->_dof.resize(sizeLambdaS,sizeLambdaT);
	this->_GCV.resize(sizeLambdaS,sizeLambdaT);
	if(regressionData_.getCovariates()->rows()!=0)
	{
		this->_beta.resize(sizeLambdaS,sizeLambdaT);
	}

	VectorXr rhs= _rightHandSide;

       const VectorXr* obsp=regressionData_.getObservations();
        //std::cout<<"Size lambdaT is "<<sizeLambdaT<<std::endl;
	for(UInt s = 0; s<sizeLambdaS; ++s)
	{
		for(UInt t = 0; t<sizeLambdaT; ++t)
		{

			Real lambdaS = (*(regressionData_.getLambdaS()))[s];
			Real lambdaT = (*(regressionData_.getLambdaT()))[t];
			_rightHandSide=rhs;
                        if (!regressionData_.isSpaceTime())
			    { //DEBUGGING, poi si elimina e si setta tramite carrier
		             set_lambda(lambdaS);
			     buildSystemMatrix(lambda_);}
			else
                             buildSystemMatrix(lambdaS, lambdaT);

			//! right-hand side correction for space varying PDEs
			if(this->isSpaceVarying)
			{
			    _rightHandSide.bottomRows(nnodes)= (-lambdaS)*rhs_ft_correction_;
			}

			//! righ-hand side correction for initial condition in parabolic case
			if(regressionData_.isSpaceTime() && regressionData_.getFlagParabolic())
			{
				for(UInt i = 0; i<regressionData_.getInitialValues()->rows(); i++)
				{
					_rightHandSide(nnodes+i) -= lambdaS*rhs_ic_correction_(i);
				}
			}
			//Applying boundary conditions if necessary
			if(regressionData_.getDirichletIndices()->size() != 0)  // if areal data NO BOUNDARY CONDITIONS
				addDirichletBC();
			//factorization of the system for woodbury decomposition
			system_factorize();

			// system solution
			_solution(s,t) = this->template system_solve(this->_rightHandSide);


			if(regressionData_.computeGCV()&&(isGAMData||regressionData_.isSpaceTime()))
			{
				if (regressionData_.computeDOF())
					computeDegreesOfFreedom(s,t,lambdaS,lambdaT);
				computeGeneralizedCrossValidation(s,t,lambdaS,lambdaT);
			}
			else
			{
				_dof(s,t) = -1;
				_GCV(s,t) = -1;
			}

			// covariates computation
			if(regressionData_.getCovariates()->rows()!=0&&(isGAMData||regressionData_.isSpaceTime()))
			{
				MatrixXr W(*(this->regressionData_.getCovariates()));
				VectorXr P(*(this->regressionData_.getWeightsMatrix()));
				VectorXr beta_rhs;
				if( P.size() !=0){
					beta_rhs = W.transpose()*P.asDiagonal()*(*obsp - psi_*_solution(s,t).topRows(psi_.cols()));
				}else{
					beta_rhs = W.transpose()*(*obsp - psi_*_solution(s,t).topRows(psi_.cols()));
				}
				_beta(s,t) = WTW_.solve(beta_rhs);
			}
		}
	}
if (!(isGAMData||regressionData_.isSpaceTime())&&lambda_!=last_lambda)
    this->last_lambda=lambda_;
  }
//std::cout<<_solution(0,0).size()<<std::endl; per il caso GCV semplice la solution è il valore in posizione 0,0
return this->_solution;
}


template<>
class MixedFERegression<RegressionData> : public MixedFERegressionBase<RegressionData>
{
public:
	MixedFERegression(const RegressionData& regressionData, const UInt& nnodes_):MixedFERegressionBase<RegressionData>(regressionData, nnodes_){};
	MixedFERegression(const std::vector<Real>& mesh_time, const RegressionData& regressionData, const UInt& nnodes_, const UInt& spline_degree):MixedFERegressionBase<RegressionData>(mesh_time, regressionData, nnodes_, spline_degree){};

	template< UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
	void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData>::preapply<ORDER,mydim,ndim, IntegratorSpace, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>(stiff, ForcingTerm(std::vector<Real>(1)), mesh);
	}
};

template<>
class MixedFERegression<RegressionDataElliptic> : public MixedFERegressionBase<RegressionDataElliptic>
{
public:
	MixedFERegression(const RegressionDataElliptic& regressionData, const UInt& nnodes_):MixedFERegressionBase<RegressionDataElliptic>(regressionData, nnodes_){};
	MixedFERegression(const std::vector<Real>& mesh_time, const RegressionDataElliptic& regressionData, const UInt& nnodes_, const UInt& spline_degree):MixedFERegressionBase<RegressionDataElliptic>(mesh_time, regressionData, nnodes_, spline_degree){};
	template<UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
	void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
	{
	if(mydim!=2 || ndim !=2){

		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = this->regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& b =  this->regressionData_.getBeta();

	    MixedFERegressionBase<RegressionDataElliptic>::preapply<ORDER,mydim,ndim, IntegratorSpace, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>(c*mass+stiff[K]+dot(b,grad), ForcingTerm(std::vector<Real>(1)), mesh);
	}
	}
};

template<>
class MixedFERegression<RegressionDataEllipticSpaceVarying> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying>
{
public:
	MixedFERegression(const RegressionDataEllipticSpaceVarying& regressionData, const UInt& nnodes_):MixedFERegressionBase<RegressionDataEllipticSpaceVarying>(regressionData, nnodes_){};
	MixedFERegression(const std::vector<Real>& mesh_time, const RegressionDataEllipticSpaceVarying& regressionData, const UInt& nnodes_, const UInt& spline_degree):MixedFERegressionBase<RegressionDataEllipticSpaceVarying>(mesh_time, regressionData, nnodes_, spline_degree){};


	template< UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
	void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
	{
	if(mydim!=2 || ndim !=2){

		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = this->regressionData_.getC();
		const Diffusivity& K = this->regressionData_.getK();
		const Advection& b = this->regressionData_.getBeta();
		const ForcingTerm& u= this->regressionData_.getU();

		this->isSpaceVarying=TRUE;

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying>::preapply<ORDER,mydim,ndim, IntegratorSpace, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>(c*mass+stiff[K]+dot(b,grad), u, mesh);
	}
	}
};


// Temporal part
template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::setPhi(){

		Spline<Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
		UInt M = spline.num_knots()-SPLINE_DEGREE-1;
		UInt m = regressionData_.getNumberofTimeObservations();

		phi_.resize(m, M);
		Real value;

    for (UInt i = 0; i < m; ++i)
        for (UInt j = 0; j < M; ++j)
        {
					value = spline.BasisFunction(SPLINE_DEGREE, j, this->regressionData_.getTimeLocations()[i]);
					if (value!=0)
					{
						phi_.coeffRef(i,j) = value;
					}
				}
    phi_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::setTimeMass(){

    using ETTimeMass = EOExpr<TimeMass>;

    Spline<Integrator, SPLINE_DEGREE, 0> spline(mesh_time_);

    TimeMass ETimeMass;
    ETTimeMass timeMass(ETimeMass);

    Assembler::operKernel(timeMass, spline, timeMass_);
}

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::smoothSecondDerivative(){

    using ETTimeMass = EOExpr<TimeMass>;

    Spline<Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);

    TimeMass ETimeMass;
    ETTimeMass timeMass(ETimeMass);

    Assembler::operKernel(timeMass, spline, Pt_);
}

// Parabolic
template<typename InputHandler>
void MixedFDRegression<InputHandler>::setDerOperator(){
	//Build the matrix of finite differences [1/dt1 0 .......... 0]
	//																			 [-1/dt1 1/dt1 0 ....0]
	//																										...
	//																			 [0 ....0 -1/dtM 1/dtM]
	UInt M = mesh_time_.size()-1;
	derOpL_.resize(M, M);

	// set the first and the last rows
	Real delta = mesh_time_[1] - mesh_time_[0];
	derOpL_.coeffRef(0,0) = 1/delta;

	delta = mesh_time_[M-1] - mesh_time_[M-2];
	derOpL_.coeffRef(M-1,M-1) = 1/delta;
	derOpL_.coeffRef(M-1,M-2) = -1/delta;

	for (UInt i = 1; i < M-1; ++i)
	{
		delta = mesh_time_[i] - mesh_time_[i-1];
		derOpL_.coeffRef(i,i-1) = -1/delta;
		derOpL_.coeffRef(i,i) 	= 1/delta;
	}

	derOpL_.makeCompressed();

}


// template specification for GAMData
template<>
class MixedFERegression<GAMDataLaplace> : public MixedFERegressionBase<RegressionData>
{
public:
	MixedFERegression(const RegressionData& regressionData, const UInt& nnodes_):MixedFERegressionBase<RegressionData> (regressionData, nnodes_){};

	template<UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
	void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData>::preapply<ORDER,mydim,ndim, IntegratorSpace, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>(stiff, ForcingTerm(std::vector<Real>(1)), mesh);
	}



};

template<>
class MixedFERegression<GAMDataElliptic> : public MixedFERegressionBase<RegressionDataElliptic>
{
public:
	MixedFERegression( const RegressionDataElliptic& regressionData, const UInt& nnodes_):MixedFERegressionBase<RegressionDataElliptic> ( regressionData, nnodes_){};

	 template< UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
	void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
	{
		if(mydim!=2 || ndim !=2){

			Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");

		}else{
			typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		    const Real& c = this->regressionData_.getC();
		    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
		    const Eigen::Matrix<Real,2,1>& b = this->regressionData_.getBeta();


		    MixedFERegressionBase<RegressionDataElliptic>::preapply<ORDER,mydim,ndim, IntegratorSpace, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>(c*mass+stiff[K]+dot(b,grad), ForcingTerm(std::vector<Real>(1)), mesh);
		}
	}


};

template<>
class MixedFERegression<GAMDataEllipticSpaceVarying> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying>
{
public:
	MixedFERegression(const RegressionDataEllipticSpaceVarying& regressionData, const UInt& nnodes_):MixedFERegressionBase<RegressionDataEllipticSpaceVarying> ( regressionData, nnodes_){};

	template<UInt ORDER, UInt mydim, UInt ndim, typename IntegratorSpace, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
	void preapply(const MeshHandler<ORDER,mydim,ndim> & mesh)
	{
		if(mydim!=2 || ndim !=2){

			Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");

		}else{
			typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

			const Reaction& c = this->regressionData_.getC();
			const Diffusivity& K = this->regressionData_.getK();
			const Advection& b = this->regressionData_.getBeta();
			const ForcingTerm& u= this->regressionData_.getU();

			this->isSpaceVarying=TRUE;

			MixedFERegressionBase<RegressionDataEllipticSpaceVarying>::preapply<ORDER,mydim,ndim, IntegratorSpace, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE>
			(c*mass+stiff[K]+dot(b,grad), u, mesh);
		}
	}



};



#endif
