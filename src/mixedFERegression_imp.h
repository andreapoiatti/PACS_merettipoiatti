#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

//#include <libseq/mpi.h>
#include "../inst/include/dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::addDirichletBC()
{
	UInt id1,id3;

	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bc_indices[i];
			id3=id1+nnodes;

			A_.coeffRef(id1,id1)=pen;
			A_.coeffRef(id3,id3)=pen;


			_b(id1)+=bc_values[i]*pen;
			_b(id3)=0;
	 }

	A_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setPsi(){

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	psi_.resize(nlocations, nnodes);
	if (regressionData_.isLocationsByNodes()){

		std::vector<coeff> tripletAll;
		auto k = regressionData_.getObservationsIndices();
		tripletAll.reserve(k.size());
		for (int i = 0; i< k.size(); ++i){
			tripletAll.push_back(coeff(i,k[i],1.0));
		}
		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		psi_.makeCompressed();
	}
	else {
		constexpr UInt Nodes = mydim==2? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,Nodes,1> coefficients;

		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else
			{
				for(UInt node = 0; node < Nodes ; ++node)
				{
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<Nodes,mydim,ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}

		psi_.makeCompressed();
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::LeftMultiplybyQ(const MatrixXr& u)
{
	if (regressionData_.getCovariates().rows() == 0){
		return u;
	}
	else{
		MatrixXr W(this->regressionData_.getCovariates());
		if (isWTWfactorized_ == false ){
			WTWinv_.compute(W.transpose()*W);
			isWTWfactorized_=true;
		}
		MatrixXr Pu= W*WTWinv_.solve(W.transpose()*u);
		return u-Pu; //_più efficiente che non usare la Q direttamente, la vedo come ortogonale alla proiezione con P
	}

}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildA(const SpMat& Psi,  const SpMat& R1,  const SpMat& R0) {

	UInt nnodes = mesh_.num_nodes();

	SpMat DMat = Psi.transpose()*Psi;

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
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

	A_.setZero();
	A_.resize(2*nnodes,2*nnodes);
	A_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	A_.makeCompressed();
	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::system_factorize() {

	UInt nnodes = mesh_.num_nodes();

	// First phase: Factorization of matrix A
	Adec_.compute(A_);

	if (regressionData_.getCovariates().rows() != 0) {
		// Second phase: factorization of matrix  G =  C + [V * A^-1 * U]

		// Definition of matrix U = [ psi * W | 0 ]^T
		MatrixXr W(this->regressionData_.getCovariates());
		U_ = MatrixXr::Zero(2*nnodes, W.cols());
		U_.topRows(nnodes) = psi_.transpose()*W;

		// D = U^T * A^-1 * U
		MatrixXr D = U_.transpose()*Adec_.solve(U_);
		// G = C + D
		MatrixXr G = -W.transpose()*W + D;
		Gdec_.compute(G);
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> &b) {

	// Resolution of the system A * x1 = b
	MatrixXr x1 = Adec_.solve(b);

	if (regressionData_.getCovariates().rows() != 0) {
		// Resolution of G * x2 = U^T * x1
		MatrixXr x2 = Gdec_.solve(U_.transpose()*x1);
		// Resolution of the system A * x3 = U * x2
		x1 -= Adec_.solve(U_*x2);
	}

	return x1;
}


 template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setQ()
 {
 	//std::cout<<"Computing Orthogonal Space Projection Matrix"<<std::endl;
 	Q_.resize(H_.rows(),H_.cols());
 	Q_ = -H_;
 	for (int i=0; i<H_.rows();++i)
 	{
 		Q_(i,i) += 1;
 	}
 }

 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim, ndim>::setH()
 {
 	//std::cout<<"Computing Projection Matrix"<<std::endl;
 	//UInt nnodes = mesh_.num_nodes();
 	UInt nlocations = regressionData_.getNumberofObservations();

 	//regressionData_.printCovariates(std::cout);
 	MatrixXr W(this->regressionData_.getCovariates());
 	//std::cout<<"W "<< W <<std::endl;
 	//total number of mesh nodes
 	//UInt nnodes = mesh_.num_nodes();
 	if(regressionData_.isLocationsByNodes())
 	{
 		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
 		for (auto i=0; i<nlocations;++i)
 		{
 			auto index_i = regressionData_.getObservationsIndices()[i];
 			for (auto j=0; j<W.cols();++j)
 			{
 				W_reduced(i,j) = W(index_i,j);
 			}
 		}
 		W = W_reduced;
 	}


 	MatrixXr WTW(W.transpose()*W);

 	H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix
 }


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
 {
 	//I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
 	//_coeffmatrix.setFromTriplets(tripletA.begin(),tripletA.end());

 	UInt nnodes = mesh_.num_nodes();

 	std::vector<coeff> tripletAll;
 	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

 	for (int k=0; k<DMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(DMat,k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
 	  }
 	for (int k=0; k<MMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(MMat,k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
 	  }
 	for (int k=0; k<AMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
 	  }
 	for (int k=0; k<AMat.outerSize(); ++k)
 	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
 	  {
 		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
 	  }

 	_coeffmatrix.setZero();
 	_coeffmatrix.resize(2*nnodes,2*nnodes);
 	_coeffmatrix.setFromTriplets(tripletAll.begin(),tripletAll.end());
 	_coeffmatrix.makeCompressed();
 	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
 }

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::getDataMatrix(SpMat& DMat)
 {
 		UInt nnodes = mesh_.num_nodes();
 		//UInt nlocations = regressionData_.getNumberofObservations();

 		DMat.resize(nnodes,nnodes);

 		if (regressionData_.getCovariates().rows() == 0)
 			DMat = psi_.transpose()*psi_;
 		else
 		{
 			DMat = (SpMat(psi_.transpose())*Q_*psi_).sparseView();
 		}
 }

 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::getDataMatrixByIndices(SpMat& DMat)
 {
 		UInt nnodes = mesh_.num_nodes();
 		UInt nlocations = regressionData_.getNumberofObservations();

 		DMat.resize(nnodes,nnodes);

 		if (regressionData_.getCovariates().rows() == 0)
 		{
 			DMat.reserve(1);
 			for (auto i = 0; i<nlocations; ++i)
 			{
 				auto index = regressionData_.getObservationsIndices()[i];
 				DMat.insert(index,index) = 1;
 			}
 		}
 		else
 		{
 			//May be inefficient
 			for (auto i = 0; i<nlocations; ++i)
 			{
 				auto index_i = regressionData_.getObservationsIndices()[i];
 				for (auto j = 0; j<nlocations; ++j)
 				{
 					auto index_j = regressionData_.getObservationsIndices()[j];
 					DMat.insert(index_i,index_j) = Q_(i,j);
 				}
 			}
 		}
 }

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData=psi_.transpose()*regressionData_.getObservations();
		}
	}
	else
	{
		rightHandData=psi_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedom(UInt output_index, Real lambda)
{
	int GCVmethod = regressionData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_index, lambda);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_index, lambda);
			break;
	}
}
//_NB degrees of freedom sono tr(S)
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedomExact(UInt output_index, Real lambda)
{

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees=0; //_è un numero reale, il valore della traccia di S

	// Case 1: MUMPS
	if (regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() == 0 )
	{
		auto k = regressionData_.getObservationsIndices();
		DMUMPS_STRUC_C id;
		int myid, ierr;
        int argc=0;
        char ** argv= NULL;
        //MPI_Init(&argc,&argv);
		//ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

		id.sym=0;
		id.par=1;
		id.job=JOB_INIT;
		id.comm_fortran=USE_COMM_WORLD;
		dmumps_c(&id);

		std::vector<int> irn;
		std::vector<int> jcn;
		std::vector<double> a;
		std::vector<int> irhs_ptr;
		std::vector<int> irhs_sparse;
		double* rhs_sparse= (double*)malloc(nlocations*sizeof(double));

		//if( myid==0){
			id.n=2*nnodes;
			for (int j=0; j<A_.outerSize(); ++j){
				for (SpMat::InnerIterator it(A_,j); it; ++it){
					irn.push_back(it.row()+1);
					jcn.push_back(it.col()+1);
					a.push_back(it.value());
				}
			}
		//}
		id.nz=irn.size();
		id.irn=irn.data();
		id.jcn=jcn.data();
		id.a=a.data();
		id.nz_rhs=nlocations;
		id.nrhs=2*nnodes;
		int j = 1;
		irhs_ptr.push_back(j);
		for (int l=0; l<k[0]-1; ++l) {
			irhs_ptr.push_back(j);
		}
		for (int i=0; i<k.size()-1; ++i) {
			++j;
			for (int l=0; l<k[i+1]-k[i]; ++l) {
				irhs_ptr.push_back(j);
			}

		}
		++j;
		for (int i=k[k.size()-1]; i < id.nrhs; ++i) {
			irhs_ptr.push_back(j);
		}
		for (int i=0; i<nlocations; ++i){
			irhs_sparse.push_back(k[i]+1);
		}
		id.irhs_sparse=irhs_sparse.data();
		id.irhs_ptr=irhs_ptr.data();
		id.rhs_sparse=rhs_sparse;

		#define ICNTL(I) icntl[(I)-1]
		//Output messages suppressed
		id.ICNTL(1)=-1;
		id.ICNTL(2)=-1;
		id.ICNTL(3)=-1;
		id.ICNTL(4)=0;
		id.ICNTL(20)=1;
		id.ICNTL(30)=1;
		id.ICNTL(14)=200;

		id.job=6;
		dmumps_c(&id);
		id.job=JOB_END;
		dmumps_c(&id);

		//if (myid==0){
			for (int i=0; i< nlocations; ++i){
				//std::cout << "rhs_sparse" << rhs_sparse[i] << std::endl;
				degrees+=rhs_sparse[i];
			}
		//}
		free(rhs_sparse);

		//MPI_Finalize();
	}
	// Case 2: Eigen
	else{
		MatrixXr X1 = psi_.transpose() * LeftMultiplybyQ(psi_);
                MatrixXr Aux=MatrixXr::Identity(nlocations,nlocations)-Q_; //_used for z_hat
		if (isRcomputed_ == false ){
			isRcomputed_ = true;
			Eigen::SparseLU<SpMat> solver;
			solver.compute(R0_);
			auto X2 = solver.solve(R1_); //_R0^-1*R1
			R_ = R1_.transpose() * X2; //_R1^T*R0^-1*R1
		}

		MatrixXr X3 = X1 + lambda * R_; //_psi^T*Q*Psi+lambda*R1^T*R0\R1
                SS_=X3; //to store for derivative of GCV
		Eigen::LDLT<MatrixXr> Dsolver(X3);

		auto k = regressionData_.getObservationsIndices(); //_cerco gl indici delle locations: pur coincidendo con le locations, i nodi possono avere indici diversi)
//NB migliorabile: nel caso isLocationsByNodes()==true, non ha senso aver usato psi, è l'identità, si possono risparmiare passaggi-> anche se le locations possono essere meno dei nodi!!!
		if(regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() != 0)
		{
		 //_number of rows of covariates, serve solo per controllo, verifico che ci siano covariate
			degrees += regressionData_.getCovariates().cols();  //_calcola q+tr(S), questo è q, numero di covariate, numero colonne della matrice delle covariate

			// Setup rhs B
			MatrixXr B; //_NB: se i nodi coincidono con i pi, la matrice Psi è l'identità, perchè diventa psi_i(p_j)=delta_ij
			B = MatrixXr::Zero(nnodes,nlocations);
			// B = I(:,k) * Q
			for (auto i=0; i<nlocations;++i) {
				VectorXr ei = VectorXr::Zero(nlocations);
				ei(i) = 1;
				VectorXr Qi = LeftMultiplybyQ(ei);
				for (int j=0; j<nlocations; ++j) {
					B(k[i], j) = Qi(j); //_riempio i nodi, dove ho anche le locations, che possono avere indici diversi rispetto ai nodi
				}
			}
			// Solve the system TX = B
			//_B=Q (PSI è l'identità)
			MatrixXr X;
			X = Dsolver.solve(B); //_X3^-1*B (psi identità)
			V_=X; //è già quella che cerco!
			// Compute trace(X)->ovviamente uso k,i come indici perchè gugarda i nodi con la location nell'ordine dei location
			for (int i = 0; i < k.size(); ++i) {
				degrees += X(k[i], i);
			}
		}
	             VectorXr z;
    	   	if(regressionData_.isLocationsByNodes())
    	   	{
    		   z=VectorXr::Zero(nlocations);
    		   for(auto i=0;i<regressionData_.getObservationsIndices().size();i++)
    			   z(regressionData_.getObservationsIndices()[i])=regressionData_.getObservationData()[i]; //_mette i valori nei nodi in cui si pongono le locations (nodi coincidono con le location, ma le locations possono avere ordini di numerazione diversi!!
    	   	}
    	   	else {
    		   	z=regressionData_.getObservationData();
    		    	}

		if (!regressionData_.isLocationsByNodes()){ //_ora psi non è l'identità!! (nb stu-hunter SAngalli pdf da pag10 in poi)
		//_usa la proprietà della traccia: tr(psi*(X3^-1)*psi^T*Q)=tr((X3^-1)*psi^T*Q*psi), quindi
		//_calcola (X3^-1)*psi^T*Q*psi=(X3^-1)*X1 così ha già calcolato il pezzo psi^T*Q*psi in calcolo precedente

			MatrixXr X, X4; //_needed the correct formulation of S, not the computation for the trace
			//X = Dsolver.solve(MatrixXr(X1)); //_X3^-1*X1
                        X4=LeftMultiplybyQ(psi_.transpose());
			X = Dsolver.solve(MatrixXr(X4)); //_X3^-1*X4
                        V_=X;
			X=psi_*X; //_ora X è la S!
			z_hat_=(Aux+LeftMultiplybyQ(X))*z;
			if (regressionData_.getCovariates().rows() != 0) {
				degrees += regressionData_.getCovariates().cols();  //_calcola q+tr(S), questo è q, numero di covariate, numero colonne della matrice delle covariate
			}
			for (int i = 0; i<nnodes; ++i) { //_???perchè nnodes e non nlocations? S è nxn!! (forse n<nnodes e quindi è lo stesso, aggiungo zeri)
				degrees += X(i,i); //_computes the trace of the matrix S (il +q è già stato calcolato)
			}
		}
	}
	_dof[output_index] = degrees;
}

//_Verificare a cosa seve la stochastic???

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedomStochastic(UInt output_index, Real lambda)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	std::default_random_engine generator;
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
//_IN Negri la Q è l'identità
	// Define the first right hand side : | I  0 |^T * psi^T * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	b.topRows(nnodes) = psi_.transpose()* LeftMultiplybyQ(u);

	// Resolution of the system
	//MatrixXr x = system_solve(b);
	Eigen::SparseLU<SpMat> solver;
	solver.compute(_coeffmatrix);
	auto x = solver.solve(b);

	MatrixXr uTpsi = u.transpose()*psi_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (regressionData_.getCovariates().rows() != 0) {
		q = regressionData_.getCovariates().cols();
	}
	// For any realization we calculate the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {
		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof[output_index] = mean;
}




template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper)
{
	UInt nnodes=mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	setPsi();

	if(!regressionData_.getCovariates().rows() == 0)
 	{
 		setH();
 		setQ();
 	}

	if(!regressionData_.isLocationsByNodes())
	{
		getDataMatrix(DMat_);
	}
	else
	{
		getDataMatrixByIndices(DMat_);
	}

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);
	VectorXr rightHandData;
	getRightHandData(rightHandData);
	this->_b = VectorXr::Zero(2*nnodes);
	this->_b.topRows(nnodes)=rightHandData;
	this->_solution.resize(regressionData_.getLambda().size());
	this->_dof.resize(regressionData_.getLambda().size());

	for(UInt i = 0; i<regressionData_.getLambda().size(); ++i) //_cacloal già per tutte le lambda e passa all'esterno, poi via R si usa la getGCV e si ottiene il calcolo
	{
		Real lambda = regressionData_.getLambda()[i];
		SpMat R1_lambda = (-lambda)*R1_;
		SpMat R0_lambda = (-lambda)*R0_;
		this->buildA(psi_, R1_lambda, R0_lambda);
		this->buildCoeffMatrix(DMat_, R1_lambda, R0_lambda);
		//std::cout << coeffmatrix_ << std::endl;
		//Applying boundary conditions if necessary
		if(regressionData_.getDirichletIndices().size() != 0)
			addDirichletBC();

		system_factorize();
		_solution[i] = this->template system_solve(this->_b);
		if(regressionData_.computeDOF())  //_NB compute dof semplicemente dice se dovrà calcolare i dof o meno, se è true, ci mette il calcolo, che è tr(S)+q
			{
				computeDegreesOfFreedom(i,lambda);
				computeGCV(i); //_to check if it's correct, stampa il valore della GCV
                                computeGCV_derivative(i); //_to check if it's correct, stampa il valore della derivata dGCV/dlambda
                                //NB ovviamente non si userà più il vettore dei dofs.
				//NB fare confronto di questi valori della GCV con quelli calcolati tramite R
			 }


		else
			_dof[i] = -1;

	}
} //_calcola già il vettore completo dei dofs( uno per ogni lambda), che sono i valori della tr(S)+q, che poi passa esternamente a R che le usa nella getGCV

//_computation of GCV (possibile implementazione)
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
Real MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeGCV(UInt output_index)
{ //_da modificare il fatto che _dof sia un vettore, serve un solo valore, non serve output index!!
  //_trovare modo efficiente di calcolare le zhat, come (I-Q+QS)*z oppure posso usare la Hat matrix...cercare metdo efficiente
	UInt s;
	//_UInt q=regressionData_.getCovariates().cols(); //_serve se si vuole usare l'articolo di stuHuntersangalli, è già implicito nel degrees of freedom
	VectorXr z;
	s= regressionData_.getNumberofObservations(); //_così ho anche il caso in uci ho meno locations dei nodi (pur coincidenti)
	if(regressionData_.isLocationsByNodes())
	{
		//s= this->mesh_.num_nodes(); //_vuol dire che il numero di locations (e quindi il numero di osservazioni, coincide col numero di nodi e le posizioni sono esattamente quelle dei nodi)
    //_s è la n dell'articolo stuHuntersangalli pdf pag.12, numero locations, dove ho le osservazioni
		z=VectorXr::Zero(s);
		for(auto i=0;i<regressionData_.getObservationsIndices().size();i++)
			z(regressionData_.getObservationsIndices()[i])=regressionData_.getObservationData()[i]; //_mette i valori nei nodi in cui si pongono le locations (nodi coincidono con le location, ma le locations possono avere ordini di numerazione diversi!!
	} else {

		z=regressionData_.getObservationData();
	}
	//_calcolo z_hat, suppongo di avere i valori come per z
	Real norm_squared=(z-z_hat_).transpose()*(z-z_hat_);
	if (s-dof_[output_index]<0) { //_dof_ non servirà, sarà un valore unico!
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->regressionData_.getLambda()[output_index]);
			#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->regressionData_.getLambda()[output_index] <<"\n";
			#endif
			}
	Real stderror=norm_squared/(s-dof_[output_index]); //così è ancora fatta sul vettore
        Real GCV_val=(s/(s-dof_[output_index]))*stderror;
	#ifdef R_VERSION_
		Rprintf("GCV=%f\n",GCV_val);
	#else
		std::cout << "GCV value="<<GCV_val<<std::endl;
	#endif

	return GCV_val; //_Calcolo della GCV come s*(z-zhat)^T*(z-zhat)/(s-(q+trS))^2

}


//_computation of GCV_derivative (possibile implementazione)
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
Real MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeGCV_derivative(UInt output_index)
{ //_da modificare il fatto che _dof sia un vettore, serve un solo valore, non serve output index!!
  //_trovare modo efficiente di calcolare le zhat, come (I-Q+QS)*z oppure posso usare la Hat matrix...cercare metdo efficiente
	UInt s;
	//_UInt q=regressionData_.getCovariates().cols(); //_serve se si vuole usare l'articolo di stuHuntersangalli, è già implicito nel degrees of freedom
	VectorXr z;
	s= regressionData_.getNumberofObservations(); //_così ho anche il caso in uci ho meno locations dei nodi (pur coincidenti)
	if(regressionData_.isLocationsByNodes())
	{
		//s= this->mesh_.num_nodes(); //_vuol dire che il numero di locations (e quindi il numero di osservazioni, coincide col numero di nodi e le posizioni sono esattamente quelle dei nodi)
    //_s è la n dell'articolo stuHuntersangalli pdf pag.12, numero locations, dove ho le osservazioni
		z=VectorXr::Zero(s);
		for(auto i=0;i<regressionData_.getObservationsIndices().size();i++)
			z(regressionData_.getObservationsIndices()[i])=regressionData_.getObservationData()[i]; //_mette i valori nei nodi in cui si pongono le locations (nodi coincidono con le location, ma le locations possono avere ordini di numerazione diversi!!
	}
	 else
	 {

		z=regressionData_.getObservationData();
	  }

	//NB _questo caso è ok se i nodi non coincidono con le location, è ridondante (migliorabile!!) se coindicono, perchè psi è l'identitò (si può fare come nel calcolo dei deg of freedom per essere più efficiente)
	Real norm_squared=(z-z_hat_).transpose()*(z-z_hat_);
        Real trace_=0.0;
if(regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() != 0) //da controllare per questione matrie identità psi che usa o meno

{       auto k = regressionData_.getObservationsIndices();
	Eigen::LDLT<MatrixXr> Dsolver( SS_ );
	MatrixXr d_S=-Dsolver.solve( R_*V_ ); //_psi è identità //_se dà errore, provare MatrixXr(R_*V_), per ricreare al più la matrice
	//_d_S=-psi*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*psi^T*Q

	for (int i = 0; i < k.size(); ++i)
	{
		trace_ += dS_(k[i], i);
 	}
}
else
       {
	Eigen::LDLT<MatrixXr> Dsolver( SS_ );
	MatrixXr d_S=-psi_*Dsolver.solve( R_*V_ ); //_se dà errore, provare MatrixXr(R_*V_), per ricreare al più la matrice
	//_d_S=-psi*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*psi^T*Q


        for (Uint i=0; i<mesh_.num_nodes(); i++) //_anche se sarebbe più corretto il numero di osservazioni, è nxn
		trace_+=d_S(i,i); //_tr(dS/dlambda)=d(tr(S))/dlambda
	}

	if(s-dof_[output_index]<0){ //_dof_ non servirà, sarà un valore unico!
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->regressionData_.getLambda()[output_index]);
			#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->regressionData_.getLambda()[output_index] <<"\n";
			#endif
			}
	Real stderror=norm_squared/(s-dof_[output_index]); //così è ancora fatta sul vettore
	Real GCV_der_val=2*(s/((s-dof_[output_index]) * (s-dof_[output_index])))*stderror*trace_;
	#ifdef R_VERSION_
		Rprintf("GCV=%f\n",GCV_der_val);
	#else
		std::cout << "GCV value="<<GCV_der_val<<std::endl;
	#endif

	return GCV_der_val; //_Calcolo della derivata della GCV

}







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
