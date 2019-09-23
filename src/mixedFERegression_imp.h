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
		Eigen::Matrix<Real ,Nodes, 1> 	coefficients;
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

 	if (regressionData_.isLocationsByNodes())
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
	UInt nnodes = mesh_.num_nodes();
	DMat.resize(nnodes, nnodes);

	// Check if Q == Identity
	if (regressionData_.getCovariates().rows() == 0)
		DMat = psi_.transpose()*psi_;
	else
	{
		DMat = (SpMat(psi_.transpose())*Q_*psi_).sparseView();
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::getDataMatrixByIndices(SpMat & DMat)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	DMat.resize(nnodes, nnodes);

	if (regressionData_.getCovariates().rows() == 0)
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
				DMat.insert(index_i, index_j) = Q_(i,j);
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

	// REMARK: the minus in top left [???]

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

	if (regressionData_.getCovariates().rows() == 0)
	{
		if (regressionData_.isLocationsByNodes())
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
	if (regressionData_.getCovariates().rows() == 0)
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
// DOF

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator, ORDER, mydim,ndim>::computeDegreesOfFreedom(UInt output_index, Real lambda)
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
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::computeDegreesOfFreedomExact(UInt output_index, Real lambda)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees = 0; 	//_è un numero reale, il valore della traccia di S
				// Noi ne dovremo restiuire uno solo (ora finisce ancora nel vettore dei dof)

	// Case 1: MUMPS [location_by_nodes && no covariates Q = Identity]
	if (regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() == 0)
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
	else
	{
		MatrixXr X1 = psi_.transpose() * LeftMultiplybyQ(psi_); // Top-left factor of A [[POSSIBILE MIGLIORAMENTO]]
									// [[ Siccome la condizioe prima dell'else è un && Q potrebbe ancora esere Identità,
									// Non sarebbe opportuno fare un check per evitare un calcolo inutile? if regressionData_.getCovariates().rows() != 0 ... else ... ]]

		// Check if R_ has already been computed
		if (isRcomputed_ == false)
		{
			isRcomputed_ = true;		// from now on it won't be computed again
			Eigen::SparseLU<SpMat> solver;	// define a factorized empty sparse LU solver
			solver.compute(R0_);		// apply it to R0 to simplify the inverse
			auto X2 = solver.solve(R1_); 	// compute _R0^{-1}*R1 with the sparse LU factorization
			R_ = R1_.transpose() * X2; 	// R == _R1^t*R0^{-1}*R1
		}

		MatrixXr X3 = X1 + lambda*R_; 		// X3 is SS_ in our code: Psi^t*Q*Psi + lambda*R1^t*R0^{-1}*R1
		SS_ = X3;				// Storing the temporary in SS_ [[ cant' we aviod the teporary X3, ci serve per riciclarlo tra diversi lambda ???]]
		Eigen::LDLT<MatrixXr> Dsolver(X3);	// Factorize also X3 [[ We might store this instead of SS_ ???]]

		// [[ TOP: Cercare un metodo iterativo pe iroslevere  il fascio di matrici (A+lambda*B)^-1]]

		auto k = regressionData_.getObservationsIndices(); //_cerco gli indici delle locations: pur coincidendo con le locations, i nodi possono avere indici diversi)
									// [[ Non ne sono molto covinto, sia per il dicscorso che Psi!=Identità molto spesso, sia perchè la condizione
									// del blocco else di cui siamo interi è un && e non un ||]]

		//NB migliorabile: nel caso isLocationsByNodes()==true, non ha senso aver usato psi, è l'identità, si possono risparmiare passaggi [[vedi nota superiore]]
		VectorXr z;
		//MatrixXr I_=MatrixXr::Identity(nlocations,nlocations);
		//if (regressionData_.isLocationsByNodes())
		//{
			//z = VectorXr::Zero(nlocations);
			//for (auto i = 0; i < regressionData_.getObservationsIndices().size(); i++)
			     // z(regressionData_.getObservationsIndices()[i]) = regressionData_.getObservationData()[i]; //_mette i valori nei nodi in cui si pongono le locations (nodi coincidono con le location, ma le locations possono avere ordini di numerazione diversi!!
		//}
		//else
		//{
		z = regressionData_.getObservations();
		//_number of rows of covariates, serve solo per controllo, verifico che ci siano covariate [[ in teoria è ovvio dall'else, controllo superfluo
		// se la prima condizione è soddisfatta la seconda lo deve essere per non violare il blocco di elese in cui ci troviamo, al limite si commenta di fianco
		// ma si evita ogni volta un controllo]]
		if (regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() != 0)
		{
			degrees += regressionData_.getCovariates().cols();  //_ alla fine voglio q+tr(S); questo è solo q, numero di covariate (numero colonne della matrice delle covariate)

			// Setup rhs B
			MatrixXr B; //_NB: se i nodi coincidono con TUTTI i p_i NELLO STESSO ORDINE, la matrice Psi è l'identità, perchè diventa psi_i(p_j)=delta_ij
			B = MatrixXr::Zero(nnodes, nlocations); // B = Psi^t*Q, Q != Identity in the block of the "if".
			// B = I(:,k) * Q


			for (auto i = 0; i < nlocations; ++i)
			{
				VectorXr ei = VectorXr::Zero(nlocations);
				ei(i) = 1;
				VectorXr Qi = LeftMultiplybyQ(ei);
				for (int j = 0; j < nlocations; ++j)
				{
					B(k[i], j) = Qi(j);  //_La B si calcola inserendo le righe i-esime di Q in posizione riga k[i] (nodo k[i]) di B, che è una matrice nnodes x nlocations
			        }
		        }
			// Solve the system TX = B
			MatrixXr X;
			S_=MatrixXr(nlocations,nlocations);

			X = Dsolver.solve(B); //_X3^-1*B
			V_ = X;
			// [[Ancora perchè creare un temporaneo se non serve a niente, solo se dovremo riciclare il vecchio valore???]]
			//_creo il temporaneo perchè non so se V sarà ancora quella giusta o dà problemi perchè è quella del ciclo precedente...per evitare problemi

			// Compute trace(X)->ovviamente uso k,i come indici perchè gugarda i nodi con la location nell'ordine dei location
			// [[ Più precisamene S = Psi*V e qui chiamiamo V con il nome X; k è dim(getObservationsIndices)==nlocations, infatti
			// Psi è nloc x nnodes e V è nnod x nloc, alla fine S è nloc x nloc e devo prendere la traccia quindi sommo
			// sum_{i=1}^{n_locations} S(i,i), ove S(i,i) = sum{j=1}^{nnodes} Psi(i,j)V(j,i), tuttavia,  per il fatto che
			// isLocationsByNodes ==  true come specificato in setPsi(), Psi(i,j) = delta_{i,k[i]} perciò S(i,i) == 1*V(k[i],i) == X(k[i],i)]]
			// Calcolo l'azione della pre-moltiplicazione per Psi: S = Psi * V (in pratica estraggo solo le nlocations righe della V (cioè della X) corrispondenti agli indici presenti in k e di queste calcolo la traccia->e.g.posizione 0 (riga 0 della matrice psi*V) corrisponde a riga k[0] in V, perciò il valore della traccia sarà V[k[0],0]), ecc...
			for (int i = 0; i < k.size(); ++i) {
				degrees += X(k[i], i);
			}
			for (int i = 0; i < k.size(); ++i)
			    for (int j = 0; j < k.size(); ++j)
			 {
				S_(i,j) = X(k[i], j);
			}

		z_hat_ 	= (H_+LeftMultiplybyQ(S_))*z;
		//z_hat_ 	= (I_-LeftMultiplybyQ(I_))*z;


		}


   		//}
		// [[???]]
		// [[Non mi è chiara la logica dietro questo passaggio: in particolare non capisco perchè facciamo quello che c'è nell'if
		// z è il vettore delle locations che sono un SOTTOINSIEME dei nodi (perchè alcuni osservaioni possono essere lette da R)
		// come NA proprio per segnalare che il nodo corispondente NON è una location, quindi n_locations è in principio minore dei nodi
		// Se noi facciamo come è critto qui staremmo prendendo l'osservazione relativa alla location i-esima e mettendola dentro z non
		// in posizione reale, ma in posizione del suo corrisponente nodo. Primo, a cosa mi serve? Secondo, siccome il numero di nodi è superiore
		// di base al numero di locations simmone io faccio un resizing di z l numero delle locations alcuni nodi potrebbero sforare la dimensione
		// dichairata del vettore (maaagari eigen fa un resizing automatico del vettore in questo caso, ma non ne ho idea..), in ogni caso non vedo la
		// ragione dietro questa procedura la z finale appartiene a R^{n_locations} e non R^{n_nodes}...]]
		if (!regressionData_.isLocationsByNodes())
		{ //_ora psi non è l'identità[[o meglio di soli 1 e 0]]!! (nb stu-hunter SAngalli pdf da pag10 in poi)
		//_usa la proprietà della traccia: tr(psi*(X3^-1)*psi^T*Q)=tr((X3^-1)*psi^T*Q*psi), quindi
		//_calcola (X3^-1)*psi^T*Q*psi=(X3^-1)*X1 così ha già calcolato il pezzo psi^T*Q*psi in calcolo precedente

			MatrixXr X, X4; //_needed the correct formulation of S, not the computation for the trace
			//X = Dsolver.solve(MatrixXr(X1)); //_X3^-1*X1
			X4 	= LeftMultiplybyQ(psi_.transpose());
			X 	= Dsolver.solve(MatrixXr(X4)); //_X3^-1*X4
			V_	= X;
			S_ 	= psi_*X;
			z_hat_ 	= (H_+LeftMultiplybyQ(S_))*z;
			//z_hat_ 	= (I_-LeftMultiplybyQ(I_))*z;


			if (regressionData_.getCovariates().rows() != 0)
			{
				degrees += regressionData_.getCovariates().cols();  //_calcola q+tr(S), questo è q, numero di covariate, numero colonne della matrice delle covariate
			}
			for (int i = 0; i<nnodes; ++i) //_???perchè nnodes e non nlocations? S è nxn!! (forse n<nnodes e quindi è lo stesso, aggiungo zeri)
							//[[Concordo penso che sia un errore]] [[POSSIBILE ERRORE]]
			{
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

//----------------------------------------------------------------------------//
// Factorizer

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::system_factorize()
{
	UInt nnodes = mesh_.num_nodes();

	// First phase: Factorization of matrix A
	Adec_.compute(A_);

	// If Q!=Identity
	if (regressionData_.getCovariates().rows() != 0)
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

	if (regressionData_.getCovariates().rows() != 0)
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

	// Set matrix Psi
	setPsi();

	// If there are covariates in the model set H and Q
	if (!regressionData_.getCovariates().rows() == 0)
 	{
 		setH();
 		setQ();
 	}

	// Set the data matrix
	if (!regressionData_.isLocationsByNodes())
	{
		getDataMatrix(DMat_);
	}
	else
	{
		getDataMatrixByIndices(DMat_);
	}

	// Assemble matrices R0-R1
	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);

	// Define Right hand data
	VectorXr rightHandData;
	getRightHandData(rightHandData);

	// Define vector b
	this->_b = VectorXr::Zero(2*nnodes);
	this->_b.topRows(nnodes) = rightHandData;

	// Resize solution to be a vector [[!!! TODO !!!]]
	this->_solution.resize(regressionData_.getLambda().size());
	this->_dof.resize(regressionData_.getLambda().size());

	// Cycle to compute solutions [[!!! TODO !!!]]
	for (UInt i = 0; i < regressionData_.getLambda().size(); ++i)
	//_cacloal già per tutte le lambda e passa all'esterno, poi via R si usa la getGCV e si ottiene il calcolo
	{
		Real  lambda 	= regressionData_.getLambda()[i];
		SpMat R1_lambda = (-lambda) * R1_;
		SpMat R0_lambda = (-lambda) * R0_;

		// Build A and coeffmatrix
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
		_solution[i] = this->template system_solve(this->_b);

		// Compute dof is necessary
		if (regressionData_.computeDOF())  //_NB compute dof semplicemente dice se dovrà calcolare i dof o meno, se è true, ci mette il calcolo, che è tr(S)+q
		{
			computeDegreesOfFreedom(i,lambda);
			Real t_=computeGCV(i); //_to check if it's correct, stampa il valore della GCV
	                Real p_=computeGCV_derivative(i); //_to check if it's correct, stampa il valore della derivata dGCV/dlambda
	                //NB ovviamente non si userà più il vettore dei dofs.
			//NB fare confronto di questi valori della GCV con quelli calcolati tramite R
		 }
		else
			_dof[i] = -1;

	}
} //_calcola già il vettore completo dei dofs( uno per ogni lambda), che sono i valori della tr(S)+q, che poi passa esternamente a R che le usa nella getGCV

//----------------------------------------------------------------------------//
// GCV [[POiatti]]

// Inserire poia
//_computation of GCV (possibile implementazione)
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::computeGCV(UInt output_index)
{ //_da modificare il fatto che _dof sia un vettore, serve un solo valore, non serve output index!!
  //_trovare modo efficiente di calcolare le zhat, come (I-Q+QS)*z oppure posso usare la Hat matrix...cercare metdo efficiente
	UInt s;
	//_UInt q=regressionData_.getCovariates().cols(); //_serve se si vuole usare l'articolo di stuHuntersangalli, è già implicito nel degrees of freedom
	VectorXr z;
	s= regressionData_.getNumberofObservations(); //_così ho anche il caso in cui ho meno locations dei nodi (pur coincidenti)
		//s= this->mesh_.num_nodes(); //_vuol dire che il numero di locations (e quindi il numero di osservazioni, coincide col numero di nodi e le posizioni sono esattamente quelle dei nodi)
    //_s è la n dell'articolo stuHuntersangalli pdf pag.12, numero locations, dove ho le osservazioni
         //ho eliminato la parte di cui abbiamo discusso, usata da Negri
	z=regressionData_.getObservations();
	//_calcolo z_hat, suppongo di avere i valori come per z
	Real norm_squared=(z-z_hat_).transpose()*(z-z_hat_);
	//Real norm_squared=(z).transpose()*(z);
	if (s-_dof[output_index]<0) { //_dof_ non servirà, sarà un valore unico!
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->regressionData_.getLambda()[output_index]);
			#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->regressionData_.getLambda()[output_index] <<"\n";
			#endif
			}
	Real stderror=norm_squared/(s-_dof[output_index]); //così è ancora fatta sul vettore
        Real GCV_val=(s/(s-_dof[output_index]))*stderror;
	#ifdef R_VERSION_
		Rprintf("GCV=%f\n",GCV_val);
	#else
		std::cout << "GCV value="<<GCV_val<<std::endl;
	#endif

	return GCV_val; //_Calcolo della GCV come s*(z-zhat)^T*(z-zhat)/(s-(q+trS))^2

}


//_computation of GCV_derivative (possibile implementazione)
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::computeGCV_derivative(UInt output_index)
{ //_da modificare il fatto che _dof sia un vettore, serve un solo valore, non serve output index!!
  //_trovare modo efficiente di calcolare le zhat, come (I-Q+QS)*z oppure posso usare la Hat matrix...cercare metdo efficiente
	UInt s;
	//_UInt q=regressionData_.getCovariates().cols(); //_serve se si vuole usare l'articolo di stuHuntersangalli, è già implicito nel degrees of freedom
	VectorXr z;
	s= regressionData_.getNumberofObservations(); //_così ho anche il caso in uci ho meno locations dei nodi (pur coincidenti)
	//if(regressionData_.isLocationsByNodes())
	//{
		//s= this->mesh_.num_nodes(); //_vuol dire che il numero di locations (e quindi il numero di osservazioni, coincide col numero di nodi e le posizioni sono esattamente quelle dei nodi)
    //_s è la n dell'articolo stuHuntersangalli pdf pag.12, numero locations, dove ho le osservazioni
		//z=VectorXr::Zero(s);
		//for(auto i=0;i<regressionData_.getObservationsIndices().size();i++)
		//	z(regressionData_.getObservationsIndices()[i])=regressionData_.getObservationData()[i]; //_mette i valori nei nodi in cui si pongono le locations (nodi coincidono con le location, ma le locations possono avere ordini di numerazione diversi!!
	//}
	// else
	// {

	z=regressionData_.getObservations();
	  //}
        MatrixXr I=MatrixXr::Identity(s,s);
	//NB _questo caso è ok se i nodi non coincidono con le location, è ridondante (migliorabile!!) se coindicono, perchè psi è l'identitò (si può fare come nel calcolo dei deg of freedom per essere più efficiente)
	Real norm_squared=(z-z_hat_).transpose()*(z-z_hat_);
        Real trace_=0.0;
	MatrixXr dS_(s,s); //S derivative dlambda
       //al più si può migliorare evitando l'ultimo prodotto per psi e ragionando con la k
	Eigen::LDLT<MatrixXr> Dsolver( SS_ );
	dS_=-psi_*Dsolver.solve( MatrixXr(R_*V_) ); //_se dà errore, provare MatrixXr(R_*V_), per ricreare al più la matrice
	//_d_S=-psi*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*psi^T*Q


        //for (UInt i=0; i<mesh_.num_nodes(); i++) //_anche se sarebbe più corretto il numero di osservazioni, è nxn

	for (UInt i=0; i<s; i++)
		trace_+=dS_(i,i); //_tr(dS/dlambda)=d(tr(S))/dlambda

	if(s-_dof[output_index]<0){ //_dof_ non servirà, sarà un valore unico!
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->regressionData_.getLambda()[output_index]);
 			#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->regressionData_.getLambda()[output_index] <<"\n";
			#endif
			}
	Real stderror=norm_squared/(s-_dof[output_index]); //così è ancora fatta sul vettore
       //uso la proprietà di simmetri e idempotenza di Q, ho Q^T*Q=Q*Q=Q
	VectorXr second_=s/((s-_dof[output_index])*(s-_dof[output_index]))*(z.transpose()*(-dS_.transpose())*LeftMultiplybyQ(MatrixXr(I-S_))*z+z.transpose()*(I-S_.transpose())*LeftMultiplybyQ(MatrixXr(-dS_))*z); //prodotto per matrici può essere un vettore, non un Real->estraggo la componente 0
	Real first_=2*(s/((s-_dof[output_index]) * (s-_dof[output_index])))*stderror*trace_;
	Real GCV_der_val=first_+second_[0];
	#ifdef R_VERSION_
		Rprintf("GCV_derivative=%f\n",GCV_der_val);
	#else
		std::cout << "GCV_derivative value="<<GCV_der_val<<std::endl;
	#endif

	return GCV_der_val; //_Calcolo della derivata della GCV

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
