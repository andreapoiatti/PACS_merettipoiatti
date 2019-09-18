#ifndef MATRIX_ASSEMBLER_IMP_H_
#define MATRIX_ASSEMBLER_IMP_H_

// 2D mesh implementation
template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper, const MeshHandler<ORDER, 2 ,2> & mesh,
	                   FiniteElement<Integrator, ORDER, 2, 2> & fe, SpMat & OpMat)
{
	// Fix assembler accuracy
	Real eps 	= 2.2204e-016;
	Real tolerance 	= 10 * eps;

	// Define an empty vector of triplets to be filled
	// Content of the triplet (i, j, value)
	std::vector<coeff> triplets;

  	for (auto t = 0; t < mesh.num_elements(); t++)
  	{
		// For every element in the global order [t], put it in the fe class (fetching)
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		// Needed to identify the pair of indices in the sparse matrix
		std::vector<UInt> identifiers;
		identifiers.resize(3*ORDER);
		for (auto q = 0; q < 3*ORDER; q++)
			identifiers[q] = mesh.getElement(t)[q].id();

		//localM = localMassMatrix(currentelem);

		// For each (i, j) couple of nodes compute the value in the operator
		// in each integration node and sum:
		// fourmula is sum_{l=0}^{NNODES} entity_ij(l) * det(J_T) * (weight_l * area_reference)
		for (int i = 0; i < 3*ORDER; i++)
		{
			for (int j = 0; j < 3*ORDER; j++)
			{
				Real s = 0;
				for (int l = 0; l < Integrator::NNODES; l++)
				{
					s += (oper(fe, i, j, l) * fe.getDet() * fe.getAreaReference() * Integrator::WEIGHTS[l]);
				}
			  	triplets.push_back(coeff(identifiers[i], identifiers[j], s));
			}
		}
	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);	// Give the matrix the nnodes*nnodes size

	// Insert the triplets: if thre are duplicates in the indices they are rightly summed
	OpMat.setFromTriplets(triplets.begin(), triplets.end());
	OpMat.prune(tolerance);   // suppresses al non zeros smaller than tolerance
}

template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER, 2, 2> & mesh,
	                    FiniteElement<Integrator, ORDER, 2, 2> & fe,
			    const ForcingTerm & u, VectorXr & forcingTerm)
{
	// Empty initialize the forcing vector with num_nodes size
	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for (auto t = 0; t < mesh.num_elements(); t++)
  	{
		// For every element in the global order [t], put it in the fe class (fetching)
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		// Needed to identify the pair of indices in the final vector
		std::vector<UInt> identifiers;
		identifiers.resize(3*ORDER);
		for (auto q = 0; q < 3*ORDER; q++)
			identifiers[q] = mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);

		// For each i in the nodes compute the value of the forcing term
		// in each integration node and sum:
		// fourmula is sum_{iq=0}^{NNODES} phi_i(iq)* u(iq)* det(J_T) * (weight_l * area_reference)
		for (int i = 0; i < 3*ORDER; i++)
		{
			Real s = 0;
			for (int iq = 0; iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i, iq) * u(globalIndex) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[iq]; //(*)
			}
			forcingTerm[identifiers[i]] += s;
		}
	}
}

//----------------------------------------------------------------------------//
// Surface mesh implementation

template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper, const MeshHandler<ORDER, 2, 3> & mesh,
	                   FiniteElement<Integrator, ORDER, 2, 3> & fe, SpMat & OpMat)
{
	// Fix assembler accuracy
	Real eps 	= 2.2204e-016;
	Real tolerance  = 10 * eps;

	// Define an empty vector of triplets to be filled
	// Content of the triplet (i, j, value)
	std::vector<coeff> triplets;

  	for (auto t = 0; t < mesh.num_elements(); t++)
  	{
		// For every element in the global order [t], put it in the fe class (fetching)
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		// Needed to identify the pair of indices in the sparse matrix
		std::vector<UInt> identifiers;
		identifiers.resize(3*ORDER);
		for (auto q = 0; q < 3*ORDER; q++)
			identifiers[q] = mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);

		// For each (i, j) couple of nodes compute the value in the operator
		// in each integration node and sum:
		// fourmula is sum_{l=0}^{NNODES} entity_ij(l) * det(J_T) * (weight_l * area_reference)
		// Here we have to remember that the det we get is of J_T*J_T^t, we have to take a sqrt
		for (int i = 0; i < 3*ORDER; i++)
		{
			for (int j = 0; j < 3*ORDER; j++)
			{
				Real s = 0;
				for (int l = 0; l < Integrator::NNODES; l++)
				{
					s += oper(fe, i, j, l) * std::sqrt(fe.getDet()) * fe.getAreaReference()* Integrator::WEIGHTS[l];
				}
			  	triplets.push_back(coeff(identifiers[i], identifiers[j], s));
			}
		}
	}

	UInt nnodes = mesh.num_nodes();
   	OpMat.resize(nnodes, nnodes);	// Give the matrix the nnodes*nnodes size

 	// Insert the triplets: if thre are duplicates in the indices they are rightly summed
 	OpMat.setFromTriplets(triplets.begin(), triplets.end());
 	OpMat.prune(tolerance);   // suppresses al non zeros smaller than tolerance
}

template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER, 2, 3> & mesh,
	                     FiniteElement<Integrator, ORDER, 2, 3> &  fe,
			     const ForcingTerm & u, VectorXr & forcingTerm)
{
	// Empty initialize the forcing vector with num_nodes size
	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for (auto t = 0; t < mesh.num_elements(); t++)
  	{
		// For every element in the global order [t], put it in the fe class (fetching)
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		// Needed to identify the pair of indices in the final vector
		std::vector<UInt> identifiers;
		identifiers.resize(3*ORDER);
		for (auto q = 0; q < 3*ORDER; q++)
			identifiers[q] = mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);

		// For each i in the nodes compute the value of the forcing term
		// in each integration node and sum:
		// fourmula is sum_{iq=0}^{NNODES} phi_i(iq)* u(iq)* det(J_T) * (weight_l * area_reference)
		// Here we have to remember that the det we get is of J_T*J_T^t, we have to take a sqrt
		for (int i = 0; i < 3*ORDER; i++)
		{
			Real s = 0;
			for (int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i, iq)* u(globalIndex) * std::sqrt(fe.getDet()) * fe.getAreaReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}
	}
}

//----------------------------------------------------------------------------//
// Volume mesh implementation

template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper, const MeshHandler<ORDER, 3, 3> & mesh,
	                   FiniteElement<Integrator, ORDER, 3, 3> & fe, SpMat & OpMat)
{
	// Fix assembler accuracy
	Real eps 	= 2.2204e-016;
	Real tolerance  = 10 * eps;

	// Define an empty vector of triplets to be filled
	// Content of the triplet (i, j, value)
	std::vector<coeff> triplets;

  	for (auto t = 0; t < mesh.num_elements(); t++)
  	{
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		// Needed to identify the pair of indices in the sparse matrix
		std::vector<UInt> identifiers;
		identifiers.resize(6*ORDER-2);
		for (auto q = 0; q < 6*ORDER-2; q++)
			identifiers[q] = mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);

		// For each (i, j) couple of nodes compute the value in the operator
		// in each integration node and sum:
		// fourmula is sum_{l=0}^{NNODES} entity_ij(l) * det(J_T) * (weight_l * volume_reference)
		// Here we have to remember that the det we get is of J_T*J_T^t, we have to take a sqrt
		for (int i = 0; i < 6*ORDER-2; i++)
		{
			for (int j = 0; j < 6*ORDER-2; j++)
			{
				Real s = 0;
				for (int l = 0; l < Integrator::NNODES; l++)
				{
					s += oper(fe, i, j, l) * std::sqrt(fe.getDet()) * fe.getVolumeReference() * Integrator::WEIGHTS[l];
				}
			  		triplets.push_back(coeff(identifiers[i], identifiers[j], s));
			}
		}
	}

	UInt nnodes = mesh.num_nodes();
    	OpMat.resize(nnodes, nnodes);	// Give the matrix the nnodes*nnodes size

  	// Insert the triplets: if thre are duplicates in the indices they are rightly summed
  	OpMat.setFromTriplets(triplets.begin(), triplets.end());
  	OpMat.prune(tolerance);   // suppresses al non zeros smaller than tolerance
}

template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER, 3, 3> & mesh,
	                     FiniteElement<Integrator, ORDER, 3, 3> & fe,
			     const ForcingTerm & u, VectorXr & forcingTerm)
{
	// Empty initialize the forcing vector with num_nodes size
	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for (auto t = 0; t < mesh.num_elements(); t++)
  	{
		// For every element in the global order [t], put it in the fe class (fetching)
		fe.updateElement(mesh.getElement(t));

		// Vector of vertices indices (link local to global indexing system)
		// Needed to identify the pair of indices in the final vector
		std::vector<UInt> identifiers;
		identifiers.resize(6*ORDER-2);
		for (auto q = 0; q < 6*ORDER-2; q++)
			identifiers[q] = mesh.getElement(t)[q].id();

		//localM=localMassMatrix(currentelem);

		// For each i in the nodes compute the value of the forcing term
		// in each integration node and sum:
		// fourmula is sum_{iq=0}^{NNODES} phi_i(iq)* u(iq)* det(J_T) * (weight_l * volume_reference)
		// Here we have to remember that the det we get is of J_T*J_T^t, we have to take a sqrt
		for (int i = 0; i < 6*ORDER-2; i++)
		{
			Real s = 0;
			for (int iq = 0; iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i, iq)* u(globalIndex) * std::sqrt(fe.getDet()) * fe.getVolumeReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}
	}
}

#endif
