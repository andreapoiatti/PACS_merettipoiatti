#ifndef __AUXILIARY_SKELETON_H__
#define __AUXILIARY_SKELETON_H__

#include "../fdaPDE.h"
#include "../finite_element.h"
#include "../matrix_assembler.h"
#include "../mesh.h"

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP get_integration_points_skeleton(SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	FiniteElement<Integrator,ORDER, mydim, ndim> fe;

	SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, 2*Integrator::NNODES*mesh.num_elements()));
	for(UInt i=0; i<mesh.num_elements(); i++)
	{
		fe.updateElement(mesh.getElement(i));
		for(UInt l = 0;l < Integrator::NNODES; l++)
		{
			Point p = fe.coorQuadPt(l);
			REAL(result)[i*Integrator::NNODES + l] = p[0];
			REAL(result)[mesh.num_elements()*Integrator::NNODES + i*Integrator::NNODES + l] = p[1];
		}
	}

	UNPROTECT(1);
	return(result);
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim, typename A>
SEXP get_FEM_Matrix_skeleton(SEXP Rmesh, EOExpr<A> oper)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	SpMat AMat;
	Assembler::operKernel(oper, mesh, fe, AMat);

	//Copy result in R memory
	SEXP result;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, AMat.nonZeros() , 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, AMat.nonZeros()));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
	Real  *rans2 = REAL(VECTOR_ELT(result, 1));
	UInt i = 0;
	for (UInt k=0; k < AMat.outerSize(); ++k)
		{
			for (SpMat::InnerIterator it(AMat,k); it; ++it)
			{
				//std::cout << "(" << it.row() <<","<< it.col() <<","<< it.value() <<")\n";
				rans[i] = 1+it.row();
				rans[i + AMat.nonZeros()] = 1+it.col();
				rans2[i] = it.value();
				i++;
			}
		}
	UNPROTECT(1);
	return(result);
}

#endif
