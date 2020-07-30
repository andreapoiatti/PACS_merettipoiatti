#ifndef __DE_INITIALIZATION_SKELETON_H__
#define __DE_INITIALIZATION_SKELETON_H__

#include "../finite_element.h"
#include "../fdaPDE.h"
#include "../mesh_objects.h"
#include "../mesh.h"
#include "../mixedFERegression.h"
#include "../matrix_assembler.h"
#include "../regressionData.h"
#include "../solverdefinitions.h"

//Density Estimation
#include "../DataProblem.h"
#include "../FunctionalProblem.h"
#include "../OptimizationAlgorithm.h"
#include "../OptimizationAlgorithm_factory.h"
#include "../FEDensityEstimation.h"

template<typename Integrator, typename Integrator_noPoly, UInt ORDER, UInt mydim, UInt ndim>
SEXP DE_init_skeleton(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals,
	SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rmesh, SEXP Rsearch, const std::string & init, UInt init_fold)
{
	// Construct data problem object
	DataProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rsearch, Rmesh);

	// Construct functional problem object
	FunctionalProblem<Integrator, Integrator_noPoly, ORDER, mydim, ndim> functionalProblem(dataProblem);

	if(init == "Heat"){

		// Construct densityInit object
		std::unique_ptr<DensityInitialization<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> densityInit = make_unique<HeatProcess<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dataProblem, functionalProblem);

		// fill fInit
		std::vector<VectorXr> fInit(dataProblem.getNlambda());
		for(UInt l = 0; l < dataProblem.getNlambda(); l++){
			fInit[l] = *(densityInit-> chooseInitialization(dataProblem.getLambda(l)));
		}

		// Copy result in R memory
		SEXP result = NILSXP;
		result = PROTECT(Rf_allocVector(VECSXP, 1));
		SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, ((fInit[0])).size(), fInit.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt j = 0; j < fInit.size(); j++)
		{
			for(UInt i = 0; i < (fInit[0]).size(); i++)
				rans[i + (fInit[0]).size()*j] = (fInit[j])[i];
		}

		UNPROTECT(1);

		return(result);
	}

	else if(init=="CV"){

		// Construct densityInit object
		std::unique_ptr<Heat_CV<Integrator, Integrator_noPoly, ORDER, mydim, ndim>> densityInit = make_unique<Heat_CV<Integrator, Integrator_noPoly, ORDER, mydim, ndim>>(dataProblem, functionalProblem, init_fold);

		// fill fInit
		VectorXr fInit;
		fInit = *(densityInit->chooseInitialization(0));

		// Copy result in R memory
		SEXP result = NILSXP;
		result = PROTECT(Rf_allocVector(VECSXP, 1));
		SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, fInit.size()));

		Real *rans = REAL(VECTOR_ELT(result, 0));
		for(UInt i = 0; i < fInit.size(); i++)
		{
			rans[i] = fInit[i];
		}

		UNPROTECT(1);

		return(result);
	}
	else{

		#ifdef R_VERSION_
		Rprintf("Invalid initialization");
		#endif

		return NILSXP;
	}

}


#endif
