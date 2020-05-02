#ifndef __FPCA_SKELETON_HPP__
#define __FPCA_SKELETON_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "solverdefinitions.h"
#include "solution_builders.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "mixedFEFPCA.h"
#include "mixedFERegression.h"
#include "mixedFEFPCAfactory.h"

template<typename Integrator,UInt ORDER, UInt mydim, UInt ndim>
SEXP FPCA_skeleton(FPCAData &fPCAData, SEXP Rmesh, std::string validation)
{

	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	std::unique_ptr<MixedFEFPCABase<Integrator, ORDER, mydim, ndim>> fpca = MixedFEFPCAfactory<Integrator, ORDER, mydim, ndim>::createFPCAsolver(validation, mesh, fPCAData);

	fpca->apply();

	const std::vector<VectorXr>& loadings = fpca->getLoadingsMat();
	const std::vector<VectorXr>& scores = fpca->getScoresMat();
	const std::vector<Real>& lambdas = fpca->getLambdaPC();
	const std::vector<Real>& variance_explained = fpca->getVarianceExplained();
	const std::vector<Real>& cumsum_percentage = fpca->getCumulativePercentage();
	const std::vector<Real>& var = fpca->getVar();

	//Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 7));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, loadings[0].size(), loadings.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, scores[0].size(), scores.size()));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, lambdas.size()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, variance_explained.size()));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, cumsum_percentage.size()));
	SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, var.size()));
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < loadings.size(); j++)
	{
		for(UInt i = 0; i < loadings[0].size(); i++)
			rans[i + loadings[0].size()*j] = loadings[j][i];
	}

	Real *rans1 = REAL(VECTOR_ELT(result, 1));
	for(UInt j = 0; j < scores.size(); j++)
	{
		for(UInt i = 0; i < scores[0].size(); i++)
			rans1[i + scores[0].size()*j] = scores[j][i];
	}

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < lambdas.size(); i++)
	{
		rans2[i] = lambdas[i];
	}

	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < variance_explained.size(); i++)
	{
		rans3[i] = variance_explained[i];
	}

	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < cumsum_percentage.size(); i++)
	{
		rans4[i] = cumsum_percentage[i];
	}
	Real *rans5 = REAL(VECTOR_ELT(result, 5));
	for(UInt i = 0; i < var.size(); i++)
	{
		rans5[i] = var[i];
	}

	UNPROTECT(1);

	return(result);
}

#endif
