#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "solverdefinitions.h"
#include "optimizationData.h"
#include "carrier.h"
#include "lambda_optimizer.h"
#include "newton.h"
#include "vector_eval.h"
#include "opt_methods_factory.h"
#include "solution_builders.h"
#include "mixedFEFPCA.h"
#include "mixedFERegression.h"
#include "mixedFEFPCAfactory.h"
#include "auxiliarySkeleton.h"

extern "C"
{
	//! A function required for anysotropic and nonstationary regression (only 2D)
	/*!
	    \return points where the PDE space-varying params are evaluated in the R code
	*/
	SEXP get_integration_points(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
	{
		//Declare pointer to access data from C++
		int order = INTEGER(Rorder)[0];

		//Get mydim and ndim
		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];
	//Not implemented for ndim==3
	    if(order == 1 && ndim ==2)
	    	return(get_integration_points_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh));
	    else if(order == 2 && ndim==2)
	    	return(get_integration_points_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh));
	    return(NILSXP);
	}


	//! A utility, not used for system solution, may be used for debugging

	SEXP get_FEM_mass_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
	{
		int order = INTEGER(Rorder)[0];

		//Get mydim and ndim
		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);

	    if(order==1 && ndim==2)
	    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, mass));
		if(order==2 && ndim==2)
			return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, mass));
		return(NILSXP);
	}

	//! A utility, not used for system solution, may be used for debugging
	SEXP get_FEM_stiff_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
	{
		int order = INTEGER(Rorder)[0];

		//Get mydim and ndim
		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		typedef EOExpr<Stiff> ETMass;   Stiff EStiff;   ETMass stiff(EStiff);

	    if(order==1 && ndim==2)
	    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, stiff));
		if(order==2 && ndim==2)
			return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, stiff));
		return(NILSXP);
	}

	/*
	//! A utility, not used for system solution, may be used for debugging
	SEXP get_FEM_PDE_matrix(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc,
					   SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations)
	{
		RegressionDataElliptic regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, DOF, RGCVmethod, Rnrealizations);

		//Get mydim and ndim
		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Real& c = regressionData.getC();
		const Eigen::Matrix<Real,2,2>& K = regressionData.getK();
		const Eigen::Matrix<Real,2,1>& beta = regressionData.getBeta();

	    if(regressionData.getOrder()==1 && ndim==2)
	    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
		if(regressionData.getOrder()==2 && ndim==2)
			return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
		return(NILSXP);
	}

	//! A utility, not used for system solution, may be used for debugging
	SEXP get_FEM_PDE_space_varying_matrix(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
			   SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations)
	{
		RegressionDataEllipticSpaceVarying regressionData(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, DOF, RGCVmethod, Rnrealizations);

		//Get mydim and ndim
		//UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = regressionData.getC();
		const Diffusivity& K = regressionData.getK();
		const Advection& beta = regressionData.getBeta();

	    if(regressionData.getOrder()==1 && ndim==2)
	    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
		if(regressionData.getOrder()==2 && ndim==2)
			return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
		return(NILSXP);
	}
	*/
}
