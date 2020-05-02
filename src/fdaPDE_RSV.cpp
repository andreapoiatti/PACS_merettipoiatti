#include "fdaPDE.h"
#include "regressionSkeleton.h"

extern "C"
{
        /*!
		This function is then called from R code.
		\param Robservations an R-vector containing the values of the observations.
		\param Rdesmat an R-matrix containing the design matrix for the regression.
		\param Rmesh an R-object containg the output mesh from Trilibrary
		\param Rorder an R-integer containing the order of the approximating basis.
		\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
		\param RK an R object representing the diffusivity tensor of the model
		\param Rbeta an R object representing the advection function of the model
		\param Rc an R object representing the reaction function of the model
		\param Ru an R object representing the forcing function of the model
		\param Rcovariates an R-matrix of covariates for the regression model
		\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
		\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
				the other are automatically considered in Neumann Condition.
		\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
		\param DOF an R boolean indicating whether dofs of the model have to be computed or not
		\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
		\param Rnrealizations the number of random points used in the stochastic computation of the dofs
		\return R-vector containg the coefficients of the solution
	*/


	SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
		SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
		SEXP ROPTmethod, SEXP Rlambda, SEXP Rinitial_lambda, SEXP Rnrealizations)
	{
	    	//Set data
		RegressionDataEllipticSpaceVarying regressionData(Rlocations, Robservations, Rorder, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues);

		// Set optimization method
		OptimizationData optimizationData(ROPTmethod, Rlambda, Rinitial_lambda, Rnrealizations);

		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		if(regressionData.getOrder() == 1 && ndim==2)
			return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh, optimizationData));
		else if(regressionData.getOrder() == 2 && ndim==2)
			return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh, optimizationData));
		else if(regressionData.getOrder() == 1 && ndim==3)
			return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh, optimizationData));
		else if(regressionData.getOrder() == 2 && ndim==3)
			return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh, optimizationData));
		return(NILSXP);
	}
}
