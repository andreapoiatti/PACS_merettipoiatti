#include "../../FdaPDE.h"
#include "../../Skeletons/Headers/Auxiliary_Skeleton.h"
#include "../../Skeletons/Headers/Regression_Skeleton.h"
#include "../../Skeletons/Headers/Regression_Skeleton_Time.h"
#include "../../Skeletons/Headers/GAM_Skeleton.h"
#include "../Headers/RegressionData.h"
#include "../../FE_Assemblers_Solvers/Headers/Integration.h"
#include "../../Lambda_Optimization/Headers/Optimization_Data.h"

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
        	\param GCV an R boolean indicating whether dofs of the model have to be computed or not
        	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when GCV is TRUE, can be either 1 (exact) or 2 (stochastic)
        	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
        	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
        	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.

        	\return R-vector containg the coefficients of the solution
        */
        SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
                SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg,
                SEXP Rsearch, SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune)
        {
                //Set data
        	RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch);
                OptimizationData optimizationData(Roptim, Rlambda, Rnrealizations, Rseed, RDOF_matrix, Rtune);

        	UInt mydim = INTEGER(Rmydim)[0];
        	UInt ndim = INTEGER(Rndim)[0];

        	if(regressionData.getOrder() == 1 && ndim==2)
        		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 2>(regressionData, optimizationData, Rmesh));
        	else if(regressionData.getOrder() == 2 && ndim==2)
        		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 2>(regressionData, optimizationData, Rmesh));
        	else if(regressionData.getOrder() == 1 && ndim==3)
        		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 3>(regressionData, optimizationData, Rmesh));
        	else if(regressionData.getOrder() == 2 && ndim==3)
        		return(regression_skeleton<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 3>(regressionData, optimizationData, Rmesh));
        	return(NILSXP);
        }

        /*!
                This function is then called from R code.
                \param Rlocations an R-matrix containing the spatial locations of the observations
                \param Rtime_locations an R-vector containing the temporal locations of the observations
                \param Robservations an R-vector containing the values of the observations.
                \param Rmesh an R-object containing the spatial mesh
                \param Rmesh_time an R-vector containing the temporal mesh
                \param Rorder an R-integer containing the order of the approximating basis in space.
                \param Rmydim an R-integer specifying if the mesh nodes lie in R^2 or R^3
                \param Rndim  an R-integer specifying if the "local dimension" is 2 or 3
                \param RlambdaS an R-double containing the penalization term of the empirical evidence respect to the prior one.
                \param RlambdaT an R-double containing the penalization term of the empirical evidence respect to the prior one.
                \param RK an R object representing the diffusivity tensor of the model
                \param Rbeta an R object representing the advection function of the model
                \param Rc an R object representing the reaction function of the model
                \param Ru an R object representing the forcing function of the model
                \param Rcovariates an R-matrix of covariates for the regression model
                \param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
                \param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
                                the other are automatically considered in Neumann Condition.
                \param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
                \param Rflag_mass an R-integer that in case of separable problem specifies whether to use mass discretization or identity discretization
                \param Rflag_parabolic an R-integer specifying if the problem is parabolic or separable
                \param Ric an R-vector containing the initial condition needed in case of parabolic problem
                \param GCV an R-integer indicating if the GCV has to be computed or not
                \param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
                \param DOF an R boolean indicating whether dofs of the model have to be computed or not
                \param RDOF_matrix a R-matrix containing the dofs (for every combination of the values in RlambdaS and RlambdaT) if they are already known from precedent computations
                \param Rnrealizations the number of random points used in the stochastic computation of the dofs
                \param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
                \param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.

                \return R-vector containg the coefficients of the solution
        */
        SEXP regression_PDE_space_varying_time(SEXP Rlocations, SEXP RbaryLocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
		SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues,  SEXP RincidenceMatrix, SEXP RarealDataAvg,
                SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP Rsearch, SEXP Roptim, SEXP Rlambda_S, SEXP Rlambda_T, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune)
        {
                //Set data
                RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Rtime_locations, Robservations, Rorder, RK, Rbeta, Rc, Ru,
                         Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rflag_mass, Rflag_parabolic, Ric, Rsearch);
                OptimizationData optimizationData(Roptim, Rlambda_S, Rlambda_T, Rnrealizations, Rseed, RDOF_matrix, Rtune);

                UInt mydim = INTEGER(Rmydim)[0];
                UInt ndim = INTEGER(Rndim)[0];

                if(regressionData.getOrder() == 1 && ndim==2)
                        return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, optimizationData, Rmesh, Rmesh_time));
                else if(regressionData.getOrder() == 2 && ndim==2)
                        return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, optimizationData, Rmesh, Rmesh_time));
                else if(regressionData.getOrder() == 1 && ndim==3)
                        return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, optimizationData, Rmesh, Rmesh_time));
                else if(regressionData.getOrder() == 2 && ndim==3)
                        return(regression_skeleton_time<RegressionDataEllipticSpaceVarying,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, optimizationData, Rmesh, Rmesh_time));
                return(NILSXP);
        }

        /*!
        	This function is then called from R code.
        	\param Rlocations an R-matrix containing the location of the observations.
        	\param Robservations an R-vector containing the values of the observations.
        	\param Rmesh an R-object containg the output mesh from Trilibrary.
        	\param Rorder an R-integer containing the order of the approximating basis.
        	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
        	\param RK an R object representing the diffusivity tensor of the model
        	\param Rbeta an R object representing the advection function of the model
        	\param Rc an R object representing the reaction function of the model
        	\param Ru an R object representing the forcing function of the model
        	\param Rcovariates an R-matrix of covariates for the regression model.
        	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data.
        	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
        			the other are automatically considered in Neumann Condition.
        	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices.
        	\param DOF an R boolean indicating whether dofs of the model have to be computed or not.
        	\param GCV an R boolean indicating whether GCV of the model have to be computed or not.
        	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic).
        	\param Rnrealizations the number of random points used in the stochastic computation of the dofs.
        	\param Rfamily Denotes the distribution of the data, within the exponential family.
        	\param Rmax_num_iteration Maximum number of steps run in the PIRLS algorithm, set to 15 by default.
        	\param Rtreshold an R-double used for arresting FPIRLS algorithm. Algorithm stops when two successive iterations lead to improvement in penalized log-likelihood smaller than threshold.
        	\param Rmu0 Initial value of the mean (natural parameter). There exists a default value for each familiy
        	\param RscaleParam If necessary and not supplied, the scale parameter \phi is estimated. See method.phi for details.
        	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
        	\param Rtune a R-double, Tuning parameter used for the estimation of GCV. called 'GCV.inflation.factor' in R code.
        	\param RarealDataAvg an R boolean indicating whether the areal data are averaged or not.

        	\return R-vector containg the outputs.
        */
          SEXP gam_PDE_space_varying(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
          	SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg,
          	SEXP Rfamily, SEXP Rmax_num_iteration, SEXP Rtreshold, SEXP Rmu0, SEXP RscaleParam, SEXP Rsearch,
                SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune)
        {
        	GAMDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch, Rmax_num_iteration, Rtreshold);
                OptimizationData optimizationData(Roptim, Rlambda, Rnrealizations, Rseed, RDOF_matrix, Rtune);

        	UInt mydim = INTEGER(Rmydim)[0];
        	UInt ndim = INTEGER(Rndim)[0];

          	std::string family = CHAR(STRING_ELT(Rfamily,0));

                if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
                	return(GAM_skeleton<GAMDataEllipticSpaceVarying,IntegratorTriangleP2, 1, 2, 2>(regressionData, optimizationData, Rmesh, Rmu0 , family, RscaleParam));
                else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
                	return(GAM_skeleton<GAMDataEllipticSpaceVarying,IntegratorTriangleP4, 2, 2, 2>(regressionData, optimizationData, Rmesh, Rmu0, family, RscaleParam));
                return(R_NilValue);
        }

        //! A utility, not used for system solution, may be used for debugging
        SEXP get_FEM_PDE_space_varying_matrix(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru,
                           SEXP Rcovariates, SEXP RBCIndices, SEXP RBCValues, SEXP RincidenceMatrix, SEXP RarealDataAvg, SEXP Rsearch)
        {
                RegressionDataEllipticSpaceVarying regressionData(Rlocations, RbaryLocations, Robservations, Rorder, RK, Rbeta, Rc, Ru, Rcovariates, RBCIndices, RBCValues, RincidenceMatrix, RarealDataAvg, Rsearch);

                //UInt mydim = INTEGER(Rmydim)[0];
                UInt ndim = INTEGER(Rndim)[0];

                typedef EOExpr<Mass>  ETMass;  Mass  EMass;  ETMass  mass(EMass);
                typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
                typedef EOExpr<Grad>  ETGrad;  Grad  EGrad;  ETGrad  grad(EGrad);

                const Reaction & c = regressionData.getC();
                const Diffusivity & K = regressionData.getK();
                const Advection & beta = regressionData.getBeta();

                if(regressionData.getOrder()==1 && ndim==2)
                        return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
                if(regressionData.getOrder()==2 && ndim==2)
                        return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
                return(NILSXP);
        }
}
