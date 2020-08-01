#include "../../FdaPDE.h"
#include "../../Skeletons/Headers/DE_Skeleton.h"
#include "../../Skeletons/Headers/DE_Initialization_Skeleton.h"
#include "../../Regression/Headers/RegressionData.h"
#include "../../Mesh/Headers/Mesh_Objects.h"
#include "../../FE_Assemblers_Solvers/Headers/Integration.h"
#include "../../Mesh/Headers/Mesh.h"
#include "../../FE_Assemblers_Solvers/Headers/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Headers/Matrix_Assembler.h"
#include "../../Global_Utilities/Headers/Solver_Definitions.h"
//#include <chrono>

#include "../../Regression/Headers/MixedFERegression.h"

//Density Estimation
#include "../Headers/DataProblem.h"
#include "../Headers/FunctionalProblem.h"
#include "../Headers/OptimizationAlgorithm.h"
#include "../Headers/OptimizationAlgorithm_factory.h"
#include "../Headers/FEDensityEstimation.h"


extern "C" {

        //! This function manages the various options for DE-PDE algorithm
        /*!
        	This function is than called from R code.
        	\param Rdata an R-matrix containing the data.
        	\param Rmesh an R-object containg the output mesh from Trilibrary
        	\param Rorder an R-integer containing the order of the approximating basis.
        	\param Rmydim an R-integer containing the dimension of the problem we are considering.
        	\param Rndim an R-integer containing the dimension of the space in which the location are.
        	\param Rfvec an R-vector containing the the initial solution coefficients given by the user.
        	\param RheatStep an R-double containing the step for the heat equation initialization.
        	\para, RheatIter an R-integer containing the number of iterations to perfrom the heat equation initialization.
        	\param Rlambda an R-vector containing the penalization terms.
        	\param Rnfolds an R-integer specifying the number of folds for cross validation.
        	\param Rnsim an R-integer specifying the number of iterations to use in the optimization algorithm.
        	\param RstepProposals an R-vector containing the step parameters useful for the descent algotihm.
        	\param Rtol1 an R-double specifying the tolerance to use for the termination criterion based on the percentage differences.
        	\param Rtol2 an R-double specifying the tolerance to use for the termination criterion based on the norm of the gradient.
        	\param Rprint and R-integer specifying if print on console.
        	\param RstepMethod an R-string containing the method to use to choose the step in the optimization algorithm.
        	\param RdirectionMethod an R-string containing the descent direction to use in the optimization algorithm.
        	\param RpreprocessMethod an R-string containing the cross-validation method to use.
        	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).

        	\return R-list containg solutions.
        */

        SEXP Density_Estimation(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
        	 SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP RstepMethod, SEXP RdirectionMethod, SEXP RpreprocessMethod, SEXP Rsearch)
        {
        	UInt order= INTEGER(Rorder)[0];
          UInt mydim=INTEGER(Rmydim)[0];
        	UInt ndim=INTEGER(Rndim)[0];

        	std::string step_method=CHAR(STRING_ELT(RstepMethod, 0));
        	std::string direction_method=CHAR(STRING_ELT(RdirectionMethod, 0));
        	std::string preprocess_method=CHAR(STRING_ELT(RpreprocessMethod, 0));

          if(order== 1 && mydim==2 && ndim==2)
        		return(DE_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order== 2 && mydim==2 && ndim==2)
        		return(DE_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order== 1 && mydim==2 && ndim==3)
        		return(DE_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order== 2 && mydim==2 && ndim==3)
        		return(DE_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	else if(order == 1 && mydim==3 && ndim==3)
        		return(DE_skeleton<IntegratorTetrahedronP2, IntegratorGaussTetra3, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));
        	// else if(order == 1 && mydim==3 && ndim==3)
        	// 	return(DE_skeleton<IntegratorTetrahedronP2, IntegratorGaussTetra3, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, step_method, direction_method, preprocess_method));

        	return(NILSXP);
        }


          SEXP Density_Initialization(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda,
        	 SEXP Rnfolds, SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch, SEXP Rinit, SEXP Rinit_fold)
        {
        	UInt order= INTEGER(Rorder)[0];
          UInt mydim=INTEGER(Rmydim)[0];
        	UInt ndim=INTEGER(Rndim)[0];

        	UInt init_fold=INTEGER(Rinit_fold)[0];

        	std::string init=CHAR(STRING_ELT(Rinit, 0));

          if(order== 1 && mydim==2 && ndim==2)
        		return(DE_init_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order== 2 && mydim==2 && ndim==2)
        		return(DE_init_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 2>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order== 1 && mydim==2 && ndim==3)
        		return(DE_init_skeleton<IntegratorTriangleP2, IntegratorGaussTriangle3, 1, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order== 2 && mydim==2 && ndim==3)
        		return(DE_init_skeleton<IntegratorTriangleP4, IntegratorGaussTriangle3, 2, 2, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	else if(order == 1 && mydim==3 && ndim==3)
        		return(DE_init_skeleton<IntegratorTetrahedronP2, IntegratorGaussTetra3, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));
        	// else if(order == 1 && mydim==3 && ndim==3)
        	// 	return(DE_init_skeleton<IntegratorTetrahedronP2, IntegratorTetrahedronP2, 1, 3, 3>(Rdata, Rorder, Rfvec, RheatStep, RheatIter, Rlambda, Rnfolds, Rnsim, RstepProposals, Rtol1, Rtol2, Rprint, Rmesh, Rsearch, init, init_fold));

        	return(NILSXP);
        }
}
