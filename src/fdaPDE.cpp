
#define R_VERSION_

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
//#include <chrono>

#include "mixedFEFPCA.h"
#include "mixedFERegression.h"
#include "mixedFEFPCAfactory.h"

template<typename CarrierType>
SEXP optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
SEXP optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler & regressionData, SEXP Rmesh, OptimizationData & optimizationData)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	// Build the mixer
	MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim> regression(mesh,regressionData,optimizationData);
	regression.preapply();

	//Build the carrier
	if(regression.check_is_space_varying())
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			Rprintf("Areal-forced\n");
			Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>, Forced, Areal>
				carrier = CarrierBuilder<InputHandler,MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>::build_forced_areal_carrier(regressionData, regression, optimizationData);
			return optimizer_method_selection<Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>, Forced, Areal>>(carrier);
		}
		else
		{
			Rprintf("Pointwise-forced\n");
			Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>, Forced>
				carrier = CarrierBuilder<InputHandler,MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>::build_forced_carrier(regressionData, regression, optimizationData);
			return optimizer_method_selection<Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>, Forced>>(carrier);
		}
	}
	else
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			Rprintf("Areal\n");
			Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>, Areal>
				carrier = CarrierBuilder<InputHandler,MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>::build_areal_carrier(regressionData, regression, optimizationData);
			return optimizer_method_selection<Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>, Areal>>(carrier);
		}
		else
		{
			Rprintf("Pointwise\n");
			Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>
				carrier = CarrierBuilder<InputHandler,MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>::build_plain_carrier(regressionData, regression, optimizationData);
			return optimizer_method_selection<Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>>(carrier);
		}
	}
}

template<typename CarrierType>
SEXP optimizer_method_selection(CarrierType & carrier)
{
	// Build the optimizer				[[ TO DO factory]]
	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_method_() == "gcv" && optr->get_evaluation_() == "exact")
	{
		Rprintf("GCV exact\n");
		GCV_Exact<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Exact<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else if(optr->get_method_() == "gcv" && optr->get_evaluation_() == "stochastic")
	{
		Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else // E.g. K-FOLD CV
	{
		Rprintf("Error, no other method implemented\n");
		return NILSXP;
	}
}

template<typename EvaluationType, typename CarrierType>
SEXP optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
	// Build wraper and newton method
	Function_Wrapper<Real, Real, Real, Real, EvaluationType> Fun(optim);
	typedef Function_Wrapper<Real, Real, Real, Real, EvaluationType> FunWr;

	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_criterion_() == "batch")
	{       timer Time_partial;
	        Time_partial.start();
		Rprintf("WARNING: start taking time\n");
		//this will be used when batch will be correctly implemented, also for return elements
		//Eval_GCV<Real, Real, GCV_Exact<Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>, 1>> eval(Fun, *(optimizationData.get_lambdas_()));
		Eval_GCV<Real, Real, EvaluationType> eval(Fun, *(optr->get_lambdas_()));  //debugging dummy trial: working
		output_Data_opt output_vec = eval.Get_optimization_vectorial();

		Rprintf("WARNING: partial time after the batch method\n");
		timespec T = Time_partial.stop();

		//to compute f and g hat
		carrier.get_tracep()->apply(output_vec.lambda_opt);

		// Get the solution
		VectorXr solution = carrier.get_tracep()->getSolution();

                return Solution_builders::GCV_batch_sol(solution, output_vec, T);
	}
	else
	{
		std::unique_ptr<Opt_methods<Real,Real,EvaluationType>> optim_p =
			Opt_method_factory<FunWr, Real, Real, EvaluationType>::create_Opt_method(optr->get_criterion_(), Fun);

		// Compute optimal lambda
		Checker ch;
		Real lambda = optr->get_initial_lambda_();
		if(lambda <=0)
		{
			// [[ TO DO]] automatic clever method to be implemented
			lambda = 0.01;
		}

		timer Time_partial;
		Time_partial.start();

		Rprintf("WARNING: start taking time\n");

		std::pair<Real, UInt> lambda_couple = optim_p->compute(lambda, 1e-5, 40, ch);

		Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		//now the last values in GCV_exact are the correct ones, related to the last iteration
		const output_Data & output = optim_p->F.get_output(lambda_couple, T); //this is why F has to be public in Opt_methods
		carrier.get_tracep()->apply(lambda_couple.first);

		// Get the solution
		//to compute f and g hat
		VectorXr solution = carrier.get_tracep()->getSolution();

		return Solution_builders::GCV_Newton_sol(solution, output);  // [[TO DO make this a template according to methods]]
	}
}



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

extern "C"
{
	void fill_optimization_data (OptimizationData & optimizationData, SEXP ROPTmethod, SEXP Rlambdas, SEXP Rinitial_lambda, SEXP Rnrealizations)
	{
		UInt criterion = INTEGER(ROPTmethod)[0];
		if(criterion == 2)
		{
			optimizationData.set_criterion_("newton_fd");
		}
		else if(criterion == 1)
		{
			optimizationData.set_criterion_("newton");
		}
		else if(criterion == 0)
			optimizationData.set_criterion_("batch");

		UInt method = INTEGER(ROPTmethod)[1];
		if(method == 0)
		{
			optimizationData.set_method_("gcv");
		}

		UInt stochastic = INTEGER(ROPTmethod)[2];
		if(stochastic == 1)
		{
			optimizationData.set_evaluation_("stochastic");
			optimizationData.set_nrealizations_(INTEGER(Rnrealizations)[0]);
		}
		else
		{
			optimizationData.set_evaluation_("exact");
		}

		if(criterion == 0)
		{       std::cout<<"Valori"<<std::endl;
			UInt n_lambdas_ = Rf_length(Rlambdas);
			std::vector<Real> lambdas_;
			lambdas_.resize(n_lambdas_);

			for(UInt i=0; i<n_lambdas_; ++i)
			{
				lambdas_[i] = REAL(Rlambdas)[i];
			}

			optimizationData.set_lambdas_(lambdas_);
		}
		else if(criterion == 1 || criterion==2)
		{
			UInt initialization = INTEGER(ROPTmethod)[3];
			if(initialization == 1)
				optimizationData.set_initial_lambda_(REAL(Rinitial_lambda)[0]);
		}
	}

	//! This function manages the various options for Spatial Regression, Sangalli et al version
	/*!
		This function is then called from R code.
		\param Robservations an R-vector containing the values of the observations.
		\param Rdesmat an R-matrix containing the design matrix for the regression.
		\param Rmesh an R-object containg the output mesh from Trilibrary
		\param Rorder an R-integer containing the order of the approximating basis.
		\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
		\param Rcovariates an R-matrix of covariates for the regression model
		\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
		\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
				the other are automatically considered in Neumann Condition.
		\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
		\param ROPTmethod an R-integer indicating the method to use to compute the dofs
		\param Rnrealizations the number of random points used in the stochastic computation of the dofs
		\return R-vector containg the coefficients of the solution
	*/

	SEXP regression_Laplace(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim,
	        SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
	        SEXP ROPTmethod, SEXP Rlambda, SEXP Rinitial_lambda, SEXP Rnrealizations)
	{
	        //Set input data
	        RegressionData regressionData(Rlocations, Robservations, Rorder, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues);

	        // Set optimization method
		OptimizationData optimizationData;
		fill_optimization_data(optimizationData, ROPTmethod, Rlambda, Rinitial_lambda, Rnrealizations);

	        UInt mydim = INTEGER(Rmydim)[0];
	        UInt ndim = INTEGER(Rndim)[0];

	        if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
	        	return(regression_skeleton<RegressionData,IntegratorTriangleP2,1,2,2>(regressionData, Rmesh, optimizationData));
	        else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
	        	return(regression_skeleton<RegressionData,IntegratorTriangleP4,2,2,2>(regressionData, Rmesh, optimizationData));
	        else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
	        	return(regression_skeleton<RegressionData,IntegratorTriangleP2,1,2,3>(regressionData, Rmesh, optimizationData));
	        else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
	        	return(regression_skeleton<RegressionData,IntegratorTriangleP4,2,2,3>(regressionData, Rmesh, optimizationData));
	        else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
	        	return(regression_skeleton<RegressionData,IntegratorTetrahedronP2,1,3,3>(regressionData, Rmesh, optimizationData));
		return(NILSXP);

	}

	/*!
		This function is then called from R code.
		\param Robservations an R-vector containing the values of the observations.
		\param Rdesmat an R-matrix containing the design matrix for the regression.
		\param Rmesh an R-object containg the output mesh from Trilibrary
		\param Rorder an R-integer containing the order of the approximating basis.
		\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
		\param RK an R-matrix representing the diffusivity matrix of the model
		\param Rbeta an R-vector representing the advection term of the model
		\param Rc an R-double representing the reaction term of the model
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

	SEXP regression_PDE(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
		SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
		SEXP ROPTmethod, SEXP Rlambda, SEXP Rinitial_lambda, SEXP Rnrealizations)
	{
		RegressionDataElliptic regressionData(Rlocations, Robservations, Rorder, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues);

		// Set optimization method
		OptimizationData optimizationData;
		fill_optimization_data(optimizationData, ROPTmethod, Rlambda, Rinitial_lambda, Rnrealizations);

		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		if(regressionData.getOrder() == 1 && ndim==2)
			return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh, optimizationData));
		else if(regressionData.getOrder() == 2 && ndim==2)
			return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh, optimizationData));
		else if(regressionData.getOrder() == 1 && ndim==3)
			return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh, optimizationData));
		else if(regressionData.getOrder() == 2 && ndim==3)
			return(regression_skeleton<RegressionDataElliptic,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh, optimizationData));
		return(NILSXP);
	}

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
		OptimizationData optimizationData;
		fill_optimization_data(optimizationData, ROPTmethod, Rlambda, Rinitial_lambda, Rnrealizations);

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

// DEBUGGING UTILITIES TO BE INTEGRATED IN THE FINAL SYNTAX
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


	//! This function manages the various options for SF-PCA
	/*!
		This function is than called from R code.
		\param Rdatamatrix an R-matrix containing the datamatrix of the problem.
		\param Rlocations an R-matrix containing the location of the observations.
		\param Rmesh an R-object containg the output mesh from Trilibrary
		\param Rorder an R-integer containing the order of the approximating basis.
		\param RincidenceMatrix an R-matrix representing the incidence matrix defining regions in the model with areal data
		\param Rmydim an R-integer containing the dimension of the problem we are considering.
		\param Rndim an R-integer containing the dimension of the space in which the location are.
		\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
		\param RnPC an R-integer specifying the number of principal components to compute.
		\param Rvalidation an R-string containing the method to use for the cross-validation of the penalization term lambda.
		\param RnFolds an R-integer specifying the number of folds to use if K-Fold cross validation method is chosen.
		\param RGCVmethod an R-integer specifying if the GCV computation has to be exact(if = 1) or stochastic (if = 2).
		\param Rnrealizations an R-integer specifying the number of realizations to use when computing the GCV stochastically.

		\return R-vector containg the coefficients of the solution
	*/
	SEXP Smooth_FPCA(SEXP Rlocations, SEXP Rdatamatrix, SEXP Rmesh, SEXP Rorder, SEXP RincidenceMatrix, SEXP Rmydim, SEXP Rndim, SEXP Rlambda, SEXP RnPC, SEXP Rvalidation, SEXP RnFolds, SEXP RGCVmethod, SEXP Rnrealizations){
	//Set data

		FPCAData fPCAdata(Rlocations, Rdatamatrix, Rorder, RincidenceMatrix, Rlambda, RnPC, RnFolds, RGCVmethod, Rnrealizations);

		UInt mydim=INTEGER(Rmydim)[0];
		UInt ndim=INTEGER(Rndim)[0];

		std::string validation=CHAR(STRING_ELT(Rvalidation,0));

		if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==2)
			return(FPCA_skeleton<IntegratorTriangleP2, 1, 2, 2>(fPCAdata, Rmesh, validation));
		else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==2)
			return(FPCA_skeleton<IntegratorTriangleP4, 2, 2, 2>(fPCAdata, Rmesh, validation));
		else if(fPCAdata.getOrder() == 1 && mydim==2 && ndim==3)
			return(FPCA_skeleton<IntegratorTriangleP2, 1, 2, 3>(fPCAdata, Rmesh, validation));
		else if(fPCAdata.getOrder() == 2 && mydim==2 && ndim==3)
			return(FPCA_skeleton<IntegratorTriangleP4, 2, 2, 3>(fPCAdata, Rmesh, validation));
		else if(fPCAdata.getOrder() == 1 && mydim==3 && ndim==3)
			return(FPCA_skeleton<IntegratorTetrahedronP2, 1, 3, 3>(fPCAdata, Rmesh, validation));
		return(NILSXP);
		 }
}
