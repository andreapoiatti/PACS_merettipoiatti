#ifndef __REGRESSION_SKELETON_HPP__
#define __REGRESSION_SKELETON_HPP__

#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "solverdefinitions.h"
#include "optimizationData.h"
#include "carrier.h"
#include "lambda_optimizer.h"
#include "newton.h"
#include "vector_eval.h"
#include "opt_methods_factory.h"
#include "solution_builders.h"
#include "timing.h"


template<typename CarrierType>
SEXP optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
SEXP optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler & regressionData, SEXP Rmesh, OptimizationData & optimizationData)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

        timer Time_partial_JJ;
	Time_partial_JJ.start();
        Rprintf("WARNING: start taking time JJ is the matter\n");
	// Build the mixer
	MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim> regression(mesh,regressionData,optimizationData);
	regression.preapply();
	Rprintf("WARNING: JJ is the matter done\n");
        timespec T_JJ = Time_partial_JJ.stop();
        



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
{         timer Time_partial_n;
		Time_partial_n.start();

		Rprintf("WARNING: start taking time to build newton method\n");
	// Build wraper and newton method
	Function_Wrapper<Real, Real, Real, Real, EvaluationType> Fun(optim);
	typedef Function_Wrapper<Real, Real, Real, Real, EvaluationType> FunWr;

        Rprintf("WARNING: partial time after the building newton method\n");
		timespec T_n = Time_partial_n.stop();
 
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

#endif
