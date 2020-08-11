#ifndef __REGRESSION_SKELETON_H__
#define __REGRESSION_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../Regression/Include/MixedFERegression.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Carrier.h"
#include "../../Lambda_Optimization/Include/Lambda_Optimizer.h"
#include "../../Lambda_Optimization/Include/Newton.h"
#include "../../Lambda_Optimization/Include/Batch_Evaluator.h"
#include "../../Lambda_Optimization/Include/Optimization_Methods_Factory.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"

template<typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler & regressionData, OptimizationData & optimizationData, SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	MixedFERegression<InputHandler> regression(regressionData, optimizationData, mesh.num_nodes());

	regression.template preapply<ORDER,mydim,ndim, Integrator, IntegratorGaussP3, 0, 0>(mesh);

        std::pair<MatrixXr, output_Data> solution_bricks;
	//Build the carrier
	if(regression.isSV())
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			Rprintf("Areal-forced\n");
			Carrier<InputHandler,Forced,Areal>
				carrier = CarrierBuilder<InputHandler>::build_forced_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler, Forced,Areal>>(carrier);
		}
		else
		{
			Rprintf("Pointwise-forced\n");
			Carrier<InputHandler,Forced>
				carrier = CarrierBuilder<InputHandler>::build_forced_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Forced>>(carrier);
		}
	}
	else
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			Rprintf("Areal\n");
			Carrier<InputHandler,Areal>
				carrier = CarrierBuilder<InputHandler>::build_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler,Areal>>(carrier);
		}
		else
		{
			Rprintf("Pointwise\n");
			Carrier<InputHandler>
				carrier = CarrierBuilder<InputHandler>::build_plain_carrier(regressionData, regression, optimizationData);
			solution_bricks = optimizer_method_selection<Carrier<InputHandler>>(carrier);
		}
	}

 	return Solution_Builders::build_solution_plain_regression<InputHandler, ORDER, mydim, ndim>(solution_bricks.first,solution_bricks.second,mesh,regressionData);
}

template<typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_method_selection(CarrierType & carrier)
{
	// Build the optimizer				[[ TO DO factory]]
	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
	{
		Rprintf("GCV exact\n");
		GCV_Exact<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Exact<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else if(optr->get_loss_function() == "GCV" && (optr->get_DOF_evaluation() == "stochastic" || optr->get_DOF_evaluation() == "not_required"))
	{
		Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
	{
		Rprintf("Pure evaluation\n");
		auto model  = carrier.get_model();
		GCV_Exact<CarrierType, 1> optim(carrier);

		timer Time_partial;
		Time_partial.start();
		Rprintf("WARNING: start taking time\n");

		// Get the solution
		output_Data output;
		output.z_hat.resize(carrier.get_psip()->rows(),carrier.get_opt_data()->get_size_S());
		output.lambda_vec = carrier.get_opt_data()->get_lambda_S();
		MatrixXr solution;

		for(UInt j=0; j<carrier.get_opt_data()->get_size_S(); j++)
		{
			if(j==0)
			{
				MatrixXr sol = carrier.apply(carrier.get_opt_data()->get_lambda_S()[j]);
				solution.resize(sol.rows(),carrier.get_opt_data()->get_size_S());
				solution.col(j) = sol;
			}
			else
			{
				solution.col(j) = carrier.apply(carrier.get_opt_data()->get_lambda_S()[j]);
			}
			optim.combine_output_prediction(solution.topRows(solution.rows()/2).col(j),output,j);
		}

		Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                //postponed after apply in order to have betas computed
                output.betas = carrier.get_model()->getBeta();

                return {solution, output};
	}
}

template<typename EvaluationType, typename CarrierType>
std::pair<MatrixXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
{
	// Build wraper and newton method
	Function_Wrapper<Real, Real, Real, Real, EvaluationType> Fun(optim);
	typedef Function_Wrapper<Real, Real, Real, Real, EvaluationType> FunWr;

	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_criterion() == "batch")
	{
		timer Time_partial;
		Time_partial.start();
		Rprintf("WARNING: start taking time\n");
		//this will be used when batch will be correctly implemented, also for return elements
		//Eval_GCV<Real, Real, GCV_Exact<Carrier<MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim>>, 1>> eval(Fun, *(optimizationData.get_lambdas_()));
		Eval_GCV<Real, Real, EvaluationType> eval(Fun, optr->get_lambda_S());  //debugging dummy trial: working
		output_Data output = eval.Get_optimization_vectorial();

		Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		MatrixXr solution = carrier.apply(output.lambda_sol);

		output.time_partial = T.tv_sec + 1e-9*T.tv_nsec;

                //postponed after apply in order to have betas computed
                output.betas = carrier.get_model()->getBeta();

                return {solution, output};
		//Solution_Builders::GCV_batch_sol(solution, output_vec);
	}
	else // 'not_required' optimization can't enter here!! [checked in R code]
	{
		std::unique_ptr<Opt_methods<Real,Real,EvaluationType>> optim_p =
			Opt_method_factory<FunWr, Real, Real, EvaluationType>::create_Opt_method(optr->get_criterion(), Fun);

                // Compute optimal lambda
		Checker ch;
		std::vector<Real> lambda_v_;
		std::vector<Real> GCV_v_;
		Real lambda = optr->get_initial_lambda_S();
		if(lambda <=0)
		{
			lambda = -1.0;
		}

		timer Time_partial;
		Time_partial.start();
		Rprintf("WARNING: start taking time\n");

		std::pair<Real, UInt> lambda_couple = optim_p->compute(lambda, 5e-2, 40, ch, GCV_v_, lambda_v_); //era 1e-5

		Rprintf("WARNING: partial time after the optimization method\n");
		timespec T = Time_partial.stop();

		// Get the solution
		//to compute f and g hat
		MatrixXr solution = carrier.apply(lambda_couple.first);

                 //postponed after apply in order to have betas computed
		//now the last values in GCV_exact are the correct ones, related to the last iteration
	         output_Data  output = optim_p->F.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which()); //this is why F has to be public in Opt_methods
		 //the copy is necessary for the bulders outside

		return {solution, output};
	}
}

#endif
