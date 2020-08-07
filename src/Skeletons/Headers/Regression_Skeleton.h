#ifndef __REGRESSION_SKELETON_H__
#define __REGRESSION_SKELETON_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Headers/Mesh.h"
#include "../../Regression/Headers/MixedFERegression.h"
#include "../../Lambda_Optimization/Headers/Optimization_Data.h"
#include "../../Lambda_Optimization/Headers/Carrier.h"
#include "../../Lambda_Optimization/Headers/Lambda_Optimizer.h"
#include "../../Lambda_Optimization/Headers/Newton.h"
#include "../../Lambda_Optimization/Headers/Batch_Evaluator.h"
#include "../../Lambda_Optimization/Headers/Optimization_Methods_Factory.h"
//#include "../../Lambda_Optimization/Headers/Solution_Builders.h"

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP build_sol(const VectorXr & solution, const output_Data & output, const MeshHandler<ORDER, mydim, ndim> & mesh , const InputHandler & regressionData )
{

        MatrixXv beta;
        if(regressionData.getCovariates()->rows()==0)
        {
        	beta.resize(1,1);
        	beta(0,0).resize(1);
        	beta(0,0)(0) = 10e20;
        }
        else
        	 beta = output.betas;
      UInt code_string;
      if (output.content=="full_optimization")
              code_string=0;
      else if (output.content=="full_dof_batch")
                code_string=1;
           else code_string=2;

      const MatrixXr & barycenters = regressionData.getBarycenters();
      const VectorXi & elementIds = regressionData.getElementIds();

       //Copy result in R memory
       SEXP result = NILSXP;
       result = PROTECT(Rf_allocVector(VECSXP, 22));
       SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, solution.size()));
       Real *rans = REAL(VECTOR_ELT(result, 0));
       for(UInt j = 0; j < solution.size(); j++)  //[TO DO ] //sono le f_hat e g_hat, si potrebbe rimuovere, cambiando la chiamata da R in  smooth.FEM.basis
       {
               rans[j] = solution[j];
       }


       UInt size_z=output.z_hat.size();
       SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, size_z));
       rans = REAL(VECTOR_ELT(result, 1));

       for(UInt j = 0; j < size_z; j++)
       {
               rans[j] = output.z_hat[j];
       }


       //Rprintf("Hey doc,  %f %f %f\n", output.z_hat[0], output.z_hat[1], output.z_hat[3]);
       SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 2));
       rans[0] = output.rmse;

       SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, 1));
       rans= REAL(VECTOR_ELT(result, 3));
       rans[0] = output.sigma_hat_sq;

       SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 4));
       rans[0] = output.lambda_sol;

       SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1));
       UInt *rans1 = INTEGER(VECTOR_ELT(result, 5));
       rans1[0] = output.lambda_pos;

       SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 6));
       rans[0] = output.GCV_opt;

       SET_VECTOR_ELT(result, 7, Rf_allocVector(INTSXP, 1));
        UInt * rans2 = INTEGER(VECTOR_ELT(result, 7));
       rans2[0] = output.n_it;

       SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 8));
       rans[0] = output.termination;

       SET_VECTOR_ELT(result, 9, Rf_allocVector(INTSXP, 1));
       UInt *rans3 = INTEGER(VECTOR_ELT(result, 9));
       rans3[0] = code_string;

       UInt size_dof=output.dof.size();
       SET_VECTOR_ELT(result, 10, Rf_allocVector(REALSXP, size_dof));
       rans = REAL(VECTOR_ELT(result, 10));






       /*for(UInt j = 0; j < size_dof; j++)
       {
               rans[j] = output.dof[j];
       }




              for(UInt j = 0; j < output.lambda_vec.size(); j++)
              {
                   std::cout<< output.lambda_vec[j]<<std::endl;
              }

       std::cout<<output.GCV_evals.size()<<std::endl;
       std::cout<<"sono qui"<<std::endl;
       UInt a;
       std::cin>>a;
*/

       UInt size_lambda=output.lambda_vec.size();
       SET_VECTOR_ELT(result, 11, Rf_allocVector(REALSXP, size_lambda));
       rans = REAL(VECTOR_ELT(result, 11));

       for(UInt j = 0; j < size_lambda; j++)
       {
               rans[j] = output.lambda_vec[j];
       }


       UInt size_vec=output.GCV_evals.size();
       SET_VECTOR_ELT(result, 12, Rf_allocVector(REALSXP, size_vec));
       rans = REAL(VECTOR_ELT(result, 12));
       for(UInt j = 0; j < size_vec; j++)
       {
               rans[j] = output.GCV_evals[j];
       }


       SET_VECTOR_ELT(result, 13, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 13));
       rans[0] = output.time_partial;




       //! Copy betas
       SET_VECTOR_ELT(result, 14, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));
       Real *rans4 = REAL(VECTOR_ELT(result, 14));
       for(UInt j = 0; j < beta.size(); j++)
       {
       	for(UInt i = 0; i < beta(0).size(); i++)
       		rans4[i + beta(0).size()*j] = beta(j)(i);
       }

       //SEND TREE INFORMATION TO R
       SET_VECTOR_ELT(result, 15, Rf_allocVector(INTSXP, 1)); //tree_header information
       int *rans5 = INTEGER(VECTOR_ELT(result, 15));
       rans5[0] = mesh.getTree().gettreeheader().gettreelev();

       SET_VECTOR_ELT(result, 16, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
       Real *rans6 = REAL(VECTOR_ELT(result, 16));
       for(UInt i = 0; i < ndim*2; i++)
       	rans6[i] = mesh.getTree().gettreeheader().domainorig(i);

       SET_VECTOR_ELT(result, 17, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
       Real *rans7 = REAL(VECTOR_ELT(result, 17));
       for(UInt i = 0; i < ndim*2; i++)
       	rans7[i] = mesh.getTree().gettreeheader().domainscal(i);


       UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
       SET_VECTOR_ELT(result, 18, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
       int *rans8 = INTEGER(VECTOR_ELT(result, 18));
       for(UInt i = 0; i < num_tree_nodes; i++)
       		rans8[i] = mesh.getTree().gettreenode(i).getid();

       for(UInt i = 0; i < num_tree_nodes; i++)
       		rans8[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

       for(UInt i = 0; i < num_tree_nodes; i++)
       		rans8[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

       SET_VECTOR_ELT(result, 19, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
       Real *rans9 = REAL(VECTOR_ELT(result, 19));
       for(UInt j = 0; j < ndim*2; j++)
       {
       	for(UInt i = 0; i < num_tree_nodes; i++)
       		rans9[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
       }

       //SEND BARYCENTER INFORMATION TO R
       SET_VECTOR_ELT(result, 20, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
       int *rans10 = INTEGER(VECTOR_ELT(result, 20));
       for(UInt i = 0; i < elementIds.rows(); i++)
       	rans10[i] = elementIds(i);

       SET_VECTOR_ELT(result, 21, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
       Real *rans11 = REAL(VECTOR_ELT(result, 21));
       for(UInt j = 0; j < barycenters.cols(); j++)
       {
       	for(UInt i = 0; i < barycenters.rows(); i++)
       		rans11[i + barycenters.rows()*j] = barycenters(i,j);
       }




       UNPROTECT(1);

       return(result);
}


template<typename CarrierType>
std::pair<VectorXr, output_Data> optimizer_method_selection(CarrierType & carrier);
template<typename EvaluationType, typename CarrierType>
std::pair<VectorXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier);

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler & regressionData, OptimizationData & optimizationData, SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	MixedFERegression<InputHandler> regression(regressionData, optimizationData, mesh.num_nodes());

	regression.template preapply<ORDER,mydim,ndim, Integrator, IntegratorGaussP3, 0, 0>(mesh);

std::pair<VectorXr, output_Data> solution_bricks;
	//Build the carrier
	if(regression.isSV())
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			Rprintf("Areal-forced\n");
			Carrier<InputHandler,Forced,Areal>
				carrier = CarrierBuilder<InputHandler>::build_forced_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks=optimizer_method_selection<Carrier<InputHandler, Forced,Areal>>(carrier);
		}
		else
		{
			Rprintf("Pointwise-forced\n");
			Carrier<InputHandler,Forced>
				carrier = CarrierBuilder<InputHandler>::build_forced_carrier(regressionData, regression, optimizationData);
			solution_bricks=optimizer_method_selection<Carrier<InputHandler,Forced>>(carrier);
		}
	}
	else
	{
		if(regressionData.getNumberOfRegions()>0)
		{
			Rprintf("Areal\n");
			Carrier<InputHandler,Areal>
				carrier = CarrierBuilder<InputHandler>::build_areal_carrier(regressionData, regression, optimizationData);
			solution_bricks=optimizer_method_selection<Carrier<InputHandler,Areal>>(carrier);
		}
		else
		{
			Rprintf("Pointwise\n");
			Carrier<InputHandler>
				carrier = CarrierBuilder<InputHandler>::build_plain_carrier(regressionData, regression, optimizationData);
			solution_bricks=optimizer_method_selection<Carrier<InputHandler>>(carrier);
		}
	}

return build_sol<InputHandler, ORDER, mydim, ndim>(solution_bricks.first, solution_bricks.second, mesh, regressionData);

}

template<typename CarrierType>
std::pair<VectorXr, output_Data> optimizer_method_selection(CarrierType & carrier)
{
	// Build the optimizer				[[ TO DO factory]]
	const OptimizationData * optr = carrier.get_opt_data();
	if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "exact")
	{
		Rprintf("GCV exact\n");
		GCV_Exact<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Exact<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "stochastic")
	{
		Rprintf("GCV stochastic\n");
		GCV_Stochastic<CarrierType, 1> optim(carrier);
		return optimizer_strategy_selection<GCV_Stochastic<CarrierType, 1>, CarrierType>(optim, carrier);
	}
	else if(optr->get_loss_function() == "GCV" && optr->get_DOF_evaluation() == "not_required")
	{
		//ADD
	}
	else if(optr->get_loss_function() == "unused" && optr->get_DOF_evaluation() == "not_required")
	{
		//ADD
	}
}

template<typename EvaluationType, typename CarrierType>
std::pair<VectorXr, output_Data> optimizer_strategy_selection(EvaluationType & optim, CarrierType & carrier)
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
		VectorXr solution = carrier.apply(output.lambda_sol);

		output.time_partial=T.tv_sec + 1e-9*T.tv_nsec;

                //postponed after apply in order to have betas computed
                output.betas=carrier.get_model()->getBeta();

                //std::cout<<"sol"<<solution<<std::endl;

                return {solution, output};
		//Solution_Builders::GCV_batch_sol(solution, output_vec);
	}
	else
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
		VectorXr solution = carrier.apply(lambda_couple.first);
                //std::cout<<"sol"<<solution<<std::endl;
                 //postponed after apply in order to have betas computed
		//now the last values in GCV_exact are the correct ones, related to the last iteration
	         output_Data  output = optim_p->F.get_output(lambda_couple, T, GCV_v_, lambda_v_, ch.which()); //this is why F has to be public in Opt_methods
		 //the copy is necessary for the bulders outside


		return {solution, output};
		//Solution_Builders::GCV_Newton_sol(solution, output);  // [[TO DO make this a template according to methods]]
	}
}
/*
const MatrixXv& solution = regression.getSolution();
const MatrixXr& dof = regression.getDOF();
const MatrixXr & GCV = regression.getGCV();
UInt bestLambda = optimizationData.get_best_lambda_S();
MatrixXv beta;
if(regressionData.getCovariates()->rows()==0)
{
	beta.resize(1,1);
	beta(0,0).resize(1);
	beta(0,0)(0) = 10e20;
}
else
	 beta = regression.getBeta();

const MatrixXr & barycenters = regression.getBarycenters();
const VectorXi & elementIds = regression.getElementIds();

//Copy result in R memory
SEXP result = NILSXP;
result = PROTECT(Rf_allocVector(VECSXP, 5+5+2));
SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution(0).size(), solution.size()));
SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));
SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, solution.size()));
SET_VECTOR_ELT(result, 3, Rf_allocVector(INTSXP, 1));
SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));

Real *rans = REAL(VECTOR_ELT(result, 0));
for(UInt j = 0; j < solution.size(); j++)
{
	for(UInt i = 0; i < solution(0).size(); i++)
		rans[i + solution(0).size()*j] = solution(j)(i);
}

Real *rans1 = REAL(VECTOR_ELT(result, 1));
for(UInt i = 0; i < solution.size(); i++)
{
	rans1[i] = dof(i);
}

//! Copy GCV vector
Real *rans2 = REAL(VECTOR_ELT(result, 2));
for(UInt i = 0; i < solution.size(); i++)
{
	rans2[i] = GCV(i);
}

//! Copy best lambda
UInt *rans3 = INTEGER(VECTOR_ELT(result, 3));
rans3[0] = bestLambda;

//! Copy betas
Real *rans4 = REAL(VECTOR_ELT(result, 4));
for(UInt j = 0; j < beta.size(); j++)
{
	for(UInt i = 0; i < beta(0).size(); i++)
		rans4[i + beta(0).size()*j] = beta(j)(i);
}

//SEND TREE INFORMATION TO R
SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1)); //tree_header information
int *rans5 = INTEGER(VECTOR_ELT(result, 5));
rans5[0] = mesh.getTree().gettreeheader().gettreelev();

SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
Real *rans6 = REAL(VECTOR_ELT(result, 6));
for(UInt i = 0; i < ndim*2; i++)
	rans6[i] = mesh.getTree().gettreeheader().domainorig(i);

SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
Real *rans7 = REAL(VECTOR_ELT(result, 7));
for(UInt i = 0; i < ndim*2; i++)
	rans7[i] = mesh.getTree().gettreeheader().domainscal(i);


UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
SET_VECTOR_ELT(result, 8, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
int *rans8 = INTEGER(VECTOR_ELT(result, 8));
for(UInt i = 0; i < num_tree_nodes; i++)
		rans8[i] = mesh.getTree().gettreenode(i).getid();

for(UInt i = 0; i < num_tree_nodes; i++)
		rans8[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

for(UInt i = 0; i < num_tree_nodes; i++)
		rans8[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

SET_VECTOR_ELT(result, 9, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
Real *rans9 = REAL(VECTOR_ELT(result, 9));
for(UInt j = 0; j < ndim*2; j++)
{
	for(UInt i = 0; i < num_tree_nodes; i++)
		rans9[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
}

//SEND BARYCENTER INFORMATION TO R
SET_VECTOR_ELT(result, 10, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
int *rans10 = INTEGER(VECTOR_ELT(result, 10));
for(UInt i = 0; i < elementIds.rows(); i++)
	rans10[i] = elementIds(i);

SET_VECTOR_ELT(result, 11, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
Real *rans11 = REAL(VECTOR_ELT(result, 11));
for(UInt j = 0; j < barycenters.cols(); j++)
{
	for(UInt i = 0; i < barycenters.rows(); i++)
		rans11[i + barycenters.rows()*j] = barycenters(i,j);
}

UNPROTECT(1);
return(result);
*/

#endif
