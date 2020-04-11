#ifndef VECTOREVAL_H
#define VECTOREVAL_H

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "fdaPDE.h"
#include "solver.h"
#include "function_variadic.h"



template <typename Tuple, typename Hessian, typename... Extensions>
class Vec_evaluation
{
        protected:
                std::vector<Real> lambda_vec;

                Vec_evaluation(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_, const std::vector<Real> &lambda_vec_): F(F_), lambda_vec(lambda_vec_) {Rprintf("Vector evaluator built\n");};

                virtual void compute_specific_parameters(void) {;}; //Computes particular parameters related to the mimizing solution. It does nothing if not implemented



        public:
                Function_Wrapper<Tuple, Tuple, Tuple, Hessian, Extensions...> & F; //needed to be public, to be able to access to other methods of the class F from outside


                std::pair<std::vector<Real>, UInt> compute_vector(void) //returns the vector of evaluations and the index of the corresponding minimum
                {
                  UInt dim=lambda_vec.size();

                  std::vector<Real> evaluations(dim);


                  UInt index_min=0; //Assume the first one is the minimum


                  for (UInt i=0;i<dim;i++)
                  {
                        evaluations[i]=this->F.evaluate_f(this->lambda_vec[i]); //only scalar functions;

                        if (evaluations[i]<evaluations[index_min])
                                {
                                    index_min=i;
                                    this->compute_specific_parameters();
                                }
                   }

//DEBUGGING
                for (UInt i=0;i<dim;i++)
                   {
                     std::cout<<"Lambda: "<<this->lambda_vec[i]<<" GCV: "<<evaluations[i]<<std::endl;
                   }

                std::cout<<"Lambda opt: "<<this->lambda_vec[index_min]<<" and GCV: "<<evaluations[index_min]<<std::endl;

                return {evaluations,index_min};

                };
};


struct output_Data_opt
{
        VectorXr        z_hat_opt;                          //!< Model predicted values in the locations
        Real            SS_res_opt;                         //!< Model predicted sum of squares of the residuals
        Real            sigma_hat_sq_opt;                   //!< Model estimated variance of errors
        Real            lambda_opt;                         //!< Lambda obtained in the solution
        std::vector<Real> GCV_evals;                        //!< GCV evaluations vector
};

// i metodi per le funzioni in F si chiamano evaluate_f, evaluate_first_derivative, evaluate_second_derivative
template <typename Tuple, typename Hessian, typename ...Extensions>
class Eval_GCV: public Vec_evaluation<Tuple, Hessian, Extensions...>
{
        protected:

                output_Data_opt output;

                void compute_specific_parameters(void) override
                {
                 this->output.z_hat_opt=this->F.get_output_partial().z_hat;
                 this->output.SS_res_opt=this->F.get_output_partial().SS_res;
                 this->output.sigma_hat_sq_opt=this->F.get_output_partial().sigma_hat_sq;
                 };

        public:
                Eval_GCV(Function_Wrapper<Tuple, Real, Tuple, Real, Extensions...> & F_, const std::vector<Real> &lambda_vec_): Vec_evaluation<Tuple, Hessian, Extensions...>(F_,lambda_vec_) {}; //esempio di possibile constructor
                // non può prendere in ingresso const ref, deve modificare l'oggetto F

                 const output_Data_opt& Get_optimization_vectorial(void)
                 {
                  std::pair<std::vector<Real>, UInt> p=this->compute_vector();
                  this->output.GCV_evals=p.first;
                  this->output.lambda_opt=this->lambda_vec.at(p.second); //safer to access, if there have been errors before
                  return this->output;
                  };


};


#endif
