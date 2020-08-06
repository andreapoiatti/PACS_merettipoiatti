#ifndef __BATCH_EVALUATOR_H__
#define __BATCH_EVALUATOR_H__

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Headers/Solver.h"
#include "Function_Variadic.h"


//! Father class for a scalar function evaluation of a given vector of lambda values computing the minimum fuction value
/*!
 * \tparam       Tuple          image type of the gradient of the function
 * \tparam       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam       Extensions    input class if the computations need members already stored in a class
 */
template <typename Tuple, typename Hessian, typename... Extensions>
class Vec_evaluation
{
        protected:
                std::vector<Real> lambda_vec;    //!< Vector of lambda to be evaluated

                Vec_evaluation(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_, const std::vector<Real> &lambda_vec_): F(F_), lambda_vec(lambda_vec_) {Rprintf("Vector evaluator built\n");}; //< Constructor

                /*!
                Function to compute particular parameters related to the mimizing solution. It does nothing if not implemented. It is not pure virtual in order to be general and leave the possibility of instantiating the object without impementing that function
                */
                virtual void compute_specific_parameters(void) {;}; //!<Computes particular parameters related to the mimizing solution. It does nothing if not implemented




        public:
                Function_Wrapper<Tuple, Tuple, Tuple, Hessian, Extensions...> & F; /*!F needed to be public, to be able to access to other methods of the class F from outside*/


                std::pair<std::vector<Real>, UInt> compute_vector(void) //!< Function which returns the vector of evaluations of GCV and the index of the corresponding minimum
                {
                  UInt dim=lambda_vec.size();

                  std::vector<Real> evaluations(dim);


                  UInt index_min=0; //Assume the first one is the minimum


                  for (UInt i=0;i<dim;i++)
                  {
                        evaluations[i]=this->F.evaluate_f(this->lambda_vec[i]); //only scalar functions;
                        if (i==0)
                           this->compute_specific_parameters(); //for the first evaluation it is needed to be computed
                        if (evaluations[i]<evaluations[index_min])
                                {
                                    index_min=i;
                                    this->compute_specific_parameters();
                                }
                   }


//DEBUGGING PURPOSE
                for (UInt i=0;i<dim;i++)
                   {
                     Rprintf("\nLambda: %f, GCV: %f\n",this->lambda_vec[i],evaluations[i]);
                   }

                Rprintf("\nLambda opt: %f and GCV: %f\n",this->lambda_vec[index_min], evaluations[index_min]);

                return {evaluations,index_min};

                };
};


struct output_Data_opt
{
        VectorXr        z_hat_opt;                          //!< Model predicted values in the locations
        Real            SS_res_opt;                         //!< Model predicted sum of squares of the residuals
        Real            sigma_hat_sq_opt;                   //!< Model estimated variance of errors
        Real            lambda_opt;                        //!< Lambda obtained in the solution
        Real            GCV_opt;                           //!<GCV optimal comptued in the vector of lambdas
        std::vector<Real> GCV_evals;                       //!< GCV evaluations vector
};


/*!
 * \tparam       Tuple          image type of the gradient of the function
 * \tparam       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam       Extensions    input class if the computations need members already stored in a class
 */
//!<Class inheriting form class Vec_evaluation, the function used is GCV evaluation
template <typename Tuple, typename Hessian, typename ...Extensions>
class Eval_GCV: public Vec_evaluation<Tuple, Hessian, Extensions...>
{
        protected:

                output_Data_opt output;

                void compute_specific_parameters(void) override //!< Computes specific parameters needed for GCV
                {

                 Rprintf("Specific parameters for GCV computed\n");

                 this->output.z_hat_opt=this->F.get_output_partial().z_hat;
                 this->output.SS_res_opt=this->F.get_output_partial().SS_res;
                 this->output.sigma_hat_sq_opt=this->F.get_output_partial().sigma_hat_sq;

                 };

        public:
                Eval_GCV(Function_Wrapper<Tuple, Real, Tuple, Real, Extensions...> & F_, const std::vector<Real> &lambda_vec_): Vec_evaluation<Tuple, Hessian, Extensions...>(F_,lambda_vec_) {}; //! Constructor

                 const output_Data_opt& Get_optimization_vectorial(void) //! Output constructor
                 {

                  std::pair<std::vector<Real>, UInt> p=this->compute_vector();
                  this->output.GCV_evals=p.first;
                  this->output.lambda_opt=this->lambda_vec.at(p.second);
                  this->output.GCV_opt=p.first.at(p.second);

                  return this->output;

                  };


};


#endif
