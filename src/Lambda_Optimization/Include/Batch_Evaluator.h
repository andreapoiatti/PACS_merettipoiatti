#ifndef __BATCH_EVALUATOR_H__
#define __BATCH_EVALUATOR_H__

// Include
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "Function_Variadic.h"
#include "Solution_Builders.h"


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

                Vec_evaluation(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_, const std::vector<Real> & lambda_vec_):
                        F(F_), lambda_vec(lambda_vec_) {//Debugging purpose//Rprintf("Vector evaluator built\n");
                        }; //!< Constructor

                /*!
                Function to compute particular parameters related to the mimizing solution. It does nothing if not implemented. It is not pure virtual in order to be general and leave the possibility of instantiating the object without impementing that function
                */
                virtual void compute_specific_parameters(void) {}; //!<Computes particular parameters related to the mimizing solution, for all the values evaluated. It does nothing if not implemented
                virtual void compute_specific_parameters_best(void) {}; //!<Computes particular parameters related to the mimizing solution, only for minimizing solutions. It does nothing if not implemented

        public:
                Function_Wrapper<Tuple, Tuple, Tuple, Hessian, Extensions...> & F; //!< F needed to be public, to be able to access to other methods of the class F from outside*/

                std::pair<std::vector<Real>, UInt> compute_vector(void) //!< Function which returns the vector of evaluations of GCV and the index of the corresponding minimum
                {
                        UInt dim = lambda_vec.size();
                        UInt index_min = 0; //Assume the first one is the minimum
                        std::vector<Real> evaluations(dim);

                        for (UInt i=0; i<dim; i++)
                        {
                                this->F.set_index(i);
                                evaluations[i] = this->F.evaluate_f(this->lambda_vec[i]); //only scalar functions;

                                this->compute_specific_parameters();
                                if (i==0)
                                     this->compute_specific_parameters_best();

                                if (evaluations[i]<evaluations[index_min])
                                {
                                        this->compute_specific_parameters_best();
                                        index_min=i;
                                }
                        }

                        //DEBUGGING PURPOSE
                        /*for (UInt i=0; i<dim; i++)
                        {
                                Rprintf("\nLambda: %f, GCV: %f\n",this->lambda_vec[i], evaluations[i]);
                        }

                        Rprintf("\nLambda opt: %f and GCV: %f\n",this->lambda_vec[index_min], evaluations[index_min]);
                         */
                        return {evaluations,index_min};
                }
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

                void compute_specific_parameters(void) override
                {
                 this->F.set_output_partial();
                }

                void compute_specific_parameters_best(void) override //!< Computes specific parameters needed for GCV
                {
                        //Debugging purpose
                        // Rprintf("Specific parameters for GCV computed\n");

                        this->F.set_output_partial_best();
                 }

        public:
                Eval_GCV(Function_Wrapper<Tuple, Real, Tuple, Real, Extensions...> & F_, const std::vector<Real> & lambda_vec_):
                        Vec_evaluation<Tuple, Hessian, Extensions...>(F_,lambda_vec_) {}; //! Constructor

                output_Data  Get_optimization_vectorial(void) //! Output constructor
                {
                        std::pair<std::vector<Real>, UInt> p = this->compute_vector();
                        output_Data output=this->F.get_output_full();
                        output.GCV_evals  = p.first;
                        output.lambda_sol = this->lambda_vec.at(p.second); //Safer use of at instead of []
                        output.lambda_pos = 1+p.second; //in R numbering
                        output.lambda_vec = this->lambda_vec;
                        output.GCV_opt    = p.first.at(p.second);

                        return output;
                }
};

#endif
