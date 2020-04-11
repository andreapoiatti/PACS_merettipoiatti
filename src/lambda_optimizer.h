#ifndef GCV_H
#define GCV_H

// Headers
#include "fdaPDE.h"
#include "mixedFERegression.h"
#include "solver.h"
#include "auxiliary_optimizer.h"
#include <algorithm>


//Output struct to be used to return values in R
struct output_Data
{
        VectorXr        z_hat;                          //!< Model predicted values in the locations
        Real            SS_res;                         //!< Model predicted sum of squares of the residuals
        Real            sigma_hat_sq;                   //!< Model estimated variance of errors
        Real            lambda_sol;                     //!<Lambda obratained in the solution
        UInt            n_it;                           //!< Number of iterations for the method
        Real            time_partial;                   //!<Time, from beginning to end of the optimization method
};


// Classes
// **** General method ***
template <typename InputCarrier, UInt size>
class Lambda_optimizer
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
*/
};

/*! /brief{father class for all optimization methods}
 * This virtual class stores the model to be used by all its children,
 * i. e. the classes actually instantiating a real optimization method and
 * performing evaluations. The class also contains a tuple to keep track
 * of the updates used by the derived classes.
 */
template <typename InputCarrier>
class Lambda_optimizer<InputCarrier, 1>
{
        protected:
                std::tuple<Real, Real, Real> last_lambda;              //!< tuple of previousy used lambdas for respectively, f, fp, fs
                //! Model containing all the information necessary for the computation of the optimal value
                const InputCarrier & the_carrier;

                // Constructors
                //! Constructor of the class given the model
                /*! \param model the structure from which to take all the data for the derived classes
                 */
                Lambda_optimizer<InputCarrier, 1>(InputCarrier & the_carrier_):
                        the_carrier(the_carrier_), last_lambda(std::make_tuple(-1, -1, -1)) {}

                Lambda_optimizer<InputCarrier, 1>(InputCarrier & the_carrier_, Real lambda0):
                        the_carrier(the_carrier_), last_lambda(std::make_tuple(lambda0, -1, -1)) {}

        virtual void update_parameters(Real lambda) = 0;
};

//----------------------------------------------------------------------------//
// *** GCV-based***

template <typename InputCarrier, UInt size>
class GCV_Family: public Lambda_optimizer<InputCarrier, size>
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
*/
};


template <typename InputCarrier>
class GCV_Family<InputCarrier, 1>: Lambda_optimizer<InputCarrier, 1>
{
        protected:
                using  Lambda_optimizer<InputCarrier, 1>::last_lambda;
                using  Lambda_optimizer<InputCarrier, 1>::the_carrier;

                //Useful common data
                VectorXr        z_hat;                          //!< Model predicted values in the locations
                VectorXr        eps_hat;                        //!< Model predicted error in the locations (residuals)
                Real            SS_res;                         //!< Model predicted sum of squares of the residuals
                UInt            s;                              //!< Model number of observations
                Real            sigma_hat_sq;                   //!< Model estimated variance of errors
                Real            aux;                            //!< Stores the value of <eps_hat, dS*z>
                output_Data     output;            //Output, needed to be used in FdaPDE.h, necessarily public

                // Utility matrices
                SpMat   	R_; 			        //!< stores the value of R1^t*R0^{-1}*R1                          [[nnodes x nnodes]]
                SpMat   	T_; 				//!< stores the value of Psi^t*Q*Psi+lambda*R                     [[nnodes x nnodes]]
                SpMat   	V_; 			        //!< stores the value of T^{-1}*Psi^t*Q                           [[nnodes x   s   ]]
                SpMat           S_;                             //!< stores the value of Psi*V [as in Stu-Hunter Sangalli]        [[   s   x   s   ]]
                Real            trS_;                           //!< stores the value of the trace of S
                SpMat           K_;                             //!< stores T^{-1}*R                                              [[nnodes x nnodes]]
                SpMat           dS_;                            //!< stores the derivative of S w.r.t. lambda                     [[   s   x   s   ]]
                Real            trdS_;                          //!< stores the value of the trace of dS
                SpMat           ddS_;                           //!< stores the second derivative of S w.r.t. lambda              [[   s   x   s   ]]
                Real            trddS_;                         //!< stores the value of the trace of ddS
                MatrixXr        US_;

                // Degrees of freedom
                Real dof;                                       //!< tr(S) + q, degrees of freedom of the model
                Real dor;                                       //!< s - dof, degrees of freedom of the residuals
                bool us = false;

                // Setters of the matrices
                void set_R_(void);
                void set_T_(Real lambda);
                void set_V_(void);
                void set_S_and_trS_(void);
                void set_dS_and_trdS_(void);
                void set_ddS_and_trddS_(void);
                void set_US_(void);


                // Setters of the common data
                void compute_s(void);
                void compute_z_hat(void);
                void compute_eps_hat(void);
                void compute_SS_res(void);
                void compute_sigma_hat_sq(void);
                void compute_aux(void);

                // Updaters
                void update_family_p1(Real lambda);  //Part 1
                void update_family_p2(void);     //Part 2, common to all
                void update_family(Real lambda);

                void zero_updater(Real lambda);
                void first_updater(Real lambda);
                void second_updater(Real lambda);

                void f_updater(Real lambda);
                void fp_updater(Real lambda);
                void fs_updater(Real lambda);

                // Utilities
                void LeftMultiplybyPsiAndTrace(Real & trace, SpMat & ret, const SpMat & mat);

                // DOF methods
        virtual void update_dof(Real lambda)    = 0;
        virtual void update_dor(Real lambda)    = 0;

                // Constructors
                //! Constructor of the class in absenceof an initial value
                GCV_Family<InputCarrier, 1>(InputCarrier & the_carrier_):
                        Lambda_optimizer<InputCarrier, 1>(the_carrier_)
                        {
                                compute_s();
                                set_R_();
                        }

                //! Initial guess about lambda
                GCV_Family<InputCarrier, 1>(InputCarrier & the_carrier_, Real lambda0):
                        Lambda_optimizer<InputCarrier, 1>(the_carrier_, lambda0)
                        {
                                compute_s();
                                set_R_();
                                update_family(lambda0);
                        }

        public:
                Real compute_f( Real lambda);
                Real compute_fp(Real lambda);
                Real compute_fs(Real lambda);
        virtual void update_parameters(Real lambda) = 0;
                 //Set and return output data, plus lambda final and number of iterations

                const output_Data & get_output(std::pair<Real, UInt> p, const timespec & T);
                const output_Data & get_output_partial(void);

};

//----------------------------------------------------------------------------//
// ** GCV_Exact **

template<typename InputCarrier, UInt size>
class GCV_Exact: public GCV_Family<InputCarrier, size>
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
*/
};


template<typename InputCarrier>
class GCV_Exact<InputCarrier, 1>: public GCV_Family<InputCarrier, 1>
{
        private:
                void update_dof(Real lambda)    override;
                void update_dor(Real lambda)    override;

        public:
                GCV_Exact<InputCarrier, 1>(InputCarrier & the_carrier_): GCV_Family<InputCarrier, 1>(the_carrier_) {}

                GCV_Exact<InputCarrier, 1>(InputCarrier & the_carrier_, Real lambda0):
                        GCV_Family<InputCarrier, 1>(the_carrier_, lambda0)
                        {
                                update_dof(lambda0);
                                update_dor(lambda0);
                        }

                void update_parameters(Real lambda) override;
};

//----------------------------------------------------------------------------//
// ** GCV_Stochastic **

template<typename InputCarrier, UInt size>
class GCV_Stochastic: public GCV_Family<InputCarrier, size>
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
*/
};
/*
template<typename InputCarrier, UInt ndim>
class GCV_Stochastic<InputCarrier, 1>: public GCV_Family<InputCarrier, 1>
{
        private:
                void compute_z_hat (Real lambda);
                void update_dof(Real lambda)    override;
                void update_dor(Real lambda)    override;

        public:
                GCV_Stochastic<InputCarrier, 1>(InputCarrier & the_carrier_): GCV_Family<InputCarrier, 1>(the_carrier_){}

                GCV_Stochastic<InputCarrier, 1>(InputCarrier & the_carrier_, Real lambda0):
                        GCV_Family<InputCarrier, 1>(the_carrier_, lambda0)
                        {
                                update_dof(lambda0);
                                update_dor(lambda0);
                        }

                void update_parameters(Real lambda) override;
};
*/

//----------------------------------------------------------------------------//
// *** K-FOLD GCV ***

template<typename InputCarrier, UInt size>
class K_fold_GCV: public Lambda_optimizer<InputCarrier, size>
{
        // [[ TO BE IMPLEMENTED ]]
};

// [[TODO]]
// look at
// VISITOR PATTERN ALEXANDRESCU
// Inheritance + overloading -> "virtual"

#include "lambda_optimizer_imp.h"

#endif
