#ifndef GCV_H
#define GCV_H

// Headers
#include "fdaPDE.h"
#include "gof_updater.h"
#include "mixedFERegression.h"
#include "solver.h"
#include "auxiliary_optimizer.h"
#include <algorithm>

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
                //! Model containing all the information necessary for the computation of the optimal value
                const InputCarrier & the_carrier;

                // Constructors
                //! Constructor of the class given the model
                /*! \param model the structure from which to take all the data for the derived classes
                 */
                Lambda_optimizer<InputCarrier, 1>(InputCarrier & the_carrier_):
                        the_carrier(the_carrier_) {}

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
                using  Lambda_optimizer<InputCarrier, 1>::the_carrier;

                // Common data
                VectorXr        z_hat;                          //!< Model predicted values in the locations
                VectorXr        eps_hat;                        //!< Model predicted error in the locations (residuals)
                Real            SS_res;                         //!< Model predicted sum of squares of the residuals
                Real            sigma_hat_sq;                   //!< Model estimated variance of errors
                UInt            s;                              //!< Model number of observations (i.e. locations)
                output_Data     output;            //Output, needed to be used in FdaPDE.h, necessarily public

                // Degrees of freedom
                Real dof;                                       //!< tr(S) + q, degrees of freedom of the model
                Real dor;                                       //!< s - dof, degrees of freedom of the residuals

                // Setters of the common data
        virtual void compute_z_hat(Real lambda) = 0;
                void compute_eps_hat(void);
                void compute_SS_res(void);
                void compute_sigma_hat_sq(void);
                void compute_s(void);

                // Updaters
                void update_errors(Real lambda);

                // DOF methods
        virtual void update_dof(Real lambda)    = 0;
        virtual void update_dor(Real lambda)    = 0;

                // Constructors
                //! Constructor of the class in absenceof an initial value
                GCV_Family<InputCarrier, 1>(InputCarrier & the_carrier_):
                        Lambda_optimizer<InputCarrier, 1>(the_carrier_)
                        {
                                this->compute_s();
                        }

        public:
                void zero_updater(Real lambda);

        virtual Real compute_f( Real lambda)        = 0;

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
                GOF_updater<GCV_Exact<InputCarrier, 1>, Real> gu;

                MatrixXr  R_; 		//!< stores the value of R1^t*R0^{-1}*R1                          [[nnodes x nnodes]]
                MatrixXr  T_; 		//!< stores the value of Psi^t*Q*Psi+lambda*R                     [[nnodes x nnodes]]
                MatrixXr  V_; 		//!< stores the value of T^{-1}*Psi^t*Q                           [[nnodes x   s   ]]
                MatrixXr  S_;           //!< stores the value of Psi*V [as in Stu-Hunter Sangalli]        [[   s   x   s   ]]
                Real      trS_;         //!< stores the value of the trace of S
                MatrixXr  dS_;          //!< stores the derivative of S w.r.t. lambda                     [[   s   x   s   ]]
                Real      trdS_;        //!< stores the value of the trace of dS
                MatrixXr  ddS_;         //!< stores the second derivative of S w.r.t. lambda              [[   s   x   s   ]]
                Real      trddS_;       //!< stores the value of the trace of ddS

                // Utility matrices [just the ones for the specific carrier]
                AuxiliaryData<InputCarrier> adt;

                void compute_z_hat (Real lambda) override;
                void update_dof(Real lambda)     override;
                void update_dor(Real lambda)     override;

                // Setters of the matrices
                void set_R_(void);
                void set_T_(Real lambda);
                void set_V_(void);
                void set_S_and_trS_(void);
                void set_dS_and_trdS_(void);
                void set_ddS_and_trddS_(void);

                // Utilities
                void LeftMultiplybyPsiAndTrace(Real & trace, MatrixXr & ret, const MatrixXr & mat);

                // Partial updaters
                void update_matrices(Real lambda);

        public:
                GCV_Exact<InputCarrier, 1>(InputCarrier & the_carrier_):
                        GCV_Family<InputCarrier, 1>(the_carrier_)
                        {
                                this->set_R_();
                        }

                void first_updater(Real lambda);
                void second_updater(Real lambda);

                Real compute_f( Real lambda) override;
                Real compute_fp(Real lambda);
                Real compute_fs(Real lambda);

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


template<typename InputCarrier>
class GCV_Stochastic<InputCarrier, 1>: public GCV_Family<InputCarrier, 1>
{
        private:
                GOF_updater<GCV_Stochastic<InputCarrier, 1>, Real> gu;

                MatrixXr US_;
                bool     us = false;

                void compute_z_hat (Real lambda) override;
                void update_dof(Real lambda)     override;
                void update_dor(Real lambda)     override;

                void set_US_(void);

        public:
                GCV_Stochastic<InputCarrier, 1>(InputCarrier & the_carrier_):
                        GCV_Family<InputCarrier, 1>(the_carrier_)
                        {
                                this->set_US_();
                        }

                void first_updater(Real lambda)  {; /*Dummy*/}
                void second_updater(Real lambda) {; /*Dummy*/}

                Real compute_f( Real lambda) override;
                Real compute_fp(Real lambda) {return 0; /*Dummy*/}
                Real compute_fs(Real lambda) {return 0; /*Dummy*/}

                void update_parameters(Real lambda) override;
};


//----------------------------------------------------------------------------//
// *** K-FOLD GCV ***

template<typename InputCarrier, UInt size>
class K_fold_GCV: public Lambda_optimizer<InputCarrier, size>
{
        // [[ TO BE IMPLEMENTED ]]
};

// Visitor pattern Alexandrescu

#include "lambda_optimizer_imp.h"

#endif
