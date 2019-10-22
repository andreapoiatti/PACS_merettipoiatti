#ifndef GCV_H
#define GCV_H

// Headers
#include "fdaPDE.h"
#include "mixedFERegression.h"
#include "solver.h"


// Classes
// **** General method ****
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim, UInt size>
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
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>
{
        protected:
                std::tuple<Real, Real, Real> last_lambda;              //!< tuple of previousy used lambdas for respectively, f, fp, fs
                //! Model containing all the information necessary for the computation of the optimal value
                MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model;

                // Constructors
                //! Constructor of the class given the model
                /*! \param model the structure from which to take all the data for the derived classes
                 */
                Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_):
                        model(model_), last_lambda(std::make_tuple(-1, -1, -1)) {}

                Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_, Real lambda0):
                        model(model_), last_lambda(std::make_tuple(lambda0, -1, -1)) {}

        virtual void update_parameters(Real lambda) = 0;
};

//----------------------------------------------------------------------------//
// *** GCV-based***

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim, UInt size>
class GCV_Family: public Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, size>
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
*/
};


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>: Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>
{
        protected:
                using  Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>::last_lambda;
                using  Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>::model;

                //Useful common data
                VectorXr        z_hat;                          //!< Model predicted values in the locations
                VectorXr        eps_hat;                        //!< Model predicted error in the locations (residuals)
                Real            SS_res;                         //!< Model predicted sum of squares of the residuals
                UInt            s;                              //!< Model number of observations
                Real            sigma_hat_sq;                   //!< Model estimated variance of errors
                Real            aux;                            //!< Stores the value of <eps_hat, dS*z>

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

                // Degrees of freedom
                Real dof;                                       //!< tr(S) + q, degrees of freedom of the model
                Real dor;                                       //!< s - dof, degrees of freedom of the residuals

                // Setters of the matrices
                void set_R_(void);
                void set_T_(Real lambda);
                void set_V_(void);
                void set_S_and_trS_(void);
                void set_dS_and_trdS_(void);
                void set_ddS_and_trddS_(void);

                // Setters of the common data
                void compute_s(void);
                void compute_z_hat(void);
                void compute_eps_hat(void);
                void compute_SS_res(void);
                void compute_sigma_hat_sq(void);
                void compute_aux(void);

                // Updaters
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
        virtual void update_dof(void)           = 0;
        virtual void update_dor(Real lambda)    = 0;

                // Constructors
                //! Constructor of the class in absenceof an initial value
                GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_):
                        Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>(model_)
                        {
                                compute_s();
                                set_R_();
                        }

                //! Initial guess about lambda
                GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_, Real lambda0):
                        Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>(model_, lambda0)
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
};

//----------------------------------------------------------------------------//
// ** GCV_Exact **

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim, UInt size>
class GCV_Exact: public GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, size>
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
        private:
                void update_dof (void) override;

        public:
                GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, size>
                        (const MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim, size> & model_):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, size>(model) {}

                GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, size>
                        (const MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim, size> & model_, const VectorXr & lambda0):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, size>(model, lambda0)
                        {
                                update_dof();
                                update_dor();
                        }

                void update_parameters(const VectorXr & lambda) override;
*/
};


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, 1>: public GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>
{
        private:
                void update_dof(void)           override;
                void update_dor(Real lambda)    override;

        public:
                GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>(model_) {}

                GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_, Real lambda0):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>(model_, lambda0)
                        {
                                update_dof();
                                update_dor();
                        }

                void update_parameters(Real lambda) override;
};

//----------------------------------------------------------------------------//
// ** GCV_Stochastic **

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim, UInt size>
class GCV_Stochastic: public GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, size>
{
/*
        [[VERSION WITH TIMES STILL TO BE IMPLEMENTED]]
        private:
                void update_dof (void) override;

        public:
                GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, size>
                        (const MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim, size> & model_):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, size>(model) {}

                GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, size>
                        (const MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim, size> & model_, const VectorXr & lambda0):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, size>(model, lambda0)
                        {
                                update_dof();
                                update_dor();
                        }

                void update_parameters(const VectorXr & lambda) override;
*/
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, 1>: public GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>
{
        private:
                void update_dof(void)           override;
                void update_dor(Real lambda)    override;

        public:
                GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>(model_){}

                GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, 1>
                        (MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim> & model_, Real lambda0):
                        GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>(model_, lambda0)
                        {
                                update_dof();
                                update_dor();
                        }

                void update_parameters(Real lambda) override;
};

//----------------------------------------------------------------------------//
// *** K-FOLD GCV ***

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim, UInt size>
class K_fold_GCV: public Lambda_optimizer<InputHandler, Integrator, ORDER, mydim, ndim, 1>
{
};

// [[TODO]]
// look at
// VISITOR PATTERN ALEXANDRESCU
// Inheritance + overloading -> "virtual"

#include "lambda_optimizer_imp.h"

#endif
