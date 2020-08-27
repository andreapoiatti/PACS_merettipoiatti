#ifndef __AUXILIARY_OPTIMIZER_H__
#define __AUXILIARY_OPTIMIZER_H__

// HEADERS
#include <functional>
#include <string>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"
#include "Carrier.h"
#include "Solution_Builders.h"

// CLASSES
//! Template class for data storing and efficient management.
/*!
 General purpose class storing data useful for fastening computation in
 Lambda_optimizer derived classes and AuxiliaryOptimizer. Its content
 are matrices, vectors and doubles useful for GCV calculations.
 \tparam InputCarrier the type of Carrier used in the optimization.
 \tparam Enaable dummy typename for SFINAE instantiation of a more refined version for problems with forcing terms.
 \sa Carrier, AuxiliaryOptimizer, Lambda_optimizer
*/
template<typename InputCarrier, typename Enable = void>
struct AuxiliaryData
{
        MatrixXr K_;                            //!< stores T^{-1}*R                            [nnodes x nnodes]
        MatrixXr F_;                            //!< stores K*v                                 [nnodes x nnodes]
        VectorXr t_;                            //!< stores dS*z;
        Real     a_;                            //!< Stores the value of <eps_hat, dS*z>
        Real     b_;                            //!< Stores <t, Q*t>
        Real     c_;                            //!< Stores <eps_hat, ddS*z>
};

//! Template class for data storing and efficient management in forcing term based problems
/*!
 General purpose class storing data useful for fastening computation in
 Lambda_optimizer derived classes and AuxiliaryOptimizer, spercialized under forcing
 term based problems. Its content are matrices, vectors and doubles
 useful for GCV calculations.
 \tparam InputCarrier the type of Carrier used in the optimization.
 \tparam typename for SFINAE instantiation of this more refined version for problems with forcing terms, distinguishing from base version.
 \sa AuxiliaryData, Carrier, AuxiliaryOptimizer, Lambda_optimizer
*/
template<typename InputCarrier>
struct AuxiliaryData<InputCarrier, typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value>::type>
{
        MatrixXr K_;                            //!< stores T^{-1}*R                            [nnodes x nnodes]
        MatrixXr F_;                            //!< stores K*v                                 [nnodes x nnodes]
        VectorXr t_;                            //!< stores dS*z;
        Real     a_;                            //!< Stores the value of <eps_hat, dS*z>
        Real     b_;                            //!< Stores <t, Q*t>
        Real     c_;                            //!< Stores <eps_hat, ddS*z>
        VectorXr f_;                            //!< Stores R1^T*R0^{-1}*u
        VectorXr g_;                            //!< Stores T^{-1}*f
        VectorXr h_;                            //!< Stores (lambda*K-I)*g
        VectorXr p_;                            //!< Stores Psi*h-t
        VectorXr r_;                            //!< Stores Q*s

        void left_multiply_by_psi(const InputCarrier & carrier, VectorXr & ret, const VectorXr & vec);
};


//! General purpose class to support efficient case-driven computation of Lambda_Optimizer
/*!
 This struct is a collection of static methods called "universal" in their name
 whose main purpose is to use SFINAE on the template InputCarrier type to provide
 correct implementation for each possible method in Lambda_Optimizer derived classes.
 Since functons are static no object of this class needs to be ever created
 \sa AuxiliaryData,  Lambda_optimizer
*/
struct AuxiliaryOptimizer
{
        static void bc_utility(MatrixXr & mat, const std::vector<UInt> * bc_idxp);
        static void bc_utility(SpMat & mat, const std::vector<UInt> * bc_idxp);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_T_setter(MatrixXr & T, InputCarrier & carrier);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_T_setter(MatrixXr & T, InputCarrier & carrier);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_E_setter(MatrixXr & E, const InputCarrier & carrier);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_E_setter(MatrixXr & E, const InputCarrier & carrier);

        static void set_E_ln_W_ptw(MatrixXr & E, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt nr, UInt s);
        static void set_E_lnn_W_ptw(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp);
        static void set_E_W_a(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap);
        static void set_E_nW_a(MatrixXr & E, const SpMat * psi_tp, const VectorXr * Ap);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_z_hat_setter(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_z_hat_setter(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda);

        template<typename InputCarrier>
        static void common_z_hat_part(VectorXr & z_hat, InputCarrier & carrier, const MatrixXr & S);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, InputCarrier & carrier, const MatrixXr & US, const UInt nnodes);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, InputCarrier & carrier, const MatrixXr & US, const UInt nnodes);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_second_updater(AuxiliaryData<InputCarrier> & adt, InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_second_updater(AuxiliaryData<InputCarrier> & adt, InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
                universal_GCV(const Real s, const Real sigma_hat_sq, const Real dor);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
                universal_GCV_d(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS);

        /* -------------------------------------------------------------------*/
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
                universal_GCV_dd(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS, const Real trddS);
};

#include "Auxiliary_Optimizer_imp.h"

#endif
