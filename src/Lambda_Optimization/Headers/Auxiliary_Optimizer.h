#ifndef __AUXILIARY_OPTIMIZER_H__
#define __AUXILIARY_OPTIMIZER_H__

#include <functional>
#include <string>
#include "../../FdaPDE.h"
#include "Carrier.h"
#include "../../FE_Assemblers_Solvers/Headers/Solver.h"
#include "../../Global_Utilities/Headers/Solver_Definitions.h"
#include "Solution_Builders.h"


template<typename InputCarrier, typename Enable = void>
struct AuxiliaryData
{
        MatrixXr K_;                            //!< stores T^{-1}*R                                              [[nnodes x nnodes]]
        MatrixXr F_;                            //!< stores K*v
        VectorXr t_;                            //!< stores dS*z;
        Real     a_;                            //!< Stores the value of <eps_hat, dS*z>
        Real     b_;                            //!< Stores <t, Q*t>
        Real     c_;                            //!< Stores <eps_hat, ddS*z>
};


template<typename InputCarrier>
struct AuxiliaryData<InputCarrier, typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value>::type>
{
        public:
                MatrixXr K_;                            //!< stores T^{-1}*R                                              [[nnodes x nnodes]]
                MatrixXr F_;                            //!< stores K*v
                VectorXr t_;                            //!< stores dS*z;
                Real     a_;                            //!< Stores the value of <eps_hat, dS*z>
                Real     b_;                            //!< Stores <t, Q*t>
                Real     c_;                            //!< Stores <eps_hat, ddS*z>
                VectorXr f_;
                VectorXr g_;
                VectorXr h_;
                VectorXr p_;
                VectorXr r_;

        void left_multiply_by_psi(const InputCarrier & carrier, VectorXr & ret, const VectorXr & vec);
};

struct AuxiliaryOptimizer
{
        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_T_setter(MatrixXr & T, const InputCarrier & carrier);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_T_setter(MatrixXr & T, const InputCarrier & carrier);

        static void set_T_nW_a(MatrixXr & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp);
        static void set_T_W_a(MatrixXr & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp);
        static void set_T_ln_nW_ptw(MatrixXr & T, const std::vector<UInt> * kp, UInt s);
        static void set_T_ln_W_ptw(MatrixXr & T, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt s);
        static void set_T_lnn_nW_ptw(MatrixXr & T, const SpMat * psip, const SpMat * psi_tp);
        static void set_T_lnn_W_ptw(MatrixXr & T, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp);
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
                universal_z_hat_setter(VectorXr & z_hat, const InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_z_hat_setter(VectorXr & z_hat, const InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda);

        template<typename InputCarrier>
        static void common_z_hat_part(VectorXr & z_hat, const InputCarrier & carrier, const MatrixXr & S);

        static void set_z_hat_W(VectorXr & z_hat, const MatrixXr * Hp, const MatrixXr * Qp, const MatrixXr & S, const VectorXr * zp);
        static void set_z_hat_nW(VectorXr & z_hat, const MatrixXr & S, const VectorXr * zp);
        /* -------------------------------------------------------------------*/

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, const UInt nnodes);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, const UInt nnodes);
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
                universal_second_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda);

        template<typename InputCarrier>
        static typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_second_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda);
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
