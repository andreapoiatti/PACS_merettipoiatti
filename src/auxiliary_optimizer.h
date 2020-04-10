#ifndef __AUXILIARYOPTIMIZER_HPP__
#define __AUXILIARYOPTIMIZER_HPP__

#include "fdaPDE.h"
#include "carrier.h"
#include "solver.h"

template <bool ... b>
struct multi_bool_type
{

};

typedef multi_bool_type<true> t_type;
typedef multi_bool_type<false> f_type;
typedef multi_bool_type<true, true> tt_type;
typedef multi_bool_type<false, true> ft_type;
typedef multi_bool_type<true, false> tf_type;
typedef multi_bool_type<false, false> ff_type;

struct AuxiliaryOptimizer
{
        template<typename InputCarrier>
        static typename std::enable_if< std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_T_setter(SpMat & T, const InputCarrier & carrier);

        template<typename InputCarrier>
        static typename std::enable_if< std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_T_setter(SpMat & T, const InputCarrier & carrier);

        static void set_T_ln_nW_a(SpMat & T, const VectorXr * Ap, const std::vector<UInt> * kp, UInt s);
        static void set_T_ln_W_a(SpMat & T, const VectorXr * Ap, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt s);
        static void set_T_lnn_nW_a(SpMat & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp);
        static void set_T_lnn_W_a(SpMat & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp);
        static void set_T_ln_nW_ptw(SpMat & T, const std::vector<UInt> * kp, UInt s);
        static void set_T_ln_W_ptw(SpMat & T, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt s);
        static void set_T_lnn_nW_ptw(SpMat & T, const SpMat * psip, const SpMat * psi_tp);
        static void set_T_lnn_W_ptw(SpMat & T, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp);

        template<typename InputCarrier>
        static typename std::enable_if< std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_E_setter(SpMat & E, const InputCarrier & carrier);

        template<typename InputCarrier>
        static typename std::enable_if< std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_E_setter(SpMat & E, const InputCarrier & carrier);

        static void set_E_ln_W_ptw(SpMat & E, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt nr, UInt s);
        static void set_E_lnn_W_ptw(SpMat & E, const SpMat * psi_tp, const MatrixXr * Qp);
        static void set_E_ln_W_a(SpMat & E, const std::vector<UInt> * kp, const MatrixXr * Qp, const VectorXr * Ap, UInt nr, UInt s);
        static void set_E_lnn_W_a(SpMat & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap);
        static void set_E_ln_nW_a(SpMat & E, const std::vector<UInt> * kp, const VectorXr * Ap, UInt nr, UInt s);
        static void set_E_lnn_nW_a(SpMat & E, const SpMat * psi_tp, const VectorXr * Ap);

        static void set_z_hat_W(VectorXr & z_hat, const MatrixXr * Hp, const MatrixXr * Qp, const SpMat & S, const VectorXr * zp);
        static void set_z_hat_nW(VectorXr & z_hat, const SpMat & S, const VectorXr * zp);
};


template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_T_setter(SpMat & T, const InputCarrier & carrier)
        {
                const VectorXr * Ap = carrier.get_Ap();
                if (carrier.loc_are_nodes())
                {
                        // Non-null entries are added directly to the structure
                        const std::vector<UInt> * kp  = carrier.get_obs_indicesp();
                        const UInt s = carrier.get_n_obs();

                        if (carrier.has_W())
                        {
                                const MatrixXr * Qp = carrier.get_Qp();
                                AuxiliaryOptimizer::set_T_ln_W_a(T, Ap, kp, Qp, s);
                        }
                        else
                        {
                                AuxiliaryOptimizer::set_T_ln_nW_a(T, Ap, kp, s);
                        }

                }
                else
                {
                        const SpMat * psip = carrier.get_psip();
                        const SpMat * psi_tp = carrier.get_psi_tp();

                        if (carrier.has_W())
                        {
                                const MatrixXr * Qp = carrier.get_Qp();
                                AuxiliaryOptimizer::set_T_lnn_W_a(T, Ap, psip, psi_tp, Qp);
                        }
                        else
                        {
                                AuxiliaryOptimizer::set_T_lnn_nW_a(T, Ap, psip, psi_tp);
                        }
                }

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_T_setter(SpMat & T, const InputCarrier & carrier)
        {
                if (carrier.loc_are_nodes())
                {
                        // Non-null entries are added directly to the structure
                        const std::vector<UInt> * kp  = carrier.get_obs_indicesp();
                        const UInt s = carrier.get_n_obs();

                        if (carrier.has_W())
                        {
                                const MatrixXr * Qp = carrier.get_Qp();
                                AuxiliaryOptimizer::set_T_ln_W_ptw(T, kp, Qp, s);
                        }
                        else
                        {
                                AuxiliaryOptimizer::set_T_ln_nW_ptw(T, kp, s);
                        }

                }
                else
                {
                        const SpMat * psip = carrier.get_psip();
                        const SpMat * psi_tp = carrier.get_psi_tp();

                        if (carrier.has_W())
                        {
                                const MatrixXr * Qp = carrier.get_Qp();
                                AuxiliaryOptimizer::set_T_lnn_W_ptw(T, psip, psi_tp, Qp);
                        }
                        else
                        {
                                AuxiliaryOptimizer::set_T_lnn_nW_ptw(T, psip, psi_tp);
                        }
                }

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_E_setter(SpMat & E, const InputCarrier & carrier)
        {
                const VectorXr * Ap = carrier.get_Ap();
                if (carrier.has_W())
                {
                        const MatrixXr * Qp = carrier.get_Qp();        // Q != I
                        if (carrier.loc_are_nodes())
                        {
                                // Psi is permutation
                                const UInt nr = carrier.get_psip()->cols();
                                const UInt s = carrier.get_n_obs();
                                const std::vector<UInt> * kp  = carrier.get_obs_indicesp();
                                AuxiliaryOptimizer::set_E_ln_W_a(E, kp, Qp, Ap, nr, s);
                        }
                        else
                        {
                                // Psi is full
                                const SpMat * psi_tp = carrier.get_psi_tp();
                                AuxiliaryOptimizer::set_E_lnn_W_a(E, psi_tp, Qp, Ap);
                        }

                }
                else
                {
                        // Q == I
                        if (carrier.loc_are_nodes())
                        {
                                // Psi is permutation
                                const UInt nr = carrier.get_psip()->cols();
                                const UInt s = carrier.get_n_obs();
                                const std::vector<UInt> * kp  = carrier.get_obs_indicesp();
                                AuxiliaryOptimizer::set_E_ln_nW_a(E, kp, Ap, nr, s);
                        }
                        else
                        {
                                // Psi is full
                                const SpMat * psi_tp = carrier.get_psi_tp();
                                AuxiliaryOptimizer::set_E_lnn_nW_a(E, psi_tp, Ap);
                        }

                }
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_E_setter(SpMat & E, const InputCarrier & carrier)
        {
                const MatrixXr * Qp = carrier.get_Qp();        // Q != I
                if (carrier.loc_are_nodes())
                {
                        // Psi is permutation
                        const UInt nr = carrier.get_psip()->cols();
                        const UInt s = carrier.get_n_obs();
                        const std::vector<UInt> * kp  = carrier.get_obs_indicesp();
                        AuxiliaryOptimizer::set_E_ln_W_ptw(E, kp, Qp, nr, s);
                }
                else
                {
                        // Psi is full
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        AuxiliaryOptimizer::set_E_lnn_W_ptw(E, psi_tp, Qp);
                }
        }

#endif
