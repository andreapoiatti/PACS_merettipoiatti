#ifndef __AUXILIARYOPTIMIZER_HPP__
#define __AUXILIARYOPTIMIZER_HPP__

#include "fdaPDE.h"
#include "carrier.h"
#include "solver.h"

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

template <typename LambdaOptim, typename T>
class GOF_updater
{
        private:
                std::vector<T> last_lambda_derivatives;

                inline void call_from_to(UInt start, UInt finish, T lambda, LambdaOptim * lopt_ptr)
                {
                        for(UInt i=start; i<=finish; ++i)
                        {
                                lopt_ptr->updaters(i, lambda);
                                last_lambda_derivatives[i] = lambda;
                        }
                }

        public:
                GOF_updater(void) = default;

                inline void initialize(const std::vector<T> & first_lambdas)
                {
                        last_lambda_derivatives = first_lambdas;
                }

                inline void call_to(UInt finish, T lambda, LambdaOptim * lopt_ptr)
                {
                        bool found = false;
                        for(UInt i = 0; i<=finish && found==false; ++i)
                                if(lambda != last_lambda_derivatives[i])
                                {
                                        call_from_to(i, finish, lambda, lopt_ptr);
                                        found = true;
                                }
                }
};

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

        static void set_T_nW_a(SpMat & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp);
        static void set_T_W_a(SpMat & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp);
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
        static void set_E_W_a(SpMat & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap);
        static void set_E_nW_a(SpMat & E, const SpMat * psi_tp, const VectorXr * Ap);

        static void set_z_hat_W(VectorXr & z_hat, const MatrixXr * Hp, const MatrixXr * Qp, const SpMat & S, const VectorXr * zp);
        static void set_z_hat_nW(VectorXr & z_hat, const SpMat & S, const VectorXr * zp);

        template<typename InputCarrier>
        static typename std::enable_if< std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, UInt nnodes);

        template<typename InputCarrier>
        static typename std::enable_if< std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
                universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, UInt nnodes);
};


template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_T_setter(SpMat & T, const InputCarrier & carrier)
        {
                const VectorXr * Ap = carrier.get_Ap();
                const SpMat * psip = carrier.get_psip();
                const SpMat * psi_tp = carrier.get_psi_tp();

                if (carrier.has_W())
                {
                        const MatrixXr * Qp = carrier.get_Qp();
                        AuxiliaryOptimizer::set_T_W_a(T, Ap, psip, psi_tp, Qp);
                }
                else
                {
                        AuxiliaryOptimizer::set_T_nW_a(T, Ap, psip, psi_tp);
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
                        // Psi is full && Q != I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        const MatrixXr * Qp = carrier.get_Qp();
                        AuxiliaryOptimizer::set_E_W_a(E, psi_tp, Qp, Ap);

                }
                else
                {
                        // Psi is full && Q == I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        AuxiliaryOptimizer::set_E_nW_a(E, psi_tp, Ap);
                }
                return 0;
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
                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, UInt nnodes)
        {
                if (carrier.has_W())
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*(carrier.get_Ap()->asDiagonal())*(*carrier.get_Qp())*US;
                }
                else
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*(carrier.get_Ap()->asDiagonal())*US;
                }
                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, UInt nnodes)
        {
                if (carrier.has_W())
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*(*carrier.get_Qp())*US;
                }
                else
                {
                        b.topRows(nnodes) = (*carrier.get_psi_tp())*US;
                }
                return 0;
        }

#endif
