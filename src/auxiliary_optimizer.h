#ifndef __AUXILIARYOPTIMIZER_HPP__
#define __AUXILIARYOPTIMIZER_HPP__

#include <functional>
#include "fdaPDE.h"
#include "carrier.h"
#include "solver.h"
#include "solverdefinitions.h"

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

template<typename InputCarrier>
void AuxiliaryData<InputCarrier, typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value>::type>::left_multiply_by_psi(const InputCarrier & carrier, VectorXr & ret, const VectorXr & vec)
{
        if (carrier.loc_are_nodes())
        {
                const UInt s = carrier.get_n_obs();
                ret = VectorXr::Zero(s);

                const std::vector<UInt> * kp = carrier.get_obs_indicesp();
                for (UInt i = 0; i < s; i++)
                        for (UInt j = 0; j < s; j++)
                                        ret.coeffRef(i) += vec.coeff((*kp)[i]);
        }
        else
        {
                // Psi is full
                ret = (*carrier.get_psip())*vec;
        }
}


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

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                const SpMat * R1p_= carrier.get_R1p();         // Get the value of matrix R1
                //Sparse_LU solver;
                Eigen::LLT<MatrixXr> factorized_R0p(MatrixXr(*(carrier.get_R0p())));	                                 // define a factorized empty sparse Cholesky solver  [[LDLT???]]
                //solver.compute(*(carrier.get_R0p()));		 // apply it to R0 to simplify the inverse
                R = (*R1p_).transpose()*factorized_R0p.solve(MatrixXr(*R1p_));            // R == _R1^t*R0^{-1}*R1
                adt.f_ = ((*R1p_).transpose())*factorized_R0p.solve((MatrixXr(*carrier.get_up())));

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_R_setter(MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                const SpMat * R1p_= carrier.get_R1p();         // Get the value of matrix R1
                //Sparse_LU solver;	                                 // define a factorized empty sparse Cholesky solver  [[LDLT???]]
                //solver.compute(*(carrier.get_R0p()));		 // apply it to R0 to simplify the inverse
                Eigen::LLT<MatrixXr> factorized_R0p(MatrixXr(*(carrier.get_R0p())));

                R = (*R1p_).transpose()*factorized_R0p.solve(MatrixXr(*R1p_));            // R == _R1^t*R0^{-1}*R1

                return 0;
        }


template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_T_setter(MatrixXr & T, const InputCarrier & carrier)
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
        AuxiliaryOptimizer::universal_T_setter(MatrixXr & T, const InputCarrier & carrier)
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
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                Eigen::LDLT<MatrixXr> factorized_T(T);	        // define a factorization

                if(!carrier.is_areal() && !carrier.has_W())
                {

                        // Q == I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        V = factorized_T.solve(MatrixXr(*psi_tp));      // find the value of V = T^{-1}*Psi^t
                }
                else
                {
                        MatrixXr E_;                // Declare an empty auxiliary matrix
                        const UInt ret =  AuxiliaryOptimizer::universal_E_setter<InputCarrier>(E_, carrier);
                        V = factorized_T.solve(E_);          // find the value of V = T^{-1}*E
                }
                adt.K_ = factorized_T.solve(R);              // K = T^{-1}*R
                adt.g_ = factorized_T.solve(adt.f_);

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_V_setter(MatrixXr & V, const MatrixXr & T, const MatrixXr & R, const InputCarrier & carrier, AuxiliaryData<InputCarrier> & adt)
        {
                Eigen::LDLT<MatrixXr> factorized_T(T);	        // define a factorization

                if(!carrier.is_areal() && !carrier.has_W())
                {

                        // Q == I
                        const SpMat * psi_tp = carrier.get_psi_tp();
                        V = factorized_T.solve(MatrixXr(*psi_tp));      // find the value of V = T^{-1}*Psi^t
                }
                else
                {
                        MatrixXr E_;                // Declare an empty auxiliary matrix
                        const UInt ret =  AuxiliaryOptimizer::universal_E_setter<InputCarrier>(E_, carrier);
                        V = factorized_T.solve(E_);          // find the value of V = T^{-1}*E
                }
                adt.K_ = factorized_T.solve(R);              // K = T^{-1}*R

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_E_setter(MatrixXr & E, const InputCarrier & carrier)
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
        AuxiliaryOptimizer::universal_E_setter(MatrixXr & E, const InputCarrier & carrier)
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
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_z_hat_setter(VectorXr & z_hat, const InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda)
        {
                common_z_hat_part(z_hat, carrier, S);

                adt.left_multiply_by_psi(carrier, adt.r_, adt.g_);

                if (carrier.has_W())
                {
                        const MatrixXr * Qp = carrier.get_Qp();
                        adt.r_ = lambda*(*Qp)*adt.r_;
                }
                else
                {
                        adt.r_ = lambda * adt.r_;
                }

                z_hat += adt.r_;

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_z_hat_setter(VectorXr & z_hat, const InputCarrier & carrier, const MatrixXr & S, AuxiliaryData<InputCarrier> & adt, const Real lambda)
        {
                common_z_hat_part(z_hat, carrier, S);

                return 0;
        }

template<typename InputCarrier>
void AuxiliaryOptimizer::common_z_hat_part(VectorXr & z_hat, const InputCarrier & carrier, const MatrixXr & S)
{
        const VectorXr * zp = carrier.get_zp();
        if(carrier.has_W())
        {
                const MatrixXr * Hp = carrier.get_Hp();
                const MatrixXr * Qp = carrier.get_Qp();
                AuxiliaryOptimizer::set_z_hat_W(z_hat, Hp, Qp, S, zp);
        }
        else
        {
                AuxiliaryOptimizer::set_z_hat_nW(z_hat, S, zp);
        }
}

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Areal, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, const UInt nnodes)
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
        AuxiliaryOptimizer::universal_b_setter(MatrixXr & b, const InputCarrier & carrier, const MatrixXr & US, const UInt nnodes)
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

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                adt.t_ = dS*(*zp);
                MatrixXr temp = lambda*adt.K_;
                for (UInt i=0; i<temp.cols(); i++)
                {
                        temp.coeffRef(i,i) -= 1;
                }
                adt.h_ = temp*adt.g_;
                adt.left_multiply_by_psi(carrier, adt.p_, adt.h_);
                adt.a_ = eps.transpose()*adt.p_;  // note different from previous case!!

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_first_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & dS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                adt.t_ = dS*(*zp);
                adt.a_ = -eps.transpose()*adt.t_;

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,t_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_second_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                if (carrier.has_W())
                        adt.b_ = adt.p_.transpose()*(*carrier.get_Qp())*adt.p_;
                else
                        adt.b_ = adt.p_.squaredNorm();

                VectorXr aux;
                adt.left_multiply_by_psi(carrier, aux, -2*adt.K_*adt.h_);

                adt.c_ = eps.transpose()*(-ddS*(*zp) + aux);

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<multi_bool_type<std::is_base_of<Forced, InputCarrier>::value>,f_type>::value, UInt>::type
        AuxiliaryOptimizer::universal_second_updater(AuxiliaryData<InputCarrier> & adt, const InputCarrier & carrier, const MatrixXr & ddS, const VectorXr & eps, const Real lambda)
        {
                const VectorXr * zp = carrier.get_zp();
                if (carrier.has_W())
                        adt.b_ = adt.t_.transpose()*(*carrier.get_Qp())*adt.t_;
                else
                        adt.b_ = adt.t_.squaredNorm();
                adt.c_ = -eps.transpose()*ddS*(*zp);

                return 0;
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
        AuxiliaryOptimizer::universal_GCV(const Real s, const Real sigma_hat_sq, const Real dor)
        {
                return s*sigma_hat_sq/Real(dor);
        }

template<typename InputCarrier>
typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
        AuxiliaryOptimizer::universal_GCV_d(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS)
        {
                return 2*s*(sigma_hat_sq*trdS + adt.a_)/Real(dor*dor);
        }


template<typename InputCarrier>
typename std::enable_if<std::is_same<t_type,t_type>::value, Real>::type
        AuxiliaryOptimizer::universal_GCV_dd(const AuxiliaryData<InputCarrier> & adt, const Real s, const Real sigma_hat_sq, const Real dor, const Real trdS, const Real trddS)
        {
                return 2*s*(trdS*(3*sigma_hat_sq*trdS+4*adt.a_)/dor + sigma_hat_sq*trddS + adt.b_ + adt.c_)/Real(dor*dor);
        }


#endif
