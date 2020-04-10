#ifndef __CARRIER_HPP__
#define __CARRIER_HPP__

#include <type_traits>
#include "fdaPDE.h"
#include "optimizationData.h"

// Declaration of classes that will be used as Extensions for Carrier
class Areal;
class Temporal;

template<typename Origin, typename... Extensions>
class Carrier: public Extensions...
{
        private:
                //Basic Data
                const Origin * trace;
                const OptimizationData * opt_data;

                bool locations_are_nodes = false;
                bool has_covariates = false;
                bool areal_data = false;
                bool temporal_data = false;

                UInt n_obs;

                const std::vector<UInt> * obs_indicesp;
                const VectorXr * zp;
                const MatrixXr * Wp;
                const MatrixXr * Hp;
                const MatrixXr * Qp;

                const SpMat * R1p;
                const SpMat * R0p;
                const SpMat * psip;
                const SpMat * psi_tp;

        public:
                Carrier() = default;

                //! Constructor taking object Extensions and initializing with it the new class
                template<typename Dummy = typename std::enable_if<sizeof...(Extensions)!=0, void>, typename... Bricks>
                Carrier(Bricks && ... ext, const Origin * trace_, const OptimizationData * opt_data_, bool locations_are_nodes_,
                                bool has_covariates_, UInt n_obs_, const std::vector<UInt> * obs_indicesp_,
                                const VectorXr * zp_, const MatrixXr * Wp_, const MatrixXr * Hp_, const MatrixXr * Qp_, const SpMat * R1p_,
                                const SpMat * R0p_, const SpMat * psip_, const SpMat * psi_tp_): trace(trace_), opt_data(opt_data_),
                                locations_are_nodes(locations_are_nodes_), has_covariates(has_covariates_), n_obs(n_obs_),
                                obs_indicesp(obs_indicesp_), zp(zp_), Wp(Wp_), Hp(Hp_), Qp(Qp_), R1p(R1p_), R0p(R0p_),
                                psip(psip_), psi_tp(psi_tp_), Extensions(std::forward<Bricks>(ext)) ...
                {
                        if(std::is_base_of<Areal, Carrier>::value)
                                areal_data = true;
                        if(std::is_base_of<Temporal, Carrier>::value)
                                temporal_data = true;
                };

                inline const Origin * get_tracep(void) const {return trace;}
                inline const OptimizationData * get_opt_data(void) const {return opt_data;}

                inline bool loc_are_nodes(void) const {return locations_are_nodes;}
                inline bool has_W(void) const {return has_covariates;}
                inline bool is_areal(void) const {return areal_data;}
                inline bool is_temporal(void) const {return temporal_data;}

                inline UInt get_n_obs(void) const {return n_obs;}

                inline const std::vector<UInt> * get_obs_indicesp(void) const {return obs_indicesp;};
                inline const VectorXr * get_zp(void) const {return zp;};
                inline const MatrixXr * get_Wp(void) const {return Wp;};
                inline const MatrixXr * get_Hp(void) const {return Hp;};
                inline const MatrixXr * get_Qp(void) const {return Qp;};

                inline const SpMat * get_R1p(void) const {return R1p;};
                inline const SpMat * get_R0p(void) const {return R0p;};
                inline const SpMat * get_psip(void) const {return psip;};
                inline const SpMat * get_psi_tp(void) const {return psi_tp;};
};

class Areal
{
        private:
                MatrixXi incidenceMatrix_;
                UInt nRegions_;
                const VectorXr * Ap;

        public:
                Areal(const VectorXr * Ap_):Ap(Ap_) {};
                inline const VectorXr * get_Ap(void) const {return Ap;};
};

class Temporal
{
        // [[ TO BE IMPLEMENTED]]
};

#endif
