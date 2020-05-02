#ifndef __CARRIER_HPP__
#define __CARRIER_HPP__

#include <type_traits>
#include "fdaPDE.h"
#include "optimizationData.h"

// Declaration of classes that will be used as Extensions for Carrier
class Areal;
class Forced;
class Temporal;

template<typename Origin, typename... Extensions>
class Carrier: public Extensions...
{
        private:
                //Basic Data
                Origin * trace;
                const OptimizationData * opt_data;

                bool locations_are_nodes = false;
                bool has_covariates = false;
                bool areal_data = false;
                bool forced_data = false;
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
                Carrier(Bricks && ... ext): Extensions(std::forward<Bricks>(ext))...{};


                inline void set_all(Origin * trace_, const OptimizationData * opt_data_,
                        bool locations_are_nodes_, bool has_covariates_, UInt n_obs_, const std::vector<UInt> * obs_indicesp_,
                        const VectorXr * zp_, const MatrixXr * Wp_, const MatrixXr * Hp_, const MatrixXr * Qp_,
                        const SpMat * R1p_, const SpMat * R0p_, const SpMat * psip_, const SpMat * psi_tp_)
                                {
                                        set_tracep(trace_);
                                        get_opt_data(opt_data_);
                                        set_loc_are_nodes(locations_are_nodes_);
                                        set_has_W(has_covariates_);
                                        set_n_obs(n_obs_);
                                        set_obs_indicesp(obs_indicesp_);
                                        set_zp(zp_);
                                        set_Wp(Wp_);
                                        set_Hp(Hp_);
                                        set_Qp(Qp_);
                                        set_R1p(R1p_);
                                        set_R0p( R0p_);
                                        set_psip(psip_);
                                        set_psi_tp(psi_tp_);

                                        if(std::is_base_of<Areal, Carrier>::value)
                                                this->areal_data = true;
                                        if(std::is_base_of<Forced, Carrier>::value)
                                                this->forced_data = true;
                                        if(std::is_base_of<Temporal, Carrier>::value)
                                                this->temporal_data = true;
                                };

                // Getters
                inline Origin * get_tracep(void) const {return this->trace;}
                inline const OptimizationData * get_opt_data(void) const {return this->opt_data;}
                inline bool loc_are_nodes(void) const {return this->locations_are_nodes;}
                inline bool has_W(void) const {return this->has_covariates;}
                inline bool is_areal(void) const {return this->areal_data;}
                inline bool is_temporal(void) const {return this->temporal_data;}
                inline UInt get_n_obs(void) const {return this->n_obs;}
                inline const std::vector<UInt> * get_obs_indicesp(void) const {return this->obs_indicesp;}
                inline const VectorXr * get_zp(void) const {return this->zp;}
                inline const MatrixXr * get_Wp(void) const {return this->Wp;}
                inline const MatrixXr * get_Hp(void) const {return this->Hp;}
                inline const MatrixXr * get_Qp(void) const {return this->Qp;}
                inline const SpMat * get_R1p(void) const {return this->R1p;}
                inline const SpMat * get_R0p(void) const {return this->R0p;}
                inline const SpMat * get_psip(void) const {return this->psip;}
                inline const SpMat * get_psi_tp(void) const {return this->psi_tp;}

                // Setters
                inline void set_tracep(Origin * trace_) {this->trace = trace_;}
                inline void get_opt_data(const OptimizationData * opt_data_) {this->opt_data = opt_data_;}
                inline void set_loc_are_nodes(const bool locations_are_nodes_) {this->locations_are_nodes = locations_are_nodes_;}
                inline void set_has_W(const bool has_covariates_) {this->has_covariates = has_covariates_;}
                inline void set_n_obs(const UInt n_obs_) {this->n_obs = n_obs_;}
                inline void set_obs_indicesp(const std::vector<UInt> * obs_indicesp_) {this->obs_indicesp = obs_indicesp_;}
                inline void set_zp(const VectorXr * zp_) {this->zp = zp_;}
                inline void set_Wp(const MatrixXr * Wp_) {this->Wp = Wp_;}
                inline void set_Hp(const MatrixXr * Hp_) {this->Hp = Hp_;}
                inline void set_Qp(const MatrixXr * Qp_) {this->Qp = Qp_;}
                inline void set_R1p(const SpMat * R1p_) {this->R1p = R1p_;}
                inline void set_R0p(const SpMat * R0p_) {this->R0p = R0p_;}
                inline void set_psip(const SpMat * psip_) {this->psip = psip_;}
                inline void set_psi_tp(const SpMat * psi_tp_) {this->psi_tp = psi_tp_;}
};

class Areal
{
        private:
                UInt n_regions;
                const VectorXr * Ap;

        public:
                Areal() = default;
                Areal(UInt n_regions_, const VectorXr * Ap_):n_regions(n_regions_), Ap(Ap_) {};

                inline void set_all_areal(UInt n_regions_, const VectorXr * Ap_)
                {
                        set_n_regions(n_regions_);
                        set_Ap(Ap_);
                }

                // Getters
                inline UInt get_n_regions(void) const {return this->n_regions;}
                inline const VectorXr * get_Ap(void) const {return this->Ap;}

                // Setters
                inline void set_n_regions(UInt n_regions_) {this->n_regions = n_regions_;}
                inline void set_Ap(const VectorXr * Ap_) {this->Ap = Ap_;}
};

class Forced
{
        private:
                const VectorXr * up;

        public:
                Forced() = default;
                Forced(const VectorXr * up_): up(up_) {};

                inline void set_all_forced(const VectorXr * up_)
                {
                        set_up(up_);
                }

                // Getters
                inline const VectorXr * get_up(void) const {return this->up;}

                // Setters
                inline void set_up(const VectorXr * up_) {this->up = up_;}
};

class Temporal
{
        // [[ TO BE IMPLEMENTED]]
};

template<typename DataHandler, typename MixedClass>
class CarrierBuilder
{
        public:
                static Carrier<MixedClass> build_plain_carrier(const DataHandler & data, MixedClass & mc, const OptimizationData & optimizationData)
                {
                        Carrier<MixedClass> car;
                        car.set_all(&mc, &optimizationData, mc.check_is_loc_by_n(),
                        	mc.checkisRegression_(), data.getNumberofObservations(), data.getObservationsIndices(),
                        	data.getObservations(), data.getCovariates(), mc.getH_(), mc.getQ_(), mc.getR1_(),
                        	mc.getR0_(), mc.getPsi_(), mc.getPsi_t());

                        return car;
                }

                static Carrier<MixedClass,Areal> build_areal_carrier(const DataHandler & data, MixedClass & mc, const OptimizationData & optimizationData)
                {
                        Carrier<MixedClass,Areal> car;
                        car.set_all(&mc, &optimizationData, mc.check_is_loc_by_n(),
                                mc.checkisRegression_(), data.getNumberofObservations(), data.getObservationsIndices(),
                                data.getObservations(), data.getCovariates(), mc.getH_(), mc.getQ_(), mc.getR1_(),
                                mc.getR0_(), mc.getPsi_(), mc.getPsi_t());
                        car.set_all_areal(data.getNumberOfRegions(), mc.getA_());

                        return car;
                }

                static Carrier<MixedClass,Forced> build_forced_carrier(const DataHandler & data, MixedClass & mc, const OptimizationData & optimizationData)
                {
                        Carrier<MixedClass,Forced> car;
                        car.set_all(&mc, &optimizationData, mc.check_is_loc_by_n(),
                                mc.checkisRegression_(), data.getNumberofObservations(), data.getObservationsIndices(),
                                data.getObservations(), data.getCovariates(), mc.getH_(), mc.getQ_(), mc.getR1_(),
                                mc.getR0_(), mc.getPsi_(), mc.getPsi_t());
                        car.set_all_forced(mc.getu_());

                        return car;
                }

                static Carrier<MixedClass,Forced,Areal> build_forced_areal_carrier(const DataHandler & data, MixedClass & mc, const OptimizationData & optimizationData)
                {
                        Carrier<MixedClass,Forced,Areal> car;
                        car.set_all(&mc, &optimizationData, mc.check_is_loc_by_n(),
                                mc.checkisRegression_(), data.getNumberofObservations(), data.getObservationsIndices(),
                                data.getObservations(), data.getCovariates(), mc.getH_(), mc.getQ_(), mc.getR1_(),
                                mc.getR0_(), mc.getPsi_(), mc.getPsi_t());
                        car.set_all_areal(data.getNumberOfRegions(), mc.getA_());
                        car.set_all_forced(mc.getu_());

                        return car;
                }
};

#endif
