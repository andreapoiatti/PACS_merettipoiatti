#ifndef __CARRIER_HPP__
#define __CARRIER_HPP__

// HEADERS
#include <memory>
#include <type_traits>
#include "fdaPDE.h"
#include "mixedFERegression.h"
#include "optimizationData.h"
#include "regressionData.h"

// Declaration of classes that will be used as Extensions for Carrier
class Areal;
class Forced;
class Temporal;
//----------------------------------------------------------------------------//

// *** CARRIER CLASS ***
//! Plain Carrier class [inherits specializations via variadic template]
/*!
 This class contains all the information needed by optimization methods and is
 conceived in a pattern-based fashion: inheriting additional structures according to
 the type of optimization problem that has to be solved. Besides the main purpose
 of being a clean and all-inclusive transfer class, this aslo provides two apply methods
 that can be used to call back the regression problem and solve the main system.
 \tparam InputHandler the typy of regression problem from which to build the Carrier
 \tparam Extensions... further information that has to be included into the Carrier
*/
template<typename InputHandler, typename... Extensions>
class Carrier: public Extensions...
{
        private:
                // --- DATA ---
                const OptimizationData * opt_data;            //!< Pointer to the optimization data needed for the method
                MixedFERegressionBase<InputHandler> * model;  //!< Pointer to the model data [not const, since it must modify members with the apply]

                // Booleans for particular computations
                bool locations_are_nodes = false;             //!< locations are the nodes boolean
                bool has_covariates = false;                  //!< regression problem has covariates boolean
                bool areal_data = false;                      //!< areal problem boolean [depends on inheritance]
                bool boundary_conditions = false;             //!< boundary conditions present boolean [depends on inheritance]
                bool forced_data = false;                     //!< presence of forcing term boolean [depends on inheritance]
                bool temporal_data = false;                   //!< presence of temporal data  [depends on inheritance]

                // General data for any problem
                UInt n_obs;                                   //!< number of locations
                UInt n_nodes;                                 //!< number of nodes
                const std::vector<UInt> * obs_indicesp;       //!< ndeices of the locations for efficient search and storage

                const VectorXr * zp;                          //!< pointer to the observations in the locations [size n_obs]
                const MatrixXr * Wp;                          //!< pointer to the matrix of covariates [size n_obs x n_covariates]
                const MatrixXr * Hp;                          //!< pointer to the hat matrix [size n_covariates x n_covariates]
                const MatrixXr * Qp;                          //!< pointer to the identity - (hat matrix) [size n_covariates x n_covariates]

                const SpMat * DMatp;                          //!< pointer to the north-west block of system matrix [size n_nodes x n_nodes]
                const SpMat * R1p;                            //!< pointer to R1 matrix [size n_nodes x n_nodes]
                const SpMat * R0p;                            //!< pointer to R0 matrix [size n_nodes x n_nodes]
                const SpMat * psip;                           //!< pointer to location-to-nodes matrix [size n_obs x n_nodes]
                const SpMat * psi_tp;                         //!< pointer to the transpose of the location-to-nodes matrix [size n_nodes x n_obs]

                const VectorXr * rhsp;                        //!< pointer to the right hand side of the system

                const std::vector<Real> * bc_valuesp;         //!< pointer to the boundary conditions vector of values
                const std::vector<UInt> * bc_indicesp;        //!< pointer to the boundary conditions vector of indices


        public:
                // CONSTRUCTORS
                //! Default constructor
                Carrier() = default;

                //! Constructor taking object Extensions and initializing with it the new class [used for minimal inheritance]
                /*!
                 \param ext eventual extensions to be added to the plain Carrier
                */
                template<typename Dummy = typename std::enable_if<sizeof...(Extensions)!=0, void>, typename... Bricks>
                Carrier(Bricks && ... ext): Extensions(std::forward<Bricks>(ext))...{};

                //! Universal setter of the class: fills all the plain parameters [all pointers are const except for the model]
                /*!
                 \param model_ pointer of the miexd object from which the carrier is derived, used for apply(s) purpose
                 \param opt_data_ stores the data related to the optimization procedure to be followed
                 \parm locations_are_nodes_ boolean to check if locations are nodes [for simplified computations]
                 \param has_covariates_ boolean to check if the problem has regressors [for simplified computations]
                 \param n_obs_ number of locations and hence of observations
                 \param n_noes_number of nodes of the mesh
                 \param obs_indicesp_ pointer collectig the indices of the getObservations
                 \param zp_ pointer to the observations in the locations
                 \param Wp_ pointer to the matrix of covariates
                 \param Hp_ pointer to the hat matrix
                 \param Qp_ pointer to identity - hat matrix
                 \param DMatp_ pointer to the north-west blockk of the system matrix
                 \param R1p_ pointer to R1 matrix
                 \param R0p_ pointer to R0 matrix
                 \param psip_ pointer to Psi matrix
                 \param psi_tp_ pointer to Psi^T matrix
                 \parma rhsp_ pointer to the right hand side of system matrix
                 \param bc_values_ pointer to the values of boundary conditions
                 \param bc_indicesp_ pointer tp the indices of the boundary conditions
                */
                inline void set_all(MixedFERegressionBase<InputHandler> * model_, const OptimizationData * opt_data_,
                        bool locations_are_nodes_, bool has_covariates_, UInt n_obs_, UInt n_nodes_, const std::vector<UInt> * obs_indicesp_,
                        const VectorXr * zp_, const MatrixXr * Wp_, const MatrixXr * Hp_, const MatrixXr * Qp_,
                        const SpMat * DMatp_, const SpMat * R1p_, const SpMat * R0p_, const SpMat * psip_, const SpMat * psi_tp_,
                        const VectorXr * rhsp_, const std::vector<Real> * bc_valuesp_, const std::vector<UInt> * bc_indicesp_)
                {
                        // Set all the data through the private setters
                        set_model(model_);
                        set_opt_data(opt_data_);
                        set_loc_are_nodes(locations_are_nodes_);
                        set_has_W(has_covariates_);
                        set_n_obs(n_obs_);
                        set_n_nodes(n_nodes_);
                        set_obs_indicesp(obs_indicesp_);
                        set_zp(zp_);
                        set_Wp(Wp_);
                        set_Hp(Hp_);
                        set_Qp(Qp_);
                        set_DMatp(DMatp_);
                        set_R1p(R1p_);
                        set_R0p(R0p_);
                        set_psip(psip_);
                        set_psi_tp(psi_tp_);
                        set_rhsp(rhsp_);
                        set_bc_valuesp(bc_valuesp_);
                        set_bc_indicesp(bc_indicesp_);

                        // Update the booleans [note some consistency constraints]
                        if (bc_indicesp_->size() > 0)
                        {
                                this->boundary_conditions = true;
                        }

                        if(std::is_base_of<Areal, Carrier>::value)
                        {
                                this->areal_data = true;
                                this->locations_are_nodes = false; // Areal data can't hav locaions by nodes
                                this->boundary_conditions = false; // Areal data can't have boundary conditions
                        }
                        if(std::is_base_of<Forced, Carrier>::value)
                                this->forced_data = true;

                        if(std::is_base_of<Temporal, Carrier>::value)
                                this->temporal_data = true;
                }

                // GETTERS
                inline const OptimizationData * get_opt_data(void) const {return this->opt_data;}               //!< Getter of opt_data
                inline bool loc_are_nodes(void) const {return this->locations_are_nodes;}                       //!< Getter of locations_are_nodes
                inline bool has_W(void) const {return this->has_covariates;}                                    //!< Getter of has_covariates
                inline bool is_areal(void) const {return this->areal_data;}                                     //!< Getter of areal_data
                inline bool is_temporal(void) const {return this->temporal_data;}                               //!< Getter of temporal_data
                inline UInt get_n_obs(void) const {return this->n_obs;}                                         //!< Getter of n_obs [# locations]
                inline UInt get_n_nodes(void) const {return this->n_nodes;}                                     //!< Getter of n_nodes [# nodes]
                inline const std::vector<UInt> * get_obs_indicesp(void) const {return this->obs_indicesp;}      //!< Getter of obs_indicesp
                inline const VectorXr * get_zp(void) const {return this->zp;}                                   //!< Getter of zp
                inline const MatrixXr * get_Wp(void) const {return this->Wp;}                                   //!< Getter of Wp
                inline const MatrixXr * get_Hp(void) const {return this->Hp;}                                   //!< Getter of Hp
                inline const MatrixXr * get_Qp(void) const {return this->Qp;}                                   //!< Getter of Qp
                inline const SpMat * get_DMatp(void) const {return this->DMatp;}                                //!< Getter of DMatp
                inline const SpMat * get_R1p(void) const {return this->R1p;}                                    //!< Getter of R1p
                inline const SpMat * get_R0p(void) const {return this->R0p;}                                    //!< Getter of R0p
                inline const SpMat * get_psip(void) const {return this->psip;}                                  //!< Getter of psip
                inline const SpMat * get_psi_tp(void) const {return this->psi_tp;}                              //!< Getter of psi_tp
                inline const VectorXr * get_rhsp(void) const {return this->rhsp;}                               //!< Getter of rhsp
                inline const std::vector<Real> * get_bc_valuesp(void) const {return this->bc_valuesp;}          //!< Getter of bc_valuesp
                inline const std::vector<UInt> * get_bc_indicesp(void) const {return this->bc_indicesp;}        //!< Getter of bc_indicesp


                // SETTERS
                inline void set_model(MixedFERegressionBase<InputHandler> * md) {this->model=md;};                                      //!< Setter of model
                inline void set_opt_data(const OptimizationData * opt_data_) {this->opt_data = opt_data_;}                              //!< Setter of opt_data
                inline void set_loc_are_nodes(const bool locations_are_nodes_) {this->locations_are_nodes = locations_are_nodes_;}      //!< Setter of locations_are_nodes
                inline void set_has_W(const bool has_covariates_) {this->has_covariates = has_covariates_;}                             //!< Setter of has_covariates
                inline void set_n_obs(const UInt n_obs_) {this->n_obs = n_obs_;}                                                        //!< Setter of n_obs
                inline void set_n_nodes(const UInt n_nodes_) {this->n_nodes = n_nodes_;}                                                //!< Setter of n_nodes
                inline void set_obs_indicesp(const std::vector<UInt> * obs_indicesp_) {this->obs_indicesp = obs_indicesp_;}             //!< Setter of obs_indicesp
                inline void set_zp(const VectorXr * zp_) {this->zp = zp_;}                                                              //!< Setter of zp
                inline void set_Wp(const MatrixXr * Wp_) {this->Wp = Wp_;}                                                              //!< Setter of Wp
                inline void set_Hp(const MatrixXr * Hp_) {this->Hp = Hp_;}                                                              //!< Setter of Hp
                inline void set_Qp(const MatrixXr * Qp_) {this->Qp = Qp_;}                                                              //!< Setter of Qp
                inline void set_DMatp(const SpMat * DMatp_) {this->DMatp = DMatp_;}                                                     //!< Setter of DMatp
                inline void set_R1p(const SpMat * R1p_) {this->R1p = R1p_;}                                                             //!< Setter of R1p
                inline void set_R0p(const SpMat * R0p_) {this->R0p = R0p_;}                                                             //!< Setter of R0p
                inline void set_psip(const SpMat * psip_) {this->psip = psip_;}                                                         //!< Setter of psip
                inline void set_psi_tp(const SpMat * psi_tp_) {this->psi_tp = psi_tp_;}                                                 //!< Setter of psi_tp
                inline void set_rhsp(const VectorXr * rhsp_) {this->rhsp = rhsp_;}                                                      //!< Setter of rhsp
                inline void set_bc_valuesp(const std::vector<Real> * bc_valuesp_) {this->bc_valuesp = bc_valuesp_;}                     //!< Setter of bc_valuesp
                inline void set_bc_indicesp(const std::vector<UInt> * bc_indicesp_) {this->bc_indicesp = bc_indicesp_;}                 //!< Setter of bc-bc_indicesp

                // APPLY FUNCTIONS
                // Apply is called here in order not to make the non const pointer public
                //! Method to the system given a lambda and a right hand side of the system
                /*!
                 \param b the right hand side of the system to be solved vi system matrix
                 \param lambda the optimization parameter with which to build the system matrix
                 \return the solution of the system
                */
                inline MatrixXr apply_to_b(const MatrixXr & b, Real lambda){
                        this->model->set_lambda(lambda); // set the lambda value
                        return (this->model->apply_to_b(b))(0,0);} //specific for spatial case

                //! Method to the system given a lambda [rith hand sde is the usual of the problem]
                /*!
                 \param lambda the optimization parameter with which to build the system matrix
                 \return the solution of the system
                */
                inline MatrixXr apply(Real lambda) {
                        this->model->set_lambda(lambda); // set the lambda value
                        return (this->model->apply())(0,0);}
};
//----------------------------------------------------------------------------//

// *** CARRIER EXTENSIONS ***
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
//----------------------------------------------------------------------------//

// *** CARRIER BUILDER ***
template<typename DataHandler>
class CarrierBuilder
{
        private:
                template<typename... Extensions>
                static void set_plain_data(Carrier<DataHandler, Extensions...> & car, const DataHandler & data, MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {
                        car.set_all(&mc, &optimizationData, data.isLocationsByNodes(), bool(data.getCovariates()->rows()>0 && data.getCovariates()->cols()>0),
                                data.getNumberofObservations(), mc.getnnodes_(), data.getObservationsIndices(),
                                data.getObservations(), data.getCovariates(), mc.getH_(), mc.getQ_(), mc.getDMat_(), mc.getR1_(),
                                mc.getR0_(), mc.getpsi_(), mc.getpsi_t_(), mc.getrhs_(), data.getDirichletValues(), data.getDirichletIndices());
                }

                template<typename... Extensions>
                static void set_areal_data(Carrier<DataHandler, Extensions...> & car, const DataHandler & data, MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {
                        car.set_all_areal(data.getNumberOfRegions(), mc.getA_());
                }

                template<typename... Extensions>
                static void set_forced_data(Carrier<DataHandler, Extensions...> & car, const DataHandler & data, MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {
                        car.set_all_forced(mc.getu_());
                }

        public:
                static Carrier<DataHandler> build_plain_carrier(const DataHandler & data, MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {       //da modificare, non serve copiarli due volte!
                        Carrier<DataHandler> car;
                        set_plain_data(car, data, mc, optimizationData);

                        return car;
                }

                static Carrier<DataHandler, Areal> build_areal_carrier(const DataHandler & data,  MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {
                        Carrier<DataHandler, Areal> car;
                        set_plain_data(car, data, mc, optimizationData);
                        set_areal_data(car, data, mc, optimizationData);

                        return car;
                }

                static Carrier<DataHandler, Forced> build_forced_carrier(const DataHandler & data, MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {
                        Carrier<DataHandler, Forced> car;
                        set_plain_data(car, data, mc, optimizationData);
                        set_forced_data(car, data, mc, optimizationData);

                        return car;
                }

                static Carrier<DataHandler, Forced, Areal> build_forced_areal_carrier(const DataHandler & data, MixedFERegressionBase<DataHandler> & mc, const OptimizationData & optimizationData)
                {
                        Carrier<DataHandler, Forced, Areal> car;
                        set_plain_data(car, data, mc, optimizationData);
                        set_areal_data(car, data, mc, optimizationData);
                        set_forced_data(car, data, mc, optimizationData);

                        return car;
                }
};

#endif
