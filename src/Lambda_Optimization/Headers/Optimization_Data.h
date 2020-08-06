#ifndef __OPTIMIZATIONDATA_HPP__
#define __OPTIMIZATIONDATA_HPP__

#include "../../FdaPDE.h"
#include <string>

//!  Class to collect data for optimization
/*!
 * This class collects all the data used in the optiization framework
*/
class  OptimizationData
{
        private:
                std::string criterion     = "batch";            //!< Batch, newton or newton_fd
                std::string DOF_evaluation = "not_required";     //!< not_required, stochastic or exact
                std::string loss_function = "unused";           //!< unused or GCV

                // For batch
                std::vector<Real> lambda_S = {-1.};
                std::vector<Real> lambda_T = {-1.};
                UInt size_S;
                UInt size_T;

                // For batch fixed method
                UInt best_lambda_S = 0;	        //!< Stores the index of the best lambdaS according to method
                UInt best_lambda_T = 0;	        //!< Stores the index of the best lambdaT according to method
                Real best_value    = std::numeric_limits<Real>::max();	//!< Stores the value of the best loss function

                // For optimized methods
                Real initial_lambda_S = 0.0;
                Real initial_lambda_T = 0.0;
                UInt seed             = 0;
                UInt nrealizations    = 0;      //!< The number of random points used in the stochastic computation of the dofs

                // To keep track of optimization
                Real last_lS_used;
                Real last_lT_used;

                // If already present
                MatrixXr DOF_matrix;

                // Tuning parameter
                Real tuning = 1.;

                // For GAM
                Real current_lambdaS;
                std::vector<Real> lambdaS_backup;


                void builder_utility(SEXP Roptim, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune);
                void fill_lambda(SEXP Rlambda, std::vector<Real> & vect, UInt & size);
                void initialize_lambda(SEXP Rlambda, Real & init);

        public:
                OptimizationData() = default;

                // [[TODO ADD SEED TO THE CONSTRUCTORS]]
                OptimizationData(SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune);
                OptimizationData(SEXP Roptim, SEXP Rlambda_S, SEXP Rlambda_T, SEXP Rflag_parabolic, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune);

                // Setters
                inline void set_criterion(const std::string && criterion_) {criterion = criterion_;}
                inline void set_DOF_evaluation(const std::string && DOF_evaluation_) {DOF_evaluation = DOF_evaluation_;}
                inline void set_loss_function(const std::string && loss_function_) {loss_function = loss_function_;}
                inline void set_lambda_S(const std::vector<Real> & lambda_S_) {lambda_S = lambda_S_;}
                inline void set_lambda_T(const std::vector<Real> & lambda_T_) {lambda_T = lambda_T_;}
                inline void set_best_lambda_S(const UInt best_lambda_S_) {best_lambda_S = best_lambda_S_;}
                inline void set_best_lambda_T(const UInt best_lambda_T_) {best_lambda_T = best_lambda_T_;}
                inline void set_best_value(const Real best_value_) {best_value = best_value_;}
                inline void set_initial_lambda_S(const Real initial_lambda_S_) {initial_lambda_S = initial_lambda_S_;}
                inline void set_initial_lambda_T(const Real initial_lambda_T_) {initial_lambda_T = initial_lambda_T_;}
                inline void set_seed(const UInt seed_){seed = seed_;}
                inline void set_nrealizations(const UInt nrealizations_) {nrealizations = nrealizations_;}
                inline void set_last_lS_used(const Real last_lS_used_) {last_lS_used = last_lS_used_;}
                inline void set_last_lT_used(const Real last_lT_used_) {last_lT_used = last_lT_used_;}
                inline void set_DOF_matrix(const MatrixXr & DOF_matrix_) {DOF_matrix = DOF_matrix_;}
                inline void set_tuning(const Real tuning_) {tuning = tuning_;}
                inline void set_current_lambdaS(const Real new_lambdaS) {current_lambdaS = new_lambdaS;}
                inline void setCurrentLambda(UInt lambda_index) {lambda_S = std::vector<Real>(1,lambdaS_backup[lambda_index]);}
                inline void set_lambdaS_backup(void) {lambdaS_backup = lambda_S;}

                // Getters
                inline std::string get_criterion(void) const {return criterion;}
                inline std::string get_DOF_evaluation(void) const {return DOF_evaluation;}
                inline std::string get_loss_function(void) const {return loss_function;}
                inline std::vector<Real> get_lambda_S(void) const {return lambda_S;}
                inline std::vector<Real> get_lambda_T(void) const {return lambda_T;}
                inline UInt get_size_S(void) const {return lambda_S.size();}
                inline UInt get_size_T(void) const {return lambda_T.size();}
                inline UInt get_best_lambda_S(void) const {return best_lambda_S;}
                inline UInt get_best_lambda_T(void) const {return best_lambda_T;}
                inline Real get_best_value(void) const {return best_value;}
                inline Real get_initial_lambda_S(void) const {return initial_lambda_S;}
                inline Real get_initial_lambda_T(void) const {return initial_lambda_T;}
                inline UInt get_seed(void) const {return seed;}
                inline UInt get_nrealizations(void) const {return nrealizations;}
                inline Real get_last_lS_used(void) const {return last_lS_used;}
                inline Real get_last_lT_used(void) const {return last_lT_used;}
                inline MatrixXr const & get_DOF_matrix(void) const {return DOF_matrix;}
                inline Real get_tuning(void) const {return tuning;}
                inline Real get_current_lambdaS(void) const {return current_lambdaS;}
                inline const std::vector<Real> * get_LambdaS_vector() const {return &lambdaS_backup;}

                // Debugging
                void print_opt_data(void) const;
};


//\param DOF an R boolean indicating whether dofs of the model have to be computed or not
//\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)

#endif
