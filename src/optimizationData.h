#ifndef __OPTIMIZATIONDATA_HPP__
#define __OPTIMIZATIONDATA_HPP__

#include "fdaPDE.h"
#include <string>
#include <vector>

//!  Class to collect data for optimization
/*!
 * This class collects all the data used in the optiization framework
*/
class  OptimizationData
{
        private:
                std::string criterion_;          // Batch or optimized
                std::string evaluation_;         // Stochastic or exact
                std::string method_;             // GCV

                std::vector<Real> lambdas_;
                Real initial_lambda_ = 0.0;
                UInt nrealizations_ = 0;      // the number of random points used in the stochastic computation of the dofs

        public:
                OptimizationData() = default;

                // Setters
                inline void set_criterion_ (const std::string && criterion) {criterion_ = criterion;}
                inline void set_evaluation_ (const std::string && evaluation) {evaluation_ = evaluation;}
                inline void set_method_ (const std::string && method) {method_ = method;}
                inline void set_lambdas_(const std::vector<Real> & lambdas) {lambdas_ = lambdas;}
                inline void set_initial_lambda_(const Real initial_lambda){initial_lambda_ = initial_lambda;}
                inline void set_nrealizations_(const UInt nrealizations){nrealizations_=nrealizations;}

                // Getters
                inline std::string get_criterion_(void) const {return criterion_;}
                inline std::string get_evaluation_(void) const {return evaluation_;}
                inline std::string get_method_(void) const {return method_;}
                inline const std::vector<Real> * get_lambdas_ (void) const {return &lambdas_;}
                inline Real get_initial_lambda_(void) const {return initial_lambda_;}
                inline UInt get_nrealizations_(void) const {return nrealizations_;}

                // Debugging
                void print_opt_data(void) const;
};



//\param DOF an R boolean indicating whether dofs of the model have to be computed or not
//\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)

#endif
