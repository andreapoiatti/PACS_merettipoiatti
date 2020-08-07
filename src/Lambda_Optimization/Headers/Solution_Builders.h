#ifndef __SOLUTION_BUILDERS_H__
#define __SOLUTION_BUILDERS_H__

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"

struct output_Data_opt
{
        VectorXr        z_hat_opt;                          //!< Model predicted values in the locations
        Real            SS_res_opt;                         //!< Model predicted sum of squares of the residuals
        Real            sigma_hat_sq_opt;                   //!< Model estimated variance of errors
        Real            lambda_opt;                        //!< Lambda obtained in the solution
        Real            GCV_opt;                           //!<GCV optimal comptued in the vector of lambdas
        std::vector<Real> GCV_evals;                       //!< GCV evaluations vector
};

//Output struct to be used to return values in R
struct output_Data
{
        std::string     content{"Empty"};          //!< Suggests what the output is containing and how it should be used
        VectorXr        z_hat;                     //!< Model predicted values in the locations
        Real            SS_res = 0.0;              //!< Model predicted sum of squares of the residuals
        Real            rmse = 0.0;                //!< Model root mean squared error
        Real            sigma_hat_sq = 0.0;        //!< Model estimated variance of errors
        Real            dof = 0.0;                 //!< tr(S) + q, degrees of freedom of the model
        Real            dor = 0.0;                 //!< s - dof, degrees of freedom of the residuals
        Real            lambda_sol = 0.0;          //!<Lambda obratained in the solution
        UInt            n_it = 0;                  //!< Number of iterations for the method
        Real            time_partial = 0;          //!<Time, from beginning to end of the optimization method
};

namespace Solution_Builders  //Unique class to manage the output
{
        SEXP GCV_Newton_sol(const VectorXr & solution, const output_Data & output);
        SEXP GCV_batch_sol(const VectorXr & solution, const output_Data_opt & output, const timespec & partial_t);
}

#endif
