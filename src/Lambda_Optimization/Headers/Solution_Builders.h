#ifndef __SOLUTION_BUILDERS_H__
#define __SOLUTION_BUILDERS_H__

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "../../Mesh/Headers/Mesh.h"
#include "../../Regression/Headers/RegressionData.h"

/*struct output_Data_opt
{
        VectorXr        z_hat_opt;                          //!< Model predicted values in the locations
        Real            SS_res_opt;                         //!< Model predicted sum of squares of the residuals
        Real            sigma_hat_sq_opt;                   //!< Model estimated variance of errors
        Real            lambda_opt;                        //!< Lambda obtained in the solution
        Real            GCV_opt;                           //!<GCV optimal comptued in the vector of lambdas
        std::vector<Real> GCV_evals;                       //!< GCV evaluations vector
};*/

//Output struct to be used to return values in R
struct output_Data
{
        std::string     content{"Empty"};          //!< Suggests what the output is containing and how it should be used
        VectorXr        z_hat;                     //!< Model predicted values in the locations
        Real            rmse = 0.0;                //!< Model root mean squared error
        Real            sigma_hat_sq = 0.0;        //!< Model estimated variance of errors
        std::vector<Real> dof = {};             //!< tr(S) + q, degrees of freedom of the model
        Real            lambda_sol = 0.0;          //!<Lambda obratained in the solution
        UInt            lambda_pos=0;              //!<Position of optimal lambda, only for batch evaluation, in R numebring starting from 1 (0 means no batch used)
        UInt            n_it = 0;                  //!< Number of iterations for the method
        Real            time_partial = 0.0;          //!<Time, from beginning to end of the optimization method
        std::vector<Real> GCV_evals={-1};          //!< GCV evaluations vector of explored lambda, with the optimization iterative method or batch
        std::vector<Real> lambda_vec={-1};          //!<Vector of explored lambda with with the optimization iterative method or batch
        Real            GCV_opt=-1;                //!<GCV optimal comptued in the vector of lambdas
        int            termination=-2;             //!< Reason of termination of the iterative optimization method (reached tolerance or max number of iterations)
        MatrixXv       betas;                      //! beta coefficients of the optimal solution
};

//namespace Solution_Builders  //Unique class to manage the output
//{
        //template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
        //SEXP build_sol(const VectorXr & solution, const output_Data & output, const MeshHandler<ORDER, mydim, ndim> & mesh, const InputHandler & regressionData );
        //SEXP GCV_batch_sol(const VectorXr & solution, const output_Data_opt & output, const timespec & partial_t);
//}

#endif
