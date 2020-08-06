#ifndef __SOLUTION_BUILDERS_H__
#define __SOLUTION_BUILDERS_H__

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "Function_Variadic.h"
#include "Newton.h"
#include "Batch_Evaluator.h"
#include "Auxiliary_Optimizer.h"

namespace Solution_Builders  //Unique class to manage the output
{
        SEXP GCV_Newton_sol(const VectorXr & solution, const output_Data & output);
        SEXP GCV_batch_sol(const VectorXr & solution, const output_Data_opt & output, const timespec & partial_t);
}


#endif
