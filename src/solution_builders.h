#ifndef __SOLUTIONBUILDERS_H__
#define __SOLUTIONBUILDERS_H__

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "fdaPDE.h"
#include "solver.h"
#include "function_variadic.h"
#include "vector_eval.h"


namespace Solution_builders  //Unique class to manage the output
{
        SEXP GCV_Newton_sol(const VectorXr & solution, const output_Data & output);
        SEXP GCV_batch_sol(const VectorXr & solution, const output_Data_opt & output, const timespec & partial_t);
}


#endif
