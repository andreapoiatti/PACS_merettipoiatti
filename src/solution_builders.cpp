#include "solution_builders.h"

SEXP Solution_builders::GCV_Newton_sol(const VectorXr & solution, const output_Data & output)
{
       //Copy result in R memory
       SEXP result = NILSXP;
       result = PROTECT(Rf_allocVector(VECSXP, 7));
       SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, solution.size()));
       Real *rans = REAL(VECTOR_ELT(result, 0));
       for(UInt j = 0; j < solution.size(); j++)  //[TO DO ] //sono le f_hat e g_hat, si potrebbe rimuovere, cambiando la chiamata da R in  smooth.FEM.basis
       {
               rans[j] = solution[j];
       }


       UInt size_z=output.z_hat.size();
       SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, size_z));
       rans = REAL(VECTOR_ELT(result, 1));

       for(UInt j = 0; j < size_z; j++)
       {
               rans[j] = output.z_hat[j];
       }

       //Rprintf("Hey doc,  %f %f %f\n", output.z_hat[0], output.z_hat[1], output.z_hat[3]);
       SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 2));
       rans[0] = output.SS_res;

       SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, 1));
       rans= REAL(VECTOR_ELT(result, 3));
       rans[0] = output.sigma_hat_sq;

       SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 4));
       rans[0] = output.lambda_sol;

       SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 5));
       rans[0] = output.n_it;

       SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 6));
       rans[0] = output.time_partial;


       UNPROTECT(1);

       return(result);
}

SEXP Solution_builders::GCV_batch_sol(const VectorXr & solution, const output_Data_opt & output, const timespec & partial_t)
{
       //Copy result in R memory
       //useless
       SEXP result = NILSXP;
       result = PROTECT(Rf_allocVector(VECSXP, 8));
       SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, solution.size()));
       Real *rans = REAL(VECTOR_ELT(result, 0));
       for(UInt j = 0; j < solution.size(); j++)  //[TO DO ] //sono le f_hat e g_hat, si potrebbe rimuovere, cambiando la chiamata da R in  smooth.FEM.basis
       {
               rans[j] = solution[j];
       }

       UInt size_z=output.z_hat_opt.size();
       SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, size_z));
       rans = REAL(VECTOR_ELT(result, 1));

       for(UInt j = 0; j < size_z; j++)
       {
               rans[j] = output.z_hat_opt[j];
       }

       //Rprintf("Hey doc,  %f %f %f\n", output.z_hat[0], output.z_hat[1], output.z_hat[3]);
       SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 2));
       rans[0] = output.SS_res_opt;

       SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, 1));
       rans= REAL(VECTOR_ELT(result, 3));
       rans[0] = output.sigma_hat_sq_opt;

       SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 4));
       rans[0] = output.lambda_opt;

       UInt size_vec=output.GCV_evals.size();
       SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, size_vec));
       rans = REAL(VECTOR_ELT(result, 5));
       for(UInt j = 0; j < size_vec; j++)
       {
               rans[j] = output.GCV_evals[j];
       }


       SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 6));
       rans[0] = output.GCV_opt;


       SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, 1));
       rans = REAL(VECTOR_ELT(result, 7));
       rans[0] = partial_t.tv_sec + 1e-9*partial_t.tv_nsec;


       UNPROTECT(1);

       return(result);
}
