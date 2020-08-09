#include "../Headers/Optimization_Data.h"

void OptimizationData::builder_utility(SEXP Roptim, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune)
{
        UInt criterion = INTEGER(Roptim)[0];
        if(criterion == 2)
        {
                this->set_criterion("newton_fd");
        }
        else if(criterion == 1)
        {
                this->set_criterion("newton");
        }
        else if(criterion == 0)
                this->set_criterion("batch");

        UInt DOF_evaluation = INTEGER(Roptim)[1];
        if(DOF_evaluation == 0)
        {
                this->set_DOF_evaluation("not_required");
        }
        else if(DOF_evaluation == 1)
        {
                this->set_DOF_evaluation("stochastic");
                this->set_nrealizations(INTEGER(Rnrealizations)[0]);
                this->set_seed(INTEGER(Rseed)[0]);
        }
        else
        {
                this->set_DOF_evaluation("exact");
        }

        UInt loss_function = INTEGER(Roptim)[2];
        if(loss_function == 0)
        {
                this->set_loss_function("unused");
        }
        else if(loss_function == 1)
        {
                this->set_loss_function("GCV");
        }

        // Tuning parameter
        this->set_tuning(REAL(Rtune)[0]);

        // DOF_matrix
        UInt n_ = INTEGER(Rf_getAttrib(RDOF_matrix, R_DimSymbol))[0];
        UInt p_ = INTEGER(Rf_getAttrib(RDOF_matrix, R_DimSymbol))[1];
        DOF_matrix.resize(n_, p_);
        for(auto i=0; i<n_; ++i)
        {
                for(auto j=0; j<p_ ; ++j)
                {
                        DOF_matrix(i,j) = REAL(RDOF_matrix)[i+ n_*j];
                }
        }

        if(n_==0 || p_==0)
        {
                n_ = 0;
                p_ = 0;
        }
        DOF_matrix.resize(n_,p_);
}
void OptimizationData::fill_lambda(SEXP Rlambda, std::vector<Real> & vect, UInt & size)
{
        size = Rf_length(Rlambda);
        vect.resize(size);

        for(UInt i=0; i<size; ++i)
        {
                vect[i] = REAL(Rlambda)[i];
        }
}

void OptimizationData::initialize_lambda(SEXP Rlambda, Real & init)
{
        if(Rf_length(Rlambda)>0)
                init = REAL(Rlambda)[0];
}

OptimizationData::OptimizationData(SEXP Roptim, SEXP Rlambda, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune)
{
        builder_utility(Roptim, Rnrealizations, Rseed, RDOF_matrix, Rtune);

        // Lambda
        if(this->criterion == "batch")
        {
                fill_lambda(Rlambda, this->lambda_S, this->size_S);
                set_lambdaS_backup();
        }
        else
        {
                initialize_lambda(Rlambda, this->initial_lambda_S);
        }
}

OptimizationData::OptimizationData(SEXP Roptim, SEXP Rlambda_S, SEXP Rlambda_T, SEXP Rflag_parabolic, SEXP Rnrealizations, SEXP Rseed, SEXP RDOF_matrix, SEXP Rtune)
{
        builder_utility(Roptim, Rnrealizations, Rseed, RDOF_matrix, Rtune);

        // Lambda
        if(this->criterion == "batch")
        {
                fill_lambda(Rlambda_S, this->lambda_S, this->size_S);
                if(INTEGER(Rflag_parabolic)[0] == 0)
                        fill_lambda(Rlambda_T, this->lambda_T, this->size_T);
                set_lambdaS_backup();
        }
        else
        {
                initialize_lambda(Rlambda_S, this->initial_lambda_S);
                if(INTEGER(Rflag_parabolic)[0] == 0)
                        initialize_lambda(Rlambda_T, this->initial_lambda_T);
        }
}

void OptimizationData::print_opt_data(void) const
{
        std::cout << "Optimization data:\n";
        std::cout << "Criterion: " << criterion << "\n";
        std::cout << "DOF valuation: " << DOF_evaluation << "\n";
        std::cout << "Loss Function: " << loss_function << "\n";
}
