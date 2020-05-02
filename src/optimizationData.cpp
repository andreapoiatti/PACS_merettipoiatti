#include "optimizationData.h"

OptimizationData::OptimizationData(SEXP ROPTmethod, SEXP Rlambdas, SEXP Rinitial_lambda, SEXP Rnrealizations)
{
        UInt criterion = INTEGER(ROPTmethod)[0];
        if(criterion == 2)
        {
                this->set_criterion_("newton_fd");
        }
        else if(criterion == 1)
        {
                this->set_criterion_("newton");
        }
        else if(criterion == 0)
                this->set_criterion_("batch");

        UInt method = INTEGER(ROPTmethod)[1];
        if(method == 0)
        {
                this->set_method_("gcv");
        }

        UInt stochastic = INTEGER(ROPTmethod)[2];
        if(stochastic == 1)
        {
                this->set_evaluation_("stochastic");
                this->set_nrealizations_(INTEGER(Rnrealizations)[0]);
        }
        else
        {
                this->set_evaluation_("exact");
        }

        if(criterion == 0)
        {       std::cout<<"Valori"<<std::endl;
                UInt n_lambdas_ = Rf_length(Rlambdas);
                std::vector<Real> lambdas_;
                lambdas_.resize(n_lambdas_);

                for(UInt i=0; i<n_lambdas_; ++i)
                {
                        lambdas_[i] = REAL(Rlambdas)[i];
                }

                this->set_lambdas_(lambdas_);
        }
        else if(criterion == 1 || criterion==2)
        {
                UInt initialization = INTEGER(ROPTmethod)[3];
                if(initialization == 1)
                        this->set_initial_lambda_(REAL(Rinitial_lambda)[0]);
        }
}

void OptimizationData::print_opt_data(void) const
{
        std::cout << "Optimization data:\n";
        std::cout << "Criterion: " << criterion_ << "\n";
        std::cout << "Evaluation: " << evaluation_ << "\n";
        std::cout << "Method: " << method_ << "\n";

        std::cout << "Lambdas:\n";
        for(UInt i=0; i<lambdas_.size(); ++i)
                std::cout << lambdas_[i] << " ";
        std::cout << "\n";
        std::cout << "Initial lambda: " << initial_lambda_ << "\n";
        std::cout << "N. realizations: " << nrealizations_ << std::endl;
}
