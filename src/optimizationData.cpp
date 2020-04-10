#include "optimizationData.h"

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
