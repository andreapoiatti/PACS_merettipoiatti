#ifndef __GOF_UPDATER_HPP__
#define __GOF_UPDATER_HPP__

#include<functional>
#include<vector>

template <typename LambdaOptim, typename T>
class GOF_updater
{
        private:
                std::vector<T> last_lambda_derivatives;
                std::vector<std::function<void(Real)>> updaters;
                LambdaOptim * start_ptr = nullptr;


                inline void call_from_to(UInt start, UInt finish, T lambda)
                {
                        for(UInt i=start; i<=finish; ++i)
                        {
                                updaters[i](lambda);
                                last_lambda_derivatives[i] = lambda;
                        }
                }

                inline void updaters_setter(LambdaOptim * lopt_ptr)
                {
                        this->updaters.reserve(3);
                        this->updaters.push_back(std::bind(&LambdaOptim::zero_updater, lopt_ptr, std::placeholders::_1));
                        this->updaters.push_back(std::bind(&LambdaOptim::first_updater, lopt_ptr, std::placeholders::_1));
                        this->updaters.push_back(std::bind(&LambdaOptim::second_updater, lopt_ptr, std::placeholders::_1));
                }

        public:
                GOF_updater(void) = default;

                inline void initialize(const std::vector<T> & first_lambdas)
                {
                        last_lambda_derivatives = first_lambdas;
                }

                inline void call_to(UInt finish, T lambda, LambdaOptim * lopt_ptr)
                {
                        if(start_ptr != lopt_ptr)
                        {
                                Rprintf("--- Set updaters ---\n");
                                initialize(std::vector<Real>{-1.,-1.,-1.});
                                updaters_setter(lopt_ptr);
                                start_ptr = lopt_ptr;
                        }

                        bool found = false;
                        for(UInt i = 0; i<=finish && found==false; ++i)
                                if(lambda != last_lambda_derivatives[i])
                                {
                                        call_from_to(i, finish, lambda);
                                        found = true;
                                }
                }
};

#endif
