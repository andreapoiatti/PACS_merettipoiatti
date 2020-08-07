#ifndef __OPTIMIZATON_METHODS_FACTORY_H__
#define __OPTIMIZATON_METHODS_FACTORY_H__

#include "Newton.h"
#include "../../Global_Utilities/Headers/Make_Unique.h"
#include <memory>

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda for each PC.
template<typename Function, typename Tuple, typename Hessian, typename EvaluationType>
class Opt_method_factory
{
	public:
        	//! A method that takes as parameter a string and builds a pointer to the right object
        	static std::unique_ptr<Opt_methods<Tuple,Hessian,EvaluationType>> create_Opt_method(const std::string & validation, Function & F)
                {
                	if(validation=="newton")
                                return make_unique<Newton_ex<Real, Real, EvaluationType>>(F);
                	if(validation=="newton_fd")
                                return make_unique<Newton_fd<Real, Real, EvaluationType>>(F);
        	}
};

#endif
