#ifndef __OPTMETHODSFACTORY_HPP__
#define __OPTMETHODSFACTORY_HPP__


#include "newton.h"


#include <memory>


template<typename T, typename... Args>  //per avere possibilità di scegliere più argomenti
std::unique_ptr<T> make_unique_opt(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda for each PC.
template<typename Function, typename Tuple, typename Hessian, typename EvaluationType>
class Opt_method_factory
{
	public:
	//! A method that takes as parameter a string and builds a pointer to the right object
	static std::unique_ptr<Opt_methods<Tuple,Hessian,EvaluationType>> create_Opt_method(const std::string & validation, Function & F){
	if(validation=="newton") return make_unique_opt<Newton_ex<Real, Real, EvaluationType>>(F);
	if(validation=="newton_fd") return make_unique_opt<Newton_fd<Real, Real, EvaluationType>>(F);
        //altri tipi
	}

};

#endif