#ifndef __OPTMETHODSFACTORY_HPP__
#define __OPTMETHODSFACTORY_HPP__


#include "Opt_methods" //per ora è Newton


#include <memory>


template<typename T, typename... Args>  //pr avere possibilità di scegliere più argomenti
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda for each PC.
template<typename Tuple, typename Hessian>
class Opt_method_factory
{
	public:
	//! A method that takes as parameter a string and builds a pointer to the right object for the cross-validation
	static std::unique_ptr<Opt_methods> create_Opt_method(const std::string &method, Function&F){
	if(validation=="Newton") return make_unique<Newton<Tuple,Hessian>>(F);
	if(validation=="Newton_mod") return make_unique<Newton_mod<Tuple,Hessian>>(F, //altro);
        //altri tipi
	}

};

#endif
