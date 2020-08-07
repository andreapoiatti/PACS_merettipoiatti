#ifndef __FUNCTION_VARIADIC_H__
#define __FUNCTION_VARIADIC_H__

// Headers
#include <functional>
#include <iostream>
#include <type_traits>
#include "../../FdaPDE.h"

// Classes

//! Function wrapper. Stores a function from R^m to R^n. Variadic template to store either from inheritance or directly
/*!
 * \tparam       Dtype          domain type of the function
 * \tparam       Ctype          image type of the function
 * \tparam       Tuple          image type of the gradient of the function
 * \tparam       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam       Extensions    input class if the computations need members already stored in a class
 */

 template<typename Dtype, typename Ctype, typename Tuple, typename Hessian, typename... Extensions>
 class Function_Wrapper : public Extensions...
 {
        private:
                 //! Actual functions (function, first derivative/gradient, second derivative/Hessian) stored inside the Function class
                 /*!
                  * \param       _g      the std:::function from which to copy
                  * \param       _dg     the derivative std:::function from which to copy
                  * \param       _ddg     the second derivative std:::function from which to copy

                  */
                  std::function<Ctype(Dtype)>   g;
                  std::function<Tuple(Dtype)>   dg; //derivative
                  std::function<Hessian(Dtype)> ddg; //second_derivative
        public:
                //  -- CONSTRUCTORS --
                //! Constructor taking object Extensions and initializing with it the new class
                template<typename D=typename std::enable_if<sizeof...(Extensions)!=0,void>, typename ...T>
                Function_Wrapper(T&& ...ext):Extensions(std::forward<T>(ext))...{};

                //! Since I have defined a constructor I need to indicate the default constructor
                Function_Wrapper() = default;

                // --SETTERS --
                //! Function to inizialize the g,dg,ddg, used only if you explicit directly the functions
                inline void set_functions(const std::function<Ctype(Dtype)> & g_,const std::function<Tuple(Dtype)> & dg_,const std::function<Hessian(Dtype)> & ddg_)
                {
                        g   = g_;
                        dg  = dg_;
                        ddg = ddg_;
                }

                // -- OPEATORS --
                //! Function version of a std::function () operator, use SFINAE to differentiate from the two possibilities of instantiation (derived or not)
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)==0 || std::is_enum<U>::value, Ctype>::type //is_enum<U> is always false for us, used to make the deduction argument and SFINAE work
                evaluate_f(U lambda)
                {
                        return g(lambda);
                }

                //! Function version of a std::function () operator, use SFINAE to differentiate from the two possibilities of instantiation (derived or not)
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)!=0 || std::is_enum<U>::value, Ctype>::type
                evaluate_f(U lambda )
                {
                        return this->compute_f(lambda);
                }

                //! Evaluation of first derivative, if derived
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)==0 || std::is_enum<U>::value, Tuple>::type //is_enum<U> is always false for us, used to make the deduction argument and SFINAE work
                evaluate_first_derivative(U lambda)
                {
                        return dg(lambda);
                }

                //! Evaluatio of first derivative, if not derived
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)!=0 || std::is_enum<U>::value, Tuple>::type
                evaluate_first_derivative(U lambda)
                {
                        return this->compute_fp(lambda);
                }

                //! Evaluation of second derivative, if derived
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)==0 || std::is_enum<U>::value, Hessian>::type //is_enum<U> is always false for us, used to make the deduction argument and SFINAE work
                evaluate_second_derivative(U lambda)
                {
                        return ddg(lambda);
                }

                //! Evaluation of second derivative, if not derived
                template <typename U>
                typename std::enable_if<sizeof...(Extensions)!=0 || std::is_enum<U>::value, Hessian>::type
                evaluate_second_derivative(U lambda )
                {
                        return this->compute_fs(lambda);
                }
 };

//! DEBUGGING PURPORSE: Dummy class for debugging purpose
/*class Prova
{       int v=0;

        public:
        // Operator
        //! Function version of a std::function () operator
        Prova() = default;

        Prova(int b): v(b){};

        inline double compute_f  ( double  lambda ) {return 4+v;};


        //! Evaluation of first derivative
        inline double compute_fp ( double lambda ) {return 34-v;};

        //! Evaluation of second derivative
        inline double compute_fs  (double  lambda ) {return 42+v; };
};
*/


#endif
