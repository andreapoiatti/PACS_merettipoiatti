#ifndef FUNCTION_H
#define FUNCTION_H

// Headers
#include <functional>
#include <type_traits>
#include <iostream>
#include "lambda_optimizer.h"
//#include "fdaPDE.h"

// Classes

//! Function wrapper. Stores a function from R^m to R^n. Variadic template, to store
/*!
 * \param       Dtype          domain type of the function
 * \param       Ctype          image type of the function
 * \param       Tuple          image type of the gradient of the function
 * \param       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \param       Extensions    input class if th e computations need members already stored in a class
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
                  std::function<Ctype(Dtype)> g;
                  std::function<Tuple(Dtype)> dg; //derivative
                  std::function<Hessian(Dtype)> ddg; //second_derivative
        public:

                //! Constructor taking object Extensions and initializing with it the new class
                   template<typename D=typename std::enable_if<sizeof...(Extensions)!=0,void>,typename ...T>  //o risulta true il primo argomento e definisce typedef, oppure non lo deinisce proprio
                   Function_Wrapper(T&& ...ext):Extensions(std::forward<T>(ext))...{};//{std::cout<<"end"<<std::endl;};
                   //! Since I have defined a constructor I need to indicate the default constructor
                   Function_Wrapper()=default;


                 //! Function to inizialize the g,dg,ddg, used only if you explicit directly the functions
                 inline void set_functions(const std::function<Ctype(Dtype)> & g_,const std::function<Tuple(Dtype)> & dg_,const std::function<Hessian(Dtype)> & ddg_)
                 {
                   g=g_;
                   dg=dg_;
                   ddg=ddg_;
                  };

                  // Operator
                  //! Function version of a std::function () operator, use SFINAE to differentiate from the two possibilities of instantiation
                  template <typename U>
                  typename std::enable_if<sizeof...(Extensions)==0 || std::is_enum<U>::value, Ctype>::type //is_enum<U> is always false for us, used to make the deduction argument and SFINAE work
                   evaluate_f  ( U  lambda) {
                                              return g(lambda);
                                             };

                 template <typename U>
                 typename std::enable_if<sizeof...(Extensions)!=0 || std::is_enum<U>::value, Ctype>::type
                  evaluate_f  ( U  lambda ) {
                                              return this->compute_f(lambda);
                                             };


                  //! Evaluation of first derivative
                  template <typename U>
                  typename std::enable_if<sizeof...(Extensions)==0 || std::is_enum<U>::value, Tuple>::type //is_enum<U> is always false for us, used to make the deduction argument and SFINAE work
                   evaluate_first_derivative  ( U  lambda) {
                                            return dg(lambda);
                                           };

                 template <typename U>
                 typename std::enable_if<sizeof...(Extensions)!=0 || std::is_enum<U>::value, Tuple>::type
                  evaluate_first_derivative  ( U  lambda ) {
                                                         return this->compute_fp(lambda);
                                                            };

                  //! Evaluation of second derivative

                  template <typename U>
                  typename std::enable_if<sizeof...(Extensions)==0 || std::is_enum<U>::value, Hessian>::type //is_enum<U> is always false for us, used to make the deduction argument and SFINAE work
                   evaluate_second_derivative  ( U  lambda) {
                                            return ddg(lambda);
                                           };

                 template <typename U>
                 typename std::enable_if<sizeof...(Extensions)!=0 || std::is_enum<U>::value, Hessian>::type
                  evaluate_second_derivative  ( U  lambda ) {
                                                         return this->compute_fs(lambda);
                                                             };

 };

//usata come prova, sarà la GCV ecc
class Prova
{       int v=0;

        public:
        // Operator
        //! Function version of a std::function () operator
        Prova()=default;
        Prova(int b): v(b){};
        inline double compute_f  ( double  lambda ) {return 4+v;};


        //! Evaluation of first derivative
        inline double compute_fp ( double lambda ) {return 34-v;};

        //! Evaluation of second derivative
        inline double compute_fs  (double  lambda ) {return 42+v; };

};


//! Checker, contains data regarding the last process
class Checker
{
        private:
                bool reached_max_iter;
                bool reached_tolerance;

        public:
                Checker(void): reached_max_iter(false), reached_tolerance(false) {}

                inline void set_max_iter(void)  {reached_max_iter  = true;}
                inline void set_tolerance(void) {reached_tolerance = true;}

                inline UInt which(void) const
                {
                        if (reached_tolerance == true)
                                return 1;
                        else if (reached_max_iter ==  true)
                                return 2;
                        else
                                return -1;
                }
};



#endif
