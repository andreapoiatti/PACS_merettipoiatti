#ifndef NEWTON_H
#define NEWTON_H

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "fdaPDE.h"
#include "solver.h"
#include "function_variadic.h"

//il .h diventerà OPT methods
//I metodi saranno tutti in questa classe!!

template <typename Tuple, typename Hessian, typename... Extensions> //o altri che si scelgono, come Ctype o Dtype
class Opt_methods
{
        protected:
                Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F;
                //virtual members
                //constructor prende in ingresso Function F e eventuali altri oggetti

                // Contructor
                Opt_methods(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): F(F_) {}
};


// Classes
template <typename Tuple>
struct Auxiliary
{
};

template<>
struct Auxiliary<Real>
{
        public:
                Auxiliary(void) {};

                static inline bool isNull(Real n)                       {return (n == 0);}
                static inline void divide(Real a, Real b, Real & x)     {x = b/a;}
                static inline Real error(Real a)                        {return std::abs(a);}
};

template<>
struct Auxiliary<VectorXr>
{
        public:
                Auxiliary(void) {};

                static inline bool isNull(MatrixXr n)
                {
                        UInt sz = n.size();
                        return (n == MatrixXr::Zero(sz,sz));
                }

                static inline void divide(const MatrixXr & A, const VectorXr & b, VectorXr & x)
                {
                        Cholesky::solve(A, b, x);
                }

                static inline Real error(VectorXr a)
                {
                        return a.norm();
                }
};

// i metodi per le funzioni in F si chiamano evaluate_f, evaluate_first_derivative, evaluate_second_derivative
template <typename Tuple, typename Hessian, typename ...Extensions>
class Newton: public Opt_methods<Tuple, Hessian, Extensions...>
{
        //eventuali metodi saranno overridden da opt methods
        public:
                Newton(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): Opt_methods<Tuple, Hessian, Extensions...>(F_) {}; //esempio di possibile constructor
                // non può prendere in ingresso const ref, deve modificare l'oggetto F

                std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch)
                {
                        // Initialize the algorithm
                        Tuple x_old;
                        Tuple x      = x0;
                        UInt  n_iter = 0;
                        Real  error  = std::numeric_limits<Real>::infinity();

                        while (n_iter < max_iter)
                        {
                                //Debugging purpose f(x)
                                Real    fp  = this->F.evaluate_f(x);
                                Tuple   fpx = this->F.evaluate_first_derivative (x);
                                Hessian fsx = this->F.evaluate_second_derivative(x);

                                if (Auxiliary<Tuple>::isNull(fsx))
                                {
                                        // Debug message
                                        // std::cout << "Division by zero" << std::endl;
                                        return {x, n_iter};
                                }

                                ++n_iter;
                                x_old = x;
                                Auxiliary<Tuple>::divide(fsx, fpx, x);
                                x = x_old - x;

                                error = Auxiliary<Tuple>::error(fpx);
                                Rprintf("error:%f\n", error);
                                if (error < tolerance)
                                {
                                        ch.set_tolerance();
                                        return {x, n_iter};
                                }
                        }

                        ch.set_max_iter();
                        return {x, n_iter};
                }
};

#endif
