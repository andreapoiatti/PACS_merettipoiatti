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
                //virtual members
                //constructor prende in ingresso Function F e eventuali altri oggetti

                // Contructor
                Opt_methods(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): F(F_) {}
        public:
                Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F; //needed to be public, to be able to access to other methods of the class F from outside
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
                        return (n == MatrixXr::Zero(sz,sz));        // DEBUGGING PURPOSE
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
class Newton_ex: public Opt_methods<Tuple, Hessian, Extensions...>        // DEBUGGING PURPOSE
{
        //eventuali metodi saranno overridden da opt methods
        public:
                Newton_ex(Function_Wrapper<Tuple, Real, Tuple, Real, Extensions...> & F_): Opt_methods<Tuple, Hessian, Extensions...>(F_) {Rprintf("Newton method built\n");}; //esempio di possibile constructor
                // non può prendere in ingresso const ref, deve modificare l'oggetto F

                std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch)
                {
                        // Initialize the algorithm
                        Tuple x_old;
                        Tuple x      = x0;
                        UInt  n_iter = 0;
                        Real  error  = std::numeric_limits<Real>::infinity();

                        //only the first time applied here
                        Real   fx  = this->F.evaluate_f(x);
                        Tuple   fpx = this->F.evaluate_first_derivative (x);
                        Hessian fsx = this->F.evaluate_second_derivative(x);

                        while (n_iter < max_iter)
                        {
                                //Debugging purpose f(x)

                                if (Auxiliary<Tuple>::isNull(fsx))
                                {
                                        // Debug message
                                        // std::cout << "Division by zero" << std::endl;
                                        return {x, n_iter};
                                }

                                ++n_iter;

                                Rprintf("\nStep number %d  of EXACT-NEWTON\n", n_iter);
                                x_old = x;
                                Auxiliary<Tuple>::divide(fsx, fpx, x);
                                x = x_old - x;

                                //put here the updates in order to compute error on the correct derivative and to have z_hat updated for the solution
                                fx = this->F.evaluate_f(x);
                                fpx = this->F.evaluate_first_derivative (x);

                                fsx = this->F.evaluate_second_derivative(x);

                                error = Auxiliary<Tuple>::error(fpx);
                                Rprintf("error: %f\n", error);
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

template <typename Tuple, typename Hessian, typename ...Extensions>
class Newton_fd: public Opt_methods<Tuple, Hessian, Extensions...>
{

};

template <typename ...Extensions>
class Newton_fd<Real, Real, Extensions...>: public Opt_methods<Real, Real, Extensions...>
{
        //eventuali metodi saranno overridden da opt methods
        public:
                Newton_fd(Function_Wrapper<Real, Real, Real, Real, Extensions...> & F_): Opt_methods<Real, Real, Extensions...>(F_) {}; //esempio di possibile constructor
                // non può prendere in ingresso const ref, deve modificare l'oggetto F

                std::pair<Real, UInt> compute (const Real & x0, const Real tolerance, const UInt max_iter, Checker & ch)
                {
                        // Initialize the algorithm
                        Real x_old;
                        Real x      = x0;
                        UInt  n_iter = 0;
                        Real  error  = std::numeric_limits<Real>::infinity();

                        Real h = 1e-4;

                        //only the first time applied here
                        Rprintf("Forward: \n");
                        Real fxph = this->F.evaluate_f(x+h);
                        Rprintf("Backward: \n");
                        Real fxmh = this->F.evaluate_f(x-h);
                        Rprintf("Center: \n");
                        Real fx  = this->F.evaluate_f(x);

                        Real fpx = (fxph-fxmh)/(2*h);
                        Rprintf("fp(x): %f\n", fpx);

                        Real fsx = (fxph+fxmh-(2*fx))/(h*h);
                        Rprintf("fs(x): %f\n", fsx);

                        while (n_iter < max_iter)
                        {
                                //Debugging purpose f(x)


                                if (Auxiliary<Real>::isNull(fsx))
                                {
                                        // Debug message
                                        // std::cout << "Division by zero" << std::endl;
                                        return {x, n_iter};
                                }

                                ++n_iter;

                                Rprintf("\nStep number %d  of FD-NEWTON\n", n_iter);
                                x_old = x;
                                Auxiliary<Real>::divide(fsx, fpx, x);
                                x = x_old - x;

                                //put here the updates in order to compute error on the correct derivative and to have z_hat updated for the solution
                                Rprintf("Forward:\n");
                                fxph = this->F.evaluate_f(x+h);
                                Rprintf("Backward: \n");
                                fxmh = this->F.evaluate_f(x-h);
                                Rprintf("Center: \n");
                                fx  = this->F.evaluate_f(x);

                                fpx = (fxph-fxmh)/(2*h);
                                Rprintf("fp(x): %f\n", fpx);

                                fsx = (fxph+fxmh-(2*fx))/(h*h);
                                Rprintf("fs(x): %f\n", fsx);

                                error = Auxiliary<Real>::error(fpx);
                                Rprintf("error: %f\n", error);
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
