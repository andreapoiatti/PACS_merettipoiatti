#ifndef __NEWTON_H__
#define __NEWTON_H__

// Headers
#include <cmath>
#include <limits>
#include <utility>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Headers/Solver.h"
#include "Function_Variadic.h"

//! Checker, contains data regarding the last process, used in optimizaion processes
class Checker
{
        private:
                bool reached_max_iter;
                bool reached_tolerance;

        public:
                Checker(void): reached_max_iter(false), reached_tolerance(false) {}

                inline void set_max_iter(void)  {reached_max_iter  = true;} //!< Sets max number of iterations
                inline void set_tolerance(void) {reached_tolerance = true;} //!< Sets the tolerance for the optimization method

                inline UInt which(void) const   //!<Returns the reason of conclusion of the iterative method
                {
                        if (reached_tolerance == true)
                                return 1;
                        else if (reached_max_iter ==  true)
                                return 2;
                        else
                                return -1;
                }
};

//! Father class to apply generic optimization methods
/*!
 * \tparam       Tuple          image type of the gradient of the function
 * \tparam       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam       Extensions    input class if the computations need members already stored in a class
 */
template <typename Tuple, typename Hessian, typename... Extensions>
class Opt_methods
{
        protected:
                //!virtual members

                //! Contructor
                Opt_methods(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): F(F_) {}
        public:
                Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F; /*! needed to be public, to be able to access to other methods of the class F from outside */
                virtual std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch) = 0; //!< Function to apply the optimization method and obtain as a result the couple (optimal lambda, optimal value of the function)
};

// Classes
template <typename Tuple>
struct Auxiliary
{
        //NOT yet implemented
};

//!< Auxiliary class to perform elementary mathematical operations and checks: specialization for 1 dimensional case
template<>
struct Auxiliary<Real>
{
        public:
                Auxiliary(void) {};

                static inline bool isNull(Real n)                       {return (n == 0);} //! Check if the input value is zero
                static inline void divide(Real a, Real b, Real & x)     {x = b/a;}  //! Apply a division
                static inline Real residual(Real a)                        {return std::abs(a);}  //! Compute the norm of the residual
};

template<>
struct Auxiliary<VectorXr>  //!< Auxiliary class to perform elementary mathematical operations and checks: specialization for n dimensional case
{
        public:
                Auxiliary(void) {};

                static inline bool isNull(MatrixXr n)  //! Check if the input value is zero
                {
                        UInt sz = n.size();
                        return (n == MatrixXr::Zero(sz,sz));
                }

                static inline void divide(const MatrixXr & A, const VectorXr & b, VectorXr & x)  //! Solve a linear system in the optimization method
                {
                        Cholesky::solve(A, b, x);
                }

                static inline Real residual(VectorXr a)    //! Compute the norm of the residual
                {
                        return a.norm();
                }
};


//! Class to apply Newton exact method, inheriting from Opt_methods
/*!
 * \tparam       Tuple          image type of the gradient of the function
 * \tparam       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam       Extensions    input class if the computations need members already stored in a class
 */
 template <typename Tuple, typename Hessian, typename ...Extensions>
 class Newton_ex: public Opt_methods<Tuple, Hessian, Extensions...>
 {
         public:
                 Newton_ex(Function_Wrapper<Tuple, Real, Tuple, Real, Extensions...> & F_): Opt_methods<Tuple, Hessian, Extensions...>(F_) {Rprintf("Newton method built\n");}; //!Constructor
                 /*! F cannot be const, it must be modified*/

                 std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch) override //!< Apply Newton's method
                 {
                        // Initialize the algorithm
                        Tuple x_old;
                        Tuple x      = x0;
                        UInt  n_iter = 0;
                        Real  error  = std::numeric_limits<Real>::infinity();

                        Rprintf("\n Starting Initializing lambda phase\n"); /*! Start from 6 lambda and find the minimum value of GCV to start from it the newton's method*/

                        Real valmin, valcur, lambda_min;
                        UInt Nm = 6;
                        std::vector<Real> vals={5.000000e-05, 1.442700e-03, 4.162766e-02, 1.201124e+00, 3.465724e+01, 1.000000e+03};
                        valcur = this->F.evaluate_f(vals[0]);
                        lambda_min = 5e-5;
                        valmin = valcur;

                        for(UInt i=1; i<Nm; i++)
                        {
                                valcur = this->F.evaluate_f(vals[i]);
                                if(valcur<valmin)
                                {
                                        valmin = valcur;
                                        lambda_min = vals[i];
                                }
                        }

                        if(x>lambda_min/10)
                        {
                                x = lambda_min/10;
                        }

                        Rprintf("\n Starting Newton's iterations: starting point lambda=%f\n",x);

                        //only the first time applied here
                        Real    fx  = this->F.evaluate_f(x);
                        Tuple   fpx = this->F.evaluate_first_derivative (x);
                        Hessian fsx = this->F.evaluate_second_derivative(x);

                        while(n_iter < max_iter)
                        {
                                //Debugging purpose f(x)

                                if(Auxiliary<Tuple>::isNull(fsx))
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
                                //if (x<1e-8) x=(1./(2*n_iter))*flesso/50;
                                if (x<1e-8)
                                {
                                        x = (1./(2*n_iter))*lambda_min/10;
                                }

                                //put here the updates in order to compute error on the correct derivative and to have z_hat updated for the solution
                                fx  = this->F.evaluate_f(x);
                                fpx = this->F.evaluate_first_derivative (x);
                                fsx = this->F.evaluate_second_derivative(x);

                                error = Auxiliary<Tuple>::residual(fpx);
                                Rprintf("Residual: %f\n", error);

                                if(error<tolerance)
                                {
                                        /* fpx=this->F.evaluate_f(x-0.5);
                                        fsx=this->F.evaluate_f(x+0.5);
                                        if (std::abs(fpx-fsx)/fsx<=0.01)
                                                Rprintf("GCV has a non standard shape");*/

                                        ch.set_tolerance();
                                        fx = this->F.evaluate_f(x);
                                        return {x, n_iter};
                                }
                        }

                        fx = this->F.evaluate_f(x);
                        ch.set_max_iter();
                        return {x, n_iter};
                }
};


template <typename Tuple, typename Hessian, typename ...Extensions>
class Newton_fd: public Opt_methods<Tuple, Hessian, Extensions...>
{
        //NOT yet implemented
};

//! Class to apply Newton method exploting finite differences to compute derivatives, inheriting from Opt_methods
/*!
 * \tparam       Tuple          image type of the gradient of the function
 * \tparam       Hessian        image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam       Extensions    input class if the computations need members already stored in a class
 */
template <typename ...Extensions>
class Newton_fd<Real, Real, Extensions...>: public Opt_methods<Real, Real, Extensions...>
{
        public:
                Newton_fd(Function_Wrapper<Real, Real, Real, Real, Extensions...> & F_): Opt_methods<Real, Real, Extensions...>(F_) {}; //! Constructor
                // NB F cannot be const

                std::pair<Real, UInt> compute (const Real & x0, const Real tolerance, const UInt max_iter, Checker & ch) override
                {
                        // Initialize the algorithm
                        Real x_old;
                        Real x      = x0;
                        UInt n_iter = 0;
                        Real error  = std::numeric_limits<Real>::infinity();
			Real h      = 4e-6;

                        Rprintf("\n Starting Initializing lambda phase"); /*! Start from 6 lambda and find the minimum value of GCV to start from it the newton's method*/

                        Real valmin, valcur, lambda_min;
                        UInt Nm = 6;
			std::vector<Real> vals = {5.000000e-05, 1.442700e-03, 4.162766e-02, 1.201124e+00, 3.465724e+01, 1.000000e+03};
			valcur=this->F.evaluate_f(vals[0]);
			lambda_min = 5e-5;
			valmin = valcur;

                        for(UInt i=1; i<Nm; i++)
                        {
                                valcur = this->F.evaluate_f(vals[i]);
                                if(valcur<valmin)
                                {
                                        valmin = valcur;
                                        lambda_min = vals[i];
                                }
                        }

			if (x>lambda_min/10)
                        {
                                x = lambda_min/10;
                        }
			Rprintf("\n Starting Newton's iterations: starting point lambda=%f\n",x);

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

                        while(n_iter < max_iter)
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

                                if (x<1e-7) {x=(1./(2*n_iter))*lambda_min/10;

           			if (x<4e-6)
                                        x=5e-6;}
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

                                error = Auxiliary<Real>::residual(fpx);
                                Rprintf("residual: %f\n", error);
                                if (error < tolerance)
                                {       /*fpx=this->F.evaluate_f(x-0.5);
                                        fsx=this->F.evaluate_f(x+0.5);
                                        if (std::abs(fpx-fsx)/fsx<=0.01)
                                           Rprintf("GCV has a non standard shape");*/
                                        ch.set_tolerance();
                                        fx  = this->F.evaluate_f(x); //eventuale miglioramento: va fatto altirmenti prende gli z:hat di quellos sbagliato.
                                        return {x, n_iter};

                                }
                        }
                        fx  = this->F.evaluate_f(x);
                        ch.set_max_iter();
                        return {x, n_iter};
                }
};

#endif
