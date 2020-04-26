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
                virtual std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch) = 0;
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

Real bisection(const Real &aa,const  Real &bb, const UInt& max_it, const Real &eval_a_, const Real &eval_b_)
{    Real a=aa, b=bb;
     Real eval_a=eval_a_;
    if( eval_a*eval_b_>= 0)
    {
        Rprintf("\n Incorrect a and b\n");
        return 1;
    }


    Real c=a;
    UInt n_it=0;
    Real eval_c;

    while (n_it<max_it) //tell about interval is sufficently small
    {
        c = std::sqrt(a*b);

        eval_c=this->F.evaluate_second_derivative(c);
        if ( eval_c== 0.0){

            break;
        }
        else if (eval_c*eval_a < 0){

                b=c;
        }
        else{

                a=c;
                eval_a=eval_c;
        }
        n_it++;
       
    }
  return c;
}


                std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch) override
                {
                        // Initialize the algorithm
                        Tuple x_old;
                        Tuple x      = x0;
                        UInt  n_iter = 0;
                        Real  error  = std::numeric_limits<Real>::infinity();



                        Real a=1e-7;
			Real b=2;
			Real eval_a=this->F.evaluate_second_derivative(a); //useful to save one evaluation in the bisection method if the extrema are already correct
			Real eval_b=this->F.evaluate_second_derivative(b);
                        while (eval_a<=0)
                                  {a*=10;
				   eval_a=this->F.evaluate_second_derivative(a);
                                   }
			while (eval_b>=0)
                                  {b*=10;
				   eval_b=this->F.evaluate_second_derivative(b);}

                        Rprintf("\n Starting interval for preprocessing: [%f;%f]\n", a,b);

                        Real flesso=bisection(a,b,4, eval_a, eval_b);
                        Rprintf("\nFlesso at %f\n", flesso);
                        
			if (x>flesso/5 || x<a)
                              {x=flesso/50;
				Rprintf("\nInitial value inserted is out of range, using default value lambda=%f\n",x);
				}	
                        Rprintf("\n Starting Newton's iteratios: starting point lambda=%f\n",x);

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


                Real bisection(const Real &aa,const  Real &bb, const UInt& max_it, const Real &eval_a_, const Real &eval_b_, const Real& h)
		{    Real a=aa, b=bb;
     		     Real eval_a=eval_a_;
    		if( eval_a*eval_b_>= 0)
                    {
                        Rprintf("\n Incorrect a and b\n");
                        return 1;
                    }


                    Real c=a;
                    UInt n_it=0;
                    Real eval_c;
		  
                    while (n_it<max_it) //tell about interval is sufficently small
                    {
                        c = std::sqrt(a*b);

                        eval_c=this->second_derivative(c,h);
                        if ( eval_c== 0.0){

                            break;
                        }
                        else if (eval_c*eval_a < 0){

                                b=c;
                        }
                        else{

                                a=c;
                                
               		        eval_a=eval_c;
                        }
                        n_it++;
                      
                    }
                  return c;
                }


                Real second_derivative(const Real& x, const Real& h)
               {
                       Rprintf("Forward: \n");
                       Real fxph = this->F.evaluate_f(x+h);
                       Rprintf("Backward: \n");
                       Real fxmh = this->F.evaluate_f(x-h);
                       Rprintf("Center: \n");
                       Real fx  = this->F.evaluate_f(x);
                       return (fxph+fxmh-(2*fx))/(h*h);

                 }

                std::pair<Real, UInt> compute (const Real & x0, const Real tolerance, const UInt max_iter, Checker & ch) override
                {
                        // Initialize the algorithm
                        Real x_old;
                        Real x      = x0;
                        UInt  n_iter = 0;
                        Real  error  = std::numeric_limits<Real>::infinity();
		        Real a=1e-7;
			Real b=2;


			Real h = 1e-5;                        

			Real eval_a=this->second_derivative(a,h); //useful to save one evaluation in the bisection method if the extrema are already correct
			Real eval_b=this->second_derivative(b,h);
                        while (eval_a<=0)
                                  {a*=10;
				   eval_a=this->second_derivative(a,h);
                                   }
			while (eval_b>=0)
                                  {b*=10;
				   eval_b=this->second_derivative(b,h);}

                        Rprintf("\n Starting interval for preprocessing: [%f;%f]\n", a,b);

                        Real flesso=bisection(a,b,4, eval_a, eval_b,h);
 			Rprintf("\nFlesso at %f\n", flesso);
                        
                      if (x>flesso/5 || x<a)
                              {x=flesso/50;
				Rprintf("\nInitial value inserted is out of range, using default value lambda=%f\n",x);
				}


			Rprintf("\n Starting Newton's iteratios: starting point lambda=%f\n",x);

                        

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
