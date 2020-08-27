#ifndef __NEWTON_H__
#define __NEWTON_H__

// HEADERS
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include "../../FdaPDE.h"
#include "../../FE_Assemblers_Solvers/Include/Solver.h"
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

                inline int which(void) const   //!< Returns the reason of conclusion of the iterative method
                {
                        if (reached_tolerance == true)
                                return 1;
                        else if (reached_max_iter ==  true)

                                return 2;
                        else
                                return -1; //error status
                }
};

//! Father class to apply generic optimization methods
/*!
 * \tparam Tuple image type of the gradient of the function
 * \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam Extensions input class if the computations need members already stored in a class
 */
template <typename Tuple, typename Hessian, typename... Extensions>
class Opt_methods
{
        protected:
                //! Contructor
                Opt_methods(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): F(F_) {}

        public:
                Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F; //!< needed to be public, to be able to access to other methods of the class F from outside

                //! Function to apply the optimization method and obtain as a result the couple (optimal lambda, optimal value of the function)
                virtual std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<Real> & lambda_v) = 0;
};

// Classes
template <typename Tuple>
struct Auxiliary
{
        // NOT yet implemented
};

//! Auxiliary class to perform elementary mathematical operations and checks: specialization for 1 dimensional case
template<>
struct Auxiliary<Real>
{
        public:
                Auxiliary(void) {};

                static inline bool isNull(Real n)                       {return (n == 0);}      //!< Check if the input value is zero
                static inline void divide(Real a, Real b, Real & x)     {x = b/a;}              //!< Apply a division
                static inline Real residual(Real a)                     {return std::abs(a);}   //!< Compute the norm of the residual
};

//! Auxiliary class to perform elementary mathematical operations and checks: specialization for n dimensional case
template<>
struct Auxiliary<VectorXr>
{
        public:
                Auxiliary(void) {};

                //! Check if the input value is zero
                static inline bool isNull(MatrixXr n)
                {
                        UInt sz = n.size();
                        return (n == MatrixXr::Zero(sz,sz));
                }

                //! Solve a linear system in the optimization method
                static inline void divide(const MatrixXr & A, const VectorXr & b, VectorXr & x)
                {
                        Cholesky::solve(A, b, x);
                }

                //! Compute the norm of the residual
                static inline Real residual(VectorXr a)
                {
                        return a.norm();
                }
};


//! Class to apply Newton exact method, inheriting from Opt_methods
/*!
 * \tparam Tuple image type of the gradient of the function
 * \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 * \tparam Extensions input class if the computations need members already stored in a class
 */
 template <typename Tuple, typename Hessian, typename ...Extensions>
 class Newton_ex: public Opt_methods<Tuple, Hessian, Extensions...>
 {
         public:
                 //! Constructor
                 /*!
                  \note F cannot be const, it must be modified
                 */
                 Newton_ex(Function_Wrapper<Tuple, Real, Tuple, Hessian, Extensions...> & F_): Opt_methods<Tuple, Hessian, Extensions...>(F_)
                 {
                         // Debugging purpose
                         // Rprintf("Newton method built\n");
                 };

                 //! Apply Newton's method
                 std::pair<Tuple, UInt> compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<Real> & lambda_v) override;
};


template <typename Tuple, typename Hessian, typename ...Extensions>
class Newton_fd: public Opt_methods<Tuple, Hessian, Extensions...>
{
        //NOT yet implemented
};

//! Class to apply Newton method exploting finite differences to compute derivatives, inheriting from Opt_methods
/*!
 \tparam Tuple image type of the gradient of the function
 \tparam Hessian image type of the Hessian of the function: if the dimension of the image is >1 (and domain >1), problems to store the hessian, it's a tensor
 \tparam Extensions input class if the computations need members already stored in a class
*/
template <typename ...Extensions>
class Newton_fd<Real, Real, Extensions...>: public Opt_methods<Real, Real, Extensions...>
{
        public:
                //! Constructor
                /*!
                 \note F cannot be const, it must be modified
                */
                Newton_fd(Function_Wrapper<Real, Real, Real, Real, Extensions...> & F_): Opt_methods<Real, Real, Extensions...>(F_) {};

                //! Apply Newton fd method
                std::pair<Real, UInt> compute (const Real & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<Real> & lambda_v) override;
};


#include "Newton_imp.h"

#endif
