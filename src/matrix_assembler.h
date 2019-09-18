#ifndef MATRIX_ASSEMBLER_H_
#define MATRIX_ASSEMBLER_H_

// Headers
#include "fdaPDE.h"
#include "finite_element.h"
#include "mesh_objects.h"
#include "param_functors.h"

// Classes

// ****  STIFFNESS ****
// *** ISOTROPIC ***
//! A Stiff class: a class for the stiffness operator.
class Stiff
{
        public:
                // Operators
        	//! A definition of operator () taking three arguments.
                /*!
                * Evaluates the stiffness operator (i,j) of the current planar finite element.
                * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER, 2, 2>, current planar finite element
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
        	template<class Integrator, UInt ORDER>
        	inline Real operator() (FiniteElement<Integrator, ORDER, 2, 2> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
        	{
                        // Evaluates partial component (i,j) of stiffness isotropic matrix
        	   	Real s = 0;
        	   	for (UInt icoor = 0; icoor < 2; icoor++)
        	   	{
        	   		s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);
        	   	}
        	   	return s;
        	}

                //! A definition of operator () taking three arguments.
                /*!
                * Evaluates the stiffness operator (i,j) of the current superficial finite element.
                * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2, 3>, current finite element
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
        	template<class Integrator, UInt ORDER>
        	inline Real operator() (FiniteElement<Integrator, ORDER, 2, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
        	{
        		Real s = 0;

        		Eigen::Matrix<Real, 2, 1> grad_phi_i;
        		Eigen::Matrix<Real, 2, 1> grad_phi_j;

        		grad_phi_i(0) = currentfe_.phiDerMaster(i, 0, iq);
        		grad_phi_i(1) = currentfe_.phiDerMaster(i, 1, iq);
        		grad_phi_j(0) = currentfe_.phiDerMaster(j, 0, iq);
        		grad_phi_j(1) = currentfe_.phiDerMaster(j, 1, iq);

        		s = grad_phi_i.dot(currentfe_.metric() * grad_phi_j);

        	   	return s;
        	}

                //! A definition of operator () taking three arguments.
                /*!
                * Evaluates the stiffness operator (i,j) of the current superficial finite element.
                * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER, 3, 3>, current finite element
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
        	template<class Integrator, UInt ORDER>
        	inline Real operator() (FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
        	{
        	   	Real s = 0;
        	   	for (UInt icoor = 0; icoor < 3; icoor++)
        	   	{
        	   		s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);
        	   	}
        	   	return s;
        	}
};

// *** ANYSOTROPIC ***
template <class Type>
class StiffAnys
{
};

// ** 2D in 2D **
template <>
class StiffAnys<Eigen::Matrix<Real, 2, 2>>
{
        private:
                //! A reference to FiniteElement<Integrator>
                /*!
                * Stores a reference to the finite element where the stiffness operator is evaluated.
                */
                //FiniteElement<Integrator, ORDER> & currentfe_;

                const Eigen::Matrix<Real, 2, 2> & K_;

        public:
                // Constructors
                //! A constructor.
                /*!
                 \param fe is a reference to FiniteElement<Integrator>
                 */
                StiffAnys(const Eigen::Matrix<Real, 2, 2> & K): K_(K){};

                // Operators
                //! A definition of operator () taking two arguments.
                /*!
                * Evaluates the stiffness operator (i,j) of the current finite elemente.
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
                //[TODO]

                //! A definition of operator () taking four arguments.
                /*!
                * Evaluates the product of: the derivative of basis(i) with respect to coordinate ic1 and the derivative of basis(j) with respect
                * to coordinate ic2 ,on current finite elemente.
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
                * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
                * returns a double.
                */
                template<class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 2> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                        // Implemeemts <J_{T^-1}^t * grad(\hat{phi}_i(p_hat)), K * J_{T^-1}^t * grad(\hat{phi}_j(p_hat))>
                   	Real s = 0;
                   	for (UInt icoor = 0; icoor < 2; icoor++)
                   	{
                   		s += (currentfe_.invTrJPhiDerMaster(i, 0, iq) * K_(0, icoor) * currentfe_.invTrJPhiDerMaster(j, icoor, iq) +
                        	      currentfe_.invTrJPhiDerMaster(i, 1, iq) * K_(1, icoor) * currentfe_.invTrJPhiDerMaster(j, icoor, iq));

                   		// Anys. version of s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);
                   	}
                   	return s;
                }

                template<class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0) {return 0;}

                template<class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0) {return 0;}
};

template <>
class StiffAnys<Diffusivity>
{
        private:
                //! A reference to FiniteElement<Integrator>
                /*!
                * Stores a reference to the finite element where the stiffness operator is evaluated.
                */
                //FiniteElement<Integrator, ORDER>& currentfe_;

                const Diffusivity & K_;

        public:
                // Constructors
                //! A constructor.
                /*!
                 \param fe is a reference to FiniteElement<Integrator>
                 */
                StiffAnys(const Diffusivity & K): K_(K){};

                // Operators
                //! A definition of operator () taking two arguments.
                /*!
                * Evaluates the stiffness operator (i,j) of the current finite elemente.
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
                //[TODO]

                //! A definition of operator () taking four arguments.
                /*!
                * Evaluates the product of: the derivative of basis(i) with respect to coordinate ic1 and the derivative of basis(j) with respect
                * to coordinate ic2 ,on current finite elemente.
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
                * \param ic1 is an unsigned int, the variable respect whom the derivative is take: ic1=0 abscissa, ic1=1 ordinata
                * returns a double.
                */
                template <class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 2> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                   	Real s = 0;
                   	for (UInt icoor = 0; icoor < 2; icoor++)
                   	{
                   		UInt globalIndex = currentfe_.getGlobalIndex(iq); // Find global index of quadrature node (K may be space-dependent)
                   		s += (currentfe_.invTrJPhiDerMaster(i, 0, iq) * K_(globalIndex)(0, icoor) * currentfe_.invTrJPhiDerMaster(j, icoor, iq) +
                		      currentfe_.invTrJPhiDerMaster(i, 1, iq) * K_(globalIndex)(1, icoor) * currentfe_.invTrJPhiDerMaster(j, icoor, iq));

                   		// Anys. version of s += currentfe_.invTrJPhiDerMaster(i, icoor, iq)*currentfe_.invTrJPhiDerMaster(j, icoor, iq);
                   	}
                   	return s;
                }

                template <class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0) {return 0;}

                template <class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0) {return 0;}
};

//----------------------------------------------------------------------------//

// **** MASS ****
//! A Mass class: a class for the mass operator. [Note term c (reaction) is absent]
class Mass
{
        public:
                // Operators
                //! A definition of operator () taking three arguments.
                /*!
                * Evaluates the mass operator (i,j) of the current finite element.
                * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,2>, current planar finite element
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
                template <class Integrator ,UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 2> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                	return currentfe_.phiMaster(i, iq) * currentfe_.phiMaster(j, iq);
                }

                //! A definition of operator () taking three arguments.
                /*!
                * Evaluates the mass operator (i,j) of the current finite element.
                * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,2,3>, current planar finite element
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
                template <class Integrator ,UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER  2, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                	return currentfe_.phiMaster(i, iq) * currentfe_.phiMaster(j, iq);
                }

                //! A definition of operator () taking three arguments.
                /*!
                * Evaluates the mass operator (i,j) of the current finite element.
                * \param currentfe_ is an object of class FiniteElement<Integrator, ORDER,3,3>, current planar finite element
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * returns a double.
                */
                template <class Integrator ,UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                	return currentfe_.phiMaster(i, iq) * currentfe_.phiMaster(j, iq);
                }
};

//----------------------------------------------------------------------------//

// **** VECTORIAL GRADIENT ****
//! A vGrad class: a class for the the vectorial Gradient operator. [Note term b (advetion) is absent]
class Grad
{
	public:
                // Operators
                //! A definition of operator () taking three arguments.
                /*!
                * Evaluates the component ic of the vGrad operator (i,j) on the current finite elemente.
                * \param i is an unsigned int, current finite element local index
                * \param j is an unsigned int, current finite element local index
                * \param ic is an unsigned int, vGrad component to be evaluated
                * returns a double.
                */
                template<class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 2> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                	 return currentfe_.phiMaster(i, iq) * currentfe_.invTrJPhiDerMaster(j, ic, iq);
                }

                // [TODO]: following are dummy versions, to be implemented

                template<class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 2, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0) {return 0;};

                template<class Integrator, UInt ORDER>
                inline Real operator() (FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0) {return 0;};
};

//----------------------------------------------------------------------------//

// **** EXPRESSION TEMPLATES FOR ASSEMBLER ****
// Generic template class wrapper
//! A ETWrapper class: Expression Template Wrapper.
/*!
 * Class that mimic the behaviour of a generic operator defined above: following
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
template<typename A>
class EOExpr
{
	private:
                //! "A" is a generic type
                A a_;

	public:
                // Constructors
                //! A constructor.
                /*!
                * \param object is a constant reference to a generic operator.
                */
                EOExpr(const A & a): a_(a){};

                // Operators
                //! A definition of operator () which takes two arguments.
                /*!
                * Masks the behaviour of the correspondent operator in the above classes.
                * \param i is an unsigned int
                * \param j is an unsigned int
                * returns a P variable.
                */
                //P operator() (UInt i, UInt j) {return a_(i,j);}

                //! A definition of operator () which takes three arguments.
                /*!
                * Masks the behaviour of the correspondent operator in the above classes.
                * \param i is an unsigned int
                * \param j is an unsigned int
                * \param ic is an unsigned int
                * returns a P variable.
                */
                //P operator() (UInt i, UInt j, UInt ic) { return a_(i,j,ic);

                //! A definition of operator () which takes four arguments.
                /*!
                * Masks the behaviour of the correspondent operator in the above classes.
                * \param i is an unsigned int
                * \param j is an unsigned int
                * \param ic1 is an unsigned int
                * \param ic2 is an unsigned int
                * returns a P variable.
                */
                //P operator() (UInt i, UInt j, UInt ic1, UInt ic2) { return a_(i,j,ic1,ic2);}

                EOExpr<StiffAnys<Eigen::Matrix<Real, 2, 2>>>    operator[] (const Eigen::Matrix<Real, 2, 2> & K)
                {
                        typedef EOExpr<StiffAnys<Eigen::Matrix<Real, 2, 2>>> ExprT;
                        StiffAnys<Eigen::Matrix<Real, 2, 2>> anys(K);
                        return ExprT(anys);
                        //StiffAnys<Eigen::Matrix<Real,2,2> > a(K);
                }

                EOExpr<StiffAnys<Diffusivity>>                 operator[] (const Diffusivity & K)
                {
                        typedef EOExpr<StiffAnys<Diffusivity> > ExprT;
                        StiffAnys<Diffusivity> anys(K);
                        return ExprT(anys);
                        //return EOExpr<P,A>(A(K));
                }

                template<typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
                Real operator() (FiniteElement<Integrator, ORDER, mydim, ndim> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                        return a_(currentfe_, i, j, iq, ic);
                }
};

// Composition of two wrappers (operator)
//! A ETWBinOp class: Expression Template Wrapper Binary Operation
/*!
 * Class that implements an abstract binary operation defined by Op between two ETWrappers, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
template<typename A, typename B, typename Op>
class EOBinOp
{
	private:
		//! "A" is a generic type.
		/*!
		 * Stores the first operand.
		 */
		A a_;

		//! "B" is a generic type.
		/*!
		 * Stores the second operand.
		 */
		B b_;

	public:
                // Constructor
                //! A constructor.
                /*!
                * \param a is a constant reference to a generic type.
                * \param b is a constant reference to a generic type.
                */
                EOBinOp(const A & a ,const B & b): a_(a),b_(b){};

                // Operators
                //! A definition of operator () taking two arguments.
                /*!
                * \param i is an unsigned int
                * \param j is an unsigned int
                * applies the generic operation defined by the type Op to the two generic objects a_, b_;
                * returns a type P variable
                */
                template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
                Real operator () (FiniteElement<Integrator, ORDER, mydim, ndim> & currentfe_, UInt i, UInt j, UInt iq, UInt ic = 0)
                {
                        return Op::apply(a_(currentfe_, i, j, iq, ic), b_(currentfe_, i, j, iq, ic));
                }
};

template<class B, class Op>
class EOBinOp<Real, B, Op>
{
        private:
        	Real   M_a;
        	B      M_b;

        public:
                // Constructor
        	EOBinOp(Real a, const B& b): M_a(a),M_b(b) {};

                // Operators
        	template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
        	inline Real operator()(FiniteElement<Integrator, ORDER, mydim, ndim> & currentfe_, int i, int j, int iq, int ic = 0)
        	{
        		return Op::apply(M_a,M_b(currentfe_, i, j, iq, ic));
        	}
};

template<class B, class Op>
class EOBinOp<Function, B, Op>
{
        private:
        	const Function & M_a;
        	B M_b;

        public:
                // Constructor
        	EOBinOp(const Function & a, const B & b): M_a(a),M_b(b) {};

                // Operators
        	template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
        	inline Real operator()(FiniteElement<Integrator, ORDER, mydim, ndim> & currentfe_, int i, int j, int iq, int ic = 0)
        	{
        		UInt globalIndex = currentfe_.getGlobalIndex(iq);
        		return Op::apply(M_a(globalIndex), M_b(currentfe_, i, j, iq, ic));
        	}
};

// Wrappers addition
//! A ETWAdd class: Expression Template Wrapper Addition
/*!
 * Class that defines Addition operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
class EOAdd
{
	public:
                /// Constructor
        	//! A constructor.
        	EOAdd(){}

                // Operators
        	//! A stastic inline method taking two arguments.
        	/*!
        	 *The actual addition operation
        	 * \param a is of P type, first addend
        	 * \param b is of P type, second addend
        	 */
        	static inline Real apply(Real a, Real b) {return (a+b);}
};

// Multiplication by real scalar
//! A ETWMult class: Expression Template Wrapper Multiplication.
/*!
 * Class that defines Multiplication operation, following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 */
class EOMult
{
	public:
                // Constructor
        	//! A constructor
        	EOMult(){}

                // Operators
        	//! A stastic inline method taking two arguments.
        	/*!
        	 * The actual multiplication operation.
        	 * \param a is of P type, first operand
        	 * \param b is a Real, second operand
        	 */
        	static inline Real apply(Real a, Real b) {return (a*b);}

                // Note:
        	// not needed since in ETRBinOp I did "Op::apply(a_(i,j),b_)"
        	// static inline P apply(Real b, const P & a) {return (a*b);}
};

// Dot
template<class A, class B>
class EODotProd
{
        private:
        	A a_;
        	B b_;

	public:
                // Constructor
        	EODotProd(const A & a, const B & b): a_(a), b_(b) {};

                // Operators
        	template<typename Integrator, UInt ORDER, UInt ndim>
        	inline Real operator()(FiniteElement<Integrator, ORDER, 2, ndim> & currentfe_, int i, int j, int iq, int ic = 0)
        	{
        		Real s = 0.;
        		for (ic = 0; ic < 2; ++ic)
        			s += a_(ic)*b_(currentfe_, i, j, iq, ic);
                        return s;
        	}

        	template<typename Integrator, UInt ORDER>
        	inline Real operator()(FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, int i, int j, int iq, int ic = 0)
        	{
        		Real s = 0.;
        		for (ic = 0; ic < 3; ++ic)
        			s += a_(ic)*b_(currentfe_, i, j, iq, ic);
                        return s;
        	}
};

//Dot
template<class B>
class EODotProd<Function, B>
{
        private:
        	const Function & M_a;
        	B                M_b;

        public:
                // Constructor
        	EODotProd(const Function & a, const B & b): M_a(a), M_b(b) {};

                // Operators
        	template<typename Integrator, UInt ORDER, UInt ndim>
        	inline Real operator()(FiniteElement<Integrator, ORDER, 2, ndim> & currentfe_, int i, int j, int iq, int ic = 0)
        	{
        		Real s = 0.;
        		UInt globalIndex = currentfe_.getGlobalIndex(iq);
        		for (ic = 0; ic < 2; ic++)
        			s += M_a(globalIndex, ic)*M_b(currentfe_, i, j, iq, ic);
                        return s;
        	}

        	template<typename Integrator, UInt ORDER>
        	inline Real operator()(FiniteElement<Integrator, ORDER, 3, 3> & currentfe_, int i, int j, int iq, int ic = 0)
        	{
        		Real s = 0.;
        		UInt globalIndex = currentfe_.getGlobalIndex(iq);
        		for (ic = 0; ic < 3; ic++)
        			s += M_a(globalIndex, ic)*M_b(currentfe_, i, j, iq, ic);
                        return s;
        	}
};

// Operator +
//! Overloading of operator +.
/*!
 * Following:
 * "Expression Templates Implementation of Continuous and DIscontinous Galerkin Methods"
 * D.A. Di Pietro, A. Veneziani
 * Takes two arguments:
 * \param a is const reference ETWrapper<P, A>
 * \param b is const reference ETWrapper<P, A>
 * \return a ETWrapper<P,ETWBinOp<P, ETWrapper<P,A>, ETWrapper<P, B>, ETWAdd<P> > which is resolved at compile time.
 */
template<typename A, typename B>
EOExpr<EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd>>
operator + (const EOExpr<A> &  a, const EOExpr<B> &  b)
{
	  typedef EOBinOp<EOExpr<A>, EOExpr<B>, EOAdd> ExprT;
	  return EOExpr<ExprT> (ExprT(a, b));
}

template<typename B>
EOExpr<EOBinOp<Function, EOExpr<B>, EOMult>>
operator * (const Function & a, const EOExpr<B> &  b)
{
	  typedef EOBinOp<Function, EOExpr<B>, EOMult> ExprT;
	  return EOExpr<ExprT> (ExprT(a, b));
}

template<typename B>
EOExpr<EOBinOp<Real, EOExpr<B>, EOMult>>
operator * (Real a, const EOExpr<B> &  b)
{
	  typedef EOBinOp<Real, EOExpr<B>, EOMult> ExprT;
	  return EOExpr<ExprT> (ExprT(a, b));
}

template<typename B>
EOExpr<EODotProd<Eigen::Matrix<Real, 2, 1>, EOExpr<B>>>
dot(const Eigen::Matrix<Real, 2, 1> & a, const EOExpr<B> & b)
{
	  typedef EODotProd<Eigen::Matrix<Real, 2, 1>, EOExpr<B>> ExprT;
	  return EOExpr<ExprT> (ExprT(a, b));
}

template<typename B>
EOExpr<EODotProd<Eigen::Matrix<Real, 3, 1>, EOExpr<B>>>
dot(const Eigen::Matrix<Real, 3, 1> & a, const EOExpr<B> & b)
{
	  typedef EODotProd<Eigen::Matrix<Real, 3, 1>, EOExpr<B>> ExprT;
	  return EOExpr<ExprT> (ExprT(a, b));
}

template<typename B>
EOExpr<EODotProd<Function, EOExpr<B>>>
dot(const Function & a, const EOExpr<B> & b)
{
	  typedef EODotProd<Function, EOExpr<B>> ExprT;
	  return EOExpr<ExprT> (ExprT(a, b));
}

//----------------------------------------------------------------------------//
// **** ASSEMBLER ****
//! A Assmbler class: discretize a generic differential operator in a sparse matrix
//! template<UInt mydim, UInt ndim>
class Assembler
{
	public:
                // Constructor
                //! A constructor
                //Assembler (){};

                // Assemblers
                // 2D -> 2D
                //! A template member taking three arguments: discretize differential operator
                /*!
                * \param oper is a template expression : the differential operator to be discretized.
                * \param mesh is const reference to a MeshHandler<ORDER,2,2>: the mesh where we want to discretize the operator.
                * \param fe is a const reference to a FiniteElement
                * stores the discretization in SPoper_mat_
                */
                //Return triplets vector
                template<UInt ORDER, typename Integrator, typename A>
                static void operKernel( EOExpr<A> oper, const MeshHandler<ORDER, 2, 2> & mesh,
                                        FiniteElement<Integrator, ORDER, 2, 2> & fe, SpMat & OpMat);

                template<UInt ORDER, typename Integrator>
                static void forcingTerm(const MeshHandler<ORDER, 2, 2> & mesh,
                                        FiniteElement<Integrator, ORDER, 2, 2> & fe,
                                        const ForcingTerm & u, VectorXr & forcingTerm);

                // 2D -> 3D
                //! A template member taking three arguments: discretize differential operator
                /*!
                * \param oper is a template expression : the differential operator to be discretized.
                * \param mesh is const reference to a MeshHandler<ORDER,2,3>: the mesh where we want to discretize the operator.
                * \param fe is a const reference to a FiniteElement
                * stores the discretization in SPoper_mat_
                */
                template<UInt ORDER, typename Integrator, typename A>
                static void operKernel( EOExpr<A> oper, const MeshHandler<ORDER, 2, 3> & mesh,
                                        FiniteElement<Integrator, ORDER, 2, 3> & fe, SpMat & OpMat);

                template<UInt ORDER, typename Integrator>
                static void forcingTerm(const MeshHandler<ORDER, 2, 3> & mesh,
                                        FiniteElement<Integrator, ORDER, 2, 3> & fe,
                                        const ForcingTerm & u, VectorXr & forcingTerm);

                // 3D -> 3D
                template<UInt ORDER, typename Integrator, typename A>
                static void operKernel( EOExpr<A> oper, const MeshHandler<ORDER, 3, 3> & mesh,
                                        FiniteElement<Integrator, ORDER, 3, 3> & fe, SpMat & OpMat);

                template<UInt ORDER, typename Integrator>
                static void forcingTerm(const MeshHandler<ORDER, 3, 3> & mesh,
                                        FiniteElement<Integrator, ORDER, 3, 3> & fe,
                                        const ForcingTerm & u, VectorXr & forcingTerm);
};

// Alternative impleentations
/*
template<>
class Assembler<2,2>{
	private:

	public:
	  //! A constructor
	  //Assembler (){};
	  //! A template member taking three arguments: discretize differential operator
	  //!
	  // \param oper is a template expression : the differential operator to be discretized.
	  // \param mesh is const reference to a MeshHandler<ORDER>: the mesh where we want to discretize the operator.
	  // \param fe is a const reference to a FiniteElement
	  // stores the discretization in SPoper_mat_


	  //Return triplets vector
	  template<UInt ORDER, typename Integrator, typename A>
	  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,2>& mesh,
	  	                     FiniteElement<Integrator, ORDER,2,2>& fe, SpMat& OpMat);

	  template<UInt ORDER, typename Integrator>
	  static void forcingTerm(const MeshHandler<ORDER,2,2>& mesh, FiniteElement<Integrator, ORDER,2,2>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

	};


template<>
class Assembler<2,3>{
	private:

	public:
	  //! A constructor
	  //Assembler (){};
	  //! A template member taking three arguments: discretize differential operator
	  //!
	  // \param oper is a template expression : the differential operator to be discretized.
	  // \param mesh is const reference to a MeshHandler<ORDER>: the mesh where we want to discretize the operator.
	  // \param fe is a const reference to a FiniteElement
	  // stores the discretization in SPoper_mat_
	  //

	  //Return triplets vector
	  template<UInt ORDER, typename Integrator, typename A>
	  static void operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,3>& mesh,
	  	                     FiniteElement<Integrator, ORDER,2,3>& fe, SpMat& OpMat);

	  template<UInt ORDER, typename Integrator>
	  static void forcingTerm(const MeshHandler<ORDER,2,3>& mesh, FiniteElement<Integrator, ORDER,2,3>& fe, const ForcingTerm& u, VectorXr& forcingTerm);

	};*/

#include "matrix_assembler_imp.h"

#endif
