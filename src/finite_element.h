#ifndef __FINITE_ELEMENT_HPP__
#define __FINITE_ELEMENT_HPP__

// HEADERS
#include "fdaPDE.h"
#include "integration.h"
#include "mesh_objects.h"

// CLASSES
//!  This class implements all properties of a Triangular or Tetrahedral Finite Element
/*!
 * This class is the most important one of the entire code
 * and implements everything needed by a triangular or tetrahedral finite elemnt
 *
 * It takes as a template parameter a class that implements the method used
 * for determining the mass, stiff and grad matrices
 *
 * the order is the polynomial order needed [r], mydim the dimension of the element
 * ndim the dimension of the space in which t is embedded [d]
 *
 * Here are also computed values for the integration matrices R0 and R1, in particular:
 * 1) the evaluation of \hat{phi} in the integration nodes
 * 2) the evaluation of the gradient of \hat{phi} in the integration NODES
 * If we consider, on a generic node, I := int_{Omega}(<grad(phi_i(p)), grad(phi_j(p))> dp)
 * Letting T: Omega_Master -> Omega s.t p = T(p_hat), (here linear) and invertible
 * Then I ==
 * int{Omega_Master}(<J_{T^-1}^t * grad(\hat{phi}_i(p_hat)), J_{T^-1}^t * grad(\hat{phi}_j(p_hat))> * det(J_{T}) dp_hat)
 * Thus we also compute:
 * 3) J_{T^-1}^t * grad(\hat{phi}_i(p_hat)) for all i in 1:mydim
 */
template <class Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FiniteElement
{
};

// *** TRIANGLE in 2D***
// Remark in a triangle in d == 2 the relationship between N_r and r is (r+1)*(r+2)/2
// that in our case r = 1,2 can be simplified with N_r == 3*r [not extendable!!]
template <class Integrator, UInt ORDER>
class FiniteElement<Integrator, ORDER, 2, 2>
{
        private:
                Element<3*ORDER, 2, 2> reference_;  // reference [master] element
                Element<3*ORDER, 2, 2> t_;

                // Matrices of evaluations on the master element:
                // of the i-th \hat{phi} on the iq-th integration node
                // ex i-row: [\hat{phi}_i(qn_1), ...,  \hat{phi}_i(qn_NNODES)]
                Eigen::Matrix<Real, 3*ORDER, Integrator::NNODES>   phiMapMaster_;
                // of the i-th (\hat{phi}_x, \hat{phi}_y) on the iq-th integration node (grad_hat(\hat{phi}))
                // ex i-row: [\hat{phi}_x_i(qn_1), \hat{phi}_y_i(qn_1), ..., \hat{phi}_x_i(qn_NNODES), \hat{phi}_y_i(qn_NNODES)]
                Eigen::Matrix<Real, 3*ORDER, Integrator::NNODES*2> phiDerMapMaster_;
                // of the i-th (J_{T^-1}^t*\hat{phi}_x, J_{T^-1}^t*\hat{phi}_y) on the iq-th integration node
                // ex i-row: [J_{T^-1}^t*\hat{phi}_x_i(qn_1), J_{T^-1}^t*\hat{phi}_y_i(qn_1), ..., J_{T^-1}^t*\hat{phi}_x_i(qn_NNODES), J_{T^-1}^t*\hat{phi}_y_i(qn_NNODES)]
                Eigen::Matrix<Real, 3*ORDER, Integrator::NNODES*2> invTrJPhiDerMapMaster_;

                // Setters for constructor
                void setPhiMaster();
                void setPhiDerMaster();
                void setInvTrJPhiDerMaster();
        public:
                // Constructor
                //! This is an empty constructor
                /*!
                * For efficiency and Expression Templates organization of the
                * code, the use of this class is based on the updateElement class
                */
                FiniteElement();

                // Updater
                //! A member updating the Finite Element properties
                /*!
                * \param t an element from which to update the finite element properties
                */
                void updateElement(const Element<3*ORDER, 2, 2> & t);

                // Getters
                Real getAreaReference()     {return reference_.getArea();}
                Real getDet()               {return t_.getDetJ();}

                //General utilities
                Point coorQuadPt(UInt iq)
                {
                return Point(t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1] + t_[0][0],
                             t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1] + t_[0][1]);
                }

                UInt getGlobalIndex(UInt iq) {return Integrator::NNODES * t_.getId() + iq;}

                // Access
                //Returns \hat{phi}
                Real phiMaster(UInt i, UInt iq)                   const;
                //Returns \nabla \hat{phi}
                Real phiDerMaster(UInt i, UInt ic, UInt iq)       const;
                //Returns J^{-T} \nabla \hat{phi}
                Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
};

// *** TRIANGLE in 3D***
// Remark in a triangle in d_master == 2 the relationship between N_r and r is (r+1)*(r+2)/2
// that in our case r = 1,2 can be simplified with N_r == 3*r [not extendable!!]
template <class Integrator ,UInt ORDER>
class FiniteElement<Integrator, ORDER, 2, 3>
{
        private:
                Element<3*ORDER, 2, 2> reference_; // reference [master] element
                Element<3*ORDER, 2, 3> t_;

                // Matrices of evaluations on the master element:
                // of the i-th \hat{phi} on the iq-th integration node
                // ex i-row: [\hat{phi}_i(qn_1), ...,  \hat{phi}_i(qn_NNODES)]
                Eigen::Matrix<Real, 3*ORDER, Integrator::NNODES>       phiMapMaster_;
                // of the i-th (\hat{phi}_x, \hat{phi}_y) on the iq-th integration node (grad_hat(\hat{phi}))
                // ex i-row: [\hat{phi}_x_i(qn_1), \hat{phi}_y_i(qn_1), ..., \hat{phi}_x_i(qn_NNODES), \hat{phi}_y_i(qn_NNODES)]
                Eigen::Matrix<Real, 3*ORDER, Integrator::NNODES*2>     phiDerMapMaster_;

                // Problem of invertibility of T... to be fixed
                // [TODO] Eigen::Matrix<Real,3*ORDER, Integrator::NNODES*2> invTrJPhiDerMapMaster_;
                Eigen::Matrix<Real,2,2> metric_;

                // Setters for constructor
                void setPhiMaster();
                void setPhiDerMaster();
                // [TODO] void setInvTrJPhiDerMaster();

        public:
                // Constructor
                //! This is an empty constructor
                /*!
                For efficiency and Expression Templates organization of the
                code, the use of this class is based on the updateElement class
                */
                FiniteElement();

                // Updater
                //! A member updating the Finite Element properties
                /*!
                * \param t an element from which to update the finite element properties
                */
                void updateElement(const Element<3*ORDER, 2, 3> & t);

                // Getters
                Real getAreaReference() {return reference_.getArea();}
                Real getDet()           {return t_.getDetJ();}

                Point coorQuadPt(UInt iq)
                {
                return Point(   t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1] + t_[0][0],
                                t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1] + t_[0][1],
                                t_.getM_J()(2,0)*Integrator::NODES[iq][0] + t_.getM_J()(2,1)*Integrator::NODES[iq][1] + t_[0][2]);
                }

                UInt getGlobalIndex(UInt iq) {return Integrator::NNODES * t_.getId() + iq;}

                // Access
                //Returns \hat{phi}
                Real phiMaster(UInt i, UInt iq)                   const;
                //Returns \nabla \hat{phi}
                Real phiDerMaster(UInt i, UInt ic, UInt iq)       const;
                //Returns J^{-1} \nabla \hat{phi}
                //[TODO]: still to be implemented
              //Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
                // Metric
                Eigen::Matrix<Real,2,2> metric()                  const {return metric_;};
};

// *** TETRAHEDRON ***
//Implementation FiniteElement with mydim == 3 & ndim == 3
template <class Integrator ,UInt ORDER>
class FiniteElement<Integrator, ORDER, 3, 3>
{
        private:
                Element<6*ORDER-2, 3, 3> reference_;
                Element<6*ORDER-2, 3, 3> t_;

                // Matrices of evaluations on the master element:
                // of the i-th \hat{phi} on the iq-th integration node
                // ex i-row: [\hat{phi}_i(qn_1), ...,  \hat{phi}_i(qn_NNODES)]
                Eigen::Matrix<Real, 6*ORDER-2, Integrator::NNODES>      phiMapMaster_;
                // of the i-th (\hat{phi}_x, \hat{phi}_y) on the iq-th integration node (grad_hat(\hat{phi}))
                // ex i-row: [\hat{phi}_x_i(qn_1), \hat{phi}_y_i(qn_1), \hat{phi}_z(qn_i) ..., \hat{phi}_x_i(qn_NNODES), \hat{phi}_y_i(qn_NNODES), \hat{phi}_z_i(qn_NNODES)]
                Eigen::Matrix<Real, 6*ORDER-2, Integrator::NNODES*3>    phiDerMapMaster_;
                // of the i-th (J_{T^-1}^t*\hat{phi}_x, J_{T^-1}^t*\hat{phi}_y,  J_{T^-1}^t*\hat{phi}_z) on the iq-th integration node
                // ex i-row: [J_{T^-1}^t*\hat{phi}_x_i(qn_1), J_{T^-1}^t*\hat{phi}_y_i(qn_1), J_{T^-1}^t*\hat{phi}_z_i(qn_1), ..., J_{T^-1}^t*\hat{phi}_x_i(qn_NNODES), J_{T^-1}^t*\hat{phi}_y_i(qn_NNODES), J_{T^-1}^t*\hat{phi}_z_i(qn_NNODES)]
                Eigen::Matrix<Real, 6*ORDER-2, Integrator::NNODES*3>    invTrJPhiDerMapMaster_;
                Eigen::Matrix<Real, 3, 3> metric_;

                // Setters for constructor
                void setPhiMaster();
                void setPhiDerMaster();
                void setInvTrJPhiDerMaster();

        public:
                // Constructor
                //! This is an empty constructor
                /*!
                * For efficiency and Expression Templates organization of the
                * code, the use of this class is based on the updateElement class
                */
                FiniteElement();

                //Updater
                //! A member updating the Finite Element properties
                /*!
                * \param t an element from which to update the finite element properties
                */
                void updateElement(const Element<6*ORDER-2, 3, 3> & t);

                // Getters
                Real getVolumeReference()       {return reference_.getVolume();}
                Real getDet()                   {return t_.getDetJ();}

                // General utilities
                Point coorQuadPt(UInt iq)
                {
                return Point(   t_.getM_J()(0,0)*Integrator::NODES[iq][0] + t_.getM_J()(0,1)*Integrator::NODES[iq][1] + t_.getM_J()(0,2)*Integrator::NODES[iq][2] + t_[0][0],
                                t_.getM_J()(1,0)*Integrator::NODES[iq][0] + t_.getM_J()(1,1)*Integrator::NODES[iq][1] + t_.getM_J()(1,2)*Integrator::NODES[iq][2] + t_[0][1],
                                t_.getM_J()(2,0)*Integrator::NODES[iq][0] + t_.getM_J()(2,1)*Integrator::NODES[iq][1] + t_.getM_J()(2,2)*Integrator::NODES[iq][2] + t_[0][2]);
                }

                UInt getGlobalIndex(UInt iq) {return Integrator::NNODES*t_.getId() + iq;}

                // Access
                //Returns \hat{phi}
                Real phiMaster(UInt i, UInt iq)                   const;
                //Returns \nabla \hat{phi}
                Real phiDerMaster(UInt i, UInt ic, UInt iq)       const;
                //Returns J^{-1} \nabla \hat{phi}
                Real invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const;
                // Metric
                Eigen::Matrix<Real, 3, 3> metric()                const {return metric_;};
};

#include "finite_element_imp.h"

#endif
