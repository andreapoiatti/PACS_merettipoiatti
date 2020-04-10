#ifndef __REGRESSIONDATA_HPP__
#define __REGRESSIONDATA_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  RegressionData
{
	private:
		std::vector<Point> locations_;			// Vector collecting the sampling locations

		VectorXr observations_;				// Vector containing sampled values in the locations
		std::vector<UInt> observations_indices_;
		bool locations_by_nodes_;			// If true locations are the nodes of the mesh

		//Design matrix
		MatrixXr covariates_;
		UInt n_;
		UInt p_;

		//Areal data
		MatrixXi incidenceMatrix_;
		UInt nRegions_;

		//Order approximating basis
		UInt order_;

		//Boundary conditions
		std::vector<Real> bc_values_;
		std::vector<UInt> bc_indices_;

		#ifdef R_VERSION_ // [[DEPRECAED??]]
		void setLocations(SEXP Rlocations);
		void setObservations(SEXP Robservations);
		void setCovariates(SEXP Rcovariates);
		void setIncidenceMatrix(SEXP RincidenceMatrix);
		#endif

	public:
		//! A basic version of the constructor.
		RegressionData(){};

		#ifdef R_VERSION_
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
		\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
		\param Rnrealizations the number of random points used in the stochastic computation of the dofs
		*/
		explicit RegressionData(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rcovariates, SEXP RincidenceMatrix,
			SEXP RBCIndices, SEXP RBCValues);
		#endif

		explicit RegressionData(std::vector<Point> & locations, VectorXr & observations, UInt order, MatrixXr& covariates,
			MatrixXi & incidenceMatrix, std::vector<UInt> & bc_indices, std::vector<Real>& bc_values);

		// Printers
		void printObservations(std::ostream & out) const;
		void printCovariates(std::ostream & out) const;
		void printLocations(std::ostream & out) const;
		void printIncidenceMatrix(std::ostream & out) const;

		// Getters
		//! A method returning a reference to the observations vector
		inline const VectorXr * getObservations() const {return &observations_;}
		//! A method returning a reference to the design matrix
		inline const MatrixXr * getCovariates() const {return &covariates_;}
		//! A method returning a reference to the incidence matrix
		inline const MatrixXi * getIncidenceMatrix() const {return &incidenceMatrix_;}
		//! A method returning the number of observations
		inline UInt getNumberofObservations() const {return observations_.size();}
		//! A method returning the locations of the observations
		inline const std::vector<Point> * getLocations() const {return &locations_;}
		//! A method returning the number of regions
		inline UInt getNumberOfRegions() const {return nRegions_;}
		inline bool isLocationsByNodes() const {return locations_by_nodes_;}
		inline const std::vector<UInt> * getObservationsIndices() const {return &observations_indices_;}
		//! A method returning the the penalization term
		inline UInt getOrder() const {return order_;}
		//! A method returning the indexes of the nodes for which is needed to apply Dirichlet Conditions
		inline std::vector<UInt>getDirichletIndices() const {return bc_indices_;}
		//! A method returning the values to apply for Dirichlet Conditions
		inline std::vector<Real> getDirichletValues() const {return bc_values_;}
};


class RegressionDataElliptic: public RegressionData
{
	private:
		Eigen::Matrix<Real,2,2> K_;
		Eigen::Matrix<Real,2,1> beta_;
		Real c_;

	public:
		#ifdef R_VERSION_
		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
		*/
		explicit RegressionDataElliptic(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP RK, SEXP Rbeta,
			SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues);
		#endif

		explicit RegressionDataElliptic(std::vector<Point> & locations, VectorXr & observations, UInt order,
			Eigen::Matrix<Real,2,2> & K, Eigen::Matrix<Real,2,1> & beta, Real c, MatrixXr & covariates,
			MatrixXi & incidenceMatrix, std::vector<UInt> & bc_indices, std::vector<Real> & bc_values);

		inline Eigen::Matrix<Real,2,2> getK() const {return K_;}
		inline Eigen::Matrix<Real,2,1> getBeta() const {return beta_;}
		inline Real getC() const {return c_;}
};

class RegressionDataEllipticSpaceVarying: public RegressionData
{
	private:
		Diffusivity K_;
		Advection beta_;
		Reaction c_;
		ForcingTerm u_;

	public:
		#ifdef R_VERSION_
		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rlocations an R-matrix containing the location of the observations.
			\param Robservations an R-vector containing the values of the observations.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param  Ru an R-double vector of length #triangles that contaiins the forcing term integrals.
			\param Rcovariates an R-matrix storing the covariates of the regression
			\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions in the model with areal data
			\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
		*/
		explicit RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP RK, SEXP Rbeta, SEXP Rc,
			 SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues);
		#endif


		explicit RegressionDataEllipticSpaceVarying(std::vector<Point> & locations, VectorXr & observations, UInt order,
			const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2>>> & K,
			const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1>>> & beta,
			const std::vector<Real> & c, const std::vector<Real> & u,
			MatrixXr & covariates, MatrixXi & incidenceMatrix, std::vector<UInt> & bc_indices, std::vector<Real> & bc_values);

		inline Diffusivity getK() const {return K_;}
		inline Advection getBeta() const {return beta_;}
		inline Reaction getC() const {return c_;}
		inline ForcingTerm getU() const {return u_;}

		void print(std::ostream & out) const;
};

#include "regressionData_imp.h"

#endif
