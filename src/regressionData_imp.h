#ifndef __REGRESSION_DATA_IMP_H__
#define __REGRESSION_DATA_IMP_H__
// create GAM constructors with right template inheritances

// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order,
									std::vector<Real> lambdaS, MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search,
									UInt max_num_iterations, Real threshold, Real tune, bool arealDataAvg):
		 RegressionData(locations, observations, order, lambdaS, covariates, incidenceMatrix, bc_indices, bc_values, DOF, false, search, tune, arealDataAvg),
		 				max_num_iterations_(max_num_iterations), threshold_(threshold), initialObservations_(observations), global_lambda_(lambdaS), GCV_GAM_(GCV)
{;}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order,
												std::vector<Real> lambdaS, Eigen::Matrix<Real,2,2>& K,
												Eigen::Matrix<Real,2,1>& beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search,
												UInt max_num_iterations, Real threshold, Real tune, bool arealDataAvg):
		 RegressionDataElliptic(locations, observations, order, lambdaS, K, beta, c, covariates, incidenceMatrix, bc_indices, bc_values, DOF, false, search, tune, arealDataAvg),
		 									max_num_iterations_(max_num_iterations), threshold_(threshold), initialObservations_(observations), global_lambda_(lambdaS), GCV_GAM_(GCV)
{;}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point>& locations,
									VectorXr& observations, UInt order, std::vector<Real> lambdaS,
									const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
									const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
									const std::vector<Real>& c, const std::vector<Real>& u,
									MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF, bool GCV, UInt search,
									UInt max_num_iterations, Real threshold, Real tune, bool arealDataAvg):
		 RegressionDataEllipticSpaceVarying(locations, observations, order, lambdaS, K, beta, c, u, covariates, incidenceMatrix, bc_indices, bc_values, DOF, false, search, tune, arealDataAvg),
		 									max_num_iterations_(max_num_iterations), threshold_(threshold), initialObservations_(observations), global_lambda_(lambdaS), GCV_GAM_(GCV)
{;}


// GAM constructors

// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS,
				SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
				SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch,
				SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rtune, SEXP RarealDataAvg):
	RegressionData(Rlocations, RbaryLocations, Robservations, Rorder, RlambdaS, Rcovariates,
					RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod,
					Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg)
{

	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];

    initialObservations_ = this->observations_;
    global_lambda_ = this->lambdaS_;
    GCV_GAM_ = this->GCV_;
    this->GCV_ = false;
}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS,
				SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
				SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch,
				SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rtune, SEXP RarealDataAvg):
	RegressionDataElliptic(Rlocations, RbaryLocations, Robservations, Rorder, RlambdaS,
				RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];

    initialObservations_ = this->observations_;
    global_lambda_ = this->lambdaS_;
    GCV_GAM_ = this->GCV_;
    this->GCV_ = false;
}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP RbaryLocations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS,
				SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
				SEXP GCV, SEXP RGCVmethod, SEXP Rnrealizations, SEXP DOF, SEXP RDOF_matrix, SEXP Rsearch,
				SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rtune, SEXP RarealDataAvg):
	RegressionDataEllipticSpaceVarying(Rlocations, RbaryLocations, Robservations, Rorder, RlambdaS,
										RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues,
										GCV, RGCVmethod, Rnrealizations, DOF, RDOF_matrix, Rsearch, Rtune, RarealDataAvg)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
    initialObservations_ = this->observations_;
    global_lambda_ = this->lambdaS_;
    GCV_GAM_ = this->GCV_;
    this->GCV_ = false;
}

#endif
