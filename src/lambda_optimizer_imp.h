#ifndef _GCV_IMP_
#define _GCV_IMP_

#include <iostream>
#include <vector>

// *** GCV-based ***
// Setters

//! Method to set the value of member SpMat R_
/*
 * /remark {R = R1^t * R0^{-1} * R1 therefore is NOT dependent on \lambda }
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::set_R_(void)
{
        SpMat R1_= this->model.getR1_();                               // Get the value of matrix R1
        Sparse_LU solver;	                                 // define a factorized empty sparse Cholesky solver
        solver.compute(this->model.getR0_());		                 // apply it to R0 to simplify the inverse
        R_ = R1_.transpose()*solver.solve(R1_);                  // R == _R1^t*R0^{-1}*R1
        R_.makeCompressed();                                     // Compress the matrix to speed up the next operations
}

//! Method to set the value of member SpMat T_
/*
 * /remark {T = D + \lambda * R where D is the top-left block of the matrix DMat}
 * /pre {set_R_ must be called before set_T_, the matrix D_ [DataMatrix] must be constructed in the model, s must be defined}
 * /sa {set_R_(void), getDataMatrix(SpMat & DMat)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::set_T_(Real lambda)
{
        // THEORETICAL REMARK:
        // It's easy to check that, since Psi is a rectangular permutation matrix,
        // multiplying Psi^t*Psi produces a [n_nodes x n_nodes] diagonal matrix:
        // there is a 1 in all the indices for which there is a permutation column
        // in Psi and 0 otherwise.
        // More formally:
        // if function k: loctions -> nodes s.t. Psi = Indicator(i,k[i]) then
        // (1) Psi^t*Psi   == Indicator(k[i],k[i])
        // (2) Psi^t*Q*Psi == q_ij*Indicator(k[i],k[j]).

        // IMPLEMENTATION OF THE REMARK:
        // When nodes and locations are cohincident is more advantageous
        // to avoid a direct sum with the full matrix D_ taking only the
        // non-null components on its diagonal and subsequently summing
        // them direcly to the second block of T: lambda*R.

        if (this->model.checkhasLocations_())
        {
                // Full this->model, no simplification allowed
                const SpMat & D_ = this->model.get_DMat_();
                T_ = D_+lambda*R_;
        }
        else
        {
                // D_ non-null entries are added directly to the structure
                T_ = lambda*R_;

                std::vector<UInt> k  = this->model.getData()->getObservationsIndices();
                if (this->model.checkisRegression_())
                {
                        // T = Psi^t*Q*Psi == q_ij*Indicator(k[i],k[j])
                        const MatrixXr &  Q_ = this->model.getQ_();
                        for (UInt i = 0; i < s ; i++)
                                for (int j = 0; j < s; j++)
                                        T_.coeffRef(k[i], k[j]) += Q_.coeff(i, j);
                }
                else
                {
                        // T = Psi^t*Psi == Indicator(k[i],k[i])
                        for (UInt i = 0; i < s ; i++)
                        {
                                T_.coeffRef(k[i], k[i]) += 1;
                        }
                }
        }

        T_.makeCompressed();    // Compressing the matrix for further computation
}

//! Method to set the value of member SpMat V_
/*
 * /remark {V = T^{-1}*Psi^t*Q }
 * /pre {set_T_ must be called before set_V_}
 * /sa {set_T_(Real lambda)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::set_V_(void)
{
        const SpMat &    Psi_   = this->model.getPsi_();      // Get matrix Psi_
        const MatrixXr & Q_     = this->model.getQ_();        // Get matrix Q_

        if (this->model.checkisRegression_())
        {
                // Q != I
                SpMat E_; // = Psi^t*Q               // Declare an empty auxiliary matrix
                if (this->model.checkhasLocations_())
                {
                        // Psi is full
                        E_ = (SpMat(Psi_.transpose())*Q_).sparseView();
                }
                else
                {
                        // Psi is permutation

                        // THEORETICAL REMARK:
                        // Since Psi is a rectangular permutation matrix, if function
                        // k: loctions -> nodes s.t. Psi = Indicator(i,k[i]) then
                        // Psi^t*Q   == Indicator(k[i],j)*q_ij

                        // IMPLEMENTATION OF THE REMARK:
                        // the number of non-null entries of E is at most s^2,
                        // we reserve a vector containing such entries and
                        // we set the final matrix from these triplets
                        E_.resize(Psi_.cols(), s);
                        std::vector<coeff> vec;
                        vec.reserve(s*s);
                        std::vector<UInt> k = this->model.getData()->getObservationsIndices();
                        for (UInt i = 0; i < s ; i++)
                                for (UInt j = 0; j < s; j++)
                                        vec.push_back(coeff(k[i], j, Q_.coeff(i, j)));
                        E_.setFromTriplets(vec.begin(), vec.end());
                }

                E_.makeCompressed();                      // compress the matrix to save space

                Sparse_LU solver;	                  // define a factorized empty sparse LU solver
                solver.compute(T_);                        // apply it to T to simplify the inverse
                V_ = solver.solve(E_);                    // find the value of V = T^{-1}*E
        }
        else
        {
                // Q == I
                Sparse_LU solver;	                  // define a factorized empty sparse LU solver
                solver.compute(T_);                        // apply it to T to simplify the inverse
                V_ = solver.solve(SpMat(Psi_.transpose()));      // find the value of V = T^{-1}*Psi^t

        }

        V_.makeCompressed();            // compress the matrix to save space
}

//! Method to set the value of member SpMat S_ and its trace trS_
/*
 * /remark {S = Psi*V }
 * /pre {set_V_ must be called before set_S_}
 * /sa {set_V_(void)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::set_S_and_trS_(void)
{
        // S = Psi*V
        trS_ = 0.0;
        this->LeftMultiplybyPsiAndTrace(trS_, S_, V_);
        S_.makeCompressed();    // Compress the matrix to save space
}

//! Method to set the value of member SpMat dS_ and its trace trdS_, also computes utility matrix K_
/*
 * /remark {dS_ = -Psi*T^{-1}*R*V}
 * /pre {definition of base matrices R_, T_, V_}
 * /sa {set_V_(void)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::set_dS_and_trdS_(void)
{
        // dS_ = -Psi*(Psi^t*Q*Psi+lambda*R1^t*R0^-1*R1)^{-1}*R1^t*R0^{-1}*R1*(Psi^t*Q*Psi+lambda*R1^t*R0^{-1}*R1)^{-1}*Psi^t*Q
        //     = -Psi*T^{-1}*R*V
        //     =  Psi*(-K*V)
        //    :=  Psi*F
        Sparse_Cholesky solver;	                // define a factorized empty sparse Cholesky solver
        solver.compute(T_);                     // apply it to T_ to simplify the inverse
        K_ = solver.solve(R_);                   // K = T^{-1}*R
        SpMat F_ = -K_*V_;                        // F = -K*V
        trdS_ = 0.0;

        this->LeftMultiplybyPsiAndTrace(trdS_, dS_, F_);
        dS_.makeCompressed();    // Compress the matrix to save space
}

//! Method to set the value of member SpMat ddS_ and its trace trddS_
/*
 * /remark {ddS_ = -2*Psi*K^2*V}
 * /pre {definition of V_ and K_}
 * /sa {set_V_(void), set_dS_and_trdS_(void)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::set_ddS_and_trddS_(void)
{
        // ddS_ = -
        //      = -
        //     := 2*Psi*K^2*V
        SpMat G_ = 2*K_*K_*V_;                        // G = -2*K^2*V
        trddS_ = 0.0;

        this->LeftMultiplybyPsiAndTrace(trddS_, ddS_, G_);
        ddS_.makeCompressed();    // Compress the matrix to save space
}

// Utilities

//! Utility to left multiply a matrix by Psi_ and compute the trace of the new matrix
/*
 * /param trace real where to store the value of the computed matrix
 * /param ret sparse matrix where to store the computed matrix
 * /param mat matrix to left multiply by Psi_
 * /sa {set_S_and_trS_(void), set_dS_and_trdS_(void), set_ddS_and_trddS_(void)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::LeftMultiplybyPsiAndTrace(Real & trace, SpMat & ret, const SpMat & mat)
{
        const SpMat & Psi_ = this->model.getPsi_();   // Get matrix Psi

        if (this->model.checkhasLocations_())
        {
                // Psi is full
                ret = Psi_*mat;
                for (int k = 0; k < s; ++k)
                        trace += ret.coeff(k, k);

        }
        else
        {
                // Psi is permutation

                // THEORETICAL REMARK:
                // Since Psi is a rectangular permutation matrix, if function
                // k: loctions -> nodes s.t. Psi = Indicator(i,k[i]) then
                // Psi*F   == Indicator(i,j)*f_{k[i]j}

                // IMPLEMENTATION OF THE REMARK:
                // the number of non-null entries of E is at most s^2,
                // we reserve a vector containing such entries and
                // we set the final matrix from these triplets

                std::vector<coeff> vec;
                vec.reserve(s*s);
                ret.resize(s, s);

                std::vector<UInt> k = this->model.getData()->getObservationsIndices();
                for (UInt i = 0; i < s; i++)
                {
                        for (UInt j = 0; j < s; j++)
                        {
                                if (i == j)
                                {
                                        Real v = mat.coeff(k[i], j);
                                        vec.push_back(coeff(i, j, v));
                                        trace += v;
                                }
                                else
                                {
                                        vec.push_back(coeff(i, j, mat.coeff(k[i], j)));
                                }
                        }
                }

                ret.setFromTriplets(vec.begin(), vec.end());
        }
}

//! Utility to compute the predicted values in the locations
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_z_hat(void)
{
        if(this->model.checkisRegression_())
                z_hat = (this->model.getH_()+this->model.getQ_()*S_)*this->model.getData()->getObservations();
        else
                z_hat = S_*this->model.getData()->getObservations();
}

//! Utility to compute the predicted residuals in the locations
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_eps_hat(void)
{
        eps_hat = this->model.getData()->getObservations()-z_hat;
}

//! Utility to compute the sum of the squares of the residuals
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_SS_res(void)
{
        SS_res = eps_hat.squaredNorm();
}

//! Utility to compute the size of the model
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_s(void)
{
        s = this->model.getData()->getNumberofObservations();
}


//! Utility to compute the estimated variance of the error
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_sigma_hat_sq(void)
{
        sigma_hat_sq = SS_res/Real(dor);
}

//! Utility to compute a term useful in the first and second derivative of the GCV
/*
 * /sa {compute_fp(Real lambda), compute_fs(Real lambda)}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_aux(void)
{
        aux =  eps_hat.transpose()*dS_*this->model.getData()->getObservations();
}

// Updaters

//! Setting all the parameters which are recursively lambda dependent
/*
 * /remark{The order in which functions are invoked is essential for the consistency of the procedure}
 */
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_family(Real lambda)
{
        this->set_T_(lambda);
        this->set_V_();
        this->set_S_and_trS_();
        this->compute_z_hat();
        this->compute_eps_hat();
        this->compute_SS_res();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::zero_updater(Real lambda)
{
        this->update_parameters(lambda);  // Update all parameters depending on lambda
        std::get<0>(this->last_lambda) = lambda;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::first_updater(Real lambda)
{
        this->set_dS_and_trdS_();
        this->compute_aux();
        std::get<1>(this->last_lambda) = lambda;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::second_updater(Real lambda)
{
        this->set_ddS_and_trddS_();
        std::get<2>(this->last_lambda) = lambda;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::f_updater(Real lambda)
{
        if (lambda != std::get<0>(this->last_lambda))
                this->zero_updater(lambda);
        // else everything is up to date
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::fp_updater(Real lambda)
{
        if (lambda != std::get<0>(this->last_lambda))
        {
                // Update all parameters depending on lambda
                this->zero_updater(lambda);
                this->first_updater(lambda);
        }
        else if (lambda != std::get<1>(this->last_lambda))
        {
                // Update only the parameters for the first derivative
                this->first_updater(lambda);
        }
        // else everything is up to date
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::fs_updater(Real lambda)
{
        if (lambda != std::get<0>(this->last_lambda))
        {
                // Update all parameters depending on lambda
                this->zero_updater(lambda);
                this->first_updater(lambda);
                this->second_updater(lambda);
        }
        else if (lambda != std::get<1>(this->last_lambda))
        {
                // Update only the parameters for first and second derivative
                this->first_updater(lambda);
                this->second_updater(lambda);
        }
        else if (lambda != std::get<2>(this->last_lambda))
        {
                // Update only the parameters for the second derivative
                this->second_updater(lambda);
        }
        // else everything is up to date
}

// GCV function and derivatives

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_f(Real lambda)
{
        // GCV = s*(z-zhat)^t*(z-zhat)/(s-(q+trS))^2
        //     = SS_res*s/(dor^2)
        //     = sigma_hat_^2*s/dor
        this->f_updater(lambda);

        Real GCV_val = s*sigma_hat_sq/dor;    // Compute the value of the GCV and print it
	#ifdef R_VERSION_
                Rprintf("LAMBDA = %f\n",lambda);
		Rprintf("GCV = %f\n",GCV_val);
	#else
		std::cout << "GCV value =" << GCV_val << std::endl;
	#endif

	return GCV_val;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_fp(Real lambda)
{
        // dGCV(lambda)/dlambda = s * (d(1/dor(lambda)^2)/dlambda * SSres + dSSres(lambda)/dlambda * 1/dor^2)
        //                                    [1]                                        [2]
        // where [1] = 2/dor^3 * d(tr(S))/dlambda = 2/dor^3 * tr(Phi*T^{-1}*R*V)
        // and   [2] = 2*eps_hat^t*d(eps_hat)/dlambda = -2*eps^hat*dS
        // summing: 2*s/(dor^2) * (sigma_hat_^2*tr(dS/dlambda) - eps_hat*dS/dlambda*z)

        this->fp_updater(lambda);

        // Compute the value of the GCV first derivative and print it
	Real GCV_der_val = 2*s* (sigma_hat_sq*trdS_ - aux)/(dor*dor) ;
	#ifdef R_VERSION_
		Rprintf("GCV_derivative = %f\n", GCV_der_val);
	#else
		std::cout << "GCV_derivative value = " << GCV_der_val << std::endl;
	#endif

	return GCV_der_val;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real GCV_Family<InputHandler, Integrator, ORDER, mydim, ndim, 1>::compute_fs(Real lambda)
{
        this->fs_updater(lambda);

        VectorXr z = this->model.getData()->getObservations();
        VectorXr t = dS_*z;
        Real aux2;
        if (this->model.checkisRegression_())
                aux2 = t.transpose()*this->model.getQ_()*t;
        else
                aux2 = t.squaredNorm();
        Real aux3 = eps_hat.transpose()*ddS_*z;

        // Compute second derivative and print it
	Real GCV_sec_der_val =
                2*s*(trdS_*(3*sigma_hat_sq*trdS_-4*aux)/dor + sigma_hat_sq*trddS_ + aux2 - aux3)/(dor*dor);
	#ifdef R_VERSION_
		Rprintf("GCV_second_derivative = %f\n", GCV_sec_der_val);
	#else
		std::cout << "GCV_second_derivative value = " << GCV_sec_der_val << std::endl;
	#endif

	return GCV_sec_der_val;
}

//----------------------------------------------------------------------------//
// ** GCV_Exact **

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_dof()
{
	this->dof = this->trS_;

        if(this->model.checkisRegression_())
                this->dof += this->model.getData()->getCovariates().cols();
        Rprintf("trS:%f\n", this->trS_);
        Rprintf("degrees:%f\n", this->dof);
}

//! Utility to compute the degrees of freedom of the residuals
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_dor(Real lambda)
{
        this->dor = this->s-this->dof;

        if (this->dor < 0)   // Just in case of bad computation
        {
                #ifdef R_VERSION_
                              Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'. Value of 'lambda' that produces an error is: %d \n", lambda);
                #else
                              std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. " <<
                                            "This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'. " <<
                                            "Value of 'lambda' that produces an error is:" << lambda <<"\n";
                #endif
        }
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Exact<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_parameters(Real lambda)
{
        this->update_family(lambda);
        this->update_dof();
        this->update_dor(lambda);
        this->compute_sigma_hat_sq();
}

//----------------------------------------------------------------------------//
// ** GCV_Stochastic **

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_dof()
{
        /*
	UInt nnodes    = R_.rows();
        UInt nr        = model.getData()->getNrealizations();

	std::default_random_engine generator;
	std::bernoulli_distribution distribution(0.5);

	// Creation of the random matrix
	MatrixXr u(s, nr);
        for (int i = 0; i < s; ++i)
        {
		for (int j = 0; j < nr; ++j)
                {
			if (distribution(generator))
                        {
				u(i, j) = 1.0;
			}
			else
                        {
				u(i, j) = -1.0;
			}
		}
	}

        // STA CALCOLANDO LA SOLUZIONE DEL SISTEMA GENERALE? con dentro il lambda?
        // METODOLOGIA TERRIBILMENTE LENTA, DA USARE SOLO SE IN CASO DI NECESSITÀ
        // CI VORREBBE UNA STRUTURA DI HERROR HANDLING, IN CUI SOLO SE IL METODO ESATTO FALLISCE GLI SI F CALCOLARE QULLO STOCASTICO
	//_IN Negri la Q è l'identità
	// Define the first right hand side : | I  0 |^T * psi^T * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	b.topRows(nnodes) = psi_.transpose()* LeftMultiplybyQ(u);

	// Resolution of the system
	//MatrixXr x = system_solve(b);
	Eigen::SparseLU<SpMat> solver;
	solver.compute(_coeffmatrix);
	auto x = solver.solve(b);

	MatrixXr uTpsi = u.transpose()*psi_;
	VectorXr edf_vect(nr);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (regressionData_.getCovariates().rows() != 0) {
		q = regressionData_.getCovariates().cols();
	}
	// For any realization we calculate the degrees of freedom
	for (int i=0; i<nr; ++i) {
		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nr;
	_dof[output_index] = mean;
        */
}

//! Utility to compute the degrees of freedom of the residuals
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_dor(Real lambda)
{
        this->dor = this->s-this->dof;

        if (this->dor < 0)   // Just in case of bad computation
        {
                #ifdef R_VERSION_
                              Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'. Value of 'lambda' that produces an error is: %d \n", lambda);
                #else
                              std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. " <<
                                            "This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'. " <<
                                            "Value of 'lambda' that produces an error is:" << lambda <<"\n";
                #endif
        }
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void GCV_Stochastic<InputHandler, Integrator, ORDER, mydim, ndim, 1>::update_parameters(Real lambda)
{
        this->update_family(lambda);
        this->update_dof();
        this->update_dor(lambda);
        this->compute_sigma_hat_sq();
}

//----------------------------------------------------------------------------//

#endif
