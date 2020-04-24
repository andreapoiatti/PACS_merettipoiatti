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

template<typename InputCarrier>
const output_Data & GCV_Family<InputCarrier, 1>::get_output(std::pair<Real,UInt> p, const timespec & T)
{
        this->output.lambda_sol         = p.first;
        this->output.n_it               = p.second;
        this->output.z_hat              = this->z_hat;
        this->output.SS_res             = this->SS_res;
        this->output.sigma_hat_sq       = this->sigma_hat_sq;
        this->output.time_partial       = T.tv_sec + 1e-9*T.tv_nsec;

        return this->output;
}


template<typename InputCarrier>
const output_Data & GCV_Family<InputCarrier, 1>::get_output_partial(void)
{
        this->output.z_hat              = this->z_hat;
        this->output.SS_res             = this->SS_res;
        this->output.sigma_hat_sq       = this->sigma_hat_sq;

        return this->output;
}


template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_R_(void)
{
        const SpMat * R1p_= this->the_carrier.get_R1p();          // Get the value of matrix R1
        Sparse_LU solver;	                                 // define a factorized empty sparse Cholesky solver
        solver.compute(*(this->the_carrier.get_R0p()));		 // apply it to R0 to simplify the inverse
        R_ = (*R1p_).transpose()*solver.solve(*R1p_);            // R == _R1^t*R0^{-1}*R1
        R_.makeCompressed();                                     // Compress the matrix to speed up the next operations
}

//! Method to set the value of member SpMat T_
/*
 * /remark {T = D + \lambda * R where D is the top-left block of the matrix DMat}
 * /pre {set_R_ must be called before set_T_, the matrix D_ [DataMatrix] must be constructed in the model, s must be defined}
 * /sa {set_R_(void), getDataMatrix(SpMat & DMat)}
 */
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_T_(Real lambda)
{
        T_ = lambda*R_;
        UInt ret =  AuxiliaryOptimizer::universal_T_setter<InputCarrier>(T_, this->the_carrier);
        T_.makeCompressed();    // Compressing the matrix for further computation
}

//! Method to set the value of member SpMat V_
/*
 * /remark {V = T^{-1}*Psi^t*Q }
 * /pre {set_T_ must be called before set_V_}
 * /sa {set_T_(Real lambda)}
 */
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_V_(void)
{
        Sparse_LU solver;	        // define a factorized empty sparse LU solver
        solver.compute(T_);             // apply it to T to simplify the inverse

        if(!this->the_carrier.is_areal() && !this->the_carrier.has_W())
        {
                // Q == I
                const SpMat * psi_tp = this->the_carrier.get_psi_tp();
                V_ = solver.solve(*psi_tp);      // find the value of V = T^{-1}*Psi^t
        }
        else
        {
                SpMat E_;                // Declare an empty auxiliary matrix
                UInt ret =  AuxiliaryOptimizer::universal_E_setter<InputCarrier>(E_, this->the_carrier);
                E_.makeCompressed();                      // compress the matrix to save space
                V_ = solver.solve(E_);          // find the value of V = T^{-1}*E
        }

        V_.makeCompressed();            // compress the matrix to save space
}

//! Method to set the value of member SpMat S_ and its trace trS_
/*
 * /remark {S = Psi*V }
 * /pre {set_V_ must be called before set_S_}
 * /sa {set_V_(void)}
 */
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_S_and_trS_(void)
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
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_dS_and_trdS_(void)
{
        // dS_ = -Psi*(Psi^t*Q*Psi+lambda*R1^t*R0^-1*R1)^{-1}*R1^t*R0^{-1}*R1*(Psi^t*Q*Psi+lambda*R1^t*R0^{-1}*R1)^{-1}*Psi^t*Q
        //     = -Psi*T^{-1}*R*V
        //     =  Psi*(-K*V)
        //    :=  Psi*F
        Sparse_Cholesky solver;	                // define a factorized empty sparse Cholesky solver [[NOT BETTER TO SOLVE THE ONE WE USE FOR T?]]
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
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_ddS_and_trddS_(void)
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
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::LeftMultiplybyPsiAndTrace(Real & trace, SpMat & ret, const SpMat & mat)
{
        if (this->the_carrier.loc_are_nodes())
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

                const std::vector<UInt> * kp = this->the_carrier.get_obs_indicesp();
                for (UInt i = 0; i < s; i++)
                        for (UInt j = 0; j < s; j++)
                        {
                                if (i == j)
                                {
                                        Real v = mat.coeff((*kp)[i], j);
                                        vec.push_back(coeff(i, j, v));
                                        trace += v;
                                }
                                else
                                {
                                        vec.push_back(coeff(i, j, mat.coeff((*kp)[i], j)));
                                }
                        }


                ret.setFromTriplets(vec.begin(), vec.end());
        }
        else
        {
                // Psi is full
                ret = (*this->the_carrier.get_psip())*mat;
                for (int i = 0; i < s; ++i)
                        trace += ret.coeff(i, i);
        }
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::set_US_(void)
{
        // Creation of the random matrix
        std::default_random_engine generator;
	std::bernoulli_distribution distribution(0.5);

        Rprintf("Set_US");

        UInt nr = this->the_carrier.get_opt_data()->get_nrealizations_();
        this->US_ = MatrixXr::Zero(this->s, nr);

        for (UInt i=0; i<this->s; ++i)
                for (UInt j=0; j<nr; ++j)
                {
                        if (distribution(generator))
                        {
                                this->US_.coeffRef(i, j) = 1.0;
                        }
                        else
                        {
                                this->US_.coeffRef(i, j) = -1.0;
                        }
                }

        this->us = true;
}

//! Utility to compute the predicted values in the locations
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_z_hat(void)
{
        // [[ TODO NOT IMPLEMENTED PRESENCE OF FORCING TERM]]
        const VectorXr * zp = this->the_carrier.get_zp();
        if(this->the_carrier.has_W())
        {
                const MatrixXr * Hp = this->the_carrier.get_Hp();
                const MatrixXr * Qp = this->the_carrier.get_Qp();
                AuxiliaryOptimizer::set_z_hat_W(z_hat, Hp, Qp, S_, zp);
        }
        else
        {
                AuxiliaryOptimizer::set_z_hat_nW(z_hat, S_, zp);
        }

        // Debugging purpose
        // Rprintf("z_hat: \n");
        // for(int i = 0; i < 20; i++)
        //        Rprintf("%f,", this->z_hat[i]);
	// Rprintf("\n");
}

//! Utility to compute the predicted residuals in the locations
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_eps_hat(void)
{
        eps_hat = (*this->the_carrier.get_zp())-z_hat;

        // Debugging purpose
        //Rprintf("Eps_hat \n");
        //for(UInt i = 0; i < this->s; i++)
        //        Rprintf("%f, ", this->eps_hat[i]);
        //Rprintf("\n");
}

//! Utility to compute the sum of the squares of the residuals
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_SS_res(void)
{
        SS_res = eps_hat.squaredNorm();

        // Debugging purpose
        //Rprintf("SS_res  = %f\n", this->SS_res);
}

//! Utility to compute the size of the model
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_s(void)
{
        s = this->the_carrier.get_n_obs();
}


//! Utility to compute the estimated variance of the error
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_sigma_hat_sq(void)
{
        sigma_hat_sq = SS_res/Real(dor);

        // Debugging purpose
        //Rprintf("sigma_hat_sq = %f\n", this->sigma_hat_sq);
}

//! Utility to compute a term useful in the first and second derivative of the GCV
/*
 * /sa {compute_fp(Real lambda), compute_fs(Real lambda)}
 */
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::compute_aux(void)
{
        aux =  eps_hat.transpose()*dS_*(*this->the_carrier.get_zp());
}

// Updaters

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::update_family_p1(Real lambda)
{
        this->set_T_(lambda);
        this->set_V_();
        this->set_S_and_trS_();
        this->compute_z_hat();
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::update_family_p2(void)
{
        this->compute_eps_hat();
        this->compute_SS_res();
}

//! Setting all the parameters which are recursively lambda dependent
/*
 * /remark{The order in which functions are invoked is essential for the consistency of the procedure}
 */
template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::update_family(Real lambda)
{
        this->update_family_p1(lambda);
        this->update_family_p2();
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::zero_updater(Real lambda)
{
        this->update_parameters(lambda);  // Update all parameters depending on lambda
        std::get<0>(this->last_lambda) = lambda;
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::first_updater(Real lambda)
{
        this->set_dS_and_trdS_();
        this->compute_aux();
        std::get<1>(this->last_lambda) = lambda;
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::second_updater(Real lambda)
{
        this->set_ddS_and_trddS_();
        std::get<2>(this->last_lambda) = lambda;
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::f_updater(Real lambda)
{
        if (lambda != std::get<0>(this->last_lambda))
                this->zero_updater(lambda);
        // else everything is up to date
}

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::fp_updater(Real lambda)
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

template<typename InputCarrier>
void GCV_Family<InputCarrier, 1>::fs_updater(Real lambda)
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

template<typename InputCarrier>
Real GCV_Family<InputCarrier, 1>::compute_f(Real lambda)
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

template<typename InputCarrier>
Real GCV_Family<InputCarrier, 1>::compute_fp(Real lambda)
{
        // dGCV(lambda)/dlambda = s * (d(1/dor(lambda)^2)/dlambda * SSres + dSSres(lambda)/dlambda * 1/dor^2)
        //                                    [1]                                        [2]
        // where [1] = 2/dor^3 * d(tr(S))/dlambda = 2/dor^3 * tr(Phi*T^{-1}*R*V)
        // and   [2] = 2*eps_hat^t*d(eps_hat)/dlambda = -2*eps^hat*dS
        // summing: 2*s/(dor^2) * (sigma_hat_^2*tr(dS/dlambda) - eps_hat*dS/dlambda*z)

        this->fp_updater(lambda);

        // Compute the value of the GCV first derivative and print it
	Real GCV_der_val = 2*s* (sigma_hat_sq*trdS_ - aux)/(dor*dor);

	#ifdef R_VERSION_
		Rprintf("GCV_derivative = %f\n", GCV_der_val);
	#else
		std::cout << "GCV_derivative value = " << GCV_der_val << std::endl;
	#endif

	return GCV_der_val;
}

template<typename InputCarrier>
Real GCV_Family<InputCarrier, 1>::compute_fs(Real lambda)
{
        this->fs_updater(lambda);

        const VectorXr * zp = this->the_carrier.get_zp();
        VectorXr t = dS_*(*zp);
        Real aux2;
        if (this->the_carrier.has_W())
                aux2 = t.transpose()*(*this->the_carrier.get_Qp())*t;
        else
                aux2 = t.squaredNorm();
        Real aux3 = eps_hat.transpose()*ddS_*(*zp);

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

template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_dof(Real lambda)
{
	this->dof = this->trS_;

        if(this->the_carrier.has_W())
                this->dof += (*this->the_carrier.get_Wp()).cols();

        // Debugging purpose
        // Rprintf("DOF: %f\n", this->dof);
}

//! Utility to compute the degrees of freedom of the residuals
template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_dor(Real lambda)
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

template<typename InputCarrier>
void GCV_Exact<InputCarrier, 1>::update_parameters(Real lambda)
{
        this->update_family(lambda);
        this->update_dof(lambda);
        this->update_dor(lambda);
        this->compute_sigma_hat_sq();
}

//----------------------------------------------------------------------------//
// ** GCV_Stochastic **

template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::update_dof(Real lambda)
{
	UInt nnodes    = this->R_.rows();
        UInt nr        = this->the_carrier.get_opt_data()->get_nrealizations_();

        if(this->us == false)
        {
                this->set_US_();
        }

	// Define the first right hand side : | I  0 |^T * psi^T * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes, this->US_.cols());
        UInt ret = AuxiliaryOptimizer::universal_b_setter(b, this->the_carrier, this->US_, nnodes);


	// Resolution of the system [[TO BE OPTIMIZED in PHI]]
        SpMat R1_lambda = (-lambda)*(*this->the_carrier.get_R1p());
	SpMat R0_lambda = (-lambda)*(*this->the_carrier.get_R0p());

	this->the_carrier.get_tracep()->buildMatrixNoCov((*this->the_carrier.get_psip()), R1_lambda, R0_lambda);

	// Factorize the system
	this->the_carrier.get_tracep()->system_factorize();

	// Solve the system
    	auto x = this->the_carrier.get_tracep()->system_solve(b);

	MatrixXr USTpsi = this->US_.transpose()*(*this->the_carrier.get_psip());
	VectorXr edf_vect(nr);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (this->the_carrier.has_W())
        {
		q = this->the_carrier.get_Wp()->cols();
	}
	// For any realization we calculate the degrees of freedom
	for (UInt i = 0; i < nr; ++i)
        {
		edf_vect(i) = USTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	this->dof = edf_vect.sum()/nr;

        // Debugging purpose
        //Rprintf("DOF:%f\n", this->dof);
}

//! Utility to compute the degrees of freedom of the residuals
template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::update_dor(Real lambda)
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

template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::compute_z_hat(Real lambda)
{
        this->the_carrier.get_tracep()->apply(lambda);

        UInt nnodes    = this->R_.rows();
        VectorXr f_hat = this->the_carrier.get_tracep()->getSolution().head(nnodes);


        //[[todo fixing ]]
        if (this->the_carrier.has_W())
        {
                this->z_hat = (*this->the_carrier.get_Hp())*(*this->the_carrier.get_zp()) + (*this->the_carrier.get_Qp())*(*this->the_carrier.get_psip())*f_hat;
        }
        else
        {
                this->z_hat = (*this->the_carrier.get_psip())*f_hat;
        }
}

template<typename InputCarrier>
void GCV_Stochastic<InputCarrier, 1>::update_parameters(Real lambda)
{
        this->compute_z_hat(lambda);
        this->update_family_p2();
        this->update_dof(lambda);
        this->update_dor(lambda);
        this->compute_sigma_hat_sq();
}

//----------------------------------------------------------------------------//

#endif
