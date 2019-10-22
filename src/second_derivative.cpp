
//----------------------------------------------------------------------------//

//_computation of GCV_second_derivative 
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>::computeGCV_second_derivative(UInt output_index)
{ //_da modificare il fatto che _dof sia un vettore, serve un solo valore, non serve output index!!
	UInt s;
	//_UInt q=regressionData_.getCovariates().cols(); //_serve se si vuole usare l'articolo di stuHuntersangalli, è già implicito nel degrees of freedom
	VectorXr z;
	s= regressionData_.getNumberofObservations(); //_così ho anche il caso in uci ho meno locations dei nodi (pur coincidenti)
	//if(regressionData_.isLocationsByNodes())
	//{
		//s= this->mesh_.num_nodes(); //_vuol dire che il numero di locations (e quindi il numero di osservazioni, coincide col numero di nodi e le posizioni sono esattamente quelle dei nodi)
    //_s è la n dell'articolo stuHuntersangalli pdf pag.12, numero locations, dove ho le osservazioni
		//z=VectorXr::Zero(s);
		//for(auto i=0;i<regressionData_.getObservationsIndices().size();i++)
		//	z(regressionData_.getObservationsIndices()[i])=regressionData_.getObservationData()[i]; //_mette i valori nei nodi in cui si pongono le locations (nodi coincidono con le location, ma le locations possono avere ordini di numerazione diversi!!
	//}
	// else
	// {

	z=regressionData_.getObservations();
	  //}
        MatrixXr I=MatrixXr::Identity(s,s);
	//NB _questo caso è ok se i nodi non coincidono con le location, è ridondante (migliorabile!!) se coindicono, perchè psi è l'identitò (si può fare come nel calcolo dei deg of freedom per essere più efficiente)
	Real norm_squared=(z-z_hat_).transpose()*(z-z_hat_);
        Real trace_=0.0, trace_der=0.0; //_trace of derivative and trace of second derivative
	MatrixXr dS_(s,s),ddS_(s,s); //S derivative dlambda e S second derivative
	MatrixXr aux; //matrix auxiliary for second derivative
       //al più si può migliorare evitando l'ultimo prodotto per psi e ragionando con la k
	Eigen::LDLT<MatrixXr> Dsolver( SS_ );
	dS_=-psi_*Dsolver.solve( MatrixXr(R_*V_) ); //_se dà errore, provare MatrixXr(R_*V_), per ricreare al più la matrice
	//_d_S=-psi*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*psi^T*Q
        aux=Dsolver.solve(R_); //(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1
	ddS_=2*psi_*aux*aux*V_; //2*psi*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*R1*R0^(-1)*R1*(psi^T*Q*psi+lambda*R1*R0^-1*R1)^(-1)*psi^T*Q
	//NB possibile da stoccare dS, anche se deve essere chiamata prima la derivata per stoccarla->per ora  la tengo così

        //for (UInt i=0; i<mesh_.num_nodes(); i++) //_anche se sarebbe più corretto il numero di osservazioni, è nxn

	for (UInt i=0; i<s; i++)
		trace_+=dS_(i,i); //_tr(dS/dlambda)=d(tr(S))/dlambda

	for (UInt i=0; i<s; i++)
		trace_der+=ddS_(i,i); //_tr(d^2S/dlambda^2)=d^2(tr(S))/dlambda^2

	if(s-_dof[output_index]<0){ //_dof_ non servirà, sarà un valore unico!
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->regressionData_.getLambda()[output_index]);
 			#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->regressionData_.getLambda()[output_index] <<"\n";
			#endif
			}
	//Real stderror=norm_squared/(s-_dof[output_index]); //così è ancora fatta sul vettore
       //uso la proprietà di simmetri e idempotenza di Q, ho Q^T*Q=Q*Q=Q
        Real c_=1/((s-_dof[output_index]) * (s-_dof[output_index])*(s-_dof[output_index]));
	Real first_=2*c_*s*trace_;
	VectorXr r_=z.transpose()*(-dS_.transpose())*LeftMultiplybyQ(MatrixXr(I-S_))*z+z.transpose()*(I-S_.transpose())*LeftMultiplybyQ(MatrixXr(-dS_))*z; //prodotto per matrici può essere un vettore, non un Real->estraggo la componente 0
        Real der_aux=c_*c_*norm_squared*(2*s)*(trace_der/c_+3*trace_*trace_*(s-_dof[output_index])*(s-_dof[output_index]));
        VectorXr second_=-z.transpose()*ddS_.transpose()*Q_*(I-S_)*z+z.transpose()*dS_.transpose()*Q_*dS_*z+z.transpose()*dS_.transpose()*Q_*dS_*z-z.transpose()*(I-S_).transpose()*Q_*ddS_*z;
	VectorXr GCV_sec_der=2*first_*r_+second_*s*c_*(s-_dof[output_index]);
	Real GCV_sec_der_val=GCV_sec_der[0]+der_aux;
	#ifdef R_VERSION_
		Rprintf("GCV_second_derivative=%f\n",GCV_sec_der_val);
	#else
		std::cout << "GCV_second_derivative value="<<GCV_sec_der_val<<std::endl;
	#endif

	return GCV_sec_der_val; //_Calcolo della derivata seconda della GCV

}


//----------------------------------------------------------------------------//

