#include "../Include/Auxiliary_Optimizer.h"

// THEORETICAL REMARK:
// Since Psi is a rectangular permutation matrix, if function
// k: loctions -> nodes s.t. Psi = Indicator(i,k[i]) then
// Psi^t*Q   == Indicator(k[i],j)*q_ij

// IMPLEMENTATION OF THE REMARK:
// the number of non-null entries of E is at most s^2,
// we reserve a vector containing such entries and
// we set the final matrix from these triplets

void AuxiliaryOptimizer::set_E_ln_W_ptw(MatrixXr & E, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt nr, UInt s)
{
        E = MatrixXr::Zero(nr, s);

        for (UInt i = 0; i < s ; i++)
                for (int j = 0; j < s; j++)
                        E.coeffRef((*kp)[i], j) += (*Qp).coeff(i, j);
}

void AuxiliaryOptimizer::set_E_lnn_W_ptw(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp)
{
        E = ((*psi_tp)*(*Qp));
}

void AuxiliaryOptimizer::set_E_W_a(MatrixXr & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap)
{
        E = ((*psi_tp)*(*Ap).asDiagonal()*(*Qp));
}

void AuxiliaryOptimizer::set_E_nW_a(MatrixXr & E, const SpMat * psi_tp, const VectorXr * Ap)
{
        E = ((*psi_tp)*(*Ap).asDiagonal());
}


void AuxiliaryOptimizer::set_z_hat_nW(VectorXr & z_hat, const MatrixXr & S, const VectorXr * zp)
{
        z_hat = S*(*zp);
}
