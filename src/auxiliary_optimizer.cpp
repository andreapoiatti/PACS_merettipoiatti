#include"auxiliary_optimizer.h"

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
// to avoid a direct sum with the full matrix, taking only the
// non-null components on its diagonal and subsequently summing
// them direcly to the second block.

void AuxiliaryOptimizer::set_T_nW_a(SpMat & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp)
{
        // Avoid using Q
        T += (*psi_tp)*(*Ap).asDiagonal()*(*psip);
}

void AuxiliaryOptimizer::set_T_W_a(SpMat & T, const VectorXr * Ap, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp)
{
        // Full model, no simplification allowed
        SpMat temp = ((*Qp)*(*psip)).sparseView();
        T += (*psi_tp)*(*Ap).asDiagonal()*temp;
}

void AuxiliaryOptimizer::set_T_ln_nW_ptw(SpMat & T, const std::vector<UInt> * kp, UInt s)
{
        // T = Psi^t*Psi == Indicator(k[i],k[i])
        for (UInt i = 0; i < s ; i++)
                T.coeffRef((*kp)[i], (*kp)[i]) += 1;
}

void AuxiliaryOptimizer::set_T_ln_W_ptw(SpMat & T, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt s)
{
        // T = Psi^t*Q*Psi == q_ij*Indicator(k[i],k[j])
        for (UInt i = 0; i < s ; i++)
                for (int j = 0; j < s; j++)
                        T.coeffRef((*kp)[i], (*kp)[j]) += (*Qp).coeff(i, j);
}

void AuxiliaryOptimizer::set_T_lnn_nW_ptw(SpMat & T, const SpMat * psip, const SpMat * psi_tp)
{
        // Avoid using Q
        T += (*psi_tp)*(*psip);
}

void AuxiliaryOptimizer::set_T_lnn_W_ptw(SpMat & T, const SpMat * psip, const SpMat * psi_tp, const MatrixXr * Qp)
{
        // Full model, no simplification allowed
        SpMat temp = ((*Qp)*(*psip)).sparseView();
        T += (*psi_tp)*temp;
}

// THEORETICAL REMARK:
// Since Psi is a rectangular permutation matrix, if function
// k: loctions -> nodes s.t. Psi = Indicator(i,k[i]) then
// Psi^t*Q   == Indicator(k[i],j)*q_ij

// IMPLEMENTATION OF THE REMARK:
// the number of non-null entries of E is at most s^2,
// we reserve a vector containing such entries and
// we set the final matrix from these triplets

void AuxiliaryOptimizer::set_E_ln_W_ptw(SpMat & E, const std::vector<UInt> * kp, const MatrixXr * Qp, UInt nr, UInt s)
{
        E.resize(nr, s);
        std::vector<coeff> vec;
        vec.reserve(s*s);

        for (UInt i = 0; i < s ; i++)
                for (UInt j = 0; j < s; j++)
                        vec.push_back(coeff((*kp)[i], j, (*Qp).coeff(i, j)));

        E.setFromTriplets(vec.begin(), vec.end());
}

void AuxiliaryOptimizer::set_E_lnn_W_ptw(SpMat & E, const SpMat * psi_tp, const MatrixXr * Qp)
{
        E = ((*psi_tp)*(*Qp)).sparseView();
}

void AuxiliaryOptimizer::set_E_W_a(SpMat & E, const SpMat * psi_tp, const MatrixXr * Qp, const VectorXr * Ap)
{
        E = ((*psi_tp)*(*Ap).asDiagonal()*(*Qp)).sparseView();
}

void AuxiliaryOptimizer::set_E_nW_a(SpMat & E, const SpMat * psi_tp, const VectorXr * Ap)
{
        E = ((*psi_tp)*(*Ap).asDiagonal());
}


void AuxiliaryOptimizer::set_z_hat_W(VectorXr & z_hat, const MatrixXr * Hp, const MatrixXr * Qp, const SpMat & S, const VectorXr * zp)
{
        z_hat = ((*Hp)+(*Qp)*S)*(*zp);
}

void AuxiliaryOptimizer::set_z_hat_nW(VectorXr & z_hat, const SpMat & S, const VectorXr * zp)
{
        z_hat = S*(*zp);
}
