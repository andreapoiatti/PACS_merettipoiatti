#ifndef __SOLUTION_BUILDERS_IMP_H__
#define __SOLUTION_BUILDERS_IMP_H__

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
SEXP Solution_Builders::build_solution_plain_regression(const VectorXr & solution, const output_Data & output, const MeshHandler<ORDER, mydim, ndim> & mesh , const InputHandler & regressionData )
{
        MatrixXv beta;
        if(regressionData.getCovariates()->rows()==0)
        {
                beta.resize(1,1);
                beta(0,0).resize(1);
                beta(0,0)(0) = 10e20;
        }
        else
        {
                beta = output.betas;
        }

        UInt code_string;
        if(output.content=="full_optimization")
        {
                code_string = 0;
        }
        else if(output.content=="full_dof_batch")
        {
                code_string = 1;
        }
        else
        {
                code_string = 2;
        }

        const MatrixXr & barycenters = regressionData.getBarycenters();
        const VectorXi & elementIds = regressionData.getElementIds();

        //Copy result in R memory
        SEXP result = NILSXP;
        result = PROTECT(Rf_allocVector(VECSXP, 22));

        SET_VECTOR_ELT(result, 0, Rf_allocVector(REALSXP, solution.size()));
        Real *rans = REAL(VECTOR_ELT(result, 0));
        for(UInt j = 0; j < solution.size(); j++)  //[TO DO ] //sono le f_hat e g_hat, si potrebbe rimuovere, cambiando la chiamata da R in  smooth.FEM.basis
        {
               rans[j] = solution[j];
        }

        UInt size_z=output.z_hat.size();
        SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, size_z));
        rans = REAL(VECTOR_ELT(result, 1));
        for(UInt j = 0; j < size_z; j++)
        {
               rans[j] = output.z_hat[j];
        }

        //Rprintf("Hey doc,  %f %f %f\n", output.z_hat[0], output.z_hat[1], output.z_hat[3]);
        SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 2));
        rans[0] = output.rmse;

        SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, 1));
        rans= REAL(VECTOR_ELT(result, 3));
        rans[0] = output.sigma_hat_sq;

        SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 4));
        rans[0] = output.lambda_sol;

        SET_VECTOR_ELT(result, 5, Rf_allocVector(INTSXP, 1));
        UInt *rans1 = INTEGER(VECTOR_ELT(result, 5));
        rans1[0] = output.lambda_pos;

        SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 6));
        rans[0] = output.GCV_opt;

        SET_VECTOR_ELT(result, 7, Rf_allocVector(INTSXP, 1));
        UInt * rans2 = INTEGER(VECTOR_ELT(result, 7));
        rans2[0] = output.n_it;

        SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 8));
        rans[0] = output.termination;

        SET_VECTOR_ELT(result, 9, Rf_allocVector(INTSXP, 1));
        UInt *rans3 = INTEGER(VECTOR_ELT(result, 9));
        rans3[0] = code_string;

        UInt size_dof=output.dof.size();
        SET_VECTOR_ELT(result, 10, Rf_allocVector(REALSXP, size_dof));
        rans = REAL(VECTOR_ELT(result, 10));
        for(UInt j = 0; j < size_dof; j++)
        {
               rans[j] = output.dof[j];
        }

        UInt size_lambda=output.lambda_vec.size();
        SET_VECTOR_ELT(result, 11, Rf_allocVector(REALSXP, size_lambda));
        rans = REAL(VECTOR_ELT(result, 11));

        for(UInt j = 0; j < size_lambda; j++)
        {
               rans[j] = output.lambda_vec[j];
        }

        UInt size_vec=output.GCV_evals.size();
        SET_VECTOR_ELT(result, 12, Rf_allocVector(REALSXP, size_vec));
        rans = REAL(VECTOR_ELT(result, 12));
        for(UInt j = 0; j < size_vec; j++)
        {
               rans[j] = output.GCV_evals[j];
        }

        SET_VECTOR_ELT(result, 13, Rf_allocVector(REALSXP, 1));
        rans = REAL(VECTOR_ELT(result, 13));
        rans[0] = output.time_partial;

        //! Copy betas
        SET_VECTOR_ELT(result, 14, Rf_allocMatrix(REALSXP, beta(0).size(), beta.size()));
        Real *rans4 = REAL(VECTOR_ELT(result, 14));
        for(UInt j = 0; j < beta.size(); j++)
        {
                for(UInt i = 0; i < beta(0).size(); i++)
                        rans4[i + beta(0).size()*j] = beta(j)(i);
        }

        //SEND TREE INFORMATION TO R
        SET_VECTOR_ELT(result, 15, Rf_allocVector(INTSXP, 1)); //tree_header information
        int *rans5 = INTEGER(VECTOR_ELT(result, 15));
        rans5[0] = mesh.getTree().gettreeheader().gettreelev();

        SET_VECTOR_ELT(result, 16, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain origin
        Real *rans6 = REAL(VECTOR_ELT(result, 16));
        for(UInt i = 0; i < ndim*2; i++)
        rans6[i] = mesh.getTree().gettreeheader().domainorig(i);

        SET_VECTOR_ELT(result, 17, Rf_allocVector(REALSXP, ndim*2)); //tree_header domain scale
        Real *rans7 = REAL(VECTOR_ELT(result, 17));
        for(UInt i = 0; i < ndim*2; i++)
        rans7[i] = mesh.getTree().gettreeheader().domainscal(i);


        UInt num_tree_nodes = mesh.num_elements()+1; //Be careful! This is not equal to number of elements
        SET_VECTOR_ELT(result, 18, Rf_allocMatrix(INTSXP, num_tree_nodes, 3)); //treenode information
        int *rans8 = INTEGER(VECTOR_ELT(result, 18));
        for(UInt i = 0; i < num_tree_nodes; i++)
                rans8[i] = mesh.getTree().gettreenode(i).getid();

        for(UInt i = 0; i < num_tree_nodes; i++)
                rans8[i + num_tree_nodes*1] = mesh.getTree().gettreenode(i).getchild(0);

        for(UInt i = 0; i < num_tree_nodes; i++)
                rans8[i + num_tree_nodes*2] = mesh.getTree().gettreenode(i).getchild(1);

        SET_VECTOR_ELT(result, 19, Rf_allocMatrix(REALSXP, num_tree_nodes, ndim*2)); //treenode box coordinate
        Real *rans9 = REAL(VECTOR_ELT(result, 19));
        for(UInt j = 0; j < ndim*2; j++)
        {
                for(UInt i = 0; i < num_tree_nodes; i++)
                        rans9[i + num_tree_nodes*j] = mesh.getTree().gettreenode(i).getbox().get()[j];
        }

        //SEND BARYCENTER INFORMATION TO R
        SET_VECTOR_ELT(result, 20, Rf_allocVector(INTSXP, elementIds.rows())); //element id of the locations point (vector)
        int *rans10 = INTEGER(VECTOR_ELT(result, 20));
        for(UInt i = 0; i < elementIds.rows(); i++)
        rans10[i] = elementIds(i);

        SET_VECTOR_ELT(result, 21, Rf_allocMatrix(REALSXP, barycenters.rows(), barycenters.cols())); //barycenter information (matrix)
        Real *rans11 = REAL(VECTOR_ELT(result, 21));
        for(UInt j = 0; j < barycenters.cols(); j++)
        {
        for(UInt i = 0; i < barycenters.rows(); i++)
                rans11[i + barycenters.rows()*j] = barycenters(i,j);
        }


        UNPROTECT(1);

        return(result);
}

#endif