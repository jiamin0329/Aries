/*!
 * \file matrix_structure.inl
 * \brief In-Line subroutines of the <i>matrix_structure.hpp</i> file.
 */

#ifndef ARIES_MATH_MATRIX_INLINE
#define ARIES_MATH_MATRIX_INLINE

namespace ARIES
{
    namespace MATH
    {

        inline void MATH_Matrix::SetValZero(void)
        {
            for (unsigned long index = 0; index < nnz*nVar*nEqn; index++)
                matrix[index] = 0.0;
        }

        inline MATH_Matrix_MatrixVectorProduct::MATH_Matrix_MatrixVectorProduct(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref)
        {
            sparse_matrix = &matrix_ref;
            geometry = geometry_ref;
            config = config_ref;
        }

        inline void MATH_Matrix_MatrixVectorProduct::operator()(const MATH_Vector & u, MATH_Vector & v) const
        {
            if (sparse_matrix == NULL) 
            {
                std::cerr << "MATH_MatrixVectorProduct::operator()(const MATH_Vector &, MATH_Vector &): " << std::endl;
                std::cerr << "pointer to sparse matrix is NULL." << std::endl;
                throw(-1);
            }
            sparse_matrix->MatrixVectorProduct(u, v, geometry, config);
        }

        inline MATH_JacobiPreconditioner::MATH_JacobiPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref) 
        {
            sparse_matrix = &matrix_ref;
            geometry = geometry_ref;
            config = config_ref;
        }

        inline void MATH_JacobiPreconditioner::operator()(const MATH_Vector & u, MATH_Vector & v) const 
        {
            if (sparse_matrix == NULL) 
            {
                std::cerr << "CJacobiPreconditioner::operator()(const MATH_Vector &, MATH_Vector &): " << std::endl;
                std::cerr << "pointer to sparse matrix is NULL." << std::endl;
                throw(-1);
            }
            sparse_matrix->ComputeJacobiPreconditioner(u, v, geometry, config);
        }

        inline MATH_ILUPreconditioner::MATH_ILUPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref) 
        {
            sparse_matrix = &matrix_ref;
            geometry = geometry_ref;
            config = config_ref;
        }

        inline void MATH_ILUPreconditioner::operator()(const MATH_Vector & u, MATH_Vector & v) const 
        {
            if (sparse_matrix == NULL) 
            {
                std::cerr << "CILUPreconditioner::operator()(const MATH_Vector &, MATH_Vector &): " << std::endl;
                std::cerr << "pointer to sparse matrix is NULL." << std::endl;
                throw(-1);
            }
            sparse_matrix->ComputeILUPreconditioner(u, v, geometry, config);
        }

        inline MATH_LUSGSPreconditioner::MATH_LUSGSPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref)
        {
            sparse_matrix = &matrix_ref;
            geometry = geometry_ref;
            config = config_ref;
        }

        inline void MATH_LUSGSPreconditioner::operator()(const MATH_Vector & u, MATH_Vector & v) const
        {
            if (sparse_matrix == NULL) 
            {
                std::cerr << "CLU_SGSPreconditioner::operator()(const MATH_Vector &, MATH_Vector &): " << std::endl;
                std::cerr << "pointer to sparse matrix is NULL." << std::endl;
                throw(-1);
            }
            sparse_matrix->ComputeLU_SGSPreconditioner(u, v, geometry, config);
        }

        inline MATH_LineletPreconditioner::MATH_LineletPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref)
        {
            sparse_matrix = &matrix_ref;
            geometry = geometry_ref;
            config = config_ref;
        }

        inline void MATH_LineletPreconditioner::operator()(const MATH_Vector & u, MATH_Vector & v) const
        {
            if (sparse_matrix == NULL) 
            {
                std::cerr << "CLineletPreconditioner::operator()(const MATH_Vector &, MATH_Vector &): " << std::endl;
                std::cerr << "pointer to sparse matrix is NULL." << std::endl;
                throw(-1);
            }
            sparse_matrix->ComputeLineletPreconditioner(u, v, geometry, config);
        }
    }
}

#endif