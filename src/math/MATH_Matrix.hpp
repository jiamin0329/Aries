/*********************************************************************************
 *                         ARIES Copyright(C), 2015.
 *
 *  \file    COMM_Vector.hpp
 *  \brief   Main class for defining sparse matrices-by-blocks
 *           with compressed row format.
 *********************************************************************************
 *      Date        Author        Version                   Reason
 *    6/11/2015    Jiamin XU        1.0                  Initial release
 *
 *
 */

#ifndef ARIES_MATH_MATRIX_HPP
#define ARIES_MATH_MATRIX_HPP

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

//ARIES headers
#include "../Common/TBOX_Config.hpp"
#include "../Geometry/GEOM_Geometry.hpp"
#include "MATH_Vector.hpp"

namespace ARIES
{
    namespace MATH
    {
        class MATH_Matrix
        {
        public:
            /*!
             * \brief Constructor of the class.
             */
            MATH_Matrix(void);

            /*!
             * \brief Destructor of the class.
             */
            ~MATH_Matrix(void);

            /*!
               * \brief Initializes space matrix system.
               * \param[in] nVar - Number of variables.
               * \param[in] nEqn - Number of equations.
               * \param[in] geometry - Geometrical definition of the problem.
               * \param[in] config - Definition of the particular problem.
               */
            void Initialize(unsigned long nPoint, unsigned long nPointDomain, unsigned short nVar, unsigned short nEqn, bool EdgeConnect, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
               * \brief Assings values to the sparse-matrix structure.
               * \param[in] val_nPoint - Number of points in the nPoint x nPoint block structure
               * \param[in] val_nVar - Number of nVar x nVar variables in each subblock of the matrix-by-block structure.
               * \param[in] val_nEq - Number of nEqn x nVar variables in each subblock of the matrix-by-block structure.
               * \param[in] val_row_ptr - Pointers to the first element in each row.
               * \param[in] val_col_ind - Column index for each of the elements in val().
               * \param[in] val_nnz - Number of possible nonzero entries in the matrix.
               * \param[in] config - Definition of the particular problem.
               */
            void SetIndexes(unsigned long val_nPoint, unsigned long val_nPointDomain, unsigned short val_nVar, unsigned short val_nEq, unsigned long* val_row_ptr, unsigned long* val_col_ind, unsigned long val_nnz, TBOX::TBOX_Config *config);

            /*!
             * \brief Sets to zero all the entries of the sparse matrix.
             */
            void SetValZero(void);

            /*!
               * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               */
            double *GetBlock(unsigned long block_i, unsigned long block_j);

            /*!
               * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               */
            double GetBlock(unsigned long block_i, unsigned long block_j, unsigned short iVar, unsigned short jVar);

            /*!
               * \brief Set the value of a block in the sparse matrix.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] **val_block - Block to set to A(i, j).
               */
            void SetBlock(unsigned long block_i, unsigned long block_j, double **val_block);

            /*!
               * \brief Set the value of a block in the sparse matrix.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] **val_block - Block to set to A(i, j).
               */
            void SetBlock(unsigned long block_i, unsigned long block_j, double *val_block);

            /*!
             * \brief Adds the specified block to the sparse matrix.
             * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
             * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
             * \param[in] **val_block - Block to add to A(i, j).
             */
            void AddBlock(unsigned long block_i, unsigned long block_j, double **val_block);

            /*!
             * \brief Subtracts the specified block to the sparse matrix.
             * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
             * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
             * \param[in] **val_block - Block to subtract to A(i, j).
             */
            void SubtractBlock(unsigned long block_i, unsigned long block_j, double **val_block);

            /*!
               * \brief Copies the block (i, j) of the matrix-by-blocks structure in the internal variable *block.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               */
            double *GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j);

            /*!
               * \brief Set the value of a block in the sparse matrix.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] **val_block - Block to set to A(i, j).
               */
            void SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, double *val_block);

            /*!
             * \brief Subtracts the specified block to the sparse matrix.
             * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
             * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
             * \param[in] **val_block - Block to subtract to A(i, j).
             */
            void SubtractBlock_ILUMatrix(unsigned long block_i, unsigned long block_j, double *val_block);

            /*!
             * \brief Adds the specified value to the diagonal of the (i, i) subblock
             *        of the matrix-by-blocks structure.
             * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
             * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
             */
            void AddVal2Diag(unsigned long block_i, double val_matrix);

            /*!
             * \brief Sets the specified value to the diagonal of the (i, i) subblock
             *        of the matrix-by-blocks structure.
             * \param[in] block_i - Index of the block in the matrix-by-blocks structure.
             * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
             */
            void SetVal2Diag(unsigned long block_i, double val_matrix);

            /*!
               * \brief Calculates the matrix-vector product
               * \param[in] matrix
               * \param[in] vector
               * \param[out] product
               */
            void MatrixVectorProduct(double *matrix, double *vector, double *product);

            /*!
             * \brief Calculates the matrix-matrix product
             * \param[in] matrix_a
             * \param[in] matrix_b
             * \param[out] product
             */
            void MatrixMatrixProduct(double *matrix_a, double *matrix_b, double *product);

            /*!
             * \brief Deletes the values of the row i of the sparse matrix.
             * \param[in] i - Index of the row.
             */
            void DeleteValsRowi(unsigned long i);

            /*!
             * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
             * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
             * \param[in] rhs - Right-hand-side of the linear system.
             * \return Solution of the linear system (overwritten on rhs).
             */
            void Gauss_Elimination(unsigned long block_i, double* rhs);

            /*!
             * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
             * \param[in] Block - matrix-by-blocks structure.
             * \param[in] rhs - Right-hand-side of the linear system.
             * \return Solution of the linear system (overwritten on rhs).
             */
            void Gauss_Elimination(double* Block, double* rhs);

            /*!
               * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i, i) subblock and rhs.
               * \param[in] block_i - Index of the (i, i) subblock in the matrix-by-blocks structure.
               * \param[in] rhs - Right-hand-side of the linear system.
               * \return Solution of the linear system (overwritten on rhs).
               */
            void Gauss_Elimination_ILUMatrix(unsigned long block_i, double* rhs);

            /*!
               * \fn void MATH_Matrix::ProdBlockVector(unsigned long block_i, unsigned long block_j, double* vec);
               * \brief Performs the product of the block (i, j) by vector vec.
               * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] block_j - Indexes of the block in the matrix-by-blocks structure.
               * \param[in] vec - Vector to be multiplied by the block (i, j) of the sparse matrix A.
               * \return Product of A(i, j) by vector *vec (stored at *prod_block_vector).
               */
            void ProdBlockVector(unsigned long block_i, unsigned long block_j, const MATH_Vector & vec);

            /*!
               * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
               * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
               * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
               * \return prod Result of the product U(A)*vec (stored at *prod_row_vector).
               */
            void UpperProduct(MATH_Vector & vec, unsigned long row_i);

            /*!
               * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
               * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
               * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
               * \return prod Result of the product L(A)*vec (stored at *prod_row_vector).
               */
            void LowerProduct(MATH_Vector & vec, unsigned long row_i);

            /*!
               * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
               * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
               * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
               * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
               */
            void DiagonalProduct(MATH_Vector & vec, unsigned long row_i);

            /*!
               * \brief Send receive the solution using MPI.
               * \param[in] x - Solution..
               * \param[in] geometry - Geometrical definition of the problem.
               * \param[in] config - Definition of the particular problem.
               */
            void SendReceive_Solution(MATH_Vector & x, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
               * \brief Performs the product of i-th row of a sparse matrix by a vector.
               * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
               * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
               * \return Result of the product (stored at *prod_row_vector).
               */
            void RowProduct(const MATH_Vector & vec, unsigned long row_i);

            /*!
               * \brief Performs the product of a sparse matrix by a vector.
               * \param[in] vec - Vector to be multiplied by the sparse matrix A.
               * \param[out] prod - Result of the product.
               * \return Result of the product A*vec.
               */
            void MatrixVectorProduct(const MATH_Vector & vec, MATH_Vector & prod);

            /*!
             * \brief Performs the product of a sparse matrix by a MATH_Vector.
             * \param[in] vec - MATH_Vector to be multiplied by the sparse matrix A.
             * \param[out] prod - Result of the product.
             */
            void MatrixVectorProduct(const MATH_Vector & vec, MATH_Vector & prod, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \brief Performs the product of two block matrices.
             */
            void GetMultBlockBlock(double *c, double *a, double *b);

            /*!
             * \brief Performs the product of a block matrices by a vector.
             */
            void GetMultBlockVector(double *c, double *a, double *b);

            /*!
             * \brief Performs the subtraction of two matrices.
             */
            void GetSubsBlock(double *c, double *a, double *b);

            /*!
             * \brief Performs the subtraction of two vectors.
             */
            void GetSubsVector(double *c, double *a, double *b);

            /*!
             * \brief Inverse diagonal block.
             * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
             * \param[out] invBlock - Inverse block.
             */
            void InverseDiagonalBlock(unsigned long block_i, double *invBlock);

            /*!
             * \brief Inverse diagonal block.
             * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
             * \param[out] invBlock - Inverse block.
             */
            void InverseDiagonalBlock_ILUMatrix(unsigned long block_i, double *invBlock);

            /*!
             * \brief Inverse a block.
             * \param[in] Block - block matrix.
             * \param[out] invBlock - Inverse block.
             */
            void InverseBlock(double *Block, double *invBlock);

            /*!
             * \brief Build the Jacobi preconditioner.
             */
            void BuildJacobiPreconditioner(void);

            /*!
             * \brief Build the Jacobi preconditioner.
             */
            void BuildILUPreconditioner(void);

            /*!
             * \brief Build the Linelet preconditioner.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             */
            unsigned short BuildLineletPreconditioner(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \brief Multiply MATH_Vector by the preconditioner
             * \param[in] vec - MATH_Vector to be multiplied by the preconditioner.
             * \param[out] prod - Result of the product A*vec.
             */
            void ComputeJacobiPreconditioner(const MATH_Vector & vec, MATH_Vector & prod, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \brief Multiply MATH_Vector by the preconditioner
             * \param[in] vec - MATH_Vector to be multiplied by the preconditioner.
             * \param[out] prod - Result of the product A*vec.
             */
            void ComputeILUPreconditioner(const MATH_Vector & vec, MATH_Vector & prod, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
               * \brief Multiply MATH_Vector by the preconditioner
               * \param[in] vec - MATH_Vector to be multiplied by the preconditioner.
               * \param[out] prod - Result of the product A*vec.
               */
            void ComputeLU_SGSPreconditioner(const MATH_Vector & vec, MATH_Vector & prod, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \brief Multiply MATH_Vector by the preconditioner
             * \param[in] vec - MATH_Vector to be multiplied by the preconditioner.
             * \param[out] prod - Result of the product A*vec.
             */
            void ComputeLineletPreconditioner(const MATH_Vector & vec, MATH_Vector & prod, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
               * \brief Compute the residual Ax-b
               * \param[in] sol - MATH_Vector to be multiplied by the preconditioner.
               * \param[in] f - Result of the product A*vec.
               * \param[out] res - Result of the product A*vec.
               */
            void ComputeResidual(const MATH_Vector & sol, const MATH_Vector & f, MATH_Vector & res);

        private:
            unsigned long nPoint,                           /*!< \brief Number of points in the grid. */
                nPointDomain,                               /*!< \brief Number of points in the grid. */
                nVar,                                       /*!< \brief Number of variables. */
                nEqn;                                       /*!< \brief Number of equations. */
            double *matrix;                                 /*!< \brief Entries of the sparse matrix. */
            double *ILU_matrix;                             /*!< \brief Entries of the ILU sparse matrix. */
            unsigned long *row_ptr;                         /*!< \brief Pointers to the first element in each row. */
            unsigned long *col_ind;                         /*!< \brief Column index for each of the elements in val(). */
            unsigned long nnz;                              /*!< \brief Number of possible nonzero entries in the matrix. */
            double *block;                                  /*!< \brief Internal array to store a subblock of the matrix. */
            double *block_inverse;                          /*!< \brief Internal array to store a subblock of the matrix. */
            double *block_weight;                           /*!< \brief Internal array to store a subblock of the matrix. */
            double *prod_block_vector;                      /*!< \brief Internal array to store the product of a subblock with a vector. */
            double *prod_row_vector;                        /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
            double *aux_vector;                             /*!< \brief Auxilar array to store intermediate results. */
            double *sum_vector;                             /*!< \brief Auxilar array to store intermediate results. */
            double *invM;                                   /*!< \brief Inverse of (Jacobi) preconditioner. */
            bool *LineletBool;                              /*!< \brief Identify if a point belong to a linelet. */
            std::vector<unsigned long> *LineletPoint;       /*!< \brief Linelet structure. */
            unsigned long nLinelet;                         /*!< \brief Number of Linelets in the system. */
            double **UBlock, **invUBlock, **LBlock,
                **yVector, **zVector, **rVector, *LFBlock,
                *LyVector, *FzVector, *AuxVector;           /*!< \brief Arrays of the Linelet preconditioner methodology. */
            unsigned long max_nElem;
        };


        /*!
        * \class COMM_MatrixVectorProduct
        * \brief specialization of matrix-vector product that uses COMM_Matrix class
        */
        class MATH_Matrix_MatrixVectorProduct : public MATH_MatrixVectorProduct
        {
        private:
            MATH_Matrix* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
            GEOM::GEOM_Geometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
            TBOX::TBOX_Config* config; /*!< \brief pointer to matrix that defines the config. */

        public:

            /*!
            * \brief constructor of the class
            * \param[in] matrix_ref - matrix reference that will be used to define the products
            */
            MATH_Matrix_MatrixVectorProduct(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref);

            /*!
            * \brief destructor of the class
            */
            ~MATH_Matrix_MatrixVectorProduct() {}

            /*!
            * \brief operator that defines the MATH_Matrix-MATH_Vector product
            * \param[in] u - MATH_Vector that is being multiplied by the sparse matrix
            * \param[out] v - MATH_Vector that is the result of the product
            */
            void operator()(const MATH_Vector & u, MATH_Vector & v) const;
        };

        /*!
        * \class CILUPreconditioner
        * \brief specialization of preconditioner that uses MATH_Matrix class
        */
        class MATH_ILUPreconditioner : public MATH_Preconditioner 
        {
        private:
            MATH_Matrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
            GEOM::GEOM_Geometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
            TBOX::TBOX_Config* config; /*!< \brief pointer to matrix that defines the config. */

        public:

            /*!
            * \brief constructor of the class
            * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
            */
            MATH_ILUPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref);

            /*!
            * \brief destructor of the class
            */
            ~MATH_ILUPreconditioner() {}

            /*!
            * \brief operator that defines the preconditioner operation
            * \param[in] u - MATH_Vector that is being preconditioned
            * \param[out] v - MATH_Vector that is the result of the preconditioning
            */
            void operator()(const MATH_Vector & u, MATH_Vector & v) const;
        };

        /*!
        * \class CJacobiPreconditioner
        * \brief specialization of preconditioner that uses MATH_Matrix class
        */
        class MATH_JacobiPreconditioner : public MATH_Preconditioner
        {
        private:
            MATH_Matrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
            GEOM::GEOM_Geometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
            TBOX::TBOX_Config* config; /*!< \brief pointer to matrix that defines the config. */

        public:

            /*!
            * \brief constructor of the class
            * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
            */
            MATH_JacobiPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref);

            /*!
            * \brief destructor of the class
            */
            ~MATH_JacobiPreconditioner() {}

            /*!
            * \brief operator that defines the preconditioner operation
            * \param[in] u - MATH_Vector that is being preconditioned
            * \param[out] v - MATH_Vector that is the result of the preconditioning
            */
            void operator()(const MATH_Vector & u, MATH_Vector & v) const;
        };


        /*!
        * \class CLineletPreconditioner
        * \brief specialization of preconditioner that uses MATH_Matrix class
        */
        class MATH_LineletPreconditioner : public MATH_Preconditioner
        {
        private:
            MATH_Matrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
            GEOM::GEOM_Geometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
            TBOX::TBOX_Config* config; /*!< \brief pointer to matrix that defines the config. */

        public:

            /*!
            * \brief constructor of the class
            * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
            */
            MATH_LineletPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref);

            /*!
            * \brief destructor of the class
            */
            ~MATH_LineletPreconditioner() {}

            /*!
            * \brief operator that defines the preconditioner operation
            * \param[in] u - MATH_Vector that is being preconditioned
            * \param[out] v - MATH_Vector that is the result of the preconditioning
            */
            void operator()(const MATH_Vector & u, MATH_Vector & v) const;

        private:
            unsigned long nPoint,   /*!< \brief Number of points in the grid. */
                nPointDomain,           /*!< \brief Number of points in the grid. */
                nVar,                   /*!< \brief Number of variables. */
                nEqn;                   /*!< \brief Number of equations. */
            double *matrix;            /*!< \brief Entries of the sparse matrix. */
            double *ILU_matrix;         /*!< \brief Entries of the ILU sparse matrix. */
            unsigned long *row_ptr;    /*!< \brief Pointers to the first element in each row. */
            unsigned long *col_ind;    /*!< \brief Column index for each of the elements in val(). */
            unsigned long nnz;         /*!< \brief Number of possible nonzero entries in the matrix. */
            double *block;             /*!< \brief Internal array to store a subblock of the matrix. */
            double *block_inverse;             /*!< \brief Internal array to store a subblock of the matrix. */
            double *block_weight;             /*!< \brief Internal array to store a subblock of the matrix. */
            double *prod_block_vector; /*!< \brief Internal array to store the product of a subblock with a vector. */
            double *prod_row_vector;   /*!< \brief Internal array to store the product of a matrix-by-blocks "row" with a vector. */
            double *aux_vector;         /*!< \brief Auxilar array to store intermediate results. */
            double *sum_vector;         /*!< \brief Auxilar array to store intermediate results. */
            double *invM;              /*!< \brief Inverse of (Jacobi) preconditioner. */
            bool *LineletBool;                          /*!< \brief Identify if a point belong to a linelet. */
            std::vector<unsigned long> *LineletPoint;        /*!< \brief Linelet structure. */
            unsigned long nLinelet;                     /*!< \brief Number of Linelets in the system. */
            double **UBlock, **invUBlock, **LBlock,
                **yVector, **zVector, **rVector, *LFBlock,
                *LyVector, *FzVector, *AuxVector;           /*!< \brief Arrays of the Linelet preconditioner methodology. */
            unsigned long max_nElem;

        };

        /*!
        * \class CLU_SGSPreconditioner
        * \brief specialization of preconditioner that uses MATH_Matrix class
        */
        class MATH_LUSGSPreconditioner : public MATH_Preconditioner 
        {
        private:
            MATH_Matrix* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
            GEOM::GEOM_Geometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
            TBOX::TBOX_Config* config; /*!< \brief pointer to matrix that defines the config. */

        public:

            /*!
            * \brief constructor of the class
            * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
            */
            MATH_LUSGSPreconditioner(MATH_Matrix & matrix_ref, GEOM::GEOM_Geometry *geometry_ref, TBOX::TBOX_Config *config_ref);

            /*!
            * \brief destructor of the class
            */
            ~MATH_LUSGSPreconditioner() {}

            /*!
            * \brief operator that defines the preconditioner operation
            * \param[in] u - MATH_Vector that is being preconditioned
            * \param[out] v - MATH_Vector that is the result of the preconditioning
            */
            void operator()(const MATH_Vector & u, MATH_Vector & v) const;
        };
    }
}

#include "MATH_Matrix.inl"

#endif