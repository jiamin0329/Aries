/*********************************************************************************
 *                         ARIES Copyright(C), 2015.
 *
 *  \file    MATH_Vector.hpp
 *  \brief   Class for holding and manipulating vectors needed by linear solvers
 *           We could use the STL vector as a base class here, but this gives us
 *           more flexibility with the underlying data (e.g. we may decide to
 *           use a block storage scheme rather than a continuous storage scheme).
 *********************************************************************************
 *      Date        Author        Version                   Reason
 *    6/11/2015    Jiamin XU        1.0                  Initial release
 *
 *
 */


#ifndef ARIES_MATH_VECTOR_HPP
#define ARIES_MATH_VECTOR_HPP

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include <climits>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

namespace ARIES
{
    namespace MATH
    {
        const double eps = std::numeric_limits<double>::epsilon(); /*!< \brief machine epsilon */

        class MATH_Vector
        {
        public:
            /*!
             * \brief default constructor of the class.
             */
            MATH_Vector(void);

            /*!
             * \brief constructor of the class.
             * \param[in] size - number of elements locally
             * \param[in] val - default value for elements
             */
            MATH_Vector(const unsigned long & size, const double & val = 0.0);

            /*!
             * \brief constructor of the class.
             * \param[in] numBlk - number of blocks locally
             * \param[in] numVar - number of variables in each block
             * \param[in] val - default value for elements
             */
            MATH_Vector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const double & val = 0.0);

            /*!
             * \brief copy constructor of the class.
             * \param[in] u -Common_Vector that is being copied
             */
            MATH_Vector(const MATH_Vector & u);

            /*!
             * \brief constructor from array
             * \param[in] size - number of elements locally
             * \param[in] u_array - vector stored as array being copied
             */
            explicit MATH_Vector(const unsigned long & size, const double* u_array);

            /*!
             * \brief constructor from array
             * \param[in] numBlk - number of blocks locally
             * \param[in] numBlkDomain - number of blocks locally (without g cells)
             * \param[in] numVar - number of variables in each block
             * \param[in] u_array - vector stored as array being copied
             */
            explicit MATH_Vector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const double* u_array);

            /*!
             * \brief class destructor
             */
            virtual ~MATH_Vector();

            /*!
             * \brief Sets to zero all the entries of the vector.
             */
            void SetValZero(void);

            /*!
             * \brief Initialize the class.
             * \param[in] numBlk - number of blocks locally
             * \param[in] numVar - number of variables in each block
             * \param[in] val - default value for elements
             */
            void Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const double & val = 0.0);

            /*!
             * \brief return the number of local elements in the Common_Vector
             */
            unsigned long GetLocSize() const;

            /*!
             * \brief return the size of the MATH_Vector (over all processors)
             */
            unsigned long GetSize() const;

            /*!
             * \brief return the number of variables at each block (typically number per node)
             */
            unsigned short GetNVar() const;

            /*!
             * \brief return the number of blocks (typically number of nodes locally)
             */
            unsigned long GetNBlk() const;

            /*!
             * \brief return the number of blocks (typically number of nodes locally)
             */
            unsigned long GetNBlkDomain() const;

            /*!
             * \brief set calling MATH_Vector to scaling of another MATH_Vector
             * \param[in] a - scalar factor for x
             * \param[in] x - MATH_Vector that is being scaled
             */
            void Equals_AX(const double & a, MATH_Vector & x);

            /*!
             * \brief adds a scaled MATH_Vector to calling MATH_Vector
             * \param[in] a - scalar factor for x
             * \param[in] x - MATH_Vector that is being scaled
             */
            void Plus_AX(const double & a, MATH_Vector & x);

            /*!
             * \brief general linear combination of two MATH_Vectors
             * \param[in] a - scalar factor for x
             * \param[in] x - first MATH_Vector in linear combination
             * \param[in] b - scalar factor for y
             * \param[in] y - second MATH_Vector in linear combination
             */
            void Equals_AX_Plus_BY(const double & a, MATH_Vector & x, const double & b, MATH_Vector & y);

            /*!
             * \brief assignment operator with deep copy
             * \param[in] u - MATH_Vector whose values are being assigned
             */
            MATH_Vector & operator=(const MATH_Vector & u);

            /*!
             * \brief MATH_Vector=double assignment operator
             * \param[in] val - value assigned to each element of MATH_Vector
             */
            MATH_Vector & operator=(const double & val);

            /*!
             * \brief addition operator
             * \param[in] u - MATH_Vector being added to *this
             */
            MATH_Vector operator+(const MATH_Vector & u) const;

            /*!
             * \brief compound addition-assignment operator
             * \param[in] u - MATH_Vector being added to calling object
             */
            MATH_Vector & operator+=(const MATH_Vector & u);

            /*!
             * \brief subtraction operator
             * \param[in] u - MATH_Vector being subtracted from *this
             */
            MATH_Vector operator-(const MATH_Vector & u) const;

            /*!
             * \brief compound subtraction-assignment operator
             * \param[in] u - MATH_Vector being subtracted from calling object
             */
            MATH_Vector & operator-=(const MATH_Vector & u);

            /*!
             * \brief vector * scalar multiplication operator
             * \param[in] val - value to multiply *this by
             */
            MATH_Vector operator*(const double & val) const;

            /*!
             * \brief scalar * vector multiplication operator
             * \param[in] val - scalar value to multiply by
             * \param[in] u - MATH_Vector having its elements scaled
             */
            friend MATH_Vector operator*(const double & val, const MATH_Vector & u);

            /*!
             * \brief compound scalar multiplication-assignment operator
             * \param[in] val - value to multiply calling object by
             */
            MATH_Vector & operator*=(const double & val);

            /*!
             * \brief vector-scalar division operator (no scalar/vector operator)
             * \param[in] val - value to divide elements of *this by
             */
            MATH_Vector operator/(const double & val) const;

            /*!
             * \brief compound scalar division-assignment operator
             * \param[in] val - value to divide elements of calling object by
             */
            MATH_Vector & operator/=(const double & val);

            /*!
             * \brief indexing operator with assignment permitted
             * \param[in] i = local index to access
             */
            double & operator[](const unsigned long & i);

            /*!
             * \brief indexing operator with assignment not permitted
             * \param[in] i = local index to access
             */
            const double & operator[](const unsigned long & i) const;

            /*!
             * \brief the L2 norm of the MATH_Vector
             * \result the L2 norm
             */
            double norm() const;

            /*!
             * \brief copies the contents of the calling MATH_Vector into an array
             * \param[out] u_array - array into which information is being copied
             * \pre u_array must be allocated and have the same size as MATH_Vector
             */
            void CopyToArray(double* u_array);

            /*!
             * \brief Subtract val_residual to the residual.
             * \param[in] val_ipoint - index of the point where subtract the residual.
             * \param[in] val_residual - Value to subtract to the residual.
             */
            void SubtractBlock(unsigned long val_ipoint, double *val_residual);

            /*!
             * \brief Add val_residual to the residual.
             * \param[in] val_ipoint - index of the point where add the residual.
             * \param[in] val_residual - Value to add to the residual.
             */
            void AddBlock(unsigned long val_ipoint, double *val_residual);

            /*!
             * \brief Set val_residual to the residual.
             * \param[in] val_ipoint - index of the point where set the residual.
             * \param[in] val_var - inde of the residual to be set.
             * \param[in] val_residual - Value to set to the residual.
             */
            void SetBlock(unsigned long val_ipoint, unsigned short val_var, double val_residual);

            /*!
             * \brief Set val_residual to the residual.
             * \param[in] val_ipoint - index of the point where set the residual.
             * \param[in] val_residual - Value to set to the residual.
             */
            void SetBlock(unsigned long val_ipoint, double *val_residual);

            /*!
             * \brief Set the residual to zero.
             * \param[in] val_ipoint - index of the point where set the residual.
             */
            void SetBlock_Zero(unsigned long val_ipoint);

            /*!
             * \brief Set the velocity residual to zero.
             * \param[in] val_ipoint - index of the point where set the residual.
             */
            void SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var);

            /*!
             * \brief Get the value of the residual.
             * \param[in] val_ipoint - index of the point where set the residual.
             * \return Pointer to the residual.
             */
            double *GetBlock(unsigned long val_ipoint);

            /*!
             * \brief Get the value of the residual.
             * \param[in] val_ipoint - index of the point where set the residual.
             * \param[in] val_var - inde of the residual to be set.
             * \return Value of the residual.
             */
            double GetBlock(unsigned long val_ipoint, unsigned short val_var);

            /*!
             * \brief dot-product between two MATH_Vectors
             * \param[in] u - first MATH_Vector in dot product
             * \param[in] v - second MATH_Vector in dot product
             */
            friend double dotProd(const MATH_Vector & u, const MATH_Vector & v);

        private:
            unsigned long d_nElm;                 /*!< \brief total number of elements (or number elements on this processor) */
            unsigned long d_nElmDomain;           /*!< \brief total number of elements (or number elements on this processor without Ghost cells) */
            unsigned long d_nElmGlobal;           /*!< \brief total number of elements over all processors */
            unsigned short d_nVar;                /*!< \brief number of elements in a block */
            unsigned long d_nBlk;                 /*!< \brief number of blocks (or number of blocks on this processor) */
            unsigned long d_nBlkDomain;           /*!< \brief number of blocks (or number of blocks on this processor without Ghost cells) */
            double* d_vec_val;                    /*!< \brief storage for the element values */
        };
    }
}


namespace ARIES
{
    namespace MATH
    {
        /*!
        * \class CMatrixVectorProduct
        * \brief abstract base class for defining matrix-vector products
        * \author J. Hicken.
        * \version 3.2.9 "eagle"
        *
        * The Krylov-subspace solvers require only matrix-vector products and
        * not the actual matrix/Jacobian.  We need some way to indicate which
        * function will perform the product.  However, sometimes the
        * functions that define the product will require different numbers
        * and types of inputs.  For example, the forward-difference
        * approximation to a Jacobian-vector product requires the vector that
        * defines the Jacobian and a perturbation parameter.  The
        * CMatrixVectorProduct class is used to derive child classes that can
        * handle the different types of matrix-vector products and still be
        * passed to a single implementation of the Krylov solvers.
        */
        class MATH_MatrixVectorProduct
        {
        public:
            virtual ~MATH_MatrixVectorProduct() = 0; ///< class destructor
            virtual void operator()(const MATH_Vector & u, MATH_Vector & v)
                const = 0; ///< matrix-vector product operation
        };
        inline MATH_MatrixVectorProduct::~MATH_MatrixVectorProduct() {}

        /*!
        * \class CPreconditioner
        * \brief abstract base class for defining preconditioning operation
        * \author J. Hicken.
        * \version 3.2.9 "eagle"
        *
        * See the remarks regarding the CMatrixVectorProduct class.  The same
        * idea applies here to the preconditioning operation.
        */
        class MATH_Preconditioner
        {
        public:
            virtual ~MATH_Preconditioner() = 0; ///< class destructor
            virtual void operator()(const MATH_Vector & u, MATH_Vector & v)
                const = 0; ///< preconditioning operation
        };
        inline MATH_Preconditioner::~MATH_Preconditioner() {}
    }//MATH
}//ARIES

#include "MATH_Vector.inl"

#endif