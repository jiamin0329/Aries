/*********************************************************************************
 *                         ARIES Copyright(C), 2015.
 *
 *  \file    MATH_Vector.cpp
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

#include "MATH_Vector.hpp"

namespace ARIES
{
    namespace MATH
    {
        MATH_Vector::MATH_Vector(void)
        {
            d_vec_val = NULL;
        }

        MATH_Vector::MATH_Vector(const unsigned long & size, const double & val)
        {
            d_nElm = size;
            d_nElmDomain = size;
            d_nBlk = d_nElm;
            d_nBlkDomain = d_nElmDomain;
            d_nVar = 1;

            /*--- Check for invalid size, then allocate memory and initialize values ---*/
            if ((d_nElm <= 0) || (d_nElm >= UINT_MAX))
            {
                std::cerr << "Common_Vector::Common_Vector(unsigned int, double): "
                    << "invalid input: size = " << size << std::endl;
                throw(-1);
            }

            d_vec_val = new double[d_nElm];
            for (unsigned int i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = val;

#ifdef HAVE_MPI
            unsigned long nElmLocal = (unsigned long)nElm;
            MPI_Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
        }

        MATH_Vector::MATH_Vector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const double & val)
        {
            d_nElm = numBlk*numVar;
            d_nElmDomain = numBlkDomain*numVar;
            d_nBlk = numBlk;
            d_nBlkDomain = numBlkDomain;
            d_nVar = numVar;

            /*--- Check for invalid input, then allocate memory and initialize values ---*/
            if ((d_nElm <= 0) || (d_nElm >= ULONG_MAX))
            {
                std::cerr << "CSysVector::CSysVector(unsigned int, unsigned int, double): "
                    << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << std::endl;
                throw(-1);
            }

            d_vec_val = new double[d_nElm];
            for (unsigned int i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = val;

#ifdef HAVE_MPI
            int myrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            unsigned long nElmLocal = (unsigned long)nElm;
            MPI_Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
        }

        MATH_Vector::MATH_Vector(const MATH_Vector & u)
        {
            /*--- Copy size information, allocate memory, and initialize values ---*/
            d_nElm = u.d_nElm;
            d_nElmDomain = u.d_nElmDomain;
            d_nBlk = u.d_nBlk;
            d_nBlkDomain = u.d_nBlkDomain;
            d_nVar = u.d_nVar;

            d_vec_val = new double[d_nElm];
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = u.d_vec_val[i];

#ifdef HAVE_MPI
            nElmGlobal = u.nElmGlobal;
#endif
        }

        MATH_Vector::MATH_Vector(const unsigned long & size, const double* u_array)
        {
            d_nElm = size;
            d_nElmDomain = size;
            d_nBlk = d_nElm;
            d_nBlkDomain = d_nElmDomain;
            d_nVar = 1;

            /*--- Check for invalid size, then allocate memory and initialize values ---*/
            if ((d_nElm <= 0) || (d_nElm >= ULONG_MAX))
            {
                std::cerr << "CSysVector::CSysVector(unsigned int, double*): "
                    << "invalid input: size = " << size << std::endl;
                throw(-1);
            }

            d_vec_val = new double[d_nElm];
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = u_array[i];

#ifdef HAVE_MPI
            int myrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            unsigned long nElmLocal = (unsigned long)nElm;
            MPI_Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
        }

        MATH_Vector::MATH_Vector(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const double* u_array)
        {
            d_nElm = numBlk*numVar;
            d_nElmDomain = numBlkDomain*numVar;
            d_nBlk = numBlk;
            d_nBlkDomain = numBlkDomain;
            d_nVar = numVar;

            /*--- check for invalid input, then allocate memory and initialize values ---*/
            if ((d_nElm <= 0) || (d_nElm >= ULONG_MAX))
            {
                std::cerr << "CSysVector::CSysVector(unsigned int, unsigned int, double*): "
                    << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << std::endl;
                throw(-1);
            }

            d_vec_val = new double[d_nElm];
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = u_array[i];

#ifdef HAVE_MPI
            unsigned long nElmLocal = (unsigned long)nElm;
            MPI_Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
        }

        MATH_Vector::~MATH_Vector()
        {
            delete[] d_vec_val;

            d_nElm = 0;
            d_nElmDomain = 0;
            d_nBlk = 0;
            d_nBlkDomain = 0;
            d_nVar = 0;
        }

        void MATH_Vector::Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const double & val)
        {
            d_nElm = numBlk*numVar;
            d_nElmDomain = numBlkDomain*numVar;
            d_nBlk = numBlk;
            d_nBlkDomain = numBlkDomain;
            d_nVar = numVar;

            /*--- Check for invalid input, then allocate memory and initialize values ---*/
            if ((d_nElm <= 0) || (d_nElm >= ULONG_MAX))
            {
                std::cerr << "CSysVector::CSysVector(unsigned int, unsigned int, double): "
                    << "invalid input: numBlk, numVar = " << numBlk << "," << numVar << std::endl;
                throw(-1);
            }

            d_vec_val = new double[d_nElm];
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = val;

#ifdef HAVE_MPI
            unsigned long nElmLocal = (unsigned long)nElm;
            MPI_Allreduce(&nElmLocal, &nElmGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

        }

        void MATH_Vector::Equals_AX(const double & a, MATH_Vector & x)
        {
            /*--- check that *this and x are compatible ---*/
            if (d_nElm != x.d_nElm)
            {
                std::cerr << "MATH_Vector::Equals_AX(): " << "sizes do not match";
                throw(-1);
            }
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = a * x.d_vec_val[i];
        }

        void MATH_Vector::Plus_AX(const double & a, MATH_Vector & x)
        {
            /*--- check that *this and x are compatible ---*/
            if (d_nElm != x.d_nElm)
            {
                std::cerr << "CSysVector::Plus_AX(): " << "sizes do not match";
                throw(-1);
            }
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] += a * x.d_vec_val[i];
        }

        void MATH_Vector::Equals_AX_Plus_BY(const double & a, MATH_Vector & x, const double & b, MATH_Vector & y)
        {
            /*--- check that *this, x and y are compatible ---*/
            if ((d_nElm != x.d_nElm) || (d_nElm != y.d_nElm))
            {
                std::cerr << "CSysVector::Equals_AX_Plus_BY(): " << "sizes do not match";
                throw(-1);
            }
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = a * x.d_vec_val[i] + b * y.d_vec_val[i];
        }

        MATH_Vector & MATH_Vector::operator=(const MATH_Vector & u)
        {
            /*--- check if self-assignment, otherwise perform deep copy ---*/
            if (this == &u) return *this;
            delete[] d_vec_val; // in case the size is different
            d_nElm = u.d_nElm;
            d_nElmDomain = u.d_nElmDomain;

            d_nBlk = u.d_nBlk;
            d_nBlkDomain = u.d_nBlkDomain;

            d_nVar = u.d_nVar;
            d_vec_val = new double[d_nElm];
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = u.d_vec_val[i];

#ifdef HAVE_MPI
            nElmGlobal = u.nElmGlobal;
#endif
            return *this;
        }

        MATH_Vector & MATH_Vector::operator=(const double & val)
        {
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] = val;
            return *this;
        }

        MATH_Vector MATH_Vector::operator+(const MATH_Vector & u) const
        {
            /*--- Use copy constructor and compound addition-assignment ---*/
            MATH_Vector sum(*this);
            sum += u;
            return sum;
        }

        MATH_Vector & MATH_Vector::operator+=(const MATH_Vector & u)
        {
            /*--- Check for consistent sizes, then add elements ---*/
            if (d_nElm != u.d_nElm)
            {
                std::cerr << "MATH_Vector::operator+=(MATH_Vector): " << "sizes do not match";
                throw(-1);
            }
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] += u.d_vec_val[i];
            return *this;
        }

        MATH_Vector MATH_Vector::operator-(const MATH_Vector & u) const
        {
            /*--- Use copy constructor and compound subtraction-assignment ---*/
            MATH_Vector diff(*this);
            diff -= u;
            return diff;
        }

        MATH_Vector & MATH_Vector::operator-=(const MATH_Vector & u)
        {
            /*--- Check for consistent sizes, then subtract elements ---*/
            if (d_nElm != u.d_nElm)
            {
                std::cerr << "MATH_Vector::operator-=(MATH_Vector): " << "sizes do not match";
                throw(-1);
            }
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] -= u.d_vec_val[i];
            return *this;
        }

        MATH_Vector MATH_Vector::operator*(const double & val) const
        {
            /*--- use copy constructor and compound scalar multiplication-assignment ---*/
            MATH_Vector prod(*this);
            prod *= val;
            return prod;
        }

        MATH_Vector operator*(const double & val, const MATH_Vector & u)
        {
            /*--- use copy constructor and compound scalar multiplication-assignment ---*/
            MATH_Vector prod(u);
            prod *= val;
            return prod;
        }

        MATH_Vector & MATH_Vector::operator*=(const double & val)
        {
            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] *= val;
            return *this;
        }

        MATH_Vector MATH_Vector::operator/(const double & val) const
        {
            /*--- use copy constructor and compound scalar division-assignment ---*/
            MATH_Vector quotient(*this);
            quotient /= val;
            return quotient;
        }

        MATH_Vector & MATH_Vector::operator/=(const double & val)
        {

            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                d_vec_val[i] /= val;
            return *this;
        }

        double MATH_Vector::norm() const
        {
            /*--- just call dotProd on this*, then sqrt ---*/
            double val = dotProd(*this, *this);
            if (val < 0.0)
            {
                std::cerr << "MATH_Vector::norm(): " << "inner product of MATH_Vector is negative";
                throw(-1);
            }
            return sqrt(val);
        }

        void MATH_Vector::CopyToArray(double* u_array)
        {

            for (unsigned long i = 0; i <= d_nElm - 1; i++)
                u_array[i] = d_vec_val[i];
        }

        void MATH_Vector::AddBlock(unsigned long val_ipoint, double *val_residual)
        {
            unsigned short iVar;

            for (iVar = 0; iVar <= d_nVar - 1; iVar++)
                d_vec_val[val_ipoint*d_nVar + iVar] += val_residual[iVar];
        }

        void MATH_Vector::SubtractBlock(unsigned long val_ipoint, double *val_residual)
        {
            unsigned short iVar;

            for (iVar = 0; iVar <= d_nVar - 1; iVar++)
                d_vec_val[val_ipoint*d_nVar + iVar] -= val_residual[iVar];
        }

        void MATH_Vector::SetBlock(unsigned long val_ipoint, double *val_residual)
        {
            unsigned short iVar;

            for (iVar = 0; iVar <= d_nVar - 1; iVar++)
                d_vec_val[val_ipoint*d_nVar + iVar] = val_residual[iVar];
        }

        void MATH_Vector::SetBlock(unsigned long val_ipoint, unsigned short val_var, double val_residual)
        {
            d_vec_val[val_ipoint*d_nVar + val_var] = val_residual;
        }

        void MATH_Vector::SetBlock_Zero(unsigned long val_ipoint)
        {
            unsigned short iVar;

            for (iVar = 0; iVar <= d_nVar - 1; iVar++)
                d_vec_val[val_ipoint*d_nVar + iVar] = 0.0;
        }

        void MATH_Vector::SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var)
        {
            d_vec_val[val_ipoint*d_nVar + val_var] = 0.0;
        }

        double MATH_Vector::GetBlock(unsigned long val_ipoint, unsigned short val_var)
        {
            return d_vec_val[val_ipoint*d_nVar + val_var];
        }

        double *MATH_Vector::GetBlock(unsigned long val_ipoint)
        {
            return &d_vec_val[val_ipoint*d_nVar];
        }

        double dotProd(const MATH_Vector & u, const MATH_Vector & v)
        {
            /*--- check for consistent sizes ---*/
            if (u.d_nElm != v.d_nElm)
            {
                std::cerr << "MATH_Vector friend dotProd(MATH_Vector, MATH_Vector): "
                    << "MATH_Vector sizes do not match";
                throw(-1);
            }

            /*--- find local inner product and, if a parallel run, sum over all processors (we use nElemDomain instead of nElem) ---*/
            double loc_prod = 0.0;
            for (unsigned long i = 0; i <= u.d_nElmDomain - 1; i++)
                loc_prod += u.d_vec_val[i] * v.d_vec_val[i];
            double prod = 0.0;

#ifdef HAVE_MPI
            MPI_Allreduce(&loc_prod, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
            prod = loc_prod;
#endif
            return prod;
        }

    }
}