/*********************************************************************************
*                         ARIES Copyright(C), 2015.
*
*  \file    MATH_Vector.inl
*  \brief   Class for holding and manipulating vectors needed by linear solvers
*           We could use the STL vector as a base class here, but this gives us
*           more flexibility with the underlying data (e.g. we may decide to
*           use a block storage scheme rather than a continuous storage scheme).
*********************************************************************************
*      Date        Author        Version                   Reason
*    06/11/2015    Jiamin XU        1.0                  Initial release
*
*
*/

#ifndef ARIES_MATH_VECTOR_INLINE
#define ARIES_MATH_VECTOR_INLINE

namespace ARIES
{
    namespace MATH
    {
        inline void MATH_Vector::SetValZero(void) 
        {
            for (unsigned long i = 0; i <= d_nElm-1; i++)
                d_vec_val[i] = 0.0;
        }

        inline unsigned long MATH_Vector::GetLocSize() const 
        { 
            return d_nElm; 
        }

        inline unsigned long MATH_Vector::GetSize() const 
        {
#ifdef HAVE_MPI
            return d_nElmGlobal;
#else
            return (unsigned long)d_nElm;
#endif
        }

        inline unsigned short MATH_Vector::GetNVar() const 
        { 
            return d_nVar; 
        }

        inline unsigned long MATH_Vector::GetNBlk() const 
        { 
            return d_nBlk; 
        }

        inline unsigned long MATH_Vector::GetNBlkDomain() const 
        { 
            return d_nBlkDomain; 
        }

        inline double & MATH_Vector::operator[](const unsigned long & i) 
        { 
            return d_vec_val[i]; 
        }

        inline const double & MATH_Vector::operator[](const unsigned long & i) const 
        { 
            return d_vec_val[i]; 
        }
    }
}

#endif