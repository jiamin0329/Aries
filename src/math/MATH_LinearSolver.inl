
#ifndef ARIES_MATH_LINEARSOLVER_INLINE
#define ARIES_MATH_LINEARSOLVER_INLINE

namespace ARIES
{
    namespace MATH
    {
        inline double MATH_LinearSolver::Sign(const double & x, const double & y) const
        {
            if (y == 0.0)
                return 0.0;
            else 
            {
                return (y < 0 ? -fabs(x) : fabs(x));
            }
        }
    }
}

#endif
