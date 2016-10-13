#ifndef ARIES_GRID_VOLUMETRICMOVEMENT_INLINE
#define ARIES_GRID_VOLUMETRICMOVEMENT_INLINE

namespace ARIES
{
    namespace GRID
    {
        inline double GRID_VolumetricMovement::Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22)
        {
            return A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20);
        }
    }
}

#endif