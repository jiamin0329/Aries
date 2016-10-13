#ifndef ARIES_GRID_SURFACEMOVEMENT_INLINE
#define ARIES_GRID_SURFACEMOVEMENT_INLINE

namespace ARIES
{
    namespace GRID
    {
        inline unsigned short GRID_SurfaceMovement::GetnLevel(void) { return nLevel; }

        inline unsigned short GRID_SurfaceMovement::GetnFFDBox(void) { return nFFDBox; }

        inline bool GRID_SurfaceMovement::GetFFDBoxDefinition(void) { return FFDBoxDefinition; }

    }
}

#endif