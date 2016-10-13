

#ifndef ARIES_GRID_GRIDMOVEMENT_HPP
#define ARIES_GRID_GRIDMOVEMENT_HPP

#include "../Geometry/GEOM_Geometry.hpp"
#include "../Common/TBOX_Config.hpp"

namespace ARIES
{
    namespace GRID
    {
        class GRID_Gridmovement 
        {
        public:
            /*!
            * \brief Constructor of the class.
            */
            GRID_Gridmovement(void);

            /*!
            * \brief Destructor of the class.
            */
            ~GRID_Gridmovement(void);


            /*!
            * \brief A pure virtual member.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            virtual void SetSurface_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

        };
    }
}

#include "GRID_Gridmovement.inl"

#endif