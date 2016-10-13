



#ifndef ARIES_GEOM_GEOMETRYPERIODIC_HPP
#define ARIES_GEOM_GEOMETRYPERIODIC_HPP

#include "GEOM_Geometry.hpp"

#include "../Common/TBOX_Config.hpp"

#include "../Grid/GRID_VertexMPI.hpp"
#include "../Grid/GRID_Line.hpp"
#include "../Grid/GRID_Triangle.hpp"
#include "../Grid/GRID_Rectangle.hpp"
#include "../Grid/GRID_Tetrahedron.hpp"
#include "../Grid/GRID_Hexahedron.hpp"
#include "../Grid/GRID_Prism.hpp"
#include "../Grid/GRID_Pyramid.hpp"

namespace ARIES
{
    namespace GEOM
    {
        /*!
        * \class GEOM_GeometryPeriodic
        * \brief Class for defining a periodic boundary condition.
        * \author T. Economon, F. Palacios
        * \version 3.2.9 "eagle"
        */
        class GEOM_GeometryPeriodic : public GEOM_Geometry
        {
            GRID::GRID_Primal*** newBoundPer;            /*!< \brief Boundary std::vector for new periodic elements (primal grid information). */
            unsigned long *nNewElem_BoundPer;			/*!< \brief Number of new periodic elements of the boundary. */

        public:
            /*!
            * \brief Constructor of the class.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            GEOM_GeometryPeriodic(GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Destructor of the class.
            */
            ~GEOM_GeometryPeriodic(void);

            /*!
            * \brief Set the periodic boundaries of the grid.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void SetPeriodicBoundary(GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Set the Tecplot file.
            * \param[in] config_filename - Name of the file where the Tecplot
            *            information is going to be stored.
            */
            void SetTecPlot(char config_filename[TBOX::MAX_STRING_SIZE], bool new_file);

            /*!
            * \brief Write the .su2 file.
            * \param[in] config - Definition of the particular problem.
            * \param[in] val_mesh_out_filename - Name of the output file.
            */
            void SetMeshFile(GEOM_Geometry *geometry, TBOX::TBOX_Config *config, std::string val_mesh_out_filename);
        };
    }
}

#include "GEOM_GeometryPeriodic.inl"

#endif