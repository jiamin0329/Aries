




/*!
* \class GEOM_GeometryMultigrid
* \brief Class for defining the multigrid geometry, the main delicated part is the
*        agglomeration stage, which is done in the declaration.
* \author F. Palacios
* \version 3.2.9 "eagle"
*/

#ifndef ARIES_GEOM_GEOMETRYMULTIGRID_HPP
#define ARIES_GEOM_GEOMETRYMULTIGRID_HPP

#include "../Common/TBOX_Config.hpp"
#include "GEOM_Geometry.hpp"

namespace ARIES
{
    namespace GEOM
    {
        class GEOM_GeometryMultigrid : public GEOM_Geometry
        {

        public:

            /*!
            * \brief Constructor of the class.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iMesh - Level of the multigrid.
            * \param[in] iZone - Current zone in the mesh.
            */
            GEOM_GeometryMultigrid(GEOM_Geometry ***geometry, TBOX::TBOX_Config **config_container, unsigned short iMesh, unsigned short iZone);

            /*!
            * \brief Destructor of the class.
            */
            ~GEOM_GeometryMultigrid(void);

            /*!
            * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
            * \param[in] CVPoint - Control volume to be agglomerated.
            * \param[in] marker_seed - Marker of the seed.
            * \param[in] fine_grid - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \return <code>TRUE</code> or <code>FALSE</code> depending if the control volume can be agglomerated.
            */
            bool SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, GEOM_Geometry *fine_grid, TBOX::TBOX_Config *config);

            /*!
            * \brief Determine if a can be agglomerated using geometrical criteria.
            * \param[in] iPoint - Seed point.
            * \param[in] fine_grid - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            bool GeometricalCheck(unsigned long iPoint, GEOM_Geometry *fine_grid, TBOX::TBOX_Config *config);

            /*!
            * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
            * \param[in] Suitable_Indirect_Neighbors - List of Indirect Neighbours that can be agglomerated.
            * \param[in] iPoint - Seed point.
            * \param[in] Index_CoarseCV - Index of agglomerated point.
            * \param[in] fine_grid - Geometrical definition of the problem.
            */
            void SetSuitableNeighbors(std::vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint, unsigned long Index_CoarseCV, GEOM_Geometry *fine_grid);

            /*!
            * \brief Set boundary vertex.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void SetVertex(GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Set points which surround a point.
            * \param[in] geometry - Geometrical definition of the problem.
            */
            void SetPoint_Connectivity(GEOM_Geometry *geometry);

            /*!
            * \brief Function declaration to avoid partially overridden classes.
            */
            void SetPoint_Connectivity(void);

            /*!
            * \brief Set the edge structure of the agglomerated control volume.
            * \param[in] config - Definition of the particular problem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] action - Allocate or not the new elements.
            */
            void SetControlVolume(TBOX::TBOX_Config *config, GEOM_Geometry *geometry, unsigned short action);

            /*!
            * \brief Mach the near field boundary condition.
            * \param[in] config - Definition of the particular problem.
            */
            void MatchNearField(TBOX::TBOX_Config *config);

            /*!
            * \brief Mach the near field boundary condition.
            * \param[in] config - Definition of the particular problem.
            */
            void MatchActuator_Disk(TBOX::TBOX_Config *config);

            /*!
            * \brief Mach the interface boundary condition.
            * \param[in] config - Definition of the particular problem.
            */
            void MatchInterface(TBOX::TBOX_Config *config);

            /*!
            * \brief Set boundary vertex structure of the agglomerated control volume.
            * \param[in] config - Definition of the particular problem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] action - Allocate or not the new elements.
            */
            void SetBoundControlVolume(TBOX::TBOX_Config *config, GEOM_Geometry *geometry, unsigned short action);

            /*!
            * \brief Set a representative coordinates of the agglomerated control volume.
            * \param[in] geometry - Geometrical definition of the problem.
            */
            void SetCoord(GEOM_Geometry *geometry);

            /*!
            * \brief Set the rotational velocity at each grid point on a coarse mesh.
            * \param[in] config - Definition of the particular problem.
            */
            void SetRotationalVelocity(TBOX::TBOX_Config *config);

            /*!
            * \brief Set the grid velocity at each node in the coarse mesh level.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iter - Current physical time step.
            */
            void SetGridVelocity(TBOX::TBOX_Config *config, unsigned long iter);

            /*!
            * \brief Set the grid velocity at each node in the coarse mesh level based
            *        on a restriction from a finer mesh.
            * \param[in] fine_mesh - Geometry container for the finer mesh level.
            * \param[in] config - Definition of the particular problem.
            */
            void SetRestricted_GridVelocity(GEOM_Geometry *fine_mesh, TBOX::TBOX_Config *config);

            /*!
            * \brief Find and store the closest neighbor to a vertex.
            * \param[in] config - Definition of the particular problem.
            */
            void FindNormal_Neighbor(TBOX::TBOX_Config *config);

            /*!
            * \brief Indentify geometrical planes in the mesh
            */
            void SetGeometryPlanes(TBOX::TBOX_Config *config);

            /*!
            * \brief Get geometrical planes in the mesh
            */
            std::vector<double> GetGeometryPlanes();

            /*!
            * \brief Get x coords of geometrical planes in the mesh
            */
            std::vector<std::vector<double> > GetXCoord();

            /*!
            * \brief Get y coords of geometrical planes in the mesh
            */
            std::vector<std::vector<double> > GetYCoord();

            /*!
            * \brief Get z coords of geometrical planes in the mesh
            */
            std::vector<std::vector<double> > GetZCoord();

            /*!
            * \brief Get all points on a geometrical plane in the mesh
            */
            std::vector<std::vector<unsigned long> > GetPlanarPoints();
        };
    }
}

#include "GEOM_GeometryMultigrid.inl"

#endif