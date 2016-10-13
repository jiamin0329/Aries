



/*!
 * \class GEOM_GeometryPhysical
 * \brief Class for reading a defining the primal grid which is read from the
 *        grid file in .su2 format.
 * \author F. Palacios
 * \version 3.2.9 "eagle"
 */


#ifndef ARIES_GEOM_GEOMETRYPHYSICAL_HPP
#define ARIES_GEOM_GEOMETRYPHYSICAL_HPP

#include "../Common/TBOX_Config.hpp"
#include "GEOM_Geometry.hpp"

namespace ARIES
{
    namespace GEOM
    {
        class GEOM_GeometryPhysical : public GEOM_Geometry
        {

        public:

            /*!
             * \brief Constructor of the class.
             */
            GEOM_GeometryPhysical(void);

            /*!
             * \overload
             * \brief Reads the geometry of the grid and adjust the boundary
             *        conditions with the configuration file.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_filename - Name of the file with the grid information.
             * \param[in] val_format - Format of the file with the grid information.
             * \param[in] val_iZone - Domain to be read from the grid file.
             * \param[in] val_nZone - Total number of domains in the grid file.
             */
            GEOM_GeometryPhysical(TBOX::TBOX_Config *config, unsigned short val_iZone, unsigned short val_nZone);

            /*!
             * \overload
             * \brief Reads the geometry of the grid and adjust the boundary
             *        conditions with the configuration file.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_filename - Name of the file with the grid information.
             * \param[in] val_format - Format of the file with the grid information.
             * \param[in] val_iZone - Domain to be read from the grid file.
             * \param[in] val_nZone - Total number of domains in the grid file.
             */
            GEOM_GeometryPhysical(GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \overload
             * \brief Reads the geometry of the grid and adjust the boundary
             *        conditions with the configuration file for parmetis version.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_filename - Name of the file with the grid information.
             * \param[in] val_format - Format of the file with the grid information.
             * \param[in] val_iZone - Domain to be read from the grid file.
             * \param[in] val_nZone - Total number of domains in the grid file.
             */
            GEOM_GeometryPhysical(GEOM_Geometry *geometry, TBOX::TBOX_Config *config, int options);

            /*!
             * \brief Destructor of the class.
             */
            ~GEOM_GeometryPhysical(void);

            /*!
             * \brief Set the send receive boundaries of the grid.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_domain - Number of domains for parallelization purposes.
             */
            void SetSendReceive(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the send receive boundaries of the grid.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_domain - Number of domains for parallelization purposes.
             */
            void SetBoundaries(TBOX::TBOX_Config *config);

            /*!
             * \brief Get the local index that correspond with the global numbering index.
             * \param[in] val_ipoint - Global point.
             * \returns Local index that correspond with the global index.
             */
            long GetGlobal_to_Local_Point(long val_ipoint);

            /*!
             * \brief Get the local marker that correspond with the global marker.
             * \param[in] val_ipoint - Global marker.
             * \returns Local marker that correspond with the global index.
             */
            unsigned short GetGlobal_to_Local_Marker(unsigned short val_imarker);

            /*!
             * \brief Reads the geometry of the grid and adjust the boundary
             *        conditions with the configuration file in parallel (for parmetis).
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_filename - Name of the file with the grid information.
             * \param[in] val_format - Format of the file with the grid information.
             * \param[in] val_iZone - Domain to be read from the grid file.
             * \param[in] val_nZone - Total number of domains in the grid file.
             */
            void Read_SU2_Format_Parallel(TBOX::TBOX_Config *config, std::string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);


            /*!
             * \brief Reads the geometry of the grid and adjust the boundary
             *        conditions with the configuration file in parallel (for parmetis).
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_filename - Name of the file with the grid information.
             * \param[in] val_format - Format of the file with the grid information.
             * \param[in] val_iZone - Domain to be read from the grid file.
             * \param[in] val_nZone - Total number of domains in the grid file.
             */
            void Read_CGNS_Format_Parallel(TBOX::TBOX_Config *config, std::string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone);

            /*!
             * \brief Find repeated nodes between two elements to identify the common face.
             * \param[in] first_elem - Identification of the first element.
             * \param[in] second_elem - Identification of the second element.
             * \param[in] face_first_elem - Index of the common face for the first element.
             * \param[in] face_second_elem - Index of the common face for the second element.
             * \return It provides 0 or 1 depending if there is a common face or not.
             */
            bool FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                unsigned short &face_second_elem);

            /*!
             * \brief Computes the distance to the nearest no-slip wall for each grid node.
             * \param[in] config - Definition of the particular problem.
             */
            void ComputeWall_Distance(TBOX::TBOX_Config *config);

            /*!
             * \brief Compute surface area (positive z-direction) for force coefficient non-dimensionalization.
             * \param[in] config - Definition of the particular problem.
             */
            void SetPositive_ZArea(TBOX::TBOX_Config *config);

            /*!
             * \brief Set points which surround a point.
             */
            void SetPoint_Connectivity(void);

            /*!
             * \brief Set a renumbering using a Reverse Cuthill-McKee Algorithm
             * \param[in] config - Definition of the particular problem.
             */
            void SetRCM_Ordering(TBOX::TBOX_Config *config);

            /*!
             * \brief Function declaration to avoid partially overridden classes.
             * \param[in] geometry - Geometrical definition of the problem.
             */
            void SetPoint_Connectivity(GEOM_Geometry *geometry);

            /*!
             * \brief Set elements which surround an element.
             */
            void SetElement_Connectivity(void);

            /*!
             * \brief Set the volume element associated to each boundary element.
             */
            void SetBoundVolume(void);

            /*!
             * \brief Set boundary vertex.
             * \param[in] config - Definition of the particular problem.
             */
            void SetVertex(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the center of gravity of the face, elements and edges.
             */
            void SetCG(void);

            /*!
             * \brief Set the edge structure of the control volume.
             * \param[in] config - Definition of the particular problem.
             * \param[in] action - Allocate or not the new elements.
             */
            void SetControlVolume(TBOX::TBOX_Config *config, unsigned short action);

            /*!
             * \brief Visualize the structure of the control volume(s).
             * \param[in] config - Definition of the particular problem.
             * \param[in] action - Allocate or not the new elements.
             */
            void VisualizeControlVolume(TBOX::TBOX_Config *config, unsigned short action);

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
             * \brief Mach the interface boundary condition.
             * \param[in] config - Definition of the particular problem.
             * \param[in] geometry_donor - Geometry of the donor zone.
             * \param[in] config_donor - Definition of the donor problem.
             */
            void MatchZone(TBOX::TBOX_Config *config, GEOM_Geometry *geometry_donor, TBOX::TBOX_Config *config_donor,
                unsigned short val_iZone, unsigned short val_nZone);

            /*!
             * \brief Set boundary vertex structure of the control volume.
             * \param[in] config - Definition of the particular problem.
             * \param[in] action - Allocate or not the new elements.
             */
            void SetBoundControlVolume(TBOX::TBOX_Config *config, unsigned short action);

            /*!
             * \brief Set the Tecplot file.
             * \param[in] config_filename - Name of the file where the Tecplot
             *            information is going to be stored.
             * \param[in] new_file - Create a new file.
             */
            void SetTecPlot(char config_filename[TBOX::MAX_STRING_SIZE], bool new_file);

            /*!
             * \brief Set the output file for boundaries in Tecplot
             * \param[in] config - Definition of the particular problem.
             * \param[in] mesh_filename - Name of the file where the Tecplot
             *            information is going to be stored.
             * \param[in] new_file - Create a new file.
             */
            void SetBoundTecPlot(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file, TBOX::TBOX_Config *config);

            /*!
             * \brief Set the output file for boundaries in STL CAD format
             * \param[in] config - Definition of the particular problem.
             * \param[in] mesh_filename - Name of the file where the STL
             *            information is going to be stored.
             * \param[in] new_file - Create a new file.
             */
            void SetBoundSTL(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file, TBOX::TBOX_Config *config);

            /*!
             * \brief Check the volume element orientation.
             * \param[in] config - Definition of the particular problem.
             */
            void Check_IntElem_Orientation(TBOX::TBOX_Config *config);

            /*!
             * \brief Check the volume element orientation.
             * \param[in] config - Definition of the particular problem.
             */
            void Check_BoundElem_Orientation(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the domains for grid grid partitioning using METIS.
             * \param[in] config - Definition of the particular problem.
             */
            void SetColorGrid(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the domains for grid grid partitioning using ParMETIS.
             * \param[in] config - Definition of the particular problem.
             */
            void SetColorGrid_Parallel(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the rotational velocity at each node.
             * \param[in] config - Definition of the particular problem.
             */
            void SetRotationalVelocity(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the grid velocity via finite differencing at each node.
             * \param[in] config - Definition of the particular problem.
             */
            void SetGridVelocity(TBOX::TBOX_Config *config, unsigned long iter);

            /*!
             * \brief Perform the MPI communication for the grid coordinates (dynamic meshes).
             * \param[in] config - Definition of the particular problem.
             */
            void Set_MPI_Coord(TBOX::TBOX_Config *config);

            /*!
             * \brief Perform the MPI communication for the grid velocities.
             * \param[in] config - Definition of the particular problem.
             */
            void Set_MPI_GridVel(TBOX::TBOX_Config *config);

            /*!
             * \brief Set the periodic boundary conditions.
             * \param[in] config - Definition of the particular problem.
             */
            void SetPeriodicBoundary(TBOX::TBOX_Config *config);

            /*!
             * \brief Do an implicit smoothing of the grid coordinates.
             * \param[in] val_nSmooth - Number of smoothing iterations.
             * \param[in] val_smooth_coeff - Relaxation factor.
             * \param[in] config - Definition of the particular problem.
             */
            void SetCoord_Smoothing(unsigned short val_nSmooth, double val_smooth_coeff, TBOX::TBOX_Config *config);

            /*!
             * \brief Write the .su2 file.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_out_filename - Name of the output file.
             */
            void SetMeshFile(TBOX::TBOX_Config *config, std::string val_mesh_out_filename);

            /*!
             * \brief Compute some parameters about the grid quality.
             * \param[out] statistics - Information about the grid quality, statistics[0] = (r/R)_min, statistics[1] = (r/R)_ave.
             */
            void GetQualityStatistics(double *statistics);

            /*!
             * \brief Find and store all vertices on a sharp corner in the geometry.
             * \param[in] config - Definition of the particular problem.
             */
            void ComputeSurf_Curvature(TBOX::TBOX_Config *config);

            /*!
             * \brief Find and store the closest neighbor to a vertex.
             * \param[in] config - Definition of the particular problem.
             */
            void FindNormal_Neighbor(TBOX::TBOX_Config *config);

            /*!
             * \brief Retrieve total number of nodes in a simulation across all processors (including halos).
             * \returns Total number of nodes in a simulation across all processors (including halos).
             */
            unsigned long GetGlobal_nPoint();

            /*!
             * \brief Retrieve total number of nodes in a simulation across all processors (excluding halos).
             * \returns Total number of nodes in a simulation across all processors (excluding halos).
             */
            unsigned long GetGlobal_nPointDomain();

            /*!
             * \brief Retrieve total number of elements in a simulation across all processors.
             * \returns Total number of elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElem();

            /*!
             * \brief Retrieve total number of triangular elements in a simulation across all processors.
             * \returns Total number of line elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemLine();

            /*!
             * \brief Retrieve total number of triangular elements in a simulation across all processors.
             * \returns Total number of triangular elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemTria();

            /*!
             * \brief Retrieve total number of quadrilateral elements in a simulation across all processors.
             * \returns Total number of quadrilateral elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemQuad();

            /*!
             * \brief Retrieve total number of tetrahedral elements in a simulation across all processors.
             * \returns Total number of tetrahedral elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemTetr();

            /*!
             * \brief Retrieve total number of hexahedral elements in a simulation across all processors.
             * \returns Total number of hexahedral elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemHexa();

            /*!
             * \brief Retrieve total number of prism elements in a simulation across all processors.
             * \returns Total number of prism elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemPris();

            /*!
             * \brief Retrieve total number of pyramid elements in a simulation across all processors.
             * \returns Total number of pyramid elements in a simulation across all processors.
             */
            unsigned long GetGlobal_nElemPyra();

            /*!
             * \brief Get number of triangular elements.
             * \return Number of line elements.
             */
            unsigned long GetnElemLine();

            /*!
             * \brief Get number of triangular elements.
             * \return Number of triangular elements.
             */
            unsigned long GetnElemTria();

            /*!
             * \brief Get number of quadrilateral elements.
             * \return Number of quadrilateral elements.
             */
            unsigned long GetnElemQuad();

            /*!
             * \brief Get number of tetrahedral elements.
             * \return Number of tetrahedral elements.
             */
            unsigned long GetnElemTetr();

            /*!
             * \brief Get number of hexahedral elements.
             * \return Number of hexahedral elements.
             */
            unsigned long GetnElemHexa();

            /*!
             * \brief Get number of prism elements.
             * \return Number of prism elements.
             */
            unsigned long GetnElemPris();

            /*!
             * \brief Get number of pyramid elements.
             * \return Number of pyramid elements.
             */
            unsigned long GetnElemPyra();

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

            /*!
             * \brief Read the sensitivity from an input file.
             * \param[in] config - Definition of the particular problem.
             */
            void SetBoundSensitivity(TBOX::TBOX_Config *config);

 

        };
    }
}

#include "GEOM_GeometryPhysical.inl"

#endif
