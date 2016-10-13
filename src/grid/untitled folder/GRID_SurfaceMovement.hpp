
#ifndef ARIES_GRID_SURFACEMOVEMENT_HPP
#define ARIES_GRID_SURFACEMOVEMENT_HPP

#include "GRID_Gridmovement.hpp"
#include "GRID_FreeFormDefBox.hpp"

namespace ARIES
{
    namespace GRID
    {
        /*!
        * \class GRID_SurfaceMovement
        * \brief Class for moving the surface numerical grid.
        * \author F. Palacios, T. Economon.
        * \version 3.2.9 "eagle"
        */
        class GRID_SurfaceMovement : public GRID_Gridmovement 
        {
        protected:
            GRID_FreeFormDefBox** FFDBox;	/*!< \brief Definition of the Free Form Deformation Box. */
            unsigned short nFFDBox;	/*!< \brief Number of FFD FFDBoxes. */
            unsigned short nLevel;	/*!< \brief Level of the FFD FFDBoxes (parent/child). */
            bool FFDBoxDefinition;	/*!< \brief If the FFD FFDBox has been defined in the input file. */
            std::vector<double> GlobalCoordX[TBOX::MAX_NUMBER_FFD];
            std::vector<double> GlobalCoordY[TBOX::MAX_NUMBER_FFD];
            std::vector<double> GlobalCoordZ[TBOX::MAX_NUMBER_FFD];
            std::vector<std::string> GlobalTag[TBOX::MAX_NUMBER_FFD];
            std::vector<unsigned long> GlobalPoint[TBOX::MAX_NUMBER_FFD];

        public:

            /*!
            * \brief Constructor of the class.
            */
            GRID_SurfaceMovement(void);

            /*!
            * \brief Destructor of the class.
            */
            ~GRID_SurfaceMovement(void);

            /*!
            * \brief Set a Hicks-Henne deformation bump functions on an airfoil.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetHicksHenne(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a NACA 4 digits airfoil family for airfoil deformation.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            */
            void SetNACA_4Digits(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config);

            /*!
            * \brief Set a parabolic family for airfoil deformation.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            */
            void SetParabolic(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config);

            /*!
            * \brief Set a obstacle in a channel.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            */
            void SetAirfoil(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config);

            /*!
            * \brief Set a rotation for surface movement.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetRotation(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set the translational/rotational velocity for a moving wall.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iZone - Zone number in the mesh.
            * \param[in] iter - Physical time iteration number.
            */
            void Moving_Walls(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter);

            /*!
            * \brief Computes the displacement of a translating surface for a dynamic mesh simulation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iter - Current physical time iteration.
            * \param[in] iZone - Zone number in the mesh.
            */
            void Surface_Translating(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
                unsigned long iter, unsigned short iZone);

            /*!
            * \brief Computes the displacement of a plunging surface for a dynamic mesh simulation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iter - Current physical time iteration.
            * \param[in] iZone - Zone number in the mesh.
            */
            void Surface_Plunging(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
                unsigned long iter, unsigned short iZone);

            /*!
            * \brief Computes the displacement of a pitching surface for a dynamic mesh simulation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iter - Current physical time iteration.
            * \param[in] iZone - Zone number in the mesh.
            */
            void Surface_Pitching(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
                unsigned long iter, unsigned short iZone);

            /*!
            * \brief Computes the displacement of a rotating surface for a dynamic mesh simulation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iter - Current physical time iteration.
            * \param[in] iZone - Zone number in the mesh.
            */
            void Surface_Rotating(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
                unsigned long iter, unsigned short iZone);

            /*!
            * \brief Unsteady aeroelastic grid movement by deforming the mesh.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] ExtIter - Physical iteration number.
            * \param[in] iMarker - Marker to deform.
            * \param[in] iMarker_Monitoring - Marker we are monitoring.
            * \param[in] displacements - solution of typical section wing model.
            */
            void AeroelasticDeform(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned long ExtIter, unsigned short iMarker, unsigned short iMarker_Monitoring, double displacements[4]);

            /*!
            * \brief Deforms a 3-D flutter/pitching surface during an unsteady simulation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iter - Current physical time iteration.
            * \param[in] iZone - Zone number in the mesh.
            */
            void SetBoundary_Flutter3D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
                GRID_FreeFormDefBox **FFDBox, unsigned long iter, unsigned short iZone);

            /*!
            * \brief Set the collective pitch for a blade surface movement.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void SetCollective_Pitch(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Set any surface deformationsbased on an input file.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iZone - Zone number in the mesh.
            * \param[in] iter - Current physical time iteration.
            */
            void SetExternal_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter);

            /*!
            * \brief Set a displacement for surface movement.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetTranslation(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a displacement for surface movement.
            * \param[in] boundary - Geometry of the boundary.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetScale(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Copy the boundary coordinates to each vertex.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void CopyBoundary(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Set the surface/boundary deformation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void SetSurface_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Compute the parametric coordinates of a grid point using a point inversion strategy
            *        in the free form FFDBox.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            */
            void SetParametricCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox);

            /*!
            * \brief Update the parametric coordinates of a grid point using a point inversion strategy
            *        in the free form FFDBox.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iFFDBox - _____________________.
            */
            void UpdateParametricCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox);

            /*!
            * \brief Check the intersections of the FFD with the surface
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iFFDBox - _____________________.
            */
            void CheckFFDIntersections(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox);

            /*!
            * \brief _____________________.
            * \param[in] geometry - _____________________.
            * \param[in] config - _____________________.
            * \param[in] FFDBoxParent - _____________________.
            * \param[in] FFDBoxChild - _____________________.
            */
            void SetParametricCoordCP(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBoxParent, GRID_FreeFormDefBox *FFDBoxChild);

            /*!
            * \brief _____________________.
            * \param[in] geometry - _____________________.
            * \param[in] config - _____________________.
            * \param[in] FFDBoxParent - _____________________.
            * \param[in] FFDBoxChild - _____________________.
            */
            void GetCartesianCoordCP(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBoxParent, GRID_FreeFormDefBox *FFDBoxChild);

            /*!
            * \brief Recompute the cartesian coordinates using the control points position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iFFDBox - _____________________.
            */
            void SetCartesianCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox);

            /*!
            * \brief Set the deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDCPChange_2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set the deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDCPChange(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a camber deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDCamber_2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a thickness deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDThickness_2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a camber deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDCamber(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a thickness deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDThickness(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a dihedral angle deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDDihedralAngle(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a twist angle deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDTwistAngle(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a rotation angle deformation of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDRotation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Set a rotation angle deformation in a control surface of the Free From box using the control point position.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] iDV - Index of the design variable.
            * \param[in] ResetDef - Reset the deformation before starting a new one.
            */
            void SetFFDControl_Surface(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iDV, bool ResetDef);

            /*!
            * \brief Read the free form information from the grid input file.
            * \note If there is no control point information, and no parametric
            *       coordinates information, the code will compute that information.
            * \param[in] config - Definition of the particular problem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            * \param[in] val_mesh_filename - Name of the grid input file.
            */
            void ReadFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox **FFDBox, std::string val_mesh_filename);

            /*!
            * \brief Read the free form information from the grid input file.
            * \note If there is no control point information, and no parametric
            *       coordinates information, the code will compute that information.
            * \param[in] config - Definition of the particular problem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] FFDBox - Array with all the free forms FFDBoxes of the computation.
            */
            void ReadFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox **FFDBox);

            /*!
            * \brief Merge the Free Form information in the SU2 file.
            * \param[in] config - Definition of the particular problem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] val_mesh_filename - Name of the grid output file.
            */
            void MergeFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Write the Free Form information in the SU2 file.
            * \param[in] config - Definition of the particular problem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] val_mesh_filename - Name of the grid output file.
            */
            void WriteFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Get information about if there is a complete FFDBox definition, or it is necessary to
            *        compute the parametric coordinates.
            * \return <code>TRUE</code> if the input grid file has a complete information; otherwise <code>FALSE</code>.
            */
            bool GetFFDBoxDefinition(void);

            /*!
            * \brief Obtain the number of FFDBoxes.
            * \return Number of FFD FFDBoxes.
            */
            unsigned short GetnFFDBox(void);

            /*!
            * \brief Obtain the number of levels.
            * \return Number of FFD levels.
            */
            unsigned short GetnLevel(void);

        };
    }
}

#include "GRID_SurfaceMovement.inl"

#endif