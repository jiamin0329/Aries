

#ifndef ARIES_GRID_VOLUMETRICMOVEMENT_HPP
#define ARIES_GRID_VOLUMETRICMOVEMENT_HPP


#include "GRID_Gridmovement.hpp"
#include "../Math/Math_Matrix.hpp"
#include "../Math/Math_Vector.hpp"
#include "../Math/Math_LinearSolver.hpp"

namespace ARIES
{
    namespace GRID
    {
        /*!
        * \class GRID_VolumetricMovement
        * \brief Class for moving the volumetric numerical grid.
        * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
        * \version 3.2.9 "eagle"
        */
        class GRID_VolumetricMovement : public GRID_Gridmovement {
        protected:

            unsigned short nDim;		/*!< \brief Number of dimensions. */
            unsigned short nVar;		/*!< \brief Number of variables. */

            unsigned long nPoint;		/*!< \brief Number of points. */
            unsigned long nPointDomain;		/*!< \brief Number of points in the domain. */

            MATH::MATH_Matrix StiffMatrix; /*!< \brief Matrix to store the point-to-point stiffness. */
            MATH::MATH_Vector LinSysSol;
            MATH::MATH_Vector LinSysRes;

        public:

            /*!
            * \brief Constructor of the class.
            */
            GRID_VolumetricMovement(GEOM::GEOM_Geometry *geometry);

            /*!
            * \brief Destructor of the class.
            */
            ~GRID_VolumetricMovement(void);

            /*!
            * \brief Update the value of the coordinates after the grid movement.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void UpdateGridCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Update the dual grid after the grid movement (edges and control volumes).
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void UpdateDualGrid(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Update the coarse multigrid levels after the grid movement.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void UpdateMultiGrid(GEOM::GEOM_Geometry **geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Compute the stiffness matrix for grid deformation using spring analogy.
            * \param[in] geometry - Geometrical definition of the problem.
            * \return Value of the length of the smallest edge of the grid.
            */
            double SetFEAMethodContributions_Elem(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Build the stiffness matrix for a 3-D hexahedron element. The result will be placed in StiffMatrix_Elem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
            * \param[in] CoordCorners - Index value for Node 1 of the current hexahedron.
            */
            void SetFEA_StiffMatrix3D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, double **StiffMatrix_Elem, unsigned long PointCorners[8], double CoordCorners[8][3], unsigned short nNodes, double scale);

            /*!
            * \brief Build the stiffness matrix for a 3-D hexahedron element. The result will be placed in StiffMatrix_Elem.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
            * \param[in] CoordCorners - Index value for Node 1 of the current hexahedron.
            */
            void SetFEA_StiffMatrix2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, double **StiffMatrix_Elem, unsigned long PointCorners[8], double CoordCorners[8][3], unsigned short nNodes, double scale);

            /*!
            * \brief Shape functions and derivative of the shape functions
            * \param[in] Xi - Local coordinates.
            * \param[in] Eta - Local coordinates.
            * \param[in] Mu - Local coordinates.
            * \param[in] CoordCorners - Coordiantes of the corners.
            * \param[in] DShapeFunction - Shape function information
            */
            double ShapeFunc_Hexa(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]);

            /*!
            * \brief Shape functions and derivative of the shape functions
            * \param[in] Xi - Local coordinates.
            * \param[in] Eta - Local coordinates.
            * \param[in] Mu - Local coordinates.
            * \param[in] CoordCorners - Coordiantes of the corners.
            * \param[in] DShapeFunction - Shape function information
            */
            double ShapeFunc_Tetra(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]);

            /*!
            * \brief Shape functions and derivative of the shape functions
            * \param[in] Xi - Local coordinates.
            * \param[in] Eta - Local coordinates.
            * \param[in] Mu - Local coordinates.
            * \param[in] CoordCorners - Coordiantes of the corners.
            * \param[in] DShapeFunction - Shape function information
            */
            double ShapeFunc_Pyram(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]);

            /*!
            * \brief Shape functions and derivative of the shape functions
            * \param[in] Xi - Local coordinates.
            * \param[in] Eta - Local coordinates.
            * \param[in] Mu - Local coordinates.
            * \param[in] CoordCorners - Coordiantes of the corners.
            * \param[in] DShapeFunction - Shape function information
            */
            double ShapeFunc_Prism(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]);

            /*!
            * \brief Shape functions and derivative of the shape functions
            * \param[in] Xi - Local coordinates.
            * \param[in] Eta - Local coordinates.
            * \param[in] CoordCorners - Coordiantes of the corners.
            * \param[in] DShapeFunction - Shape function information
            */
            double ShapeFunc_Triangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]);

            /*!
            * \brief Shape functions and derivative of the shape functions
            * \param[in] Xi - Local coordinates.
            * \param[in] Eta - Local coordinates.
            * \param[in] CoordCorners - Coordiantes of the corners.
            * \param[in] DShapeFunction - Shape function information
            */
            double ShapeFunc_Rectangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]);

            /*!
            * \brief Compute the shape functions for hexahedron
            * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
            */
            double GetHexa_Volume(double CoordCorners[8][3]);

            /*!
            * \brief Compute the shape functions for hexahedron
            * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
            */
            double GetTetra_Volume(double CoordCorners[8][3]);

            /*!
            * \brief Compute the shape functions for hexahedron
            * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
            */
            double GetPrism_Volume(double CoordCorners[8][3]);

            /*!
            * \brief Compute the shape functions for hexahedron
            * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
            */
            double GetPyram_Volume(double CoordCorners[8][3]);

            /*!
            * \brief Compute the shape functions for hexahedron
            * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
            */
            double GetTriangle_Area(double CoordCorners[8][3]);

            /*!
            * \brief Compute the shape functions for hexahedron
            * \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
            */
            double GetRectangle_Area(double CoordCorners[8][3]);

            /*!
            * \brief Add the stiffness matrix for a 2-D triangular element to the global stiffness matrix for the entire mesh (node-based).
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
            * \param[in] PointCorners
            * \param[in] nNodes
            */
            void AddFEA_StiffMatrix(GEOM::GEOM_Geometry *geometry, double **StiffMatrix_Elem, unsigned long PointCorners[8], unsigned short nNodes);

            /*!
            * \brief Check for negative volumes (all elements) after performing grid deformation.
            * \param[in] geometry - Geometrical definition of the problem.
            */
            double Check_Grid(GEOM::GEOM_Geometry *geometry);

            /*!
            * \brief Compute the minimum distance to the nearest deforming surface.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void ComputeDeforming_Wall_Distance(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Check the boundary vertex that are going to be moved.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void SetBoundaryDisplacements(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Check the domain points vertex that are going to be moved.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            */
            void SetDomainDisplacements(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
            * \brief Unsteady grid movement using rigid mesh rotation.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iZone - Zone number in the mesh.
            * \param[in] iter - Physical time iteration number.
            */
            void Rigid_Rotation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter);

            /*!
            * \brief Unsteady pitching grid movement using rigid mesh motion.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iZone - Zone number in the mesh.
            * \param[in] iter - Physical time iteration number.
            */
            void Rigid_Pitching(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter);

            /*!
            * \brief Unsteady plunging grid movement using rigid mesh motion.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iZone - Zone number in the mesh.
            * \param[in] iter - Physical time iteration number.
            */
            void Rigid_Plunging(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter);

            /*!
            * \brief Unsteady translational grid movement using rigid mesh motion.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iZone - Zone number in the mesh.
            * \param[in] iter - Physical time iteration number.
            */
            void Rigid_Translation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter);

            /*!
            * \brief Grid deformation using the spring analogy method.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] UpdateGeo - Update geometry.
            */
            void SetVolume_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, bool UpdateGeo);

            /*!
            * \brief Compute the determinant of a 3 by 3 matrix.
            * 3 by 3 matrix elements
            * \param[in] A00
            * \param[in] A01
            * \param[in] A02
            * \param[in] A10
            * \param[in] A11
            * \param[in] A12
            * \param[in] A20
            * \param[in] A21
            * \param[in] A22
            * \result Determinant of the matrix
            */
            double Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22);

        };
    }
}

#include "GRID_VolumetricMovement.inl"

#endif