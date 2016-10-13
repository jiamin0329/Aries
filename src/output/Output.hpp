/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    class for data output
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_OUTPUT_HPP
#define ARIES_OUTPUT_HPP

/* external includes */
#ifdef HAVE_MPI
#include "mpi.h"
#endif
#ifdef HAVE_CGNS
#include "cgnslib.h"
#endif
#ifdef HAVE_TECIO
#include "TECIO.h"
#endif

/* c++ includes */
#include <fstream>
#include <cmath>
#include <time.h>

/* Aries includes */
#include "IMesh.hpp"
#include "IProcData.hpp"
#include "ISolution.hpp"

/*
 *================================================================================
 *                            Forward Declaration
 *================================================================================
 */


using namespace std;

/*
 *================================================================================
 *                             Namespace Class
 *================================================================================
 */

namespace ARIES
{
    class Output
    {
    public:
        /*!
         * \brief Constructor of the class.
         */
        Output();
		  
        /*!
         * \brief Destructor of the class.
         */
        ~Output();

        /*!
         * \brief Writes and organizes the all the output files, except the history one, for serial computations.
         * \param[in] sol - Container vector with all the solutions.
         * \param[in] mesh - Geometrical definition of the problem.
         * \param[in] procData - processor data 
         * \param[in] iExtIter - Current external (time) iteration.
         * \param[in] val_nZone - Total number of domains in the grid file.
         */
        void SetResultFiles(ISolution* sol, IMesh* mesh, IProcData* procData, unsigned long iExtIter, unsigned short val_nZone);

        /*!
         * \brief Create and write the file with the flow coefficient on the surface.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] FlowSolution - Flow solution.
         * \param[in] iExtIter - Current external (time) iteration.
         * \param[in] val_iZone - Current zone number in the grid file.
         */
        bool SetSurfaceCsvFlow(ISolution* sol, IMesh* mesh, IProcData* procData, unsigned long iExtIter, unsigned short val_iZone);

        /*!
         * \brief Create and write the file with the adjoint coefficients on the surface for serial computations.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] AdjSolution - Adjoint solution.
         * \param[in] FlowSolution - Flow solution.
         * \param[in] iExtIter - Current external (time) iteration.
         * \param[in] val_iZone - Current zone number in the grid file.
         */
        bool SetSurfaceCsvAdjoint(ISolution* sol, IMesh* mesh, IProcData* procData, unsigned long iExtIter, unsigned short val_iZone);

        /*!
         * \brief Create and write the file with linearized coefficient on the surface for serial computations
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] LinSolution - Linearized solution.
         * \param[in] val_filename - Name of the output file.
         * \param[in] iExtIter - Current external (time) iteration.
         */
        bool SetSurfaceCsvLinearized(ISolution* sol, IMesh* mesh, IProcData* procData, string val_filename, unsigned long iExtIter); 

        /*!
         * \brief Merge the geometry into a data structure used for output file writing.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] val_nZone - iZone index.
         */
        void MergeConnectivity(IProcData* procData, IMesh* mesh, unsigned short val_iZone);

        /*!
         * \brief Merge the connectivity for a single element type from all processors.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] Elem_Type - VTK index of the element type being merged.
         */
        void MergeVolumetricConnectivity(IProcData* procData, IMesh* mesh, unsigned short Elem_Type);

        /*!
         * \brief Merge the connectivity for a single element type from all processors.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] Elem_Type - VTK index of the element type being merged.
         */
        void MergeSurfaceConnectivity(IProcData* procData, IMesh* mesh, unsigned short Elem_Type);


        
        /*!
         * \brief Merge the node coordinates from all processors.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         */
        void MergeCoordinates(IProcData* procData, IMesh* mesh);
        
        
        /*!
         * \brief Merge the solution into a data structure used for output file writing.
         * \param[in] config - Definition of the particular problem.
         * \param[in] geometry - Geometrical definition of the problem.
         * \param[in] solution - Flow, adjoint or linearized solution.
         * \param[in] val_nZone - iZone index.
         */
        void MergeSolution(IProcData* procData, IMesh* mesh, ISolution* sol, unsigned short val_iZone);
        
        
            
        //            /*!
        //             * \brief Writes and organizes the all the output files, except the history one, for serial computations.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             * \param[in] val_iZone - Total number of domains in the grid file.
        //             * \param[in] val_nZone - Total number of domains in the grid file.
        //             */
        //            void SetBaselineResult_Files(SOLV::SOLV_Solver **solver, GEOM::GEOM_Geometry **geometry, TBOX::TBOX_Config **config,
        //                                         unsigned long iExtIter, unsigned short val_nZone);
        //
        //            /*!
        //             * \brief Writes and organizes the all the output files, except the history one, for serial computations.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] val_nZone - Total number of domains in the grid file.
        //             */
        //            void SetMesh_Files(GEOM::GEOM_Geometry **geometry, TBOX::TBOX_Config **config, unsigned short val_nZone, bool new_file, bool su2_file);
        //
        //            /*!
        //             * \brief Writes equivalent area.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             */
        //            void SetEquivalentArea(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
        //                                   unsigned long iExtIter);
        //
        //            /*!
        //             * \brief Writes inverse design.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             */
        //            void SetCp_InverseDesign(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
        //                                     unsigned long iExtIter);
        //
        //            /*!
        //             * \brief Writes inverse design.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             */
        //            void SetHeat_InverseDesign(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
        //                                       unsigned long iExtIter);
        //
        //            /*!
        //             * \brief Writes forces at different sections.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             */
        //            void SetForceSections(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
        //                                  unsigned long iExtIter);
        //
        //            /*!
        //             * \brief Writes one dimensional output.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             */
        //            void OneDimensionalOutput(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);
        //
        //            /*!
        //             * \brief Writes mass flow rate output at monitored marker.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             */
        //            void SetMassFlowRate(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);





        //            /*!
        //             * \brief Merge the solution into a data structure used for output file writing.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] solution - Flow, adjoint or linearized solution.
        //             * \param[in] val_nZone - iZone index.
        //             */
        //            void MergeBaselineSolution(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, SOLV::SOLV_Solver *solver, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write a native SU2 restart file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetRestart(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, SOLV::SOLV_Solver **solver, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write the x, y, & z coordinates to a CGNS output file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetCGNS_Coordinates(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write the element connectivity to a CGNS output file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetCGNS_Connectivity(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write solution data to a CGNS output file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetCGNS_Solution(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write a Paraview ASCII solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - Current zone.
        //             * \param[in] val_nZone - Total number of zones.
        //             */
        //            void SetParaview_ASCII(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
        //
        //            /*!
        //             * \brief Write a Paraview ASCII solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - Current zone.
        //             * \param[in] val_nZone - Total number of zones.
        //             */
        //            void SetParaview_MeshASCII(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol, bool new_file);
        //
        //            /*!
        //             * \brief Write a Tecplot ASCII solution file.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             */
        //            void SetTecplotASCII_LowMemory(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, SOLV::SOLV_Solver **solver, char mesh_filename[TBOX::MAX_STRING_SIZE], bool surf_sol);
        //
        //            /*!
        //             * \brief Write a Tecplot ASCII solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - Current zone.
        //             * \param[in] val_nZone - Total number of zones.
        //             */
        //            void SetTecplotASCII(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, SOLV::SOLV_Solver **solver, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetTecplotASCII_Mesh(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, bool surf_sol, bool new_file);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            std::string AssembleVariableNames(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short nVar_Consv, unsigned short *NVar);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetSU2_MeshASCII(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetSU2_MeshBinary(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetTecplotBinary_DomainMesh(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write the coordinates and connectivity to a Tecplot binary surface mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetTecplotBinary_SurfaceMesh(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write solution data to a Tecplot binary volume solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetTecplotBinary_DomainSolution(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write solution data to a Tecplot binary surface solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetTecplotBinary_SurfaceSolution(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write a Tecplot ASCII solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - Current zone.
        //             * \param[in] val_nZone - Total number of zones.
        //             */
        //            void SetFieldViewASCII(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetFieldViewASCII_Mesh(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry);
        //
        //            /*!
        //             * \brief Write the nodal coordinates and connectivity to a Tecplot binary mesh file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetFieldViewBinary_Mesh(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry);
        //
        //            /*!
        //             * \brief Write solution data to a Tecplot binary volume solution file.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] val_iZone - iZone index.
        //             */
        //            void SetFieldViewBinary(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone);
        //
        //            /*!
        //             * \brief Deallocate temporary memory needed for merging and writing coordinates.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             */
        //            void DeallocateCoordinates(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry);
        //
        //            /*!
        //             * \brief Deallocate temporary memory needed for merging and writing connectivity.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             */
        //            void DeallocateConnectivity(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry, bool surf_sol);
        //
        //            /*!
        //             * \brief Deallocate temporary memory needed for merging and writing solution variables.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             */
        //            void DeallocateSolution(TBOX::TBOX_Config *config, GEOM::GEOM_Geometry *geometry);
        //
        //            /*!
        //             * \brief Write the header of the history file.
        //             * \param[in] ConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
        //             * \param[in] config - Definition of the particular problem.
        //             */
        //            void SetConvHistory_Header(ofstream *ConvHist_file, TBOX::TBOX_Config *config);
        //
        //            /*!
        //             * \brief Write the history file and the convergence on the screen for serial computations.
        //             * \param[in] ConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] integration - Generic subroutines for space integration, time integration, and monitoring.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             * \param[in] timeused - Current number of clock tick in the computation (related with total time).
        //             * \param[in] val_nZone - iZone index.
        //             */
        //            void SetConvHistory_Body(ofstream *ConvHist_file, GEOM::GEOM_Geometry ***geometry, SOLV::SOLV_Solver ****solver_container, TBOX::TBOX_Config **config,
        //                                     INTE::INTE_Integration ***integration, bool DualTime, double timeused, unsigned short val_iZone);
        //
        //            /*!
        //             * \brief Write the history file and the convergence on the screen for serial computations.
        //             * \param[in] ConvHist_file - Pointer to the convergence history file (which is defined in the main subroutine).
        //             * \param[in] geometry - Geometrical definition of the problem.
        //             * \param[in] solver_container - Container vector with all the solutions.
        //             * \param[in] config - Definition of the particular problem.
        //             * \param[in] integration - Generic subroutines for space integration, time integration, and monitoring.
        //             * \param[in] iExtIter - Current external (time) iteration.
        //             * \param[in] timeused - Current number of clock tick in the computation (related with total time).
        //             * \param[in] val_nZone - iZone index.
        //             */
        //            void SetForces_Breakdown(GEOM::GEOM_Geometry ***geometry, SOLV::SOLV_Solver ****solver_container, TBOX::TBOX_Config **config,
        //                                     INTE::INTE_Integration ***integration, unsigned short val_iZone);
        //
        //            void SetCFL_Number(SOLV::SOLV_Solver ****solver_container, TBOX::TBOX_Config **config, unsigned short val_iZone);

    private:
        unsigned short **nOutput_Vars;
        double ****data_container;
            
        unsigned long d_numGlobalPointDomain;          // Global number of domains
        unsigned long d_numGlobalPoints;           // Global number of nodes with halos
        unsigned long d_numSurfPoints;             // Global number of nodes of the surface
        unsigned long d_numSurfElems;              // Global number of surface elems without halos
        unsigned long d_numGlobalElems;            // Global number of elems without halos
        unsigned long d_numGlobalBoundLines;
        unsigned long d_numGlobalBoundTrias;
        unsigned long d_numGlobalBoundQuads;
        unsigned long d_numGlobalTrias;
        unsigned long d_numGlobalQuads;
        unsigned long d_numGlobalTetrs;
        unsigned long d_numGlobalHexas;
        unsigned long d_numGlobalPriss;
        unsigned long d_numGlobalPyras;
            
        double **d_coords; // node i (x, y, z) = (Coords[0][i], Coords[1][i], Coords[2][i])
            
        int *d_connBoundLine;
        int *d_connBoundTria;
        int *d_connBoundQuad;
        int *d_connTria;	// triangle 1 = Conn_Tria[0], Conn_Tria[1], Conn_Tria[3]
        int *d_connQuad;
        int *d_connTetr;
        int *d_connHexa;
        int *d_connPris;
        int *d_connPyra;
            
        double *Volume;
        double **Data;
        double **residuals, **consv_vars;					// placeholders
        double *p, *rho, *M, *Cp, *Cf, *Ch, *h, *yplus;		// placeholders 
        unsigned short nVar_Consv, nVar_Total, nVar_Extra, nZones;

        // output flags
        bool d_isBaseOutput;
        bool d_isSurfOutput;
        bool d_isCgnsOutput;
        bool d_isTecplotOutput;
        bool d_isParaviewOutput;

            
        double d_rhoResNew;
        double d_rhoResOld;

        /*
         *  cgns flags
         */
        int cgns_base, cgns_zone, cgns_base_results, cgns_zone_results;
    }; // end Output class
} // end ARIES namespace

#endif










