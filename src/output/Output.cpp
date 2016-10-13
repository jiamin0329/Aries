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


#include "Output.hpp"

#include "const_def.h"

#include <iostream>
#include <vector>

using namespace std;

namespace ARIES
{
    /*
     *  constructor
     */ 
    Output::Output()
    {
        // Initialize point and connectivity counters to zero
        d_numGlobalPointDomain = 0;
        d_numGlobalPoints = 0;
        d_numSurfPoints = 0;
        d_numGlobalElems = 0;
        d_numSurfElems = 0;
        d_numGlobalTrias = 0;
        d_numGlobalQuads = 0;
        d_numGlobalTetrs = 0;
        d_numGlobalHexas = 0;
        d_numGlobalPriss = 0;
        d_numGlobalPyras = 0;
        d_numGlobalBoundLines = 0;
        d_numGlobalBoundTrias = 0;
        d_numGlobalBoundQuads = 0;
	  
        // Initialize flags
        d_isBaseOutput = false;
        d_isSurfOutput = false;
        d_isCgnsOutput = false;
        d_isParaviewOutput = false;
        d_isTecplotOutput = false;

        // Initialize residual
        d_rhoResNew = EPS;
        d_rhoResOld = EPS;
    }

    Output::~Output(void)
    {
    }

    void Output::SetResultFiles(ISolution* sol, IMesh* mesh, IProcData* procData, unsigned long iExtIter, unsigned short val_nZone)
    {
        int rank = MASTER_NODE;

#ifdef ARIES_HAVE_MPI
        int size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
            
        unsigned short iZone;

        for (iZone = 0; iZone < val_nZone; iZone++)
        {
            /*
             *  Flags identifying the types of files to be written.
             */
            bool isWrtVol  = procData->IsWrtVolSol(iZone);
            bool isWrtSurf = procData->IsWrtSurfSol(iZone);
            bool isWrtCsv  = procData->IsWrtCsvSol(iZone);

            SolverType solverType = procData->GetSolverType(iZone);
            OutputType outputType = procData->GetOutputType(iZone);
                
#ifdef ARIES_HAVE_MPI
            /*
             *  Do not merge the volume solutions if we are running in parallel.
             *  Force the use of SU2_SOL to merge the volume sols in this case.
             */
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            if (size > SINGLE_NODE)
            {
                isWrtVol = false;
                isWrtSurf = false;
            }
#endif
            if (rank == MASTER_NODE)
                cout << endl << "Writing comma-separated values (CSV) surface files." << endl;

            if (isWrtCsv)
            {
                switch (solverType)
                {
                    case EULER:
                    case NAVIER_STOKES:
                    case RANS:
                    case FLUID_STRUCTURE_EULER:
                    case FLUID_STRUCTURE_NAVIER_STOKES:
                    case FLUID_STRUCTURE_RANS:
                    case TNE2_EULER:
                    case TNE2_NAVIER_STOKES:
                        SetSurfaceCsvFlow(sol, mesh, procData, iExtIter, iZone);
                        break;
                    case ADJ_EULER:
                    case ADJ_NAVIER_STOKES:
                    case ADJ_RANS:
                    case ADJ_TNE2_EULER:
                    case ADJ_TNE2_NAVIER_STOKES:
                        SetSurfaceCsvAdjoint(sol, mesh, procData, iExtIter, iZone);
                        break;
                    
                    case LIN_EULER:
                    case LIN_NAVIER_STOKES:
                        SetSurfaceCsvLinearized(sol, mesh, procData, procData->GetSurfLinCoeffFileName(), iExtIter);
                        break;
                    default:
                        SetSurfaceCsvFlow(sol, mesh, procData, iExtIter, iZone);
                        break;
                }
            }

            /*
             *  Merge the node coordinates and connectivity, if necessary. This
             *  is only performed if a volume solution file is requested, and it
             *  is active by default.
             */
            if (isWrtVol || isWrtSurf)
            {
                if (rank == MASTER_NODE)
                    cout << "Merging connectivities in the Master node." << endl;
                MergeConnectivity(procData, mesh, iZone);
            }

            /*
             *  Merge coordinates of all grid nodes (excluding ghost points).
             *  The grid coordinates are always merged and included first in the
             *  restart files.
             */
            if (rank == MASTER_NODE)
                cout << "Merging coordinates in the Master node." << endl;
            MergeCoordinates(procData, mesh);

            if ((rank == MASTER_NODE) && (isWrtVol || isWrtSurf))
            {
                if (outputType == TECPLOT_BINARY)
                {
                    if (rank == MASTER_NODE)
                        cout << "Writing Tecplot binary volume and surface mesh files." << endl;
                    //SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][MESH_0], iZone);
                    //SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][MESH_0], iZone);
                    //if (!d_isBaseOutput)
                    //DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
                    //if (!d_isSurfOutput)
                    //DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], d_isSurfOutput);
                }
            }

            /*
             *  Merge the solution data needed for volume solutions and restarts
             */
            if (rank == MASTER_NODE) 
                cout << "Merging solution in the Master node." << endl;
            MergeSolution(procData, mesh, sol, iZone);

            /*
             *  Write restart, or Tecplot files using the merged data.
             *  This data lives only on the master, and these routines are currently
             *  executed by the master proc alone (as if in serial).
             */
            if (rank == MASTER_NODE)
            {
                /*
                 *  Write a native restart file
                 */
                if (rank == MASTER_NODE) 
                    cout << "Writing SU2 native restart file." << endl;
                //SetRestart(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone);

                if (isWrtVol)
                {
                    switch (outputType)
                    {
						case TECPLOT:		
							if (rank == MASTER_NODE)
								cout << "Writing Tecplot ASCII file volume solution file." << endl;
							//SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone, val_nZone, false);
							//DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
							break;

						case FIELDVIEW:
							if (rank == MASTER_NODE)
								cout << "Writing FieldView ASCII file volume solution file." << endl;
							//SetFieldViewASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
							//DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
							break;

						case TECPLOT_BINARY:
							if (rank == MASTER_NODE)
								cout << "Writing Tecplot binary volume solution file." << endl;
							//SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][MESH_0], iZone);
							break;

						case FIELDVIEW_BINARY:
							if (rank == MASTER_NODE)
								cout << "Writing FieldView binary file volume solution file." << endl;
							//SetFieldViewBinary(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
							//DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
							break;

						case PARAVIEW:
							if (rank == MASTER_NODE)
								cout << "Writing Paraview ASCII volume solution file." << endl;
							//SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
							//DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
							break;

						default:
							break;
                    }

                }

                if (isWrtSurf)
                {
                    switch (outputType)
                    {
						case TECPLOT:
							if (rank == MASTER_NODE)
								cout << "Writing Tecplot ASCII surface solution file." << endl;
							//SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone, val_nZone, true);
							//DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
							break;

						case TECPLOT_BINARY:
							if (rank == MASTER_NODE)
								cout << "Writing Tecplot binary surface solution file." << endl;
							//SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][MESH_0], iZone);
							break;

						case PARAVIEW:
							if (rank == MASTER_NODE)
								cout << "Writing Paraview ASCII surface solution file." << endl;
							//SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
							//DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
							break;

						default:
							break;
                    }
                }

                /*--- Release memory needed for merging the solution data. ---*/
                //DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
                //DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);
            }

            /*
             *  Final broadcast (informing other procs that the base output file was written).
             */
#ifdef ARIES_HAVE_MPI
            MPI_Bcast(&d_isBaseOutput, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Bcast(&d_isSurfOutput, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
        } // end iZone loop 
    } // end SetResultFiles
} // end ARIES namespace

