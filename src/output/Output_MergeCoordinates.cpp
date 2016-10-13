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
            MergeSolution(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone);

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

        
    /*
     * \brief Create and write the file with the flow coefficient on the surface.
     */
    bool Output::SetSurfaceCsvFlow(ISolution *sol, IMesh *mesh, IProcData *procData, unsigned long iExtIter, unsigned short iZone)
    {
        unsigned short iMarker;
        unsigned long iPoint, iVertex, globalIndex;
            
        double pressCoeff = 0.0, skinFrictionCoeff, heatFlux;
        double xCoord = 0.0, yCoord = 0.0,  zCoord = 0.0, mach, pressure;
        char cstr[200];

        SolverType solverType = procData->GetSolverType(iZone);
        UnsteadyType unsteadyType = procData->GetUnsteadyType();
        unsigned short nDim = mesh->GetDim();

        bool isWrtUnsteady = procData->GetWriteUnsteady();
#ifndef AREIS_HAVE_MPI
        char buffer[50];
        ofstream surfFlowFile;

        /*--- Write file name with extension if unsteady ---*/
        strcpy(cstr, procData->GetSurfFlowCoeffFileName().c_str());
		
        if (unsteadyType == TIME_SPECTRAL)
        {
            if (int(iZone) < 10)
                sprintf(buffer, "_0000%d.csv", int(iZone));
            if ((int(iZone) >= 10) && (int(iZone) < 100))
                sprintf(buffer, "_000%d.csv", int(iZone));
            if ((int(iZone) >= 100) && (int(iZone) < 1000))
                sprintf(buffer, "_00%d.csv", int(iZone));
            if ((int(iZone) >= 1000) && (int(iZone) < 10000))
                sprintf(buffer, "_0%d.csv", int(iZone));
            if (int(iZone) >= 10000) sprintf(buffer, "_%d.csv", int(iZone));

        }
        else if (unsteadyType && isWrtUnsteady)
        {
            if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))
                sprintf(buffer, "_0000%d.csv", int(iExtIter));
            if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))
                sprintf(buffer, "_000%d.csv", int(iExtIter));
            if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))
                sprintf(buffer, "_00%d.csv", int(iExtIter));
            if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000))
                sprintf(buffer, "_0%d.csv", int(iExtIter));
            if (int(iExtIter) >= 10000)
                sprintf(buffer, "_%d.csv", int(iExtIter));
        }
        else
            sprintf(buffer, ".csv");

        strcat(cstr, buffer);
        surfFlowFile.precision(15);
        surfFlowFile.open(cstr, ios::out);

        surfFlowFile << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
        if (nDim == 3)
            surfFlowFile << "\"z_coord\", ";
        surfFlowFile << "\"Pressure\", \"Pressure_Coefficient\", ";

        switch (solverType)
        {
            case EULER:
                surfFlowFile << "\"Mach_Number\"" << endl;
                break;
            case NAVIER_STOKES:
            case RANS:
                surfFlowFile << "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl;
                break;
            case TNE2_EULER:
                surfFlowFile << "\"Mach_Number\"" << endl;
                break;
            case TNE2_NAVIER_STOKES:
                surfFlowFile << "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl;
                break;
            default:
                break;
        }

        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->IsMarkerPlotted(iMarker))
            {
                for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                {
                    iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                        
                    globalIndex = mesh->GetNode(iPoint)->GetGlobalIndex();
                    xCoord      = mesh->GetNode(iPoint)->GetCoord(0);
                    yCoord      = mesh->GetNode(iPoint)->GetCoord(1);
                    if (nDim == 3)
                        zCoord = mesh->GetNode(iPoint)->GetCoord(2);

                    /*
                     *  The output should be in inches
                     */
                    if (procData->GetSystemMeasurements() == US)
                    {
                        xCoord *= 12.0;
                        yCoord *= 12.0;
                        if (nDim == 3)
                            zCoord *= 12.0;
                    }
                        
                    pressure = sol->GetVariable(iPoint)->GetPressure();
                    pressCoeff = sol->GetCoeffPressure(iMarker, iVertex);
                    surfFlowFile << scientific << globalIndex << ", " << xCoord << ", " << yCoord << ", ";
                    if (nDim == 3)
                        surfFlowFile << scientific << zCoord << ", ";
                    surfFlowFile << scientific << pressure << ", " << pressCoeff << ", ";

                    switch (solverType)
                    {
                        case EULER:
                            mach = sqrt(sol->GetVariable(iPoint)->GetVelocity2()) / sol->GetVariable(iPoint)->GetSoundSpeed();
                            surfFlowFile << scientific << mach << endl;
                            break;
                        case NAVIER_STOKES:
                        case RANS:
                            skinFrictionCoeff = sol->GetCoeffSkinFriction(iMarker, iVertex);
                            heatFlux = sol->GetHeatFlux(iMarker, iVertex);
                            surfFlowFile << scientific << skinFrictionCoeff << ", " << heatFlux << endl;
                            break;
                        case TNE2_EULER:
                            mach = sqrt(sol->GetVariable(iPoint)->GetVelocity2()) / sol->GetVariable(iPoint)->GetSoundSpeed();
                            surfFlowFile << scientific << mach << endl;
                            break;
                        case TNE2_NAVIER_STOKES:
                            skinFrictionCoeff = sol->GetCoeffSkinFriction(iMarker, iVertex);
                            heatFlux = sol->GetHeatFlux(iMarker, iVertex);
                            surfFlowFile << scientific << skinFrictionCoeff << ", " << heatFlux << endl;
                            break;
                        default:
                            break;
                    }
                }
            }
        }

        surfFlowFile.close();
#else
        int rank, iProcessor, nProcessor;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

        unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
        unsigned long nVertex_Surface = 0, nLocalVertex_Surface = 0;
        unsigned long MaxLocalVertex_Surface = 0;

        /*--- Find the max number of surface vertices among all
          partitions and set up buffers. The master node will handle the
          writing of the CSV file after gathering all of the data. ---*/

        nLocalVertex_Surface = 0;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Plotting(iMarker) == YES)
                for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
                }

        /*--- Communicate the number of local vertices on each partition
          to the master node ---*/

        Buffer_Send_nVertex[0] = nLocalVertex_Surface;
        if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long[nProcessor];

        MPI_Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*--- Send and Recv buffers ---*/

        double *Buffer_Send_Coord_x = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_Coord_x = NULL;

        double *Buffer_Send_Coord_y = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_Coord_y = NULL;

        double *Buffer_Send_Coord_z = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_Coord_z = NULL;

        double *Buffer_Send_Press = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_Press = NULL;

        double *Buffer_Send_CPress = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_CPress = NULL;

        double *Buffer_Send_Mach = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_Mach = NULL;

        double *Buffer_Send_SkinFriction = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_SkinFriction = NULL;

        double *Buffer_Send_HeatTransfer = new double[MaxLocalVertex_Surface];
        double *Buffer_Recv_HeatTransfer = NULL;

        unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalVertex_Surface];
        unsigned long *Buffer_Recv_GlobalIndex = NULL;

        /*--- Prepare the receive buffers on the master node only. ---*/

        if (rank == MASTER_NODE) {
            Buffer_Recv_Coord_x = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_Coord_y = new double[nProcessor*MaxLocalVertex_Surface];
            if (nDim == 3) Buffer_Recv_Coord_z = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_Press = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_CPress = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_Mach = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_SkinFriction = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_HeatTransfer = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalVertex_Surface];
        }

        /*--- Loop over all vertices in this partition and load the
          data of the specified type into the buffer to be sent to
          the master node. ---*/

        nVertex_Surface = 0;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Plotting(iMarker) == YES)
                for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    if (geometry->node[iPoint]->GetDomain()) {
                        Buffer_Send_Press[nVertex_Surface] = sol->GetVariable(iPoint)->GetPressure();
                        Buffer_Send_CPress[nVertex_Surface] = FlowSolver->GetCPressure(iMarker, iVertex);
                        Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
                        Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
                        if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2); }

                        /*--- If US system, the output should be in inches ---*/

                        if (config->GetSystemMeasurements() == US) {
                            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
                            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
                            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
                        }

                        Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();

                        if (solver == EULER)
                            Buffer_Send_Mach[nVertex_Surface] = sqrt(sol->GetVariable(iPoint)->GetVelocity2()) / sol->GetVariable(iPoint)->GetSoundSpeed();
                        if ((solver == NAVIER_STOKES) || (solver == RANS))
                            Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker, iVertex);
                        if (solver == TNE2_EULER)
                            Buffer_Send_Mach[nVertex_Surface] = sqrt(sol->GetVariable(iPoint)->GetVelocity2()) / sol->GetVariable(iPoint)->GetSoundSpeed();
                        if (solver == TNE2_NAVIER_STOKES) {
                            Buffer_Send_SkinFriction[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker, iVertex);
                            Buffer_Send_HeatTransfer[nVertex_Surface] = FlowSolver->GetHeatFlux(iMarker, iVertex);
                        }
                        nVertex_Surface++;
                    }
                }

        /*--- Send the information to the master node ---*/

        MPI_Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (nDim == 3) MPI_Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Press, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Press, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (solver == EULER) MPI_Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if ((solver == NAVIER_STOKES) || (solver == RANS)) MPI_Gather(Buffer_Send_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (solver == TNE2_EULER) MPI_Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (solver == TNE2_NAVIER_STOKES) {
            MPI_Gather(Buffer_Send_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_HeatTransfer, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_HeatTransfer, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        }
        MPI_Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*--- The master node unpacks the data and writes the surface CSV file ---*/

        if (rank == MASTER_NODE) {

            /*--- Write file name with extension if unsteady ---*/
            char buffer[50];
            string filename = config->GetSurfFlowCoeff_FileName();
            ofstream surfFlowFile;

            /*--- Write file name with extension if unsteady ---*/
            strcpy(cstr, filename.c_str());
            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
                if (int(iZone) < 10) sprintf(buffer, "_0000%d.csv", int(iZone));
                if ((int(iZone) >= 10) && (int(iZone) < 100)) sprintf(buffer, "_000%d.csv", int(iZone));
                if ((int(iZone) >= 100) && (int(iZone) < 1000)) sprintf(buffer, "_00%d.csv", int(iZone));
                if ((int(iZone) >= 1000) && (int(iZone) < 10000)) sprintf(buffer, "_0%d.csv", int(iZone));
                if (int(iZone) >= 10000) sprintf(buffer, "_%d.csv", int(iZone));

            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
                if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))    sprintf(buffer, "_0000%d.csv", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))   sprintf(buffer, "_000%d.csv", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))  sprintf(buffer, "_00%d.csv", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.csv", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.csv", int(iExtIter));
            }
            else
                sprintf(buffer, ".csv");

            strcat(cstr, buffer);
            surfFlowFile.precision(15);
            surfFlowFile.open(cstr, ios::out);

            surfFlowFile << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
            if (nDim == 3) surfFlowFile << "\"z_coord\", ";
            surfFlowFile << "\"Pressure\", \"Pressure_Coefficient\", ";

            switch (solver) {
                case EULER: surfFlowFile << "\"Mach_Number\"" << endl; break;
                case NAVIER_STOKES: case RANS: surfFlowFile << "\"Skin_Friction_Coefficient\"" << endl; break;
                case TNE2_EULER: surfFlowFile << "\"Mach_Number\"" << endl; break;
                case TNE2_NAVIER_STOKES: surfFlowFile << "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
            }

            /*--- Loop through all of the collected data and write each node's values ---*/

            unsigned long Total_Index;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

                    /*--- Current index position and global index ---*/
                    Total_Index = iProcessor*MaxLocalVertex_Surface + iVertex;
                    Global_Index = Buffer_Recv_GlobalIndex[Total_Index];

                    /*--- Retrieve the merged data for this node ---*/
                    xCoord = Buffer_Recv_Coord_x[Total_Index];
                    yCoord = Buffer_Recv_Coord_y[Total_Index];
                    if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
                    Pressure = Buffer_Recv_Press[Total_Index];
                    PressCoeff = Buffer_Recv_CPress[Total_Index];

                    /*--- Write the first part of the data ---*/
                    surfFlowFile << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
                    if (nDim == 3) surfFlowFile << scientific << zCoord << ", ";
                    surfFlowFile << scientific << Pressure << ", " << PressCoeff << ", ";

                    /*--- Write the solver-dependent part of the data ---*/
                    switch (solver) {
                        case EULER:
                            Mach = Buffer_Recv_Mach[Total_Index];
                            surfFlowFile << scientific << Mach << endl;
                            break;
                        case NAVIER_STOKES: case RANS:
                            SkinFrictionCoeff = Buffer_Recv_SkinFriction[Total_Index];
                            surfFlowFile << scientific << SkinFrictionCoeff << endl;
                            break;
                        case TNE2_EULER:
                            Mach = Buffer_Recv_Mach[Total_Index];
                            surfFlowFile << scientific << Mach << endl;
                            break;
                        case TNE2_NAVIER_STOKES:
                            SkinFrictionCoeff = Buffer_Recv_SkinFriction[Total_Index];
                            surfFlowFile << scientific << SkinFrictionCoeff << endl;
                            HeatFlux = Buffer_Recv_HeatTransfer[Total_Index];
                            surfFlowFile << scientific << HeatFlux << endl;
                            break;
                    }
                }
            }

            /*--- Close the CSV file ---*/
            surfFlowFile.close();

            /*--- Release the recv buffers on the master node ---*/

            delete[] Buffer_Recv_Coord_x;
            delete[] Buffer_Recv_Coord_y;
            if (nDim == 3) delete[] Buffer_Recv_Coord_z;
            delete[] Buffer_Recv_Press;
            delete[] Buffer_Recv_CPress;
            delete[] Buffer_Recv_Mach;
            delete[] Buffer_Recv_SkinFriction;
            delete[] Buffer_Recv_HeatTransfer;
            delete[] Buffer_Recv_GlobalIndex;

            delete[] Buffer_Recv_nVertex;

        }

        /*--- Release the memory for the remaining buffers and exit ---*/

        delete[] Buffer_Send_Coord_x;
        delete[] Buffer_Send_Coord_y;
        delete[] Buffer_Send_Coord_z;
        delete[] Buffer_Send_Press;
        delete[] Buffer_Send_CPress;
        delete[] Buffer_Send_Mach;
        delete[] Buffer_Send_SkinFriction;
        delete[] Buffer_Send_HeatTransfer;
        delete[] Buffer_Send_GlobalIndex;
#endif
        return true;
    }

    
    /*
     * \brief Create and write the file with the adjoint coefficients on the surface.
     */
    bool Output::SetSurfaceCsvAdjoint(ISolution* sol, IMesh* mesh, IProcData* procData, unsigned long iExtIter, unsigned short val_iZone)
    {
#ifndef ARIES_HAVE_MPI
        unsigned long iPoint, iVertex;
        double *solution;
        double xCoord, yCoord, zCoord;
        unsigned short iMarker;
        char cstr[200], buffer[50];
        ofstream surfAdjFile;

        /*
         *  Write file name with extension if unsteady
         */
        strcpy(cstr, procData->GetSurfAdjCoeffFileName().c_str());

        if (procData->GetUnsteadyType() == TIME_SPECTRAL)
        {
            if ( int(val_iZone) <     10)                              sprintf(buffer, "_0000%d.csv", int(val_iZone));
            if ((int(val_iZone) >=    10) && (int(val_iZone) <   100)) sprintf(buffer, "_000%d.csv", int(val_iZone));
            if ((int(val_iZone) >=   100) && (int(val_iZone) <  1000)) sprintf(buffer, "_00%d.csv", int(val_iZone));
            if ((int(val_iZone) >=  1000) && (int(val_iZone) < 10000)) sprintf(buffer, "_0%d.csv", int(val_iZone));
            if ( int(val_iZone) >= 10000)                              sprintf(buffer, "_%d.csv", int(val_iZone));
        }
        else if (procData->GetUnsteadyType() && procData->GetWriteUnsteady())
        {
            if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))
                sprintf(buffer, "_0000%d.csv", int(iExtIter));
            if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))
                sprintf(buffer, "_000%d.csv", int(iExtIter));
            if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))
                sprintf(buffer, "_00%d.csv", int(iExtIter));
            if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000))
                sprintf(buffer, "_0%d.csv", int(iExtIter));
            if (int(iExtIter) >= 10000)
                sprintf(buffer, "_%d.csv", int(iExtIter));
        }
        else
            sprintf(buffer, ".csv");

        strcat(cstr, buffer);
        surfAdjFile.precision(15);
        surfAdjFile.open(cstr, ios::out);

        if (mesh->GetDim() == 2)
        {
            surfAdjFile << "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
            for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
            {
                if (procData->IsMarkerPlotted(iMarker))
                {
                    for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                    {
                        iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                            
                        solution = sol->GetVariableAdj(iPoint)->GetSolution();
                        xCoord      = mesh->GetNode(iPoint)->GetCoord(0);
                        yCoord      = mesh->GetNode(iPoint)->GetCoord(1);
                        /*--- If US system, the output should be in inches ---*/
                        if (procData->GetSystemMeasurements() == US)
                        {
                            xCoord *= 12.0;
                            yCoord *= 12.0;
                        }

                        surfAdjFile << scientific << iPoint << ", " << sol->GetCoeffSensitivity(iMarker, iVertex) << ", "
                                    << solution[0] << ", " << solution[1] << ", " << solution[2] << ", " << solution[3] << ", "
                                    << xCoord << ", " << yCoord << endl;
                    }
                }
            }
        }

        if (mesh->GetDim() == 3)
        {
            surfAdjFile << "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
            for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
            {
                if (procData->IsMarkerPlotted(iMarker))
                {
                    for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                    {
                        iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                            
                        solution = sol->GetVariableAdj(iPoint)->GetSolution();
                        xCoord      = mesh->GetNode(iPoint)->GetCoord(0);
                        yCoord      = mesh->GetNode(iPoint)->GetCoord(1);
                        zCoord      = mesh->GetNode(iPoint)->GetCoord(3);
                        /*--- If US system, the output should be in inches ---*/

                        if (procData->GetSystemMeasurements() == US)
                        {
                            xCoord *= 12.0;
                            yCoord *= 12.0;
                            zCoord *= 12.0;
                        }

                        surfAdjFile << scientific << iPoint << ", " << sol->GetCoeffSensitivity(iMarker, iVertex) << ", "
                                    << solution[0] << ", " << solution[1] << ", " << solution[2] << ", " << solution[3] << ", " << solution[4] << ", "
                                    << xCoord << ", " << yCoord << ", " << zCoord << endl;
                    }
                }
            }
        }
        surfAdjFile.close();
#else
        int rank, iProcessor, nProcessor;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

        unsigned short nDim = geometry->GetnDim(), iMarker;
        double *Solution, *Normal, *d, *Coord;
        unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
                MaxLocalVertex_Surface = 0, nBuffer_Scalar;
        unsigned long *Buffer_Receive_nVertex = NULL;
        ofstream surfAdjFile;

        /*--- Write the surface .csv file ---*/
        nLocalVertex_Surface = 0;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Plotting(iMarker) == YES)
                for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
                }

        if (rank == MASTER_NODE)
            Buffer_Receive_nVertex = new unsigned long[nProcessor];

        Buffer_Send_nVertex[0] = nLocalVertex_Surface;

        MPI_Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

        double *Buffer_Send_Coord_x = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_Coord_y = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_Coord_z = new double[MaxLocalVertex_Surface];
        unsigned long *Buffer_Send_GlobalPoint = new unsigned long[MaxLocalVertex_Surface];
        double *Buffer_Send_Sensitivity = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_PsiRho = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_Phi_x = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_Phi_y = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_Phi_z = new double[MaxLocalVertex_Surface];
        double *Buffer_Send_PsiE = new double[MaxLocalVertex_Surface];

        nVertex_Surface = 0;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Plotting(iMarker) == YES)
                for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    if (geometry->node[iPoint]->GetDomain()) {
                        Solution = AdjSolver->node[iPoint]->GetSolution();
                        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                        Coord = geometry->node[iPoint]->GetCoord();
                        d = AdjSolver->node[iPoint]->GetForceProj_Vector();
                        Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
                        Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
                        Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
                        Buffer_Send_Sensitivity[nVertex_Surface] = AdjSolver->GetCSensitivity(iMarker, iVertex);
                        Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
                        Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
                        Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
                        if (nDim == 2) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
                        if (nDim == 3) {
                            Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
                            Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
                            Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
                        }

                        /*--- If US system, the output should be in inches ---*/

                        if (config->GetSystemMeasurements() == US) {
                            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
                            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
                            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
                        }

                        nVertex_Surface++;
                    }
                }

        double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
                *Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
                *Buffer_Receive_PsiE = NULL;
        unsigned long *Buffer_Receive_GlobalPoint = NULL;

        if (rank == MASTER_NODE) {
            Buffer_Receive_Coord_x = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_Coord_y = new double[nProcessor*MaxLocalVertex_Surface];
            if (nDim == 3) Buffer_Receive_Coord_z = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_GlobalPoint = new unsigned long[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_Sensitivity = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_PsiRho = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_Phi_x = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_Phi_y = new double[nProcessor*MaxLocalVertex_Surface];
            if (nDim == 3) Buffer_Receive_Phi_z = new double[nProcessor*MaxLocalVertex_Surface];
            Buffer_Receive_PsiE = new double[nProcessor*MaxLocalVertex_Surface];
        }

        nBuffer_Scalar = MaxLocalVertex_Surface;

        /*--- Send the information to the Master node ---*/
        MPI_Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (nDim == 3) MPI_Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (nDim == 3) MPI_Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

        /*--- The master node is the one who writes the surface files ---*/
        if (rank == MASTER_NODE) {
            unsigned long iVertex, GlobalPoint, position;
            char cstr[200], buffer[50];
            ofstream surfAdjFile;
            string filename = config->GetSurfAdjCoeff_FileName();

            /*--- Write file name with extension if unsteady ---*/
            strcpy(cstr, filename.c_str());

            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
                if (int(val_iZone) < 10) sprintf(buffer, "_0000%d.csv", int(val_iZone));
                if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf(buffer, "_000%d.csv", int(val_iZone));
                if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf(buffer, "_00%d.csv", int(val_iZone));
                if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf(buffer, "_0%d.csv", int(val_iZone));
                if (int(val_iZone) >= 10000) sprintf(buffer, "_%d.csv", int(val_iZone));

            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
                if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf(buffer, "_0000%d.csv", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf(buffer, "_000%d.csv", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf(buffer, "_00%d.csv", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.csv", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.csv", int(iExtIter));
            }
            else
                sprintf(buffer, ".csv");

            strcat(cstr, buffer);
            surfAdjFile.open(cstr, ios::out);
            surfAdjFile.precision(15);

            /*--- Write the 2D surface flow coefficient file ---*/
            if (geometry->GetnDim() == 2) {

                surfAdjFile << "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;

                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {

                        position = iProcessor*MaxLocalVertex_Surface + iVertex;
                        GlobalPoint = Buffer_Receive_GlobalPoint[position];

                        surfAdjFile << scientific << GlobalPoint <<
                                ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
                                ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
                                ", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
                                ", " << Buffer_Receive_Coord_y[position] << endl;
                    }
            }

            /*--- Write the 3D surface flow coefficient file ---*/
            if (geometry->GetnDim() == 3) {

                surfAdjFile << "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;

                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
                        position = iProcessor*MaxLocalVertex_Surface + iVertex;
                        GlobalPoint = Buffer_Receive_GlobalPoint[position];

                        surfAdjFile << scientific << GlobalPoint <<
                                ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
                                ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
                                ", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
                                ", " << Buffer_Receive_Coord_y[position] << ", " << Buffer_Receive_Coord_z[position] << endl;
                    }
            }

        }

        if (rank == MASTER_NODE) {
            delete[] Buffer_Receive_nVertex;
            delete[] Buffer_Receive_Coord_x;
            delete[] Buffer_Receive_Coord_y;
            if (nDim == 3) delete[] Buffer_Receive_Coord_z;
            delete[] Buffer_Receive_Sensitivity;
            delete[] Buffer_Receive_PsiRho;
            delete[] Buffer_Receive_Phi_x;
            delete[] Buffer_Receive_Phi_y;
            if (nDim == 3) delete[] Buffer_Receive_Phi_z;
            delete[] Buffer_Receive_PsiE;
            delete[] Buffer_Receive_GlobalPoint;
        }

        delete[] Buffer_Send_Coord_x;
        delete[] Buffer_Send_Coord_y;
        delete[] Buffer_Send_Coord_z;
        delete[] Buffer_Send_GlobalPoint;
        delete[] Buffer_Send_Sensitivity;
        delete[] Buffer_Send_PsiRho;
        delete[] Buffer_Send_Phi_x;
        delete[] Buffer_Send_Phi_y;
        delete[] Buffer_Send_Phi_z;
        delete[] Buffer_Send_PsiE;

        surfAdjFile.close();

#endif
        return true;
    }

    bool Output::SetSurfaceCsvLinearized(ISolution* sol, IMesh* mesh, IProcData* procData, string val_filename, unsigned long iExtIter)
    {

        return true;
    }


    /*!
     * \brief Merge the geometry into a data structure used for output file writing.
     */
    void Output::MergeConnectivity(IProcData* procData, IMesh* mesh, unsigned short iZone)
    {
        int rank = MASTER_NODE;
        int size = SINGLE_NODE;

#ifdef ARIES_HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

        /*
         *  Flags identifying the types of files to be written.
         */
        bool isWrtVol  = procData->IsWrtVolSol(iZone);
        bool isWrtSurf = procData->IsWrtSurfSol(iZone);

        /*
         *  Merge connectivity for each type of element (excluding halos). Note
         *  that we only need to merge the connectivity once, as it does not change
         *  during computation. Check whether the base file has been written.
         */

        /*
         *  Merge volumetric grid.
         */
        if (isWrtVol)
        {
            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalTrias != 0))
                cout << "Merging volumetric triangle grid connectivity." << endl;
            MergeVolumetricConnectivity(procData, mesh, TRIANGLE);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalQuads != 0))
                cout << "Merging volumetric rectangle grid connectivity." << endl;
            MergeVolumetricConnectivity(procData, mesh, RECTANGLE);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalTetrs != 0))
                cout << "Merging volumetric tetrahedron grid connectivity." << endl;
            MergeVolumetricConnectivity(procData, mesh, TETRAHEDRON);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalHexas != 0))
                cout << "Merging volumetric hexahedron grid connectivity." << endl;
            MergeVolumetricConnectivity(procData, mesh, HEXAHEDRON);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalPriss != 0))
                cout << "Merging volumetric prism grid connectivity." << endl;
            MergeVolumetricConnectivity(procData, mesh, PRISM);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalPyras != 0))
                cout << "Merging volumetric pyramid grid connectivity." << endl;
            MergeVolumetricConnectivity(procData, mesh, PYRAMID);
        }

        /*
         *  Merge surface grid.
         */
        if (isWrtSurf)
        {
            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalBoundLines != 0))
                cout << "Merging surface line grid connectivity." << endl;
            MergeSurfaceConnectivity(procData, mesh, LINE);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalBoundTrias != 0))
                cout << "Merging surface triangle grid connectivity." << endl;
            MergeSurfaceConnectivity(procData, mesh, TRIANGLE);

            if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (d_numGlobalBoundQuads != 0))
                cout << "Merging surface rectangle grid connectivity." << endl;
            MergeSurfaceConnectivity(procData, mesh, RECTANGLE);
        }

        /*
         *  Update total number of volume elements after merge.
         */
        d_numGlobalElems = d_numGlobalTrias + d_numGlobalQuads + d_numGlobalTetrs + d_numGlobalHexas + d_numGlobalPyras + d_numGlobalPriss;

        /*
         *  Update total number of surface elements after merge.
         */
        d_numSurfElems = d_numGlobalBoundLines + d_numGlobalBoundTrias + d_numGlobalBoundQuads;
    }


    /*!
     * \brief Merge the connectivity for a single element type from all processors.
     */
    void Output::MergeVolumetricConnectivity(IProcData* procData, IMesh* mesh, unsigned short elemType) 
    {
        int iProcessor;
        unsigned short NODES_PER_ELEMENT;
        unsigned long iPoint, iNode, jNode;
        unsigned long iElem = 0;
        unsigned long nLocalElem = 0, nElemTotal = 0;

        unsigned long iVertex, iMarker;
        unsigned long jElem;
        int sendRecv, recvFrom;

        unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
        unsigned long nBuffer_Scalar = 0;
        unsigned long kNode = 0, kElem = 0;
        unsigned long maxLocalElem = 0, iGlobalIndex, jPoint, kPoint;

        bool isWrtHalo = procData->IsWrtHalo();
        bool *isWriteElem = NULL, notPeriodic, notHalo, isAddedPeriodic;

        int *connElem = NULL;

        int rank = MASTER_NODE;
        int size = SINGLE_NODE;

#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
        /*
         *  Store the local number of this element type and the number of nodes
         *  per this element type. In serial, this will be the total number of this
         *  element type in the entire mesh. In parallel, it is the number on only
         *  the current partition.
         */
        switch (elemType) 
        {
            case TRIANGLE:
                nLocalElem = mesh->GetNumElemTria();
                NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
                break;
            case RECTANGLE:
                nLocalElem = mesh->GetNumElemQuad();
                NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
                break;
            case TETRAHEDRON:
                nLocalElem = mesh->GetNumElemTetr();
                NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
                break;
            case HEXAHEDRON:
                nLocalElem = mesh->GetNumElemHexa();
                NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
                break;
            case PRISM:
                nLocalElem = mesh->GetNumElemPris();
                NODES_PER_ELEMENT = N_POINTS_PRISM;
                break;
            case PYRAMID:
                nLocalElem = mesh->GetNumElemPyra();
                NODES_PER_ELEMENT = N_POINTS_PYRAMID;
                break;
            default:
                cout << "Error: Unrecognized element type \n";
                exit(EXIT_FAILURE);
                break;
        }

        /*
         *  Find the max number of this element type among all
         *  partitions and set up buffers.
         */
        Buffer_Send_nElem[0] = nLocalElem;
        if (rank == MASTER_NODE) 
            Buffer_Recv_nElem = new unsigned long[size];

#ifdef ARIES_HAVE_MPI
        MPI_Allreduce(&nLocalElem, &maxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
        maxLocalElem = nLocalElem;
        Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

        nBuffer_Scalar = maxLocalElem*NODES_PER_ELEMENT;

        /*
         *  Send and Recv buffers
         */
        unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
        unsigned long *Buffer_Recv_Elem = NULL;

        unsigned short *Buffer_Send_Halo = new unsigned short[maxLocalElem];
        unsigned short *Buffer_Recv_Halo = NULL;

        /*
         *  Prepare the receive buffers on the master node only.
         */
        if (rank == MASTER_NODE) 
        {
            Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
            Buffer_Recv_Halo = new unsigned short[size*maxLocalElem];
            connElem = new int[size*maxLocalElem*NODES_PER_ELEMENT];
        }

        /*
         *  Force the removal of all added periodic elements (use global index).
         *  First, we isolate and create a list of all added periodic points, excluding
         *  those that we part of the original domain (we want these to be in the
         *  output files).
         */
        vector<unsigned long> addedPeriodic;
        addedPeriodic.clear();
        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->GetMarkerBCType(iMarker) == SEND_RECEIVE)
            {
                sendRecv = procData->GetMarkerSendRecv(iMarker);
                for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                {
                    iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                    if ((mesh->GetVertex(iMarker, iVertex)->GetRotationType() > 0) &&
                        (mesh->GetVertex(iMarker, iVertex)->GetRotationType() % 2 == 0) &&
                        (sendRecv < 0))
                    {
                        addedPeriodic.push_back(mesh->GetNode(iPoint)->GetGlobalIndex());
                    }
                }
            }
        }

        /*
         *  Now we communicate this information to all processors, so that they
         *  can force the removal of these particular nodes by flagging them as halo
         *  points. In general, this should be a small percentage of the total mesh,
         *  so the communication/storage costs here shouldn't be prohibitive.
         */
            
        /*
         *  First communicate the number of points that each rank has found
         */
        unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
        unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
        Buffer_Recv_nAddedPeriodic = new unsigned long[size];

        nAddedPeriodic = addedPeriodic.size();
        Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef ARIES_HAVE_MPI
        MPI_Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nAddedPeriodic, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
        maxAddedPeriodic = nAddedPeriodic;
        Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

        /*
         *  Communicate the global index values of all added periodic nodes.
         */
        unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
        unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

        for (iPoint = 0; iPoint < addedPeriodic.size(); iPoint++) 
        {
            Buffer_Send_AddedPeriodic[iPoint] = addedPeriodic[iPoint];
        }

        /*
         *  Gather the element connectivity information. All processors will now
         *  have a copy of the global index values for all added periodic points.
         */

#ifdef ARIES_HAVE_MPI
        MPI_Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG, Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
            Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

        /*
         *  Search all send/recv boundaries on this partition for halo cells. In
         *  particular, consider only the recv conditions (these are the true halo
         *  nodes). Check the ranks of the processors that are communicating and
         *  choose to keep only the halo cells from the higher rank processor. Here,
         *  we are also choosing to keep periodic nodes that were part of the original
         *  domain. We will check the communicated list of added periodic points.
         */
        int *localHalo = new int[mesh->GetNumPoints()];
        for (iPoint = 0; iPoint < mesh->GetNumPoints(); iPoint++)
            localHalo[iPoint] = !mesh->GetNode(iPoint)->GetDomain();

        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->GetMarkerBCType(iMarker) == SEND_RECEIVE)
            {
                sendRecv = procData->GetMarkerSendRecv(iMarker);
                recvFrom = abs(sendRecv) - 1;

                for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++) 
                {
                    iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                    iGlobalIndex = mesh->GetNode(iPoint)->GetGlobalIndex();

                    /*
                     *  We need to keep one copy of overlapping halo cells.
                     */
                    notHalo = ((mesh->GetVertex(iMarker, iVertex)->GetRotationType() == 0) &&
                               (sendRecv < 0) && (rank > recvFrom));

                    /*
                     *  We want to keep the periodic nodes that were part of the original domain
                     */
                    notPeriodic = ((mesh->GetVertex(iMarker, iVertex)->GetRotationType()  > 0) &&
                                   (mesh->GetVertex(iMarker, iVertex)->GetRotationType()  % 2 == 1) &&
                                   (sendRecv < 0));

                    /*
                     *  Lastly, check that this isn't an added periodic point that
                     *  we will forcibly remove. Use the communicated list of these points.
                     */
                    isAddedPeriodic = false; kPoint = 0;
                    for (iProcessor = 0; iProcessor < size; iProcessor++) 
                    {
                        for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) 
                        {
                            if (iGlobalIndex == Buffer_Recv_AddedPeriodic[kPoint + jPoint])
                                isAddedPeriodic = true;
                        }
                        /*
                         *  Adjust jNode to index of next proc's data in the buffers.
                         */
                        kPoint = (iProcessor + 1)*maxAddedPeriodic;
                    }

                    /*
                     *  If we found either of these types of nodes, flag them to be kept.
                     */
                    if ((notHalo || notPeriodic) && !isAddedPeriodic) 
                    {
                        localHalo[iPoint] = false;
                    }
                }
            }
        }

        /*
         *  Loop over all elements in this partition and load the
         *  elements of the current type into the buffer to be sent to
         *  the master node.
         */
        jNode = 0; jElem = 0;
        for (iElem = 0; iElem < mesh->GetNumElements(); iElem++) 
        {
            if (mesh->GetElement(iElem)->GetVtkType() == elemType) 
            {
                /*
                 *  Loop over all nodes in this element and load the
                 *  connectivity into the send buffer.
                 */
                Buffer_Send_Halo[jElem] = false;
                for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) 
                {
                    /*
                     *  Store the global index values directly.
                     */
                    iPoint = mesh->GetElement(iElem)->GetNode(iNode);
                    Buffer_Send_Elem[jNode] = mesh->GetNode(iPoint)->GetGlobalIndex();

                    /*
                     *  Check if this is a halo node. If so, flag this element
                     *  as a halo cell. We will use this later to sort and remove
                     *  any duplicates from the connectivity list.
                     */
                    if (localHalo[iPoint]) 
                    {
                        Buffer_Send_Halo[jElem] = true;
                    }

                    /*
                     *  Increment jNode as the counter. We need this because iElem
                     *  may include other elements that we skip over during this loop.
                     */
                    jNode++;
                }
                jElem++;
            }
        }

        /*
         *  Gather the element connectivity information.
         */
#ifdef ARIES_HAVE_MPI
        MPI_Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Halo, maxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, maxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) 
            Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
        for (iPoint = 0; iPoint < maxLocalElem; iPoint++) 
            Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

        /*
         *  The master node unpacks and sorts the connectivity.
         */
        if (rank == MASTER_NODE) 
        {
            /*
             *  We need to remove any duplicate elements (halo cells) that
             *  exist on multiple partitions. Start by initializing all elements
             *  to the "write" state by using a boolean array.
             */
            isWriteElem = new bool[size*maxLocalElem];
            for (iElem = 0; iElem < size*maxLocalElem; iElem++)
            {
                isWriteElem[iElem] = true;
            }

            /*
             *  Remove the rind layer from the solution only if requested
             */
            if (!isWrtHalo) 
            {
                /*
                 *  Loop for flagging duplicate elements so that they are not
                 *  included in the final connectivity list.
                 */
                kElem = 0;
                for (iProcessor = 0; iProcessor < size; iProcessor++) 
                {
                    for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) 
                    {
                        /*
                         *  Check if this element was marked as a halo.
                         */
                        if (Buffer_Recv_Halo[kElem + iElem])
                            isWriteElem[kElem + iElem] = false;
                    }
                    kElem = (iProcessor + 1)*maxLocalElem;
                }
            }

            /*
             *  Store the unique connectivity list for this element type.
             */
            jNode = 0; kNode = 0; jElem = 0; nElemTotal = 0;
            for (iProcessor = 0; iProcessor < size; iProcessor++) 
            {
                for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) 
                {
                    /*
                     *  Only write the elements that were flagged for it.
                     */
                    if (isWriteElem[jElem + iElem]) 
                    {
                        /*
                         *  Increment total count for this element type
                         */
                        nElemTotal++;

                        /*
                         *  Get global index, then loop over each variable and store.
                         *  Note that we are adding one to the index value because CGNS/Tecplot
                         *  use 1-based indexing.
                         */
                        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++)
                        {
                            connElem[kNode] = (int)Buffer_Recv_Elem[jNode + iElem*NODES_PER_ELEMENT + iNode] + 1;
                            kNode++;
                        }
                    }
                }
                /*
                 *  Adjust jNode to index of next proc's data in the buffers.
                 */
                jElem = (iProcessor + 1)*maxLocalElem;
                jNode = (iProcessor + 1)*nBuffer_Scalar;
            }
        }

        /*
         *  Immediately release the temporary buffers.
         */
        delete[] Buffer_Send_Elem;
        delete[] Buffer_Send_Halo;
        delete[] Buffer_Recv_nAddedPeriodic;
        delete[] Buffer_Send_AddedPeriodic;
        delete[] Buffer_Recv_AddedPeriodic;
        delete[] localHalo;
        if (rank == MASTER_NODE) 
        {
            delete[] Buffer_Recv_nElem;
            delete[] Buffer_Recv_Elem;
            delete[] Buffer_Recv_Halo;
            delete[] isWriteElem;
        }

        /*
         *  Store the particular global element count in the class data,
         *  and set the class data pointer to the connectivity array.
         */
        if (rank == MASTER_NODE)
        {
            switch (elemType)
            {
                case TRIANGLE:
                    d_numGlobalTrias = nElemTotal;
                    if (d_numGlobalTrias > 0) d_connTria = connElem;
                    break;
                case RECTANGLE:
                    d_numGlobalQuads = nElemTotal;
                    if (d_numGlobalQuads > 0) d_connQuad = connElem;
                    break;
                case TETRAHEDRON:
                    d_numGlobalTetrs = nElemTotal;
                    if (d_numGlobalTetrs > 0) d_connTetr = connElem;
                    break;
                case HEXAHEDRON:
                    d_numGlobalHexas = nElemTotal;
                    if (d_numGlobalHexas > 0) d_connHexa = connElem;
                    break;
                case PRISM:
                    d_numGlobalPriss = nElemTotal;
                    if (d_numGlobalPriss > 0) d_connPris = connElem;
                    break;
                case PYRAMID:
                    d_numGlobalPyras = nElemTotal;
                    if (d_numGlobalPyras > 0) d_connPyra = connElem;
                    break;
                default:
                    cout << "Error: Unrecognized element type \n";
                    exit(EXIT_FAILURE); break;
            }
        }
    }


    /*!
     * \brief Merge the connectivity for a single element type from all processors.
     */
    void Output::MergeSurfaceConnectivity(IProcData* procData, IMesh* mesh, unsigned short elemType)
    {
        unsigned short NODES_PER_ELEMENT;

        unsigned short iMarker;
        unsigned long iPoint, iNode, jNode;
        unsigned long iElem = 0;
        unsigned long nLocalElem = 0, nElemTotal = 0;

        int iProcessor;
        unsigned long jElem;
        unsigned long iVertex;

        int sendRecv, recvFrom;

        unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
        unsigned long nBufferScalar = 0;
        unsigned long kNode = 0, kElem = 0;
        unsigned long maxLocalElem = 0, iGlobalIndex, jPoint, kPoint;

        bool isWrtHalo = procData->IsWrtHalo();
        bool *writeElem = NULL, notPeriodic, notHalo, isAddedPeriodic;


        int *connElem = NULL;

        int rank = MASTER_NODE;
        int size = SINGLE_NODE;

#ifdef ARIES_HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

        /*
         *  Store the local number of this element type and the number of nodes
         *  per this element type. In serial, this will be the total number of this
         *  element type in the entire mesh. In parallel, it is the number on only
         *  the current partition.
         */
        nLocalElem = 0;

        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->IsMarkerPlotted(iMarker))
            {
                for (iElem = 0; iElem < mesh->GetNumElemBoundary(iMarker); iElem++)
                {
                    if (mesh->GetElementBoundary(iMarker, iElem)->GetVtkType() == elemType)
                    {
                        nLocalElem++;
                    }
                }
            }
        }

        switch (elemType)
        {
            case LINE:
                NODES_PER_ELEMENT = N_POINTS_LINE;
                break;
            case TRIANGLE:
                NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
                break;
            case RECTANGLE:
                NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
                break;
            default:
                cout << "Error: Unrecognized element type \n";
                exit(EXIT_FAILURE);
                break;
        }

        /*
         *  Find the max number of this element type among all
         *  partitions and set up buffers.
         */
        Buffer_Send_nElem[0] = nLocalElem;
        if (rank == MASTER_NODE)
            Buffer_Recv_nElem = new unsigned long[size];

#ifdef ARIES_HAVE_MPI
        MPI_Allreduce(&nLocalElem, &maxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
        maxLocalElem = nLocalElem;
        Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

        nBufferScalar = maxLocalElem*NODES_PER_ELEMENT;

        /*
         *  Send and Recv buffers
         */
        unsigned long *Buffer_Send_Elem = new unsigned long[nBufferScalar];
        unsigned long *Buffer_Recv_Elem = NULL;

        unsigned short *Buffer_Send_Halo = new unsigned short[maxLocalElem];
        unsigned short *Buffer_Recv_Halo = NULL;

        /*
         *  Prepare the receive buffers on the master node only.
         */
        if (rank == MASTER_NODE)
        {
            Buffer_Recv_Elem = new unsigned long[size*nBufferScalar];
            Buffer_Recv_Halo = new unsigned short[size*maxLocalElem];
            connElem = new int[size*maxLocalElem*NODES_PER_ELEMENT];
        }

        /*
         *  Force the removal of all added periodic elements (use global index).
         *  First, we isolate and create a list of all added periodic points, excluding
         *  those that we part of the original domain (we want these to be in the
         *  output files).
         */
        vector<unsigned long> addedPeriodic;
        addedPeriodic.clear();
        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->GetMarkerBCType(iMarker) == SEND_RECEIVE)
            {
                sendRecv = procData->GetMarkerSendRecv(iMarker);
                for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                {
                    iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                    if (( mesh->GetVertex(iMarker, iVertex)->GetRotationType() > 0) &&
                        ( mesh->GetVertex(iMarker, iVertex)->GetRotationType() % 2 == 0) &&
                        ( sendRecv < 0))
                    {
                        addedPeriodic.push_back(mesh->GetNode(iPoint)->GetGlobalIndex());
                    }
                }
            }
        }

        /*
         *  Now we communicate this information to all processors, so that they
         *  can force the removal of these particular nodes by flagging them as halo
         *  points. In general, this should be a small percentage of the total mesh,
         *  so the communication/storage costs here shouldn't be prohibitive.
         */

        /*
         *  First communicate the number of points that each rank has found
         */
        unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
        unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
        Buffer_Recv_nAddedPeriodic = new unsigned long[size];

        nAddedPeriodic = addedPeriodic.size();
        Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef ARIES_HAVE_MPI
        MPI_Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nAddedPeriodic, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
        maxAddedPeriodic = nAddedPeriodic;
        Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

        /*
         *  Communicate the global index values of all added periodic nodes.
         */
        unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
        unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

        for (iPoint = 0; iPoint < addedPeriodic.size(); iPoint++)
        {
            Buffer_Send_AddedPeriodic[iPoint] = addedPeriodic[iPoint];
        }

        /*
         *  Gather the element connectivity information. All processors will now
         *  have a copy of the global index values for all added periodic points.
         */
#ifdef ARIES_HAVE_MPI
        MPI_Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                      Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                      MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
            Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

        /*
         *  Search all send/recv boundaries on this partition for halo cells. In
         *  particular, consider only the recv conditions (these are the true halo
         *  nodes). Check the ranks of the processors that are communicating and
         *  choose to keep only the halo cells from the higher rank processor. Here,
         *  we are also choosing to keep periodic nodes that were part of the original
         *  domain. We will check the communicated list of added periodic points.
         */
        int *localHalo = new int[mesh->GetNumPoints()];
        for (iPoint = 0; iPoint < mesh->GetNumPoints(); iPoint++)
            localHalo[iPoint] = !mesh->GetNode(iPoint)->GetDomain();

        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->GetMarkerBCType(iMarker) == SEND_RECEIVE)
            {
                sendRecv = procData->GetMarkerSendRecv(iMarker);
                recvFrom = abs(sendRecv) - 1;

                for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                {
                    iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                    iGlobalIndex = mesh->GetNode(iPoint)->GetGlobalIndex();

                    /*
                     *  We need to keep one copy of overlapping halo cells.
                     */
                    notHalo = ((mesh->GetVertex(iMarker, iVertex)->GetRotationType() == 0) &&
                               (sendRecv < 0) &&
                               (rank > recvFrom));

                    /*
                     *  We want to keep the periodic nodes that were part of the original domain
                     */
                    notPeriodic = ((mesh->GetVertex(iMarker, iVertex)->GetRotationType() > 0) &&
                                   (mesh->GetVertex(iMarker, iVertex)->GetRotationType() % 2 == 1) &&
                                   (sendRecv < 0));

                    /*
                     *  Lastly, check that this isn't an added periodic point that
                     *  we will forcibly remove. Use the communicated list of these points.
                     */
                    isAddedPeriodic = false; kPoint = 0;
                    for (iProcessor = 0; iProcessor < size; iProcessor++)
                    {
                        for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++)
                        {
                            if (iGlobalIndex == Buffer_Recv_AddedPeriodic[kPoint + jPoint])
                                isAddedPeriodic = true;
                        }
                        /*
                         *  Adjust jNode to index of next proc's data in the buffers.
                         */
                        kPoint = (iProcessor + 1)*maxAddedPeriodic;
                    }

                    /*
                     *  If we found either of these types of nodes, flag them to be kept.
                     */
                    if ((notHalo || notPeriodic) && !isAddedPeriodic)
                    {
                        localHalo[iPoint] = false;
                    }
                }
            }
        }

        /*
         *  Loop over all elements in this partition and load the
         *  elements of the current type into the buffer to be sent to
         *  the master node.
         */
        jNode = 0; jElem = 0;
        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->IsMarkerPlotted(iMarker))
            {
                for (iElem = 0; iElem < mesh->GetNumElemBoundary(iMarker); iElem++)
                {
                    if (mesh->GetElementBoundary(iMarker, iElem)->GetVtkType() == elemType)
                    {
                        /*
                         *  Loop over all nodes in this element and load the
                         *  connectivity into the send buffer.
                         */
                        Buffer_Send_Halo[jElem] = false;
                        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++)
                        {
                            /*
                             *  Store the global index values directly.
                             */
                            iPoint = mesh->GetElementBoundary(iMarker, iElem)->GetNode(iNode);
                            Buffer_Send_Elem[jNode] = mesh->GetNode(iPoint)->GetGlobalIndex();

                            /*
                             *  Check if this is a halo node. If so, flag this element
                             *  as a halo cell. We will use this later to sort and remove
                             *  any duplicates from the connectivity list.
                             */
                            if (localHalo[iPoint])
                                Buffer_Send_Halo[jElem] = true;

                            /*
                             *  Increment jNode as the counter. We need this because iElem
                             *  may include other elements that we skip over during this loop.
                             */
                            jNode++;
                        }
                        jElem++;
                    }
                }
            }
        }

        /*
         *  Gather the element connectivity information.
         */
#ifdef HAVE_MPI
        MPI_Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Halo, maxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, maxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBufferScalar; iPoint++)
            Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];

        for (iPoint = 0; iPoint < maxLocalElem; iPoint++)
            Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

        /*
         *  The master node unpacks and sorts the connectivity.
         */
        if (rank == MASTER_NODE)
        {
            /*
             *  We need to remove any duplicate elements (halo cells) that
             *  exist on multiple partitions. Start by initializing all elements
             *  to the "write" state by using a boolean array.
             */
            writeElem = new bool[size*maxLocalElem];
            for (iElem = 0; iElem < size*maxLocalElem; iElem++)
            {
                writeElem[iElem] = true;
            }

            /*
             *  Remove the rind layer from the solution only if requested
             */
            if (!isWrtHalo)
            {
                /*
                 *  Loop for flagging duplicate elements so that they are not
                 *  included in the final connectivity list.
                 */
                kElem = 0;
                for (iProcessor = 0; iProcessor < size; iProcessor++)
                {
                    for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++)
                    {
                        /*
                         *  Check if this element was marked as a halo.
                         */
                        if (Buffer_Recv_Halo[kElem + iElem])
                            writeElem[kElem + iElem] = false;
                    }
                    kElem = (iProcessor + 1)*maxLocalElem;
                }
            }

            /*
             *  Store the unique connectivity list for this element type.
             */
            jNode = 0; kNode = 0; jElem = 0; nElemTotal = 0;
            for (iProcessor = 0; iProcessor < size; iProcessor++)
            {
                for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++)
                {
                    /*
                     *  Only write the elements that were flagged for it.
                     */
                    if (writeElem[jElem + iElem])
                    {
                        /*
                         *  Increment total count for this element type
                         */
                        nElemTotal++;

                        /*
                         *  Get global index, then loop over each variable and store.
                         *  Note that we are adding one to the index value because CGNS/Tecplot
                         *  use 1-based indexing.
                         */
                        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++)
                        {
                            connElem[kNode] = (int)Buffer_Recv_Elem[jNode + iElem*NODES_PER_ELEMENT + iNode] + 1;
                            kNode++;
                        }
                    }
                }
                /*
                 *  Adjust jNode to index of next proc's data in the buffers.
                 */
                jElem = (iProcessor + 1)*maxLocalElem;
                jNode = (iProcessor + 1)*nBufferScalar;
            }
        }

        /*
         *  mmediately release the temporary buffers.
         */
        delete[] Buffer_Send_Elem;
        delete[] Buffer_Send_Halo;
        delete[] Buffer_Recv_nAddedPeriodic;
        delete[] Buffer_Send_AddedPeriodic;
        delete[] Buffer_Recv_AddedPeriodic;
        delete[] localHalo;
        if (rank == MASTER_NODE)
        {
            delete[] Buffer_Recv_nElem;
            delete[] Buffer_Recv_Elem;
            delete[] Buffer_Recv_Halo;
            delete[] writeElem;
        }

        /*
         *  Store the particular global element count in the class data,
         *  and set the class data pointer to the connectivity array.
         */
        if (rank == MASTER_NODE)
        {
            switch (elemType)
            {
                case LINE:
                    d_numGlobalBoundLines = nElemTotal;
                    if (d_numGlobalBoundLines > 0)
                        d_connBoundLine = connElem;
                    break;
                case TRIANGLE:
                    d_numGlobalBoundTrias = nElemTotal;
                    if (d_numGlobalBoundTrias > 0)
                        d_connBoundTria = connElem;
                    break;
                case RECTANGLE:
                    d_numGlobalBoundQuads = nElemTotal;
                    if (d_numGlobalBoundQuads > 0)
                        d_connBoundQuad = connElem;
                    break;
                default:
                    cout << "Error: Unrecognized element type \n";
                    exit(EXIT_FAILURE);
                    break;
            }
        }
    }

    /*!
     * \brief Merge the node coordinates from all processors.
     */
    void Output::MergeCoordinates(IProcData* procData, IMesh* mesh)
    {
        /*
         *  Local variables needed on all processors
         */
        unsigned short iDim, nDim = mesh->GetDim();
        unsigned long iPoint;

#ifndef ARIES_HAVE_MPI
        /*
         *  In serial, the single process has access to all geometry, so simply
         *  load the coordinates into the data structure.
         */
        unsigned short iMarker;
        unsigned long iVertex, nTotalPoints = 0;
        int sendRecv;
        
        /*
         *  First, create a structure to locate any periodic halo nodes
         */
        int *localHalo = new int[mesh->GetNumPoints()];
        for (iPoint = 0; iPoint < mesh->GetNumPoints(); iPoint++)
            localHalo[iPoint] = !mesh->GetNode(iPoint)->GetDomain();

        for (iMarker = 0; iMarker < procData->GetNumMarkers(); iMarker++)
        {
            if (procData->GetMarkerBCType(iMarker) == SEND_RECEIVE)
            {
                sendRecv = procData->GetMarkerSendRecv(iMarker);
                for (iVertex = 0; iVertex < mesh->GetNumVerticesOnMarker(iMarker); iVertex++)
                {
                    iPoint = mesh->GetVertex(iMarker, iVertex)->GetNode();
                    if ((mesh->GetVertex(iMarker, iVertex)->GetRotationType() > 0) &&
                        (mesh->GetVertex(iMarker, iVertex)->GetRotationType() % 2 == 1) &&
                        (sendRecv < 0))
                    {
                        localHalo[iPoint] = false;
                    }
                }
            }
        }
        
        /*
         *  Total number of points in the mesh (this might include periodic points).
         */
        for (iPoint = 0; iPoint < mesh->GetNumPoints(); iPoint++)
            if (!localHalo[iPoint]) 
                nTotalPoints++;

        d_numGlobalPoints = nTotalPoints;
        d_numGlobalPointDomain = mesh->GenNumPointDomain();

        /*
         *  Allocate the coordinates data structure.
         */
        d_coords = new double* [nDim];
        for (iDim = 0; iDim < nDim; iDim++)
        {
            d_coords[iDim] = new double[d_numGlobalPoints];
        }

        /*
         *  Loop over the mesh to collect the coords of the local points
         */
        for (iPoint = 0; iPoint < mesh->GetNumPoints(); iPoint++)
        {
            /*
             *  Check if the node belongs to the domain (i.e, not a halo node).
             *  Sort by the global index, even in serial there is a renumbering (e.g. RCM).
             */
            if (!localHalo[iPoint])
            {
                /*
                 *  Retrieve the current coordinates at this node.
                 */
                unsigned long iGlobalIndex = mesh->GetNode(iPoint)->GetGlobalIndex();

                for (iDim = 0; iDim < nDim; iDim++)
                {
                    d_coords[iDim][iGlobalIndex] = mesh->GetNode(iPoint)->GetCoord(iDim);
                    /*
                     *  If US system, the output should be in inches
                     */
                    if ( procData->GetSystemMeasurements() == US )
                    {
                        d_coords[iDim][iGlobalIndex] *= 12.0;
                    }
                }
            }
        }

        delete[] localHalo;
#else
        /*--- MPI preprocessing ---*/
        int iProcessor, nProcessor, rank;
        unsigned long jPoint;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

        bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

        /*--- Local variables needed for merging the geometry with MPI. ---*/

        unsigned long iVertex, iMarker;

        int SendRecv, RecvFrom;

        unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
        unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
        unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0, periodicNodes = 0;

        if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

        int *Local_Halo = new int[geometry->GetnPoint()];
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
            Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

        /*--- Search all send/recv boundaries on this partition for any periodic
          nodes that were part of the original domain. We want to recover these
          for visualization purposes. ---*/

        if (Wrt_Halo) {
            nLocalPoint = geometry->GetnPoint();
        }
        else {
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
                    SendRecv = config->GetMarker_All_SendRecv(iMarker);
                    RecvFrom = abs(SendRecv) - 1;

                    /*--- Checking for less than or equal to the rank, because there may
                      be some periodic halo nodes that send info to the same rank. ---*/

                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
                        if (isPeriodic) Local_Halo[iPoint] = false;
                    }
                }
            }

            /*--- Sum total number of nodes that belong to the domain ---*/

            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                if (Local_Halo[iPoint] == false)
                    nLocalPoint++;
        }
        Buffer_Send_nPoin[0] = nLocalPoint;

        /*--- Communicate the total number of nodes on this domain. ---*/

        MPI_Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                   Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

        if (rank == MASTER_NODE) {
            d_numGlobalDomains = 0;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                d_numGlobalDomains += Buffer_Recv_nPoin[iProcessor];
            }
        }
        nBuffer_Scalar = MaxLocalPoint;

        /*--- Send and Recv buffers. ---*/

        double *Buffer_Send_X = new double[MaxLocalPoint];
        double *Buffer_Recv_X = NULL;

        double *Buffer_Send_Y = new double[MaxLocalPoint];
        double *Buffer_Recv_Y = NULL;

        double *Buffer_Send_Z, *Buffer_Recv_Z = NULL;
        if (nDim == 3) Buffer_Send_Z = new double[MaxLocalPoint];

        unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
        unsigned long *Buffer_Recv_GlobalIndex = NULL;

        /*--- Prepare the receive buffers in the master node only. ---*/

        if (rank == MASTER_NODE) {

            Buffer_Recv_X = new double[nProcessor*MaxLocalPoint];
            Buffer_Recv_Y = new double[nProcessor*MaxLocalPoint];
            if (nDim == 3) Buffer_Recv_Z = new double[nProcessor*MaxLocalPoint];
            Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

            /*--- Sum total number of nodes to be written and allocate arrays ---*/
            d_numGlobalPoints = 0;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                d_numGlobalPoints += Buffer_Recv_nPoin[iProcessor];
            }
            d_coords = new double*[nDim];
            for (iDim = 0; iDim < nDim; iDim++) {
                d_coords[iDim] = new double[d_numGlobalPoints];
            }
        }

        /*--- Main communication routine. Loop over each coordinate and perform
          the MPI comm. Temporary 1-D buffers are used to send the coordinates at
          all nodes on each partition to the master node. These are then unpacked
          by the master and sorted by global index in one large n-dim. array. ---*/

        /*--- Loop over this partition to collect the coords of the local points. ---*/
        double *d_coords_Local; jPoint = 0;
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

            /*--- Check for halos and write only if requested ---*/
            if (!Local_Halo[iPoint] || Wrt_Halo) {

                /*--- Retrieve local coordinates at this node. ---*/
                d_coords_Local = geometry->node[iPoint]->GetCoord();

                /*--- Load local coords into the temporary send buffer. ---*/
                Buffer_Send_X[jPoint] = d_coords_Local[0];
                Buffer_Send_Y[jPoint] = d_coords_Local[1];
                if (nDim == 3) Buffer_Send_Z[jPoint] = d_coords_Local[2];

                /*--- If US system, the output should be in inches ---*/

                if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
                    Buffer_Send_X[jPoint] *= 12.0;
                    Buffer_Send_Y[jPoint] *= 12.0;
                    if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
                }

                /*--- Store the global index for this local node. ---*/
                Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();

                /*--- Increment jPoint as the counter. We need this because iPoint
                  may include halo nodes that we skip over during this loop. ---*/
                jPoint++;
            }
        }
        
        /*
         *  Gather the coordinate data on the master node using MPI.
         */
        MPI_Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        if (nDim == 3)
            MPI_Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        MPI_Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

        /*
         *  The master node unpacks and sorts this variable by global index
         */
        if (rank == MASTER_NODE)
        {
            jPoint = 0;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
            {
                for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++)
                {
                    /*
                     *  Get global index, then loop over each variable and store
                     */
                    iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                    d_coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
                    d_coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
                    if (nDim == 3)
                        d_coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
                    jPoint++;
                }
                /*
                 *  Adjust jPoint to index of next proc's data in the buffers.
                 */
                jPoint = (iProcessor + 1)*nBuffer_Scalar;
            }
        }

        /*
         *  Immediately release the temporary data buffers.
         */
        delete[] Local_Halo;
        delete[] Buffer_Send_X;
        delete[] Buffer_Send_Y;
        if (nDim == 3)
            delete[] Buffer_Send_Z;
        delete[] Buffer_Send_GlobalIndex;
        if (rank == MASTER_NODE)
        {
            delete[] Buffer_Recv_X;
            delete[] Buffer_Recv_Y;
            if (nDim == 3)
                delete[] Buffer_Recv_Z;
            delete[] Buffer_Recv_GlobalIndex;
            delete[] Buffer_Recv_nPoin;
        }
#endif

    }

    
    /*!
     * \brief Merge solution into a data structure used for output file writing.
     */
    void Output::MergeSolution(IProcData* procData, IMesh* mesh, ISolution* sol, unsigned short val_iZone)
    {
        unsigned short Kind_Solver = config->GetKind_Solver();
        unsigned short iVar = 0, jVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
        unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0;
        unsigned short iVar_GridVel = 0, iVar_PressCp = 0, iVar_Density = 0, iVar_Lam = 0, iVar_MachMean = 0,
                iVar_Tempv = 0, iVar_EF = 0, iVar_Temp = 0, iVar_Mach = 0, iVar_Press = 0, iVar_TempLam = 0,
                iVar_ViscCoeffs = 0, iVar_Sens = 0, iVar_FEA = 0, iVar_Extra = 0, iVar_Eddy = 0, iVar_Sharp = 0;

        unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
        double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;

        double *Aux_Frict = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL, *Aux_Sens = NULL;

        unsigned short CurrentIndex;
        int SendRecv, RecvFrom, *Local_Halo;
        unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
        unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
        unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
        bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

        int iProcessor;
        int rank = MASTER_NODE;
        int size = SINGLE_NODE;

#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

        bool grid_movement = (config->GetGrid_Movement());
        bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
        bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
        bool freesurface = (config->GetKind_Regime() == FREESURFACE);
        bool transition = (config->GetKind_Trans_Model() == LM);
        bool flow = ((config->GetKind_Solver() == EULER) ||
                     (config->GetKind_Solver() == NAVIER_STOKES) ||
                     (config->GetKind_Solver() == RANS) ||
                     (config->GetKind_Solver() == ADJ_EULER) ||
                     (config->GetKind_Solver() == ADJ_NAVIER_STOKES) ||
                     (config->GetKind_Solver() == ADJ_RANS));

        unsigned short iDim;
        unsigned short nDim = geometry->GetnDim();
        double RefAreaCoeff = config->GetRefAreaCoeff();
        double Gamma = config->GetGamma();
        double RefVel2, *Normal, Area;

        /*--- Set the non-dimensionalization ---*/
        if (flow)
        {
            if (grid_movement)
            {
                Gas_Constant = config->GetGas_ConstantND();
                Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
                Mach_Motion = config->GetMach_Motion();
                RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
            }
            else
            {
                RefVel2 = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                    RefVel2 += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
            }
            RefDensity = solver[FLOW_SOL]->GetDensity_Inf();
            RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
            factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
        }
        
        /*--- Prepare send buffers for the conservative variables. Need to
          find the total number of conservative variables and also the
          index for their particular solution container. ---*/

        switch (Kind_Solver)
        {
            case EULER:
            case NAVIER_STOKES:
                FirstIndex = FLOW_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case RANS:
                FirstIndex = FLOW_SOL;
                SecondIndex = TURB_SOL;
                if (transition)
                    ThirdIndex = TRANS_SOL;
                else
                    ThirdIndex = NONE;
                break;
            case TNE2_EULER:
            case TNE2_NAVIER_STOKES:
                FirstIndex = TNE2_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case POISSON_EQUATION:
                FirstIndex = POISSON_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case WAVE_EQUATION:
                FirstIndex = WAVE_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case HEAT_EQUATION:
                FirstIndex = HEAT_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case LINEAR_ELASTICITY:
                FirstIndex = FEA_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case ADJ_EULER:
            case ADJ_NAVIER_STOKES:
                FirstIndex = ADJFLOW_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case ADJ_TNE2_EULER:
            case ADJ_TNE2_NAVIER_STOKES:
                FirstIndex = ADJTNE2_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            case ADJ_RANS:
                FirstIndex = ADJFLOW_SOL;
                if (config->GetFrozen_Visc())
                    SecondIndex = NONE;
                else
                    SecondIndex = ADJTURB_SOL;
                ThirdIndex = NONE;
                break;
            case LIN_EULER:
            case LIN_NAVIER_STOKES:
                FirstIndex = LINFLOW_SOL;
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
            default:
                SecondIndex = NONE;
                ThirdIndex = NONE;
                break;
        }

        nVar_First = solver[FirstIndex]->GetnVar();
        if (SecondIndex != NONE)
            nVar_Second = solver[SecondIndex]->GetnVar();
        if (ThirdIndex != NONE)
            nVar_Third = solver[ThirdIndex]->GetnVar();
        nVar_Consv = nVar_First + nVar_Second + nVar_Third;
        nVar_Total = nVar_Consv;

        if (!config->GetLow_MemoryOutput())
        {
            /*--- Add the limiters ---*/
            if (config->GetWrt_Limiters())
                nVar_Total += nVar_Consv;

            /*--- Add the residuals ---*/
            if (config->GetWrt_Residuals())
                nVar_Total += nVar_Consv;

            /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
            if (grid_movement)
            {
                iVar_GridVel = nVar_Total;
                if (geometry->GetnDim() == 2)
                    nVar_Total += 2;
                else if (geometry->GetnDim() == 3)
                    nVar_Total += 3;
            }

            /*--- Add density to the restart file ---*/
            if ((config->GetKind_Regime() == FREESURFACE))
            {
                iVar_Density = nVar_Total; nVar_Total += 1;
            }

            /*--- Add Pressure, Temperature, Cp, Mach to the restart file ---*/
            if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
            {
                iVar_PressCp = nVar_Total;
                nVar_Total += 3;
                iVar_MachMean = nVar_Total;
                nVar_Total += 1;
            }

            /*--- Add Laminar Viscosity, Skin Friction, Heat Flux, & yPlus to the restart file ---*/
            if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
            {
                iVar_Lam = nVar_Total;
                nVar_Total += 1;
                iVar_ViscCoeffs = nVar_Total;
                nVar_Total += 3;
            }

            /*--- Add Eddy Viscosity to the restart file ---*/
            if (Kind_Solver == RANS)
            {
                iVar_Eddy = nVar_Total;
                nVar_Total += 1;
            }

            /*--- Add Sharp edges to the restart file ---*/
            if (config->GetWrt_SharpEdges())
            {
                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                {
                    iVar_Sharp = nVar_Total;
                    nVar_Total += 1;
                }
            }

            if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES))
            {
                iVar_Mach = nVar_Total;
                nVar_Total++;
                iVar_Press = nVar_Total;
                nVar_Total++;
                iVar_Temp = nVar_Total;
                nVar_Total++;
                iVar_Tempv = nVar_Total;
                nVar_Total++;
            }

            if (Kind_Solver == TNE2_NAVIER_STOKES)
            {
                iVar_TempLam = nVar_Total; nVar_Total += config->GetnSpecies() + 3;
            }

            if (Kind_Solver == POISSON_EQUATION)
            {
                iVar_EF = nVar_Total; nVar_Total += geometry->GetnDim();
            }

            if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
                (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_TNE2_EULER) ||
                (Kind_Solver == ADJ_TNE2_NAVIER_STOKES))
            {
                iVar_Sens = nVar_Total;
                nVar_Total += 2;
            }

            if (Kind_Solver == LINEAR_ELASTICITY)
            {
                iVar_FEA = nVar_Total;
                nVar_Total += 2;
            }

            if (config->GetExtraOutput())
            {
                if (Kind_Solver == RANS)
                {
                    iVar_Extra = nVar_Total;
                    nVar_Extra = solver[TURB_SOL]->GetnOutputVariables();
                    nVar_Total += nVar_Extra;
                }
                if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES))
                {
                    iVar_Extra = nVar_Total;
                    nVar_Extra = solver[TNE2_SOL]->GetnVar();
                    nVar_Total += nVar_Extra;
                }
            }
        }

        Local_Halo = new int[geometry->GetnPoint()];
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
            Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

        /*--- Search all send/recv boundaries on this partition for any periodic
          nodes that were part of the original domain. We want to recover these
          for visualization purposes. ---*/
        if (Wrt_Halo)
        {
            nLocalPoint = geometry->GetnPoint();
        }
        else
        {
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            {
                if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)
                {
                    SendRecv = config->GetMarker_All_SendRecv(iMarker);
                    RecvFrom = abs(SendRecv) - 1;

                    /*--- Checking for less than or equal to the rank, because there may
                      be some periodic halo nodes that send info to the same rank. ---*/
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
                    {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) && (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
                        if (isPeriodic)
                            Local_Halo[iPoint] = false;
                    }
                }
            }

            /*--- Sum total number of nodes that belong to the domain ---*/
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
            {
                if (Local_Halo[iPoint] == false)
                    nLocalPoint++;
            }
        }
        Buffer_Send_nPoint[0] = nLocalPoint;

        /*--- Each processor sends its local number of nodes to the master. ---*/
        if (rank == MASTER_NODE)
            Buffer_Recv_nPoint = new unsigned long[size];

#ifdef HAVE_MPI
        MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
        MaxLocalPoint = nLocalPoint;
        Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif

        nBuffer_Scalar = MaxLocalPoint;

        /*--- Send and Recv buffers. ---*/
        double *Buffer_Send_Var = new double[MaxLocalPoint];
        double *Buffer_Recv_Var = NULL;

        double *Buffer_Send_Res = new double[MaxLocalPoint];
        double *Buffer_Recv_Res = NULL;

        double *Buffer_Send_Vol = new double[MaxLocalPoint];
        double *Buffer_Recv_Vol = NULL;

        unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
        unsigned long *Buffer_Recv_GlobalIndex = NULL;

        /*--- Auxiliary vectors for surface coefficients ---*/
        if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
        {
            Aux_Frict = new double[geometry->GetnPoint()];
            Aux_Heat = new double[geometry->GetnPoint()];
            Aux_yPlus = new double[geometry->GetnPoint()];
        }

        if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
            (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_TNE2_EULER) ||
            (Kind_Solver == ADJ_TNE2_NAVIER_STOKES))
        {
            Aux_Sens = new double[geometry->GetnPoint()];
        }

        /*--- Prepare the receive buffers in the master node only. ---*/
        if (rank == MASTER_NODE)
        {
            Buffer_Recv_Var = new double[size*MaxLocalPoint];
            Buffer_Recv_Res = new double[size*MaxLocalPoint];
            Buffer_Recv_Vol = new double[size*MaxLocalPoint];
            Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];

            /*--- Sum total number of nodes to be written and allocate arrays ---*/
            d_numGlobalPoints = 0;
            for (iProcessor = 0; iProcessor < size; iProcessor++)
            {
                d_numGlobalPoints += Buffer_Recv_nPoint[iProcessor];
            }

            Data = new double*[nVar_Total];
            for (iVar = 0; iVar < nVar_Total; iVar++)
            {
                Data[iVar] = new double[d_numGlobalPoints];
            }
        }

        /*--- Main communication routine. Loop over each variable that has
          been requested by the user and perform the MPI comm. Temporary
          1-D buffers are used to send the solution for each variable at all
          nodes on each partition to the master node. These are then unpacked
          by the master and sorted by global index in one large n-dim. array. ---*/
        for (iVar = 0; iVar < nVar_Consv; iVar++)
        {
            /*--- Logic for which solution class to draw from. ---*/
            jVar = iVar;
            CurrentIndex = FirstIndex;
            if ((SecondIndex != NONE) && (iVar > nVar_First - 1))
            {
                jVar = iVar - nVar_First;
                CurrentIndex = SecondIndex;
            }
            if ((SecondIndex != NONE) && (ThirdIndex != NONE) && (iVar > (nVar_First + nVar_Second - 1)))
            {
                jVar = iVar - nVar_First - nVar_Second;
                CurrentIndex = ThirdIndex;
            }

            /*--- Loop over this partition to collect the current variable ---*/
            jPoint = 0;
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
            {
                /*--- Check for halos & write only if requested ---*/
                if (!Local_Halo[iPoint] || Wrt_Halo)
                {
                    /*--- Get this variable into the temporary send buffer. ---*/
                    Buffer_Send_Var[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar);

                    if (!config->GetLow_MemoryOutput())
                    {
                        if (config->GetWrt_Limiters())
                        {
                            Buffer_Send_Vol[jPoint] = solver[CurrentIndex]->node[iPoint]->GetLimiter_Primitive(jVar);
                        }

                        if (config->GetWrt_Residuals())
                        {
                            Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes.GetBlock(iPoint, jVar);
                        }
                    }

                    /*--- Only send/recv the volumes & global indices during the first loop ---*/
                    if (iVar == 0)
                    {
                        Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
                    }

                    jPoint++;
                }
            }

            /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
            MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
            for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
            if (!config->GetLow_MemoryOutput())
            {
                if (config->GetWrt_Limiters())
                {
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
                }

                if (config->GetWrt_Residuals())
                {
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
                }
            }

            if (iVar == 0)
            {
#ifdef HAVE_MPI
                MPI_Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
                for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
            }

                /*--- The master node unpacks and sorts this variable by global index ---*/
                if (rank == MASTER_NODE)
                {
                    jPoint = 0;
                    for (iProcessor = 0; iProcessor < size; iProcessor++)
                    {
                        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                        {
                            /*--- Get global index, then loop over each variable and store ---*/
                            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];

                            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];

                            if (!config->GetLow_MemoryOutput())
                            {
                                if (config->GetWrt_Limiters())
                                {
                                    Data[iVar + nVar_Consv][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
                                }

                                if (config->GetWrt_Residuals())
                                {
                                    unsigned short ExtraIndex;
                                    ExtraIndex = nVar_Consv;
                                    if (config->GetWrt_Limiters()) ExtraIndex = 2 * nVar_Consv;
                                    Data[iVar + ExtraIndex][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                }
                            }
                            jPoint++;
                        }
                        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                        jPoint = (iProcessor + 1)*nBuffer_Scalar;
                    }
                }
            }

            if (!config->GetLow_MemoryOutput())
            {
                /*--- Additional communication routine for the grid velocity. Note that
                we are reusing the same temporary buffers from above for efficiency.
                Also, in the future more routines like this could be used to write
                an arbitrary number of additional variables to the file. ---*/

                if (grid_movement)
                {
                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0; double *Grid_Vel;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the three grid velocity components. ---*/

                            Grid_Vel = geometry->node[iPoint]->GetGridVel();
                            Buffer_Send_Var[jPoint] = Grid_Vel[0];
                            Buffer_Send_Res[jPoint] = Grid_Vel[1];
                            if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    if (geometry->GetnDim() == 3) 
                    {
                        MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    }
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
                    if (geometry->GetnDim() == 3)
                    {
                        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                            Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
                    }
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_GridVel;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                if (geometry->GetnDim() == 3)
                                    Data[iVar + 2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate the Density in Free-surface problems ---*/
                if (config->GetKind_Regime() == FREESURFACE)
                {
                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the pressure and mach variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetDensityInc();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_Density;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate Pressure, Cp, and Mach ---*/
                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                {
                    /*--- First, loop through the mesh in order to find and store the
                    value of the coefficient of pressure at any surface nodes. They
                    will be placed in an auxiliary vector and then communicated like
                    all other volumetric variables. ---*/

                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
                            if (compressible)
                            {
                                Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure();
                                Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature();
                                Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff;
                            }
                            if (incompressible || freesurface)
                            {
                                Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressureInc();
                                Buffer_Send_Res[jPoint] = 0.0;
                                Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff;
                            }
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_PressCp;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                Data[iVar + 2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate Mach---*/
                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                {
                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
                            if (compressible)
                            {
                                Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2()) / solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
                            }
                            if (incompressible || freesurface)
                            {
                                Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref() /
                                    sqrt(config->GetBulk_Modulus() / (solver[FLOW_SOL]->node[iPoint]->GetDensityInc()*config->GetDensity_Ref()));
                            }
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_MachMean;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Laminar Viscosity ---*/
                if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                {
                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
                            if (compressible)
                            {
                                Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
                            }
                            if (incompressible || freesurface)
                            {
                                Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc();
                            }
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_Lam;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }

                    /*--- Communicate skin friction, heat transfer, y+ ---*/

                    /*--- First, loop through the mesh in order to find and store the
                    value of the viscous coefficients at any surface nodes. They
                    will be placed in an auxiliary vector and then communicated like
                    all other volumetric variables. ---*/
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        Aux_Frict[iPoint] = 0.0;
                        Aux_Heat[iPoint] = 0.0;
                        Aux_yPlus[iPoint] = 0.0;
                    }
                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    {
                        if (config->GetMarker_All_Plotting(iMarker) == YES)
                        {
                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
                            {
                                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                                Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex);
                                Aux_Heat[iPoint] = solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex);
                                Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker, iVertex);
                            }
                        }
                    }

                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
                            if (compressible)
                            {
                                Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
                                Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
                                Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
                            }
                            if (incompressible || freesurface)
                            {
                                Buffer_Send_Var[jPoint] = Aux_Frict[iPoint];
                                Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
                                Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
                            }
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    MPI_Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_ViscCoeffs;

                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar + 0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                Data[iVar + 2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate the Eddy Viscosity ---*/
                if (Kind_Solver == RANS)
                {
                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the pressure and mach variables. ---*/
                            if (compressible)
                            {
                                Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();
                            }
                            if (incompressible || freesurface)
                            {
                                Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosityInc();
                            }
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_Eddy;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate the Sharp Edges ---*/
                if (config->GetWrt_SharpEdges())
                {
                    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                    {
                        /*--- Loop over this partition to collect the current variable ---*/
                        jPoint = 0;
                        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                        {
                            /*--- Check for halos & write only if requested ---*/
                            if (!Local_Halo[iPoint] || Wrt_Halo)
                            {
                                /*--- Load buffers with the pressure and mach variables. ---*/
                                Buffer_Send_Var[jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance();
                                jPoint++;
                            }
                        }

                        /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                        MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                            Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                        /*--- The master node unpacks and sorts this variable by global index ---*/
                        if (rank == MASTER_NODE)
                        {
                            jPoint = 0; iVar = iVar_Sharp;
                            for (iProcessor = 0; iProcessor < size; iProcessor++)
                            {
                                for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                                {
                                    /*--- Get global index, then loop over each variable and store ---*/
                                    iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                    Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                    jPoint++;
                                }

                                /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                                jPoint = (iProcessor + 1)*nBuffer_Scalar;
                            }
                        }
                    }
                }

                /*--- Communicate additional variables for the two-temperature solvers ---*/
                if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES))
                {
                    /*--- Mach number ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = sqrt(solver[TNE2_SOL]->node[iPoint]->GetVelocity2()) / solver[TNE2_SOL]->node[iPoint]->GetSoundSpeed();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_Mach;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }

                    /*--- Pressure ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetPressure();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_Press;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }

                    /*--- Temperature ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_Temp;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }

                    /*--- Vib-el Temperature ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetTemperature_ve();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/
#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_Tempv;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                if (Kind_Solver == TNE2_NAVIER_STOKES)
                {
                    /*--- Species diffusion coefficients ---*/
                    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
                    {
                        jPoint = 0;
                        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                        {
                            /*--- Check for halos & write only if requested ---*/
                            if (!Local_Halo[iPoint] || Wrt_Halo)
                            {
                                /*--- Load buffers with the Mach number variables. ---*/
                                Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetDiffusionCoeff()[iSpecies];
                                jPoint++;
                            }
                        }

                        /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                        MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                            Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                        /*--- The master node unpacks and sorts this variable by global index ---*/
                        if (rank == MASTER_NODE)
                        {
                            jPoint = 0;
                            iVar = iVar_TempLam + iSpecies;
                            for (iProcessor = 0; iProcessor < size; iProcessor++)
                            {
                                for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                                {
                                    /*--- Get global index, then loop over each variable and store ---*/
                                    iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                    Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                    jPoint++;
                                }

                                /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                                jPoint = (iProcessor + 1)*nBuffer_Scalar;
                            }
                        }
                    }

                    /*--- Laminar viscosity ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetLaminarViscosity();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_TempLam + config->GetnSpecies();
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }
                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }

                    /*--- Thermal conductivity ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_TempLam + config->GetnSpecies() + 1;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }

                    /*--- Vib-el Thermal conductivity ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the Mach number variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->node[iPoint]->GetThermalConductivity_ve();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0;
                        iVar = iVar_TempLam + config->GetnSpecies() + 2;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate the surface sensitivity ---*/
                if ((Kind_Solver == ADJ_EULER) ||
                    (Kind_Solver == ADJ_NAVIER_STOKES) ||
                    (Kind_Solver == ADJ_RANS) ||
                    (Kind_Solver == ADJ_TNE2_EULER) ||
                    (Kind_Solver == ADJ_TNE2_NAVIER_STOKES))
                {
                    /*--- First, loop through the mesh in order to find and store the
                    value of the surface sensitivity at any surface nodes. They
                    will be placed in an auxiliary vector and then communicated like
                    all other volumetric variables. ---*/

                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                        Aux_Sens[iPoint] = 0.0;
                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    {
                        if (config->GetMarker_All_Plotting(iMarker) == YES)
                        {
                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
                            {
                                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                                Area = 0.0;
                                for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim] * Normal[iDim];
                                Area = sqrt(Area);
                                Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex) / Area;
                            }
                        }
                    }

                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
                            Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
                            if (config->GetKind_ConvNumScheme() == SPACE_CENTERED)
                                Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(iPoint);
                            if (config->GetKind_ConvNumScheme() == SPACE_UPWIND)
                                Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0);

                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                    MPI_Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_Sens;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar + 0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                /*--- Communicate the Linear elasticity ---*/
                if (Kind_Solver == LINEAR_ELASTICITY)
                {
                    /*--- Loop over this partition to collect the current variable ---*/
                    jPoint = 0;
                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    {
                        /*--- Check for halos & write only if requested ---*/
                        if (!Local_Halo[iPoint] || Wrt_Halo)
                        {
                            /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
                            Buffer_Send_Var[jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
                            Buffer_Send_Res[jPoint] = solver[FEA_SOL]->node[iPoint]->GetFlow_Pressure();
                            jPoint++;
                        }
                    }

                    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                    MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                    /*--- The master node unpacks and sorts this variable by global index ---*/
                    if (rank == MASTER_NODE)
                    {
                        jPoint = 0; iVar = iVar_FEA;
                        for (iProcessor = 0; iProcessor < size; iProcessor++)
                        {
                            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                            {
                                /*--- Get global index, then loop over each variable and store ---*/
                                iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
                                jPoint++;
                            }

                            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                            jPoint = (iProcessor + 1)*nBuffer_Scalar;
                        }
                    }
                }

                if (config->GetExtraOutput())
                {
                    for (jVar = 0; jVar < nVar_Extra; jVar++)
                    {
                        /*--- Loop over this partition to collect the current variable ---*/
                        jPoint = 0;
                        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                        {
                            /*--- Check for halos & write only if requested ---*/
                            if (!Local_Halo[iPoint] || Wrt_Halo)
                            {
                                /*--- Get this variable into the temporary send buffer. ---*/
                                if (Kind_Solver == RANS)
                                {
                                    Buffer_Send_Var[jPoint] = solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra + jVar];
                                }
                                if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES))
                                {
                                    Buffer_Send_Var[jPoint] = solver[TNE2_SOL]->OutputVariables[iPoint*nVar_Extra + jVar];
                                }
                                jPoint++;
                            }
                        }

                        /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
                        MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
                        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
                            Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

                        /*--- The master node unpacks and sorts this variable by global index ---*/
                        if (rank == MASTER_NODE)
                        {
                            jPoint = 0; iVar = iVar_Extra;
                            for (iProcessor = 0; iProcessor < size; iProcessor++)
                            {
                                for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++)
                                {
                                    /*--- Get global index, then loop over each variable and store ---*/
                                    iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
                                    Data[iVar + jVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
                                    jPoint++;
                                }

                                /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
                                jPoint = (iProcessor + 1)*nBuffer_Scalar;
                            }
                        }
                    }
                }
            }

            /*--- Immediately release the temporary buffers. ---*/
            delete[] Buffer_Send_Var;
            delete[] Buffer_Send_Res;
            delete[] Buffer_Send_Vol;
            delete[] Buffer_Send_GlobalIndex;
            if (rank == MASTER_NODE)
            {
                delete[] Buffer_Recv_Var;
                delete[] Buffer_Recv_Res;
                delete[] Buffer_Recv_Vol;
                delete[] Buffer_Recv_GlobalIndex;
            }

            /*--- Release memory needed for surface coefficients ---*/
            delete[] Local_Halo;

            if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
            {
                delete[] Aux_Frict;
                delete[] Aux_Heat;
                delete[] Aux_yPlus;
            }

            if ((Kind_Solver == ADJ_EULER) ||
                (Kind_Solver == ADJ_NAVIER_STOKES) ||
                (Kind_Solver == ADJ_RANS) ||
                (Kind_Solver == ADJ_TNE2_EULER) ||
                (Kind_Solver == ADJ_TNE2_NAVIER_STOKES))
            {
                delete[] Aux_Sens;
            }
        }


} // end ARIES namespace

