
   
        /*!
         * \brief Merge the solution into a data structure used for output file writing.
         */
        void Cae_Output::MergeBaselineSolution(Config_p config, Geom_p geometry, Solver_p solver, unsigned short val_iZone)
        {
            /*--- Local variables needed on all processors ---*/
            unsigned short iVar;
            unsigned long iPoint = 0, jPoint = 0;

            nVar_Total = config->fields.size() - 1;

            /*--- Merge the solution either in serial or parallel. ---*/

#ifndef HAVE_MPI

            /*--- In serial, the single process has access to all solution data,
            so it is simple to retrieve and store inside Solution_Data. ---*/
            unsigned short iMarker;
            unsigned long iVertex, nTotalPoints = 0;
            int SendRecv;

            /*--- First, create a structure to locate any periodic halo nodes ---*/
            int *Local_Halo = new int[geometry->GetnPoint()];
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            {
                if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) 
                {
                    SendRecv = config->GetMarker_All_SendRecv(iMarker);
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) 
                    {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                            (SendRecv < 0))
                        {
                            Local_Halo[iPoint] = false;
                        }
                    }
                }
            }

            /*--- Total number of points in the mesh (this might include periodic points). ---*/
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                if (!Local_Halo[iPoint]) nTotalPoints++;

            d_numGlobalPoints = nTotalPoints;
            Data = new double*[nVar_Total];
            for (iVar = 0; iVar < nVar_Total; iVar++) 
            {
                Data[iVar] = new double[d_numGlobalPoints];
            }

            /*--- Loop over all points in the mesh, but only write data
            for nodes in the domain (ignore periodic halo nodes). ---*/
            jPoint = 0;
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
            {
                if (!Local_Halo[iPoint]) 
                {
                    /*--- Solution (first, and second system of equations) ---*/
                    unsigned short jVar = 0;
                    for (iVar = 0; iVar < nVar_Total; iVar++) 
                    {
                        Data[jVar][jPoint] = solver->node[iPoint]->GetSolution(iVar);
                        jVar++;
                    }
                }

                /*--- Increment jPoint as the counter. We need this because iPoint
                may include halo nodes that we skip over during this loop. ---*/
                jPoint++;
            }
#else

            /*--- MPI preprocessing ---*/
            int rank, nProcessor, iProcessor;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            /*--- Local variables needed for merging with MPI ---*/
            unsigned long iVertex, iMarker;

            int SendRecv, RecvFrom;

            unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
            unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
            unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

            int *Local_Halo = new int[geometry->GetnPoint()];
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

            bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

            /*--- Search all send/recv boundaries on this partition for any periodic
            nodes that were part of the original domain. We want to recover these
            for visualization purposes. ---*/
            if (Wrt_Halo) 
            {
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
            Buffer_Send_nPoint[0] = nLocalPoint;

            if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];

            MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
            MPI_Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

            nBuffer_Scalar = MaxLocalPoint;

            /*--- Send and Recv buffers. ---*/

            double *Buffer_Send_Var = new double[MaxLocalPoint];
            double *Buffer_Recv_Var = NULL;

            unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
            unsigned long *Buffer_Recv_GlobalIndex = NULL;

            /*--- Prepare the receive buffers in the master node only. ---*/
            if (rank == MASTER_NODE) {

                Buffer_Recv_Var = new double[nProcessor*MaxLocalPoint];
                Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

                /*--- Sum total number of nodes to be written and allocate arrays ---*/
                d_numGlobalPoints = 0;
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                    d_numGlobalPoints += Buffer_Recv_nPoint[iProcessor];
                }
                Data = new double*[nVar_Total];
                for (iVar = 0; iVar < nVar_Total; iVar++) {
                    Data[iVar] = new double[d_numGlobalPoints];
                }

            }

            /*--- Main communication routine. Loop over each variable that has
            been requested by the user and perform the MPI comm. Temporary
            1-D buffers are used to send the solution for each variable at all
            nodes on each partition to the master node. These are then unpacked
            by the master and sorted by global index in one large n-dim. array. ---*/

            for (iVar = 0; iVar < nVar_Total; iVar++) {

                /*--- Loop over this partition to collect the current variable ---*/
                jPoint = 0;
                for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

                    /*--- Check for halos and write only if requested ---*/
                    if (!Local_Halo[iPoint] || Wrt_Halo) {

                        /*--- Get this variable into the temporary send buffer. ---*/
                        Buffer_Send_Var[jPoint] = solver->node[iPoint]->GetSolution(iVar);

                        /*--- Only send/recv the volumes & global indices during the first loop ---*/
                        if (iVar == 0) {
                            Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
                        }
                        jPoint++;
                    }
                }

                /*--- Gather the data on the master node. ---*/

                MPI_Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
                if (iVar == 0) {
                    MPI_Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
                }

                /*--- The master node unpacks and sorts this variable by global index ---*/
                if (rank == MASTER_NODE) {
                    jPoint = 0;
                    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

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

            /*--- Immediately release the temporary buffers. ---*/

            delete[] Buffer_Send_Var;
            delete[] Buffer_Send_GlobalIndex;
            if (rank == MASTER_NODE) {
                delete[] Buffer_Recv_Var;
                delete[] Buffer_Recv_GlobalIndex;
            }

#endif
            delete[] Local_Halo;
        }

        /*!
         * \brief Write a native ARIES restart file.
         */
        void Cae_Output::SetRestart(Config_p config, Geom_p geometry, Solver_p *solver, unsigned short val_iZone) 
        {
            /*--- Local variables ---*/
            unsigned short Kind_Solver = config->GetKind_Solver();
            unsigned short iVar, iDim, nDim = geometry->GetnDim();
            unsigned long iPoint, iExtIter = config->GetExtIter();
            bool grid_movement = config->GetGrid_Movement();
            ofstream restart_file;
            string filename;

            /*--- Retrieve filename from config ---*/
            if (config->GetAdjoint()) 
            {
                filename = config->GetRestart_AdjFileName();
                filename = config->GetObjFunc_Extension(filename);
            }
            else 
            {
                filename = config->GetRestart_FlowFileName();
            }

            /*--- Unsteady problems require an iteration number to be appended. ---*/
            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) 
            {
                filename = config->GetUnsteady_FileName(filename, int(val_iZone));
            }
            else if (config->GetWrt_Unsteady()) 
            {
                filename = config->GetUnsteady_FileName(filename, int(iExtIter));
            }

            /*--- Open the restart file and write the solution. ---*/
            restart_file.open(filename.c_str(), ios::out);
            restart_file.precision(15);

            /*--- Write the header line based on the particular solver ----*/
            restart_file << "\"PointID\"";

            /*--- Mesh coordinates are always written to the restart first ---*/

            if (nDim == 2) 
            {
                restart_file << "\t\"x\"\t\"y\"";
            }
            else 
            {
                restart_file << "\t\"x\"\t\"y\"\t\"z\"";
            }

            for (iVar = 0; iVar < nVar_Consv; iVar++) 
            {
                restart_file << "\t\"Conservative_" << iVar + 1 << "\"";
            }

            if (!config->GetLow_MemoryOutput()) 
            {
                if (config->GetWrt_Limiters()) 
                {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) 
                    {
                        restart_file << "\t\"Limiter_" << iVar + 1 << "\"";
                    }
                }
                if (config->GetWrt_Residuals()) 
                {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) 
                    {
                        restart_file << "\t\"Residual_" << iVar + 1 << "\"";
                    }
                }

                /*--- Mesh velocities for dynamic mesh cases ---*/

                if (grid_movement) 
                {
                    if (nDim == 2) 
                    {
                        restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"";
                    }
                    else 
                    {
                        restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"\t\"Grid_Velz\"";
                    }
                }

                /*--- Solver specific output variables ---*/
                if (config->GetKind_Regime() == FREESURFACE) 
                {
                    restart_file << "\t\"Density\"";
                }

                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) 
                {
                    restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"Pressure_Coefficient\"\t\"Mach\"";
                }

                if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) 
                {
                    restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient\"\t\"Heat_Flux\"\t\"Y_Plus\"";
                }

                if (Kind_Solver == RANS) 
                {
                    restart_file << "\t\"Eddy_Viscosity\"";
                }

                if (config->GetWrt_SharpEdges()) 
                {
                    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                    {
                        restart_file << "\t\"Sharp_Edge_Dist\"";
                    }
                }

                if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) 
                {
                    restart_file << "\t\"Mach\"\t\"Pressure\"\t\"Temperature\"\t\"Temperature_ve\"";
                }

                if (Kind_Solver == TNE2_NAVIER_STOKES) 
                {
                    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
                        restart_file << "\t\"DiffusionCoeff_" << iSpecies << "\"";
                    restart_file << "\t\"Laminar_Viscosity\"\t\"ThermConductivity\"\t\"ThermConductivity_ve\"";
                }

                if (Kind_Solver == POISSON_EQUATION) 
                {
                    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                        restart_file << "\t\"poissonField_" << iDim + 1 << "\"";
                }

                if ((Kind_Solver == ADJ_EULER) ||
                    (Kind_Solver == ADJ_NAVIER_STOKES) ||
                    (Kind_Solver == ADJ_RANS) ||
                    (Kind_Solver == ADJ_TNE2_EULER) ||
                    (Kind_Solver == ADJ_TNE2_NAVIER_STOKES))
                {
                    restart_file << "\t\"Surface_Sensitivity\"\t\"Solution_Sensor\"";
                }

                if (Kind_Solver == LINEAR_ELASTICITY) 
                {
                    restart_file << "\t\"Von_Mises_Stress\"\t\"Flow_Pressure\"";
                }

                if (config->GetExtraOutput()) 
                {
                    string *headings = NULL;
                    //if (Kind_Solver == RANS) {
                    headings = solver[TURB_SOL]->OutputHeadingNames;
                    //}

                    for (iVar = 0; iVar < nVar_Extra; iVar++)
                    {
                        if (headings == NULL)
                        {
                            restart_file << "\t\"ExtraOutput_" << iVar + 1 << "\"";
                        }
                        else
                        {
                            restart_file << "\t\"" << headings[iVar] << "\"";
                        }
                    }
                }
            }

            restart_file << endl;

            /*--- Write the restart file ---*/

            for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) 
            {
                /*--- Index of the point ---*/
                restart_file << iPoint << "\t";

                /*--- Write the grid coordinates first ---*/
                for (iDim = 0; iDim < nDim; iDim++) 
                {
                    restart_file << scientific << Coords[iDim][iPoint] << "\t";
                }

                /*--- Loop over the variables and write the values to file ---*/
                for (iVar = 0; iVar < nVar_Total; iVar++) 
                {
                    restart_file << scientific << Data[iVar][iPoint] << "\t";
                }
                restart_file << endl;
            }

            restart_file.close();
        }

        /*!
         * \brief Deallocate temporary memory needed for merging and writing coordinates.
         */
        void Cae_Output::DeallocateCoordinates(Config_p config, Geom_p geometry) 
        {
            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables and initialization ---*/
            unsigned short iDim, nDim = geometry->GetnDim();

            /*--- The master node alone owns all data found in this routine. ---*/
            if (rank == MASTER_NODE) 
            {
                /*--- Deallocate memory for coordinate data ---*/
                for (iDim = 0; iDim < nDim; iDim++) 
                {
                    delete[] Coords[iDim];
                }
                delete[] Coords;
            }
        }

        /*!
         * \brief Deallocate temporary memory needed for merging and writing connectivity.
         */
        void Cae_Output::DeallocateConnectivity(Config_p config, Geom_p geometry, bool surf_sol) 
        {
            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
            /*--- The master node alone owns all data found in this routine. ---*/
            if (rank == MASTER_NODE) 
            {
                /*--- Deallocate memory for connectivity data ---*/
                if (surf_sol) 
                {
                    if (d_numGlobalBoundLines > 0) 
                        delete[] Conn_Line;
                    if (d_numGlobalBoundTrias > 0) 
                        delete[] Conn_BoundTria;
                    if (d_numGlobalBoundQuads > 0) 
                        delete[] Conn_BoundQuad;
                }
                else 
                {
                    if (d_numGlobalTrias > 0) 
                        delete[] Conn_Tria;
                    if (d_numGlobalQuads > 0) 
                        delete[] Conn_Quad;
                    if (d_numGlobalTetrs > 0) 
                        delete[] Conn_Tetr;
                    if (d_numGlobalHexas > 0) 
                        delete[] Conn_Hexa;
                    if (d_numGlobalPriss > 0) 
                        delete[] Conn_Pris;
                    if (d_numGlobalPyras > 0) 
                        delete[] Conn_Pyra;
                }
            }
        }

        /*!
         * \brief Deallocate temporary memory needed for merging and writing solution variables.
         */
        void Cae_Output::DeallocateSolution(Config_p config, Geom_p geometry) 
        {

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- The master node alone owns all data found in this routine. ---*/
            if (rank == MASTER_NODE) 
            {
                /*--- Deallocate memory for solution data ---*/
                for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) 
                {
                    delete[] Data[iVar];
                }
                delete[] Data;
            }
        }

        /*!
         * \brief Write the header of the history file.
         */
        void Cae_Output::SetConvHistory_Header(ofstream *ConvHist_file, Config_p config)
        {
            char cstr[200], buffer[50], turb_resid[1000];
            unsigned short iMarker, iMarker_Monitoring, iSpecies;
            string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff;

            bool rotating_frame = config->GetRotating_Frame();
            bool aeroelastic = config->GetAeroelastic_Simulation();
            bool equiv_area = config->GetEquivArea();
            bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS));
            bool frozen_turb = config->GetFrozen_Visc();
            bool freesurface = (config->GetKind_Regime() == FREESURFACE);
            bool inv_design = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());
            bool output_1d = config->GetWrt_1D_Output();
            bool output_per_surface = false;
            bool output_massflow = (config->GetKind_ObjFunc() == MASS_FLOW_RATE);
            if (config->GetnMarker_Monitoring() > 1) 
                output_per_surface = true;

            bool isothermal = false;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if ((config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC))
                    isothermal = true;

            /*--- Write file name with extension ---*/
            string filename = config->GetConv_FileName();
            strcpy(cstr, filename.data());

            if (config->GetWrt_Unsteady() && config->GetRestart())
            {
                long iExtIter = config->GetUnst_RestartIter();
                if (int(iExtIter) < 10) 
                    sprintf(buffer, "_0000%d", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) 
                    sprintf(buffer, "_000%d", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) 
                    sprintf(buffer, "_00%d", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) 
                    sprintf(buffer, "_0%d", int(iExtIter));
                if (int(iExtIter) >= 10000) 
                    sprintf(buffer, "_%d", int(iExtIter));
                strcat(cstr, buffer);
            }

            if ((config->GetOutput_FileFormat() == TECPLOT) ||
                (config->GetOutput_FileFormat() == FIELDVIEW)) 
                sprintf(buffer, ".dat");
            else if ((config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
                (config->GetOutput_FileFormat() == FIELDVIEW_BINARY))  
                sprintf(buffer, ".plt");
            else if (config->GetOutput_FileFormat() == PARAVIEW)  
                sprintf(buffer, ".csv");
            strcat(cstr, buffer);

            ConvHist_file->open(cstr, ios::out);
            ConvHist_file->precision(15);

            /*--- Begin of the header ---*/
            char begin[] = "\"Iteration\"";

            /*--- Header for the coefficients ---*/
            char flow_coeff[] = ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
            char heat_coeff[] = ",\"HeatFlux_Total\",\"HeatFlux_Maximum\"";
            char equivalent_area_coeff[] = ",\"CEquivArea\",\"CNearFieldOF\"";
            char rotating_frame_coeff[] = ",\"CMerit\",\"CT\",\"CQ\"";
            char free_surface_coeff[] = ",\"CFreeSurface\"";
            char wave_coeff[] = ",\"CWave\"";
            char fea_coeff[] = ",\"CFEA\"";
            char adj_coeff[] = ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
            char oneD_outputs[] = ",\"Avg_TotalPress\",\"Avg_Mach\",\"Avg_Temperature\",\"MassFlowRate\",\"FluxAvg_Pressure\",\"FluxAvg_Density\",\"FluxAvg_Velocity\",\"FluxAvg_Enthalpy\"";
            char Cp_inverse_design[] = ",\"Cp_Diff\"";
            char Heat_inverse_design[] = ",\"HeatFlux_Diff\"";
            char mass_flow_rate[] = ",\"MassFlowRate\"";


            /* Find the markers being monitored and create a header for them */
            for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) 
            {
                Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
                monitoring_coeff += ",\"CLift_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CDrag_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CSideForce_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CL/CD_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CFx_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CFy_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CFz_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CMx_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CMy_" + Monitoring_Tag + "\"";
                monitoring_coeff += ",\"CMz_" + Monitoring_Tag + "\"";
                aeroelastic_coeff += ",\"plunge_" + Monitoring_Tag + "\"";
                aeroelastic_coeff += ",\"pitch_" + Monitoring_Tag + "\"";
            }

            /*--- Header for the residuals ---*/
            char flow_resid[] = ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
            char adj_flow_resid[] = ",\"Res_AdjFlow[0]\",\"Res_AdjFlow[1]\",\"Res_AdjFlow[2]\",\"Res_AdjFlow[3]\",\"Res_AdjFlow[4]\"";
            switch (config->GetKind_Turb_Model()) 
            {
            case SA:	   
                sprintf(turb_resid, ",\"Res_Turb[0]\""); break;
            case SA_NEG: 
                sprintf(turb_resid, ",\"Res_Turb[0]\""); break;
            case ML:	  
                sprintf(turb_resid, ",\"Res_Turb[0]\""); break;
            case SST:   	
                sprintf(turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
            }
            char adj_turb_resid[] = ",\"Res_AdjTurb[0]\"";
            char levelset_resid[] = ",\"Res_LevelSet\"";
            char adj_levelset_resid[] = ",\"Res_AdjLevelSet\"";
            char wave_resid[] = ",\"Res_Wave[0]\",\"Res_Wave[1]\"";
            char fea_resid[] = ",\"Res_FEA\"";
            char heat_resid[] = ",\"Res_Heat\"";

            /*--- End of the header ---*/
            char end[] = ",\"Linear_Solver_Iterations\",\"CFL_Number\",\"Time(min)\"\n";

            if ((config->GetOutput_FileFormat() == TECPLOT) ||
                (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
                (config->GetOutput_FileFormat() == FIELDVIEW) ||
                (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) 
            {
                ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
                ConvHist_file[0] << "VARIABLES = ";
            }

            /*--- Write the header, case depending ---*/
            switch (config->GetKind_Solver()) 
            {
            case EULER: 
            case NAVIER_STOKES: 
            case RANS:
                ConvHist_file[0] << begin << flow_coeff;
                if (isothermal) 
                    ConvHist_file[0] << heat_coeff;
                if (equiv_area) 
                    ConvHist_file[0] << equivalent_area_coeff;
                if (inv_design)
                {
                    ConvHist_file[0] << Cp_inverse_design;
                    if (isothermal) 
                        ConvHist_file[0] << Heat_inverse_design;
                }
                if (rotating_frame) 
                    ConvHist_file[0] << rotating_frame_coeff;
                ConvHist_file[0] << flow_resid;
                if (turbulent) 
                    ConvHist_file[0] << turb_resid;
                if (aeroelastic)
                    ConvHist_file[0] << aeroelastic_coeff;
                if (output_per_surface)
                    ConvHist_file[0] << monitoring_coeff;
                if (output_1d)
                    ConvHist_file[0] << oneD_outputs;
                if (output_massflow) 
                    ConvHist_file[0] << mass_flow_rate;
                ConvHist_file[0] << end;
                if (freesurface) 
                {
                    ConvHist_file[0] << begin << flow_coeff << free_surface_coeff;
                    ConvHist_file[0] << flow_resid << levelset_resid << end;
                }
                break;
            case TNE2_EULER: 
            case TNE2_NAVIER_STOKES:
                ConvHist_file[0] << begin << flow_coeff;
                if (config->GetKind_Solver() == TNE2_NAVIER_STOKES) 
                {
                    ConvHist_file[0] << heat_coeff;
                    if (inv_design) 
                        ConvHist_file[0] << Heat_inverse_design;
                }
                for (iSpecies = 0; iSpecies < config->GetnSpecies() + 5; iSpecies++)
                    ConvHist_file[0] << ",\"Residual[" << iSpecies << "]\"";
                ConvHist_file[0] << end;
                break;
            case ADJ_EULER: 
            case ADJ_NAVIER_STOKES: 
            case ADJ_RANS:
            case ADJ_TNE2_EULER: 
            case ADJ_TNE2_NAVIER_STOKES:
                ConvHist_file[0] << begin << adj_coeff << adj_flow_resid;
                if ((turbulent) && (!frozen_turb)) 
                    ConvHist_file[0] << adj_turb_resid;
                ConvHist_file[0] << end;
                if (freesurface) 
                {
                    ConvHist_file[0] << begin << adj_coeff << adj_flow_resid << adj_levelset_resid << end;
                }
                break;
            case WAVE_EQUATION:
                ConvHist_file[0] << begin << wave_coeff;
                ConvHist_file[0] << wave_resid << end;
                break;

            case HEAT_EQUATION:
                ConvHist_file[0] << begin << heat_coeff;
                ConvHist_file[0] << heat_resid << end;
                break;

            case LINEAR_ELASTICITY:
                ConvHist_file[0] << begin << fea_coeff;
                ConvHist_file[0] << fea_resid << end;
                break;
            }

            if (config->GetOutput_FileFormat() == TECPLOT ||
                config->GetOutput_FileFormat() == TECPLOT_BINARY ||
                config->GetOutput_FileFormat() == FIELDVIEW ||
                config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
                ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
            }
        }

        /*!
         * \brief Write the history file and the convergence on the screen for serial computations.
         */
        void Cae_Output::SetConvHistory_Body(ofstream *ConvHist_file, Geom_p **geometry, Solver_p ***solver_container, Config_p *config, INTE::INTE_Integration ***integration, bool DualTime_Iteration, double timeused, unsigned short val_iZone) 
        {
            bool output_1d = config[val_iZone]->GetWrt_1D_Output();
            bool output_massflow = (config[val_iZone]->GetKind_ObjFunc() == MASS_FLOW_RATE);
            unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();

            int rank;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = MASTER_NODE;
#endif

            /*--- If 1-D outputs requested, calculated them. Requires info from all nodes,
            Get area-averaged and flux-averaged values at the specified surface ---*/
            if (output_1d) 
            {
                switch (config[val_iZone]->GetKind_Solver()) 
                {
                case EULER:                   
                case NAVIER_STOKES:                  
                case RANS:
                case FLUID_STRUCTURE_EULER:   
                case FLUID_STRUCTURE_NAVIER_STOKES:  
                case FLUID_STRUCTURE_RANS:
                case ADJ_EULER:               
                case ADJ_NAVIER_STOKES:          
                case ADJ_RANS:
                    OneDimensionalOutput(solver_container[val_iZone][FinestMesh][FLOW_SOL], geometry[val_iZone][FinestMesh], config[val_iZone]);
                    break;
                }
            }
            if (output_massflow) 
            {
                switch (config[val_iZone]->GetKind_Solver()) 
                {
                case EULER:                
                case NAVIER_STOKES:              
                case RANS:
                case FLUID_STRUCTURE_EULER: 
                case FLUID_STRUCTURE_NAVIER_STOKES: 
                case FLUID_STRUCTURE_RANS:
                case ADJ_EULER:           
                case ADJ_NAVIER_STOKES:   
                case ADJ_RANS:
                    SetMassFlowRate(solver_container[val_iZone][FinestMesh][FLOW_SOL], geometry[val_iZone][FinestMesh], config[val_iZone]);
                    break;
                }
            }

            /*--- Output using only the master node ---*/
            if (rank == MASTER_NODE)
            {
                unsigned long iIntIter = config[val_iZone]->GetIntIter();
                unsigned long iExtIter = config[val_iZone]->GetExtIter();

                /*--- WARNING: These buffers have hard-coded lengths. Note that you
                may have to adjust them to be larger if adding more entries. ---*/
                char begin[1000], direct_coeff[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000],
                    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
                    adj_turb_resid[1000], resid_aux[1000], levelset_resid[1000], adj_levelset_resid[1000], wave_coeff[1000],
                    heat_coeff[1000], fea_coeff[1000], wave_resid[1000], heat_resid[1000], fea_resid[1000], end[1000],
                    oneD_outputs[1000], massflow_outputs[1000];

                double dummy = 0.0, *Coord;
                unsigned short iVar, iMarker, iMarker_Monitoring;

                unsigned long LinSolvIter = 0, iPointMaxResid;
                double timeiter = timeused / double(iExtIter + 1);

                unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
                unsigned short nSpecies = config[val_iZone]->GetnSpecies();

                bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
                bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
                bool freesurface = (config[val_iZone]->GetKind_Regime() == FREESURFACE);

                bool rotating_frame = config[val_iZone]->GetRotating_Frame();
                bool aeroelastic = config[val_iZone]->GetAeroelastic_Simulation();
                bool equiv_area = config[val_iZone]->GetEquivArea();
                bool inv_design = (config[val_iZone]->GetInvDesign_Cp() || config[val_iZone]->GetInvDesign_HeatFlux());
                bool transition = (config[val_iZone]->GetKind_Trans_Model() == LM);
                bool isothermal = false;
                for (iMarker = 0; iMarker < config[val_iZone]->GetnMarker_All(); iMarker++)
                {
                    if ((config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
                        (config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC) ||
                        (config[val_iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC))
                        isothermal = true;
                }
                bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                    (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
                bool adjoint = config[val_iZone]->GetAdjoint();
                bool fluid_structure = ((config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) ||
                    (config[val_iZone]->GetKind_Solver() == FLUID_STRUCTURE_RANS));
                bool wave = (config[val_iZone]->GetKind_Solver() == WAVE_EQUATION);
                bool heat = (config[val_iZone]->GetKind_Solver() == HEAT_EQUATION);
                bool fea = (config[val_iZone]->GetKind_Solver() == LINEAR_ELASTICITY);
                bool TNE2 = ((config[val_iZone]->GetKind_Solver() == TNE2_EULER) || (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) ||
                    (config[val_iZone]->GetKind_Solver() == ADJ_TNE2_EULER) || (config[val_iZone]->GetKind_Solver() == ADJ_TNE2_NAVIER_STOKES));
                bool flow = (config[val_iZone]->GetKind_Solver() == EULER) || (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
                    (config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_EULER) ||
                    (config[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS);

                bool output_per_surface = false;
                if (config[val_iZone]->GetnMarker_Monitoring() > 1) output_per_surface = true;

                /*--- Initialize variables to store information from all domains (direct solution) ---*/
                double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
                    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
                    Total_CT = 0.0, Total_CQ = 0.0, Total_CFreeSurface = 0.0, Total_CWave = 0.0, Total_CHeat = 0.0, Total_CpDiff = 0.0, Total_HeatFluxDiff = 0.0,
                    Total_CFEA = 0.0, Total_Heat = 0.0, Total_MaxHeat = 0.0, Total_Mdot = 0.0;
                double OneD_AvgStagPress = 0.0, OneD_AvgMach = 0.0, OneD_AvgTemp = 0.0, OneD_MassFlowRate = 0.0,
                    OneD_FluxAvgPress = 0.0, OneD_FluxAvgDensity = 0.0, OneD_FluxAvgVelocity = 0.0, OneD_FluxAvgEntalpy = 0.0;

                /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
                double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
                double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;

                /*--- Residual arrays ---*/
                double *residual_flow = NULL,
                    *residual_turbulent = NULL,
                    *residual_transition = NULL,
                    *residual_TNE2 = NULL,
                    *residual_levelset = NULL;
                double *residual_adjflow = NULL,
                    *residual_adjturbulent = NULL,
                    *residual_adjTNE2 = NULL,
                    *residual_adjlevelset = NULL;
                double *residual_wave = NULL;
                double *residual_fea = NULL;
                double *residual_heat = NULL;

                /*--- Coefficients Monitored arrays ---*/
                double *aeroelastic_plunge = NULL,
                    *aeroelastic_pitch = NULL,
                    *Surface_CLift = NULL,
                    *Surface_CDrag = NULL,
                    *Surface_CSideForce = NULL,
                    *Surface_CEff = NULL,
                    *Surface_CFx = NULL,
                    *Surface_CFy = NULL,
                    *Surface_CFz = NULL,
                    *Surface_CMx = NULL,
                    *Surface_CMy = NULL,
                    *Surface_CMz = NULL;

                /*--- Initialize number of variables ---*/
                unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0,
                    nVar_Trans = 0, nVar_TNE2 = 0, nVar_Wave = 0, nVar_Heat = 0, nVar_FEA = 0,
                    nVar_AdjFlow = 0, nVar_AdjTNE2 = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0;

                /*--- Direct problem variables ---*/
                if (compressible) nVar_Flow = nDim + 2; else nVar_Flow = nDim + 1;
                if (turbulent) 
                {
                    switch (config[val_iZone]->GetKind_Turb_Model()) 
                    {
                    case SA:	   
                        nVar_Turb = 1; 
                        break;
                    case SA_NEG: 
                        nVar_Turb = 1; 
                        break;
                    case ML:	  
                        nVar_Turb = 1;
                        break;
                    case SST:   
                        nVar_Turb = 2; 
                        break;
                    }
                }
                if (transition) 
                    nVar_Trans = 2;
                if (TNE2) 
                    nVar_TNE2 = config[val_iZone]->GetnSpecies() + nDim + 2;
                if (wave) 
                    nVar_Wave = 2;
                if (fea) 
                    nVar_FEA = nDim;
                if (heat) 
                    nVar_Heat = 1;
                if (freesurface) 
                    nVar_LevelSet = 1;

                /*--- Adjoint problem variables ---*/
                if (compressible) 
                    nVar_AdjFlow = nDim + 2; 
                else 
                    nVar_AdjFlow = nDim + 1;

                if (turbulent) 
                {
                    switch (config[val_iZone]->GetKind_Turb_Model()) 
                    {
                    case SA:	  
                        nVar_AdjTurb = 1; 
                        break;
                    case SA_NEG: 
                        nVar_AdjTurb = 1;
                        break;
                    case ML:     
                        nVar_AdjTurb = 1; 
                        break;
                    case SST:    
                        nVar_AdjTurb = 2; 
                        break;
                    }
                }

                if (TNE2) 
                    nVar_AdjTNE2 = config[val_iZone]->GetnSpecies() + nDim + 2;
                if (freesurface) 
                    nVar_AdjLevelSet = 1;

                /*--- Allocate memory for the residual ---*/
                residual_flow = new double[nVar_Flow];
                residual_turbulent = new double[nVar_Turb];
                residual_transition = new double[nVar_Trans];
                residual_TNE2 = new double[nVar_TNE2];
                residual_levelset = new double[nVar_LevelSet];
                residual_wave = new double[nVar_Wave];
                residual_fea = new double[nVar_FEA];
                residual_heat = new double[nVar_Heat];

                residual_adjflow = new double[nVar_AdjFlow];
                residual_adjturbulent = new double[nVar_AdjTurb];
                residual_adjTNE2 = new double[nVar_AdjTNE2];
                residual_adjlevelset = new double[nVar_AdjLevelSet];

                /*--- Allocate memory for the coefficients being monitored ---*/
                aeroelastic_plunge = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                aeroelastic_pitch = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CLift = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CDrag = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CSideForce = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CEff = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFx = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFy = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFz = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMx = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMy = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMz = new double[config[ZONE_0]->GetnMarker_Monitoring()];

                /*--- Write information from nodes ---*/
                switch (config[val_iZone]->GetKind_Solver()) 
                {
                case EULER:                   
                case NAVIER_STOKES:                 
                case RANS:
                case FLUID_STRUCTURE_EULER: 
                case FLUID_STRUCTURE_NAVIER_STOKES:  
                case FLUID_STRUCTURE_RANS:
                case ADJ_EULER:              
                case ADJ_NAVIER_STOKES:      
                case ADJ_RANS:
                    /*--- Flow solution coefficients ---*/
                    Total_CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
                    Total_CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
                    Total_CSideForce = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
                    Total_CEff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEff();
                    Total_CMx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
                    Total_CMy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
                    Total_CMz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
                    Total_CFx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
                    Total_CFy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
                    Total_CFz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();

                    if (freesurface) 
                    {
                        Total_CFreeSurface = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFreeSurface();
                    }

                    if (isothermal) 
                    {
                        Total_Heat = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_HeatFlux();
                        Total_MaxHeat = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_MaxHeatFlux();
                    }

                    if (equiv_area) 
                    {
                        Total_CEquivArea = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
                        Total_CNearFieldOF = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();

                        /*--- Note that there is a redefinition of the nearfield based functionals ---*/
                        Total_CEquivArea = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0 - config[val_iZone]->GetWeightCd())*Total_CEquivArea;
                        Total_CNearFieldOF = config[val_iZone]->GetWeightCd()*Total_CDrag + (1.0 - config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
                    }

                    if (inv_design)
                    {
                        Total_CpDiff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
                        if (isothermal) 
                        {
                            Total_HeatFluxDiff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_HeatFluxDiff();
                        }
                    }

                    if (rotating_frame) 
                    {
                        Total_CT = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CT();
                        Total_CQ = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CQ();
                        Total_CMerit = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMerit();
                    }

                    if (aeroelastic) 
                    {
                        /*--- Look over the markers being monitored and get the desired values ---*/
                        for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) 
                        {
                            aeroelastic_plunge[iMarker_Monitoring] = config[val_iZone]->GetAeroelastic_plunge(iMarker_Monitoring);
                            aeroelastic_pitch[iMarker_Monitoring] = config[val_iZone]->GetAeroelastic_pitch(iMarker_Monitoring);
                        }
                    }

                    if (output_per_surface) 
                    {
                        /*--- Look over the markers being monitored and get the desired values ---*/
                        for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) 
                        {
                            Surface_CLift[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
                            Surface_CDrag[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
                            Surface_CSideForce[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce(iMarker_Monitoring);
                            Surface_CEff[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CEff(iMarker_Monitoring);
                            Surface_CFx[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring);
                            Surface_CFy[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring);
                            Surface_CFz[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring);
                            Surface_CMx[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
                            Surface_CMy[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
                            Surface_CMz[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
                        }
                    }

                    if (fluid_structure) 
                    {
                        Total_CFEA = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetTotal_CFEA();
                    }

                    if (output_1d) 
                    {
                        /*--- Get area-averaged and flux-averaged values at the specified surface ---*/
                        OneD_AvgStagPress = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_TotalPress();
                        OneD_AvgMach = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_Mach();
                        OneD_AvgTemp = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_Temp();
                        OneD_MassFlowRate = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_MassFlowRate();

                        OneD_FluxAvgPress = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgPress();
                        OneD_FluxAvgDensity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgDensity();
                        OneD_FluxAvgVelocity = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgVelocity();
                        OneD_FluxAvgEntalpy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_FluxAvgEntalpy();
                    }

                    /*--- Get Mass Flow at the Monitored Markers ---*/
                    if (output_massflow) 
                    {
                        Total_Mdot = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOneD_MassFlowRate();
                    }

                    /*--- Flow Residuals ---*/
                    for (iVar = 0; iVar < nVar_Flow; iVar++)
                        residual_flow[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);

                    /*--- Turbulent residual ---*/
                    if (turbulent) 
                    {
                        for (iVar = 0; iVar < nVar_Turb; iVar++)
                            residual_turbulent[iVar] = solver_container[val_iZone][FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
                    }

                    /*--- Transition residual ---*/
                    if (transition) 
                    {
                        for (iVar = 0; iVar < nVar_Trans; iVar++)
                            residual_transition[iVar] = solver_container[val_iZone][FinestMesh][TRANS_SOL]->GetRes_RMS(iVar);
                    }

                    /*--- Free Surface residual ---*/
                    if (freesurface) 
                    {
                        for (iVar = 0; iVar < nVar_LevelSet; iVar++)
                            residual_levelset[iVar] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(nDim + 1);
                    }

                    /*--- FEA residual ---*/
                    if (fluid_structure) 
                    {
                        for (iVar = 0; iVar < nVar_FEA; iVar++)
                            residual_fea[iVar] = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
                    }

                    /*--- Iterations of the linear solver ---*/
                    LinSolvIter = (unsigned long)solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetIterLinSolver();

                    /*--- Adjoint solver ---*/
                    if (adjoint)
                    {
                        /*--- Adjoint solution coefficients ---*/
                        Total_Sens_Geo = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Geo();
                        Total_Sens_Mach = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Mach();
                        Total_Sens_AoA = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_AoA();
                        Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Press();
                        Total_Sens_Temp = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_Temp();

                        /*--- Adjoint flow residuals ---*/

                        for (iVar = 0; iVar < nVar_AdjFlow; iVar++) 
                        {
                            residual_adjflow[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(iVar);
                        }

                        /*--- Adjoint turbulent residuals ---*/
                        if (turbulent)
                        {
                            if (!config[val_iZone]->GetFrozen_Visc()) 
                            {
                                for (iVar = 0; iVar < nVar_AdjTurb; iVar++)
                                    residual_adjturbulent[iVar] = solver_container[val_iZone][FinestMesh][ADJTURB_SOL]->GetRes_RMS(iVar);
                            }
                        }

                        /*--- Adjoint level set residuals ---*/
                        if (freesurface) 
                        {
                            for (iVar = 0; iVar < nVar_AdjLevelSet; iVar++)
                                residual_adjlevelset[iVar] = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(nDim + 1);
                        }
                    }
                    break;
                case TNE2_EULER:    
                case TNE2_NAVIER_STOKES:
                case ADJ_TNE2_EULER:
                case ADJ_TNE2_NAVIER_STOKES:
                    /*--- Coefficients ---*/
                    Total_CLift = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CLift();
                    Total_CDrag = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CDrag();
                    Total_CSideForce = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CSideForce();
                    Total_CEff = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CEff();
                    Total_CMx = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMx();
                    Total_CMy = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMy();
                    Total_CMz = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CMz();
                    Total_CFx = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFx();
                    Total_CFy = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFy();
                    Total_CFz = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_CFz();

                    if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) 
                    {
                        Total_Heat = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_HeatFlux();
                        Total_MaxHeat = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_MaxHeatFlux();
                        if (inv_design) 
                        {
                            Total_HeatFluxDiff = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetTotal_HeatFluxDiff();
                        }
                    }

                    /*--- Residuals ---*/
                    for (iVar = 0; iVar < nVar_TNE2; iVar++)
                        residual_TNE2[iVar] = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetRes_RMS(iVar);

                    /*--- Iterations of the linear solver ---*/
                    LinSolvIter = (unsigned long)solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetIterLinSolver();

                    /*--- Adjoint solver ---*/
                    if (adjoint) 
                    {
                        /*--- Adjoint solution coefficients ---*/
                        Total_Sens_Geo = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Geo();
                        Total_Sens_Mach = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Mach();
                        Total_Sens_AoA = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_AoA();
                        Total_Sens_Press = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Press();
                        Total_Sens_Temp = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetTotal_Sens_Temp();

                        /*--- Adjoint flow residuals ---*/
                        for (iVar = 0; iVar < nVar_AdjTNE2; iVar++)
                        {
                            residual_adjTNE2[iVar] = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetRes_RMS(iVar);
                        }
                    }
                    break;
                case WAVE_EQUATION:
                    /*--- Wave coefficients  ---*/
                    Total_CWave = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetTotal_CWave();

                    /*--- Wave Residuals ---*/
                    for (iVar = 0; iVar < nVar_Wave; iVar++) 
                    {
                        residual_wave[iVar] = solver_container[val_iZone][FinestMesh][WAVE_SOL]->GetRes_RMS(iVar);
                    }
                    break;
                case HEAT_EQUATION:
                    /*--- Heat coefficients  ---*/
                    Total_CHeat = solver_container[val_iZone][FinestMesh][HEAT_SOL]->GetTotal_CHeat();

                    /*--- Wave Residuals ---*/
                    for (iVar = 0; iVar < nVar_Heat; iVar++) 
                    {
                        residual_heat[iVar] = solver_container[val_iZone][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
                    }
                    break;
                case LINEAR_ELASTICITY:
                    /*--- FEA coefficients ---*/
                    Total_CFEA = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetTotal_CFEA();

                    /*--- Plasma Residuals ---*/
                    for (iVar = 0; iVar < nVar_FEA; iVar++) 
                    {
                        residual_fea[iVar] = solver_container[val_iZone][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
                    }
                    break;
                }

                /*--- Header frequency ---*/
                bool Unsteady = ((config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND));
                bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
                bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
                bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
                bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
                bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));

                bool write_heads;
                if (Unsteady) 
                    write_heads = (iIntIter == 0);
                else 
                    write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq() * 40)) == 0));

                if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) 
                {
                    /*--- Prepare the history file output, note that the dual
                    time output don't write to the history file ---*/
                    if (!DualTime_Iteration) 
                    {
                        /*--- Write the begining of the history file ---*/
                        sprintf(begin, "%12d", int(iExtIter));

                        /*--- Write the end of the history file ---*/
                        sprintf(end, ", %12.10f, %12.10f, %12.10f\n", double(LinSolvIter), config[val_iZone]->GetCFL(MESH_0), timeused / 60.0);

                        /*--- Write the solution and residual of the history file ---*/
                        switch (config[val_iZone]->GetKind_Solver()) 
                        {
                        case EULER: 
                        case NAVIER_STOKES: 
                        case RANS:
                        case FLUID_STRUCTURE_EULER:
                        case FLUID_STRUCTURE_NAVIER_STOKES:
                        case FLUID_STRUCTURE_RANS:
                        case ADJ_EULER: 
                        case ADJ_NAVIER_STOKES: 
                        case ADJ_RANS:
                            /*--- Direct coefficients ---*/
                            sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                                Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff);
                            if (isothermal)
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                                Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat);
                            if (equiv_area)
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
                            if (inv_design) 
                            {
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CpDiff);
                                Total_CpDiff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
                                if (isothermal) 
                                {
                                    sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat, Total_CpDiff, Total_HeatFluxDiff);
                                }
                            }
                            if (rotating_frame)
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                                Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CMerit, Total_CT, Total_CQ);

                            if (freesurface) 
                            {
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                                    Total_CFz, Total_CEff, Total_CFreeSurface);
                            }
                            if (fluid_structure)
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
                                Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);

                            if (aeroelastic) 
                            {
                                for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) 
                                {
                                    //Append one by one the surface coeff to aeroelastic coeff. (Think better way do this, maybe use string)
                                    if (iMarker_Monitoring == 0) 
                                    {
                                        sprintf(aeroelastic_coeff, ", %12.10f", aeroelastic_plunge[iMarker_Monitoring]);
                                    }
                                    else 
                                    {
                                        sprintf(surface_coeff, ", %12.10f", aeroelastic_plunge[iMarker_Monitoring]);
                                        strcat(aeroelastic_coeff, surface_coeff);
                                    }
                                    sprintf(surface_coeff, ", %12.10f", aeroelastic_pitch[iMarker_Monitoring]);
                                    strcat(aeroelastic_coeff, surface_coeff);
                                }
                            }

                            if (output_per_surface) 
                            {
                                for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                                    //Append one by one the surface coeff to monitoring coeff. (Think better way do this, maybe use string)
                                    if (iMarker_Monitoring == 0) 
                                    {
                                        sprintf(monitoring_coeff, ", %12.10f", Surface_CLift[iMarker_Monitoring]);
                                    }
                                    else 
                                    {
                                        sprintf(surface_coeff, ", %12.10f", Surface_CLift[iMarker_Monitoring]);
                                        strcat(monitoring_coeff, surface_coeff);
                                    }
                                    sprintf(surface_coeff, ", %12.10f", Surface_CDrag[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CSideForce[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CEff[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CFx[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CFy[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CFz[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CMx[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CMy[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                    sprintf(surface_coeff, ", %12.10f", Surface_CMz[iMarker_Monitoring]);
                                    strcat(monitoring_coeff, surface_coeff);
                                }
                            }

                            /*--- Flow residual ---*/
                            if (nDim == 2) 
                            {
                                if (compressible) 
                                    sprintf(flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_flow[0]), log10(residual_flow[1]), log10(residual_flow[2]), log10(residual_flow[3]), dummy);
                                if (incompressible || freesurface)
                                    sprintf(flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_flow[0]), log10(residual_flow[1]), log10(residual_flow[2]), dummy, dummy);
                            }
                            else
                            {
                                if (compressible) 
                                    sprintf(flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_flow[0]), log10(residual_flow[1]), log10(residual_flow[2]), log10(residual_flow[3]), log10(residual_flow[4]));
                                if (incompressible || freesurface) 
                                    sprintf(flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_flow[0]), log10(residual_flow[1]), log10(residual_flow[2]), log10(residual_flow[3]), dummy);
                            }

                            /*--- Turbulent residual ---*/
                            if (turbulent) 
                            {
                                switch (nVar_Turb) 
                                {
                                case 1: 
                                    sprintf(turb_resid, ", %12.10f", log10(residual_turbulent[0]));
                                    break;
                                case 2: 
                                    sprintf(turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); 
                                    break;
                                }
                            }

                            /*---- Averaged stagnation pressure at an exit ---- */
                            if (output_1d) 
                            {
                                sprintf(oneD_outputs, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", OneD_AvgStagPress, OneD_AvgMach, OneD_AvgTemp, OneD_MassFlowRate, OneD_FluxAvgPress, OneD_FluxAvgDensity, OneD_FluxAvgVelocity, OneD_FluxAvgEntalpy);
                            }
                            if (output_massflow) 
                            {
                                sprintf(massflow_outputs, ", %12.10f", Total_Mdot);
                            }

                            /*--- Transition residual ---*/
                            if (transition) 
                            {
                                sprintf(trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
                            }

                            /*--- Free surface residual ---*/
                            if (freesurface) 
                            {
                                sprintf(levelset_resid, ", %12.10f", log10(residual_levelset[0]));
                            }

                            /*--- Fluid structure residual ---*/
                            if (fluid_structure) 
                            {
                                if (nDim == 2)
                                    sprintf(levelset_resid, ", %12.10f, %12.10f, 0.0", log10(residual_fea[0]), log10(residual_fea[1]));
                                else 
                                    sprintf(levelset_resid, ", %12.10f, %12.10f, %12.10f", log10(residual_fea[0]), log10(residual_fea[1]), log10(residual_fea[2]));
                            }

                            if (adjoint) 
                            {
                                /*--- Adjoint coefficients ---*/
                                sprintf(adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);

                                /*--- Adjoint flow residuals ---*/
                                if (nDim == 2) 
                                {
                                    if (compressible) 
                                        sprintf(adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10(residual_adjflow[0]), log10(residual_adjflow[1]), log10(residual_adjflow[2]), log10(residual_adjflow[3]));
                                    if (incompressible || freesurface) 
                                        sprintf(adj_flow_resid, ", %12.10f, %12.10f, %12.10f, 0.0, 0.0", log10(residual_adjflow[0]), log10(residual_adjflow[1]), log10(residual_adjflow[2]));
                                }
                                else 
                                {
                                    if (compressible) 
                                        sprintf(adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_adjflow[0]), log10(residual_adjflow[1]), log10(residual_adjflow[2]), log10(residual_adjflow[3]), log10(residual_adjflow[4]));
                                    if (incompressible || freesurface) 
                                        sprintf(adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10(residual_adjflow[0]), log10(residual_adjflow[1]), log10(residual_adjflow[2]), log10(residual_adjflow[3]));
                                }

                                /*--- Adjoint turbulent residuals ---*/
                                if (turbulent)
                                    if (!config[val_iZone]->GetFrozen_Visc())
                                        sprintf(adj_turb_resid, ", %12.10f", log10(residual_adjturbulent[0]));

                                /*--- Adjoint free surface residuals ---*/
                                if (freesurface) sprintf(adj_levelset_resid, ", %12.10f", log10(residual_adjlevelset[0]));
                            }
                            break;
                        case TNE2_EULER:    
                        case TNE2_NAVIER_STOKES:
                        case ADJ_TNE2_EULER: 
                        case ADJ_TNE2_NAVIER_STOKES:
                            /*--- Direct coefficients ---*/
                            if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) 
                            {
                                if (!(inv_design))
                                    sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                                    Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                                    Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                                    Total_CEff, Total_Heat, Total_MaxHeat);
                                else
                                    sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                                    Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                                    Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                                    Total_CEff, Total_Heat, Total_MaxHeat, Total_HeatFluxDiff);
                            }
                            else
                                sprintf(direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                                Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx,
                                Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz,
                                Total_CEff);

                            /*--- Direct problem residual ---*/
                            for (iVar = 0; iVar < nSpecies + nDim + 2; iVar++) 
                            {
                                sprintf(resid_aux, ", %12.10f", log10(residual_TNE2[iVar]));
                                if (iVar == 0) 
                                    strcpy(flow_resid, resid_aux);
                                else 
                                    strcat(flow_resid, resid_aux);
                            }

                            if (adjoint) 
                            {
                                /*--- Adjoint coefficients ---*/
                                sprintf(adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);

                                /*--- Adjoint flow residuals ---*/
                                for (iVar = 0; iVar < nSpecies + nDim + 2; iVar++)
                                {
                                    sprintf(resid_aux, ", %12.10f", log10(residual_adjTNE2[iVar]));
                                    if (iVar == 0) 
                                        strcpy(adj_flow_resid, resid_aux);
                                    else 
                                        strcat(adj_flow_resid, resid_aux);
                                }
                            }
                            break;
                        case WAVE_EQUATION:
                            sprintf(direct_coeff, ", %12.10f", Total_CWave);
                            sprintf(wave_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_wave[0]), log10(residual_wave[1]), dummy, dummy, dummy);
                            break;
                        case HEAT_EQUATION:
                            sprintf(direct_coeff, ", %12.10f", Total_CHeat);
                            sprintf(heat_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_heat[0]), dummy, dummy, dummy, dummy);
                            break;
                        case LINEAR_ELASTICITY:
                            sprintf(direct_coeff, ", %12.10f", Total_CFEA);
                            sprintf(fea_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10(residual_fea[0]), dummy, dummy, dummy, dummy);
                            break;
                        }
                    }

                    /*--- Write the screen header---*/
                    if ((write_heads) && !(!DualTime_Iteration && Unsteady)) 
                    {
                        if (!Unsteady) 
                        {
                            switch (config[val_iZone]->GetKind_Solver()) 
                            {
                            case EULER:
                            case NAVIER_STOKES: 
                            case RANS:
                            case ADJ_EULER: 
                            case ADJ_NAVIER_STOKES: 
                            case ADJ_RANS:
                            case FLUID_STRUCTURE_EULER:  
                            case FLUID_STRUCTURE_NAVIER_STOKES:
                            case FLUID_STRUCTURE_RANS:
                                cout << endl << "---------------------- Local Time Stepping Summary ----------------------" << endl;
                                for (unsigned short iMesh = FinestMesh; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
                                    cout << "MG level: " << iMesh << " -> Min. DT: " << solver_container[val_iZone][iMesh][FLOW_SOL]->GetMin_Delta_Time() <<
                                    ". Max. DT: " << solver_container[val_iZone][iMesh][FLOW_SOL]->GetMax_Delta_Time() <<
                                    ". CFL: " << config[val_iZone]->GetCFL(iMesh) << "." << endl;
                                cout << "-------------------------------------------------------------------------" << endl;
                                break;
                            case TNE2_EULER:
                            case TNE2_NAVIER_STOKES:
                            case ADJ_TNE2_EULER:
                            case ADJ_TNE2_NAVIER_STOKES:
                                cout << endl << "Min Delta Time: " << solver_container[val_iZone][MESH_0][TNE2_SOL]->GetMin_Delta_Time() << ". Max Delta Time: " << solver_container[val_iZone][MESH_0][TNE2_SOL]->GetMax_Delta_Time() << ".";
                                break;
                            }
                        }
                        else 
                        {
                            if (flow) 
                            {
                                cout << endl << "Min DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMin_Delta_Time() <<
                                    ".Max DT: " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
                                    ".Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                            }
                            else 
                            {
                                cout << endl << "Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                            }
                        }

                        switch (config[val_iZone]->GetKind_Solver()) 
                        {
                        case EULER:                  
                        case NAVIER_STOKES:
                        case FLUID_STRUCTURE_EULER:  
                        case FLUID_STRUCTURE_NAVIER_STOKES:
                            /*--- Visualize the maximum residual ---*/
                            iPointMaxResid = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
                            Coord = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);

                            cout << endl << "----------------------- Residual Evolution Summary ----------------------" << endl;
                            cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;

                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                                if (nDim == 3) 
                                    cout << ", " << Coord[2];
                                cout << ")." << endl;
                            }
                            else 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] * 12.0 << ", " << Coord[1] * 12.0;
                                if (nDim == 3)
                                    cout << ", " << Coord[2] * 12.0;
                                cout << ")." << endl;
                            }

                            /*--- Print out the number of non-physical points and reconstructions ---*/
                            if (config[val_iZone]->GetNonphysical_Points() > 0)
                                cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
                            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
                                cout << "There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;

                            cout << "-------------------------------------------------------------------------" << endl;

                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << " ExtIter";

                            if (!fluid_structure)
                            {
                                if (incompressible) 
                                    cout << "   Res[Press]" << "     Res[Velx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
                                else if (freesurface) 
                                    cout << "   Res[Press]" << "     Res[Dist]" << "   CLift(Total)" << "     CLevelSet" << endl;
                                else if (rotating_frame && nDim == 3)
                                    cout << "     Res[Rho]" << "     Res[RhoE]" << " CThrust(Total)" << " CTorque(Total)" << endl;
                                else if (aeroelastic) 
                                    cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
                                else if (equiv_area) 
                                    cout << "     Res[Rho]" << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
                                else 
                                    cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
                            }
                            else if (fluid_structure)
                                cout << "     Res[Rho]" << "   Res[Displx]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
                            break;
                        case TNE2_EULER:  
                        case TNE2_NAVIER_STOKES:
                            /*--- Visualize the maximum residual ---*/
                            iPointMaxResid = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetPoint_Max(0);
                            Coord = solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetPoint_Max_Coord(0);
                            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][TNE2_SOL]->GetRes_Max(0)) << "." << endl;
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                                if (nDim == 3) 
                                    cout << ", " << Coord[2];
                                cout << ")." << endl;
                            }
                            else 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] * 12.0 << ", " << Coord[1] * 12.0;
                                if (nDim == 3) 
                                    cout << ", " << Coord[2] * 12.0;
                                cout << ")." << endl;
                            }

                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << " ExtIter";

                            cout << "     Res[Rho]" << "     Res[RhoE]" << "   Res[RhoEve]" << "   CDrag(Total)";
                            if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES)
                                cout << "   HeatLoad(Total)" << endl;
                            else 
                                cout << endl;
                            break;
                        case RANS: 
                        case FLUID_STRUCTURE_RANS:
                            /*--- Visualize the maximum residual ---*/
                            iPointMaxResid = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max(0);
                            Coord = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetPoint_Max_Coord(0);

                            cout << endl << "----------------------- Residual Evolution Summary ----------------------" << endl;

                            cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_Max(0)) << "." << endl;
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                                if (nDim == 3)
                                    cout << ", " << Coord[2];
                                cout << ")." << endl;
                            }
                            else 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] * 12.0 << ", " << Coord[1] * 12.0;
                                if (nDim == 3)
                                    cout << ", " << Coord[2] * 12.0;
                                cout << ")." << endl;
                            }
                            cout << "Maximum Omega " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetOmega_Max() << ", maximum Strain Rate " << solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetStrainMag_Max() << "." << endl;

                            /*--- Print out the number of non-physical points and reconstructions ---*/
                            if (config[val_iZone]->GetNonphysical_Points() > 0)
                                cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
                            if (config[val_iZone]->GetNonphysical_Reconstr() > 0)
                                cout << "There are " << config[val_iZone]->GetNonphysical_Reconstr() << " non-physical states in the upwind reconstruction." << endl;

                            cout << "-------------------------------------------------------------------------" << endl;

                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << " ExtIter";

                            if (incompressible || freesurface)
                                cout << "   Res[Press]";
                            else
                                cout << "      Res[Rho]";//, cout << "     Res[RhoE]";

                            switch (config[val_iZone]->GetKind_Turb_Model())
                            {
                            case SA:	   
                                cout << "       Res[nu]";
                                break;
                            case SA_NEG:
                                cout << "       Res[nu]";
                                break;
                            case ML:	  
                                cout << "       Res[nu]"; 
                                break;
                            case SST:	   
                                cout << "     Res[kine]" << "     Res[omega]"; 
                                break;
                            }

                            if (transition)
                            { 
                                cout << "      Res[Int]" << "       Res[Re]"; }
                            else if  (rotating_frame && nDim == 3) 
                                cout << "   CThrust(Total)" << "   CTorque(Total)" << endl;
                            else if (aeroelastic) 
                                cout << "   CLift(Total)" << "   CDrag(Total)" << "         plunge" << "          pitch" << endl;
                            else if (equiv_area)
                                cout << "   CLift(Total)" << "   CDrag(Total)" << "    CPress(N-F)" << endl;
                            else 
                                cout << "   CLift(Total)" << "   CDrag(Total)" << endl;
                            break;
                        case WAVE_EQUATION:
                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << "  ExtIter";
                            cout << "      Res[Wave]" << "   CWave(Total)" << endl;
                            break;
                        case HEAT_EQUATION:
                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else
                                cout << endl << " IntIter" << "  ExtIter";
                            cout << "      Res[Heat]" << "   CHeat(Total)" << endl;
                            break;
                        case LINEAR_ELASTICITY:
                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << "  ExtIter";

                            if (nDim == 2)
                                cout << "    Res[Displx]" << "    Res[Disply]" << "   CFEA(Total)" << endl;
                            if (nDim == 3)
                                cout << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "   CFEA(Total)" << endl;
                            break;
                        case ADJ_EULER:      
                        case ADJ_NAVIER_STOKES:
                            /*--- Visualize the maximum residual ---*/
                            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
                            Coord = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
                            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                                if (nDim == 3)
                                    cout << ", " << Coord[2];
                                cout << ")." << endl;
                            }
                            else 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] * 12.0 << ", " << Coord[1] * 12.0;
                                if (nDim == 3) 
                                    cout << ", " << Coord[2] * 12.0;
                                cout << ")." << endl;
                            }

                            /*--- Print out the number of non-physical points and reconstructions ---*/
                            if (config[val_iZone]->GetNonphysical_Points() > 0)
                                cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;

                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << "  ExtIter";

                            if (incompressible || freesurface) 
                                cout << "   Res[Psi_Press]" << "   Res[Psi_Velx]";
                            else 
                                cout << "   Res[Psi_Rho]" << "     Res[Psi_E]";
                            cout << "      Sens_Geo" << "     Sens_Mach" << endl;

                            if (freesurface) 
                            {
                                cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
                            }
                            break;
                        case ADJ_RANS:
                            /*--- Visualize the maximum residual ---*/
                            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max(0);
                            Coord = solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetPoint_Max_Coord(0);
                            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJFLOW_SOL]->GetRes_Max(0)) << "." << endl;
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                                if (nDim == 3) 
                                    cout << ", " << Coord[2];
                                cout << ")." << endl;
                            }
                            else 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] * 12.0 << ", " << Coord[1] * 12.0;
                                if (nDim == 3) 
                                    cout << ", " << Coord[2] * 12.0;
                                cout << ")." << endl;
                            }

                            /*--- Print out the number of non-physical points and reconstructions ---*/
                            if (config[val_iZone]->GetNonphysical_Points() > 0)
                                cout << "There are " << config[val_iZone]->GetNonphysical_Points() << " non-physical points in the solution." << endl;

                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << "  ExtIter";

                            if (incompressible || freesurface) 
                                cout << "     Res[Psi_Press]";
                            else 
                                cout << "     Res[Psi_Rho]";

                            if (!config[val_iZone]->GetFrozen_Visc()) 
                            {
                                cout << "      Res[Psi_nu]";
                            }
                            else 
                            {
                                if (incompressible || freesurface) 
                                    cout << "   Res[Psi_Velx]";
                                else 
                                    cout << "     Res[Psi_E]";
                            }
                            cout << "     Sens_Geo" << "    Sens_Mach" << endl;

                            if (freesurface) 
                            {
                                cout << "   Res[Psi_Press]" << "   Res[Psi_Dist]" << "    Sens_Geo" << "   Sens_Mach" << endl;
                            }
                            break;
                        case ADJ_TNE2_EULER:           
                        case ADJ_TNE2_NAVIER_STOKES:
                            /*--- Visualize the maximum residual ---*/
                            iPointMaxResid = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetPoint_Max(0);
                            Coord = solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetPoint_Max_Coord(0);
                            cout << endl << "log10[Maximum residual]: " << log10(solver_container[val_iZone][FinestMesh][ADJTNE2_SOL]->GetRes_Max(0)) << "." << endl;
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] << ", " << Coord[1];
                                if (nDim == 3) 
                                    cout << ", " << Coord[2];
                                cout << ")." << endl;
                            }
                            else 
                            {
                                cout << "Maximum residual point " << iPointMaxResid << ", located at (" << Coord[0] * 12.0 << ", " << Coord[1] * 12.0;
                                if (nDim == 3) 
                                    cout << ", " << Coord[2] * 12.0;
                                cout << ")." << endl;
                            }

                            if (!Unsteady) 
                                cout << endl << " Iter" << "    Time(s)";
                            else 
                                cout << endl << " IntIter" << "  ExtIter";

                            cout << "   Res[Psi_Rho]" << "     Res[Psi_E]" << "   Res[Psi_Eve]" << "     Sens_Geo" << "    Sens_Mach" << endl;
                            break;
                        }
                    }

                    /*--- Write the solution on the screen and history file ---*/
                    cout.precision(6);
                    cout.setf(ios::fixed, ios::floatfield);

                    if (!Unsteady) 
                    {
                        cout.width(5); cout << iExtIter;
                        cout.width(11); cout << timeiter;
                    }
                    else 
                    {
                        cout.width(8); cout << iIntIter;
                        cout.width(8); cout << iExtIter;
                    }

                    switch (config[val_iZone]->GetKind_Solver()) 
                    {
                    case EULER: 
                    case NAVIER_STOKES:
                    case FLUID_STRUCTURE_EULER: 
                    case FLUID_STRUCTURE_NAVIER_STOKES:
                        if (!DualTime_Iteration) 
                        {
                            if (compressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
                            if (incompressible) ConvHist_file[0] << begin << direct_coeff << flow_resid;
                            if (freesurface) ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
                            if (fluid_structure) ConvHist_file[0] << fea_resid;
                            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
                            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
                            if (output_1d) ConvHist_file[0] << oneD_outputs;
                            if (output_massflow) ConvHist_file[0] << massflow_outputs;
                            ConvHist_file[0] << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(13); cout << log10(residual_flow[0]);
                        if (!fluid_structure && !equiv_area) 
                        {
                            if (compressible) 
                            {
                                if (nDim == 2) 
                                { 
                                    cout.width(14); 
                                    cout << log10(residual_flow[3]);
                                }
                                else 
                                { 
                                    cout.width(14); 
                                    cout << log10(residual_flow[4]); 
                                }
                            }
                            if (incompressible) 
                            { 
                                cout.width(14); 
                                cout << log10(residual_flow[1]); 
                            }
                            if (freesurface) 
                            { 
                                cout.width(14);
                                cout << log10(residual_levelset[0]); 
                            }
                        }
                        else if (fluid_structure) 
                        { 
                            cout.width(14); 
                            cout << log10(residual_fea[0]); 
                        }

                        if (rotating_frame && nDim == 3) 
                        {
                            cout.setf(ios::scientific, ios::floatfield);
                            cout.width(15);
                            cout << Total_CT;
                            cout.width(15); 
                            cout << Total_CQ;
                            cout.unsetf(ios_base::floatfield);
                        }
                        else if (equiv_area) 
                        {
                            cout.width(15); cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); cout.width(15);
                            cout.precision(4);
                            cout.setf(ios::scientific, ios::floatfield);
                            cout << Total_CNearFieldOF;
                        }
                        else if (freesurface) 
                        { 
                            cout.width(15); 
                            cout << Total_CLift;
                            cout.width(15); 
                            cout << Total_CFreeSurface; 
                        }
                        else 
                        { 
                            cout.width(15);
                            cout << min(10000.0, max(-10000.0, Total_CLift)); 
                            cout.width(15); 
                            cout << min(10000.0, max(-10000.0, Total_CDrag)); 
                        }
                        if (aeroelastic) 
                        {
                            cout.setf(ios::scientific, ios::floatfield);
                            cout.width(15); 
                            cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
                            cout.width(15); 
                            cout << aeroelastic_pitch[0];
                            cout.unsetf(ios_base::floatfield);
                        }
                        cout << endl;
                        break;
                    case RANS:

                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid;
                            if (aeroelastic)
                                ConvHist_file[0] << aeroelastic_coeff;
                            if (output_per_surface) 
                                ConvHist_file[0] << monitoring_coeff;
                            if (output_1d)
                                ConvHist_file[0] << oneD_outputs;
                            if (output_massflow) 
                                ConvHist_file[0] << massflow_outputs;
                            ConvHist_file[0] << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);

                        if (incompressible || freesurface)
                            cout.width(13);
                        else  
                            cout.width(14);
                        cout << log10(residual_flow[0]);
                        //          else  cout.width(14),
                        //                 cout << log10(residual_flow[0]),
                        //                 cout.width(14);
                        //          if ( nDim==2 ) cout << log10(residual_flow[3]);
                        //          if ( nDim==3 ) cout << log10(residual_flow[4]);

                        switch (nVar_Turb) 
                        {
                        case 1: 
                            cout.width(14); cout << log10(residual_turbulent[0]); 
                            break;
                        case 2: 
                            cout.width(14); cout << log10(residual_turbulent[0]);
                            cout.width(15); cout << log10(residual_turbulent[1]); 
                            break;
                        }

                        if (transition) 
                        { 
                            cout.width(14); 
                            cout << log10(residual_transition[0]); 
                            cout.width(14); 
                            cout << log10(residual_transition[1]); 
                        }

                        if (rotating_frame && nDim == 3) 
                        {
                            cout.setf(ios::scientific, ios::floatfield);
                            cout.width(15); cout << Total_CT; cout.width(15);
                            cout << Total_CQ;
                            cout.unsetf(ios_base::floatfield);
                        }
                        else if (equiv_area)
                        {
                            cout.width(15); 
                            cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); cout.width(15);
                            cout.precision(4);
                            cout.setf(ios::scientific, ios::floatfield);
                            cout << Total_CNearFieldOF;
                        }
                        else 
                        { cout.width(15); cout << min(10000.0, max(-10000.0, Total_CLift)); cout.width(15); cout << min(10000.0, max(-10000.0, Total_CDrag)); }

                        if (aeroelastic) 
                        {
                            cout.setf(ios::scientific, ios::floatfield);
                            cout.width(15); cout << aeroelastic_plunge[0]; //Only output the first marker being monitored to the console.
                            cout.width(15); cout << aeroelastic_pitch[0];
                            cout.unsetf(ios_base::floatfield);
                        }
                        cout << endl;

                        if (freesurface) 
                        {
                            if (!DualTime_Iteration) 
                            {
                                ConvHist_file[0] << begin << direct_coeff << flow_resid << levelset_resid << end;
                                ConvHist_file[0].flush();
                            }

                            cout.precision(6);
                            cout.setf(ios::fixed, ios::floatfield);
                            cout.width(13); cout << log10(residual_flow[0]);
                            cout.width(14); cout << log10(residual_levelset[0]);
                            cout.width(15); cout << Total_CLift;
                            cout.width(14); cout << Total_CFreeSurface;

                            cout << endl;
                        }

                        break;

                    case TNE2_EULER: 
                    case TNE2_NAVIER_STOKES:
                        if (!DualTime_Iteration)
                        {
                            ConvHist_file[0] << begin << direct_coeff << flow_resid;
                            ConvHist_file[0] << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(13); cout << log10(residual_TNE2[0]);
                        cout.width(14); cout << log10(residual_TNE2[nSpecies + nDim]);
                        cout.width(14); cout << log10(residual_TNE2[nSpecies + nDim + 1]);
                        cout.width(15); cout << Total_CDrag;
                        if (config[val_iZone]->GetKind_Solver() == TNE2_NAVIER_STOKES) 
                        {
                            cout.precision(1);
                            cout.width(11); cout << Total_MaxHeat;
                        }
                        cout << endl;
                        break;

                    case WAVE_EQUATION:
                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << wave_coeff << wave_resid << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(14); cout << log10(residual_wave[0]);
                        cout.width(14); cout << Total_CWave;
                        cout << endl;
                        break;

                    case HEAT_EQUATION:
                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << heat_coeff << heat_resid << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(14); cout << log10(residual_heat[0]);
                        cout.width(14); cout << Total_CHeat;
                        cout << endl;
                        break;
                    case LINEAR_ELASTICITY:
                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << fea_coeff << fea_resid << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(15); cout << log10(residual_fea[0]);
                        cout.width(15); cout << log10(residual_fea[1]);
                        if (nDim == 3)
                        { 
                            cout.width(15); 
                            cout << log10(residual_fea[2]); 
                        }
                        cout.precision(4);
                        cout.setf(ios::scientific, ios::floatfield);
                        cout.width(14); cout << Total_CFEA;
                        cout << endl;
                        break;
                    case ADJ_EULER:              
                    case ADJ_NAVIER_STOKES:
                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        if (compressible)
                        {
                            cout.width(15); cout << log10(residual_adjflow[0]);
                            cout.width(15); cout << log10(residual_adjflow[nDim + 1]);
                        }
                        if (incompressible || freesurface) 
                        {
                            cout.width(17); cout << log10(residual_adjflow[0]);
                            cout.width(16); cout << log10(residual_adjflow[1]);
                        }
                        cout.precision(4);
                        cout.setf(ios::scientific, ios::floatfield);
                        cout.width(14); cout << Total_Sens_Geo;
                        cout.width(14); cout << Total_Sens_Mach;
                        cout << endl;
                        cout.unsetf(ios_base::floatfield);

                        if (freesurface) 
                        {
                            if (!DualTime_Iteration) 
                            {
                                ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid << end;
                                ConvHist_file[0].flush();
                            }

                            cout.precision(6);
                            cout.setf(ios::fixed, ios::floatfield);
                            cout.width(17); cout << log10(residual_adjflow[0]);
                            cout.width(16); cout << log10(residual_adjlevelset[0]);
                            cout.precision(3);
                            cout.setf(ios::scientific, ios::floatfield);
                            cout.width(12); cout << Total_Sens_Geo;
                            cout.width(12); cout << Total_Sens_Mach;
                            cout.unsetf(ios_base::floatfield);
                            cout << endl;
                        }
                        break;
                    case ADJ_RANS:
                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid;
                            if (!config[val_iZone]->GetFrozen_Visc())
                                ConvHist_file[0] << adj_turb_resid;
                            ConvHist_file[0] << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(17); cout << log10(residual_adjflow[0]);
                        if (!config[val_iZone]->GetFrozen_Visc()) 
                        {
                            cout.width(17); cout << log10(residual_adjturbulent[0]);
                        }
                        else 
                        {
                            if (compressible) 
                            {
                                if (geometry[val_iZone][FinestMesh]->GetnDim() == 2) 
                                { 
                                    cout.width(15); 
                                    cout << log10(residual_adjflow[3]); 
                                }
                                else 
                                { 
                                    cout.width(15); 
                                    cout << log10(residual_adjflow[4]); 
                                }
                            }
                            if (incompressible || freesurface) 
                            {
                                cout.width(15); cout << log10(residual_adjflow[1]);
                            }
                        }
                        cout.precision(4);
                        cout.setf(ios::scientific, ios::floatfield);
                        cout.width(14); cout << Total_Sens_Geo;
                        cout.width(14); cout << Total_Sens_Mach;
                        cout << endl;
                        cout.unsetf(ios_base::floatfield);
                        if (freesurface) 
                        {
                            if (!DualTime_Iteration) 
                            {
                                ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << adj_levelset_resid;
                                ConvHist_file[0] << end;
                                ConvHist_file[0].flush();
                            }

                            cout.precision(6);
                            cout.setf(ios::fixed, ios::floatfield);
                            cout.width(17); cout << log10(residual_adjflow[0]);
                            cout.width(16); cout << log10(residual_adjlevelset[0]);

                            cout.precision(4);
                            cout.setf(ios::scientific, ios::floatfield);
                            cout.width(12); cout << Total_Sens_Geo;
                            cout.width(12); cout << Total_Sens_Mach;
                            cout << endl;
                            cout.unsetf(ios_base::floatfield);
                        }
                        break;
                    case ADJ_TNE2_EULER:             
                    case ADJ_TNE2_NAVIER_STOKES:
                        if (!DualTime_Iteration) 
                        {
                            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
                            ConvHist_file[0].flush();
                        }

                        cout.precision(6);
                        cout.setf(ios::fixed, ios::floatfield);
                        cout.width(15); cout << log10(residual_adjTNE2[0]);
                        cout.width(15); cout << log10(residual_adjTNE2[nSpecies + nDim]);
                        cout.width(15); cout << log10(residual_adjTNE2[nSpecies + nDim + 1]);

                        cout.precision(4);
                        cout.setf(ios::scientific, ios::floatfield);
                        cout.width(14); cout << Total_Sens_Geo;
                        cout.width(14); cout << Total_Sens_Mach;
                        cout << endl;
                        cout.unsetf(ios_base::floatfield);
                        break;
                    }
                    cout.unsetf(ios::fixed);

                    delete[] residual_flow;
                    delete[] residual_turbulent;
                    delete[] residual_transition;
                    delete[] residual_TNE2;
                    delete[] residual_levelset;
                    delete[] residual_wave;
                    delete[] residual_fea;
                    delete[] residual_heat;

                    delete[] residual_adjflow;
                    delete[] residual_adjturbulent;
                    delete[] residual_adjTNE2;
                    delete[] residual_adjlevelset;

                    delete[] Surface_CLift;
                    delete[] Surface_CDrag;
                    delete[] Surface_CSideForce;
                    delete[] Surface_CEff;
                    delete[] Surface_CFx;
                    delete[] Surface_CFy;
                    delete[] Surface_CFz;
                    delete[] Surface_CMx;
                    delete[] Surface_CMy;
                    delete[] Surface_CMz;
                    delete[] aeroelastic_pitch;
                    delete[] aeroelastic_plunge;
                }
            }
        }

        /*!
         * \brief Set CFL numbers
         */
        void Cae_Output::SetCFL_Number(Solver_p ***solver_container, Config_p *config, unsigned short val_iZone)
        {
            double CFLFactor = 1.0, power = 1.0, CFL = 0.0, CFLMin = 0.0, CFLMax = 0.0, Div = 1.0, Diff = 0.0, MGFactor[100];
            unsigned short iMesh;

            unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
            unsigned long ExtIter = config[val_iZone]->GetExtIter();

            d_rhoResNew = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(0);

            if (d_rhoResNew < EPS)
                d_rhoResNew = EPS;
            if (d_rhoResOld < EPS)
                d_rhoResOld = d_rhoResNew;

            Div = d_rhoResOld / d_rhoResNew;
            Diff = d_rhoResNew - d_rhoResOld;

            /*--- Compute MG factor ---*/
            MGFactor[MESH_0] = 1.0;
            for (iMesh = 1; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
            {
                MGFactor[iMesh] = MGFactor[iMesh - 1] * config[val_iZone]->GetCFL(iMesh) / config[val_iZone]->GetCFL(iMesh - 1);
            }

            if (Div < 1.0)
                power = config[val_iZone]->GetCFL_AdaptParam(0);
            else
                power = config[val_iZone]->GetCFL_AdaptParam(1);

            /*--- Detect a stall in the residual ---*/
            if ((fabs(Diff) <= d_rhoResNew*1E-8) && (ExtIter != 0))
            {
                Div = 0.1;
                power = config[val_iZone]->GetCFL_AdaptParam(1);
            }

            CFLMin = config[val_iZone]->GetCFL_AdaptParam(2);
            CFLMax = config[val_iZone]->GetCFL_AdaptParam(3);

            CFLFactor = pow(Div, power);

            for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
            {
                CFL = config[val_iZone]->GetCFL(iMesh);
                CFL *= CFLFactor;

                if ((iMesh == MESH_0) && (CFL <= CFLMin))
                {
                    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
                    {
                        config[val_iZone]->SetCFL(iMesh, 1.001*CFLMin*MGFactor[iMesh]);
                    }
                    break;
                }
                if ((iMesh == MESH_0) && (CFL >= CFLMax))
                {
                    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
                        config[val_iZone]->SetCFL(iMesh, 0.999*CFLMax*MGFactor[iMesh]);
                    break;
                }

                config[val_iZone]->SetCFL(iMesh, CFL);
            }
            d_rhoResOld = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
        }

        void Cae_Output::SetForces_Breakdown(Geom_p **geometry, Solver_p ***solver_container, Config_p *config, INTE::INTE_Integration ***integration, unsigned short val_iZone)
        {

            char cstr[200];
            unsigned short iMarker_Monitoring;
            ofstream Breakdown_file;

            int rank = MASTER_NODE;
            bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
            bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
            bool freesurface = (config[val_iZone]->GetKind_Regime() == FREESURFACE);
            bool unsteady = (config[val_iZone]->GetUnsteady_Simulation() != NO);
            bool viscous = config[val_iZone]->GetViscous();
            bool grid_movement = config[val_iZone]->GetGrid_Movement();
            bool gravity = config[val_iZone]->GetGravityForce();
            bool turbulent = config[val_iZone]->GetKind_Solver() == RANS;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
            unsigned short nDim = geometry[val_iZone][FinestMesh]->GetnDim();
            bool flow = ((config[val_iZone]->GetKind_Solver() == EULER) || (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
                (config[val_iZone]->GetKind_Solver() == RANS));

            /*--- Output the mean flow solution using only the master node ---*/

            if ((rank == MASTER_NODE) && (flow)) 
            {
                cout << "Writing the forces breakdown file." << endl;
                /*--- Initialize variables to store information from all domains (direct solution) ---*/

                double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0,
                    Inv_CLift = 0.0, Inv_CDrag = 0.0, Inv_CSideForce = 0.0, Inv_CMx = 0.0, Inv_CMy = 0.0, Inv_CMz = 0.0, Inv_CEff = 0.0, Inv_CFx = 0.0, Inv_CFy = 0.0, Inv_CFz = 0.0,
                    *Surface_CLift = NULL, *Surface_CDrag = NULL, *Surface_CSideForce = NULL, *Surface_CEff = NULL, *Surface_CFx = NULL, *Surface_CFy = NULL, *Surface_CFz = NULL, *Surface_CMx = NULL, *Surface_CMy = NULL, *Surface_CMz = NULL,
                    *Surface_CLift_Inv = NULL, *Surface_CDrag_Inv = NULL, *Surface_CSideForce_Inv = NULL, *Surface_CEff_Inv = NULL, *Surface_CFx_Inv = NULL, *Surface_CFy_Inv = NULL, *Surface_CFz_Inv = NULL, *Surface_CMx_Inv = NULL, *Surface_CMy_Inv = NULL, *Surface_CMz_Inv = NULL;
                time_t now = time(0);
                string dt = ctime(&now); dt[24] = '.';

                /*--- Allocate memory for the coefficients being monitored ---*/

                Surface_CLift = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CDrag = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CSideForce = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CEff = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFx = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFy = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFz = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMx = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMy = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMz = new double[config[ZONE_0]->GetnMarker_Monitoring()];

                Surface_CLift_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CDrag_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CSideForce_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CEff_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFx_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFy_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CFz_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMx_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMy_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];
                Surface_CMz_Inv = new double[config[ZONE_0]->GetnMarker_Monitoring()];

                /*--- Flow solution coefficients ---*/

                Total_CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CLift();
                Total_CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CDrag();
                Total_CSideForce = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
                Total_CEff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CEff();
                Total_CMx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMx();
                Total_CMy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMy();
                Total_CMz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CMz();
                Total_CFx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFx();
                Total_CFy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFy();
                Total_CFz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetTotal_CFz();

                /*--- Flow inviscid solution coefficients ---*/

                Inv_CLift = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CLift_Inv();
                Inv_CDrag = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CDrag_Inv();
                Inv_CSideForce = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CSideForce_Inv();
                Inv_CEff = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CEff_Inv();
                Inv_CMx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CMx_Inv();
                Inv_CMy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CMy_Inv();
                Inv_CMz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CMz_Inv();
                Inv_CFx = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CFx_Inv();
                Inv_CFy = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CFy_Inv();
                Inv_CFz = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetAllBound_CFz_Inv();

                /*--- Look over the markers being monitored and get the desired values ---*/
                for (iMarker_Monitoring = 0; iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring(); iMarker_Monitoring++) 
                {
                    Surface_CLift[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
                    Surface_CDrag[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
                    Surface_CSideForce[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce(iMarker_Monitoring);
                    Surface_CEff[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CEff(iMarker_Monitoring);
                    Surface_CFx[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx(iMarker_Monitoring);
                    Surface_CFy[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy(iMarker_Monitoring);
                    Surface_CFz[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz(iMarker_Monitoring);
                    Surface_CMx[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
                    Surface_CMy[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
                    Surface_CMz[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);

                    Surface_CLift_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CLift_Inv(iMarker_Monitoring);
                    Surface_CDrag_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CDrag_Inv(iMarker_Monitoring);
                    Surface_CSideForce_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CSideForce_Inv(iMarker_Monitoring);
                    Surface_CEff_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CEff_Inv(iMarker_Monitoring);
                    Surface_CFx_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFx_Inv(iMarker_Monitoring);
                    Surface_CFy_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFy_Inv(iMarker_Monitoring);
                    Surface_CFz_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CFz_Inv(iMarker_Monitoring);
                    Surface_CMx_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMx_Inv(iMarker_Monitoring);
                    Surface_CMy_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMy_Inv(iMarker_Monitoring);
                    Surface_CMz_Inv[iMarker_Monitoring] = solver_container[val_iZone][FinestMesh][FLOW_SOL]->GetSurface_CMz_Inv(iMarker_Monitoring);
                }

                /*--- Write file name with extension ---*/
                string filename = config[val_iZone]->GetBreakdown_FileName();
                strcpy(cstr, filename.data());

                Breakdown_file.open(cstr, ios::out);

                Breakdown_file << endl << "-------------------------------------------------------------------------" << endl;
                Breakdown_file << "|   Local date and time: " << dt << "                      |" << endl;
                Breakdown_file << "-------------------------------------------------------------------------" << endl;

                Breakdown_file.precision(6);

                Breakdown_file << endl << endl << "Problem definition:" << endl << endl;

                if (compressible) 
                {
                    if (viscous) 
                    {
                        Breakdown_file << "Viscous flow: Computing pressure using the ideal gas law" << endl;
                        Breakdown_file << "based on the free-stream temperature and a density computed" << endl;
                        Breakdown_file << "from the Reynolds number." << endl;
                    }
                    else 
                    {
                        Breakdown_file << "Inviscid flow: Computing density based on free-stream" << endl;
                        Breakdown_file << "temperature and pressure using the ideal gas law." << endl;
                    }
                }

                if (grid_movement) 
                    Breakdown_file << "Force coefficients computed using MACH_MOTION." << endl;
                else
                    Breakdown_file << "Force coefficients computed using free-stream values." << endl;

                if (incompressible || freesurface) 
                {
                    Breakdown_file << "Viscous and Inviscid flow: rho_ref, and vel_ref" << endl;
                    Breakdown_file << "are based on the free-stream values, p_ref = rho_ref*vel_ref^2." << endl;
                    Breakdown_file << "The free-stream value of the pressure is 0." << endl;
                    Breakdown_file << "Mach number: " << config[val_iZone]->GetMach() << ", computed using the Bulk modulus." << endl;
                    Breakdown_file << "Angle of attack (deg): " << config[val_iZone]->GetAoA() << ", computed using the the free-stream velocity." << endl;
                    Breakdown_file << "Side slip angle (deg): " << config[val_iZone]->GetAoS() << ", computed using the the free-stream velocity." << endl;
                    if (viscous) 
                        Breakdown_file << "Reynolds number: " << config[val_iZone]->GetReynolds() << ", computed using free-stream values." << endl;
                    Breakdown_file << "Only dimensional computation, the grid should be dimensional." << endl;
                }

                Breakdown_file << "-- Input conditions:" << endl;

                if (compressible) 
                {
                    switch (config[val_iZone]->GetKind_FluidModel()) 
                    {
                    case STANDARD_AIR:
                        Breakdown_file << "Fluid Model: STANDARD_AIR " << endl;
                        Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant();
                        if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            Breakdown_file << " N.m/kg.K." << endl;
                        else if (config[val_iZone]->GetSystemMeasurements() == US)
                            Breakdown_file << " lbf.ft/slug.R." << endl;
                        Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
                        Breakdown_file << "Specific Heat Ratio: 1.4000 " << endl;
                        break;
                    case IDEAL_GAS:
                        Breakdown_file << "Fluid Model: IDEAL_GAS " << endl;
                        Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
                        Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
                        Breakdown_file << "Specific Heat Ratio: " << config[val_iZone]->GetGamma() << endl;
                        break;
                    case VW_GAS:
                        Breakdown_file << "Fluid Model: Van der Waals " << endl;
                        Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
                        Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
                        Breakdown_file << "Specific Heat Ratio: " << config[val_iZone]->GetGamma() << endl;
                        Breakdown_file << "Critical Pressure:   " << config[val_iZone]->GetPressure_Critical() << " Pa." << endl;
                        Breakdown_file << "Critical Temperature:  " << config[val_iZone]->GetTemperature_Critical() << " K." << endl;
                        Breakdown_file << "Critical Pressure (non-dim):   " << config[val_iZone]->GetPressure_Critical() / config[val_iZone]->GetPressure_Ref() << endl;
                        Breakdown_file << "Critical Temperature (non-dim) :  " << config[val_iZone]->GetTemperature_Critical() / config[val_iZone]->GetTemperature_Ref() << endl;
                        break;
                    case PR_GAS:
                        Breakdown_file << "Fluid Model: Peng-Robinson " << endl;
                        Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
                        Breakdown_file << "Specific gas constant(non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
                        Breakdown_file << "Specific Heat Ratio: " << config[val_iZone]->GetGamma() << endl;
                        Breakdown_file << "Critical Pressure:   " << config[val_iZone]->GetPressure_Critical() << " Pa." << endl;
                        Breakdown_file << "Critical Temperature:  " << config[val_iZone]->GetTemperature_Critical() << " K." << endl;
                        Breakdown_file << "Critical Pressure (non-dim):   " << config[val_iZone]->GetPressure_Critical() / config[val_iZone]->GetPressure_Ref() << endl;
                        Breakdown_file << "Critical Temperature (non-dim) :  " << config[val_iZone]->GetTemperature_Critical() / config[val_iZone]->GetTemperature_Ref() << endl;
                        break;
                    }

                    if (viscous) 
                    {
                        switch (config[val_iZone]->GetKind_ViscosityModel()) 
                        {
                        case CONSTANT_VISCOSITY:
                            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  " << endl;
                            Breakdown_file << "Laminar Viscosity: " << config[val_iZone]->GetMu_ConstantND()*config[val_iZone]->GetViscosity_Ref();
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                                Breakdown_file << " N.s/m^2." << endl;
                            else if (config[val_iZone]->GetSystemMeasurements() == US) 
                                Breakdown_file << " lbf.s/ft^2." << endl;
                            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND() << endl;
                            break;
                        case SUTHERLAND:
                            Breakdown_file << "Viscosity Model: SUTHERLAND " << endl;
                            Breakdown_file << "Ref. Laminar Viscosity: " << config[val_iZone]->GetMu_RefND()*config[val_iZone]->GetViscosity_Ref();
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                                Breakdown_file << " N.s/m^2." << endl;
                            else if (config[val_iZone]->GetSystemMeasurements() == US) 
                                Breakdown_file << " lbf.s/ft^2." << endl;
                            Breakdown_file << "Ref. Temperature: " << config[val_iZone]->GetMu_Temperature_RefND()*config[val_iZone]->GetTemperature_Ref();
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                                Breakdown_file << " K." << endl;
                            else if (config[val_iZone]->GetSystemMeasurements() == US) 
                                Breakdown_file << " R." << endl;
                            Breakdown_file << "Sutherland Constant: " << config[val_iZone]->GetMu_SND()*config[val_iZone]->GetTemperature_Ref();
                            if (config[val_iZone]->GetSystemMeasurements() == SI) 
                                Breakdown_file << " K." << endl;
                            else if (config[val_iZone]->GetSystemMeasurements() == US) 
                                Breakdown_file << " R." << endl;
                            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND() << endl;
                            Breakdown_file << "Ref. Temperature (non-dim): " << config[val_iZone]->GetMu_Temperature_RefND() << endl;
                            Breakdown_file << "Sutherland constant (non-dim): " << config[val_iZone]->GetMu_SND() << endl;
                            break;
                        }
                        switch (config[val_iZone]->GetKind_ConductivityModel()) 
                        {
                        case CONSTANT_PRANDTL:
                            Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  " << endl;
                            Breakdown_file << "Prandtl: " << config[val_iZone]->GetPrandtl_Lam() << endl;
                            break;
                        case CONSTANT_CONDUCTIVITY:
                            Breakdown_file << "Conductivity Model: CONSTANT_CONDUCTIVITY " << endl;
                            Breakdown_file << "Molecular Conductivity: " << config[val_iZone]->GetKt_ConstantND()*config[val_iZone]->GetConductivity_Ref() << " W/m^2.K." << endl;
                            Breakdown_file << "Molecular Conductivity (non-dim): " << config[val_iZone]->GetKt_ConstantND() << endl;
                            break;
                        }
                    }
                }

                if (incompressible || freesurface) 
                {
                    Breakdown_file << "Bulk modulus: " << config[val_iZone]->GetBulk_Modulus();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " Pa." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " psf." << endl;
                    Breakdown_file << "Artificial compressibility factor: " << config[val_iZone]->GetArtComp_Factor();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " Pa." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US)
                        Breakdown_file << " psf." << endl;
                }

                Breakdown_file << "Free-stream static pressure: " << config[val_iZone]->GetPressure_FreeStream();
                if (config[val_iZone]->GetSystemMeasurements() == SI) 
                    Breakdown_file << " Pa." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US) 
                    Breakdown_file << " psf." << endl;

                Breakdown_file << "Free-stream total pressure: " << config[val_iZone]->GetPressure_FreeStream() * pow(1.0 + config[val_iZone]->GetMach()*config[val_iZone]->GetMach()*0.5*(config[val_iZone]->GetGamma() - 1.0), config[val_iZone]->GetGamma() / (config[val_iZone]->GetGamma() - 1.0));
                if (config[val_iZone]->GetSystemMeasurements() == SI) 
                    Breakdown_file << " Pa." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US)
                    Breakdown_file << " psf." << endl;

                if (compressible) 
                {
                    Breakdown_file << "Free-stream temperature: " << config[val_iZone]->GetTemperature_FreeStream();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " K." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " R." << endl;
                }

                Breakdown_file << "Free-stream density: " << config[val_iZone]->GetDensity_FreeStream();
                if (config[val_iZone]->GetSystemMeasurements() == SI) 
                    Breakdown_file << " kg/m^3." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US) 
                    Breakdown_file << " slug/ft^3." << endl;

                if (nDim == 2) 
                {
                    Breakdown_file << "Free-stream velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
                    Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ")";
                }
                if (nDim == 3) 
                {
                    Breakdown_file << "Free-stream velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
                    Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ", " << config[val_iZone]->GetVelocity_FreeStream()[2] << ")";
                }
                if (config[val_iZone]->GetSystemMeasurements() == SI) 
                    Breakdown_file << " m/s. ";
                else if (config[val_iZone]->GetSystemMeasurements() == US) 
                    Breakdown_file << " ft/s. ";

                Breakdown_file << "Magnitude: " << config[val_iZone]->GetModVel_FreeStream();
                if (config[val_iZone]->GetSystemMeasurements() == SI) 
                    Breakdown_file << " m/s." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US)
                    Breakdown_file << " ft/s." << endl;

                if (compressible) 
                {
                    Breakdown_file << "Free-stream total energy per unit mass: " << config[val_iZone]->GetEnergy_FreeStream();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " m^2/s^2." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " ft^2/s^2." << endl;
                }

                if (viscous)
                {
                    Breakdown_file << "Free-stream viscosity: " << config[val_iZone]->GetViscosity_FreeStream();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " N.s/m^2." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US)
                        Breakdown_file << " lbf.s/ft^2." << endl;
                    if (turbulent) 
                    {
                        Breakdown_file << "Free-stream turb. kinetic energy per unit mass: " << config[val_iZone]->GetTke_FreeStream();
                        if (config[val_iZone]->GetSystemMeasurements() == SI) 
                            Breakdown_file << " m^2/s^2." << endl;
                        else if (config[val_iZone]->GetSystemMeasurements() == US) 
                            Breakdown_file << " ft^2/s^2." << endl;
                        Breakdown_file << "Free-stream specific dissipation: " << config[val_iZone]->GetOmega_FreeStream();
                        if (config[val_iZone]->GetSystemMeasurements() == SI)
                            Breakdown_file << " 1/s." << endl;
                        else if (config[val_iZone]->GetSystemMeasurements() == US) 
                            Breakdown_file << " 1/s." << endl;
                    }
                }

                if (unsteady) 
                { 
                    Breakdown_file << "Total time: " << config[val_iZone]->GetTotal_UnstTime() << " s. Time step: " << config[val_iZone]->GetDelta_UnstTime() << " s." << endl; 
                }

                /*--- Print out reference values. ---*/
                Breakdown_file << "-- Reference values:" << endl;

                if (compressible) 
                {
                    Breakdown_file << "Reference specific gas constant: " << config[val_iZone]->GetGas_Constant_Ref();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " N.m/kg.K." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " lbf.ft/slug.R." << endl;
                }

                Breakdown_file << "Reference pressure: " << config[val_iZone]->GetPressure_Ref();
                if (config[val_iZone]->GetSystemMeasurements() == SI)
                    Breakdown_file << " Pa." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US) 
                    Breakdown_file << " psf." << endl;

                if (compressible) 
                {
                    Breakdown_file << "Reference temperature: " << config[val_iZone]->GetTemperature_Ref();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " K." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " R." << endl;
                }

                Breakdown_file << "Reference density: " << config[val_iZone]->GetDensity_Ref();
                if (config[val_iZone]->GetSystemMeasurements() == SI)
                    Breakdown_file << " kg/m^3." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US)
                    Breakdown_file << " slug/ft^3." << endl;

                Breakdown_file << "Reference velocity: " << config[val_iZone]->GetVelocity_Ref();
                if (config[val_iZone]->GetSystemMeasurements() == SI) 
                    Breakdown_file << " m/s." << endl;
                else if (config[val_iZone]->GetSystemMeasurements() == US) 
                    Breakdown_file << " ft/s." << endl;

                if (compressible) 
                {
                    Breakdown_file << "Reference energy per unit mass: " << config[val_iZone]->GetEnergy_Ref();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " m^2/s^2." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " ft^2/s^2." << endl;
                }

                if (incompressible || freesurface) 
                {
                    Breakdown_file << "Reference length: " << config[val_iZone]->GetLength_Ref();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " m." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " in." << endl;
                }

                if (viscous) 
                {
                    Breakdown_file << "Reference viscosity: " << config[val_iZone]->GetViscosity_Ref();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " N.s/m^2." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US)
                        Breakdown_file << " lbf.s/ft^2." << endl;
                    Breakdown_file << "Reference conductivity: " << config[val_iZone]->GetConductivity_Ref();
                    if (config[val_iZone]->GetSystemMeasurements() == SI)
                        Breakdown_file << " W/m^2.K." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US) 
                        Breakdown_file << " lbf/ft.s.R." << endl;
                }


                if (unsteady) 
                    Breakdown_file << "Reference time: " << config[val_iZone]->GetTime_Ref() << " s." << endl;

                /*--- Print out resulting non-dim values here. ---*/
                Breakdown_file << "-- Resulting non-dimensional state:" << endl;
                Breakdown_file << "Mach number (non-dim): " << config[val_iZone]->GetMach() << endl;
                if (viscous) 
                {
                    Breakdown_file << "Reynolds number (non-dim): " << config[val_iZone]->GetReynolds() << ". Re length: " << config[val_iZone]->GetLength_Reynolds();
                    if (config[val_iZone]->GetSystemMeasurements() == SI) 
                        Breakdown_file << " m." << endl;
                    else if (config[val_iZone]->GetSystemMeasurements() == US)
                        Breakdown_file << " ft." << endl;
                }

                if (gravity) 
                {
                    Breakdown_file << "Froude number (non-dim): " << config[val_iZone]->GetFroude() << endl;
                    Breakdown_file << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*config[val_iZone]->GetFroude()*config[val_iZone]->GetFroude() << endl;
                }

                if (compressible) 
                {
                    Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
                    Breakdown_file << "Free-stream temperature (non-dim): " << config[val_iZone]->GetTemperature_FreeStreamND() << endl;
                }

                Breakdown_file << "Free-stream pressure (non-dim): " << config[val_iZone]->GetPressure_FreeStreamND() << endl;

                Breakdown_file << "Free-stream density (non-dim): " << config[val_iZone]->GetDensity_FreeStreamND() << endl;

                if (nDim == 2) 
                {
                    Breakdown_file << "Free-stream velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
                    Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << "). ";
                }
                else 
                {
                    Breakdown_file << "Free-stream velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
                    Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << ", " << config[val_iZone]->GetVelocity_FreeStreamND()[2] << "). ";
                }
                Breakdown_file << "Magnitude: " << config[val_iZone]->GetModVel_FreeStreamND() << endl;

                if (compressible)
                    Breakdown_file << "Free-stream total energy per unit mass (non-dim): " << config[val_iZone]->GetEnergy_FreeStreamND() << endl;

                if (viscous)
                {
                    Breakdown_file << "Free-stream viscosity (non-dim): " << config[val_iZone]->GetViscosity_FreeStreamND() << endl;
                    if (turbulent) 
                    {
                        Breakdown_file << "Free-stream turb. kinetic energy (non-dim): " << config[val_iZone]->GetTke_FreeStreamND() << endl;
                        Breakdown_file << "Free-stream specific dissipation (non-dim): " << config[val_iZone]->GetOmega_FreeStreamND() << endl;
                    }
                }

                if (unsteady) 
                {
                    Breakdown_file << "Total time (non-dim): " << config[val_iZone]->GetTotal_UnstTimeND() << endl;
                    Breakdown_file << "Time step (non-dim): " << config[val_iZone]->GetDelta_UnstTimeND() << endl;
                }

                Breakdown_file << endl << endl << "Forces breakdown:" << endl << endl;

                Breakdown_file << "Total CL:    ";
                Breakdown_file.width(11); Breakdown_file << Total_CLift;
                Breakdown_file << " | Pressure Component (";
                Breakdown_file.width(5); Breakdown_file << int((Inv_CLift*100.0) / (Total_CLift + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Inv_CLift;
                Breakdown_file << " | Friction Component (";
                Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CLift*100.0) / (Total_CLift + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Total_CLift - Inv_CLift << endl;

                Breakdown_file << "Total CD:    ";
                Breakdown_file.width(11); Breakdown_file << Total_CDrag;
                Breakdown_file << " | Pressure Component (";
                Breakdown_file.width(5); Breakdown_file << int((Inv_CDrag*100.0) / (Total_CDrag + EPS)) << "%): ";;
                Breakdown_file.width(11); Breakdown_file << Inv_CDrag;
                Breakdown_file << " | Friction Component (";
                Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CDrag*100.0) / (Total_CDrag + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Total_CDrag - Inv_CDrag << endl;

                if (nDim == 3) 
                {
                    Breakdown_file << "Total CSF:   ";
                    Breakdown_file.width(11); Breakdown_file << Total_CSideForce;
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Inv_CSideForce*100.0) / (Total_CSideForce + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Inv_CSideForce;
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CSideForce*100.0) / (Total_CSideForce + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Total_CSideForce - Inv_CSideForce << endl;
                }

                Breakdown_file << "Total CL/CD: ";
                Breakdown_file.width(11); Breakdown_file << Total_CEff;
                Breakdown_file << " | Pressure Component (";
                Breakdown_file.width(5); Breakdown_file << int((Inv_CEff*100.0) / (Total_CEff + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Inv_CEff;
                Breakdown_file << " | Friction Component (";
                Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CEff*100.0) / (Total_CEff + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Total_CEff - Inv_CEff << endl;

                if (nDim == 3) 
                {
                    Breakdown_file << "Total CMx:   ";
                    Breakdown_file.width(11); Breakdown_file << Total_CMx;
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Inv_CMx*100.0) / (Total_CMx + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Inv_CMx;
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CMx*100.0) / (Total_CMx + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Total_CMx - Inv_CMx << endl;

                    Breakdown_file << "Total CMy:   ";
                    Breakdown_file.width(11); Breakdown_file << Total_CMy;
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Inv_CMy*100.0) / (Total_CMy + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Inv_CMy;
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CMy*100.0) / (Total_CMy + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Total_CMy - Inv_CMy << endl;
                }

                Breakdown_file << "Total CMz:   ";
                Breakdown_file.width(11); Breakdown_file << Total_CMz;
                Breakdown_file << " | Pressure Component (";
                Breakdown_file.width(5); Breakdown_file << int((Inv_CMz*100.0) / (Total_CMz + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Inv_CMz;
                Breakdown_file << " | Friction Component (";
                Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CMz*100.0) / (Total_CMz + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Total_CMz - Inv_CMz << endl;

                Breakdown_file << "Total CFx:   ";
                Breakdown_file.width(11); Breakdown_file << Total_CFx;
                Breakdown_file << " | Pressure Component (";
                Breakdown_file.width(5); Breakdown_file << int((Inv_CFx*100.0) / (Total_CFx + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Inv_CFx;
                Breakdown_file << " | Friction Component (";
                Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CFx*100.0) / (Total_CFx + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Total_CFx - Inv_CFx << endl;

                Breakdown_file << "Total CFy:   ";
                Breakdown_file.width(11); Breakdown_file << Total_CFy;
                Breakdown_file << " | Pressure Component (";
                Breakdown_file.width(5); Breakdown_file << int((Inv_CFy*100.0) / (Total_CFy + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Inv_CFy;
                Breakdown_file << " | Friction Component (";
                Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CFy*100.0) / (Total_CFy + EPS));
                Breakdown_file << "%): ";
                Breakdown_file.width(11); Breakdown_file << Total_CFy - Inv_CFy << endl;

                if (nDim == 3) 
                {
                    Breakdown_file << "Total CFz:   ";
                    Breakdown_file.width(11); Breakdown_file << Total_CFz;
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Inv_CFz*100.0) / (Total_CFz + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Inv_CFz;
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Inv_CFz*100.0) / (Total_CFz + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Total_CFz - Inv_CFz << endl;
                }

                Breakdown_file << endl << endl;

                for (iMarker_Monitoring = 0; iMarker_Monitoring < config[val_iZone]->GetnMarker_Monitoring(); iMarker_Monitoring++)
                {
                    Breakdown_file << "Surface name: " << config[val_iZone]->GetMarker_Monitoring(iMarker_Monitoring) << endl << endl;

                    Breakdown_file << "Total CL    (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CLift[iMarker_Monitoring] * 100.0) / (Total_CLift + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CLift[iMarker_Monitoring];
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CLift_Inv[iMarker_Monitoring] * 100.0) / (Surface_CLift[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CLift_Inv[iMarker_Monitoring];
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CLift_Inv[iMarker_Monitoring] * 100.0) / (Surface_CLift[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CLift[iMarker_Monitoring] - Surface_CLift_Inv[iMarker_Monitoring] << endl;

                    Breakdown_file << "Total CD    (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CDrag[iMarker_Monitoring] * 100.0) / (Total_CDrag + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CDrag[iMarker_Monitoring];
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CDrag_Inv[iMarker_Monitoring] * 100.0) / (Surface_CDrag[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CDrag_Inv[iMarker_Monitoring];
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CDrag_Inv[iMarker_Monitoring] * 100.0) / (Surface_CDrag[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CDrag[iMarker_Monitoring] - Surface_CDrag_Inv[iMarker_Monitoring] << endl;

                    if (nDim == 3) 
                    {
                        Breakdown_file << "Total CSF   (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CSideForce[iMarker_Monitoring] * 100.0) / (Total_CSideForce + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CSideForce[iMarker_Monitoring];
                        Breakdown_file << " | Pressure Component (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CSideForce_Inv[iMarker_Monitoring] * 100.0) / (Surface_CSideForce[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CSideForce_Inv[iMarker_Monitoring];
                        Breakdown_file << " | Friction Component (";
                        Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CSideForce_Inv[iMarker_Monitoring] * 100.0) / (Surface_CSideForce[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CSideForce[iMarker_Monitoring] - Surface_CSideForce_Inv[iMarker_Monitoring] << endl;
                    }

                    Breakdown_file << "Total CL/CD (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CEff[iMarker_Monitoring] * 100.0) / (Total_CEff + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CEff[iMarker_Monitoring];
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CEff_Inv[iMarker_Monitoring] * 100.0) / (Surface_CEff[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CEff_Inv[iMarker_Monitoring];
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CEff_Inv[iMarker_Monitoring] * 100.0) / (Surface_CEff[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";

                    Breakdown_file.width(11); Breakdown_file << Surface_CEff[iMarker_Monitoring] - Surface_CEff_Inv[iMarker_Monitoring] << endl;

                    if (nDim == 3) 
                    {
                        Breakdown_file << "Total CMx   (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CMx[iMarker_Monitoring] * 100.0) / (Total_CMx + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CMx[iMarker_Monitoring];
                        Breakdown_file << " | Pressure Component (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CMx_Inv[iMarker_Monitoring] * 100.0) / (Surface_CMx[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CMx_Inv[iMarker_Monitoring];
                        Breakdown_file << " | Friction Component (";
                        Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CMx_Inv[iMarker_Monitoring] * 100.0) / (Surface_CMx[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CMx[iMarker_Monitoring] - Surface_CMx_Inv[iMarker_Monitoring] << endl;

                        Breakdown_file << "Total CMy   (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CMy[iMarker_Monitoring] * 100.0) / (Total_CMy + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CMy[iMarker_Monitoring];
                        Breakdown_file << " | Pressure Component (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CMy_Inv[iMarker_Monitoring] * 100.0) / (Surface_CMy[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CMy_Inv[iMarker_Monitoring];
                        Breakdown_file << " | Friction Component (";
                        Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CMy_Inv[iMarker_Monitoring] * 100.0) / (Surface_CMy[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CMy[iMarker_Monitoring] - Surface_CMy_Inv[iMarker_Monitoring] << endl;
                    }

                    Breakdown_file << "Total CMz   (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CMz[iMarker_Monitoring] * 100.0) / (Total_CMz + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CMz[iMarker_Monitoring];
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CMz_Inv[iMarker_Monitoring] * 100.0) / (Surface_CMz[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CMz_Inv[iMarker_Monitoring];
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CMz_Inv[iMarker_Monitoring] * 100.0) / (Surface_CMz[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CMz[iMarker_Monitoring] - Surface_CMz_Inv[iMarker_Monitoring] << endl;

                    Breakdown_file << "Total CFx   (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CFx[iMarker_Monitoring] * 100.0) / (Total_CFx + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CFx[iMarker_Monitoring];
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CFx_Inv[iMarker_Monitoring] * 100.0) / (Surface_CFx[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CFx_Inv[iMarker_Monitoring];
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CFx_Inv[iMarker_Monitoring] * 100.0) / (Surface_CFx[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CFx[iMarker_Monitoring] - Surface_CFx_Inv[iMarker_Monitoring] << endl;

                    Breakdown_file << "Total CFy   (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CFy[iMarker_Monitoring] * 100.0) / (Total_CFy + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CFy[iMarker_Monitoring];
                    Breakdown_file << " | Pressure Component (";
                    Breakdown_file.width(5); Breakdown_file << int((Surface_CFy_Inv[iMarker_Monitoring] * 100.0) / (Surface_CFy[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CFy_Inv[iMarker_Monitoring];
                    Breakdown_file << " | Friction Component (";
                    Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CFy_Inv[iMarker_Monitoring] * 100.0) / (Surface_CFy[iMarker_Monitoring] + EPS));
                    Breakdown_file << "%): ";
                    Breakdown_file.width(11); Breakdown_file << Surface_CFy[iMarker_Monitoring] - Surface_CFy_Inv[iMarker_Monitoring] << endl;

                    if (nDim == 3) 
                    {
                        Breakdown_file << "Total CFz   (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CFz[iMarker_Monitoring] * 100.0) / (Total_CFz + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CFz[iMarker_Monitoring];
                        Breakdown_file << " | Pressure Component (";
                        Breakdown_file.width(5); Breakdown_file << int((Surface_CFz_Inv[iMarker_Monitoring] * 100.0) / (Surface_CFz[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CFz_Inv[iMarker_Monitoring];
                        Breakdown_file << " | Friction Component (";
                        Breakdown_file.width(5); Breakdown_file << int(100.0 - (Surface_CFz_Inv[iMarker_Monitoring] * 100.0) / (Surface_CFz[iMarker_Monitoring] + EPS));
                        Breakdown_file << "%): ";
                        Breakdown_file.width(11); Breakdown_file << Surface_CFz[iMarker_Monitoring] - Surface_CFz_Inv[iMarker_Monitoring] << endl;
                    }

                    Breakdown_file << endl;
                }

                delete[] Surface_CLift;
                delete[] Surface_CDrag;
                delete[] Surface_CSideForce;
                delete[] Surface_CEff;
                delete[] Surface_CFx;
                delete[] Surface_CFy;
                delete[] Surface_CFz;
                delete[] Surface_CMx;
                delete[] Surface_CMy;
                delete[] Surface_CMz;

                delete[] Surface_CLift_Inv;
                delete[] Surface_CDrag_Inv;
                delete[] Surface_CSideForce_Inv;
                delete[] Surface_CEff_Inv;
                delete[] Surface_CFx_Inv;
                delete[] Surface_CFy_Inv;
                delete[] Surface_CFz_Inv;
                delete[] Surface_CMx_Inv;
                delete[] Surface_CMy_Inv;
                delete[] Surface_CMz_Inv;

                Breakdown_file.close();
            }
        }



        void Cae_Output::SetBaselineResult_Files(SOLV::SOLV_Solver **solver, GEOM::GEOM_Geometry **geometry, Config_p *config,
            unsigned long iExtIter, unsigned short val_nZone) {

            int rank = MASTER_NODE;
            int size = SINGLE_NODE;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

            unsigned short iZone;

            for (iZone = 0; iZone < val_nZone; iZone++) {

                /*--- Flags identifying the types of files to be written. ---*/

                bool Low_MemoryOutput = config[iZone]->GetLow_MemoryOutput();
                bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
                bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

                /*--- Get the file output format ---*/

                unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();

                /*--- Merge the node coordinates and connectivity if necessary. This
                is only performed if a volume solution file is requested, and it
                is active by default. ---*/

                if ((Wrt_Vol || Wrt_Srf) && (!Low_MemoryOutput)) {
                    if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
                    MergeConnectivity(config[iZone], geometry[iZone], iZone);
                }

                /*--- Merge the solution data needed for volume solutions and restarts ---*/

                if ((Wrt_Vol || Wrt_Srf) && (!Low_MemoryOutput)) {
                    if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
                    MergeBaselineSolution(config[iZone], geometry[iZone], solver[iZone], iZone);
                }

                /*--- Write restart, Tecplot or Paraview files using the merged data.
                This data lives only on the master, and these routines are currently
                executed by the master proc alone (as if in serial). ---*/

                if (!Low_MemoryOutput) {

                    if (rank == MASTER_NODE) {

                        if (Wrt_Vol) {

                            switch (FileFormat) {

                            case TECPLOT:

                                /*--- Write a Tecplot ASCII file ---*/

                                if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
                                SetTecplotASCII(config[iZone], geometry[iZone], solver, iZone, val_nZone, false);
                                DeallocateConnectivity(config[iZone], geometry[iZone], false);
                                break;

                            case FIELDVIEW:

                                /*--- Write a FieldView ASCII file ---*/

                                if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
                                SetFieldViewASCII(config[iZone], geometry[iZone], iZone, val_nZone);
                                DeallocateConnectivity(config[iZone], geometry[iZone], false);
                                break;

                            case TECPLOT_BINARY:

                                /*--- Write a Tecplot binary solution file ---*/

                                if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (volume grid)." << endl;
                                SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone], iZone);
                                SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone], iZone);
                                break;

                            case FIELDVIEW_BINARY:

                                /*--- Write a binary binary file ---*/

                                if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
                                SetFieldViewBinary(config[iZone], geometry[iZone], iZone, val_nZone);
                                DeallocateConnectivity(config[iZone], geometry[iZone], false);
                                break;

                            case PARAVIEW:

                                /*--- Write a Paraview ASCII file ---*/

                                if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (volume grid)." << endl;
                                SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, false);
                                DeallocateConnectivity(config[iZone], geometry[iZone], false);
                                break;

                            default:
                                break;
                            }

                        }

                        if (Wrt_Srf) {

                            switch (FileFormat) {

                            case TECPLOT:

                                /*--- Write a Tecplot ASCII file ---*/

                                if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
                                SetTecplotASCII(config[iZone], geometry[iZone], solver, iZone, val_nZone, true);
                                DeallocateConnectivity(config[iZone], geometry[iZone], true);
                                break;

                            case TECPLOT_BINARY:

                                /*--- Write a Tecplot binary solution file ---*/

                                if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (surface grid)." << endl;
                                SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone], iZone);
                                SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone], iZone);
                                break;

                            case PARAVIEW:

                                /*--- Write a Paraview ASCII file ---*/

                                if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (surface grid)." << endl;
                                SetParaview_ASCII(config[iZone], geometry[iZone], iZone, val_nZone, true);
                                DeallocateConnectivity(config[iZone], geometry[iZone], true);
                                break;

                            default:
                                break;
                            }
                        }

                        if (FileFormat == TECPLOT_BINARY) {
                            if (!d_isBaseOutput)
                                DeallocateConnectivity(config[iZone], geometry[iZone], false);
                            if (!d_isSurfOutput)
                                DeallocateConnectivity(config[iZone], geometry[iZone], d_isSurfOutput);
                        }

                        if (Wrt_Vol || Wrt_Srf)
                            DeallocateSolution(config[iZone], geometry[iZone]);
                    }

                }

                else {

                    if (Wrt_Vol) {

                        if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
                        char buffer_char[50], out_file[MAX_STRING_SIZE];

                        string filename;
                        if (!config[iZone]->GetAdjoint()) filename = config[iZone]->GetFlow_FileName();
                        else filename = config[iZone]->GetAdj_FileName();

                        if (size > 1) {
                            sprintf(buffer_char, "_%d", int(rank + 1));
                            filename = filename + buffer_char;
                        }

                        sprintf(buffer_char, ".dat");
                        strcpy(out_file, filename.c_str()); strcat(out_file, buffer_char);
                        SetTecplotASCII_LowMemory(config[iZone], geometry[iZone], solver, out_file, false);
                    }

                    if (Wrt_Srf) {

                        if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
                        char buffer_char[50], out_file[MAX_STRING_SIZE];

                        string filename;
                        if (!config[iZone]->GetAdjoint()) filename = config[iZone]->GetSurfFlowCoeff_FileName();
                        else filename = config[iZone]->GetSurfAdjCoeff_FileName();

                        if (size > 1) {
                            sprintf(buffer_char, "_%d", int(rank + 1));
                            filename = filename + buffer_char;
                        }

                        sprintf(buffer_char, ".dat");
                        strcpy(out_file, filename.c_str()); strcat(out_file, buffer_char);
                        SetTecplotASCII_LowMemory(config[iZone], geometry[iZone], solver, out_file, true);
                    }

                }

                /*--- Final broadcast (informing other procs that the base output
                file was written). ---*/

#ifdef HAVE_MPI
                MPI_Bcast(&d_isBaseOutput, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

            }
        }

        void Cae_Output::SetMesh_Files(GEOM::GEOM_Geometry **geometry, Config_p *config, unsigned short val_nZone, bool new_file, bool su2_file) {

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            unsigned short iZone;

            for (iZone = 0; iZone < val_nZone; iZone++) {

                /*--- Flags identifying the types of files to be written. ---*/

                bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol() && config[iZone]->GetVisualize_Deformation();
                bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol() && config[iZone]->GetVisualize_Deformation();;

                /*--- Merge the node coordinates and connectivity if necessary. This
                is only performed if a volume solution file is requested, and it
                is active by default. ---*/

                if (rank == MASTER_NODE) cout << "Merging grid connectivity." << endl;
                MergeConnectivity(config[iZone], geometry[iZone], iZone);

                /*--- Merge coordinates of all grid nodes (excluding ghost points).
                The grid coordinates are always merged and included first in the
                restart files. ---*/

                if (rank == MASTER_NODE) cout << "Merging grid coordinates." << endl;
                MergeCoordinates(config[iZone], geometry[iZone]);

                /*--- Write restart, Tecplot or Paraview files using the merged data.
                This data lives only on the master, and these routines are currently
                executed by the master proc alone (as if in serial). ---*/

                if (rank == MASTER_NODE) {

                    if (Wrt_Vol) {

                        if (rank == MASTER_NODE) cout << "Writing volume mesh file." << endl;

                        /*--- Write a Tecplot ASCII file ---*/
                        if (config[iZone]->GetOutput_FileFormat() == PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone, val_nZone, false, new_file);
                        else SetTecplotASCII_Mesh(config[iZone], geometry[iZone], false, new_file);

                    }

                    if (Wrt_Srf) {

                        if (rank == MASTER_NODE) cout << "Writing surface mesh file." << endl;

                        /*--- Write a Tecplot ASCII file ---*/
                        if (config[iZone]->GetOutput_FileFormat() == PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone, val_nZone, true, new_file);
                        else SetTecplotASCII_Mesh(config[iZone], geometry[iZone], true, new_file);


                    }

                    if (rank == MASTER_NODE) cout << "Writing .su2 file." << endl;

                    /*--- Write a .su2 ASCII file ---*/

                    if (su2_file) SetSU2_MeshASCII(config[iZone], geometry[iZone]);

                    /*--- Deallocate connectivity ---*/

                    DeallocateConnectivity(config[iZone], geometry[iZone], true);

                }

                /*--- Final broadcast (informing other procs that the base output
                file was written). ---*/

#ifdef HAVE_MPI
                MPI_Bcast(&d_isBaseOutput, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

            }
        }

        void Cae_Output::SetMassFlowRate(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, Config_p config) {
            unsigned short iDim, iMarker_monitor, iMarker;
            unsigned long iVertex, iPoint;
            double Vector[3], Total_Mdot = 0.0;
            unsigned short nDim = geometry->GetnDim();

            for (iMarker = 0; iMarker < config->GetnMarker_Monitoring(); iMarker++) {
                iMarker_monitor = config->GetMarker_All_Monitoring(iMarker);

                for (iVertex = 0; iVertex < geometry->nVertex[iMarker_monitor]; iVertex++) {
                    iPoint = geometry->vertex[iMarker_monitor][iVertex]->GetNode();

                    if (geometry->node[iPoint]->GetDomain()) {
                        geometry->vertex[iMarker_monitor][iVertex]->GetNormal(Vector);
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Total_Mdot += Vector[iDim] * (solver_container->node[iPoint]->GetSolution(iDim + 1));
                        }
                    }
                }
            }

#ifdef HAVE_MPI
            /*--- Add AllBound information using all the nodes ---*/
            double My_Total_Mdot = Total_Mdot;    Total_Mdot = 0.0;
            MPI_Allreduce(&My_Total_Mdot, &Total_Mdot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
            /*--- Set the output: reusing same variable from OneDimensionalOutput code ---*/
            solver_container->SetOneD_MassFlowRate(Total_Mdot);
        }

        void Cae_Output::OneDimensionalOutput(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, Config_p config) {

            unsigned long iVertex, iPoint;
            unsigned short iDim, iMarker, Out1D;
            double *Normal = NULL, Area = 0.0, OverArea = 0.0, UnitaryNormal[3],
                Stag_Pressure, Mach, Temperature, Pressure = 0.0, Density = 0.0, Velocity2, Enthalpy, RhoU, U,// local values at each node (Velocity2 = V^2). U = normal velocity
                SumPressure = 0.0, SumStagPressure = 0.0, SumArea = 0.0, SumMach = 0.0, SumTemperature = 0.0, SumForUref = 0.0, SumRhoU = 0.0, SumEnthalpy = 0.0,// sum of (local value ) * (dA) (integral)
                AveragePressure = 0.0, AverageMach = 0.0, AverageTemperature = 0.0, MassFlowRate = 0.0, // Area Averaged value ( sum / A )
                VelocityRef = 0.0, EnthalpyRef = 0.0, DensityRef = 0.0, PressureRef = 0.0; // Flux conserved values. TemperatureRef follows ideal gas

            bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
            bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
            bool freesurface = (config->GetKind_Regime() == FREESURFACE);
            double Gamma = config->GetGamma();
            unsigned short nDim = geometry->GetnDim();


            /*--- Loop over the markers ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

                Out1D = config->GetMarker_All_Out_1D(iMarker);

                /*--- Loop over the vertices to compute the output ---*/


                if (Out1D == YES) {

                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

                        /*--- Find the normal direction ---*/

                        if (geometry->node[iPoint]->GetDomain()) {


                            /*--- Compute area, and unitary normal ---*/
                            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim] * Normal[iDim]; Area = sqrt(Area);
                            for (iDim = 0; iDim < nDim; iDim++) UnitaryNormal[iDim] = -Normal[iDim] / Area;

                            if (compressible) {
                                Pressure = solver_container->node[iPoint]->GetPressure();
                                Density = solver_container->node[iPoint]->GetDensity();
                            }
                            if (incompressible || freesurface) {
                                Pressure = solver_container->node[iPoint]->GetPressureInc();
                                Density = solver_container->node[iPoint]->GetDensityInc();
                            }

                            /*-- Find velocity normal to the marked surface/opening --*/

                            U = 0.0;
                            for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                                U += UnitaryNormal[iDim] * solver_container->node[iPoint]->GetVelocity(iDim);
                            }

                            Enthalpy = solver_container->node[iPoint]->GetEnthalpy();
                            Velocity2 = solver_container->node[iPoint]->GetVelocity2();
                            Temperature = solver_container->node[iPoint]->GetTemperature();

                            Mach = (sqrt(Velocity2)) / solver_container->node[iPoint]->GetSoundSpeed();
                            Stag_Pressure = Pressure*pow((1.0 + ((Gamma - 1.0) / 2.0)*pow(Mach, 2.0)), (Gamma / (Gamma - 1.0)));

                            RhoU = U*Density;
                            SumStagPressure += Stag_Pressure * Area;
                            SumArea += Area;
                            SumMach += Mach*Area;
                            SumPressure += Pressure * Area;
                            SumTemperature += Temperature*Area;
                            SumRhoU += RhoU*Area;
                            SumForUref += RhoU*U*U*Area;
                            SumEnthalpy += RhoU*Enthalpy*Area;

                        }
                    }

                    if (SumRhoU != 0.0) { // To avoid division by 0

                        OverArea = 1.0 / SumArea;
                        AveragePressure += abs(SumStagPressure*OverArea);
                        AverageMach += abs(SumMach*OverArea);
                        AverageTemperature += abs(SumTemperature*OverArea);
                        MassFlowRate += SumRhoU;
                        PressureRef += abs(SumPressure*OverArea);
                        VelocityRef += abs(sqrt(abs(SumForUref / SumRhoU)));
                        EnthalpyRef += abs(SumEnthalpy / SumRhoU);
                        DensityRef += abs(PressureRef*Gamma / (Gamma - 1) / (EnthalpyRef - 0.5*VelocityRef*VelocityRef));

                    }

                }

            }

#ifdef HAVE_MPI

            /*--- Add AllBound information using all the nodes ---*/

            double My_AveragePressure = AveragePressure;    AveragePressure = 0.0;
            double My_AverageMach = AverageMach;        AverageMach = 0.0;
            double My_AverageTemperature = AverageTemperature; AverageTemperature = 0.0;
            double My_MassFlowRate = MassFlowRate;       MassFlowRate = 0.0;
            double My_PressureRef = PressureRef;        PressureRef = 0.0;
            double My_VelocityRef = VelocityRef;        VelocityRef = 0.0;
            double My_EnthalpyRef = EnthalpyRef;        EnthalpyRef = 0.0;
            double My_DensityRef = DensityRef;         DensityRef = 0.0;

            MPI_Allreduce(&My_AveragePressure, &AveragePressure, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_AverageMach, &AverageMach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_AverageTemperature, &AverageTemperature, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_MassFlowRate, &MassFlowRate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_PressureRef, &PressureRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_VelocityRef, &VelocityRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_EnthalpyRef, &EnthalpyRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&My_DensityRef, &DensityRef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#endif

            /*--- Set the 1D output ---*/

            solver_container->SetOneD_TotalPress(AveragePressure);
            solver_container->SetOneD_Mach(AverageMach);
            solver_container->SetOneD_Temp(AverageTemperature);
            solver_container->SetOneD_MassFlowRate(MassFlowRate);

            solver_container->SetOneD_FluxAvgPress(PressureRef);
            solver_container->SetOneD_FluxAvgDensity(DensityRef);
            solver_container->SetOneD_FluxAvgVelocity(VelocityRef);
            solver_container->SetOneD_FluxAvgEntalpy(EnthalpyRef);

        }

        void Cae_Output::SetForceSections(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, Config_p config, unsigned long iExtIter) {

            short iSection, nSection;
            unsigned long iVertex, iPoint;
            double *Plane_P0, *Plane_Normal, MinPlane, MaxPlane, *CPressure, MinXCoord, MaxXCoord, Force[3], ForceInviscid[3],
                MomentInviscid[3] = { 0.0, 0.0, 0.0 }, MomentDist[3] = { 0.0, 0.0, 0.0 }, RefDensity, RefPressure, RefAreaCoeff, *Velocity_Inf, Gas_Constant, Mach2Vel, Mach_Motion, Gamma, RefVel2 = 0.0, factor, NDPressure, *Origin, RefLengthMoment, Alpha, Beta, CDrag_Inv, CLift_Inv, CMy_Inv;
            vector<double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil, Pressure_Airfoil;
            string Marker_Tag, Slice_Filename, Slice_Ext;
            ofstream Cp_File;
            unsigned short iDim;

            bool grid_movement = config->GetGrid_Movement();
            bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
            bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
            bool freesurface = (config->GetKind_Regime() == FREESURFACE);

            Plane_P0 = new double[3];
            Plane_Normal = new double[3];
            CPressure = new double[geometry->GetnPoint()];

            /*--- Compute some reference quantities and necessary values ---*/
            RefDensity = solver_container->GetDensity_Inf();
            RefPressure = solver_container->GetPressure_Inf();
            RefAreaCoeff = config->GetRefAreaCoeff();
            Velocity_Inf = solver_container->GetVelocity_Inf();
            Gamma = config->GetGamma();
            Origin = config->GetRefOriginMoment(0);
            RefLengthMoment = config->GetRefLengthMoment();
            Alpha = config->GetAoA()*PI_NUMBER / 180.0;
            Beta = config->GetAoS()*PI_NUMBER / 180.0;

            if (grid_movement) {
                Gas_Constant = config->GetGas_ConstantND();
                Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
                Mach_Motion = config->GetMach_Motion();
                RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
            }
            else {
                RefVel2 = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                    RefVel2 += Velocity_Inf[iDim] * Velocity_Inf[iDim];
            }
            factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            if (geometry->GetnDim() == 3) {

                /*--- Copy the pressure to an auxiliar structure ---*/

                for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
                    if (compressible) {
                        CPressure[iPoint] = (solver_container->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff;
                    }
                    if (incompressible || freesurface) {
                        CPressure[iPoint] = (solver_container->node[iPoint]->GetPressureInc() - RefPressure)*factor*RefAreaCoeff;
                    }
                }

                nSection = config->GetnSections();

                for (iSection = 0; iSection < nSection; iSection++) {

                    /*--- Read the values from the config file ---*/

                    MinPlane = config->GetSection_Location(0); MaxPlane = config->GetSection_Location(1);
                    MinXCoord = -1E6; MaxXCoord = 1E6;

                    Plane_Normal[0] = 0.0;    Plane_P0[0] = 0.0;
                    Plane_Normal[1] = 0.0;    Plane_P0[1] = 0.0;
                    Plane_Normal[2] = 0.0;    Plane_P0[2] = 0.0;

                    Plane_Normal[config->GetAxis_Orientation()] = 1.0;
                    Plane_P0[config->GetAxis_Orientation()] = MinPlane + iSection*(MaxPlane - MinPlane) / double(nSection - 1);

                    /*--- Compute the airfoil sections (note that we feed in the Cp) ---*/

                    geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal,
                        MinXCoord, MaxXCoord, CPressure,
                        Xcoord_Airfoil, Ycoord_Airfoil,
                        Zcoord_Airfoil, Pressure_Airfoil, true,
                        config);

                    if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() == 0)) {
                        cout << "Please check the config file, the section " << iSection + 1 << " has not been detected." << endl;
                    }

                    /*--- Output the pressure on each section (tecplot format) ---*/

                    if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() != 0)) {

                        /*--- Write Cp at each section ---*/

                        ofstream Cp_File;
                        if (iSection == 0) {
                            Cp_File.open("cp_sections.dat", ios::out);
                            Cp_File << "TITLE = \"Airfoil sections\"" << endl;
                            Cp_File << "VARIABLES = \"X\",\"Y\",\"Z\",\"Cp\"" << endl;
                        }
                        else Cp_File.open("cp_sections.dat", ios::app);

                        Cp_File << "ZONE T=\"SECTION_" << (iSection + 1) << "\", NODES= " << Xcoord_Airfoil.size() << ", ELEMENTS= " << Xcoord_Airfoil.size() - 1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;

                        /*--- Coordinates and pressure value ---*/

                        if (config->GetSystemMeasurements() == SI) {
                            for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
                                Cp_File << Xcoord_Airfoil[iVertex] << " " << Ycoord_Airfoil[iVertex] << " " << Zcoord_Airfoil[iVertex] << " " << Pressure_Airfoil[iVertex] << endl;
                            }
                        }
                        if (config->GetSystemMeasurements() == US) {
                            for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
                                Cp_File << Xcoord_Airfoil[iVertex] * 12.0 << " " << Ycoord_Airfoil[iVertex] * 12.0 << " " << Zcoord_Airfoil[iVertex] * 12.0 << " " << Pressure_Airfoil[iVertex] << endl;
                            }
                        }

                        /*--- Basic conectivity ---*/

                        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
                            Cp_File << iVertex << "\t" << iVertex + 1 << "\n";
                        }

                        Cp_File.close();


                        /*--- Compute load distribution ---*/

                        ForceInviscid[0] = 0.0; ForceInviscid[1] = 0.0; ForceInviscid[2] = 0.0; MomentInviscid[1] = 0.0;

                        for (iVertex = 0; iVertex < Xcoord_Airfoil.size() - 1; iVertex++) {

                            NDPressure = 0.5*(Pressure_Airfoil[iVertex] + Pressure_Airfoil[iVertex + 1]);

                            Force[0] = -(Zcoord_Airfoil[iVertex + 1] - Zcoord_Airfoil[iVertex])*NDPressure;
                            Force[1] = 0.0;
                            Force[2] = (Xcoord_Airfoil[iVertex + 1] - Xcoord_Airfoil[iVertex])*NDPressure;

                            ForceInviscid[0] += Force[0];
                            ForceInviscid[1] += Force[1];
                            ForceInviscid[2] += Force[2];

                            MomentDist[0] = 0.5*(Xcoord_Airfoil[iVertex] + Xcoord_Airfoil[iVertex + 1]) - Origin[0];
                            MomentDist[1] = 0.5*(Ycoord_Airfoil[iVertex] + Ycoord_Airfoil[iVertex + 1]) - Origin[1];
                            MomentDist[2] = 0.5*(Zcoord_Airfoil[iVertex] + Zcoord_Airfoil[iVertex + 1]) - Origin[3];

                            MomentInviscid[1] += (Force[0] * MomentDist[2] - Force[2] * MomentDist[0]) / RefLengthMoment;

                        }

                        CLift_Inv = fabs(-ForceInviscid[0] * sin(Alpha) + ForceInviscid[2] * cos(Alpha));
                        CDrag_Inv = fabs(ForceInviscid[0] * cos(Alpha)*cos(Beta) + ForceInviscid[1] * sin(Beta) + ForceInviscid[2] * sin(Alpha)*cos(Beta));
                        CMy_Inv = MomentInviscid[1];


                        /*--- Write load distribution ---*/

                        ofstream Load_File;
                        if (iSection == 0) {
                            Load_File.open("load_distribution.dat", ios::out);
                            Load_File << "TITLE = \"Load distribution\"" << endl;
                            Load_File << "VARIABLES = \"Y\",\"C<sub>L</sub>\",\"C<sub>D</sub>\",\"C<supb>My</sub>\"" << endl;
                            Load_File << "ZONE T=\"Wing load distribution\", NODES= " << nSection << ", ELEMENTS= " << nSection - 1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
                        }
                        else Load_File.open("load_distribution.dat", ios::app);

                        /*--- Coordinates and pressure value ---*/

                        Load_File << Ycoord_Airfoil[0] << " " << CLift_Inv << " " << CDrag_Inv << " " << CMy_Inv << endl;

                        /*--- Basic conectivity ---*/

                        if (iSection == nSection - 1) {
                            for (iSection = 1; iSection < nSection; iSection++) {
                                Load_File << iSection << "\t" << iSection + 1 << "\n";
                            }
                        }

                        Load_File.close();


                    }

                }


            }

            /*--- Delete dynamically allocated memory ---*/

            delete[] Plane_P0;
            delete[] Plane_Normal;
            delete[] CPressure;

        }

        void Cae_Output::SetCp_InverseDesign(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, Config_p config, unsigned long iExtIter) {

            unsigned short iMarker, icommas, Boundary, iDim;
            unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
            double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = NULL, Area, PressDiff;
            bool *PointInDomain;
            string text_line, surfCp_filename;
            ifstream Surface_file;
            char buffer[50], cstr[200];


            nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
            MPI_Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
            nPointGlobal = nPointLocal;
#endif

            Point2Vertex = new unsigned long[nPointGlobal][2];
            PointInDomain = new bool[nPointGlobal];

            for (iPoint = 0; iPoint < nPointGlobal; iPoint++)
                PointInDomain[iPoint] = false;

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Boundary = config->GetMarker_All_KindBC(iMarker);

                if ((Boundary == EULER_WALL) ||
                    (Boundary == HEAT_FLUX) ||
                    (Boundary == HEAT_FLUX_CATALYTIC) ||
                    (Boundary == HEAT_FLUX_NONCATALYTIC) ||
                    (Boundary == ISOTHERMAL) ||
                    (Boundary == ISOTHERMAL_CATALYTIC) ||
                    (Boundary == ISOTHERMAL_NONCATALYTIC) ||
                    (Boundary == NEARFIELD_BOUNDARY)) {
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

                        /*--- The Pressure file uses the global numbering ---*/

#ifndef HAVE_MPI
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
                        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif

                        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
                            Point2Vertex[iPoint][0] = iMarker;
                            Point2Vertex[iPoint][1] = iVertex;
                            PointInDomain[iPoint] = true;
                            solver_container->SetCPressureTarget(iMarker, iVertex, 0.0);
                        }

                    }
                }
            }

            /*--- Prepare to read the surface pressure files (CSV) ---*/

            surfCp_filename = "TargetCp";
            strcpy(cstr, surfCp_filename.c_str());

            /*--- Write file name with extension if unsteady or steady ---*/

            if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
                (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
                if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))    sprintf(buffer, "_0000%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))   sprintf(buffer, "_000%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))  sprintf(buffer, "_00%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.dat", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.dat", int(iExtIter));
            }
            else
                sprintf(buffer, ".dat");

            strcat(cstr, buffer);

            /*--- Read the surface pressure file ---*/

            string::size_type position;

            Surface_file.open(cstr, ios::in);

            if (!(Surface_file.fail())) {

                getline(Surface_file, text_line);

                while (getline(Surface_file, text_line)) {
                    for (icommas = 0; icommas < 50; icommas++) {
                        position = text_line.find(",", 0);
                        if (position != string::npos) text_line.erase(position, 1);
                    }
                    stringstream  point_line(text_line);

                    if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff;
                    if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff;

                    if (PointInDomain[iPoint]) {

                        /*--- Find the vertex for the Point and Marker ---*/

                        iMarker = Point2Vertex[iPoint][0];
                        iVertex = Point2Vertex[iPoint][1];

                        solver_container->SetCPressureTarget(iMarker, iVertex, PressureCoeff);

                    }

                }

                Surface_file.close();

            }

            /*--- Compute the pressure difference ---*/

            PressDiff = 0.0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Boundary = config->GetMarker_All_KindBC(iMarker);

                if ((Boundary == EULER_WALL) ||
                    (Boundary == HEAT_FLUX) ||
                    (Boundary == HEAT_FLUX_CATALYTIC) ||
                    (Boundary == HEAT_FLUX_NONCATALYTIC) ||
                    (Boundary == ISOTHERMAL) ||
                    (Boundary == ISOTHERMAL_CATALYTIC) ||
                    (Boundary == ISOTHERMAL_NONCATALYTIC) ||
                    (Boundary == NEARFIELD_BOUNDARY)) {
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

                        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

                        Cp = solver_container->GetCPressure(iMarker, iVertex);
                        CpTarget = solver_container->GetCPressureTarget(iMarker, iVertex);

                        Area = 0.0;
                        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                            Area += Normal[iDim] * Normal[iDim];
                        Area = sqrt(Area);

                        PressDiff += Area * (CpTarget - Cp) * (CpTarget - Cp);
                    }

                }
            }

#ifdef HAVE_MPI
            double MyPressDiff = PressDiff;   PressDiff = 0.0;
            MPI_Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

            /*--- Update the total Cp difference coeffient ---*/

            solver_container->SetTotal_CpDiff(PressDiff);

            delete[] Point2Vertex;

        }

        void Cae_Output::SetHeat_InverseDesign(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, Config_p config, unsigned long iExtIter) {

            unsigned short iMarker, icommas, Boundary, iDim;
            unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
            double XCoord, YCoord, ZCoord, PressureCoeff, HeatFlux = 0.0, HeatFluxDiff, HeatFluxTarget, *Normal = NULL, Area,
                Pressure, Cf;
            bool *PointInDomain;
            string text_line, surfHeatFlux_filename;
            ifstream Surface_file;
            char buffer[50], cstr[200];


            nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
            MPI_Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
            nPointGlobal = nPointLocal;
#endif

            Point2Vertex = new unsigned long[nPointGlobal][2];
            PointInDomain = new bool[nPointGlobal];

            for (iPoint = 0; iPoint < nPointGlobal; iPoint++)
                PointInDomain[iPoint] = false;

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Boundary = config->GetMarker_All_KindBC(iMarker);

                if ((Boundary == EULER_WALL) ||
                    (Boundary == HEAT_FLUX) ||
                    (Boundary == HEAT_FLUX_CATALYTIC) ||
                    (Boundary == HEAT_FLUX_NONCATALYTIC) ||
                    (Boundary == ISOTHERMAL) ||
                    (Boundary == ISOTHERMAL_CATALYTIC) ||
                    (Boundary == ISOTHERMAL_NONCATALYTIC) ||
                    (Boundary == NEARFIELD_BOUNDARY)) {
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

                        /*--- The Pressure file uses the global numbering ---*/

#ifndef HAVE_MPI
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
                        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif

                        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
                            Point2Vertex[iPoint][0] = iMarker;
                            Point2Vertex[iPoint][1] = iVertex;
                            PointInDomain[iPoint] = true;
                            solver_container->SetHeatFluxTarget(iMarker, iVertex, 0.0);
                        }
                    }
                }
            }

            /*--- Prepare to read the surface pressure files (CSV) ---*/

            surfHeatFlux_filename = "TargetHeatFlux";
            strcpy(cstr, surfHeatFlux_filename.c_str());

            /*--- Write file name with extension if unsteady or steady ---*/

            if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
                (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
                if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))    sprintf(buffer, "_0000%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))   sprintf(buffer, "_000%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))  sprintf(buffer, "_00%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.dat", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.dat", int(iExtIter));
            }
            else
                sprintf(buffer, ".dat");

            strcat(cstr, buffer);

            /*--- Read the surface pressure file ---*/

            string::size_type position;

            Surface_file.open(cstr, ios::in);

            if (!(Surface_file.fail())) {

                getline(Surface_file, text_line);

                while (getline(Surface_file, text_line)) {
                    for (icommas = 0; icommas < 50; icommas++) {
                        position = text_line.find(",", 0);
                        if (position != string::npos) text_line.erase(position, 1);
                    }
                    stringstream  point_line(text_line);

                    if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;
                    if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;

                    if (PointInDomain[iPoint]) {

                        /*--- Find the vertex for the Point and Marker ---*/

                        iMarker = Point2Vertex[iPoint][0];
                        iVertex = Point2Vertex[iPoint][1];

                        solver_container->SetHeatFluxTarget(iMarker, iVertex, HeatFlux);

                    }

                }

                Surface_file.close();
            }

            /*--- Compute the pressure difference ---*/

            HeatFluxDiff = 0.0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Boundary = config->GetMarker_All_KindBC(iMarker);

                if ((Boundary == EULER_WALL) ||
                    (Boundary == HEAT_FLUX) ||
                    (Boundary == HEAT_FLUX_CATALYTIC) ||
                    (Boundary == HEAT_FLUX_NONCATALYTIC) ||
                    (Boundary == ISOTHERMAL) ||
                    (Boundary == ISOTHERMAL_CATALYTIC) ||
                    (Boundary == ISOTHERMAL_NONCATALYTIC) ||
                    (Boundary == NEARFIELD_BOUNDARY)) {
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

                        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

                        HeatFlux = solver_container->GetHeatFlux(iMarker, iVertex);
                        HeatFluxTarget = solver_container->GetHeatFluxTarget(iMarker, iVertex);

                        Area = 0.0;
                        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                            Area += Normal[iDim] * Normal[iDim];
                        Area = sqrt(Area);

                        HeatFluxDiff += Area * (HeatFluxTarget - HeatFlux) * (HeatFluxTarget - HeatFlux);

                    }

                }
            }

#ifdef HAVE_MPI
            double MyHeatFluxDiff = HeatFluxDiff;   HeatFluxDiff = 0.0;
            MPI_Allreduce(&MyHeatFluxDiff, &HeatFluxDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

            /*--- Update the total HeatFlux difference coeffient ---*/

            solver_container->SetTotal_HeatFluxDiff(HeatFluxDiff);

            delete[] Point2Vertex;

        }

        void Cae_Output::SetEquivalentArea(SOLV::SOLV_Solver *solver_container, GEOM::GEOM_Geometry *geometry, Config_p config, unsigned long iExtIter) {

            ofstream EquivArea_file, FuncGrad_file;
            unsigned short iMarker = 0, iDim;
            short *AzimuthalAngle = NULL;
            double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign = 0.0, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion,
                *Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf,
                ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
                *Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
                *Weight = NULL, jFunction, jp1Function;
            unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint,
                *IdPoint = NULL, *IdDomain = NULL, auxDomain;
            unsigned short iPhiAngle;
            ofstream NearFieldEA_file; ifstream TargetEA_file;

            double XCoordBegin_OF = config->GetEA_IntLimit(0);
            double XCoordEnd_OF = config->GetEA_IntLimit(1);

            unsigned short nDim = geometry->GetnDim();
            double AoA = -(config->GetAoA()*PI_NUMBER / 180.0);
            double EAScaleFactor = config->GetEA_ScaleFactor(); // The EA Obj. Func. should be ~ force based Obj. Func.

            int rank = MESH_0;

            Mach = config->GetMach();
            Gamma = config->GetGamma();
            Beta = sqrt(Mach*Mach - 1.0);
            R_Plane = fabs(config->GetEA_IntLimit(2));
            Pressure_Inf = config->GetPressure_FreeStreamND();
            Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
            Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
            Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
            ModVelocity_Inf = 0;
            for (iDim = 0; iDim < 3; iDim++)
                ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];

            factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);

#ifndef HAVE_MPI

            /*--- Compute the total number of points on the near-field ---*/

            nVertex_NearField = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                        Coord = geometry->node[iPoint]->GetCoord();

                        /*--- Using Face_Normal(z), and Coord(z) we identify only a surface,
                        note that there are 2 NEARFIELD_BOUNDARY surfaces ---*/

                        if ((Face_Normal[nDim - 1] > 0.0) && (Coord[nDim - 1] < 0.0)) nVertex_NearField++;
                    }

            /*--- Create an array with all the coordinates, points, pressures, face area,
            equivalent area, and nearfield weight ---*/

            Xcoord = new double[nVertex_NearField];
            Ycoord = new double[nVertex_NearField];
            Zcoord = new double[nVertex_NearField];
            AzimuthalAngle = new short[nVertex_NearField];
            IdPoint = new unsigned long[nVertex_NearField];
            IdDomain = new unsigned long[nVertex_NearField];
            Pressure = new double[nVertex_NearField];
            FaceArea = new double[nVertex_NearField];
            EquivArea = new double[nVertex_NearField];
            TargetArea = new double[nVertex_NearField];
            NearFieldWeight = new double[nVertex_NearField];
            Weight = new double[nVertex_NearField];

            /*--- Copy the boundary information to an array ---*/

            nVertex_NearField = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                        Coord = geometry->node[iPoint]->GetCoord();

                        if ((Face_Normal[nDim - 1] > 0.0) && (Coord[nDim - 1] < 0.0)) {

                            IdPoint[nVertex_NearField] = iPoint;
                            Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
                            Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);

                            if (nDim == 2) {
                                AzimuthalAngle[nVertex_NearField] = 0;
                            }

                            if (nDim == 3) {
                                Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);

                                /*--- Rotate the nearfield cylinder (AoA) only 3D ---*/

                                double YcoordRot = Ycoord[nVertex_NearField];
                                double ZcoordRot = Xcoord[nVertex_NearField] * sin(AoA) + Zcoord[nVertex_NearField] * cos(AoA);

                                /*--- Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/

                                double AngleDouble; short AngleInt;
                                AngleDouble = fabs(atan(-YcoordRot / ZcoordRot)*180.0 / PI_NUMBER);

                                /*--- Fix an azimuthal line due to misalignments of the near-field ---*/

                                double FixAzimuthalLine = config->GetFixAzimuthalLine();

                                if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1)) AngleDouble = FixAzimuthalLine - 0.1;

                                AngleInt = (short)floor(AngleDouble + 0.5);
                                if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
                                else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
                            }

                            if (AzimuthalAngle[nVertex_NearField] <= 60) {
                                Pressure[nVertex_NearField] = solver_container->node[iPoint]->GetPressure();
                                FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim - 1]);
                                nVertex_NearField++;
                            }

                        }
                    }

#else

            int nProcessor;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
            int iProcessor;

            unsigned long *Buffer_Receive_nVertex = NULL;
            if (rank == MASTER_NODE) {
                Buffer_Receive_nVertex = new unsigned long[nProcessor];
            }

            /*--- Compute the total number of points of the near-field ghost nodes ---*/

            nLocalVertex_NearField = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                        Coord = geometry->node[iPoint]->GetCoord();

                        if (geometry->node[iPoint]->GetDomain())
                            if ((Face_Normal[nDim - 1] > 0.0) && (Coord[nDim - 1] < 0.0))
                                nLocalVertex_NearField++;
                    }

            unsigned long *Buffer_Send_nVertex = new unsigned long[1];
            Buffer_Send_nVertex[0] = nLocalVertex_NearField;

            /*--- Send Near-Field vertex information --*/

            MPI_Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
            delete[] Buffer_Send_nVertex;

            double *Buffer_Send_Xcoord = new double[MaxLocalVertex_NearField];
            double *Buffer_Send_Ycoord = new double[MaxLocalVertex_NearField];
            double *Buffer_Send_Zcoord = new double[MaxLocalVertex_NearField];
            unsigned long *Buffer_Send_IdPoint = new unsigned long[MaxLocalVertex_NearField];
            double *Buffer_Send_Pressure = new double[MaxLocalVertex_NearField];
            double *Buffer_Send_FaceArea = new double[MaxLocalVertex_NearField];

            double *Buffer_Receive_Xcoord = NULL;
            double *Buffer_Receive_Ycoord = NULL;
            double *Buffer_Receive_Zcoord = NULL;
            unsigned long *Buffer_Receive_IdPoint = NULL;
            double *Buffer_Receive_Pressure = NULL;
            double *Buffer_Receive_FaceArea = NULL;

            if (rank == MASTER_NODE) {
                Buffer_Receive_Xcoord = new double[nProcessor*MaxLocalVertex_NearField];
                Buffer_Receive_Ycoord = new double[nProcessor*MaxLocalVertex_NearField];
                Buffer_Receive_Zcoord = new double[nProcessor*MaxLocalVertex_NearField];
                Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
                Buffer_Receive_Pressure = new double[nProcessor*MaxLocalVertex_NearField];
                Buffer_Receive_FaceArea = new double[nProcessor*MaxLocalVertex_NearField];
            }

            unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;
            unsigned long nBuffer_Ycoord = MaxLocalVertex_NearField;
            unsigned long nBuffer_Zcoord = MaxLocalVertex_NearField;
            unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
            unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
            unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;

            for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
                Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
                Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
                Buffer_Send_Ycoord[iVertex] = 0.0; Buffer_Send_Zcoord[iVertex] = 0.0;
            }

            /*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/

            nLocalVertex_NearField = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                        Coord = geometry->node[iPoint]->GetCoord();

                        if (geometry->node[iPoint]->GetDomain())
                            if ((Face_Normal[nDim - 1] > 0.0) && (Coord[nDim - 1] < 0.0)) {
                                Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
                                Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
                                Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
                                Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
                                Buffer_Send_Pressure[nLocalVertex_NearField] = solver_container->node[iPoint]->GetPressure();
                                Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim - 1]);
                                nLocalVertex_NearField++;
                            }
                    }

            /*--- Send all the information --*/

            MPI_Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI_DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
            delete[] Buffer_Send_Xcoord;
            delete[] Buffer_Send_Ycoord;
            delete[] Buffer_Send_Zcoord;
            delete[] Buffer_Send_IdPoint;
            delete[] Buffer_Send_Pressure;
            delete[] Buffer_Send_FaceArea;

            if (rank == MASTER_NODE) {

                Xcoord = new double[nVertex_NearField];
                Ycoord = new double[nVertex_NearField];
                Zcoord = new double[nVertex_NearField];
                AzimuthalAngle = new short[nVertex_NearField];
                IdPoint = new unsigned long[nVertex_NearField];
                IdDomain = new unsigned long[nVertex_NearField];
                Pressure = new double[nVertex_NearField];
                FaceArea = new double[nVertex_NearField];
                EquivArea = new double[nVertex_NearField];
                TargetArea = new double[nVertex_NearField];
                NearFieldWeight = new double[nVertex_NearField];
                Weight = new double[nVertex_NearField];

                nVertex_NearField = 0;
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
                        Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField + iVertex];
                        Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField + iVertex];

                        if (nDim == 2) {
                            AzimuthalAngle[nVertex_NearField] = 0;
                        }

                        if (nDim == 3) {
                            Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField + iVertex];

                            /*--- Rotate the nearfield cylinder  ---*/

                            double YcoordRot = Ycoord[nVertex_NearField];
                            double ZcoordRot = Xcoord[nVertex_NearField] * sin(AoA) + Zcoord[nVertex_NearField] * cos(AoA);

                            /*--- Compute the Azimuthal angle ---*/

                            double AngleDouble; short AngleInt;
                            AngleDouble = fabs(atan(-YcoordRot / ZcoordRot)*180.0 / PI_NUMBER);

                            /*--- Fix an azimuthal line due to misalignments of the near-field ---*/

                            double FixAzimuthalLine = config->GetFixAzimuthalLine();

                            if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1))
                                AngleDouble = FixAzimuthalLine - 0.1;

                            AngleInt = (short)floor(AngleDouble + 0.5);

                            if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
                            else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
                        }

                        if (AzimuthalAngle[nVertex_NearField] <= 60) {
                            IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField + iVertex];
                            Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField + iVertex];
                            FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField + iVertex];
                            IdDomain[nVertex_NearField] = iProcessor;
                            nVertex_NearField++;
                        }

                    }

                delete[] Buffer_Receive_nVertex;

                delete[] Buffer_Receive_Xcoord;
                delete[] Buffer_Receive_Ycoord;
                delete[] Buffer_Receive_Zcoord;
                delete[] Buffer_Receive_IdPoint;
                delete[] Buffer_Receive_Pressure;
                delete[] Buffer_Receive_FaceArea;

            }

#endif

            if (rank == MASTER_NODE) {

                vector<short> PhiAngleList;
                vector<short>::iterator IterPhiAngleList;

                for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
                    PhiAngleList.push_back(AzimuthalAngle[iVertex]);

                sort(PhiAngleList.begin(), PhiAngleList.end());
                IterPhiAngleList = unique(PhiAngleList.begin(), PhiAngleList.end());
                PhiAngleList.resize(IterPhiAngleList - PhiAngleList.begin());

                /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/

                vector<vector<double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size());
                vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size());
                vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size());
                vector<vector<double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());

                /*--- Distribute the values among the different PhiAngles ---*/

                for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
                    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                        if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
                            Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
                            Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
                            Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
                            IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
                            IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
                            Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
                            FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
                            EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
                            TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
                            NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
                            Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
                        }

                /*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/

                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                    for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
                        for (jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
                            if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex + 1]) {
                                auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex + 1]; Xcoord_PhiAngle[iPhiAngle][jVertex + 1] = auxXCoord;
                                auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex + 1]; Ycoord_PhiAngle[iPhiAngle][jVertex + 1] = auxYCoord;
                                auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex + 1]; Zcoord_PhiAngle[iPhiAngle][jVertex + 1] = auxZCoord;
                                auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex + 1]; Pressure_PhiAngle[iPhiAngle][jVertex + 1] = auxPress;
                                auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex + 1]; FaceArea_PhiAngle[iPhiAngle][jVertex + 1] = auxArea;
                                auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex + 1]; IdPoint_PhiAngle[iPhiAngle][jVertex + 1] = auxPoint;
                                auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex + 1]; IdDomain_PhiAngle[iPhiAngle][jVertex + 1] = auxDomain;
                            }


                /*--- Check that all the azimuth lists have the same size ---*/

                unsigned short nVertex = Xcoord_PhiAngle[0].size();
                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
                    unsigned short nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
                    if (nVertex_aux != nVertex) cout << "Be careful!!! one azimuth list is shorter than the other" << endl;
                    nVertex = min(nVertex, nVertex_aux);
                }

                /*--- Compute equivalent area distribution at each azimuth angle ---*/

                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
                    EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
                    for (iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
                        EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;

                        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex] * cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex] * sin(AoA);

                        for (jVertex = 0; jVertex < iVertex - 1; jVertex++) {

                            Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex] * cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex] * sin(AoA);
                            jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex + 1] * cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex + 1] * sin(AoA);

                            jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i - Coord_j);
                            jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex + 1] - Pressure_Inf)*sqrt(Coord_i - jp1Coord);

                            DeltaX = (jp1Coord - Coord_j);
                            MeanFuntion = 0.5*(jp1Function + jFunction);
                            EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFuntion;
                        }
                    }
                }

                /*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/

                NearFieldEA_file.precision(15);
                NearFieldEA_file.open("Equivalent_Area.dat", ios::out);
                NearFieldEA_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << endl;

                if (config->GetSystemMeasurements() == US)
                    NearFieldEA_file << "VARIABLES = \"Height (in) at r=" << R_Plane*12.0 << " in. (cyl. coord. system)\"";
                else
                    NearFieldEA_file << "VARIABLES = \"Height (m) at r=" << R_Plane << " m. (cylindrical coordinate system)\"";

                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
                    if (config->GetSystemMeasurements() == US)
                        NearFieldEA_file << ", \"Equivalent Area (ft<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
                    else
                        NearFieldEA_file << ", \"Equivalent Area (m<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
                }

                NearFieldEA_file << endl;
                for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {

                    double XcoordRot = Xcoord_PhiAngle[0][iVertex] * cos(AoA) - Zcoord_PhiAngle[0][iVertex] * sin(AoA);
                    double XcoordRot_init = Xcoord_PhiAngle[0][0] * cos(AoA) - Zcoord_PhiAngle[0][0] * sin(AoA);

                    if (config->GetSystemMeasurements() == US)
                        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
                    else
                        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init);

                    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
                        NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
                    }

                    NearFieldEA_file << endl;

                }
                NearFieldEA_file.close();

                /*--- Read target equivalent area from the configuration file,
                this first implementation requires a complete table (same as the original
                EA table). so... no interpolation. ---*/

                vector<vector<double> > TargetArea_PhiAngle_Trans;
                TargetEA_file.open("TargetEA.dat", ios::in);

                if (TargetEA_file.fail()) {
                    if (iExtIter == 0) {
                        cout << "There is no Target Equivalent Area file (TargetEA.dat)!!" << endl;
                        cout << "Using default parameters (Target Equiv Area = 0.0)" << endl;
                    }
                    /*--- Set the table to 0 ---*/
                    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                        for (iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
                            TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
                }
                else {

                    /*--- skip header lines ---*/

                    string line;
                    getline(TargetEA_file, line);
                    getline(TargetEA_file, line);

                    while (TargetEA_file) {

                        string line;
                        getline(TargetEA_file, line);
                        istringstream is(line);
                        vector<double> row;
                        unsigned short iter = 0;

                        while (is.good()) {
                            string token;
                            getline(is, token, ',');

                            istringstream js(token);

                            double data;
                            js >> data;

                            /*--- The first element in the table is the coordinate (in or m)---*/

                            if (iter != 0) row.push_back(data);
                            iter++;

                        }
                        TargetArea_PhiAngle_Trans.push_back(row);
                    }

                    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                        for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
                            TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];

                }

                /*--- Divide by the number of Phi angles in the nearfield ---*/

                double PhiFactor = 1.0 / double(PhiAngleList.size());

                /*--- Evaluate the objective function ---*/

                InverseDesign = 0;
                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                    for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
                        Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
                        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];

                        double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex] - TargetArea_PhiAngle[iPhiAngle][iVertex];
                        double percentage = fabs(Difference) * 100 / fabs(TargetArea_PhiAngle[iPhiAngle][iVertex]);

                        if ((percentage < 0.1) || (Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;

                        InverseDesign += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex] * Difference*Difference;

                    }

                /*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/

                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                    for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
                        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
                        NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
                        for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
                            Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
                            Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;

                            double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex] - TargetArea_PhiAngle[iPhiAngle][jVertex];
                            double percentage = fabs(Difference) * 100 / fabs(TargetArea_PhiAngle[iPhiAngle][jVertex]);

                            if ((percentage < 0.1) || (Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;

                            NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex] * 2.0*Difference*factor*sqrt(Coord_j - Coord_i);
                        }
                    }

                /*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/

                EquivArea_file.precision(15);
                EquivArea_file.open("nearfield_flow.dat", ios::out);
                EquivArea_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << endl;

                if (config->GetSystemMeasurements() == US)
                    EquivArea_file << "VARIABLES = \"Height (in) at r=" << R_Plane*12.0 << " in. (cyl. coord. system)\",\"Equivalent Area (ft<sup>2</sup>)\",\"Target Equivalent Area (ft<sup>2</sup>)\",\"Cp\"" << endl;
                else
                    EquivArea_file << "VARIABLES = \"Height (m) at r=" << R_Plane << " m. (cylindrical coordinate system)\",\"Equivalent Area (m<sup>2</sup>)\",\"Target Equivalent Area (m<sup>2</sup>)\",\"Cp\"" << endl;

                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
                    EquivArea_file << fixed << "ZONE T= \"<greek>F</greek>=" << PhiAngleList[iPhiAngle] << " deg.\"" << endl;
                    for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {

                        double XcoordRot = Xcoord_PhiAngle[0][iVertex] * cos(AoA) - Zcoord_PhiAngle[0][iVertex] * sin(AoA);
                        double XcoordRot_init = Xcoord_PhiAngle[0][0] * cos(AoA) - Zcoord_PhiAngle[0][0] * sin(AoA);

                        if (config->GetSystemMeasurements() == US)
                            EquivArea_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
                        else
                            EquivArea_file << scientific << (XcoordRot - XcoordRot_init);

                        EquivArea_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
                            << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << (Pressure_PhiAngle[iPhiAngle][iVertex] - Pressure_Inf) / Pressure_Inf << endl;
                    }
                }

                EquivArea_file.close();

                /*--- Write Weight file for adjoint computation ---*/

                FuncGrad_file.precision(15);
                FuncGrad_file.open("WeightNF.dat", ios::out);

                FuncGrad_file << scientific << "-1.0";
                for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                    FuncGrad_file << scientific << "\t" << PhiAngleList[iPhiAngle];
                FuncGrad_file << endl;

                for (iVertex = 0; iVertex < NearFieldWeight_PhiAngle[0].size(); iVertex++) {
                    double XcoordRot = Xcoord_PhiAngle[0][iVertex] * cos(AoA) - Zcoord_PhiAngle[0][iVertex] * sin(AoA);
                    FuncGrad_file << scientific << XcoordRot;
                    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
                        FuncGrad_file << scientific << "\t" << NearFieldWeight_PhiAngle[iPhiAngle][iVertex];
                    FuncGrad_file << endl;
                }
                FuncGrad_file.close();

                /*--- Delete structures ---*/

                delete[] Xcoord; delete[] Ycoord; delete[] Zcoord;
                delete[] AzimuthalAngle; delete[] IdPoint; delete[] IdDomain;
                delete[] Pressure; delete[] FaceArea;
                delete[] EquivArea; delete[] TargetArea;
                delete[] NearFieldWeight; delete[] Weight;

            }

#ifndef HAVE_MPI

            /*--- Store the value of the NearField coefficient ---*/

            solver_container->SetTotal_CEquivArea(InverseDesign);

#else

            /*--- Send the value of the NearField coefficient to all the processors ---*/

            MPI_Bcast(&InverseDesign, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

            /*--- Store the value of the NearField coefficient ---*/

            solver_container->SetTotal_CEquivArea(InverseDesign);

#endif

        }

        void Cae_Output::SetCGNS_Coordinates(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short iZone) {

#ifdef HAVE_CGNS

            /*--- local CGNS variables ---*/
            int cgns_file, cgns_coord, element_dims, physical_dims, cgns_err;
            unsigned long iExtIter = config->GetExtIter();
            string base_file, buffer, elements_name;
            stringstream name, results_file;
            bool unsteady = config->GetUnsteady_Simulation();
            cgsize_t isize[3][1];

            /*--- Create CGNS base file name ---*/
            base_file = config->GetFlow_FileName();

            /*--- Add CGNS extension. ---*/
            base_file = base_file.append(".cgns");

            /*--- Create CGNS results file name ---*/
            if (unsteady) {

                buffer = config->GetFlow_FileName();

                results_file.str(string()); results_file << buffer;
                if (((int)iExtIter >= 0) && ((int)iExtIter < 10))			results_file << "_0000" << iExtIter;
                if (((int)iExtIter >= 10) && ((int)iExtIter < 100))		results_file << "_000" << iExtIter;
                if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))		results_file << "_00" << iExtIter;
                if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))	results_file << "_0" << iExtIter;
                if ((int)iExtIter >= 10000)							results_file << iExtIter;
                results_file << ".cgns";
            }

            /*--- Write base file if not already done ---*/
            if (!d_isBaseOutput) {

                /*--- Write base file ---*/
                cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_MODIFY, &cgns_file);
                if (cgns_err) cg_error_print();

                element_dims = geometry->GetnDim();		// Currently (release 2.0) only all-2D or all-3D zones permitted
                physical_dims = element_dims;

                isize[0][0] = (cgsize_t)d_numGlobalPoints;				// vertex size
                isize[1][0] = (cgsize_t)nGlobal_Elem;				// cell size
                isize[2][0] = 0;						// boundary vertex size (zero if elements not sorted)

                cgns_err = cg_goto(cgns_file, cgns_base, "Zone_t", cgns_zone, "end");
                if (cgns_err) cg_error_print();
                //    
                //    cgns_err = cg_goto(cgns_file, cgns_base, cgns_zone,"end");
                //		if (cgns_err) cg_error_print();

                /*--- write CGNS node coordinates ---*/
                cgns_err = cg_coord_write(cgns_file, cgns_base, cgns_zone, RealDouble, "x", Coords[0], &cgns_coord);
                if (cgns_err) cg_error_print();
                cgns_err = cg_coord_write(cgns_file, cgns_base, cgns_zone, RealDouble, "y", Coords[1], &cgns_coord);
                if (cgns_err) cg_error_print();
                if (geometry->GetnDim() == 3) {
                    cgns_err = cg_coord_write(cgns_file, cgns_base, cgns_zone, RealDouble, "z", Coords[2], &cgns_coord);
                    if (cgns_err) cg_error_print();
                }

                cgns_err = cg_close(cgns_file);
                if (cgns_err) cg_error_print();

                d_isBaseOutput = true;

            }

            /*--- Set up results file for this time step if necessary ---*/
            if (unsteady) {

                cgns_err = cg_open((char *)results_file.str().c_str(), CG_MODE_WRITE, &cgns_file);

                element_dims = geometry->GetnDim();		// Currently (release 3.2.9 "eagle") only all-2D or all-3D zones permitted
                physical_dims = element_dims;

                /*--- write CGNS base data (one base assumed as of version 3.2.9 "eagle") ---*/
                cgns_err = cg_base_write(cgns_file, "SU2 Base", element_dims, physical_dims, &cgns_base_results);
                if (cgns_err) cg_error_print();

                isize[0][0] = (cgsize_t)geometry->GetGlobal_nPointDomain();				// vertex size
                isize[1][0] = (cgsize_t)nGlobal_Elem;				// cell size
                isize[2][0] = 0;						// boundary vertex size (zero if elements not sorted)

                /*--- write CGNS zone data ---*/
                cgns_err = cg_zone_write(cgns_file, cgns_base_results, "SU2 Zone", isize[0], Unstructured, &cgns_zone_results);
                if (cgns_err) cg_error_print();

                cgns_err = cg_goto(cgns_file, cgns_base_results, "Zone_t", cgns_zone_results, "end");
                if (cgns_err) cg_error_print();

                /*--- Write CGNS node coordinates, if appliciable ---*/
                if (config->GetGrid_Movement()) {

                    /*--- write CGNS node coordinates ---*/
                    cgns_err = cg_coord_write(cgns_file, cgns_base_results, cgns_zone_results, RealDouble, "x", Coords[0], &cgns_coord);
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_coord_write(cgns_file, cgns_base_results, cgns_zone_results, RealDouble, "y", Coords[1], &cgns_coord);
                    if (cgns_err) cg_error_print();
                    if (geometry->GetnDim() == 3) {
                        cgns_err = cg_coord_write(cgns_file, cgns_base_results, cgns_zone_results, RealDouble, "z", Coords[2], &cgns_coord);
                        if (cgns_err) cg_error_print();
                    }
                }
                else {
                    /*--- Write a CGNS link for the node coordinates ---*/
                    cgns_err = cg_link_write("GridCoordinates", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/GridCoordinates");
                    if (cgns_err) cg_error_print();
                }

                /*--- Write a CGNS link for each element type connectivity ---*/
                if (d_numGlobalTrias > 0) cgns_err = cg_link_write("Triangle Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Triangle Elements");
                if (d_numGlobalQuads > 0) cgns_err = cg_link_write("Quadrilateral Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Quadrilateral Elements");
                if (d_numGlobalTetrs > 0) cgns_err = cg_link_write("Tetrahedral Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Tetrahedral Elements");
                if (d_numGlobalHexas > 0) cgns_err = cg_link_write("Hexahedral Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Hexahedral Elements");
                if (d_numGlobalPyras > 0) cgns_err = cg_link_write("Pyramid Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Pyramid Elements");
                if (d_numGlobalPriss > 0) cgns_err = cg_link_write("Prism Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Prism Elements");
                if (d_numGlobalBoundLines > 0) cgns_err = cg_link_write("Line Elements", (char *)base_file.c_str(), "/SU2 Base/SU2 Zone/Line Elements");
                if (cgns_err) cg_error_print();


                /*--- Close CGNS file ---*/
                cgns_err = cg_close(cgns_file);
                if (cgns_err) cg_error_print();

            }



#else // Not built with CGNS support

            cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n";

#endif

        }

        void Cae_Output::SetCGNS_Connectivity(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short iZone) {

#ifdef HAVE_CGNS

            /*--- local CGNS variables ---*/
            int cgns_file, element_dims, physical_dims, cgns_err;
            int cgns_section;
            unsigned long iExtIter = config->GetExtIter();
            string base_file, buffer, elements_name;
            stringstream name, results_file;
            bool unsteady = config->GetUnsteady_Simulation();
            cgsize_t isize[3][1], elem_start, elem_end, N;

            /*--- Create CGNS base file name ---*/
            base_file = config->GetFlow_FileName();

            /*--- Add CGNS extension. ---*/
            base_file = base_file.append(".cgns");

            /*--- Create CGNS results file name ---*/
            if (unsteady) {

                buffer = config->GetFlow_FileName();

                results_file.str(string()); results_file << buffer;
                if (((int)iExtIter >= 0) && ((int)iExtIter < 10))			results_file << "_0000" << iExtIter;
                if (((int)iExtIter >= 10) && ((int)iExtIter < 100))		results_file << "_000" << iExtIter;
                if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))		results_file << "_00" << iExtIter;
                if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))	results_file << "_0" << iExtIter;
                if ((int)iExtIter >= 10000)							results_file << iExtIter;
                results_file << ".cgns";
            }

            /*--- Write base file if not already done ---*/
            if (!d_isBaseOutput) {

                /*--- Write base file ---*/
                cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_WRITE, &cgns_file);
                if (cgns_err) cg_error_print();

                element_dims = geometry->GetnDim();		// Currently (release 2.0) only all-2D or all-3D zones permitted
                physical_dims = element_dims;

                /*--- write CGNS base data (one base assumed as of version 3.2.9 "eagle") ---*/
                cgns_err = cg_base_write(cgns_file, "SU2 Base", element_dims, physical_dims, &cgns_base);
                if (cgns_err) cg_error_print();

                /*--- write CGNS descriptor data ---*/
                cgns_err = cg_goto(cgns_file, cgns_base, "end");
                if (cgns_err) cg_error_print();

                cgns_err = cg_equationset_write(physical_dims);
                if (cgns_err) cg_error_print();

                /*--- Write governing equations to CGNS file ---*/
                cgns_err = cg_goto(cgns_file, cgns_base, "FlowEquationSet_t", 1, "end");
                if (cgns_err) cg_error_print();
                switch (config->GetKind_Solver()) {
                case EULER:
                    cgns_err = cg_governing_write(Euler); break;
                case NAVIER_STOKES:
                    cgns_err = cg_governing_write(NSLaminar); break;
                case RANS:
                    cgns_err = cg_governing_write(NSTurbulent); break;
                default:
                    break; // cgns_err = cg_governing_write(CG_UserDefined);
                }
                if (cgns_err) cg_error_print();

                if (unsteady) cgns_err = cg_simulation_type_write(cgns_file, cgns_base, TimeAccurate);
                else cgns_err = cg_simulation_type_write(cgns_file, cgns_base, NonTimeAccurate);
                if (cgns_err) cg_error_print();

                cgns_err = cg_descriptor_write("Solver Information", "SU2 version 3.2.9 \"eagle\"");
                if (cgns_err) cg_error_print();

                isize[0][0] = (cgsize_t)geometry->GetGlobal_nPointDomain(); //;				// vertex size
                isize[1][0] = (cgsize_t)nGlobal_Elem;				// cell size
                isize[2][0] = 0;						// boundary vertex size (zero if elements not sorted)

                /*--- write CGNS zone data ---*/
                cgns_err = cg_zone_write(cgns_file, cgns_base, "SU2 Zone", isize[0], Unstructured, &cgns_zone);
                if (cgns_err) cg_error_print();

                cgns_err = cg_goto(cgns_file, cgns_base, "Zone_t", cgns_zone, "end");
                if (cgns_err) cg_error_print();

                /*--- Reference Note: CGNS element type list:
                NODE, BAR_2, BAR_3, TRI_3, TRI_6, QUAD_4, QUAD_8, QUAD_9, TETRA_4, TETRA_10, PYRA_5,
                PYRA_14, PENTA_6, PENTA_15, PENTA_18, HEXA_8, HEXA_20, HEXA_27, MIXED, PYRA_13, NGON_n, NFACE_n ---*/

                /*--- Write a CGNS section for each element type ---*/
                // ier = cg_section_write(int fn, int B, int Z, char *ElementSectionName, ElementType_t type,
                // cgsize_t start, cgsize_t end, int nbndry, cgsize_t *Elements, int *S);

                if (d_numGlobalTrias > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalTrias;
                    N = (int)d_numGlobalTrias*N_POINTS_TRIANGLE;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone,
                        "Triangle Elements", TRI_3, elem_start, elem_end,
                        0, (cgsize_t *)Conn_Tria, &cgns_section);
                }
                if (d_numGlobalQuads > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalQuads; N = (int)d_numGlobalQuads*N_POINTS_QUADRILATERAL;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone, "Quadrilateral Elements", QUAD_4,
                        elem_start, elem_end, 0, (cgsize_t *)Conn_Quad, &cgns_section);
                }
                if (d_numGlobalTetrs > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalTetrs; N = (int)d_numGlobalTetrs*N_POINTS_TETRAHEDRON;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone, "Tetrahedral Elements", TETRA_4,
                        elem_start, elem_end, 0, (cgsize_t *)Conn_Tetr, &cgns_section);
                }
                if (d_numGlobalHexas > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalHexas; N = (int)d_numGlobalHexas*N_POINTS_HEXAHEDRON;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone, "Hexahedral Elements", HEXA_8,
                        elem_start, elem_end, 0, (cgsize_t *)Conn_Hexa, &cgns_section);
                }
                if (d_numGlobalPyras > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalPyras; N = (int)d_numGlobalPyras*N_POINTS_PYRAMID;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone, "Pyramid Elements", PYRA_5,
                        elem_start, elem_end, 0, (cgsize_t *)Conn_Pyra, &cgns_section);
                }
                if (d_numGlobalPriss > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalPriss; N = (int)d_numGlobalPriss*N_POINTS_PRISM;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone, "Prism Elements", PENTA_6,
                        elem_start, elem_end, 0, (cgsize_t *)Conn_Pris, &cgns_section);
                }
                if (d_numGlobalBoundLines > 0) {
                    elem_start = 1; elem_end = (int)d_numGlobalBoundLines; N = (int)d_numGlobalBoundLines*N_POINTS_LINE;
                    cgns_err = cg_section_write(cgns_file, cgns_base, cgns_zone, "Line Elements", BAR_2,
                        elem_start, elem_end, 0, (cgsize_t *)Conn_Line, &cgns_section);
                }
                if (cgns_err) cg_error_print();


                cgns_err = cg_close(cgns_file);
                if (cgns_err) cg_error_print();

            }

#else // Not built with CGNS support

            cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n";

#endif

        }

        void Cae_Output::SetCGNS_Solution(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short iZone) {

#ifdef HAVE_CGNS

            /*--- local CGNS variables ---*/
            int cgns_file, cgns_flow, cgns_field, element_dims, physical_dims, cgns_err;
            unsigned long jVar, iVar, iExtIter = config->GetExtIter();
            string base_file, buffer, elements_name;
            stringstream name, results_file;
            bool unsteady = config->GetUnsteady_Simulation();
            cgsize_t isize[3][1];

            bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);

            /*--- Create CGNS base file name ---*/
            base_file = config->GetFlow_FileName();

            /*--- Add CGNS extension. ---*/
            base_file = base_file.append(".cgns");

            /*--- Create CGNS results file name ---*/
            if (unsteady) {

                buffer = config->GetFlow_FileName();

                results_file.str(string()); results_file << buffer;
                if (((int)iExtIter >= 0) && ((int)iExtIter < 10))			results_file << "_0000" << iExtIter;
                if (((int)iExtIter >= 10) && ((int)iExtIter < 100))		results_file << "_000" << iExtIter;
                if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))		results_file << "_00" << iExtIter;
                if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))	results_file << "_0" << iExtIter;
                if ((int)iExtIter >= 10000)							results_file << iExtIter;
                results_file << ".cgns";
            }

            isize[0][0] = (cgsize_t)d_numGlobalPoints;				// vertex size
            isize[1][0] = (cgsize_t)nGlobal_Elem;				// cell size
            isize[2][0] = 0;						// boundary vertex size (zero if elements not sorted)


            if (!unsteady) {

                /*--- Write base file ---*/
                cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_MODIFY, &cgns_file);
                if (cgns_err) cg_error_print();

                element_dims = geometry->GetnDim();		// Currently (release 2.0) only all-2D or all-3D zones permitted
                physical_dims = element_dims;

                /*--- write CGNS descriptor data ---*/
                cgns_err = cg_goto(cgns_file, cgns_base, "end");
                if (cgns_err) cg_error_print();

                /*--- Create a CGNS solution node ---*/
                cgns_err = cg_sol_write(cgns_file, cgns_base, cgns_zone, "Solution", Vertex, &cgns_flow);
                if (cgns_err) cg_error_print();

                cgns_err = cg_goto(cgns_file, cgns_base, "Zone_t", cgns_zone, "FlowSolution_t", cgns_flow, "end");
                if (cgns_err) cg_error_print();

                cgns_err = cg_gridlocation_write(Vertex);
                if (cgns_err) cg_error_print();
            }


            //d_isCgnsOutput = true;
            else {

                /*--- Set up results file for this time step if necessary ---*/

                cgns_err = cg_open((char *)results_file.str().c_str(), CG_MODE_MODIFY, &cgns_file);

                element_dims = geometry->GetnDim();		// Currently (release 2.0) only all-2D or all-3D zones permitted
                physical_dims = element_dims;

                //		/*--- write CGNS base data (one base assumed as of version 3.2.9 "eagle") ---*/
                //		cgns_err = cg_base_write(cgns_file,"SU2 Base", element_dims, physical_dims, &cgns_base);
                //		if (cgns_err) cg_error_print();

                isize[0][0] = (cgsize_t)d_numGlobalPoints;				// vertex size
                isize[1][0] = (cgsize_t)nGlobal_Elem;				// cell size
                isize[2][0] = 0;						// boundary vertex size (zero if elements not sorted)

                //		/*--- write CGNS zone data ---*/
                //		cgns_err = cg_zone_write(cgns_file, cgns_base,"SU2 Zone", isize[0],Unstructured, &cgns_zone);
                //		if (cgns_err) cg_error_print();

                cgns_err = cg_goto(cgns_file, cgns_base_results, "Zone_t", cgns_zone_results, "end");
                if (cgns_err) cg_error_print();

                /*--- Write a CGNS solution node for this time step ---*/
                cgns_err = cg_sol_write(cgns_file, cgns_base_results, cgns_zone_results, "Solution", Vertex, &cgns_flow);
                if (cgns_err) cg_error_print();

                cgns_err = cg_goto(cgns_file, cgns_base_results, "Zone_t", cgns_zone_results, "FlowSolution_t", cgns_flow, "end");
                if (cgns_err) cg_error_print();

                cgns_err = cg_gridlocation_write(Vertex);
                if (cgns_err) cg_error_print();

                cgns_base = cgns_base_results;
                cgns_zone = cgns_zone_results;
            }
            //	else {
            //    
            //		/*--- Open CGNS file for soltuion writing ---*/
            //		cgns_err = cg_open((char *)base_file.c_str(), CG_MODE_MODIFY, &cgns_file);
            //		cgns_base = 1; cgns_zone = 1; cgns_flow = 1;	// fix for multiple zones
            //    
            //	}

            /*	Reference Note on solution variables:
            index 0 --> (nVar_Consv-1)			= Conservative Variables
            nVar_Consv --> (2*nVar_Consv-1)		= Conservative Variable Residuals
            (2*nVar_Consv-1)+					= Additional p, M, T, laminar, eddy depending on solver used */

            /*--- Write conservative variables to CGNS file ---*/
            for (iVar = 0; iVar < nVar_Consv; iVar++) {
                name.str(string()); name << "Conservative Variable " << iVar + 1;
                cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, (char *)name.str().c_str(), Data[iVar], &cgns_field);
                if (cgns_err) cg_error_print();
            }

            /*--- Write primitive variable residuals to CGNS file ---*/
            if (config->GetWrt_Limiters()) {
                for (jVar = 0; jVar < nVar_Consv; jVar++) {
                    name.str(string()); name << "Primitive Limiter " << jVar + 1;
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, (char *)name.str().c_str(), Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                }
            }

            /*--- Write conservative variable residuals to CGNS file ---*/
            if (config->GetWrt_Residuals()) {
                for (jVar = 0; jVar < nVar_Consv; jVar++) {
                    name.str(string()); name << "Conservative Residual " << jVar + 1;
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, (char *)name.str().c_str(), Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                }
            }

            /*--- Write grid velocities to CGNS file, if applicable ---*/
            if (config->GetGrid_Movement()) {
                cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Grid Velocity X", Data[iVar], &cgns_field); iVar++;
                if (cgns_err) cg_error_print();
                cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Grid Velocity Y", Data[iVar], &cgns_field); iVar++;
                if (cgns_err) cg_error_print();
                if (geometry->GetnDim() == 3) {
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Grid Velocity Z", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                }
            }

            if (compressible) {
                switch (config->GetKind_Solver()) {

                    /*--- Write pressure and Mach data to CGNS file ---*/
                case EULER:
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Pressure", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Mach", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    break;

                    /*--- Write temperature and laminar viscosity to CGNS file, if applicable ---*/
                case NAVIER_STOKES:
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Pressure", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Mach", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Temperature", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Viscosity", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    break;

                    /*--- Write eddy viscosity to CGNS file, if applicable ---*/
                case RANS:
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Pressure", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Mach", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Temperature", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Viscosity", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    cgns_err = cg_field_write(cgns_file, cgns_base, cgns_zone, cgns_flow, RealDouble, "Eddy Viscosity", Data[iVar], &cgns_field); iVar++;
                    if (cgns_err) cg_error_print();
                    break;

                default:
                    cout << "Error: Unrecognized equation type \n";
                    exit(EXIT_FAILURE); break;
                }
            }

            /*--- Close CGNS file ---*/
            cgns_err = cg_close(cgns_file);
            if (cgns_err) cg_error_print();

#else // Not built with CGNS support

            cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n";

#endif

        }

        void Cae_Output::SetFieldViewASCII(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone) {

            unsigned short iDim, iVar, nDim = geometry->GetnDim(), ngrids = 1, nbvars, nvars;
            unsigned short Kind_Solver = config->GetKind_Solver();

            unsigned long iPoint, iElem, iNode, nbfaces;

            unsigned long iExtIter = config->GetExtIter();
            bool adjoint = config->GetAdjoint();
            bool grid_movement = config->GetGrid_Movement();

            char cstr[200], buffer[50];
            string filename, FieldName;

            /*--- Write file name with extension ---*/

            if (adjoint) filename = config->GetAdj_FileName();
            else filename = config->GetFlow_FileName();

            if (Kind_Solver == LINEAR_ELASTICITY)
                filename = config->GetStructure_FileName().c_str();

            if (Kind_Solver == WAVE_EQUATION)
                filename = config->GetWave_FileName().c_str();

            if (Kind_Solver == HEAT_EQUATION)
                filename = config->GetHeat_FileName().c_str();

            if (Kind_Solver == POISSON_EQUATION)
                filename = config->GetStructure_FileName().c_str();

            strcpy(cstr, filename.c_str());

            /*--- Special cases where a number needs to be appended to the file name. ---*/

            if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS ||
                Kind_Solver == ADJ_EULER || Kind_Solver == ADJ_NAVIER_STOKES || Kind_Solver == ADJ_RANS) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {

                if (config->GetKind_SU2() == SU2_SOL) { val_iZone = iExtIter; }

                if (int(val_iZone) < 10) sprintf(buffer, "_0000%d.uns", int(val_iZone));
                if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf(buffer, "_000%d.uns", int(val_iZone));
                if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf(buffer, "_00%d.uns", int(val_iZone));
                if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf(buffer, "_0%d.uns", int(val_iZone));
                if (int(val_iZone) >= 10000) sprintf(buffer, "_%d.uns", int(val_iZone));

            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
                if (int(iExtIter) < 10) sprintf(buffer, "_0000%d.uns", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf(buffer, "_000%d.uns", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf(buffer, "_00%d.uns", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.uns", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.uns", int(iExtIter));
            }
            else { sprintf(buffer, ".uns"); }

            strcat(cstr, buffer);

            /*--- Open FieldView ASCII file and write the header ---*/

            ofstream FieldView_File;
            FieldView_File.open(cstr, ios::out);
            FieldView_File.precision(6);

            FieldView_File << "FIELDVIEW 3 0" << endl;

            /*--- Output constants for time, fsmach, alpha and re. ---*/

            FieldView_File << "Constants" << endl;
            FieldView_File << config->GetExtIter() << "\t" << config->GetMach() << "\t" << config->GetAoA() << "\t" << config->GetReynolds() << endl;

            /*--- Output the number of grids. ---*/

            FieldView_File << "Grids\t" << ngrids << endl;

            /*--- Output the table of boundary types, starting with the number of types.
            Note that this differs from the binary/unformatted specification.
            Each boundary type name is preceded by 3 integer flags.
            The first flag indicates whether this boundary type is a wall.
            A flag value of 1 indicates a wall, and a value of 0 indicates
            a non-wall.  Walls are significant for streamline calculation.
            The second flag indicates whether the boundary type has surface
            results.  A value of 1 means surface results will be present for
            this boundary type (if any boundary variables are specified in the
            Boundary Variable Names section below).  A value of 0 means no surface
            results will be present.
            The third flag indicates whether boundary faces of this type have
            consistent "clockness" for the purpose of calculating a surface
            normal.  A value of 1 means that all faces of this type are
            written following a "right hand rule" for clockness.  In other
            words, if the vertices are written on counter-clockwise:
            4 --- 3
            |     |
            1 --- 2
            then the normal to the face is pointing towards you (not away
            from you).  A value of 0 means that the faces do not have any
            consistent clockness.  The "clockness" of surface normals is
            only used for calculating certain special surface integrals
            that involve surface normals.  If the surface normals flag
            is 0, these special integrals will not be available. ---*/

            FieldView_File << "Boundary Table\t1" << endl;
            FieldView_File << "1\t0\t1\tMARKER_PLOTTING" << endl;

            /*--- Output the table of variable names, starting with the number of
            variables.  The number of variables can be zero.
            Note that vector variables are specified by a ';' and vector name
            following the first scalar name of 3 scalar components of the
            vector.  If writing 2-D results, the third component must still
            be provided here, and its values must be written in the variables
            section below (typically padded with zeros.) ---*/

            if (config->GetKind_SU2() == SU2_SOL) {

                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class. ---*/

                nvars = nVar_Total - nDim;

                FieldView_File << "Variable Names\t" << nvars << endl;

                for (unsigned short iField = 1 + nDim; iField < config->fields.size(); iField++) {

                    /*--- Remove all double-quote characters ---*/

                    FieldName = config->fields[iField];

                    FieldName.erase(
                        remove(FieldName.begin(), FieldName.end(), '\"'),
                        FieldName.end()
                        );

                    FieldView_File << FieldName << endl;
                }

                /*--- SU2 does not generate boundary variables ---*/

                nbvars = 0;
                FieldView_File << "Boundary Variable Names\t" << nbvars << endl;

            }

            else {

                nvars = nVar_Total;

                FieldView_File << "Variable Names\t" << nvars << endl;

                for (iVar = 0; iVar < nVar_Consv; iVar++) {
                    FieldView_File << "Conservative_" << iVar + 1 << endl;
                }

                /*--- Add names for any extra variables (this will need to be adjusted). ---*/

                if (config->GetWrt_Limiters()) {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) {
                        FieldView_File << "Limiter_" << iVar + 1 << endl;
                    }
                }

                if (config->GetWrt_Residuals()) {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) {
                        FieldView_File << "Residual_" << iVar + 1 << endl;
                    }
                }

                if (grid_movement) {
                    if (nDim == 2) FieldView_File << "Grid_Velx\nGrid_Vely" << endl;
                    else FieldView_File << "Grid_Velx\nGrid_Vely\nGrid_Velz" << endl;
                }

                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
                    FieldView_File << "Pressure\nTemperature\nPressure_Coefficient\nMach" << endl;
                }

                if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
                    FieldView_File << "Laminar_Viscosity\nSkin_Friction_Coefficient\nHeat_Flux\nY_Plus" << endl;
                }

                if (Kind_Solver == RANS) {
                    FieldView_File << "Eddy_Viscosity" << endl;
                }

                if (config->GetWrt_SharpEdges()) {
                    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
                        FieldView_File << "Sharp_Edge_Dist" << endl;
                    }
                }

                if ((Kind_Solver == ADJ_EULER) ||
                    (Kind_Solver == ADJ_NAVIER_STOKES) ||
                    (Kind_Solver == ADJ_RANS)) {
                    FieldView_File << "Surface_Sensitivity\nSolution_Sensor" << endl;
                }

                /*--- SU2 does not generate boundary variables ---*/

                nbvars = 0;
                FieldView_File << "Boundary Variable Names\t" << nbvars << endl;

            }

            /*--- Output the node definition section for this grid
            Output the X, Y, Z coordinates of successive nodes.
            Note that this differs from the binary/unformatted specification. ---*/

            if (nDim == 3) {

                FieldView_File << "Nodes\t" << d_numGlobalPoints << endl;

                for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                    if (config->GetKind_SU2() != SU2_SOL) {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FieldView_File << scientific << Coords[iDim][iPoint] << "\t";
                    }
                    else {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FieldView_File << scientific << Data[iDim][iPoint] << "\t";
                    }
                    FieldView_File << endl;
                }

            }

            else {

                FieldView_File << "Nodes\t" << d_numGlobalPoints * 2 << endl;

                for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                    if (config->GetKind_SU2() != SU2_SOL) {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FieldView_File << scientific << Coords[iDim][iPoint] << "\t";
                    }
                    else {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FieldView_File << scientific << Data[iDim][iPoint] << "\t";
                    }
                    FieldView_File << scientific << "0.0" << endl;
                }
                for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                    if (config->GetKind_SU2() != SU2_SOL) {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FieldView_File << scientific << Coords[iDim][iPoint] << "\t";
                    }
                    else {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FieldView_File << scientific << Data[iDim][iPoint] << "\t";
                    }
                    FieldView_File << scientific << "-1E-10" << endl;
                }

            }

            /*--- Output the boundary face definitions.
            Note that this differs from the binary/unformatted specification.
            Each face is preceded by an index into the boundary table at the
            top of the file and the number of face vertices, 3 or 4.
            All faces here have 4 vertices.  If the face is triangular,
            the last vertex should be zero.
            TIP: FIELDVIEW assumes that boundary faces are not in random
            order.  It assumes that faces of the same type tend to occur
            in groups.  If your boundary faces are in random order, you
            may want to output them one boundary type at a time.  This
            will give you better performance (less memory, greater speed)
            in FIELDVIEW. ---*/


            if (nDim == 2) {

                nbfaces = d_numGlobalTrias + d_numGlobalQuads;

                FieldView_File << "Boundary Faces\t" << nbfaces << endl;

                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    FieldView_File << "1\t3\t" << Conn_Tria[iNode + 0] << "\t";
                    FieldView_File << Conn_Tria[iNode + 1] << "\t";
                    FieldView_File << Conn_Tria[iNode + 2] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    FieldView_File << "1\t4\t" << Conn_Quad[iNode + 0] << "\t";
                    FieldView_File << Conn_Quad[iNode + 1] << "\t";
                    FieldView_File << Conn_Quad[iNode + 2] << "\t";
                    FieldView_File << Conn_Quad[iNode + 3] << "\n";
                }

            }

            if (nDim == 3) {

                nbfaces = d_numGlobalBoundTrias + d_numGlobalBoundQuads;

                FieldView_File << "Boundary Faces\t" << nbfaces << endl;

                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    FieldView_File << "1\t3\t" << Conn_BoundTria[iNode + 0] << "\t";
                    FieldView_File << Conn_BoundTria[iNode + 1] << "\t";
                    FieldView_File << Conn_BoundTria[iNode + 2] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    FieldView_File << "1\t4\t" << Conn_BoundQuad[iNode + 0] << "\t";
                    FieldView_File << Conn_BoundQuad[iNode + 1] << "\t";
                    FieldView_File << Conn_BoundQuad[iNode + 2] << "\t";
                    FieldView_File << Conn_BoundQuad[iNode + 3] << "\n";
                }

            }


            /*--- Output the elements section for this grid.
            Note that this differs from the binary/unformatted specification.
            It contains the headers and node definitions of all elements.
            In this example, each element starts with 2 for type 'hex',
            with a subtype of 1 (the only subtype currently supported).
            This is followed by the node indices for the element. ---*/


            FieldView_File << "Elements" << endl;

            for (iElem = 0; iElem < d_numGlobalTrias; iElem++) {
                iNode = iElem*N_POINTS_TRIANGLE;
                FieldView_File << "3\t1\t" << Conn_Tria[iNode + 0] << "\t";
                FieldView_File << Conn_Tria[iNode + 1] << "\t";
                FieldView_File << Conn_Tria[iNode + 2] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Tria[iNode + 0] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Tria[iNode + 1] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Tria[iNode + 2] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalQuads; iElem++) {
                iNode = iElem*N_POINTS_QUADRILATERAL;
                FieldView_File << "2\t1\t" << Conn_Quad[iNode + 0] << "\t";
                FieldView_File << Conn_Quad[iNode + 1] << "\t";
                FieldView_File << Conn_Quad[iNode + 2] << "\t";
                FieldView_File << Conn_Quad[iNode + 3] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Quad[iNode + 0] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Quad[iNode + 1] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Quad[iNode + 2] << "\t";
                FieldView_File << d_numGlobalPoints + Conn_Quad[iNode + 3] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) {
                iNode = iElem*N_POINTS_TETRAHEDRON;
                FieldView_File << "1\t1\t" << Conn_Tetr[iNode + 0] << "\t" << Conn_Tetr[iNode + 1] << "\t";
                FieldView_File << Conn_Tetr[iNode + 2] << "\t" << Conn_Tetr[iNode + 3] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalHexas; iElem++) {
                iNode = iElem*N_POINTS_HEXAHEDRON;
                FieldView_File << "2\t1\t" << Conn_Hexa[iNode + 0] << "\t" << Conn_Hexa[iNode + 1] << "\t";
                FieldView_File << Conn_Hexa[iNode + 2] << "\t" << Conn_Hexa[iNode + 3] << "\t";
                FieldView_File << Conn_Hexa[iNode + 4] << "\t" << Conn_Hexa[iNode + 5] << "\t";
                FieldView_File << Conn_Hexa[iNode + 6] << "\t" << Conn_Hexa[iNode + 7] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalPriss; iElem++) {
                iNode = iElem*N_POINTS_PRISM;
                FieldView_File << "3\t1\t" << Conn_Pris[iNode + 0] << "\t" << Conn_Pris[iNode + 1] << "\t";
                FieldView_File << Conn_Pris[iNode + 2] << "\t" << Conn_Pris[iNode + 3] << "\t";
                FieldView_File << Conn_Pris[iNode + 4] << "\t" << Conn_Pris[iNode + 5] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalPyras; iElem++) {
                iNode = iElem*N_POINTS_PYRAMID;
                FieldView_File << "4\t1\t" << Conn_Pyra[iNode + 0] << "\t" << Conn_Pyra[iNode + 1] << "\t";
                FieldView_File << Conn_Pyra[iNode + 2] << "\t" << Conn_Pyra[iNode + 3] << "\t";
                FieldView_File << Conn_Pyra[iNode + 4] << "\n";
            }

            /*--- Output the variables data for this grid.
            Note that all of the data for the first variable is output
            before any of the data for the second variable.
            You should skip this section if the number of variables is zero.
            The variables must be in the same order as the "Variable Names"
            section. ---*/

            FieldView_File << "Variables" << endl;

            /*--- Loop over the vars/residuals and write the values to file ---*/

            if (config->GetKind_SU2() != SU2_SOL) {
                for (iVar = 0; iVar < nvars; iVar++) {
                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        FieldView_File << scientific << Data[iVar][iPoint] << endl;
                    }
                    if (nDim == 2) {
                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            FieldView_File << scientific << Data[iVar][iPoint] << endl;
                        }
                    }
                }
            }
            else {
                for (iVar = 0; iVar < nvars; iVar++) {
                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        FieldView_File << scientific << Data[iVar + nDim][iPoint] << endl;
                    }
                    if (nDim == 2) {
                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            FieldView_File << scientific << Data[iVar + nDim][iPoint] << endl;
                        }
                    }
                }
            }

            /*--- Output the boundary variables data for this grid.
            Note that all of the data for the first variable is output
            before any of the data for the second variable.
            Remember that the Boundary Table above has a "surface results
            flag" indicating which boundary types have surface results.
            The data should be written in the same order as the faces in
            the Boundary Faces section, skipping over faces whose boundary
            type has a surface results flag of zero (false).
            For each variable, you should write one number per boundary face.
            You should skip this section if the number of boundary
            variables is zero. ---*/

            FieldView_File << "Boundary Variables" << endl;


            FieldView_File.close();

        }

        void Cae_Output::SetFieldViewASCII_Mesh(Config_p config, GEOM::GEOM_Geometry *geometry) { }

        void Cae_Output::SetFieldViewBinary(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone) {

            unsigned short iDim, iVar, nDim = geometry->GetnDim(), ngrids = 1, nbvars, nvars;
            unsigned short Kind_Solver = config->GetKind_Solver();

            unsigned long iPoint, iElem, iNode, nbfaces;
            unsigned long iExtIter = config->GetExtIter();
            bool adjoint = config->GetAdjoint();

            char cstr[200], buffer[50];
            string filename;

            /*--- Write file name with extension ---*/

            if (adjoint) filename = config->GetAdj_FileName();
            else filename = config->GetFlow_FileName();

            if (Kind_Solver == LINEAR_ELASTICITY)
                filename = config->GetStructure_FileName().c_str();

            if (Kind_Solver == WAVE_EQUATION)
                filename = config->GetWave_FileName().c_str();

            if (Kind_Solver == HEAT_EQUATION)
                filename = config->GetHeat_FileName().c_str();

            if (Kind_Solver == POISSON_EQUATION)
                filename = config->GetStructure_FileName().c_str();

            strcpy(cstr, filename.c_str());

            /*--- Special cases where a number needs to be appended to the file name. ---*/

            if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS ||
                Kind_Solver == ADJ_EULER || Kind_Solver == ADJ_NAVIER_STOKES || Kind_Solver == ADJ_RANS) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {

                if (config->GetKind_SU2() == SU2_SOL) { val_iZone = iExtIter; }

                if (int(val_iZone) < 10) sprintf(buffer, "_0000%d.uns", int(val_iZone));
                if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf(buffer, "_000%d.uns", int(val_iZone));
                if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf(buffer, "_00%d.uns", int(val_iZone));
                if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf(buffer, "_0%d.uns", int(val_iZone));
                if (int(val_iZone) >= 10000) sprintf(buffer, "_%d.uns", int(val_iZone));

            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
                if (int(iExtIter) < 10) sprintf(buffer, "_0000%d.uns", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf(buffer, "_000%d.uns", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf(buffer, "_00%d.uns", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.uns", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.uns", int(iExtIter));
            }
            else { sprintf(buffer, ".uns"); }

            strcat(cstr, buffer);

            /*--- Open FieldView ASCII file and write the header ---*/

            ofstream FieldView_File;
            FieldView_File.open(cstr, ios::out);
            FieldView_File.precision(6);

            FieldView_File << "FIELDVIEW 3 0" << endl;

            /*--- Output constants for time, fsmach, alpha and re. ---*/

            FieldView_File << "Constants" << endl;
            FieldView_File << config->GetExtIter() << "\t" << config->GetMach() << "\t" << config->GetAoA() << "\t" << config->GetReynolds() << endl;

            /*--- Output the number of grids. ---*/

            FieldView_File << "Grids\t" << ngrids << endl;

            /*--- Output the table of boundary types, starting with the number of types.
            Note that this differs from the binary/unformatted specification.
            Each boundary type name is preceded by 3 integer flags.
            The first flag indicates whether this boundary type is a wall.
            A flag value of 1 indicates a wall, and a value of 0 indicates
            a non-wall.  Walls are significant for streamline calculation.
            The second flag indicates whether the boundary type has surface
            results.  A value of 1 means surface results will be present for
            this boundary type (if any boundary variables are specified in the
            Boundary Variable Names section below).  A value of 0 means no surface
            results will be present.
            The third flag indicates whether boundary faces of this type have
            consistent "clockness" for the purpose of calculating a surface
            normal.  A value of 1 means that all faces of this type are
            written following a "right hand rule" for clockness.  In other
            words, if the vertices are written on counter-clockwise:
            4 --- 3
            |     |
            1 --- 2
            then the normal to the face is pointing towards you (not away
            from you).  A value of 0 means that the faces do not have any
            consistent clockness.  The "clockness" of surface normals is
            only used for calculating certain special surface integrals
            that involve surface normals.  If the surface normals flag
            is 0, these special integrals will not be available. ---*/

            FieldView_File << "Boundary Table\t1" << endl;
            FieldView_File << "1\t0\t1\tMARKER_PLOTTING" << endl;

            /*--- Output the table of variable names, starting with the number of
            variables.  The number of variables can be zero.
            Note that vector variables are specified by a ';' and vector name
            following the first scalar name of 3 scalar components of the
            vector.  If writing 2-D results, the third component must still
            be provided here, and its values must be written in the variables
            section below (typically padded with zeros.) ---*/

            if (config->GetKind_SU2() == SU2_SOL) {

                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class. ---*/

                nvars = config->fields.size() - 1 - nDim;

                FieldView_File << "Variable Names\t" << nvars << endl;

                for (unsigned short iField = 1 + nDim; iField < config->fields.size(); iField++) {
                    FieldView_File << config->fields[iField] << endl;
                }

                /*--- SU2 does not generate boundary variables ---*/

                nbvars = 0;
                FieldView_File << "Boundary Variable Names\t" << nbvars << endl;

            }

            else {

                nvars = nVar_Consv + nDim;

                FieldView_File << "Variable Names\t" << nvars - nDim << endl;

                for (iVar = 0; iVar < nVar_Consv; iVar++) {
                    FieldView_File << "Conservative_" << iVar + 1 << endl;
                }

                /*--- SU2 does not generate boundary variables ---*/

                nbvars = 0;
                FieldView_File << "Boundary Variable Names\t" << nbvars << endl;

            }

            /*--- Output the node definition section for this grid ---*/

            FieldView_File << "Nodes\t" << d_numGlobalPoints << endl;

            /*--- Output the X, Y, Z coordinates of successive nodes.
            Note that this differs from the binary/unformatted specification. ---*/

            for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {

                if (config->GetKind_SU2() != SU2_SOL) {
                    for (iDim = 0; iDim < nDim; iDim++)
                        FieldView_File << scientific << Coords[iDim][iPoint] << "\t";
                }
                else {
                    for (iVar = 0; iVar < nVar_Total; iVar++)
                        FieldView_File << scientific << Data[iVar][iPoint] << "\t";
                }
                FieldView_File << endl;

            }

            /*--- Output the boundary face definitions.
            Note that this differs from the binary/unformatted specification.
            Each face is preceded by an index into the boundary table at the
            top of the file and the number of face vertices, 3 or 4.
            All faces here have 4 vertices.  If the face is triangular,
            the last vertex should be zero.
            TIP: FIELDVIEW assumes that boundary faces are not in random
            order.  It assumes that faces of the same type tend to occur
            in groups.  If your boundary faces are in random order, you
            may want to output them one boundary type at a time.  This
            will give you better performance (less memory, greater speed)
            in FIELDVIEW. ---*/

            nbfaces = d_numGlobalBoundLines + d_numGlobalBoundTrias + d_numGlobalBoundQuads;

            FieldView_File << "Boundary Faces\t" << nbfaces << endl;

            for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                iNode = iElem*N_POINTS_LINE;
                FieldView_File << "1\t2\t" << Conn_Line[iNode + 0] << "\t";
                FieldView_File << "1\t2\t" << Conn_Line[iNode + 1] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                iNode = iElem*N_POINTS_TRIANGLE;
                FieldView_File << "1\t3\t" << Conn_BoundTria[iNode + 0] << "\t";
                FieldView_File << Conn_BoundTria[iNode + 1] << "\t";
                FieldView_File << Conn_BoundTria[iNode + 2] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                iNode = iElem*N_POINTS_QUADRILATERAL;
                FieldView_File << "1\t4\t" << Conn_BoundQuad[iNode + 0] << "\t";
                FieldView_File << Conn_BoundQuad[iNode + 1] << "\t";
                FieldView_File << Conn_BoundQuad[iNode + 2] << "\t";
                FieldView_File << Conn_BoundQuad[iNode + 3] << "\n";
            }

            /*--- Output the elements section for this grid.
            Note that this differs from the binary/unformatted specification.
            It contains the headers and node definitions of all elements.
            In this example, each element starts with 2 for type 'hex',
            with a subtype of 1 (the only subtype currently supported).
            This is followed by the node indices for the element. ---*/


            FieldView_File << "Elements" << endl;

            for (iElem = 0; iElem < d_numGlobalTrias; iElem++) {
                iNode = iElem*N_POINTS_TRIANGLE;
                FieldView_File << "2\t1\t" << Conn_Tria[iNode + 0] << "\t";
                FieldView_File << Conn_Tria[iNode + 1] << "\t";
                FieldView_File << Conn_Tria[iNode + 2] << "\t";
                FieldView_File << Conn_Tria[iNode + 2] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalQuads; iElem++) {
                iNode = iElem*N_POINTS_QUADRILATERAL;
                FieldView_File << "2\t1\t" << Conn_Quad[iNode + 0] << "\t";
                FieldView_File << Conn_Quad[iNode + 1] << "\t";
                FieldView_File << Conn_Quad[iNode + 2] << "\t";
                FieldView_File << Conn_Quad[iNode + 3] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) {
                iNode = iElem*N_POINTS_TETRAHEDRON;
                FieldView_File << "1\t1\t" << Conn_Tetr[iNode + 0] << "\t" << Conn_Tetr[iNode + 1] << "\t";
                FieldView_File << Conn_Tetr[iNode + 2] << "\t" << Conn_Tetr[iNode + 2] << "\t";
                FieldView_File << Conn_Tetr[iNode + 3] << "\t" << Conn_Tetr[iNode + 3] << "\t";
                FieldView_File << Conn_Tetr[iNode + 3] << "\t" << Conn_Tetr[iNode + 3] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalHexas; iElem++) {
                iNode = iElem*N_POINTS_HEXAHEDRON;
                FieldView_File << "2\t1\t" << Conn_Hexa[iNode + 0] << "\t" << Conn_Hexa[iNode + 1] << "\t";
                FieldView_File << Conn_Hexa[iNode + 2] << "\t" << Conn_Hexa[iNode + 3] << "\t";
                FieldView_File << Conn_Hexa[iNode + 4] << "\t" << Conn_Hexa[iNode + 5] << "\t";
                FieldView_File << Conn_Hexa[iNode + 6] << "\t" << Conn_Hexa[iNode + 7] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalPriss; iElem++) {
                iNode = iElem*N_POINTS_PRISM;
                FieldView_File << "3\t1\t" << Conn_Pris[iNode + 0] << "\t" << Conn_Pris[iNode + 1] << "\t";
                FieldView_File << Conn_Pris[iNode + 1] << "\t" << Conn_Pris[iNode + 2] << "\t";
                FieldView_File << Conn_Pris[iNode + 3] << "\t" << Conn_Pris[iNode + 4] << "\t";
                FieldView_File << Conn_Pris[iNode + 4] << "\t" << Conn_Pris[iNode + 5] << "\n";
            }

            for (iElem = 0; iElem < d_numGlobalPyras; iElem++) {
                iNode = iElem*N_POINTS_PYRAMID;
                FieldView_File << "4\t1\t" << Conn_Pyra[iNode + 0] << "\t" << Conn_Pyra[iNode + 1] << "\t";
                FieldView_File << Conn_Pyra[iNode + 2] << "\t" << Conn_Pyra[iNode + 3] << "\t";
                FieldView_File << Conn_Pyra[iNode + 4] << "\t" << Conn_Pyra[iNode + 4] << "\t";
                FieldView_File << Conn_Pyra[iNode + 4] << "\t" << Conn_Pyra[iNode + 4] << "\n";
            }

            /*--- Output the variables data for this grid.
            Note that all of the data for the first variable is output
            before any of the data for the second variable.
            You should skip this section if the number of variables is zero.
            The variables must be in the same order as the "Variable Names"
            section. ---*/

            FieldView_File << "Variables" << endl;

            /*--- Loop over the vars/residuals and write the values to file ---*/
            for (iVar = nDim; iVar < nVar_Total; iVar++) {
                for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                    FieldView_File << scientific << Data[iVar][iPoint] << endl;
                }
            }

            /*--- Output the boundary variables data for this grid.
            Note that all of the data for the first variable is output
            before any of the data for the second variable.
            Remember that the Boundary Table above has a "surface results
            flag" indicating which boundary types have surface results.
            The data should be written in the same order as the faces in
            the Boundary Faces section, skipping over faces whose boundary
            type has a surface results flag of zero (false).
            For each variable, you should write one number per boundary face.
            You should skip this section if the number of boundary
            variables is zero. ---*/

            FieldView_File << "Boundary Variables" << endl;


            FieldView_File.close();

        }

        void Cae_Output::SetFieldViewBinary_Mesh(Config_p config, GEOM::GEOM_Geometry *geometry) { }

        void Cae_Output::SetParaview_ASCII(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {

            /*--- Local variables and initialization ---*/
            unsigned short iDim, iVar, nDim = geometry->GetnDim();
            unsigned short Kind_Solver = config->GetKind_Solver();

            unsigned long iPoint, iElem, iNode;
            unsigned long iExtIter = config->GetExtIter();
            unsigned long *LocalIndex = NULL;
            bool *SurfacePoint = NULL;

            unsigned long d_numSurfElems_Storage;
            unsigned long nGlobal_Elem_Storage;

            bool grid_movement = config->GetGrid_Movement();
            bool adjoint = config->GetAdjoint();

            char cstr[200], buffer[50];
            string filename;

            /*--- Write file name with extension ---*/
            if (surf_sol) {
                if (adjoint)
                    filename = config->GetSurfAdjCoeff_FileName();
                else
                    filename = config->GetSurfFlowCoeff_FileName();
            }
            else {
                if (adjoint)
                    filename = config->GetAdj_FileName();
                else
                    filename = config->GetFlow_FileName();
            }

            if (Kind_Solver == LINEAR_ELASTICITY)
                filename = config->GetStructure_FileName().c_str();

            if (Kind_Solver == WAVE_EQUATION)
                filename = config->GetWave_FileName().c_str();

            if (Kind_Solver == POISSON_EQUATION)
                filename = config->GetStructure_FileName().c_str();

            if (Kind_Solver == HEAT_EQUATION)
                filename = config->GetHeat_FileName().c_str();

            strcpy(cstr, filename.c_str());
            if (Kind_Solver == POISSON_EQUATION) strcpy(cstr, config->GetStructure_FileName().c_str());

            /*--- Special cases where a number needs to be appended to the file name. ---*/
            if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            /*--- Special cases where a number needs to be appended to the file name. ---*/
            if (((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
                if (int(val_iZone) < 10) sprintf(buffer, "_0000%d.vtk", int(val_iZone));
                if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf(buffer, "_000%d.vtk", int(val_iZone));
                if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf(buffer, "_00%d.vtk", int(val_iZone));
                if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf(buffer, "_0%d.vtk", int(val_iZone));
                if (int(val_iZone) >= 10000) sprintf(buffer, "_%d.vtk", int(val_iZone));

            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
                if (int(iExtIter) < 10) sprintf(buffer, "_0000%d.vtk", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf(buffer, "_000%d.vtk", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf(buffer, "_00%d.vtk", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.vtk", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.vtk", int(iExtIter));
            }
            else {
                sprintf(buffer, ".vtk");
            }

            strcat(cstr, buffer);

            /*--- Open Paraview ASCII file and write the header. ---*/
            ofstream Paraview_File;
            Paraview_File.open(cstr, ios::out);
            Paraview_File.precision(6);
            Paraview_File << "# vtk DataFile Version 3.0\n";
            Paraview_File << "vtk output\n";
            Paraview_File << "ASCII\n";
            Paraview_File << "DATASET UNSTRUCTURED_GRID\n";

            /*--- If it's a surface output, print only the points
            that are in the element list, change the numbering ---*/

            if (surf_sol) {

                LocalIndex = new unsigned long[d_numGlobalPoints + 1];
                SurfacePoint = new bool[d_numGlobalPoints + 1];

                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) SurfacePoint[iPoint] = false;

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    SurfacePoint[Conn_Line[iNode + 0]] = true;
                    SurfacePoint[Conn_Line[iNode + 1]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    SurfacePoint[Conn_BoundTria[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 2]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    SurfacePoint[Conn_BoundQuad[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 2]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 3]] = true;
                }

                d_numSurfPoints = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                    LocalIndex[iPoint] = 0;
                    if (SurfacePoint[iPoint]) { d_numSurfPoints++; LocalIndex[iPoint] = d_numSurfPoints; }
                }

            }

            /*--- Write the header ---*/
            if (surf_sol) Paraview_File << "POINTS " << d_numSurfPoints << " float\n";
            else Paraview_File << "POINTS " << d_numGlobalPoints << " float\n";

            /*--- Write surface and volumetric solution data. ---*/
            for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {

                if (surf_sol) {

                    if (LocalIndex[iPoint + 1] != 0) {

                        /*--- Write the node coordinates ---*/
                        if (config->GetKind_SU2() != SU2_SOL) {
                            for (iDim = 0; iDim < nDim; iDim++)
                                Paraview_File << scientific << Coords[iDim][iPoint] << "\t";
                            if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                        }
                        else {
                            for (iDim = 0; iDim < nDim; iDim++)
                                Paraview_File << scientific << Data[iDim][iPoint] << "\t";
                            if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                        }

                    }

                }
                else {

                    if (config->GetKind_SU2() != SU2_SOL) {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Paraview_File << scientific << Coords[iDim][iPoint] << "\t";
                        if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                    }
                    else {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Paraview_File << scientific << Data[iDim][iPoint] << "\t";
                        if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                    }

                }
            }

            /*--- Write the header ---*/
            d_numSurfElems_Storage = d_numGlobalBoundLines * 3 + d_numGlobalBoundTrias * 4 + d_numGlobalBoundQuads * 5;
            nGlobal_Elem_Storage = d_numGlobalTrias * 4 + d_numGlobalQuads * 5 + d_numGlobalTetrs * 5 + d_numGlobalHexas * 9 + d_numGlobalPriss * 7 + d_numGlobalPyras * 6;

            if (surf_sol) Paraview_File << "\nCELLS " << d_numSurfElems << "\t" << d_numSurfElems_Storage << "\n";
            else Paraview_File << "\nCELLS " << nGlobal_Elem << "\t" << nGlobal_Elem_Storage << "\n";

            if (surf_sol) {

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    Paraview_File << N_POINTS_LINE << "\t";
                    Paraview_File << LocalIndex[Conn_Line[iNode + 0]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_Line[iNode + 1]] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Paraview_File << N_POINTS_TRIANGLE << "\t";
                    Paraview_File << LocalIndex[Conn_BoundTria[iNode + 0]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundTria[iNode + 1]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundTria[iNode + 2]] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Paraview_File << N_POINTS_QUADRILATERAL << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 0]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 1]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 2]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 3]] - 1 << "\t";
                }

            }
            else {

                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Paraview_File << N_POINTS_TRIANGLE << "\t";
                    Paraview_File << Conn_Tria[iNode + 0] - 1 << "\t";
                    Paraview_File << Conn_Tria[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Tria[iNode + 2] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Paraview_File << N_POINTS_QUADRILATERAL << "\t";
                    Paraview_File << Conn_Quad[iNode + 0] - 1 << "\t";
                    Paraview_File << Conn_Quad[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Quad[iNode + 2] - 1 << "\t";
                    Paraview_File << Conn_Quad[iNode + 3] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) {
                    iNode = iElem*N_POINTS_TETRAHEDRON;
                    Paraview_File << N_POINTS_TETRAHEDRON << "\t";
                    Paraview_File << Conn_Tetr[iNode + 0] - 1 << "\t" << Conn_Tetr[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Tetr[iNode + 2] - 1 << "\t" << Conn_Tetr[iNode + 3] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalHexas; iElem++) {
                    iNode = iElem*N_POINTS_HEXAHEDRON;
                    Paraview_File << N_POINTS_HEXAHEDRON << "\t";
                    Paraview_File << Conn_Hexa[iNode + 0] - 1 << "\t" << Conn_Hexa[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Hexa[iNode + 2] - 1 << "\t" << Conn_Hexa[iNode + 3] - 1 << "\t";
                    Paraview_File << Conn_Hexa[iNode + 4] - 1 << "\t" << Conn_Hexa[iNode + 5] - 1 << "\t";
                    Paraview_File << Conn_Hexa[iNode + 6] - 1 << "\t" << Conn_Hexa[iNode + 7] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalPriss; iElem++) {
                    iNode = iElem*N_POINTS_PRISM;
                    Paraview_File << N_POINTS_PRISM << "\t";
                    Paraview_File << Conn_Pris[iNode + 0] - 1 << "\t" << Conn_Pris[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Pris[iNode + 2] - 1 << "\t" << Conn_Pris[iNode + 3] - 1 << "\t";
                    Paraview_File << Conn_Pris[iNode + 4] - 1 << "\t" << Conn_Pris[iNode + 5] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalPyras; iElem++) {
                    iNode = iElem*N_POINTS_PYRAMID;
                    Paraview_File << N_POINTS_PYRAMID << "\t";
                    Paraview_File << Conn_Pyra[iNode + 0] - 1 << "\t" << Conn_Pyra[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Pyra[iNode + 2] - 1 << "\t" << Conn_Pyra[iNode + 3] - 1 << "\t";
                    Paraview_File << Conn_Pyra[iNode + 4] - 1 << "\t";
                }
            }

            /*--- Write the header ---*/
            if (surf_sol) Paraview_File << "\nCELL_TYPES " << d_numSurfElems << "\n";
            else Paraview_File << "\nCELL_TYPES " << nGlobal_Elem << "\n";

            if (surf_sol) {
                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) Paraview_File << "3\t";
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) Paraview_File << "5\t";
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) Paraview_File << "9\t";

            }
            else {
                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) Paraview_File << "5\t";
                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) Paraview_File << "9\t";
                for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) Paraview_File << "10\t";
                for (iElem = 0; iElem < d_numGlobalHexas; iElem++) Paraview_File << "12\t";
                for (iElem = 0; iElem < d_numGlobalPriss; iElem++) Paraview_File << "13\t";
                for (iElem = 0; iElem < d_numGlobalPyras; iElem++) Paraview_File << "14\t";
            }



            /*--- Write the header ---*/
            if (surf_sol) Paraview_File << "\nPOINT_DATA " << d_numSurfPoints << "\n";
            else Paraview_File << "\nPOINT_DATA " << d_numGlobalPoints << "\n";

            unsigned short VarCounter = 0;

            if (config->GetKind_SU2() == SU2_SOL) {

                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class. ---*/
                for (unsigned short iField = 1; iField < config->fields.size(); iField++) {

                    double output_variable = true;
                    size_t found = config->fields[iField].find("\"x\"");
                    if (found != string::npos) output_variable = false;
                    found = config->fields[iField].find("\"y\"");
                    if (found != string::npos) output_variable = false;
                    found = config->fields[iField].find("\"z\"");
                    if (found != string::npos) output_variable = false;

                    if (output_variable) {
                        Paraview_File << "\nSCALARS " << config->fields[iField] << " float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                    }

                    VarCounter++;


                }

            }

            else {

                for (iVar = 0; iVar < nVar_Consv; iVar++) {

                    Paraview_File << "\nSCALARS Conservative_" << iVar + 1 << " float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;
                }

                if (config->GetWrt_Limiters()) {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) {

                        Paraview_File << "\nSCALARS Limiter_" << iVar + 1 << " float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        VarCounter++;
                    }
                }

                if (config->GetWrt_Residuals()) {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) {

                        Paraview_File << "\nSCALARS Residual_" << iVar + 1 << " float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        VarCounter++;
                    }
                }

                /*--- Add names for any extra variables (this will need to be adjusted). ---*/
                if (grid_movement) {

                    Paraview_File << "\nSCALARS Grid_Velx float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Grid_Vely float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    if (nDim == 3) {

                        Paraview_File << "\nSCALARS Grid_Velz float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        VarCounter++;

                    }
                }

                if (config->GetKind_Regime() == FREESURFACE) {

                    Paraview_File << "\nSCALARS Density float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {

                    Paraview_File << "\nSCALARS Pressure float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Temperature float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Pressure_Coefficient float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Mach float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {

                    Paraview_File << "\nSCALARS Laminar_Viscosity float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Skin_Friction_Coefficient float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Heat_Flux float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Y_Plus float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if (Kind_Solver == RANS) {

                    Paraview_File << "\nSCALARS Eddy_Viscosity float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if ((Kind_Solver == ADJ_EULER) ||
                    (Kind_Solver == ADJ_NAVIER_STOKES) ||
                    (Kind_Solver == ADJ_RANS)) {

                    Paraview_File << "\nSCALARS Surface_Sensitivity float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

            }

            Paraview_File.close();

            if (surf_sol) delete[] LocalIndex;

        }

        void Cae_Output::SetParaview_MeshASCII(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol, bool new_file) {

            /*--- Local variables and initialization ---*/
            unsigned short iDim, iVar, nDim = geometry->GetnDim();
            unsigned short Kind_Solver = config->GetKind_Solver();

            unsigned long iPoint, iElem, iNode;
            unsigned long iExtIter = config->GetExtIter();
            unsigned long *LocalIndex = NULL;
            bool *SurfacePoint = NULL;

            unsigned long d_numSurfElems_Storage;
            unsigned long nGlobal_Elem_Storage;

            bool grid_movement = config->GetGrid_Movement();
            bool adjoint = config->GetAdjoint();

            char cstr[200], buffer[50];
            string filename;

            /*--- Write file name with extension ---*/
            if (surf_sol) {
                if (adjoint)
                    filename = config->GetSurfAdjCoeff_FileName();
                else
                    filename = config->GetSurfFlowCoeff_FileName();
            }
            else {
                if (adjoint)
                    filename = config->GetAdj_FileName();
                else
                    filename = config->GetFlow_FileName();
            }
            if (config->GetKind_SU2() == SU2_DEF){
                if (new_file){
                    if (surf_sol) filename = "surface_grid";
                    else filename = "volumetric_grid";
                }
                else{
                    if (surf_sol) filename = "surface_deformed_grid";
                    else filename = "volumetric_deformed_grid";
                }
            }

            if (Kind_Solver == LINEAR_ELASTICITY)
                filename = config->GetStructure_FileName().c_str();

            if (Kind_Solver == WAVE_EQUATION)
                filename = config->GetWave_FileName().c_str();

            if (Kind_Solver == POISSON_EQUATION)
                filename = config->GetStructure_FileName().c_str();

            if (Kind_Solver == HEAT_EQUATION)
                filename = config->GetHeat_FileName().c_str();

            strcpy(cstr, filename.c_str());
            if (Kind_Solver == POISSON_EQUATION) strcpy(cstr, config->GetStructure_FileName().c_str());

            /*--- Special cases where a number needs to be appended to the file name. ---*/
            if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            /*--- Special cases where a number needs to be appended to the file name. ---*/
            if (((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
                if (int(val_iZone) < 10) sprintf(buffer, "_0000%d.vtk", int(val_iZone));
                if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf(buffer, "_000%d.vtk", int(val_iZone));
                if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf(buffer, "_00%d.vtk", int(val_iZone));
                if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf(buffer, "_0%d.vtk", int(val_iZone));
                if (int(val_iZone) >= 10000) sprintf(buffer, "_%d.vtk", int(val_iZone));

            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
                if (int(iExtIter) < 10) sprintf(buffer, "_0000%d.vtk", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf(buffer, "_000%d.vtk", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf(buffer, "_00%d.vtk", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.vtk", int(iExtIter));
                if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.vtk", int(iExtIter));
            }
            else {
                sprintf(buffer, ".vtk");
            }

            strcat(cstr, buffer);

            /*--- Open Paraview ASCII file and write the header. ---*/
            ofstream Paraview_File;
            Paraview_File.open(cstr, ios::out);
            Paraview_File.precision(6);
            Paraview_File << "# vtk DataFile Version 3.0\n";
            Paraview_File << "vtk output\n";
            Paraview_File << "ASCII\n";
            if (config->GetKind_SU2() != SU2_DEF) Paraview_File << "DATASET UNSTRUCTURED_GRID\n";
            else Paraview_File << "DATASET UNSTRUCTURED_GRID\n";


            /*--- If it's a surface output, print only the points
            that are in the element list, change the numbering ---*/

            if (surf_sol) {

                LocalIndex = new unsigned long[d_numGlobalPoints + 1];
                SurfacePoint = new bool[d_numGlobalPoints + 1];

                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) SurfacePoint[iPoint] = false;

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    SurfacePoint[Conn_Line[iNode + 0]] = true;
                    SurfacePoint[Conn_Line[iNode + 1]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    SurfacePoint[Conn_BoundTria[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 2]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    SurfacePoint[Conn_BoundQuad[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 2]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 3]] = true;
                }

                d_numSurfPoints = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                    LocalIndex[iPoint] = 0;
                    if (SurfacePoint[iPoint]) { d_numSurfPoints++; LocalIndex[iPoint] = d_numSurfPoints; }
                }

            }

            /*--- Write the header ---*/
            if (surf_sol) Paraview_File << "POINTS " << d_numSurfPoints << " float\n";
            else Paraview_File << "POINTS " << d_numGlobalPoints << " float\n";

            /*--- Write surface and volumetric solution data. ---*/
            for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {

                if (surf_sol) {

                    if (LocalIndex[iPoint + 1] != 0) {

                        /*--- Write the node coordinates ---*/
                        if (config->GetKind_SU2() != SU2_SOL) {
                            for (iDim = 0; iDim < nDim; iDim++)
                                Paraview_File << scientific << Coords[iDim][iPoint] << "\t";
                            if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                        }
                        else {
                            for (iDim = 0; iDim < nDim; iDim++)
                                Paraview_File << scientific << Data[iDim][iPoint] << "\t";
                            if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                        }

                    }

                }
                else {

                    if (config->GetKind_SU2() != SU2_SOL) {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Paraview_File << scientific << Coords[iDim][iPoint] << "\t";
                        if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                    }
                    else {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Paraview_File << scientific << Data[iDim][iPoint] << "\t";
                        if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
                    }

                }
            }

            /*--- Write the header ---*/
            d_numSurfElems_Storage = d_numGlobalBoundLines * 3 + d_numGlobalBoundTrias * 4 + d_numGlobalBoundQuads * 5;
            nGlobal_Elem_Storage = d_numGlobalTrias * 4 + d_numGlobalQuads * 5 + d_numGlobalTetrs * 5 + d_numGlobalHexas * 9 + d_numGlobalPriss * 7 + d_numGlobalPyras * 6;

            if (surf_sol) Paraview_File << "\nCELLS " << d_numSurfElems << "\t" << d_numSurfElems_Storage << "\n";
            else Paraview_File << "\nCELLS " << nGlobal_Elem << "\t" << nGlobal_Elem_Storage << "\n";

            if (surf_sol) {

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    Paraview_File << N_POINTS_LINE << "\t";
                    Paraview_File << LocalIndex[Conn_Line[iNode + 0]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_Line[iNode + 1]] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Paraview_File << N_POINTS_TRIANGLE << "\t";
                    Paraview_File << LocalIndex[Conn_BoundTria[iNode + 0]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundTria[iNode + 1]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundTria[iNode + 2]] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Paraview_File << N_POINTS_QUADRILATERAL << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 0]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 1]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 2]] - 1 << "\t";
                    Paraview_File << LocalIndex[Conn_BoundQuad[iNode + 3]] - 1 << "\t";
                }

            }
            else {

                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Paraview_File << N_POINTS_TRIANGLE << "\t";
                    Paraview_File << Conn_Tria[iNode + 0] - 1 << "\t";
                    Paraview_File << Conn_Tria[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Tria[iNode + 2] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Paraview_File << N_POINTS_QUADRILATERAL << "\t";
                    Paraview_File << Conn_Quad[iNode + 0] - 1 << "\t";
                    Paraview_File << Conn_Quad[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Quad[iNode + 2] - 1 << "\t";
                    Paraview_File << Conn_Quad[iNode + 3] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) {
                    iNode = iElem*N_POINTS_TETRAHEDRON;
                    Paraview_File << N_POINTS_TETRAHEDRON << "\t";
                    Paraview_File << Conn_Tetr[iNode + 0] - 1 << "\t" << Conn_Tetr[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Tetr[iNode + 2] - 1 << "\t" << Conn_Tetr[iNode + 3] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalHexas; iElem++) {
                    iNode = iElem*N_POINTS_HEXAHEDRON;
                    Paraview_File << N_POINTS_HEXAHEDRON << "\t";
                    Paraview_File << Conn_Hexa[iNode + 0] - 1 << "\t" << Conn_Hexa[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Hexa[iNode + 2] - 1 << "\t" << Conn_Hexa[iNode + 3] - 1 << "\t";
                    Paraview_File << Conn_Hexa[iNode + 4] - 1 << "\t" << Conn_Hexa[iNode + 5] - 1 << "\t";
                    Paraview_File << Conn_Hexa[iNode + 6] - 1 << "\t" << Conn_Hexa[iNode + 7] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalPriss; iElem++) {
                    iNode = iElem*N_POINTS_PRISM;
                    Paraview_File << N_POINTS_PRISM << "\t";
                    Paraview_File << Conn_Pris[iNode + 0] - 1 << "\t" << Conn_Pris[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Pris[iNode + 2] - 1 << "\t" << Conn_Pris[iNode + 3] - 1 << "\t";
                    Paraview_File << Conn_Pris[iNode + 4] - 1 << "\t" << Conn_Pris[iNode + 5] - 1 << "\t";
                }

                for (iElem = 0; iElem < d_numGlobalPyras; iElem++) {
                    iNode = iElem*N_POINTS_PYRAMID;
                    Paraview_File << N_POINTS_PYRAMID << "\t";
                    Paraview_File << Conn_Pyra[iNode + 0] - 1 << "\t" << Conn_Pyra[iNode + 1] - 1 << "\t";
                    Paraview_File << Conn_Pyra[iNode + 2] - 1 << "\t" << Conn_Pyra[iNode + 3] - 1 << "\t";
                    Paraview_File << Conn_Pyra[iNode + 4] - 1 << "\t";
                }
            }

            /*--- Write the header ---*/
            if (surf_sol) Paraview_File << "\nCELL_TYPES " << d_numSurfElems << "\n";
            else Paraview_File << "\nCELL_TYPES " << nGlobal_Elem << "\n";

            if (surf_sol) {
                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) Paraview_File << "3\t";
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) Paraview_File << "5\t";
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) Paraview_File << "9\t";

            }
            else {
                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) Paraview_File << "5\t";
                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) Paraview_File << "9\t";
                for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) Paraview_File << "10\t";
                for (iElem = 0; iElem < d_numGlobalHexas; iElem++) Paraview_File << "12\t";
                for (iElem = 0; iElem < d_numGlobalPriss; iElem++) Paraview_File << "13\t";
                for (iElem = 0; iElem < d_numGlobalPyras; iElem++) Paraview_File << "14\t";
            }



            /*--- Write the header ---*/
            if (config->GetKind_SU2() != SU2_DEF) {
                if (surf_sol) Paraview_File << "\nPOINT_DATA " << d_numSurfPoints << "\n";
                else Paraview_File << "\nPOINT_DATA " << d_numGlobalPoints << "\n";
            }

            unsigned short VarCounter = 0;

            if (config->GetKind_SU2() == SU2_SOL) {

                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class. ---*/
                for (unsigned short iField = 1; iField < config->fields.size(); iField++) {

                    double output_variable = true;
                    size_t found = config->fields[iField].find("\"x\"");
                    if (found != string::npos) output_variable = false;
                    found = config->fields[iField].find("\"y\"");
                    if (found != string::npos) output_variable = false;
                    found = config->fields[iField].find("\"z\"");
                    if (found != string::npos) output_variable = false;

                    if (output_variable) {
                        Paraview_File << "\nSCALARS " << config->fields[iField] << " float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                    }

                    VarCounter++;


                }

            }

            else if (config->GetKind_SU2() != SU2_DEF){

                for (iVar = 0; iVar < nVar_Consv; iVar++) {

                    Paraview_File << "\nSCALARS Conservative_" << iVar + 1 << " float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;
                }

                if (config->GetWrt_Limiters()) {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) {

                        Paraview_File << "\nSCALARS Limiter_" << iVar + 1 << " float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        VarCounter++;
                    }
                }

                if (config->GetWrt_Residuals()) {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) {

                        Paraview_File << "\nSCALARS Residual_" << iVar + 1 << " float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        VarCounter++;
                    }
                }

                /*--- Add names for any extra variables (this will need to be adjusted). ---*/
                if (grid_movement) {

                    Paraview_File << "\nSCALARS Grid_Velx float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Grid_Vely float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    if (nDim == 3) {

                        Paraview_File << "\nSCALARS Grid_Velz float 1\n";
                        Paraview_File << "LOOKUP_TABLE default\n";

                        for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                            if (surf_sol) {
                                if (LocalIndex[iPoint + 1] != 0) {
                                    /*--- Loop over the vars/residuals and write the values to file ---*/
                                    Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                                }
                            }
                            else {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        VarCounter++;

                    }
                }

                if (config->GetKind_Regime() == FREESURFACE) {

                    Paraview_File << "\nSCALARS Density float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {

                    Paraview_File << "\nSCALARS Pressure float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Temperature float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Pressure_Coefficient float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Mach float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {

                    Paraview_File << "\nSCALARS Laminar_Viscosity float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Skin_Friction_Coefficient float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Heat_Flux float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                    Paraview_File << "\nSCALARS Y_Plus float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if (Kind_Solver == RANS) {

                    Paraview_File << "\nSCALARS Eddy_Viscosity float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

                if ((Kind_Solver == ADJ_EULER) ||
                    (Kind_Solver == ADJ_NAVIER_STOKES) ||
                    (Kind_Solver == ADJ_RANS)) {

                    Paraview_File << "\nSCALARS Surface_Sensitivity float 1\n";
                    Paraview_File << "LOOKUP_TABLE default\n";

                    for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {
                        if (surf_sol) {
                            if (LocalIndex[iPoint + 1] != 0) {
                                /*--- Loop over the vars/residuals and write the values to file ---*/
                                Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                            }
                        }
                        else {
                            /*--- Loop over the vars/residuals and write the values to file ---*/
                            Paraview_File << scientific << Data[VarCounter][iPoint] << "\t";
                        }
                    }
                    VarCounter++;

                }

            }

            Paraview_File.close();

            if (surf_sol)  delete[] LocalIndex;
            if (SurfacePoint != NULL) delete[] SurfacePoint;

        }

        void Cae_Output::SetTecplotASCII(Config_p config, GEOM::GEOM_Geometry *geometry, SOLV::SOLV_Solver **solver, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol)
        {
            /*--- Local variables and initialization ---*/
            unsigned short iDim, iVar;
            unsigned short nDim = geometry->GetnDim();
            unsigned short Kind_Solver = config->GetKind_Solver();

            unsigned long iPoint, iElem, iNode;
            unsigned long iExtIter = config->GetExtIter();
            unsigned long *LocalIndex = NULL;
            bool *SurfacePoint = NULL;

            bool grid_movement = config->GetGrid_Movement();
            bool adjoint = config->GetAdjoint();

            char cstr[200], buffer[50];
            string filename;

            /*--- Write file name with extension ---*/
            if (surf_sol)
            {
                if (adjoint) filename = config->GetSurfAdjCoeff_FileName();
                else filename = config->GetSurfFlowCoeff_FileName();
            }
            else
            {
                if (adjoint)
                    filename = config->GetAdj_FileName();
                else
                    filename = config->GetFlow_FileName();
            }

            if (Kind_Solver == LINEAR_ELASTICITY)
            {
                if (surf_sol)
                    filename = config->GetSurfStructure_FileName().c_str();
                else
                    filename = config->GetStructure_FileName().c_str();
            }

            if (Kind_Solver == WAVE_EQUATION)
            {
                if (surf_sol)
                    filename = config->GetSurfWave_FileName().c_str();
                else
                    filename = config->GetWave_FileName().c_str();
            }

            if (Kind_Solver == HEAT_EQUATION)
            {
                if (surf_sol)
                    filename = config->GetSurfHeat_FileName().c_str();
                else
                    filename = config->GetHeat_FileName().c_str();
            }

            if (Kind_Solver == POISSON_EQUATION)
            {
                if (surf_sol)
                    filename = config->GetSurfStructure_FileName().c_str();
                else
                    filename = config->GetStructure_FileName().c_str();
            }

            strcpy(cstr, filename.c_str());

            /*--- Special cases where a number needs to be appended to the file name. ---*/
            if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS ||
                Kind_Solver == ADJ_EULER || Kind_Solver == ADJ_NAVIER_STOKES || Kind_Solver == ADJ_RANS) &&
                (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL))
            {
                sprintf(buffer, "_%d", int(val_iZone));
                strcat(cstr, buffer);
            }

            if (config->GetUnsteady_Simulation() == TIME_SPECTRAL)
            {
                if (config->GetKind_SU2() == SU2_SOL)
                {
                    val_iZone = iExtIter;
                }

                if (int(val_iZone) < 10) 
                    sprintf(buffer, "_0000%d.dat", int(val_iZone));
                if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) 
                    sprintf(buffer, "_000%d.dat", int(val_iZone));
                if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) 
                    sprintf(buffer, "_00%d.dat", int(val_iZone));
                if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) 
                    sprintf(buffer, "_0%d.dat", int(val_iZone));
                if (int(val_iZone) >= 10000) 
                    sprintf(buffer, "_%d.dat", int(val_iZone));
            }
            else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
            {
                if (int(iExtIter) < 10)
                    sprintf(buffer, "_0000%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))
                    sprintf(buffer, "_000%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))
                    sprintf(buffer, "_00%d.dat", int(iExtIter));
                if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000))
                    sprintf(buffer, "_0%d.dat", int(iExtIter));
                if (int(iExtIter) >= 10000)
                    sprintf(buffer, "_%d.dat", int(iExtIter));
            }
            else
            {
                sprintf(buffer, ".dat");
            }

            strcat(cstr, buffer);

            /*--- Open Tecplot ASCII file and write the header. ---*/
            ofstream Tecplot_File;
            Tecplot_File.open(cstr, ios::out);
            Tecplot_File.precision(6);
            if (surf_sol) 
                Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
            else 
                Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;

            /*--- Prepare the variable lists. ---*/

            /*--- Write the list of the fields in the restart file.
            Without including the PointID---*/
            if (config->GetKind_SU2() == SU2_SOL)
            {
                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class. ---*/
                Tecplot_File << "VARIABLES = ";
                nVar_Total = config->fields.size() - 1;
                for (unsigned short iField = 1; iField < config->fields.size(); iField++) 
                {
                    Tecplot_File << config->fields[iField];
                }
                Tecplot_File << endl;
            }
            else
            {
                if (nDim == 2)
                {
                    Tecplot_File << "VARIABLES = \"x\",\"y\"";
                }
                else
                {
                    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
                }

                /*--- Add names for conservative and residual variables ---*/
                for (iVar = 0; iVar < nVar_Consv; iVar++)
                {
                    Tecplot_File << ",\"Conservative_" << iVar + 1 << "\"";
                }

                if (!config->GetLow_MemoryOutput())
                {
                    if (config->GetWrt_Limiters())
                    {
                        for (iVar = 0; iVar < nVar_Consv; iVar++)
                        {
                            Tecplot_File << ",\"Limiter_" << iVar + 1 << "\"";
                        }
                    }

                    if (config->GetWrt_Residuals())
                    {
                        for (iVar = 0; iVar < nVar_Consv; iVar++) 
						{
                            Tecplot_File << ",\"Residual_" << iVar + 1 << "\"";
                        }
                    }

                    /*--- Add names for any extra variables (this will need to be adjusted). ---*/
                    if (grid_movement)
                    {
                        if (nDim == 2)
                        {
                            Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\"";
                        }
                        else
                        {
                            Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\",\"Grid_Velz\"";
                        }
                    }

                    if (config->GetKind_Regime() == FREESURFACE)
                    {
                        Tecplot_File << ",\"Density\"";
                    }

                    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                    {
                        Tecplot_File << ",\"Pressure\",\"Temperature\",\"Pressure_Coefficient\",\"Mach\"";
                    }

                    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                    {
                        Tecplot_File << ",\"Laminar_Viscosity\", \"Skin_Friction_Coefficient\", \"Heat_Flux\", \"Y_Plus\"";
                    }

                    if (Kind_Solver == RANS)
                    {
                        Tecplot_File << ", \"Eddy_Viscosity\"";
                    }

                    if (config->GetWrt_SharpEdges())
                    {
                        if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                        {
                            Tecplot_File << ", \"Sharp_Edge_Dist\"";
                        }
                    }

                    if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES))
                    {
                        Tecplot_File << ",\"Mach\",\"Pressure\",\"Temperature\",\"Temperature_ve\"";
                    }

                    if (Kind_Solver == TNE2_NAVIER_STOKES)
                    {
                        for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
                            Tecplot_File << ",\"DiffusionCoeff_" << iSpecies << "\"";
                        Tecplot_File << ",\"Laminar_Viscosity\",\"ThermConductivity\",\"ThermConductivity_ve\"";
                    }

                    if (Kind_Solver == POISSON_EQUATION)
                    {
                        unsigned short iDim;
                        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                            Tecplot_File << ",\"poissonField_" << iDim + 1 << "\"";
                    }

                    if ((Kind_Solver == ADJ_EULER) ||
                        (Kind_Solver == ADJ_NAVIER_STOKES) ||
                        (Kind_Solver == ADJ_RANS) ||
                        (Kind_Solver == ADJ_TNE2_EULER) ||
                        (Kind_Solver == ADJ_TNE2_NAVIER_STOKES))
                    {
                        Tecplot_File << ", \"Surface_Sensitivity\", \"Solution_Sensor\"";
                    }

                    if (Kind_Solver == LINEAR_ELASTICITY)
                    {
                        Tecplot_File << ", \"Von_Mises_Stress\", \"Flow_Pressure\"";
                    }

                    if (config->GetExtraOutput())
                    {
                        string *headings = NULL;
                        //if (Kind_Solver == RANS) {
                        headings = solver[TURB_SOL]->OutputHeadingNames;
                        //}
                        for (iVar = 0; iVar < nVar_Extra; iVar++)
                        {
                            //Tecplot_File << ", \"ExtraOutput_" << iVar+1<<"\"";
                            if (headings == NULL)
                            {
                                Tecplot_File << ", \"ExtraOutput_" << iVar + 1 << "\"";
                            }
                            else
                            {
                                Tecplot_File << ", \"" << headings[iVar] << "\"";
                            }
                        }
                    }
                }

                Tecplot_File << endl;
            }

            /*--- If it's a surface output, print only the points
            that are in the element list, change the numbering ---*/

            if (surf_sol)
            {
                LocalIndex = new unsigned long[d_numGlobalPoints + 1];
                SurfacePoint = new bool[d_numGlobalPoints + 1];

                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) 
					SurfacePoint[iPoint] = false;

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++)
                {
                    iNode = iElem*N_POINTS_LINE;
                    SurfacePoint[Conn_Line[iNode + 0]] = true;
                    SurfacePoint[Conn_Line[iNode + 1]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++)
                {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    SurfacePoint[Conn_BoundTria[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 2]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++)
                {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    SurfacePoint[Conn_BoundQuad[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 2]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 3]] = true;
                }

                d_numSurfPoints = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++)
                {
                    LocalIndex[iPoint] = 0;
                    if (SurfacePoint[iPoint])
                    {
                        d_numSurfPoints++;
                        LocalIndex[iPoint] = d_numSurfPoints;
                    }
                }
            }

            /*--- Write the header ---*/
            Tecplot_File << "ZONE ";
            if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady())
            {
                Tecplot_File << "STRANDID=" << int(iExtIter + 1) << ", SOLUTIONTIME=" << config->GetDelta_UnstTime()*iExtIter << ", ";
            }
            else if (config->GetUnsteady_Simulation() == TIME_SPECTRAL)
            {
                /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
                double period = config->GetTimeSpectral_Period();
                double deltaT = period / (double)(config->GetnTimeInstances());
                Tecplot_File << "STRANDID=" << int(iExtIter + 1) << ", SOLUTIONTIME=" << deltaT*iExtIter << ", ";
            }

            if (nDim == 2)
            {
                if (surf_sol) Tecplot_File << "NODES= " << d_numSurfPoints << ", ELEMENTS= " << d_numSurfElems << ", DATAPACKING=POINT, ZONETYPE=FELINESEG" << endl;
                else Tecplot_File << "NODES= " << d_numGlobalPoints << ", ELEMENTS= " << nGlobal_Elem << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
            }
            else
            {
                if (surf_sol) Tecplot_File << "NODES= " << d_numSurfPoints << ", ELEMENTS= " << d_numSurfElems << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
                else Tecplot_File << "NODES= " << d_numGlobalPoints << ", ELEMENTS= " << nGlobal_Elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << endl;
            }

            /*--- Write surface and volumetric solution data. ---*/
            for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++)
            {
                if (surf_sol)
                {
                    if (LocalIndex[iPoint + 1] != 0)
                    {
                        /*--- Write the node coordinates ---*/
                        if (config->GetKind_SU2() != SU2_SOL)
                        {
                            for (iDim = 0; iDim < nDim; iDim++)
                                Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
                        }

                        /*--- Loop over the vars/residuals and write the values to file ---*/
                        for (iVar = 0; iVar < nVar_Total; iVar++)
                            Tecplot_File << scientific << Data[iVar][iPoint] << "\t";

                        Tecplot_File << endl;
                    }
                }
                else
                {
                    /*--- Write the node coordinates ---*/
                    if (config->GetKind_SU2() != SU2_SOL)
                    {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
                    }

                    /*--- Loop over the vars/residuals and write the values to file ---*/
                    for (iVar = 0; iVar < nVar_Total; iVar++)
                        Tecplot_File << scientific << Data[iVar][iPoint] << "\t";

                    Tecplot_File << endl;
                }
            }


            /*--- Write connectivity data. ---*/
            if (surf_sol) 
            {
                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) 
                {
                    iNode = iElem*N_POINTS_LINE;
                    Tecplot_File << LocalIndex[Conn_Line[iNode + 0]] << "\t";
                    Tecplot_File << LocalIndex[Conn_Line[iNode + 1]] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) 
                {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 0]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 1]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 2]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 2]] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++)
                {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 0]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 1]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 2]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 3]] << "\n";
                }

            }
            else 
            {
                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) 
                {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Tecplot_File << Conn_Tria[iNode + 0] << "\t";
                    Tecplot_File << Conn_Tria[iNode + 1] << "\t";
                    Tecplot_File << Conn_Tria[iNode + 2] << "\t";
                    Tecplot_File << Conn_Tria[iNode + 2] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) 
                {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Tecplot_File << Conn_Quad[iNode + 0] << "\t";
                    Tecplot_File << Conn_Quad[iNode + 1] << "\t";
                    Tecplot_File << Conn_Quad[iNode + 2] << "\t";
                    Tecplot_File << Conn_Quad[iNode + 3] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) 
                {
                    iNode = iElem*N_POINTS_TETRAHEDRON;
                    Tecplot_File << Conn_Tetr[iNode + 0] << "\t" << Conn_Tetr[iNode + 1] << "\t";
                    Tecplot_File << Conn_Tetr[iNode + 2] << "\t" << Conn_Tetr[iNode + 2] << "\t";
                    Tecplot_File << Conn_Tetr[iNode + 3] << "\t" << Conn_Tetr[iNode + 3] << "\t";
                    Tecplot_File << Conn_Tetr[iNode + 3] << "\t" << Conn_Tetr[iNode + 3] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalHexas; iElem++) 
                {
                    iNode = iElem*N_POINTS_HEXAHEDRON;
                    Tecplot_File << Conn_Hexa[iNode + 0] << "\t" << Conn_Hexa[iNode + 1] << "\t";
                    Tecplot_File << Conn_Hexa[iNode + 2] << "\t" << Conn_Hexa[iNode + 3] << "\t";
                    Tecplot_File << Conn_Hexa[iNode + 4] << "\t" << Conn_Hexa[iNode + 5] << "\t";
                    Tecplot_File << Conn_Hexa[iNode + 6] << "\t" << Conn_Hexa[iNode + 7] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalPriss; iElem++) 
                {
                    iNode = iElem*N_POINTS_PRISM;
                    Tecplot_File << Conn_Pris[iNode + 0] << "\t" << Conn_Pris[iNode + 1] << "\t";
                    Tecplot_File << Conn_Pris[iNode + 1] << "\t" << Conn_Pris[iNode + 2] << "\t";
                    Tecplot_File << Conn_Pris[iNode + 3] << "\t" << Conn_Pris[iNode + 4] << "\t";
                    Tecplot_File << Conn_Pris[iNode + 4] << "\t" << Conn_Pris[iNode + 5] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalPyras; iElem++) 
                {
                    iNode = iElem*N_POINTS_PYRAMID;
                    Tecplot_File << Conn_Pyra[iNode + 0] << "\t" << Conn_Pyra[iNode + 1] << "\t";
                    Tecplot_File << Conn_Pyra[iNode + 2] << "\t" << Conn_Pyra[iNode + 3] << "\t";
                    Tecplot_File << Conn_Pyra[iNode + 4] << "\t" << Conn_Pyra[iNode + 4] << "\t";
                    Tecplot_File << Conn_Pyra[iNode + 4] << "\t" << Conn_Pyra[iNode + 4] << "\n";
                }
            }

            Tecplot_File.close();

            if (surf_sol) delete[] LocalIndex;
        }
        void Cae_Output::SetTecplotASCII_LowMemory(Config_p config, GEOM::GEOM_Geometry *geometry, SOLV::SOLV_Solver **solver, char mesh_filename[MAX_STRING_SIZE], bool surf_sol)
        {
            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            int size = SINGLE_NODE;
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            unsigned long iElem, iPoint;
            unsigned short iVar;
            bool grid_movement = config->GetGrid_Movement();
            unsigned short Kind_Solver = config->GetKind_Solver();
            ofstream Tecplot_File;
            unsigned long Total_nElem_Bound, *PointSurface = NULL, nPointSurface = 0;
            unsigned short iMarker;

            /*--- Open Tecplot ASCII file and write the header. ---*/

            Tecplot_File.open(mesh_filename, ios::out);
            Tecplot_File.precision(6);
            if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
            else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;

            /*--- Write the list of the fields in the restart file.
            Without including the PointID---*/
            if (config->GetKind_SU2() == SU2_SOL) {

                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class. ---*/
                Tecplot_File << "VARIABLES = ";
                nVar_Total = config->fields.size() - 1;
                for (unsigned short iField = 1; iField < config->fields.size(); iField++) {
                    Tecplot_File << config->fields[iField];
                }

                Tecplot_File << endl;

            }
            else {

                if (geometry->GetnDim() == 2) { Tecplot_File << "VARIABLES = \"x\",\"y\""; }
                else { Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\""; }

                /*--- Add names for conservative, limiters, and residual variables ---*/

                for (iVar = 0; iVar < nVar_Consv; iVar++) {
                    Tecplot_File << ",\"Conservative_" << iVar + 1 << "\"";
                }

                if (!config->GetLow_MemoryOutput()) {

                    if (config->GetWrt_Limiters()) {
                        for (iVar = 0; iVar < nVar_Consv; iVar++) {
                            Tecplot_File << ",\"Limiter_" << iVar + 1 << "\"";
                        }
                    }
                    if (config->GetWrt_Residuals()) {
                        for (iVar = 0; iVar < nVar_Consv; iVar++) {
                            Tecplot_File << ",\"Residual_" << iVar + 1 << "\"";
                        }
                    }

                    /*--- Add names for any extra variables (this will need to be adjusted). ---*/

                    if (grid_movement) {
                        if (geometry->GetnDim() == 2) {
                            Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\"";
                        }
                        else {
                            Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\",\"Grid_Velz\"";
                        }
                    }

                    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
                        Tecplot_File << ",\"Pressure\",\"Temperature\",\"Pressure_Coefficient\",\"Mach\"";
                        if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
                            Tecplot_File << ",\"Laminar_Viscosity\", \"Skin_Friction_Coefficient\", \"Heat_Flux\", \"Y_Plus\"";
                            if (Kind_Solver == RANS) { Tecplot_File << ", \"Eddy_Viscosity\""; }
                        }
                        if (config->GetWrt_SharpEdges()) { Tecplot_File << ", \"Sharp_Edge_Dist\""; }
                    }

                    if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) {
                        Tecplot_File << ", \"Surface_Sensitivity\", \"Solution_Sensor\"";
                    }

                }

                Tecplot_File << endl;

            }


            if (surf_sol) {

                /*--- It is important to do a renumbering to don't add points
                that do not belong to the surfaces ---*/

                PointSurface = new unsigned long[geometry->GetnPoint()];
                for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    if (geometry->node[iPoint]->GetBoundary()) {
                        PointSurface[iPoint] = nPointSurface;
                        nPointSurface++;
                    }

                /*--- Compute the total number of elements ---*/

                Total_nElem_Bound = 0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                    if (config->GetMarker_All_Plotting(iMarker) == YES) {
                        Total_nElem_Bound += geometry->GetnElem_Bound(iMarker);
                    }
                }

                if (Total_nElem_Bound != 0) {

                    /*--- Write the header of the file ---*/

                    Tecplot_File << "ZONE T= \"MPI rank: " << rank << "\", ";
                    Tecplot_File << "NODES= " << nPointSurface << ", ELEMENTS= " << Total_nElem_Bound << ", DATAPACKING= POINT";
                    if (geometry->GetnDim() == 2) Tecplot_File << ", ZONETYPE= FELINESEG" << endl;
                    if (geometry->GetnDim() == 3) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL" << endl;

                    /*--- Only write the coordiantes of the points that are on the surfaces ---*/

                    if (geometry->GetnDim() == 3) {
                        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                            if (geometry->node[iPoint]->GetBoundary()) {
                                for (iVar = 0; iVar < solver[FLOW_SOL]->GetnVar(); iVar++)
                                    Tecplot_File << scientific << solver[FLOW_SOL]->node[iPoint]->GetSolution(iVar) << "\t";
                                Tecplot_File << "\n";
                            }
                    }
                    else {
                        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                            if (geometry->node[iPoint]->GetBoundary()) {
                                for (iVar = 0; iVar < solver[FLOW_SOL]->GetnVar(); iVar++)
                                    Tecplot_File << scientific << solver[FLOW_SOL]->node[iPoint]->GetSolution(iVar) << "\t";
                                Tecplot_File << "\n";
                            }
                    }

                    /*--- Write the cells using the new numbering ---*/

                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                        if (config->GetMarker_All_Plotting(iMarker) == YES)
                            for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
                                if (geometry->GetnDim() == 2) {
                                    Tecplot_File << PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)] + 1 << " "
                                        << PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)] + 1 << endl;
                                }
                                if (geometry->GetnDim() == 3) {
                                    if (geometry->bound[iMarker][iElem]->GetnNodes() == 3) {
                                        Tecplot_File << PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)] + 1 << " "
                                            << PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)] + 1 << " "
                                            << PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)] + 1 << " "
                                            << PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)] + 1 << endl;
                                    }
                                    if (geometry->bound[iMarker][iElem]->GetnNodes() == 4) {
                                        Tecplot_File << PointSurface[geometry->bound[iMarker][iElem]->GetNode(0)] + 1 << " "
                                            << PointSurface[geometry->bound[iMarker][iElem]->GetNode(1)] + 1 << " "
                                            << PointSurface[geometry->bound[iMarker][iElem]->GetNode(2)] + 1 << " "
                                            << PointSurface[geometry->bound[iMarker][iElem]->GetNode(3)] + 1 << endl;
                                    }
                                }
                            }
                }
                else {

                    /*--- No elements in the surface ---*/

                    if (geometry->GetnDim() == 2) {
                        Tecplot_File << "ZONE ";
                        Tecplot_File << "T= \"MPI rank: " << rank << "\", ";
                        Tecplot_File << "NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG" << endl;
                        for (iVar = 0; iVar < solver[FLOW_SOL]->GetnVar(); iVar++)
                            Tecplot_File << scientific << "0.0\t";
                        Tecplot_File << "\n";
                        Tecplot_File << "1 1" << endl;
                    }
                    if (geometry->GetnDim() == 3) {
                        Tecplot_File << "ZONE ";
                        Tecplot_File << "T= \"MPI rank: " << rank << "\", ";
                        Tecplot_File << "NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
                        for (iVar = 0; iVar < solver[FLOW_SOL]->GetnVar(); iVar++)
                            Tecplot_File << scientific << "0.0\t";
                        Tecplot_File << "\n";
                        Tecplot_File << "1 1 1 1" << endl;
                    }
                }

                /*--- Dealocate memory and close the file ---*/

                delete[] PointSurface;
                Tecplot_File.close();

            }

            else {

                Tecplot_File << "ZONE ";
                Tecplot_File << "T= \"MPI rank: " << rank << "\", ";
                Tecplot_File << "NODES= " << geometry->GetnPoint() << ", ELEMENTS= " << geometry->GetnElem() << ", DATAPACKING= POINT";
                if (geometry->GetnDim() == 2) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL" << endl;
                if (geometry->GetnDim() == 3) Tecplot_File << ", ZONETYPE= FEBRICK" << endl;

                /*--- Adding coordinates ---*/

                for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
                    for (iVar = 0; iVar < solver[FLOW_SOL]->GetnVar(); iVar++)
                        Tecplot_File << scientific << solver[FLOW_SOL]->node[iPoint]->GetSolution(iVar) << "\t";
                    Tecplot_File << "\n";
                }

                /*--- Adding conectivity ---*/

                for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
                    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
                        Tecplot_File <<
                            geometry->elem[iElem]->GetNode(0) + 1 << " " << geometry->elem[iElem]->GetNode(1) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(2) + 1 << " " << geometry->elem[iElem]->GetNode(2) + 1 << endl;
                    }
                    if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {
                        Tecplot_File <<
                            geometry->elem[iElem]->GetNode(0) + 1 << " " << geometry->elem[iElem]->GetNode(1) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(2) + 1 << " " << geometry->elem[iElem]->GetNode(3) + 1 << endl;
                    }
                    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
                        Tecplot_File <<
                            geometry->elem[iElem]->GetNode(0) + 1 << " " << geometry->elem[iElem]->GetNode(1) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(2) + 1 << " " << geometry->elem[iElem]->GetNode(2) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(3) + 1 << " " << geometry->elem[iElem]->GetNode(3) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(3) + 1 << " " << geometry->elem[iElem]->GetNode(3) + 1 << endl;
                    }
                    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
                        Tecplot_File <<
                            geometry->elem[iElem]->GetNode(0) + 1 << " " << geometry->elem[iElem]->GetNode(1) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(2) + 1 << " " << geometry->elem[iElem]->GetNode(3) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(4) + 1 << " " << geometry->elem[iElem]->GetNode(5) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(6) + 1 << " " << geometry->elem[iElem]->GetNode(7) + 1 << endl;
                    }
                    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID) {
                        Tecplot_File <<
                            geometry->elem[iElem]->GetNode(0) + 1 << " " << geometry->elem[iElem]->GetNode(1) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(2) + 1 << " " << geometry->elem[iElem]->GetNode(3) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(4) + 1 << " " << geometry->elem[iElem]->GetNode(4) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(4) + 1 << " " << geometry->elem[iElem]->GetNode(4) + 1 << endl;
                    }
                    if (geometry->elem[iElem]->GetVTK_Type() == PRISM) {
                        Tecplot_File <<
                            geometry->elem[iElem]->GetNode(0) + 1 << " " << geometry->elem[iElem]->GetNode(1) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(1) + 1 << " " << geometry->elem[iElem]->GetNode(2) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(3) + 1 << " " << geometry->elem[iElem]->GetNode(4) + 1 << " " <<
                            geometry->elem[iElem]->GetNode(4) + 1 << " " << geometry->elem[iElem]->GetNode(5) + 1 << endl;
                    }
                }

                Tecplot_File.close();

            }


#ifdef HAVE_MPI

            /*--- Add solution files to a single file ---*/

            if (rank == MASTER_NODE) {

                ofstream Tecplot_File;
                string filename, text_line;
                char buffer_char[50], out_file[MAX_STRING_SIZE];

                if (!config->GetAdjoint()) {
                    if (surf_sol) filename = config->GetSurfFlowCoeff_FileName();
                    else filename = config->GetFlow_FileName();
                }
                else {
                    if (surf_sol) filename = config->GetSurfAdjCoeff_FileName();
                    else filename = config->GetAdj_FileName();
                }

                strcpy(mesh_filename, filename.c_str());
                sprintf(buffer_char, ".dat");
                strcat(mesh_filename, buffer_char);

                Tecplot_File.open(mesh_filename, ios::out);

                for (int iRank = 0; iRank < size; iRank++) {

                    if (!config->GetAdjoint()) {
                        if (surf_sol) filename = config->GetSurfFlowCoeff_FileName();
                        else filename = config->GetFlow_FileName();
                    }
                    else {
                        if (surf_sol) filename = config->GetSurfAdjCoeff_FileName();
                        else filename = config->GetAdj_FileName();
                    }

                    strcpy(out_file, filename.c_str());
                    sprintf(buffer_char, "_%i.dat", iRank + 1);
                    strcat(out_file, buffer_char);
                    ifstream Tecplot_File_;
                    Tecplot_File_.open(out_file, ios::in);
                    while (getline(Tecplot_File_, text_line)) {
                        Tecplot_File << text_line << endl;
                    }
                    Tecplot_File_.close();
                    remove(out_file);
                }

                Tecplot_File.close();

            }

#endif

        }
        void Cae_Output::SetTecplotASCII_Mesh(Config_p config, GEOM::GEOM_Geometry *geometry, bool surf_sol, bool new_file) 
		{
            unsigned short iDim, nDim = geometry->GetnDim();
            unsigned long iPoint, iElem, iNode;
            unsigned long *LocalIndex = NULL;
            bool *SurfacePoint = NULL;
            char cstr[200];
            ofstream Tecplot_File;

            if (surf_sol) 
				strcpy(cstr, "surface_grid.dat");
            else 
				strcpy(cstr, "volumetric_grid.dat");

            /*--- Open Tecplot ASCII file and write the header. ---*/

            if (new_file) {
                Tecplot_File.open(cstr, ios::out);
                Tecplot_File.precision(6);
                if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
                else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;

                if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\"";
                else Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
            }
            else Tecplot_File.open(cstr, ios::out | ios::app);
            Tecplot_File << endl;

            /*--- If it's a surface output, print only the points
            that are in the element list, change the numbering ---*/

            if (surf_sol) {

                LocalIndex = new unsigned long[d_numGlobalPoints + 1];
                SurfacePoint = new bool[d_numGlobalPoints + 1];

                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) SurfacePoint[iPoint] = false;

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    SurfacePoint[Conn_Line[iNode + 0]] = true;
                    SurfacePoint[Conn_Line[iNode + 1]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    SurfacePoint[Conn_BoundTria[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 2]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    SurfacePoint[Conn_BoundQuad[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 2]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 3]] = true;
                }

                d_numSurfPoints = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                    LocalIndex[iPoint] = 0;
                    if (SurfacePoint[iPoint]) { d_numSurfPoints++; LocalIndex[iPoint] = d_numSurfPoints; }
                }

            }

            /*--- Write the header ---*/

            Tecplot_File << "ZONE T= ";
            if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
            else Tecplot_File << "\"Deformed grid\", C=RED, ";

            if (nDim == 2) {
                if (surf_sol) Tecplot_File << "NODES= " << d_numSurfPoints << ", ELEMENTS= " << d_numSurfElems << ", DATAPACKING=POINT, ZONETYPE=FELINESEG" << endl;
                else Tecplot_File << "NODES= " << d_numGlobalPoints << ", ELEMENTS= " << d_numGlobalElems << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
            }
            else {
                if (surf_sol) Tecplot_File << "NODES= " << d_numSurfPoints << ", ELEMENTS= " << d_numSurfElems << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
                else Tecplot_File << "NODES= " << d_numGlobalPoints << ", ELEMENTS= " << d_numGlobalElems << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << endl;
            }

            /*--- Write surface and volumetric solution data. ---*/

            for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) {

                if (surf_sol) {

                    if (LocalIndex[iPoint + 1] != 0) {

                        /*--- Write the node coordinates ---*/
                        for (iDim = 0; iDim < nDim; iDim++)
                            Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";

                        Tecplot_File << endl;

                    }

                }
                else {

                    /*--- Write the node coordinates ---*/

                    for (iDim = 0; iDim < nDim; iDim++)
                        Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";


                    Tecplot_File << endl;

                }

            }


            /*--- Write connectivity data. ---*/

            if (surf_sol) {

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    Tecplot_File << LocalIndex[Conn_Line[iNode + 0]] << "\t";
                    Tecplot_File << LocalIndex[Conn_Line[iNode + 1]] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 0]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 1]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 2]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundTria[iNode + 2]] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 0]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 1]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 2]] << "\t";
                    Tecplot_File << LocalIndex[Conn_BoundQuad[iNode + 3]] << "\n";
                }

            }
            else {

                for (iElem = 0; iElem < d_numGlobalTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    Tecplot_File << Conn_Tria[iNode + 0] << "\t";
                    Tecplot_File << Conn_Tria[iNode + 1] << "\t";
                    Tecplot_File << Conn_Tria[iNode + 2] << "\t";
                    Tecplot_File << Conn_Tria[iNode + 2] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    Tecplot_File << Conn_Quad[iNode + 0] << "\t";
                    Tecplot_File << Conn_Quad[iNode + 1] << "\t";
                    Tecplot_File << Conn_Quad[iNode + 2] << "\t";
                    Tecplot_File << Conn_Quad[iNode + 3] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalTetrs; iElem++) {
                    iNode = iElem*N_POINTS_TETRAHEDRON;
                    Tecplot_File << Conn_Tetr[iNode + 0] << "\t" << Conn_Tetr[iNode + 1] << "\t";
                    Tecplot_File << Conn_Tetr[iNode + 2] << "\t" << Conn_Tetr[iNode + 2] << "\t";
                    Tecplot_File << Conn_Tetr[iNode + 3] << "\t" << Conn_Tetr[iNode + 3] << "\t";
                    Tecplot_File << Conn_Tetr[iNode + 3] << "\t" << Conn_Tetr[iNode + 3] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalHexas; iElem++) {
                    iNode = iElem*N_POINTS_HEXAHEDRON;
                    Tecplot_File << Conn_Hexa[iNode + 0] << "\t" << Conn_Hexa[iNode + 1] << "\t";
                    Tecplot_File << Conn_Hexa[iNode + 2] << "\t" << Conn_Hexa[iNode + 3] << "\t";
                    Tecplot_File << Conn_Hexa[iNode + 4] << "\t" << Conn_Hexa[iNode + 5] << "\t";
                    Tecplot_File << Conn_Hexa[iNode + 6] << "\t" << Conn_Hexa[iNode + 7] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalPriss; iElem++) {
                    iNode = iElem*N_POINTS_PRISM;
                    Tecplot_File << Conn_Pris[iNode + 0] << "\t" << Conn_Pris[iNode + 1] << "\t";
                    Tecplot_File << Conn_Pris[iNode + 1] << "\t" << Conn_Pris[iNode + 2] << "\t";
                    Tecplot_File << Conn_Pris[iNode + 3] << "\t" << Conn_Pris[iNode + 4] << "\t";
                    Tecplot_File << Conn_Pris[iNode + 4] << "\t" << Conn_Pris[iNode + 5] << "\n";
                }

                for (iElem = 0; iElem < d_numGlobalPyras; iElem++) {
                    iNode = iElem*N_POINTS_PYRAMID;
                    Tecplot_File << Conn_Pyra[iNode + 0] << "\t" << Conn_Pyra[iNode + 1] << "\t";
                    Tecplot_File << Conn_Pyra[iNode + 2] << "\t" << Conn_Pyra[iNode + 3] << "\t";
                    Tecplot_File << Conn_Pyra[iNode + 4] << "\t" << Conn_Pyra[iNode + 4] << "\t";
                    Tecplot_File << Conn_Pyra[iNode + 4] << "\t" << Conn_Pyra[iNode + 4] << "\n";
                }
            }

            Tecplot_File.close();

            if (surf_sol) delete[] LocalIndex;

        }

        void Cae_Output::SetTecplotBinary_DomainMesh(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone) {

#ifdef HAVE_TECIO

            double   t;
            INTEGER4 i, N, err, Debug, NPts, NElm, IsDouble, KMax;
            INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
            INTEGER4 *ShareFromZone = NULL, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
            string buffer, variables;
            stringstream file;
            bool first_zone = true;
            unsigned short dims = geometry->GetnDim();
            enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
            enum	 ZoneType { ORDERED = 0, FELINESEG = 1, FETRIANGLE = 2, FEQUADRILATERAL = 3, FETETRAHEDRON = 4, FEBRICK = 5, FEPOLYGON = 6, FEPOLYHEDRON = 7 };

            /*--- Consistent data for Tecplot zones ---*/

            Debug = 0;
            IsDouble = 1;
            NPts = (INTEGER4)d_numGlobalPoints;
            t = 0.0;//iExtIter*config->GetDelta_UnstTimeND();
            KMax = 0;
            ICellMax = 0;
            JCellMax = 0;
            KCellMax = 0;
            StrandID = 0;//(INTEGER4)iExtIter;
            ParentZn = 0;
            IsBlock = 1;
            NumFaceConnections = 0;
            FaceNeighborMode = 0;
            ShareConnectivityFromZone = 0;

            /*--- Write Tecplot solution file ---*/

            if (!d_isBaseOutput) {

                file.str(string());
                buffer = config->GetFlow_FileName();

                file << buffer << ".mesh.plt";
                FileType = GRID;

                if (dims == 2) variables = "x y";
                else if (dims == 3) variables = "x y z";
                else cout << "Error: wrong number of dimensions: " << dims << endl;

                /*--- Open Tecplot file ---*/

                err = TECINI112((char *)config->GetFlow_FileName().c_str(),
                    (char *)variables.c_str(),
                    (char *)file.str().c_str(),
                    (char *)".",
                    &FileType,
                    &Debug,
                    &IsDouble);
                if (err) cout << "Error in opening Tecplot file" << endl;

                first_zone = true;
                //    ShareFromZone = new INTEGER4[dims];
                //    for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                if (d_numGlobalTrias > 0) {

                    /*--- Write the zone header information ---*/
                    ZoneType = FETRIANGLE; NElm = (INTEGER4)d_numGlobalTrias; N = NElm*N_POINTS_TRIANGLE;

                    err = TECZNE112((char*)"Triangle Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/

                    if (first_zone) {

                        ShareFromZone = new INTEGER4[dims];
                        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                        if (config->GetKind_SU2() == SU2_SOL) {
                            err = TECDAT112(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Data[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        else {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_Tria);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                }
                if (d_numGlobalQuads > 0) {

                    /*--- Write the zone header information ---*/

                    ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)d_numGlobalQuads; N = NElm*N_POINTS_QUADRILATERAL;

                    err = TECZNE112((char*)"Quadrilateral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/

                    if (first_zone) {

                        ShareFromZone = new INTEGER4[dims];
                        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                        if (config->GetKind_SU2() == SU2_SOL) {
                            err = TECDAT112(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Data[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        else {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_Quad);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                }
                if (d_numGlobalTetrs > 0) {

                    /*--- Write the zone header information ---*/

                    ZoneType = FETETRAHEDRON; NElm = (INTEGER4)d_numGlobalTetrs; N = NElm*N_POINTS_TETRAHEDRON;

                    err = TECZNE112((char*)"Tetrahedral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/

                    if (first_zone) {

                        ShareFromZone = new INTEGER4[dims];
                        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                        if (config->GetKind_SU2() == SU2_SOL) {
                            err = TECDAT112(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Data[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        else {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_Tetr);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                }
                if (d_numGlobalHexas > 0) {

                    /*--- Write the zone header information ---*/

                    ZoneType = FEBRICK; NElm = (INTEGER4)d_numGlobalHexas; N = NElm*N_POINTS_HEXAHEDRON;

                    err = TECZNE112((char*)"Hexahedral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/

                    if (first_zone) {

                        ShareFromZone = new INTEGER4[dims];
                        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                        if (config->GetKind_SU2() == SU2_SOL) {
                            err = TECDAT112(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Data[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        else {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_Hexa);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                }

                if (d_numGlobalPyras > 0) {

                    /*--- Here, we reuse the hex implementation to write pyramid elements.
                    Write the zone header information. ---*/
                    ZoneType = FEBRICK; NElm = (INTEGER4)d_numGlobalPyras; N = NElm*N_POINTS_HEXAHEDRON;

                    err = TECZNE112((char*)"Pyramid Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/

                    if (first_zone) {

                        ShareFromZone = new INTEGER4[dims];
                        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                        if (config->GetKind_SU2() == SU2_SOL) {
                            err = TECDAT112(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Data[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        else {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        if (err) cout << "Error writing grid coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    /*--- Convert the pyramid connectivity from 5 nodes to 8 nodes for FEBRICK ---*/
                    int *Conn_Pyra_Mod = new int[d_numGlobalPyras*N_POINTS_HEXAHEDRON];
                    unsigned long iNode_Pyra, iNode_Hexa;
                    for (unsigned long iElem = 0; iElem < d_numGlobalPyras; iElem++) {
                        iNode_Pyra = iElem*N_POINTS_PYRAMID;
                        iNode_Hexa = iElem*N_POINTS_HEXAHEDRON;
                        Conn_Pyra_Mod[iNode_Hexa + 0] = Conn_Pyra[iNode_Pyra + 4];
                        Conn_Pyra_Mod[iNode_Hexa + 1] = Conn_Pyra[iNode_Pyra + 4];
                        Conn_Pyra_Mod[iNode_Hexa + 2] = Conn_Pyra[iNode_Pyra + 4];
                        Conn_Pyra_Mod[iNode_Hexa + 3] = Conn_Pyra[iNode_Pyra + 4];
                        Conn_Pyra_Mod[iNode_Hexa + 4] = Conn_Pyra[iNode_Pyra + 0];
                        Conn_Pyra_Mod[iNode_Hexa + 5] = Conn_Pyra[iNode_Pyra + 1];
                        Conn_Pyra_Mod[iNode_Hexa + 6] = Conn_Pyra[iNode_Pyra + 2];
                        Conn_Pyra_Mod[iNode_Hexa + 7] = Conn_Pyra[iNode_Pyra + 3];
                    }
                    err = TECNOD112(Conn_Pyra_Mod);
                    if (err) cout << "Error writing pyramid connectivity to Tecplot file" << endl;
                    delete[] Conn_Pyra_Mod;

                }

                if (d_numGlobalPriss > 0) {

                    /*--- Here, we reuse the hex implementation to write prism elements.
                    Write the zone header information ---*/
                    ZoneType = FEBRICK; NElm = (INTEGER4)d_numGlobalPriss; N = NElm*N_POINTS_HEXAHEDRON;

                    err = TECZNE112((char*)"Prism Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/

                    if (first_zone) {

                        ShareFromZone = new INTEGER4[dims];
                        for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                        if (config->GetKind_SU2() == SU2_SOL) {
                            err = TECDAT112(&NPts, Data[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Data[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Data[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        else {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
                            if (geometry->GetnDim() == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                ShareFromZone[2] = 1;
                            }
                        }
                        if (err) cout << "Error writing grid coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    /*--- Convert the prism connectivity from 6 nodes to 8 nodes for FEBRICK ---*/
                    int *Conn_Pris_Mod = new int[d_numGlobalPriss*N_POINTS_HEXAHEDRON];
                    unsigned long iNode_Pris, iNode_Hexa;
                    for (unsigned long iElem = 0; iElem < d_numGlobalPriss; iElem++) {
                        iNode_Pris = iElem*N_POINTS_PRISM;
                        iNode_Hexa = iElem*N_POINTS_HEXAHEDRON;
                        Conn_Pris_Mod[iNode_Hexa + 0] = Conn_Pris[iNode_Pris + 0];
                        Conn_Pris_Mod[iNode_Hexa + 1] = Conn_Pris[iNode_Pris + 0];
                        Conn_Pris_Mod[iNode_Hexa + 2] = Conn_Pris[iNode_Pris + 1];
                        Conn_Pris_Mod[iNode_Hexa + 3] = Conn_Pris[iNode_Pris + 2];
                        Conn_Pris_Mod[iNode_Hexa + 4] = Conn_Pris[iNode_Pris + 3];
                        Conn_Pris_Mod[iNode_Hexa + 5] = Conn_Pris[iNode_Pris + 3];
                        Conn_Pris_Mod[iNode_Hexa + 6] = Conn_Pris[iNode_Pris + 4];
                        Conn_Pris_Mod[iNode_Hexa + 7] = Conn_Pris[iNode_Pris + 5];
                    }
                    err = TECNOD112(Conn_Pris_Mod);
                    if (err) cout << "Error writing prism connectivity to Tecplot file" << endl;
                    delete[] Conn_Pris_Mod;

                }

                delete[] ShareFromZone;
                d_isBaseOutput = true;

                err = TECEND112();
                if (err) cout << "Error in closing Tecplot file" << endl;

            }

#endif

        }
        void Cae_Output::SetTecplotBinary_DomainSolution(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone) {

#ifdef HAVE_TECIO

            double   t;
            INTEGER4 i, N, iVar, err, Debug, NPts, NElm, IsDouble, KMax;
            INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
            INTEGER4 *ShareFromZone = NULL, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
            string buffer, variables;
            stringstream file;
            bool first_zone = true, unsteady = config->GetUnsteady_Simulation(), GridMovement = config->GetGrid_Movement();
            bool Wrt_Unsteady = config->GetWrt_Unsteady();
            unsigned long iExtIter = config->GetExtIter();
            unsigned short NVar, dims = geometry->GetnDim();
            enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
            enum	 ZoneType { ORDERED = 0, FELINESEG = 1, FETRIANGLE = 2, FEQUADRILATERAL = 3, FETETRAHEDRON = 4, FEBRICK = 5, FEPOLYGON = 6, FEPOLYHEDRON = 7 };

            /*--- Consistent data for Tecplot zones ---*/
            Debug = 0;
            IsDouble = 1;
            NPts = (INTEGER4)d_numGlobalPoints;
            t = iExtIter*config->GetDelta_UnstTime();
            KMax = 0;
            ICellMax = 0;
            JCellMax = 0;
            KCellMax = 0;
            StrandID = (INTEGER4)iExtIter + 1;
            ParentZn = 0;
            IsBlock = 1;
            NumFaceConnections = 0;
            FaceNeighborMode = 0;
            ShareConnectivityFromZone = 0;

            file.str(string());
            buffer = config->GetFlow_FileName();

            file << buffer;

            if (unsteady) {
                if (((int)iExtIter >= 0) && ((int)iExtIter < 10))			file << "_0000" << iExtIter;
                if (((int)iExtIter >= 10) && ((int)iExtIter < 100))		file << "_000" << iExtIter;
                if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))		file << "_00" << iExtIter;
                if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))	file << "_0" << iExtIter;
                if ((int)iExtIter >= 10000)							file << iExtIter;
            }
            file << ".sol.plt";
            FileType = SOLUTION;
            variables = AssembleVariableNames(geometry, config, nVar_Consv, &NVar);
            if (config->GetKind_SU2() == SU2_SOL) {
                if (Wrt_Unsteady && GridMovement) nVar_Total = NVar;
                else nVar_Total = NVar + dims;
            }

            /*--- Open Tecplot file ---*/
            err = TECINI112((char *)config->GetFlow_FileName().c_str(),
                (char *)variables.c_str(),
                (char *)file.str().c_str(),
                (char *)".",
                &FileType,
                &Debug,
                &IsDouble);
            if (err) cout << "Error in opening Tecplot file" << endl;

            //  first_zone = true;
            //  ShareFromZone = new INTEGER4[NVar];
            //  for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

            if (d_numGlobalTrias > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FETRIANGLE; NElm = (INTEGER4)d_numGlobalTrias; N = NElm*N_POINTS_TRIANGLE;

                err = TECZNE112((char*)"Triangle Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    ShareFromZone = new INTEGER4[NVar];
                    for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iVar = 0; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        else {
                            for (iVar = dims; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                    }
                    else {

                        if (Wrt_Unsteady && GridMovement) {

                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            if (dims == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                                ShareFromZone[i++] = 1;
                            }
                        }

                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    first_zone = false;
                }

            }
            if (d_numGlobalQuads > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)d_numGlobalQuads; N = NElm*N_POINTS_QUADRILATERAL;

                err = TECZNE112((char*)"Quadrilateral Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    ShareFromZone = new INTEGER4[NVar];
                    for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iVar = 0; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        else {
                            for (iVar = dims; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            if (dims == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                                ShareFromZone[i++] = 1;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }

                    first_zone = false;
                }

            }
            if (d_numGlobalTetrs > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FETETRAHEDRON; NElm = (INTEGER4)d_numGlobalTetrs; N = NElm*N_POINTS_TETRAHEDRON;

                err = TECZNE112((char*)"Tetrahedral Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    ShareFromZone = new INTEGER4[NVar];
                    for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iVar = 0; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        else {
                            for (iVar = dims; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            if (dims == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                                ShareFromZone[i++] = 1;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }

                    first_zone = false;
                }

            }
            if (d_numGlobalHexas > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FEBRICK; NElm = (INTEGER4)d_numGlobalHexas; N = NElm*N_POINTS_HEXAHEDRON;

                err = TECZNE112((char*)"Hexahedral Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    ShareFromZone = new INTEGER4[NVar];
                    for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iVar = 0; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        else {
                            for (iVar = dims; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            if (dims == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                                ShareFromZone[i++] = 1;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }

                    first_zone = false;
                }

            }
            if (d_numGlobalPyras > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FEBRICK; NElm = (INTEGER4)d_numGlobalPyras; N = NElm*N_POINTS_HEXAHEDRON;

                err = TECZNE112((char*)"Pyramid Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    ShareFromZone = new INTEGER4[NVar];
                    for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iVar = 0; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        else {
                            for (iVar = dims; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            if (dims == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                                ShareFromZone[i++] = 1;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }

                    first_zone = false;
                }

            }
            if (d_numGlobalPriss > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FEBRICK; NElm = (INTEGER4)d_numGlobalPriss; N = NElm*N_POINTS_HEXAHEDRON;

                err = TECZNE112((char*)"Prism Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    ShareFromZone = new INTEGER4[NVar];
                    for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iVar = 0; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        else {
                            for (iVar = dims; iVar < nVar_Total; iVar++) {
                                err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                            if (dims == 3) {
                                err = TECDAT112(&NPts, Coords[2], &IsDouble);
                                if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                                ShareFromZone[i++] = 1;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }

                    first_zone = false;
                }
            }

            delete[] ShareFromZone;

            err = TECEND112();
            if (err) cout << "Error in closing Tecplot file" << endl;

#endif

        }

        void Cae_Output::SetTecplotBinary_SurfaceMesh(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone) {

#ifdef HAVE_TECIO

            double   t;
            INTEGER4 i, N, err, Debug, NPts, NElm, IsDouble, KMax;
            INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
            INTEGER4 *ShareFromZone, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
            string buffer, variables;
            stringstream file;
            bool first_zone = true;
            unsigned short iDim, dims = geometry->GetnDim();
            unsigned long iPoint, iElem, iNode;
            enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
            enum	 ZoneType { ORDERED = 0, FELINESEG = 1, FETRIANGLE = 2, FEQUADRILATERAL = 3, FETETRAHEDRON = 4, FEBRICK = 5, FEPOLYGON = 6, FEPOLYHEDRON = 7 };

            /*--- Write Tecplot solution file ---*/
            if (!d_isSurfOutput) {

                file.str(string());
                buffer = config->GetSurfFlowCoeff_FileName();

                file << buffer << ".mesh.plt";
                FileType = GRID;

                if (dims == 2) variables = "x y";
                else if (dims == 3) variables = "x y z";
                else cout << "Error: wrong number of dimensions: " << dims << endl;

                first_zone = true;
                ShareFromZone = new INTEGER4[dims];
                for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

                /*--- Perform a renumbering for the surface points/elements ---*/
                unsigned long *LocalIndex = new unsigned long[d_numGlobalPoints + 1];
                bool *SurfacePoint = new bool[d_numGlobalPoints + 1];

                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) SurfacePoint[iPoint] = false;

                for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                    iNode = iElem*N_POINTS_LINE;
                    SurfacePoint[Conn_Line[iNode + 0]] = true;
                    SurfacePoint[Conn_Line[iNode + 1]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                    iNode = iElem*N_POINTS_TRIANGLE;
                    SurfacePoint[Conn_BoundTria[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundTria[iNode + 2]] = true;
                }
                for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                    iNode = iElem*N_POINTS_QUADRILATERAL;
                    SurfacePoint[Conn_BoundQuad[iNode + 0]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 1]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 2]] = true;
                    SurfacePoint[Conn_BoundQuad[iNode + 3]] = true;
                }

                unsigned long d_numSurfPoints = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                    LocalIndex[iPoint] = 0;
                    if (SurfacePoint[iPoint]) { d_numSurfPoints++; LocalIndex[iPoint] = d_numSurfPoints; }
                }

                /*--- Collect surface coordinates into one array as well ---*/
                /*--- Note the -1 in the Coords/Data array in order to undo the 1-based indexing ---*/
                double **Surf_Coords = new double*[dims];
                for (iDim = 0; iDim < dims; iDim++)
                    Surf_Coords[iDim] = new double[d_numSurfPoints];

                unsigned long iSurf_Poin = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                    if (SurfacePoint[iPoint]) {
                        for (iDim = 0; iDim < dims; iDim++) {
                            if (config->GetKind_SU2() == SU2_SOL)
                                Surf_Coords[iDim][iSurf_Poin] = Data[iDim][iPoint - 1];
                            else
                                Surf_Coords[iDim][iSurf_Poin] = Coords[iDim][iPoint - 1];
                        }
                        iSurf_Poin++;
                    }
                }

                /*--- Consistent data for Tecplot zones ---*/
                Debug = 0;
                IsDouble = 1;
                NPts = (INTEGER4)d_numSurfPoints;
                t = 0.0;//iExtIter*config->GetDelta_UnstTimeND();
                KMax = 0;
                ICellMax = 0;
                JCellMax = 0;
                KCellMax = 0;
                StrandID = 0;//(INTEGER4)iExtIter;
                ParentZn = 0;
                IsBlock = 1;
                NumFaceConnections = 0;
                FaceNeighborMode = 0;
                ShareConnectivityFromZone = 0;

                /*--- Open Tecplot file ---*/
                err = TECINI112((char *)config->GetSurfFlowCoeff_FileName().c_str(),
                    (char *)variables.c_str(),
                    (char *)file.str().c_str(),
                    (char *)".",
                    &FileType,
                    &Debug,
                    &IsDouble);
                if (err) cout << "Error in opening Tecplot file" << endl;


                if (d_numGlobalBoundLines > 0) {

                    /*--- Put the connectivity into a single array for writing ---*/
                    int *Conn_Line_New = new int[d_numGlobalBoundLines*N_POINTS_LINE];
                    iNode = 0;
                    for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                        iNode = iElem*N_POINTS_LINE;
                        Conn_Line_New[iNode + 0] = LocalIndex[Conn_Line[iNode + 0]];
                        Conn_Line_New[iNode + 1] = LocalIndex[Conn_Line[iNode + 1]];
                    }

                    /*--- Write the zone header information ---*/
                    ZoneType = FELINESEG; NElm = (INTEGER4)d_numGlobalBoundLines; N = NElm*N_POINTS_LINE;

                    err = TECZNE112((char*)"Line Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        NULL,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/
                    if (first_zone) {

                        err = TECDAT112(&NPts, Surf_Coords[0], &IsDouble); ShareFromZone[0] = 1;
                        err = TECDAT112(&NPts, Surf_Coords[1], &IsDouble); ShareFromZone[1] = 1;
                        if (geometry->GetnDim() == 3) {
                            err = TECDAT112(&NPts, Surf_Coords[2], &IsDouble);
                            ShareFromZone[2] = 1;
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_Line_New);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                    delete[] Conn_Line_New;
                }

                if (d_numGlobalBoundTrias > 0) {

                    /*--- Put the connectivity into a single array for writing ---*/
                    int *Conn_BoundTria_New = new int[d_numGlobalBoundTrias*N_POINTS_TRIANGLE];

                    iNode = 0;
                    for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                        iNode = iElem*N_POINTS_TRIANGLE;
                        Conn_BoundTria_New[iNode + 0] = LocalIndex[Conn_BoundTria[iNode + 0]];
                        Conn_BoundTria_New[iNode + 1] = LocalIndex[Conn_BoundTria[iNode + 1]];
                        Conn_BoundTria_New[iNode + 2] = LocalIndex[Conn_BoundTria[iNode + 2]];
                    }

                    /*--- Write the zone header information ---*/
                    ZoneType = FETRIANGLE; NElm = (INTEGER4)d_numGlobalBoundTrias; N = NElm*N_POINTS_TRIANGLE;

                    err = TECZNE112((char*)"Triangle Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/
                    if (first_zone) {

                        err = TECDAT112(&NPts, Surf_Coords[0], &IsDouble); ShareFromZone[0] = 1;
                        err = TECDAT112(&NPts, Surf_Coords[1], &IsDouble); ShareFromZone[1] = 1;
                        if (geometry->GetnDim() == 3) {
                            err = TECDAT112(&NPts, Surf_Coords[2], &IsDouble);
                            ShareFromZone[2] = 1;
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_BoundTria_New);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                    delete[] Conn_BoundTria_New;
                }

                if (d_numGlobalBoundQuads > 0) {


                    /*--- Put the connectivity into a single array for writing ---*/
                    int *Conn_BoundQuad_New = new int[d_numGlobalBoundQuads*N_POINTS_QUADRILATERAL];
                    iNode = 0;
                    for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                        iNode = iElem*N_POINTS_QUADRILATERAL;
                        Conn_BoundQuad_New[iNode + 0] = LocalIndex[Conn_BoundQuad[iNode + 0]];
                        Conn_BoundQuad_New[iNode + 1] = LocalIndex[Conn_BoundQuad[iNode + 1]];
                        Conn_BoundQuad_New[iNode + 2] = LocalIndex[Conn_BoundQuad[iNode + 2]];
                        Conn_BoundQuad_New[iNode + 3] = LocalIndex[Conn_BoundQuad[iNode + 3]];
                    }

                    /*--- Write the zone header information ---*/
                    ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)d_numGlobalBoundQuads; N = NElm*N_POINTS_QUADRILATERAL;

                    err = TECZNE112((char*)"Quadrilateral Elements",
                        &ZoneType,
                        &NPts,
                        &NElm,
                        &KMax,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &t,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        0,         /* TotalNumFaceNodes */
                        0,         /* NumConnectedBoundaryFaces */
                        0,         /* TotalNumBoundaryConnections */
                        NULL,      /* PassiveVarList */
                        NULL,      /* ValueLocation */
                        ShareFromZone,      /* ShareVarFromZone */
                        &ShareConnectivityFromZone);
                    if (err) cout << "Error writing Tecplot zone data" << endl;

                    /*--- write node coordinates and data if not done already---*/
                    if (first_zone) {

                        err = TECDAT112(&NPts, Surf_Coords[0], &IsDouble); ShareFromZone[0] = 1;
                        err = TECDAT112(&NPts, Surf_Coords[1], &IsDouble); ShareFromZone[1] = 1;
                        if (geometry->GetnDim() == 3) {
                            err = TECDAT112(&NPts, Surf_Coords[2], &IsDouble);
                            ShareFromZone[2] = 1;
                        }
                        if (err) cout << "Error writing coordinates to Tecplot file" << endl;
                        first_zone = false;
                    }

                    err = TECNOD112(Conn_BoundQuad_New);
                    if (err) cout << "Error writing connectivity to Tecplot file" << endl;

                    delete[] Conn_BoundQuad_New;
                }

                for (iDim = 0; iDim < dims; iDim++)
                    delete[] Surf_Coords[iDim];
                delete[] Surf_Coords;
                delete[] ShareFromZone;
                delete[] LocalIndex;
                delete[] SurfacePoint;
                d_isSurfOutput = true;

                err = TECEND112();
                if (err) cout << "Error in closing Tecplot file" << endl;

            }

#endif

        }
        void Cae_Output::SetTecplotBinary_SurfaceSolution(Config_p config, GEOM::GEOM_Geometry *geometry, unsigned short val_iZone) 
        {
#ifdef HAVE_TECIO

            double   t;
            INTEGER4 i, N, iVar, err, Debug, NPts, NElm, IsDouble, KMax;
            INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
            INTEGER4 *ShareFromZone, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
            string buffer, variables;
            stringstream file;
            bool first_zone = true, unsteady = config->GetUnsteady_Simulation(), GridMovement = config->GetGrid_Movement();
            bool Wrt_Unsteady = config->GetWrt_Unsteady();
            unsigned long iPoint, iElem, iNode, iSurf_Poin, iExtIter = config->GetExtIter();
            unsigned short iDim, NVar, dims = geometry->GetnDim();
            enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
            enum	 ZoneType { ORDERED = 0, FELINESEG = 1, FETRIANGLE = 2, FEQUADRILATERAL = 3, FETETRAHEDRON = 4, FEBRICK = 5, FEPOLYGON = 6, FEPOLYHEDRON = 7 };


            file.str(string());
            buffer = config->GetSurfFlowCoeff_FileName();

            file << buffer;

            if (unsteady) {
                if (((int)iExtIter >= 0) && ((int)iExtIter < 10))			file << "_0000" << iExtIter;
                if (((int)iExtIter >= 10) && ((int)iExtIter < 100))		file << "_000" << iExtIter;
                if (((int)iExtIter >= 100) && ((int)iExtIter < 1000))		file << "_00" << iExtIter;
                if (((int)iExtIter >= 1000) && ((int)iExtIter < 10000))	file << "_0" << iExtIter;
                if ((int)iExtIter >= 10000)							file << iExtIter;
            }
            file << ".sol.plt";
            FileType = SOLUTION;
            variables = AssembleVariableNames(geometry, config, nVar_Consv, &NVar);
            if (config->GetKind_SU2() == SU2_SOL) {
                if (Wrt_Unsteady && GridMovement) nVar_Total = NVar;
                else nVar_Total = NVar + dims;
            }

            first_zone = true;
            ShareFromZone = new INTEGER4[NVar];
            for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;


            /*--- Perform a renumbering for the surface points/elements ---*/
            unsigned long *LocalIndex = new unsigned long[d_numGlobalPoints + 1];
            bool *SurfacePoint = new bool[d_numGlobalPoints + 1];

            for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) SurfacePoint[iPoint] = false;

            for (iElem = 0; iElem < d_numGlobalBoundLines; iElem++) {
                iNode = iElem*N_POINTS_LINE;
                SurfacePoint[Conn_Line[iNode + 0]] = true;
                SurfacePoint[Conn_Line[iNode + 1]] = true;
            }
            for (iElem = 0; iElem < d_numGlobalBoundTrias; iElem++) {
                iNode = iElem*N_POINTS_TRIANGLE;
                SurfacePoint[Conn_BoundTria[iNode + 0]] = true;
                SurfacePoint[Conn_BoundTria[iNode + 1]] = true;
                SurfacePoint[Conn_BoundTria[iNode + 2]] = true;
            }
            for (iElem = 0; iElem < d_numGlobalBoundQuads; iElem++) {
                iNode = iElem*N_POINTS_QUADRILATERAL;
                SurfacePoint[Conn_BoundQuad[iNode + 0]] = true;
                SurfacePoint[Conn_BoundQuad[iNode + 1]] = true;
                SurfacePoint[Conn_BoundQuad[iNode + 2]] = true;
                SurfacePoint[Conn_BoundQuad[iNode + 3]] = true;
            }

            unsigned long d_numSurfPoints = 0;
            for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                LocalIndex[iPoint] = 0;
                if (SurfacePoint[iPoint]) { d_numSurfPoints++; LocalIndex[iPoint] = d_numSurfPoints; }
            }

            /*--- Collect surface coordinates into one array as well ---*/
            /*--- Note the -1 in the Coords/Data array in order to undo the 1-based indexing ---*/
            double **Surf_Coords;
            if (Wrt_Unsteady && GridMovement) {
                Surf_Coords = new double*[dims];
                for (iDim = 0; iDim < dims; iDim++)
                    Surf_Coords[iDim] = new double[d_numSurfPoints];

                iSurf_Poin = 0;
                for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                    if (SurfacePoint[iPoint]) {
                        for (iDim = 0; iDim < dims; iDim++) {
                            if (config->GetKind_SU2() == SU2_SOL)
                                Surf_Coords[iDim][iSurf_Poin] = Data[iDim][iPoint - 1];
                            else
                                Surf_Coords[iDim][iSurf_Poin] = Coords[iDim][iPoint - 1];
                        }
                        iSurf_Poin++;
                    }
                }
            }

            /*--- Collect surface data into one array for the surface as well ---*/
            /*--- Note the -1 in the Coords/Data array in order to undo the 1-based indexing ---*/
            double **Surf_Data = new double*[nVar_Total];
            for (iVar = 0; iVar < nVar_Total; iVar++)
                Surf_Data[iVar] = new double[d_numSurfPoints];

            iSurf_Poin = 0;
            for (iPoint = 0; iPoint < d_numGlobalPoints + 1; iPoint++) {
                if (SurfacePoint[iPoint]) {
                    for (iVar = 0; iVar < nVar_Total; iVar++) {
                        if (config->GetKind_SU2() == SU2_SOL) {
                            if (Wrt_Unsteady && GridMovement)
                                Surf_Data[iVar][iSurf_Poin] = Data[iVar][iPoint - 1];
                            else
                                Surf_Data[iVar][iSurf_Poin] = Data[iVar][iPoint - 1];
                        }
                        else
                            Surf_Data[iVar][iSurf_Poin] = Data[iVar][iPoint - 1];
                    }
                    iSurf_Poin++;
                }
            }

            /*--- Consistent data for Tecplot zones ---*/
            Debug = 0;
            IsDouble = 1;
            NPts = (INTEGER4)d_numSurfPoints;
            t = iExtIter*config->GetDelta_UnstTime();
            KMax = 0;
            ICellMax = 0;
            JCellMax = 0;
            KCellMax = 0;
            StrandID = (INTEGER4)iExtIter + 1;
            ParentZn = 0;
            IsBlock = 1;
            NumFaceConnections = 0;
            FaceNeighborMode = 0;
            ShareConnectivityFromZone = 0;


            /*--- Open Tecplot file ---*/
            err = TECINI112((char *)config->GetFlow_FileName().c_str(),
                (char *)variables.c_str(),
                (char *)file.str().c_str(),
                (char *)".",
                &FileType,
                &Debug,
                &IsDouble);
            if (err) cout << "Error in opening Tecplot file" << endl;


            if (d_numGlobalBoundLines > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FELINESEG; NElm = (INTEGER4)d_numGlobalBoundLines; N = NElm*N_POINTS_LINE;

                err = TECZNE112((char*)"Line Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iDim = 0; iDim < dims; iDim++) {
                                err = TECDAT112(&NPts, Surf_Data[iDim], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        for (iVar = dims; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iDim = 0; iDim < dims; iDim++) {
                                err = TECDAT112(&NPts, Surf_Coords[iDim], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    first_zone = false;
                }

            }

            if (d_numGlobalBoundTrias > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FETRIANGLE; NElm = (INTEGER4)d_numGlobalBoundTrias; N = NElm*N_POINTS_TRIANGLE;

                err = TECZNE112((char*)"Triangle Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iDim = 0; iDim < dims; iDim++) {
                                err = TECDAT112(&NPts, Surf_Data[iDim], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        for (iVar = dims; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iDim = 0; iDim < dims; iDim++) {
                                err = TECDAT112(&NPts, Surf_Coords[iDim], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    first_zone = false;
                }

            }

            if (d_numGlobalBoundQuads > 0) {

                /*--- Write the zone header information ---*/
                ZoneType = FEQUADRILATERAL; NElm = (INTEGER4)d_numGlobalBoundQuads; N = NElm*N_POINTS_QUADRILATERAL;

                err = TECZNE112((char*)"Quadrilateral Elements",
                    &ZoneType,
                    &NPts,
                    &NElm,
                    &KMax,
                    &ICellMax,
                    &JCellMax,
                    &KCellMax,
                    &t,
                    &StrandID,
                    &ParentZn,
                    &IsBlock,
                    &NumFaceConnections,
                    &FaceNeighborMode,
                    0,         /* TotalNumFaceNodes */
                    0,         /* NumConnectedBoundaryFaces */
                    0,         /* TotalNumBoundaryConnections */
                    NULL,      /* PassiveVarList */
                    NULL,      /* ValueLocation */
                    ShareFromZone,      /* ShareVarFromZone */
                    &ShareConnectivityFromZone);
                if (err) cout << "Error writing Tecplot zone data" << endl;

                /*--- write node coordinates and data if not done already---*/
                if (first_zone) {

                    i = 0;
                    if (config->GetKind_SU2() == SU2_SOL) {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iDim = 0; iDim < dims; iDim++) {
                                err = TECDAT112(&NPts, Surf_Data[iDim], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        for (iVar = dims; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    else {
                        if (Wrt_Unsteady && GridMovement) {
                            for (iDim = 0; iDim < dims; iDim++) {
                                err = TECDAT112(&NPts, Surf_Coords[iDim], &IsDouble); ShareFromZone[i++] = 1;
                                if (err) cout << "Error writing data to Tecplot file" << endl;
                            }
                        }
                        for (iVar = 0; iVar < nVar_Total; iVar++) {
                            err = TECDAT112(&NPts, Surf_Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
                            if (err) cout << "Error writing data to Tecplot file" << endl;
                        }
                    }
                    first_zone = false;
                }

            }

            for (iVar = 0; iVar < nVar_Total; iVar++)
                delete[] Surf_Data[iVar];
            delete[] Surf_Data;

            if (Wrt_Unsteady && GridMovement) {
                for (iDim = 0; iDim < dims; iDim++)
                    delete[] Surf_Coords[iDim];
                delete[] Surf_Coords;
            }
            delete[] LocalIndex;
            delete[] SurfacePoint;
            delete[] ShareFromZone;

            err = TECEND112();
            if (err) cout << "Error in closing Tecplot file" << endl;

#endif

        }

        std::string Cae_Output::AssembleVariableNames(Geom_p geometry, Config_p config, unsigned short nVar_Consv, unsigned short *NVar) 
        {
            /*--- Local variables ---*/
            std::stringstream variables; variables.str(string());
            unsigned short iVar;
            *NVar = 0;
            unsigned short iDim, nDim = geometry->GetnDim();
            unsigned short Kind_Solver = config->GetKind_Solver();
            bool grid_movement = config->GetGrid_Movement();
            bool Wrt_Unsteady = config->GetWrt_Unsteady();


            /*--- Write the basic variable header based on the particular solution ----*/

            /*--- Write the list of the fields in the restart file.
            Without including the PointID---*/
            if (config->GetKind_SU2() == SU2_SOL) 
            {
                /*--- If SU2_SOL called this routine, we already have a set of output
                variables with the appropriate string tags stored in the config class.
                We simply read in and remove the quotation marks from the var names. ---*/

                /*--- Set the number of variables to be written. Subtract off an index for
                the PointID as well as each coordinate (x, y, z). ---*/
                string varname;

                if (Wrt_Unsteady && grid_movement) 
                {
                    *NVar = config->fields.size() - 1;
                    for (unsigned short iField = 1; iField < config->fields.size(); iField++)
                    {
                        varname = config->fields[iField];
                        varname.erase(varname.begin(), varname.begin() + 1);
                        varname.erase(varname.end() - 1, varname.end());
                        variables << varname << " ";
                    }
                }
                else 
                {
                    *NVar = config->fields.size() - 1 - nDim;
                    for (unsigned short iField = 1 + nDim; iField < config->fields.size(); iField++) 
                    {
                        varname = config->fields[iField];
                        varname.erase(varname.begin(), varname.begin() + 1);
                        varname.erase(varname.end() - 1, varname.end());
                        variables << varname << " ";
                    }
                }
            }
            else 
            {
                if (Wrt_Unsteady && grid_movement) 
                {
                    if (nDim == 2) 
                    {
                        variables << "x y "; 
                        *NVar += 2;
                    }
                    else 
                    {
                        variables << "x y z "; 
                        *NVar += 3;
                    }
                }

                for (iVar = 0; iVar < nVar_Consv; iVar++) 
                {
                    variables << "Conservative_" << iVar + 1 << " "; 
                    *NVar += 1;
                }

                if (config->GetWrt_Limiters()) 
                {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) 
                    {
                        variables << "Limiter_" << iVar + 1 << " ";
                        *NVar += 1;
                    }
                }

                if (config->GetWrt_Residuals()) 
                {
                    for (iVar = 0; iVar < nVar_Consv; iVar++) 
                    {
                        variables << "Residual_" << iVar + 1 << " ";
                        *NVar += 1;
                    }
                }

                /*--- Add names for any extra variables (this will need to be adjusted). ---*/
                if (grid_movement) 
                {
                    if (nDim == 2)
                    {
                        variables << "Grid_Velx Grid_Vely "; 
                        *NVar += 2;
                    }
                    else 
                    {
                        variables << "Grid_Velx Grid_Vely Grid_Velz "; 
                        *NVar += 3;
                    }
                }

                if (config->GetKind_Regime() == FREESURFACE) 
                {
                    variables << "Density ";
                    *NVar += 1;
                }

                if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
                {
                    variables << "Pressure Temperature Pressure_Coefficient Mach ";
                    *NVar += 4;
                }

                if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) 
                {
                    variables << "Laminar_Viscosity Skin_Friction_Coefficient Heat_Flux Y_Plus ";
                    *NVar += 4;
                }

                if (Kind_Solver == RANS) 
                {
                    variables << "Eddy_Viscosity ";
                    *NVar += 1;
                }

                if (config->GetWrt_SharpEdges())
                {
                    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) 
                    {
                        variables << "Sharp_Edge_Dist ";
                        *NVar += 1;
                    }
                }

                if ((Kind_Solver == TNE2_EULER) || (Kind_Solver == TNE2_NAVIER_STOKES)) 
                {
                    variables << "Mach Pressure Temperature Temperature_ve ";
                    *NVar += 4;
                }

                if (Kind_Solver == TNE2_NAVIER_STOKES) 
                {
                    for (iVar = 0; iVar < config->GetnSpecies(); iVar++)
                        variables << "DiffusionCoeff_" << iVar << " ";
                    variables << "Laminar_Viscosity ThermConductivity ThermConductivity_ve";
                    *NVar += 4;
                }

                if (Kind_Solver == POISSON_EQUATION) 
                {
                    for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
                    {
                        variables << "poissonField_" << iDim + 1 << " ";
                        *NVar += 1;
                    }
                }

                if ((Kind_Solver == ADJ_EULER) ||
                    (Kind_Solver == ADJ_NAVIER_STOKES) ||
                    (Kind_Solver == ADJ_RANS) ||
                    (Kind_Solver == ADJ_TNE2_EULER) ||
                    (Kind_Solver == ADJ_TNE2_NAVIER_STOKES)) 
                {
                    variables << "Surface_Sensitivity Solution_Sensor ";
                    *NVar += 2;
                }
            }

            return variables.str();
        }

        void Cae_Output::SetSU2_MeshASCII(Config_p config, Geom_p geometry)
        {
            char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];
            unsigned long iElem, iPoint, iElem_Bound, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], iNode, nElem;
            unsigned short iMarker, iDim, nDim = geometry->GetnDim(), iChar, iPeriodic, nPeriodic = 0, VTK_Type, nMarker_;
            double *center, *angles, *transl;
            ofstream output_file;
            ifstream input_file;
            string Grid_Marker, text_line, Marker_Tag, str;
            string::size_type position;

            /*--- Read the name of the output and input file ---*/
            str = config->GetMesh_Out_FileName();
            strcpy(out_file, str.c_str());
            strcpy(cstr, out_file);
            output_file.precision(15);
            output_file.open(cstr, ios::out);

            /*--- Write dimensions data. ---*/
            output_file << "NDIME= " << nDim << endl;

            /*--- Write connectivity data. ---*/
            nElem = d_numGlobalTrias + d_numGlobalQuads + d_numGlobalTetrs + d_numGlobalHexas + d_numGlobalPriss + d_numGlobalPyras;

            output_file << "NELEM= " << nElem << endl;

            nElem = 0;

            for (iElem = 0; iElem < d_numGlobalTrias; iElem++) 
            {
                iNode = iElem*N_POINTS_TRIANGLE;
                output_file << "5\t";
                output_file << Conn_Tria[iNode + 0] - 1 << "\t"; output_file << Conn_Tria[iNode + 1] - 1 << "\t";
                output_file << Conn_Tria[iNode + 2] - 1 << "\t";
                output_file << nElem << "\n"; nElem++;
            }

            for (iElem = 0; iElem < d_numGlobalQuads; iElem++) 
            {
                iNode = iElem*N_POINTS_QUADRILATERAL;
                output_file << "9\t";
                output_file << Conn_Quad[iNode + 0] - 1 << "\t"; output_file << Conn_Quad[iNode + 1] - 1 << "\t";
                output_file << Conn_Quad[iNode + 2] - 1 << "\t"; output_file << Conn_Quad[iNode + 3] - 1 << "\t";
                output_file << nElem << "\n"; nElem++;
            }

            for (iElem = 0; iElem < d_numGlobalTetrs; iElem++)
            {
                iNode = iElem*N_POINTS_TETRAHEDRON;
                output_file << "10\t";
                output_file << Conn_Tetr[iNode + 0] - 1 << "\t" << Conn_Tetr[iNode + 1] - 1 << "\t";
                output_file << Conn_Tetr[iNode + 2] - 1 << "\t" << Conn_Tetr[iNode + 3] - 1 << "\t";
                output_file << nElem << "\n"; nElem++;
            }

            for (iElem = 0; iElem < d_numGlobalHexas; iElem++)
            {
                iNode = iElem*N_POINTS_HEXAHEDRON;
                output_file << "12\t";
                output_file << Conn_Hexa[iNode + 0] - 1 << "\t" << Conn_Hexa[iNode + 1] - 1 << "\t";
                output_file << Conn_Hexa[iNode + 2] - 1 << "\t" << Conn_Hexa[iNode + 3] - 1 << "\t";
                output_file << Conn_Hexa[iNode + 4] - 1 << "\t" << Conn_Hexa[iNode + 5] - 1 << "\t";
                output_file << Conn_Hexa[iNode + 6] - 1 << "\t" << Conn_Hexa[iNode + 7] - 1 << "\t";
                output_file << nElem << "\n"; nElem++;
            }

            for (iElem = 0; iElem < d_numGlobalPriss; iElem++) 
            {
                iNode = iElem*N_POINTS_PRISM;
                output_file << "13\t";
                output_file << Conn_Pris[iNode + 0] - 1 << "\t" << Conn_Pris[iNode + 1] - 1 << "\t";
                output_file << Conn_Pris[iNode + 2] - 1 << "\t" << Conn_Pris[iNode + 3] - 1 << "\t";
                output_file << Conn_Pris[iNode + 4] - 1 << "\t" << Conn_Pris[iNode + 5] - 1 << "\t";
                output_file << nElem << "\n"; nElem++;
            }

            for (iElem = 0; iElem < d_numGlobalPyras; iElem++) 
            {
                iNode = iElem*N_POINTS_PYRAMID;
                output_file << "14\t";
                output_file << Conn_Pyra[iNode + 0] - 1 << "\t" << Conn_Pyra[iNode + 1] - 1 << "\t";
                output_file << Conn_Pyra[iNode + 2] - 1 << "\t" << Conn_Pyra[iNode + 3] - 1 << "\t";
                output_file << Conn_Pyra[iNode + 4] - 1 << "\t";
                output_file << nElem << "\n"; nElem++;
            }

            /*--- Write the node coordinates ---*/
            output_file << "NPOIN= " << d_numGlobalPoints << endl;

            for (iPoint = 0; iPoint < d_numGlobalPoints; iPoint++) 
            {
                for (iDim = 0; iDim < nDim; iDim++)
                    output_file << scientific << Coords[iDim][iPoint] << "\t";
                output_file << iPoint << endl;
            }

            /*--- Read the boundary information ---*/
            input_file.open("boundary.su2", ios::out);

            /*--- Read grid file with format SU2 ---*/
            while (getline(input_file, text_line)) 
            {
                /*--- Write the physical boundaries ---*/
                position = text_line.find("NMARK=", 0);
                if (position != string::npos) 
                {
                    text_line.erase(0, 6); nMarker_ = atoi(text_line.c_str());
                    output_file << "NMARK= " << nMarker_ << endl;

                    for (iMarker = 0; iMarker < nMarker_; iMarker++) 
                    {
                        getline(input_file, text_line);
                        text_line.erase(0, 11);
                        string::size_type position;
                        for (iChar = 0; iChar < 20; iChar++) 
                        {
                            position = text_line.find(" ", 0);
                            if (position != string::npos) 
                                text_line.erase(position, 1);
                            position = text_line.find("\r", 0);
                            if (position != string::npos) 
                                text_line.erase(position, 1);
                            position = text_line.find("\n", 0);
                            if (position != string::npos) 
                                text_line.erase(position, 1);
                        }
                        Marker_Tag = text_line.c_str();

                        /*--- Standart physical boundary ---*/
                        getline(input_file, text_line);

                        text_line.erase(0, 13); 
                        nElem_Bound_ = atoi(text_line.c_str());
                        output_file << "MARKER_TAG= " << Marker_Tag << endl;
                        output_file << "MARKER_ELEMS= " << nElem_Bound_ << endl;

                        for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) 
                        {
                            getline(input_file, text_line);
                            istringstream bound_line(text_line);

                            bound_line >> VTK_Type;
                            output_file << VTK_Type;

                            switch (VTK_Type) 
                            {
                            case LINE:
                                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                                output_file << "\t" << vnodes_edge[0] << "\t" << vnodes_edge[1] << endl;
                                break;
                            case TRIANGLE:
                                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                                output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
                                break;
                            case RECTANGLE:
                                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                                output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
                                break;
                            }
                        }
                    }
                }
            }

            input_file.close();

            remove("boundary.su2");

            /*--- Get the total number of periodic transformations ---*/
            nPeriodic = config->GetnPeriodicIndex();
            output_file << "NPERIODIC= " << nPeriodic << endl;

            /*--- From iPeriodic obtain the iMarker ---*/
            for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++)
            {
                /*--- Retrieve the supplied periodic information. ---*/
                center = config->GetPeriodicCenter(iPeriodic);
                angles = config->GetPeriodicRotation(iPeriodic);
                transl = config->GetPeriodicTranslate(iPeriodic);

                output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
                output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
                output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
                output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;

            }

            output_file.close();

        }

        void Cae_Output::SetSU2_MeshBinary(Config_p config, Geom_p geometry)
        { 
        }
    }
}

