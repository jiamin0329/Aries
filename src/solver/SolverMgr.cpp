/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Solver manager
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    01-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#inlucde "SolverMgr.hpp"

namespace ARIES
{
    SolverMgr::SolverMgr()
    {
    }

    SolverMgr::~SovlerMgr()
    {
    }


    SolverMgr::DoSolve()
    {

        int rank = MASTER_NODE;

#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

        /*--- Main external loop of the solver. Within this loop, each iteration ---*/

        if (rank == MASTER_NODE)
            cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

        /*--- This is temporal and just to check. It will have to be added to the regular history file ---*/

        ofstream historyFile_FSI;
        bool writeHistFSI = config_container[ZONE_0]->GetWrite_Conv_FSI();
        if (writeHistFSI && (rank == MASTER_NODE)){
            char cstrFSI[200];
            string filenameHistFSI = config_container[ZONE_0]->GetConv_FileName_FSI();
            strcpy (cstrFSI, filenameHistFSI.data());
            historyFile_FSI.open (cstrFSI);
            historyFile_FSI << "Time,Iteration,Aitken,URes,logResidual,orderMagnResidual" << endl;
            historyFile_FSI.close();
        }

        while (ExtIter < config_container[ZONE_0]->GetnExtIter()) {

            /*--- Perform some external iteration preprocessing. ---*/

            PreprocessExtIter(ExtIter);

            /*--- Perform a single iteration of the chosen PDE solver. ---*/

            if (!fsi){

                /*--- Perform a dynamic mesh update if required. ---*/
                DynamicMeshUpdate(ExtIter);

                /*--- Run a single iteration of the problem (mean flow, wave, heat, ...). ---*/
                Run();

                /*--- Update the solution for dual time stepping strategy ---*/
                Update();
            }
            else{
                Run();      //In the FSIDriver case, mesh and solution updates are already included into the Run function
            }

            /*--- Monitor the computations after each iteration. ---*/
            Monitor(ExtIter);

            /*--- Output the solution in files. ---*/
            Output(ExtIter);

            /*--- If the convergence criteria has been met, terminate the simulation. ---*/

            if (StopCalc) break;

            ExtIter++;

  

        }


    
    }
