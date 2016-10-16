/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Main entrance for ARIES suite
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    13-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

/* Aries includes */
#include "const_def.h"
#include "IProcData.hpp"
#include "ProcData.hpp"
#include "AriesMPI.hpp"
#include "AriesManager.hpp"
#include "Logger.hpp"
#include "TimerManager.hpp"
/* C++ includes */
#include <string>


int main(int argc, char *argv[])
{
    //ARIES::AriesMPI::Init(&argc, &argv); 
    ARIES::AriesManager::Initialize(argc, argv);
    ARIES::AriesManager::Startup();
    
    ARIES::Logger::GetInstance()->Startup(true,true,"./hello.log");
    ARIES::TimerManager::GetInstance()->Print(std::cout);

    //unsigned short nZone, nDim;
    //bool isFsi;
    //char configFileName[MAX_STRING_SIZE];
    //
    ///*
    // *  Load in the number of zones and spatial dimensions in the mesh
    // *  file (If no config file is specified, defau lt.cfg is used)
    // */
    //if (argc == 2)
    //    strcpy(configFileName, argv[1]);
    //else
    //    strcpy(configFileName, "default.cfg");
    //
    ///*
    // *  Read the name and format of the input mesh file to get from the mesh
    // *  file the number of zones and dimensions from the numerical grid (required
    // *  for variables allocation)
    // */
    //ARIES::IProcData* procData = new ARIES::ProcData();
    //nDim  = procData->GetNumDim();
    //nZone = procData->GetNumZone();
    //isFsi = procData->IsFsiSimulation();



    ///*--- First, given the basic information about the number of zones and the
    //  solver types from the config, instantiate the appropriate driver for the problem
    //  and perform all the preprocessing. ---*/
    //if (nZone == SINGLE_ZONE)
    //{
    //    /*--- Single zone problem: instantiate the single zone driver class. ---*/
    //    driver = new CSingleZoneDriver(config_file_name, nZone, nDim);
    //
    //}
    //else if (config->GetUnsteady_Simulation() == TIME_SPECTRAL)
    //{
    //    /*--- Use the spectral method driver. ---*/
    //    driver = new CSpectralDriver(config_file_name, nZone, nDim);
    //}
    //else if ((nZone == 2) && fsi)
    //{
    //    /*--- FSI problem: instantiate the FSI driver class. ---*/
    //    driver = new CFSIDriver(config_file_name, nZone, nDim);
    //
    //}
    //else
    //{
    //    /*--- Multi-zone problem: instantiate the multi-zone driver class by default
    //      or a specialized driver class for a particular multi-physics problem. ---*/
    //    driver = new CMultiZoneDriver(config_file_name, nZone, nDim);
    //    /*--- Future multi-zone drivers instatiated here. ---*/
    //}
    //
    //delete config;
    //config = NULL;
    //
    ///*--- Launch the main external loop of the solver ---*/
    //driver->StartSolver();
    //
    ///*--- Postprocess all the containers, close history file, exit SU2 ---*/
    //driver->Postprocessing();
    //
    //if (driver != NULL)
    //    delete driver;
    //driver = NULL;
    
    ARIES::Logger::GetInstance()->Shutdown();
    ARIES::AriesManager::Shutdown();
    ARIES::AriesManager::Finalize();
    
    return 0;
}
