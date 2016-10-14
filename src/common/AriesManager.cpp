/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Aries class to manage package startup and shutdown
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *    13-Oct-2016     Jiamin Xu               Add AriesMPI Initialize/Finalize
 *================================================================================
 */

#include "AriesManager.hpp"

/* Aries includes */
#include "AriesMPI.hpp"
#include "StartupShutdownManager.hpp"
#include "Utilities.hpp"

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    bool AriesManager::d_initialized = false;
    bool AriesManager::d_started = false;

    void AriesManager::Initialize(int argc, char *argv[])
    {
        ARIES_ASSERT(!d_initialized);
        AriesMPI::Init(&argc, &argv);
        StartupShutdownManager::Initialize();
        d_initialized = true;
    }

    void AriesManager::Startup()
    {
        ARIES_ASSERT(d_initialized);
        ARIES_ASSERT(!d_started);
        StartupShutdownManager::Startup();
        d_started = true;
    }

    void AriesManager::Shutdown()
    {
        ARIES_ASSERT(d_initialized);
        ARIES_ASSERT(d_started);
        StartupShutdownManager::Shutdown();
        d_started = false;
    }

    void AriesManager::Finalize()
    {
        ARIES_ASSERT(d_initialized);
        StartupShutdownManager::Finalize();
        ARIES::AriesMPI::Finalize(); 
        d_initialized = false;
    }
    
} // end namespace ARIES

