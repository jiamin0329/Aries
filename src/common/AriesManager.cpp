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
 *================================================================================
 */

#include "AriesManager.hpp"

// Aries head files
#include "StartupShutdownManager.hpp"
#include "Utilities.hpp"
//*

// c++ head files

//*

namespace ARIES
{
    bool AriesManager::d_initialized = false;
    bool AriesManager::d_started = false;

    void AriesManager::Initialize()
    {
        ARIES_ASSERT(!d_initialized);
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
        d_initialized = false;
    }
    
} // end namespace ARIES

