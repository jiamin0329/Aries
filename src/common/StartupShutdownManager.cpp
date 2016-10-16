/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Manager for startup and shutdown routines to be called at program
 *    start and exit
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "StartupShutdownManager.hpp"

/* Aries includes */
#include "Utilities.hpp"
/* C++ includes */
#include <cstdlib>

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    bool StartupShutdownManager::d_singletonInitialized = false;

    bool StartupShutdownManager::d_inInitialize = false;
    bool StartupShutdownManager::d_inStartup = false;
    bool StartupShutdownManager::d_inShutdown = false;
    bool StartupShutdownManager::d_inFinalize = false;

    bool StartupShutdownManager::d_initialized = false;
    bool StartupShutdownManager::d_startuped = false;
    bool StartupShutdownManager::d_shutdowned = false;
    bool StartupShutdownManager::d_finalized = false;

    StartupShutdownManager::ListElement* StartupShutdownManager::d_managerList[d_numberOfPriorities];
    StartupShutdownManager::ListElement* StartupShutdownManager::d_managerListLast[d_numberOfPriorities];
    int StartupShutdownManager::d_numManagerItems[d_numberOfPriorities];

    StartupShutdownManager::Handler::Handler( void(*initialize)(), void(*startup)(), void(*shutdown)(), void(*finalize)(), unsigned char priority):
            d_initialize(initialize),
            d_startup(startup),
            d_shutdown(shutdown),
            d_finalize(finalize),
            d_priority(priority)
    {
        StartupShutdownManager::RegisterHandler(this);
    }

    StartupShutdownManager::Handler::~Handler()
    {
    }

    void StartupShutdownManager::Handler::Initialize()
    {
        if (d_initialize)
            (*d_initialize)();
    }

    void StartupShutdownManager::Handler::Startup()
    {
        if (d_startup)
            (*d_startup)();
    }

    void StartupShutdownManager::Handler::Shutdown()
    {
        if (d_shutdown)
            (*d_shutdown)();
    }

    void StartupShutdownManager::Handler::Finalize()
    {
        if (d_finalize)
            (*d_finalize)();
    }

    unsigned char StartupShutdownManager::Handler::GetPriority()
    {
        return d_priority;
    }

    bool StartupShutdownManager::Handler::HasInitialize()
    {
        return d_initialize != 0;
    }

    bool StartupShutdownManager::Handler::HasStartup()
    {
        return d_startup != 0;
    }

    bool StartupShutdownManager::Handler::HasShutdown()
    {
        return d_shutdown != 0;
    }

    bool StartupShutdownManager::Handler::HasFinalize()
    {
        return d_finalize != 0;
    }

    void StartupShutdownManager::RegisterHandler(IHandler* handler)
    {
        ARIES_ASSERT(handler);
        
        ARIES_ASSERT(!(d_inInitialize && handler->HasInitialize()));
        ARIES_ASSERT(!(d_inStartup && handler->HasStartup()));
        ARIES_ASSERT(!(d_inShutdown && handler->HasShutdown()));
        ARIES_ASSERT(!d_inFinalize);

        if (!d_singletonInitialized)
            SetupSingleton();
        
        ListElement* item = new ListElement;
        item->handler = handler;

        unsigned char priority = handler->GetPriority();

        item->next = 0;
        if (d_numManagerItems[priority] == 0)
            d_managerList[priority] = item;
        else
            d_managerListLast[priority]->next = item;
        
        d_managerListLast[priority] = item;
        ++d_numManagerItems[priority];
    }

    void StartupShutdownManager::Initialize()
    {
        ARIES_ASSERT(!d_initialized);

        d_initialized = true;
        // only shutdown if something was registered
        if (d_singletonInitialized)
        {
            d_inInitialize = true;

            for (int priority = 0; priority < d_numberOfPriorities; ++priority)
            {
                ListElement* item = d_managerList[priority];
                while (item)
                {
                    if (item->handler)
                    {
                        item->handler->Initialize();
                    }
                    item = item->next;
                }
            }
            d_inInitialize = false;
        }
    }

    void StartupShutdownManager::Startup()
    {
        // If this is thrown you need to make sure SAMRAIManger::initialize
        // is called before startup.
        ARIES_ASSERT(d_initialized);
        ARIES_ASSERT(!d_startuped);

        d_startuped = true;

        // only shutdown if something was registered
        if (d_singletonInitialized)
        {
            d_inStartup = true;

            for (int priority = 0; priority < d_numberOfPriorities; ++priority)
            {
                ListElement* item = d_managerList[priority];
                while (item)
                {
                    if (item->handler)
                    {
                        item->handler->Startup();
                    }
                    item = item->next;
                }
            }
            d_inStartup = false;
        }
        d_shutdowned = false;
    }

    void StartupShutdownManager::Shutdown()
    {
        ARIES_ASSERT(d_initialized);
        ARIES_ASSERT(d_startuped);
        ARIES_ASSERT(!d_shutdowned);

        d_shutdowned = true;

        // only shutdown if something was registered
        if (d_singletonInitialized)
        {
            d_inShutdown = true;

            for (int priority = d_numberOfPriorities - 1; priority > -1; --priority)
            {
                ListElement* item = d_managerList[priority];
                while (item)
                {
                    if (item->handler)
                    {
                        item->handler->Shutdown();
                    }
                    item = item->next;
                }
            }
            d_inShutdown = false;
        }
        d_startuped = false;
    }

    void StartupShutdownManager::Finalize()
    {
        ARIES_ASSERT(d_initialized);
        ARIES_ASSERT(d_shutdowned);
        ARIES_ASSERT(!d_finalized);

        d_finalized = true;

        // only finalize if something was registered
        if (d_singletonInitialized)
        {
            d_inFinalize = true;

            for (int priority = d_numberOfPriorities - 1; priority > -1; --priority)
            {
                ListElement* item = d_managerList[priority];
                while (item)
                {
                    if (item->handler)
                    {
                        item->handler->Finalize();
                    }
                    item = item->next;
                }
            }

            for (int priority = 0; priority < d_numberOfPriorities; ++priority)
            {
                ListElement* item = d_managerList[priority];
                while (item)
                {
                    ListElement* to_delete = item;
                    item = item->next;
                    delete to_delete;
                }
            }
            d_inFinalize = false;
        }
        d_initialized = false;
    }
    
    void StartupShutdownManager::SetupSingleton()
    {
        for (int priority = d_numberOfPriorities - 1; priority > -1; --priority)
        {
            d_managerList[priority] = 0;
            d_managerListLast[priority] = 0;
            d_numManagerItems[priority] = 0;
        }

        d_singletonInitialized = true;
    }

    StartupShutdownManager::ListElement::ListElement():
            handler(0),
            next(0)
    {
    }

    StartupShutdownManager::ListElement::~ListElement()
    {
    }

} // end namespace ARIES

