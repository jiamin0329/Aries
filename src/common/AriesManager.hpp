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
#ifndef ARIES_ARIESMANAGER_HPP
#define ARIES_ARIESMANAGER_HPP

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    /*!
     *  @brief Utility for managing startup and shutdown of Aries objects and
     *         for managing the maximum number of patch data components allowed.
     *
     *  The startup/shutdown mechanism in AriesManager is used to manage the
     *  allocation and deallocation of certain memory, particularly static data
     *  that must exist during the full extent of a run or for the full extent
     *  of a single problem within a run. The specific data that is controlled
     *  by this mechanism is managed using the StartupShutdownManager.
     *
     *  The startup and shutdown functions may be called multiple times within a
     *  run, in order to allow for the execution of more than one problem within one
     *  program. Initialize and finalize must be called exactly once in a single
     *  run.
     *
     *  @see AriesMPI
     *  @see StartupShutdownManager
     */
    class AriesManager
    {
    public:
        /*!
         *  @brief Initial setup of the Aries package.
         *
         *  This function should be invoked ONLY ONCE at the start of a process
         *  to initialize Aries.
         *
         *  This function initializes AriesMPI, as well as data for any classes that
         *  implement the initialize callback interface through StartupShutdownManager.
         */
        static void Initialize(int argc, char *argv[]);

        /*!
         *  @brief Startup of the Aries package.
         *
         *  This function invokes startup for any classes that implement the
         *  startup callback interface through StartupShutdownManager.
         *  This function may be invoked more than once in a process. 
         */
        static void Startup();

        /*!
         *  @brief Shutdown the Aries package.
         *
         *  This function invokes shutdown for any classes that implement the
         *  startup callback interface through StartupShutdownManager.
         *  This function may be invoked more than once in an process if
         *  solving multiple Aries problems.
         */
        static void Shutdown();

        /*!
         *  @brief Final cleanup of the Aries package.
         *
         *  This function should be invoked ONLY ONCE at the end of a process
         *  to complete the cleanup of Aries memory allocations and
         *  any other cleanup tasks. Data for any classes that implement
         *  the finalize callback interface will be finalized through StartupShutdownManager.
         */
        static void Finalize();
 
        static bool IsInitialized() { return d_initialized; }      /**< @brief Acessor of d_initialized. */
        static bool IsStarted() { return d_started; }              /**< @brief Acessor of d_started. */     

    private:
        AriesManager();                                            /**< @brief Unimplemented default constructor. */
        AriesManager(const AriesManager& other);                   /**< @brief Unimplemented copy constructor. */
        AriesManager& operator = (const AriesManager& rhs);        /**< @brief Unimplemented assignment operator. */
        
        static bool d_initialized;                                 /**< @brief Flag indicating AriesManager has been initialized. */
        static bool d_started;                                     /**< @brief Flag indicating Startup has occured. */
        
    }; // end class AriesManager
} // end namespace ARIES

#endif

