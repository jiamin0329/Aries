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

namespace ARIES
{
    /*
     * Utility for managing startup and shutdown of Aries objects
     * and for managing the maximum number of patch data components allowed.
     *
     * The startup/shutdown mechanism in AriesManager is used to manage the
     * allocation and deallocation of certain memory, particularly static data
     * that must exist during the full extent of a run or for the full extent
     * of a single problem within a run.  The specific data that is controlled
     * by this mechanism is managed using the StartupShutdownManager.
     *
     * The four steps of the startup/shutdown mechanism are:
     *
     * <ul>
     * <li> initialize -- called at the start of a program after MPI is
     *      initialized but befor any other Aries objects are used.
     * <li> startup -- called to begin a problem-specific segment of the code.
     * <li> shutdown -- called at the end of a problem-specific segment of the
     *                  code.  Shuts down and deallocates everything that was
     *                  started and allocated by startup.
     * <li> finalize -- called at the end of a program right before MPI is
     *                  finalized.
     * </ul>
     *
     * The startup and shutdown functions may be called multiple times within a
     * run, in order to allow for the execution of more than one problem within one
     * program. Initialize and finalize must be called exactly once in a single
     * run.
     *
     * @see StartupShutdownManager
     * @see ParallelIO
     */
    class AriesManager
    {
    public:
        /*!
         * @brief Initial setup of the Aries package.
         *
         * This function should be invoked ONLY ONCE at the start of a process
         * to initialize Aries and AFTER MPI is initialized (if used) by a
         * call to one of the AriesMPI init routines.
         *
         * This function initializes Aries I/O, as well as data for any classes
         * that implement the initialize callback interface through StartupShutdownManager.
         *
         * @pre !s_initialized
         */
        static void Initialize();

        /*!
         * @brief Startup of the Aries package.
         *
         * This function invokes startup for any classes that implement the
         * startup callback interface through StartupShutdownManager.
         * This function may be invoked more than once in a process 
         */
        static void Startup();

        /*!
         * @brief Shutdown the Aries package.
         *
         * This function invokes shutdown for any classes that implement the
         * startup callback interface through StartupShutdownManager.
         * This function may be invoked more than once in an process if
         * solving multiple Aries problems.
         */
        static void Shutdown();

        /*!
         * @brief Final cleanup of the Aries package.
         *
         * This function should be invoked ONLY ONCE at the end of a process
         * to complete the cleanup of Aries memory allocations and
         * any other cleanup tasks.  Aries I/O will be finalized, as well as data for any classes that implement
         * the finalize callback interface through StartupShutdownManager.
         *
         * After this function is called, the only thing that should occur before
         * exiting the program is a call to AriesMPI::Finalize().
         *
         * This function should be invoked only once.
         */
        static void Finalize();

        static bool IsInitialized() { return d_initialized; }
        static bool IsStarted() { return d_started; }

    private:
        // Unimplemented default constructor.
        AriesManager();

        // Unimplemented copy constructor.
        AriesManager(const AriesManager& other);

        // Unimplemented assignment operator.
        AriesManager& operator = (const AriesManager& rhs);

        static bool d_initialized;  // Flag indicating AriesManager has been initialized.
        static bool d_started;      // Flag indicating startup has occured.
    }; // end class AriesManager
} // end namespace ARIES

#endif




