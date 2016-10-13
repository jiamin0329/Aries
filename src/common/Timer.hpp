/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Timer class to track elapsed time in portions of a program.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_TIMER_HPP
#define ARIES_TIMER_HPP

/* Aries includes */
#include "SAMRAI/tbox/Clock.h"
/* C++ includes */
#include <string>
#include <vector>

/*
 *================================================================================
 *    Forward declarations
 *================================================================================
 */
namespace ARIES
{
    class TimerManager;
}

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    /**
     * Class Timer holds the exclusive and total start, stop, and elapsed
     * time for timers instrumented in SAMRAI.  Total time is simply the time
     * between calls to the start() and stop() functions.  Exclusive time is
     * applicable if there are nested timers called.
     *
     * System and user start and end times are stored as variables of type
     * clock_t, defined in the sys/times.h include file.  A detailed explanation
     * of the structures used to store system and user times is given in the
     * header for the Clock class. This routine simply accesses the functions
     * specified in that class.
     *
     * Wallclock time may be computed by the systems internal clocks which require
     * an object of type clock_t, or by SAMRAI_MPI::Wtime() if the code is linked
     * to MPI libraries.
     *
     * In addition to running or not running, a timer may be active or inactive.
     * An inactive timer is one that is created within a program but will never
     * be turned on or off because it is either not specified as active in
     * an input file or it was not explicitly made active by the user.  When
     * a timer is created, it is active by default.
     *
     * Note that the constructor is protected so that timer objects can only
     * be created by the TimerManager class.
     *
     * @see TimerManager
     */
    class Timer
    {
        friend class TimerManager;
    public:
        /**
         * Empty destructor for Timer class.
         */
        ~Timer();

        /**
         * Return string name for timer.
         */
        const std::string& getName() const
        {
            return d_name;
        }

        void start();
        void stop();

        //If active, SAMRAI_MPI::getSAMRAIWorld().Barrier() then start the timer.
        void barrierAndStart();

        // If active, SAMRAI_MPI::getSAMRAIWorld().Barrier() then stop the timer.
        void barrierAndStop();

        /**
         * Start exclusive time.
         */
        void startExclusive();

        /**
         * Stop exclusive time.
         */
        void stopExclusive();

        /**
         * Reset the state of the timing information.
         */
        void reset();

        /**
         * Return total system time (between starts and stops)
         */
        double getTotalSystemTime() const
        {
            return d_system_total / Clock::getClockCycle();
        }


        double getTotalUserTime() const
        {
            return d_user_total / Clock::getClockCycle();
        }

        double getTotalWallclockTime() const
        {
            return d_wallclock_total;
        }

        double getMaxWallclockTime() const
        {
            return d_max_wallclock;
        }

        double getExclusiveSystemTime() const
        {
            return d_system_exclusive / Clock::getClockCycle();
        }

        /**
         * Return exclusive user time.
         */
        double getExclusiveUserTime() const
        {
            return d_user_exclusive / Clock::getClockCycle();
        }
        
        double getExclusiveWallclockTime() const
        {
            return d_wallclock_exclusive;
        }

        bool isActive() const
        {
            return d_is_active;
        }

        bool isRunning() const
        {

            return d_is_running;
        }

        /**
         * Return number of accesses to start()-stop() functions for the
         * timer.
         */
        int getNumberAccesses() const
        {
            return d_accesses;
        }

        /**
         * Compute load balance efficiency based on wallclock (non-exclusive)
         * time.
         */
        double computeLoadBalanceEfficiency();

        /**
         * Compute max wallclock time based on total (non-exclusive) time.
         */
        void computeMaxWallclock();

        /**
         * Write timer data members to restart database.
         *
         * @pre restart_db
         */
        void putToRestart(const boost::shared_ptr<Database>& restart_db) const;

        /**
         * Read restarted times from restart database.
         *
         * @pre restart_db
         */
        void getFromRestart(const boost::shared_ptr<Database>& restart_db);

    protected:
        /**
         * The constructor for the Timer class sets timer name string
         * and integer identifiers, and initializes the timer state.
         */
        explicit Timer(const std::string& name);

        /*
         * Set this timer object to be a active or inactive.
         */
        void setActive(bool is_active)
        {
            d_is_active = is_active;
        }

        /**
         * Add Timer that running concurrently with this one.
         */
        void addConcurrentTimer(const Timer& timer)
        {
            if (!isConcurrentTimer(timer))
            {
                d_concurrent_timers.push_back(&timer);
            }
        }

        bool isConcurrentTimer(const Timer& timer) const;

    private:
        // Unimplemented default constructor.
        Timer();

        // Unimplemented copy constructor.
        Timer(const Timer& other);

        // Unimplemented assignment operator.
        Timer& operator = (const Timer& rhs);

        /*
         * Class name, id, and concurrent timer flag.
         */
        std::string d_name;
        std::vector<const Timer *> d_concurrent_timers;

        bool d_is_running;
        bool d_is_active;

        /*
         *  Total times (non-exclusive)
         */
        double d_user_total;
        double d_system_total;
        double d_wallclock_total;

        /*
         *  Exclusive times
         */
        double d_user_exclusive;
        double d_system_exclusive;
        double d_wallclock_exclusive;

        /*
         *  Cross processor times (i.e. determined across processors)
         */
        double d_max_wallclock;

        /*
         *  Timestamps.  User and system times are stored as type clock_t.
         *  Wallclock time is also stored as clock_t unless the library has
         * been compiled with MPI.  In this case, the wall time is stored
         * as type double.
         */
        clock_t d_user_start_total;
        clock_t d_user_stop_total;
        clock_t d_system_start_total;
        clock_t d_system_stop_total;
        clock_t d_user_start_exclusive;
        clock_t d_user_stop_exclusive;
        clock_t d_system_start_exclusive;
        clock_t d_system_stop_exclusive;
        double d_wallclock_start_total;
        double d_wallclock_stop_total;
        double d_wallclock_start_exclusive;
        double d_wallclock_stop_exclusive;

        /*
         * Counter of number of times timers start/stop
         * are accessed.
         */
        int d_accesses;

        static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT;
        
        /*
         * Static integer constant describing this class's version number.
         */
        static const int ARIES_TIMER_VERSION;
    }; // end class Timer
} // end namespace ARIES

#endif
