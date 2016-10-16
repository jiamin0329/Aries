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
 *    16-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_TIMER_HPP
#define ARIES_TIMER_HPP

/* Aries includes */
#include "Clock.hpp"
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
    /*!
     *  Class Timer holds the exclusive and total start, stop, and elapsed
     *  time for timers instrumented in Aries.  Total time is simply the time
     *  between calls to the start() and stop() functions.  Exclusive time is
     *  applicable if there are nested timers called.
     *  
     *  System and user start and end times are stored as variables of type
     *  clock_t, defined in the sys/times.h include file.  A detailed explanation
     *  of the structures used to store system and user times is given in the
     *  header for the Clock class. This routine simply accesses the functions
     *  specified in that class.
     *  
     *  Wallclock time may be computed by the systems internal clocks which require
     *  an object of type clock_t, or by Aries_MPI::Wtime() if the code is linked
     *  to MPI libraries.
     *  
     *  In addition to running or not running, a timer may be active or inactive.
     *  An inactive timer is one that is created within a program but will never
     *  be turned on or off because it is either not specified as active in
     *  an input file or it was not explicitly made active by the user.  When
     *  a timer is created, it is active by default.
     *  
     *  Note that the constructor is protected so that timer objects can only
     *  be created by the TimerManager class.
     *  
     *  @see TimerManager
     */
    class Timer
    {
        friend class TimerManager;
    public:
        ~Timer();                                                        /**< @brief Empty destructor for Timer class. */
        const std::string& GetName() const { return d_name; }            /**< @brief Return string name for timer. */

        void Start();
        void Stop();

        void BarrierAndStart();  /**< @brief If active, Aries_MPI::getAriesWorld().Barrier() then start the timer. */
        void BarrierAndStop();   /**< @brief If active, Aries_MPI::getAriesWorld().Barrier() then stop the timer. */

        void StartExclusive();   /**< @brief Start exclusive time.. */
        void StopExclusive();    /**< @brief Stop exclusive time. */

        void Reset();            /**< @brief Reset the state of the timing information. */

        //@{
        //! @name Return total system time (between starts and stops)
        double GetTotalSystemTime() const { return d_systemTotal / Clock::GetClockCycle(); }
        double GetTotalUserTime() const { return d_userTotal / Clock::GetClockCycle(); }
        double GetTotalWallclockTime() const { return d_wallclockTotal; }
        //@}

        double GetMaxWallclockTime() const { return d_maxWallclock; }
        
        //@{
        //! @name Return exclusive system time (between starts and stops)
        double GetExclusiveSystemTime() const { return d_systemExclusive / Clock::GetClockCycle(); }
        double GetExclusiveUserTime() const { return d_userExclusive / Clock::GetClockCycle(); }
        double GetExclusiveWallclockTime() const { return d_wallclockExclusive; }
        //@}

        bool IsActive() const { return d_isActive; }           /**< @brief Acessor for d_isActive. */
        bool IsRunning() const { return d_isRunning; }         /**< @brief Acessor for d_isRunning. */

        /*!
         *  @brief Return number of accesses to start()-stop() functions for the timer.
         */
        int GetNumberAccesses() const { return d_accesses; }      

        /*!
         *  @brief Compute load balance efficiency based on wallclock (non-exclusive) time.
         */  
        double ComputeLoadBalanceEfficiency();

        /*!
         *  @brief Compute max wallclock time based on total (non-exclusive) time.
         */
        void ComputeMaxWallclock();

        /*!
         *  @brief Write timer data members to restart database.
         */
        void PutToRestart() const;

        /*!
         *  @brief Read restarted times from restart database.
         */
        void GetFromRestart();

    protected:
        /*!
         *  @brief The constructor for the Timer class sets timer name string
         *  and integer identifiers, and initializes the timer state.
         */
        explicit Timer(const std::string& name);

        /*!
         *  @brief Set this timer object to be a active or inactive.
         */
        void SetActive(bool is_active) {  d_isActive = is_active; }

        /*!
         *  @brief Add Timer that running concurrently with this one.
         */
        void AddConcurrentTimer(const Timer& timer)
        {
            if (!IsConcurrentTimer(timer))
            {
                d_concurrentTimers.push_back(&timer);
            }
        }

        bool IsConcurrentTimer(const Timer& timer) const;

    private:
        Timer();                               /**< @brief Unimplemented default constructor. */
        Timer(const Timer& other);             /**< @brief Unimplemented copy constructor. */
        Timer& operator = (const Timer& rhs);  /**< @brief Unimplemented assignment operator. */

        std::string d_name;                                /**< @brief Timer name. */
        std::vector<const Timer *> d_concurrentTimers;     /**< @brief Concurrent timer flag. */

        bool d_isRunning;
        bool d_isActive;

        double d_userTotal;           /**< @brief Total times. */
        double d_systemTotal;         /**< @brief Total times. */
        double d_wallclockTotal;      /**< @brief Total times. */

        double d_userExclusive;           /**< @brief Exclusive times. */
        double d_systemExclusive;         /**< @brief Exclusive times. */
        double d_wallclockExclusive;      /**< @brief Exclusive times. */

        double d_maxWallclock;     /**< @brief Cross processor times. */

        //@{
        //! @name Timestamps.
        clock_t d_userStartTotal ;
        clock_t d_userStopTotal;
        clock_t d_userStartExclusive;
        clock_t d_userStopExclusive;
        
        clock_t d_systemStartTotal;
        clock_t d_systemStopTotal;
        clock_t d_systemStartExclusive;
        clock_t d_systemStopExclusive;
        
        double d_wallclockStartTotal;
        double d_wallclockStopTotal;
        double d_wallclockStartExclusive;
        double d_wallclockStopExclusive;
        //@}

        int d_accesses;  /**< @brief Counter of number of times timers start/stop are accessed. */

        static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT;
        static const int ARIES_TIMER_VERSION;
        
    }; // end class Timer
} // end namespace ARIES

#endif

