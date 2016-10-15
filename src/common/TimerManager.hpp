/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for processor data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    11-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_TIMERMANAGER_HPP
#define ARIES_TIMERMANAGER_HPP

/* Aries includes */
#include "Timer.hpp"
/* C++ includes */
#include <iostream>
#include <string>
#include <vector>
#include <list>

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    class TimerManager
    {
        friend class Timer;
    public:
        static TimerManager* GetInstance();

        Timer* GetTimer(const std::string& name, bool ignore_timer_input = false);

        bool CheckTimerExists(Timer* timer, const std::string& name) const;
        bool CheckTimerRunning(const std::string& name) const;

        void ResetAllTimers();
        void Print(std::ostream& os = plog);  /**< @brief Print timer information from each processor. */

    protected:
        explicit TimerManager();
        ~TimerManager();

        void RegisterSingletonSubclassInstance(TimerManager* subclass_instance);

        void StartTime(Timer* timer);
        void StopTime(Timer* timer);

    private:
        TimerManager();                                                                /* Default constructor. */
        TimerManager(const TimerManager& other);                                       /* Copy constructor. */
        TimerManager& operator = (const TimerManager& rhs);                            /* Assignment operator. */

        void ActivateExistingTimers();                                                 /* Activate any existing timers that have already been registered. */
        
        bool CheckTimerExistsInArray(Timer* timer,
                                     const std::string& name,
                                     const std::vector<Timer*>& timer_array) const;    /* Add timer to either the active or inactive timer array. */


        void PrintTable(const std::string& table_title,
                        const std::vector<std::string>& column_titles,
                        const std::vector<std::string>& timer_names,
                        const int column_ids[],
                        const double timer_values[][18], std::ostream& os);            /* Print a table of values, using values specified in the timer_values array. */
                                                                                          
        void PrintTable(const std::string& table_title,
                        const std::vector<std::string>& column_titles,
                        const std::vector<std::string>& timer_names,
                        const int max_processor_id[][2],
                        const int max_array_id,
                        const int column_ids[],
                        const double timer_values[][18], std::ostream& os);            /* Print a table of values, using values specified in the timer_values array. */         

        void PrintOverhead(const std::vector<std::string>& timer_names,
                           const double timer_values[][18], std::ostream& os);         /* Output overhead stats for Timers. */

        void PrintConcurrent(std::ostream& os);                                        /* Output concurrent tree of Timers. */           

        void BuildTimerArrays(double timer_values[][18],
                              int max_processor_id[][2],
                              std::vector<std::string>&timer_names);                   /* Build the timer_names, timer_values, and max_processor_id arrays. */

        void BuildOrderedList(const double timer_values[][18],
                              const int column,
                              int index[],
                              const int array_size);                                   /* Build an ordered list array, organizing timers largest to smallest. */
 
        bool CheckTimerInNameLists(const std::string& name);                           /* Checks timer name to determine if it is specified to be turned on. */

        void CheckConsistencyAcrossProcessors();                                       /* Evaluate consistency of timer database across processors. */

        void GetFromInput(const boost::shared_ptr<Database>& input_db);                /* parse input data for managing timers. */

        void AddTimerToNameLists(const std::string& name);                             /* add a timer name to the d_package, d_class, or d_class_method lists. */

        static void Quicksort(const std::vector<double>&a, int index[],
                              int lo, int hi);                                         /* Quicksort algorithm specialized for timer array. */

        static int    ComputePercentageInt   (const double frac, const double tot);    /* Simple methods to compute percentages */
        static double ComputePercentageDouble(const double frac, const double tot);    /* Simple methods to compute percentages */

        void ComputeOverheadConstants();                                               /* Compute the overhead costs of the timing routines for active and non-active timers. */
        double ComputeOverheadConstantActiveOrInactive(bool active);

        void ClearArrays();                                                            /* Clear the registered timers. */

        static void FinalizeCallback();                                                /* Deallocate the TimerManager instance. */

        static TimerManager* d_instance;                                               /* Static data members to manage the singleton timer manager instance. */

        static int d_mainTimerIdentifier;                                              /* Static constants used by timer manager. */
        static int d_inactiveTimerIdentifier;                                          /* Static constants used by timer manager. */

        double d_timerActiveAccessTime;                                                /* Timer accesss overheads. */
        double d_timerInactiveAccessTime;                                              /* Timer accesss overheads. */

        Timer* d_mainTimer;                                                            /* Main timer used to time overall run time. */
        std::vector<Timer*> d_timers;                                                  /* An array of pointers to those timers. */
        std::vector<Timer*> d_inactiveTimers;                                          /* An array of dummy inactive timers is used to record number of accesses to non-active timers. */

        std::list<Timer *> d_exclusiveTimerStack;

        std::list<std::string> d_packageNames;
        std::list<std::string> d_classNames;
        std::list<std::string> d_classMethodNames;

        int d_lengthPackageNames;
        int d_lengthClassNames;
        int d_lengthClassMethodNames;

        double d_printThreshold;

        bool d_printExclusive;
        bool d_printTotal;

        bool d_printProcessor;
        bool d_printMax;
        bool d_printSummed;

        bool d_printUser;
        bool d_printSys;
        bool d_printWall;

        bool d_printPercentage;
        bool d_printConcurrent;
        bool d_printTimerOverhead;

        static const int DEFAULT_NUMBER_OF_TIMERS_INCREMENT = 128;
        static StartupShutdownManager::Handler d_finalizeHandler;
        
    }; // end class TimerManager
} // end namespace ARIES

#endif

