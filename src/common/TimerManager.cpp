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

#include "TimerManager.hpp"

/* Aries includes */
#include "AriesMPI.hpp"
#include "AriesManager.hpp"
#include "StartupShutdownManager.hpp"
#include "Utilities.hpp"
#include "Logger.hpp"

/* C++ includes */
#include <string>
#include <iomanip>


#ifdef __INTEL_COMPILER
#pragma warning (disable:869)
#endif

namespace ARIES
{
    TimerManager* TimerManager::d_instance = NULL;

    int TimerManager::d_mainTimerIdentifier = -1;
    int TimerManager::d_inactiveTimerIdentifier = -9999;

    StartupShutdownManager::Handler TimerManager::d_finalizeHandler(0, 0, 0, TimerManager::FinalizeCallback, StartupShutdownManager::priorityTimerManger);

    TimerManager* TimerManager::GetInstance()
    {
        if (!d_instance)
        {
            d_instance = new TimerManager();
            d_instance->ComputeOverheadConstants();
            d_instance->d_mainTimer->Start();
        }
        return d_instance;
    }

    Timer* TimerManager::GetTimer(const std::string& name, bool ignore_timer_input)
    {
        Timer* timer;

        ARIES_ASSERT(!name.empty());
        bool timerActive = true;
        if (!ignore_timer_input)
        {
            timerActive = CheckTimerInNameLists(name);
        }

        /*
         * Add the timer to the appropriate array, if necessary.
         */
        if (timerActive)
        {
            if (!CheckTimerExistsInArray(timer, name, d_timers))
            {
                if (d_timers.size() == d_timers.capacity())
                {
                    d_timers.reserve(d_timers.size() + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
                }
                delete timer;
                timer = new Timer(name);
                
                d_timers.push_back(timer);
            }
        }
        else
        {
            if (!CheckTimerExistsInArray(timer, name, d_inactiveTimers))
            {
                if (d_inactiveTimers.size() == d_inactiveTimers.capacity())
                {
                    d_inactiveTimers.reserve(d_inactiveTimers.size() + DEFAULT_NUMBER_OF_TIMERS_INCREMENT);
                }
                
                delete timer;
                timer = new Timer(name);
                
                timer->SetActive(false);
                d_inactiveTimers.push_back(timer);
            }
        }
        return timer;
    }

    bool TimerManager::CheckTimerExists(Timer* timer, const std::string& name) const
    {
        ARIES_ASSERT(!name.empty());

        bool timer_found = CheckTimerExistsInArray(timer, name, d_timers);
        if (!timer_found)
            timer_found = CheckTimerExistsInArray(timer, name, d_inactiveTimers);
        
        return timer_found;
    }
    
    bool TimerManager::CheckTimerRunning(const std::string& name) const
    {
        bool is_running = false;
        ARIES_ASSERT(!name.empty());
        
        Timer* timer;
        if (CheckTimerExistsInArray(timer, name, d_timers))
            is_running = timer->IsRunning();
        
        return is_running;
    }
    
    void TimerManager::ResetAllTimers()
    {
        d_mainTimer->Stop();
        d_mainTimer->Reset();
        d_mainTimer->Start();

        for (size_t i = 0; i < d_timers.size(); ++i)
        {
            d_timers[i]->Reset();
        }
        for (size_t j = 0; j < d_inactiveTimers.size(); ++j)
        {
            d_inactiveTimers[j]->Reset();
        }
    }
    
    void TimerManager::Print(std::ostream& os)
    {
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        /*
         * There are 18 possible timer values that users may wish to look at.
         * (i.e. User/sys/wallclock time, Total or Exclusive, for individual
         * processor, max across all procs, or summed across all procs).
         * This method builds an array that holds these values and outputs
         * them in column format.
         */

        /*
         * First, stop the main_timer so we have an accurate measure of
         * overall run time so far.
         */
        d_mainTimer->Stop();

        /*
         * If we are doing max or sum operations, make sure timers are
         * consistent across processors.
         */
        if (d_printSummed || d_printMax)
        {
            CheckConsistencyAcrossProcessors();
        }

        /*
         * Invoke arrays used to hold timer names, timer_values, and
         * max_processor_ids (i.e. processor holding the maximum time).
         */
        double(*timer_values)[18] = new double[d_timers.size() + 1][18];
        int(*max_processor_id)[2] = new int[d_timers.size() + 1][2];
        std::vector<std::string> timer_names(static_cast<int>(d_timers.size()) + 1);

        /*
         * Fill in timer_values and timer_names arrays, based on values of
         * d_printTotal, d_printExclusive,
         * d_printUser, d_printSys, d_printWallclock,
         * d_printProcessor, d_printSummed, d_printMax.
         */
        BuildTimerArrays(timer_values, max_processor_id, timer_names);

        /*
         * Now that we have built up the array, select output options.
         *
         * 1) If user has requested any two (or more) of {user, system, and
         *    walltime} then use the format:
         *    [Exclusive,Total]
         *    [Processor,Summed,Max]
         *       Name     User    System   Wall  (Max Processor)
         *
         * 2) If user chose just one of {User, system, or walltime}, then use the
         *    format:
         *    [Exclusive,Total]
         *       Name     Processor   Summed   Max  Max Processor
         *
         * 3) If user chose just one of {user, system, or walltime} and just one
         *    of {processor, summed, max} then use the format:
         *       Name     [Processor,Summed,Max]  Total   Exclusive
         *
         * If the user wants overhead stats, print those as:
         *   Timer Overhead:
         *     Timer Name    number calls     Estimated overhead
         *
         * If they want output of a concurrent timer tree, print this as:
         *   Concurrent Tree:
         *     Timer Name    names of timers called by it.
         */
        
        /*
         * Determine which case we are doing - #1, #2, or #3
         * #1 case - user requested any two (or more) of [user,system,wallclock]
         * #2 case - user requested one of [user,system,wallclock]
         * #3 case - user requested one of [user,system,wallclock] and one of [processor, summed, max]
         */
        bool case1 = false;
        bool case2 = false;
        bool case3 = false;
        if ((d_printUser && d_printSys)  ||
            (d_printSys  && d_printWall) ||
            (d_printWall && d_printUser))
        {
            case1 = true;
        }
        else
        {
            case2 = true;
            if (( d_printProcessor && !d_printSummed && !d_printMax) ||
                (!d_printProcessor &&  d_printSummed && !d_printMax) ||
                (!d_printProcessor && !d_printSummed &&  d_printMax))
            {
                case2 = false;
                case3 = true;
            }
        }

        std::string table_title;
        std::vector<std::string> column_titles(4);
        int column_ids[3] = { 0, 0, 0 };
        int j, k;

        /*
         * Now print out the data
         */
        if (case1)
        {
            column_titles[0] = "";
            column_titles[1] = "";
            column_titles[2] = "";
            column_titles[3] = "Proc";

            if (d_printUser)
                column_titles[0] = "User Time";
            else 
                column_titles[0] = "";
      
            if (d_printSys) 
                column_titles[1] = "Sys Time";
            else 
                column_titles[1] = "";
      
            if (d_printWall) 
                column_titles[2] = "Wall Time";
            else 
                column_titles[2] = "";
      

            for (k = 0; k < 2; ++k)
            {
                if ((k == 0 && d_printExclusive) ||
                    (k == 1 && d_printTotal))
                {
                    for (j = 0; j < 3; ++j)
                    {
                        if ((j == 0 && d_printProcessor) ||
                            (j == 1 && d_printSummed) ||
                            (j == 2 && d_printMax))
                        {
                            if (j == 0)
                            {
                                std::ostringstream out;
                                if (k == 0)
                                {
                                    out << "EXCLUSIVE TIME \nPROCESSOR:"
                                        << mpi.GetRank();
                                    table_title = out.str();

                                    column_ids[0] = 0;
                                    column_ids[1] = 1;
                                    column_ids[2] = 2;
                                }
                                else if (k == 1)
                                {
                                    out << "TOTAL TIME \nPROCESSOR:"
                                        << mpi.GetRank();
                                    table_title = out.str();

                                    column_ids[0] = 9;
                                    column_ids[1] = 10;
                                    column_ids[2] = 11;
                                }
                                PrintTable(table_title,
                                           column_titles,
                                           timer_names,
                                           column_ids,
                                           timer_values,
                                           os);
                            }
                            else if (j == 1)
                            {
                                if (k == 0)
                                {
                                    table_title = "EXCLUSIVE TIME \nSUMMED ACROSS ALL PROCESSORS";
                                    column_ids[0] = 3;
                                    column_ids[1] = 4;
                                    column_ids[2] = 5;
                                }
                                else if (k == 1)
                                {
                                    table_title = "TOTAL TIME \nSUMMED ACROSS ALL PROCESSORS:";
                                    column_ids[0] = 12;
                                    column_ids[1] = 13;
                                    column_ids[2] = 14;
                                }
                                PrintTable(table_title,
                                           column_titles,
                                           timer_names,
                                           column_ids,
                                           timer_values,
                                           os);
                                
                            }
                            else if (j == 2)
                            {
                                int max_array_id = 0; // identifies which of the two
                                // max_processor_id values to print
                                if (k == 0)
                                {
                                    table_title = "EXCLUSIVE TIME \nMAX ACROSS ALL PROCESSORS";
                                    column_ids[0] = 6;
                                    column_ids[1] = 7;
                                    column_ids[2] = 8;
                                    max_array_id = 0;
                                }
                                else if (k == 1)
                                {
                                    table_title = "TOTAL TIME \nMAX ACROSS ALL PROCESSORS";
                                    column_ids[0] = 15;
                                    column_ids[1] = 16;
                                    column_ids[2] = 17;
                                    max_array_id = 1;
                                }
                                PrintTable(table_title,
                                           column_titles,
                                           timer_names,
                                           max_processor_id,
                                           max_array_id,
                                           column_ids,
                                           timer_values,
                                           os);
                            }
                        } // if j
                    } // for j
                } // if k
            } // for k
        } // if case 1

        if (case2)
        {
            for (k = 0; k < 2; ++k)
            {
                if ((k == 0 && d_printExclusive) ||
                    (k == 1 && d_printTotal))
                {
                    int max_array_id = 0;
                    std::string table_title_line_1;
                    std::string table_title_line_2;
                    if (k == 0)
                    {
                        table_title_line_1 = "EXCLUSVE \n";
                        max_array_id = 0;
                    }
                    else if (k == 1)
                    {
                        table_title_line_1 = "TOTAL \n";
                        max_array_id = 1;
                    }
                    if (d_printUser)
                    {
                        table_title_line_2 = "USER TIME";
                    }
                    else if (d_printSys)
                    {
                        table_title_line_2 = "SYSTEM TIME";
                    }
                    else if (d_printWall)
                    {
                        table_title_line_2 = "WALLCLOCK TIME";
                    }
                    table_title = table_title_line_1;
                    table_title += table_title_line_2;

                    column_titles[0] = "";
                    column_titles[1] = "";
                    column_titles[2] = "";
                    column_titles[3] = "";
                    if (d_printProcessor)
                    {
                        std::ostringstream out;
                        out << "Proc: " << mpi.GetRank();
                        column_titles[0] = out.str();
                    }
                    if (d_printSummed)
                    {
                        column_titles[1] = "Summed";
                    }
                    if (d_printMax)
                    {
                        column_titles[2] = "Max";
                        column_titles[3] = "Proc";
                    }

                    if (d_printUser)
                    {
                        if (k == 0)
                        {
                            column_ids[0] = 0;
                            column_ids[1] = 3;
                            column_ids[2] = 6;
                        }
                        else if (k == 1)
                        {
                            column_ids[0] = 9;
                            column_ids[1] = 12;
                            column_ids[2] = 15;
                        }
                    }
                    else if (d_printSys)
                    {
                        if (k == 0)
                        {
                            column_ids[0] = 1;
                            column_ids[1] = 4;
                            column_ids[2] = 7;
                        }
                        else if (k == 1)
                        {
                            column_ids[0] = 10;
                            column_ids[1] = 13;
                            column_ids[2] = 16;
                        }
                    }
                    else if (d_printWall)
                    {
                        if (k == 0)
                        {
                            column_ids[0] = 2;
                            column_ids[1] = 5;
                            column_ids[2] = 8;
                        }
                        else if (k == 1)
                        {
                            column_ids[0] = 11;
                            column_ids[1] = 14;
                            column_ids[2] = 17;
                        }
                    }

                    PrintTable(table_title,
                               column_titles,
                               timer_names,
                               max_processor_id,
                               max_array_id,
                               column_ids,
                               timer_values,
                               os);
                } // if k
            }  // for k
        } // if case2
        
        if (case3)
        {
            if (d_printExclusive && !d_printTotal)
            {
                column_titles[0] = "Exclusive";
                column_titles[1] = "";
            }
            else if (!d_printExclusive && d_printTotal)
            {
                column_titles[0] = "";
                column_titles[1] = "Total";
            }
            else if (d_printExclusive && d_printTotal)
            {
                column_titles[0] = "Exclusive";
                column_titles[1] = "Total";
            }
            column_titles[3] = "";

            column_ids[2] = 0;
            if (d_printUser)
            {
                if (d_printProcessor)
                {
                    std::ostringstream out;
                    out << "USER TIME \nPROCESSOR: " << mpi.GetRank();
                    table_title = out.str();

                    column_ids[0] = 0;
                    column_ids[1] = 9;
                }
                else if (d_printSummed)
                {
                    table_title = "USER TIME \nSUMMED ACROSS ALL PROCESSORS";
                    column_ids[0] = 3;
                    column_ids[1] = 12;
                }
                else if (d_printMax)
                {
                    table_title = "USER TIME \nMAX ACROSS ALL PROCESSORS";
                    column_ids[0] = 6;
                    column_ids[1] = 15;
                }
            }
            else if (d_printSys)
            {
                if (d_printProcessor)
                {
                    std::ostringstream out;
                    out << "SYSTEM TIME \nPROCESSOR: " << mpi.GetRank();
                    table_title = out.str();

                    column_ids[0] = 1;
                    column_ids[1] = 10;
                }
                else if (d_printSummed)
                {
                    table_title = "SYSTEM TIME \nSUMMED ACROSS ALL PROCESSORS";
                    column_ids[0] = 4;
                    column_ids[1] = 13;
                }
                else if (d_printMax)
                {
                    table_title = "SYSTEM TIME \nMAX ACROSS ALL PROCESSORS";
                    column_ids[0] = 7;
                    column_ids[1] = 16;
                }
            }
            else if (d_printWall)
            {
                if (d_printProcessor)
                {
                    std::ostringstream out;
                    out << "WALLCLOCK TIME \nPROCESSOR: " << mpi.GetRank();
                    table_title = out.str();

                    column_ids[0] = 2;
                    column_ids[1] = 11;
                }
                else if (d_printSummed)
                {
                    table_title = "WALLCLOCK TIME \nSUMMED ACROSS ALL PROCESSORS";
                    column_ids[0] = 5;
                    column_ids[1] = 14;
                }
                else if (d_printMax)
                {
                    table_title = "WALLCLOCK TIME \nMAX ACROSS ALL PROCESSORS";
                    column_ids[0] = 8;
                    column_ids[1] = 17;
                }
            }
            PrintTable(table_title,
                       column_titles,
                       timer_names,
                       column_ids,
                       timer_values,
                       os);
        }

        /*
         * Print overhead stats - number of accesses and estimated cost
         * (estimated cost computed based on the number of accesses and
         * a fixed d_timerAciveAccessTime value).
         * Store the number of accesses in max_processor_id[0] and the estimated
         * cost in timer_values[0] and use the PrintTable method.
         */
        if (d_printTimerOverhead)
        {
            PrintOverhead(timer_names,
                          timer_values,
                          os);
        }

        /*
         * Print tree of concurrent timers.
         */
        if (d_printConcurrent)
        {
            PrintConcurrent(os);
        }

        delete[] timer_values;
        delete[] max_processor_id;
        /*
         * Lastly, restart the main_timer that we stopped at the beginning of
         * this routine
         */
        d_mainTimer->Start();
    }

    TimerManager::TimerManager():
            d_timerActiveAccessTime(-9999.0),
            d_timerInactiveAccessTime(-9999.0),
            d_mainTimer(new Timer("TOTAL RUN TIME")),
            d_lengthPackageNames(0),
            d_lengthClassNames(0),
            d_lengthClassMethodNames(0),
            d_printThreshold(0.25),
            d_printExclusive(false),
            d_printTotal(true),
            d_printProcessor(true),
            d_printMax(false),
            d_printSummed(false),
            d_printUser(false),
            d_printSys(false),
            d_printWall(true),
            d_printPercentage(true),
            d_printConcurrent(false),
            d_printTimerOverhead(false)
    {
        GetFromInput();
    }

    TimerManager::~TimerManager()
    {
        d_mainTimer->Stop();
        d_mainTimer->Reset();

        d_timers.clear();
        d_inactiveTimers.clear();

        d_exclusiveTimerStack.clear();

        d_packageNames.clear();
        d_classNames.clear();
        d_classMethodNames.clear();
    }

    void TimerManager::RegisterSingletonSubclassInstance(TimerManager* subclassInstance)
    {
        if (!d_instance)
        {
            d_instance = subclassInstance;
        }
        else
        {
            ARIES_ERROR("TimerManager internal error...\n"
                        << "Attemptng to set Singleton instance to subclass instance,"
                        << "\n but Singleton instance already set." << std::endl);
        }
    }
    
    void TimerManager::StartTime(Timer* timer)
    {
        ARIES_ASSERT(timer != 0);

        if (timer->IsActive())
        {
        }

        if (d_printExclusive)
        {
            if (!d_exclusiveTimerStack.empty())
            {
                ((Timer*)d_exclusiveTimerStack.front())->StopExclusive();
            }
            Timer* stack_timer = timer;
            d_exclusiveTimerStack.push_front(stack_timer);
            stack_timer->StartExclusive();
        }

        if (d_printConcurrent)
        {
            for (size_t i = 0; i < d_timers.size(); ++i)
            {
                if ((d_timers[i] != timer) && d_timers[i]->IsRunning())
                {
                    d_timers[i]->AddConcurrentTimer(*d_timers[i]);
                }
            }
        }
    }

    void TimerManager::StopTime(Timer* timer)
    {
        ARIES_ASSERT(timer != 0);

        if (d_printExclusive)
        {
            timer->StopExclusive();
            if (!d_exclusiveTimerStack.empty())
            {
                d_exclusiveTimerStack.pop_front();
                if (!d_exclusiveTimerStack.empty())
                {
                    ((Timer *)d_exclusiveTimerStack.front())->StartExclusive();
                }
            }
        }
    }

    void TimerManager::ActivateExistingTimers()
    {
        std::vector<Timer* >::iterator it = d_inactiveTimers.begin();
        while (it != d_inactiveTimers.end())
        {
            bool timer_active = CheckTimerInNameLists((*it)->GetName());
            if (timer_active)
            {
                (*it)->SetActive(true);
                d_timers.push_back((*it));
                it = d_inactiveTimers.erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    bool TimerManager::CheckTimerExistsInArray(Timer* timer, const std::string& name, const std::vector<Timer* >& timer_array) const
    {
        bool timer_found = false;

        timer->Reset();
        if (!name.empty())
        {
            for (size_t i = 0; i < timer_array.size(); ++i)
            {
                if (timer_array[i]->GetName() == name)
                {
                    timer_found = true;
                    timer = timer_array[i];
                    break;
                }
            }
        }
        return timer_found;
    }

    void TimerManager::PrintTable(const std::string& table_title,
                                  const std::vector<std::string>& column_titles,
                                  const std::vector<std::string>& timer_names,
                                  const int column_ids[],
                                  const double timer_values[][18],
                                  std::ostream& os)
    {
        std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
        std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
        std::string ascii_line = ascii_line1;
        ascii_line += ascii_line2;
        
        /*
         * By default, output in C++ is right justified with the setw()
         * option e.g. cout << "[" << setw(5) << 1 << "]" will output
         * [   1].  The line below makes it left justified, so the same line
         * will generate [1   ].  We us left justification because it is
         * more convenient to output columns of tables.
         */
        os.setf(std::ios::left);

        os << ascii_line << table_title << "\n";

        int i;
        /*
         * Determine maximum name length for formatting
         */
        int maxlen = 10;
        for (unsigned int n = 0; n < d_timers.size() + 1; ++n)
        {
            i = static_cast<int>(timer_names[n].size());
            if (i > maxlen)
                maxlen = i;
        }
        
        /*
         * Print table header.  If we are only printing the overall program
         * timer (i.e. d_num_timers = 0) with only d_printProcessor,
         * d_printTotal, and d_printWall options being true (which
         * is the default case if the user doesn't add a "TimerManager"
         * section to the input file) then don't bother to print header as
         * it just clutters up the output.  Also, turn off percentages since
         * this doesn't mean anything with just one timer.
         */
        bool default_case =
                !d_printExclusive &&
                !d_printSummed &&
                !d_printMax &&
                !d_printUser &&
                !d_printSys &&
                (d_timers.size() == 0);
        
        if (default_case)
        {
            d_printPercentage = false;
        }
        else
        {
            os << ascii_line
               << std::setw(maxlen + 3) << "Timer Name" << ' ';
            for (i = 0; i < 3; ++i)
            {
                if (!column_titles[i].empty())
                {
                    os << std::setw(15) << column_titles[i].c_str() << "  ";
                }
            }
            os << std::endl;
        }

        /*
         * Organize timers largest to smallest.  Apply this to the LAST NONZERO
         * column entry for the table by forming an ordering array - ordered_list
         * - that orders these values.
         */
        int last_nonzero_column = 0;
        for (i = 0; i < 3; ++i)
        {
            if (!column_titles[i].empty())
            {
                last_nonzero_column = column_ids[i];
            }
        }
        int* ordered_list = new int[d_timers.size() + 1];
        BuildOrderedList(timer_values,
                         last_nonzero_column,
                         ordered_list,
                         static_cast<int>(d_timers.size()));

        /*
         * Tack on TOTAL TIME to end of ordered list
         */
        ordered_list[static_cast<int>(d_timers.size())] =  static_cast<int>(d_timers.size());

        /*
         * Now output the rows of the table.
         */
        for (size_t k = 0; k < d_timers.size() + 1; ++k)
        {
            int n = ordered_list[k];
            
            /*
             * Check the print threshold to see if we should print this timer.
             */
            double frac = ComputePercentageDouble(timer_values[n][last_nonzero_column],
                                                  timer_values[d_timers.size()][last_nonzero_column]);

            if (frac > d_printThreshold)
            {
                os << std::setw(maxlen + 3) << timer_names[n].c_str() << ' ';

                /*
                 * Print column values
                 */
                for (i = 0; i < 3; ++i)
                {
                    /*
                     * Print column values only title is non-null (i.e. not "")
                     */
                    if (!column_titles[i].empty())
                    {
                        /*
                         * Print percentages if requested.
                         */
                        int j = column_ids[i];

                        if (d_printPercentage)
                        {
                            int perc = ComputePercentageInt(timer_values[n][j],
                                                            timer_values[d_timers.size()][j]);

                            std::ostringstream out;
                            out << timer_values[n][j] << " (" << perc << "%)";
                            os << std::setw(15) << out.str().c_str() << "  ";

                        }
                        else
                        {
                            os << std::setw(15) << timer_values[n][j] << "  ";
                        }
                    } // if title is non-null
                } // loop over columns
                os << std::endl;
            } // if meets d_printThreshold condition
        } // loop over timers
        delete[] ordered_list;

        os << ascii_line << std::endl;
        os.setf(std::ios::right);
    }

    void TimerManager::PrintTable(const std::string& table_title,
                                  const std::vector<std::string>& column_titles,
                                  const std::vector<std::string>& timer_names,
                                  const int max_processor_id[][2],
                                  const int max_array_id,
                                  const int column_ids[],
                                  const double timer_values[][18],
                                  std::ostream& os)
    {
        std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
        std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
        std::string ascii_line = ascii_line1;
        ascii_line += ascii_line2;

        /*
         * Left-justify all output in this method.
         */
        os.setf(std::ios::left);

        os << ascii_line
           << table_title << "\n"
           << ascii_line;

        int i;

        /*
         * Determine maximum name length for formatting
         */
        int maxlen = 10;
        for (unsigned int n = 0; n < d_timers.size() + 1; ++n)
        {
            i = static_cast<int>(timer_names[n].size());
            if (i > maxlen) maxlen = i;
        }

        /*
         * Print table header
         */
        os << std::setw(maxlen + 3) << "Timer Name" << ' ';
        for (i = 0; i < 4; ++i)
        {
            if (!column_titles[i].empty())
            {
                os << std::setw(15) << column_titles[i].c_str() << "  ";
            }
        }
        os << std::endl;

        /*
         * Organize timers largest to smallest.  Apply this to the LAST NONZERO
         * column entry for the table by forming an ordering array - ordered_list
         * - that orders these values.
         */
        int last_nonzero_column = 0;
        for (i = 0; i < 3; ++i)
        {
            if (!column_titles[i].empty())
            {
                last_nonzero_column = column_ids[i];
            }
        }
        int* ordered_list = new int[d_timers.size() + 1];
        BuildOrderedList(timer_values, last_nonzero_column, ordered_list, static_cast<int>(d_timers.size()));

        /*
         * Tack on TOTAL TIME to end of ordered list
         */
        ordered_list[static_cast<int>(d_timers.size())] = static_cast<int>(d_timers.size());

        /*
         * Now output the rows of the table.
         */
        for (size_t j = 0; j < d_timers.size() + 1; ++j)
        {
            unsigned int n = ordered_list[j];

            /*
             * Check the print threshold to see if we should print this timer.
             */
            double frac = ComputePercentageDouble(timer_values[n][last_nonzero_column], timer_values[d_timers.size()][last_nonzero_column]);

            if (frac > d_printThreshold)
            {
                os << std::setw(maxlen + 3) << timer_names[n].c_str() << ' ';

                /*
                 * Print columns.
                 */
                for (i = 0; i < 4; ++i)
                {
                    /*
                     * Print column values only title is non-null (i.e. not "")
                     */
                    if (!column_titles[i].empty())
                    {
                        /*
                         * Print percentages for columns 0-2
                         */
                        if (i < 3)
                        {
                            int k = column_ids[i];

                            if (d_printPercentage)
                            {
                                int perc = ComputePercentageInt(timer_values[n][k],
                                                                timer_values[d_timers.size()][k]);
                                std::ostringstream out;
                                out << timer_values[n][k] << " (" << perc << "%)";
                                os << std::setw(15) << out.str().c_str() << "  ";
                            }
                            else
                            {
                                os << std::setw(15) << timer_values[n][k] << "  ";
                            }

                        }
                        else
                        {

                            /*
                             * Print column 3 - the processor holding processor ID
                             * with max times (don't do for TOTAL TIME - this is
                             * meaningless since all processors are synchronized
                             * before and after this call).
                             */
                            if (n < d_timers.size())
                            {
                                os << std::setw(15) << max_processor_id[n][max_array_id];
                            }
                        } // column 3
                    }  // if column title is non-null
                } // loop over columns
                os << std::endl;
            } //  matches d_printThreshold conditions
        } // loop over timers

        delete[] ordered_list;

        os << ascii_line << std::endl;
        os.setf(std::ios::right);
    }
    

    void TimerManager::PrintOverhead(const std::vector<std::string>& timer_names,
                                     const double timer_values[][18],
                                     std::ostream& os)
    {
        std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
        std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
        std::string ascii_line = ascii_line1;
        ascii_line += ascii_line2;

        /*
         * Left-justify all output in this method.
         */
        os.setf(std::ios::left);

        os << ascii_line
           << "TIMER OVERHEAD STATISTICS \n"
           << ascii_line;

        /*
         * Determine maximum name length for formatting
         */
        int maxlen = 10;
        for (unsigned int n = 0; n < d_timers.size(); ++n)
        {
            int i = static_cast<int>(timer_names[n].size());
            if (i > maxlen) maxlen = i;
        }

        /*
         * Print table header
         */
        os << std::setw(maxlen + 3) << "Timer Name"
           << std::setw(25) << "Number Accesses" << "  "
           << std::setw(25) << "Estimated Cost"
           << std::endl;

        /*
         * Compute totals: total number of REGISTERED accesses and total cost.
         * Total cost includes inactive timer costs.
         */
        int total_inactive_accesses = 0;
        for (size_t i = 0; i < d_inactiveTimers.size(); ++i)
        {
            total_inactive_accesses += d_inactiveTimers[i]->GetNumberAccesses();
        }

        double est_cost = d_timerInactiveAccessTime * total_inactive_accesses;
        double total_est_cost = est_cost;

        int total_accesses = 0;
        for (size_t n = 0; n < d_timers.size(); ++n)
        {
            total_accesses += d_timers[n]->GetNumberAccesses();
        }
        est_cost = d_timerActiveAccessTime * total_accesses;

        /*
         * If we are keeping exclusive or concurrent times, each access costs
         * roughly four times as much.  Make this correction here...
         */
        if (d_printExclusive || d_printConcurrent)
        {
            est_cost *= 4.;
        }
        total_est_cost += est_cost;
        
        /*
         * Output the rows of the table.  Start first with the inactive timers...
         */
        int num_accesses = total_inactive_accesses;
        est_cost = d_timerInactiveAccessTime * num_accesses;
        int perc = ComputePercentageInt(est_cost, total_est_cost);

        os << std::setw(maxlen + 3) << "inactive timers"
           << std::setw(25) << num_accesses << "  ";
        std::ostringstream out;
        out << est_cost << " (" << perc << "%)";
        os << std::setw(25) << out.str().c_str();
        os << std::endl;

        /*
         * Now print the rest of the timers.  While we are cycling through them,
         * add up the total cost and print it at the end...
         */
        for (unsigned int n = 0; n < d_timers.size(); ++n)
        {
            num_accesses = d_timers[n]->GetNumberAccesses();
            est_cost = d_timerActiveAccessTime * num_accesses;

            /*
             * If we are keeping exclusive or concurrent times, each access costs
             * roughly four times as much.  Make this correction here...
             */
            if (d_printExclusive || d_printConcurrent)
            {
                est_cost *= 4.;
            }

            perc = ComputePercentageInt(est_cost, total_est_cost);

            os << std::setw(maxlen + 3) << timer_names[n].c_str()
               << std::setw(25) << num_accesses << "  ";
            std::ostringstream out2;
            out2 << est_cost << " (" << perc << "%)";
            os << std::setw(25) << out2.str().c_str();
            
            os << std::endl;
        }

        /*
         * Output the totals.
         */
        os << std::setw(maxlen + 3) << "TOTAL:"
           << std::setw(25) << total_accesses << "  "
           << std::setw(25) << total_est_cost
           << "\n" << std::endl;

        /*
         * Compare the total estimated cost with overall program wallclock time.
         * If it is a significant percentage (> 5%) print a warning
         */
        double perc_dbl = ComputePercentageDouble(total_est_cost,
                                                  timer_values[d_timers.size()][11]);

        os << "Estimated Timer Costs as a percentage of overall Wallclock Time: "
           << perc_dbl << "% \n";
        if (perc_dbl > 5.)
        {
            os << "WARNING:  TIMERS ARE USING A SIGNIFICANT FRACTION OF RUN TIME"
               << std::endl;
        }

        os << ascii_line << std::endl;
        os.setf(std::ios::right);
    }
    
    void TimerManager::PrintConcurrent(std::ostream& os)
    {
        std::string ascii_line1 = "++++++++++++++++++++++++++++++++++++++++";
        std::string ascii_line2 = "++++++++++++++++++++++++++++++++++++++++\n";
        std::string ascii_line = ascii_line1;
        ascii_line += ascii_line2;

        os << ascii_line
           << "CONCURRENT TIMERS\n"
           << ascii_line;

        /*
         * Determine maximum name length for formatting
         */
        int maxlen = 10;
        for (size_t n = 0; n < d_timers.size(); ++n)
        {
            int i = int((d_timers[n]->GetName()).size());
            if (i > maxlen) maxlen = i;
        }

        /*
         * Print table header
         */
        os << std::setw(maxlen + 3) << "Timer Name"
           << std::setw(25) << "Nested Timers"
           << std::endl;

        /*
         * Output the rows of the table.
         */
        for (size_t n = 0; n < d_timers.size(); ++n)
        {
            os << std::setw(maxlen + 3) << d_timers[n]->GetName().c_str();

            int count = 0;
            for (size_t i = 0; i < d_timers.size(); ++i)
            {
                if (d_timers[n]->IsConcurrentTimer(*d_timers[i]))
                {
                    ++count;
                }
            }
            if (count == 0)
            {
                os << std::setw(25) << "none " << std::endl;
            }
            else
            {
                /*
                 * Format it like:    Timer Name      Concurrent Timer #1
                 *                                    Concurrent Timer #2
                 *                                    ...
                 * Use "count" variable defined above to identify the first
                 * line or subsequent lines.
                 */
                count = 0;
                for (size_t j = 0; j < d_timers.size(); ++j)
                {
                    if (d_timers[n]->IsConcurrentTimer(*d_timers[j]))
                    {
                        if (count == 0)
                        {
                            os << std::setw(25) << d_timers[j]->GetName().c_str()
                               << std::endl;
                        }
                        else
                        {
                            os << std::setw(maxlen + 3) << " "
                               << d_timers[j]->GetName().c_str() << std::endl;
                        }
                        ++count;
                    }
                }
            }

        }
        os << ascii_line << std::endl;
    }

    void TimerManager::BuildTimerArrays(double timer_values[][18], int max_processor_id[][2], std::vector<std::string>& timer_names)
    {
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        /*
         * timer_values - 2D array dimensioned [d_timers.size()][18]
         *     For each timer, there are 18 potential values which may be of
         *     interest.  This array stores them if they are requested.
         * max_processor_id - 2D array dimensioned [d_timers.size()][2]
         *     Holds the value of the processor that used the maximum amount
         *     of time.  [0] is for exclusive time, while [1] is for total time.
         */

        /*
         * Initialize arrays
         */
        for (unsigned int n = 0; n < d_timers.size() + 1; ++n)
        {
            timer_names[n] = "";
            max_processor_id[n][0] = 0;
            max_processor_id[n][1] = 0;
            for (int i = 0; i < 18; ++i)
            {
                timer_values[n][i] = 0.;
            }
        }

        /*
         * Build arrays.
         */
        for (unsigned int n = 0; n < d_timers.size(); ++n)
        {
            timer_names[n] = d_timers[n]->GetName();

            /*
             *  Build timer_values[n][m] array:
             *    m = 0 :  processor exclusive user time
             *    m = 1 :  processor exclusive sys time
             *    m = 2 :  processor exclusive wall time
             *    m = 3 :  summed exclusive user time
             *    m = 4 :  summed exclusive sys time
             *    m = 5 :  summed exclusive wall time
             *    m = 6 :  max exclusive user time
             *    m = 7 :  max exclusive sys time
             *    m = 8 :  max exclusive wall time
             *    m = 9 :  processor total user time
             *    m = 10 :  processor total sys time
             *    m = 11 :  processor total wall time
             *    m = 12 :  summed total user time
             *    m = 13 :  summed total sys time
             *    m = 14 :  summed total wall time
             *    m = 15 :  max total user time
             *    m = 16 :  max total sys time
             *    m = 17 :  max total wall time
             */
            for (int k = 0; k < 2; ++k)
            {
                for (int j = 0; j < 3; ++j)
                {
                    if ((k == 0 && d_printExclusive) || (k == 1 && d_printTotal))
                    {
                        if ((j == 0 && d_printProcessor) ||
                            (j == 1 && d_printSummed) ||
                            (j == 2 && d_printMax))
                        {
                            if (k == 0 && j == 0)
                            {
                                if (d_printUser)
                                    timer_values[n][0] = d_timers[n]->GetExclusiveUserTime();
                                
                                if (d_printSys)
                                    timer_values[n][1] = d_timers[n]->GetExclusiveSystemTime();
                                
                                if (d_printWall)
                                    timer_values[n][2] = d_timers[n]->GetExclusiveWallclockTime();
                                
                            }
                            else if (k == 0 && j == 1)
                            {
                                if (d_printUser)
                                {
                                    timer_values[n][3] = d_timers[n]->GetExclusiveUserTime();
                                    if (mpi.GetSize() > 1)
                                        mpi.AllReduce(&timer_values[n][3], 1, MPI_SUM); 
                                }
                                if (d_printSys)
                                {
                                    timer_values[n][3] = d_timers[n]->GetExclusiveSystemTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][4], 1, MPI_SUM);
                                    }
                                }
                                if (d_printWall)
                                {
                                    timer_values[n][3] = d_timers[n]->GetExclusiveWallclockTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][5], 1, MPI_SUM);
                                    }
                                }
                            }
                            else if (k == 0 && j == 2)
                            {
                                if (d_printUser)
                                {
                                    double user_time = d_timers[n]->GetExclusiveUserTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.Allreduce(&user_time, &timer_values[n][6], 1, MPI_DOUBLE, MPI_MAX);
                                    }
                                }
                                if (d_printSys)
                                {
                                    double sys_time = d_timers[n]->GetExclusiveSystemTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.Allreduce(&sys_time, &timer_values[n][7], 1, MPI_DOUBLE, MPI_MAX);
                                    }
                                }
                                if (d_printWall)
                                {
                                    timer_values[n][8] = d_timers[n]->GetExclusiveWallclockTime();
                                    max_processor_id[n][0] = mpi.GetRank();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][8], 1, MPI_MAXLOC, &max_processor_id[n][0]);
                                    }
                                }

                            }
                            else if (k == 1 && j == 0)
                            {
                                if (d_printUser)
                                {
                                    timer_values[n][9] = d_timers[n]->GetTotalUserTime();
                                }
                                if (d_printSys)
                                {
                                    timer_values[n][10] = d_timers[n]->GetTotalSystemTime();
                                }
                                if (d_printWall)
                                {
                                    timer_values[n][11] = d_timers[n]->GetTotalWallclockTime();
                                }
                            }
                            else if (k == 1 && j == 1)
                            {
                                if (d_printUser)
                                {
                                    timer_values[n][12] = d_timers[n]->GetTotalUserTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][12], 1, MPI_SUM);
                                    }
                                }
                                if (d_printSys)
                                {
                                    timer_values[n][13] = d_timers[n]->GetTotalSystemTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][13], 1, MPI_SUM);
                                    }
                                }
                                if (d_printWall)
                                {
                                    timer_values[n][14] = d_timers[n]->GetTotalWallclockTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][14], 1, MPI_SUM);
                                    }
                                }
                            }
                            else if (k == 1 && j == 2)
                            {
                                if (d_printUser)
                                {
                                    double user_time = d_timers[n]->GetTotalUserTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.Allreduce(&user_time, &timer_values[n][15], 1, MPI_DOUBLE, MPI_MAX);
                                    }
                                }
                                if (d_printSys)
                                {
                                    double sys_time = d_timers[n]->GetTotalSystemTime();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.Allreduce(&sys_time, &timer_values[n][16], 1, MPI_DOUBLE, MPI_MAX);
                                    }
                                }
                                if (d_printWall)
                                {
                                    timer_values[n][17] = d_timers[n]->GetTotalWallclockTime();
                                    max_processor_id[n][1] = mpi.GetRank();
                                    if (mpi.GetSize() > 1)
                                    {
                                        mpi.AllReduce(&timer_values[n][17], 1, MPI_MAXLOC, &max_processor_id[n][1]);
                                    }
                                }
                            }
                        }   // if j
                    }   // if k
                } // loop over j
            } // loop over k
        } // loop over n
        
        /*
         * Store main_timer data in timer_values[d_timers.size()][] location.  Max
         * time and exclusive time are not determined since these don't really
         * mean anything for an overall measurement of run time.
         */
        timer_names[static_cast<int>(d_timers.size())] = "TOTAL RUN TIME:";
        if (d_printUser)
        {
            double main_time = d_mainTimer->GetTotalUserTime();
            timer_values[d_timers.size()][0] = main_time;
            timer_values[d_timers.size()][3] = main_time;
            timer_values[d_timers.size()][6] = main_time;
            timer_values[d_timers.size()][9] = main_time;
            timer_values[d_timers.size()][12] = main_time;
            timer_values[d_timers.size()][15] = main_time;
            if (mpi.GetSize() > 1)
            {
                mpi.Allreduce(&main_time, &timer_values[d_timers.size()][3], 1, MPI_DOUBLE, MPI_SUM);
                mpi.Allreduce(&main_time, &timer_values[d_timers.size()][12], 1, MPI_DOUBLE, MPI_SUM);
            }
        }
        if (d_printSys)
        {
            double main_time = d_mainTimer->GetTotalSystemTime();
            timer_values[d_timers.size()][1] = main_time;
            timer_values[d_timers.size()][4] = main_time;
            timer_values[d_timers.size()][7] = main_time;
            timer_values[d_timers.size()][10] = main_time;
            timer_values[d_timers.size()][13] = main_time;
            timer_values[d_timers.size()][16] = main_time;
            if (mpi.GetSize() > 1)
            {
                mpi.Allreduce(&main_time, &timer_values[d_timers.size()][4], 1, MPI_DOUBLE, MPI_SUM);
                mpi.Allreduce(&main_time, &timer_values[d_timers.size()][13], 1, MPI_DOUBLE, MPI_SUM);
            }
        }
        if (d_printWall)
        {
            double main_time = d_mainTimer->GetTotalWallclockTime();
            timer_values[d_timers.size()][2] = main_time;
            timer_values[d_timers.size()][5] = main_time;
            timer_values[d_timers.size()][8] = main_time;
            timer_values[d_timers.size()][11] = main_time;
            timer_values[d_timers.size()][14] = main_time;
            timer_values[d_timers.size()][17] = main_time;
            if (mpi.GetSize() > 1)
            {
                mpi.Allreduce(&main_time, &timer_values[d_timers.size()][5], 1, MPI_DOUBLE, MPI_SUM);
                mpi.Allreduce(&main_time, &timer_values[d_timers.size()][14], 1, MPI_DOUBLE, MPI_SUM);
            }
        }
    }

    void TimerManager::BuildOrderedList(const double timer_values[][18], const int column, int index[], const int array_size)
    {
        /*
         * initialize the arrays
         */
        std::vector<double> timer_vals(array_size);
        for (int i = 0; i < array_size; ++i)
        {
            index[i] = i;
            timer_vals[i] = timer_values[i][column];
        }

        /*
         * Do a quicksort on timer_values array to build index array
         * ordered_list.
         */
        Quicksort(timer_vals, index, 0, array_size - 1);
    }
    
    bool TimerManager::CheckTimerInNameLists(const std::string& copy)
    {
        std::string name = copy;
        /*
         * string::size_type is generally an int, but it may depend on vendor's
         * implementation of string class.  Just to be safe, we use the definition
         * from the string class (string::size_type).
         */
        std::string::size_type string_length, list_entry_length, position;
   
        /*
         * Specification of whether we will use the timer after comparing to
         * package, class, and class::method lists.
         */
        bool will_use_timer = false;

        /*
         * Start by evaluating the Timer's package.  If the list of packages
         * has length zero, we can skip this step. We can also skip this step
         * if the timer does not have two "::" in its name.
         *
         * The manager will assume the following with regard to the specified timer
         * name:  The name has...
         *  1) Two "::" - it is of the form Package::Class::method.  A user can
         *     turn on the timer by specifying the timers Package, Class, or
         *     Class::method combination in the input file.  This is the most
         *     versatile.
         *  2) One "::" - it has the form Class::method.  The timer is assumed not
         *     to be in a package.  A user can turn on the timer by specifying its
         *     Class or class::method combination.
         *  3) No "::" - it has the form Class.  The timer is can be turned on by
         *     entering the class name in the input file.
         */
        if (d_lengthPackageNames > 0)
        {
            /*
             * See how many instances of "::" there are in the name.  If there
             * are at least two, the timer has a package which might be in the
             * package list.
             */
            int occurrences = 0;
            position = name.find("::");
            if (position < name.size())
            {
                ++occurrences;
                std::string substring = name.substr(position + 2);
                position = substring.find("::");
                if (position < substring.size())
                {
                    ++occurrences;
                }
            }

            if (occurrences >= 2)
            {
                /*
                 * The timer may be in the package list.  Parse to get its package
                 * name and compare to the list entries.
                 */
                position = name.find("::");
                std::string package = name.substr(0, position);

                /*
                 * Iterate through package list and see if the timer's package
                 * name is there.
                 */
                bool package_exists = false;
                string_length = package.size();
                for (std::list<std::string>::iterator i = d_packageNames.begin();
                     i != d_packageNames.end(); ++i)
                {
                    list_entry_length = i->size();
                    if (string_length == list_entry_length)
                    {
                        package_exists = (*i == package);
                    }
                    if (package_exists)
                    {
                        break;
                    }
                }
                will_use_timer = package_exists;
            }
        }

        if (!will_use_timer)
        {
            /*
             * The timer's package has already been compared to the package list,
             * so we can parse it off to ease further evaluations of the name.
             */
            int occurrences = 0;
            position = name.find("::");
            if (position < name.size())
            {
                ++occurrences;
                std::string substring = name.substr(position + 2);
                position = substring.find("::");
                if (position < substring.size())
                {
                    ++occurrences;
                }
            }
            if (occurrences >= 2)
            {
                position = name.find("::");
                name = name.substr(position + 2);
            }

            /*
             * See if Timer's class is in d_classNames.  If the list of classes
             * has length zero, we can skip this step.
             */
            if (d_lengthClassNames > 0)
            {
                /*
                 * If a "::" exists in the timer name, parse off the part before
                 * before it.  This is the class name.  If no "::" exists, assume
                 * the timer name is a class name (don't do any parsing) and compare
                 * it directly to the class list.
                 */
                position = name.find("::");
                std::string class_name;
                if (position < name.size())
                {
                    class_name = name.substr(0, position);
                }
                else
                {
                    class_name = name;
                }

                /*
                 * Is class name dimensional?
                 */
                string_length = class_name.size();
                std::string dim = class_name.substr(string_length - 1, 1);
                bool is_dimensional = false;
                std::string nondim_class_name;
                if (dim == "1" || dim == "2" || dim == "3")
                {
                    is_dimensional = true;
                    nondim_class_name = class_name.substr(0, string_length - 1);
                }
                else
                {
                    nondim_class_name = class_name;
                }

                /*
                 * See if class name is in class list.  Accelerate this by comparing
                 * the length of the entries.  Do non-dimensional comparison first.
                 */
                string_length = nondim_class_name.size();
                bool class_exists = false;
                for (std::list<std::string>::iterator i = d_classNames.begin();
                     i != d_classNames.end(); ++i)
                {
                    list_entry_length = i->size();
                    if (string_length == list_entry_length)
                    {
                        class_exists = (*i == nondim_class_name);
                    }
                    if (class_exists)
                    {
                        break;
                    }
                }

                /*
                 * Now do dimensional comparison.
                 */
                string_length = class_name.size();
                if (is_dimensional && !class_exists)
                {
                    for (std::list<std::string>::iterator i = d_classNames.begin(); i != d_classNames.end(); ++i)
                    {
                        list_entry_length = i->size();
                        if (string_length == list_entry_length)
                        {
                            class_exists = (*i == class_name);
                        }
                        if (class_exists)
                        {
                            break;
                        }
                    }
                }
                will_use_timer = class_exists;
            }

            /*
             * See if Timer's class::method name is in d_classMethodNames list.
             *
             * If the list of class_method_names has length zero, we can skip
             * this step.  Also, if no "::" exists in the timer name then it
             * cannot be in the timer list so lets avoid the comparison.
             */
            position = name.find("::");
            occurrences = 0;
            if (position < name.size())
            {
                occurrences = 1;
            }

            if (!will_use_timer && d_lengthClassMethodNames > 0 && occurrences > 0)
            {
                /*
                 * Parse name before "::" - this is the class name.
                 */
                position = name.find("::");
                std::string class_name = name.substr(0, position);

                /*
                 * Is class name dimensional?
                 */
                string_length = class_name.size();
                std::string dim = class_name.substr(string_length - 1, 1);
                bool is_dimensional = false;
                std::string nondim_name;
                if (dim == "1" || dim == "2" || dim == "3")
                {
                    is_dimensional = true;
                    std::string nondim_class_name = class_name.substr(0,
                                                                      string_length - 1);
                    std::string method_name = name.substr(position);
                    nondim_name = nondim_class_name;
                    nondim_name += method_name;

                }
                else
                {
                    nondim_name = name;
                }

                /*
                 * See if name is in class_method_names list.  Accelerate this by
                 * comparing the length of the entries.  Do non-dimensional
                 * comparison first.
                 */
                bool class_method_exists = false;
                string_length = nondim_name.size();
                for (std::list<std::string>::iterator i = d_classMethodNames.begin(); i != d_classMethodNames.end(); ++i)
                {
                    list_entry_length = i->size();
                    if (string_length == list_entry_length)
                    {
                        class_method_exists = (*i == nondim_name);
                    }
                    if (class_method_exists)
                    {
                        break;
                    }
                }

                /*
                 * Now do dimensional comparison.
                 */
                if (is_dimensional && !class_method_exists)
                {
                    string_length = name.size();
                    for (std::list<std::string>::iterator i = d_classMethodNames.begin(); i != d_classMethodNames.end(); ++i)
                    {
                        list_entry_length = i->size();
                        if (string_length == list_entry_length)
                        {
                            class_method_exists = (*i == name);
                        }
                        if (class_method_exists)
                        {
                            break;
                        }
                    }
                }

                will_use_timer = class_method_exists;
            }
        }
        return will_use_timer;
    }

    void TimerManager::CheckConsistencyAcrossProcessors()
    {
        /*
         * Due to the difficulty of comparing strings using MPI calls,
         * we do a rough consistency check of
         * 1. the number of timers and
         * 2. the length of each timer name.
         *
         * Steps:
         * 1. Do global reductions to get the max number of timers
         *    and the max lengths of each timer name.
         * 2. Issue a warning if the number of timers is inconsistent.
         *    This inconsistency would be found on all processes
         *    except those with the biggest number of timers.
         * 3. Issue a warning for each individual timer if
         *    its name length is less than the max length of
         *    all timers at the same index in the timer manager.
         *    Even if the number of timers are consistent, this
         *    would find wrong timer orderings or inconsistent
         *    timer names, unless the errors are for timer names
         *    with identical lengths.
         * 4. Go global reductions to get the number of inconsistencies
         *    of other processes.  Turn off printing of sum and max
         *    if any processes has inconsistencies.
         *
         * In the future, we may want to convert the strings into
         * their MD5 signatures and compare those as integers.
         */

        const AriesMPI& mpi(AriesMPI::GetAriesWorld());

        unsigned int max_num_timers = static_cast<unsigned int>(d_timers.size());
        if (mpi.GetSize() > 1)
        {
            int i = static_cast<int>(d_timers.size());
            mpi.Allreduce(&i, &max_num_timers, 1, MPI_INT, MPI_MAX);
        }

        std::vector<int> max_timer_lengths(max_num_timers);
        std::vector<int> rank_of_max(max_num_timers, mpi.GetRank());

        for (unsigned int i = 0; i < max_num_timers; ++i)
        {
            max_timer_lengths[i] = i < static_cast<unsigned int>(d_timers.size()) ? static_cast<int>(d_timers[i]->GetName().size()) : 0;
        }

        if (mpi.GetSize() > 1)
        {
            mpi.AllReduce(&max_timer_lengths[0], max_num_timers, MPI_MAXLOC, &rank_of_max[0]);
        }

        int inconsistency_count = 0;

        if (max_num_timers > d_timers.size())
        {
            ARIES_WARNING("Timer selections across processors were determined to be"
                          << "\ninconsistent.  This processor has only "
                          << d_timers.size() << " while some has " << max_num_timers
                          << ".\nThe consistency check"
                          << "\nwill continue for this process, but checking only\n"
                          << d_timers.size() << " timers."
                          << "\nIt is not possible to print global"
                          << "\nsummed or max timer information." << std::endl);
            ++inconsistency_count;
        }

        for (unsigned int i = 0; i < d_timers.size(); ++i)
        {
            if (max_timer_lengths[i] != int(d_timers[i]->GetName().size()))
            {
                ARIES_WARNING("Timer[" << i << "]: " << d_timers[i]->GetName()
                              << "\nis not consistent across all processors."
                              << "\nOther timer[" << i << "] has up to "
                              << max_timer_lengths[i] << " characters in their names."
                              << "\nIt is not possible to print global"
                              << "\nsummed or max timer information."
                              << std::endl);
                ++inconsistency_count;
            }
        }

        int max_inconsistency_count = inconsistency_count;
        if (mpi.GetSize() > 1)
        {
            mpi.Allreduce(&inconsistency_count, &max_inconsistency_count, 1,  MPI_INT, MPI_MAX);
        }
        if (max_inconsistency_count > 0)
        {
            d_printSummed = false;
            d_printMax = false;
            if (inconsistency_count == 0)
            {
                ARIES_WARNING("Though this process found no timer inconsistencies,"
                              << "\nother processes did.  It is not possible to print"
                              << "\nglobal summed or max timer information." << std::endl);
            }
        }

        /*
         * NOTE:  It might be nice to someday add the capability to remove the
         * inconsistent timers and print the max/summed values of the
         * consistent ones.   Unfortunately, this is tough to implement.  If it
         * were just a matter of comparing timer names across processors it would be
         * easy. But with MPI, only ints and doubles can be exchanged across
         * processors so it is difficult to make string comparisons.
         * It is possible to compare the MD5 sum of the strings,
         * but that may make SAMRAI dependent on the MD5 library.
         */
    }
    
    void TimerManager::GetFromInput()
    {
        //if ()
        {
            //d_printExclusive = input_db->getBoolWithDefault("print_exclusive", false);
            //d_printTotal = input_db->getBoolWithDefault("print_total", true);
            //d_printProcessor = input_db->getBoolWithDefault("print_processor", true);
            //d_printMax = input_db->getBoolWithDefault("print_max", false);
            //d_printSummed = input_db->getBoolWithDefault("print_summed", false);
            //d_printUser = input_db->getBoolWithDefault("print_user", false);
            //d_printSys = input_db->getBoolWithDefault("print_sys", false);
            //d_printWall = input_db->getBoolWithDefault("print_wall", true);
            //d_printPercentage = input_db->getBoolWithDefault("print_percentage", true);
            //d_printConcurrent = input_db->getBoolWithDefault("print_concurrent", false);
            //d_printTimerOverhead = input_db->getBoolWithDefault("print_timer_overhead", false);
            //d_printThreshold = input_db->getDoubleWithDefault("print_threshold", 0.25);

            std::vector<std::string> timer_list;
            //if (input_db->keyExists("timer_list"))
            //{
            //    timer_list = input_db->getStringVector("timer_list");
            //}

            /*
             *  Step thru the input list and call addTimerToNameLists to add
             *  the input file entry to the d_packageNames,
             *  d_classNames, and d_classMethodNames lists.
             */
            for (int i = 0; i < static_cast<int>(timer_list.size()); ++i)
            {
                std::string entry = timer_list[i];
                AddTimerToNameLists(entry);
            }
            d_lengthPackageNames = static_cast<int>(d_packageNames.size());
            d_lengthClassNames = static_cast<int>(d_classNames.size());
            d_lengthClassMethodNames = static_cast<int>(d_classMethodNames.size());
        }
    }

    void TimerManager::AddTimerToNameLists(const std::string& name)
    {
        /*
         * Evaluate whether the name is a package, class, or class::method
         * combination.  This parser supports inputs of the form:
         *
         *    *::*::*         - ALL timers added
         *    Package::*::*   - "Package" added to package list.
         *    Class           - "Class" added to class list.
         *    *::Class        - "Class" added to class list.
         *    Class::*        - "Class" added to class list.
         *    *::Class::*     - "Class" added to class list.
         *    Package::Class::method  - "Class::method" put to class_method list
         *    Class::method   - "Class::method" put to class_method list
         */
        std::string::size_type position, string_length;

        /*
         *  Step thru the input list and form the d_packageNames,
         *  d_classNames, and d_classMethodNames lists.
         */
        if (!name.empty())
        {
            // Nested if #1
            std::string entry = name;
            /*
             *  Once we have determined whether the entry is a package,
             *  class, or class::method, use this bool to jump to the next
             *  loop entry.
             */
            bool determined_entry = false;

            /*
             * Check if its a wildcard entry - "*::*::*".  If so, add all package
             * names to the package name list.
             */
            position = entry.find("*::*::*");  // if not found, "position" runs off
            // end of entry so pos > entry.size()
            if (position < entry.size())
            {
                d_packageNames.push_front("algs");
                d_packageNames.push_front("apps");
                d_packageNames.push_front("appu");
                d_packageNames.push_front("geom");
                d_packageNames.push_front("hier");
                d_packageNames.push_front("math");
                d_packageNames.push_front("mesh");
                d_packageNames.push_front("pdat");
                d_packageNames.push_front("solv");
                d_packageNames.push_front("tbox");
                d_packageNames.push_front("xfer");
                determined_entry = true;
            }

            /*
             * Is it a package?  Look for "::*::*" string.  If its there,
             * parse it off and add the package to the package list.
             */
            if (!determined_entry)
            {
                position = entry.find("::*::*");
                if (position < entry.size())
                {
                    entry = entry.substr(0, position);
                    d_packageNames.push_front(entry);
                    determined_entry = true;
                }
            }

            if (!determined_entry)
            {
                /*
                 * Is it a class?  If it doesn't have any "::", it must be a class.
                 */
                position = entry.find("::");
                if (position > entry.size())
                {
                    d_classNames.push_front(entry);
                    determined_entry = true;
                }
                if (!determined_entry)
                {
                    /*
                     * At this point, we know the entry has a "::" but wasn't
                     * identified as a package.  There are several options that
                     * might make Foo a class entry:
                     *  1) Foo::
                     *  2) *::Foo::
                     *  3) Package::Foo::
                     *  4) *::Foo
                     * Parse these as follows:  First, look for existence of "::*"
                     * at the end of the entry.  This will identify the first 3
                     * options.  Next look for existence of "*::" at front of the
                     * string.  This will identify the fourth choice.
                     *
                     * Check for first three options...
                     */
                    string_length = entry.size();
                    std::string substring = entry.substr(string_length - 3, string_length);
                    if (substring == "::*")
                    {
                        entry = entry.substr(0, string_length - 3);

                        /*
                         * If a preceeding "::" exists at the front of the entry
                         * (i.e. option 2 and 3), parse off anything before it.
                         */
                        position = entry.find("::");
                        if (position < entry.size())
                        {
                            entry = entry.substr(position + 2);
                        }
                        d_classNames.push_front(entry);
                        determined_entry = true;
                    }

                    if (!determined_entry)
                    { 
                        /*
                         * Check for option 4.  The entry has a preceeding *::. Do not
                         * accept case where there is a second "::" followed by anything
                         * but "*", since this is a class::method combination.
                         *
                         */
                        substring = entry.substr(0, 3);
                        if (substring == "*::")
                        {
                            entry = entry.substr(3);
                            position = entry.find("::");

                            /*
                             * There is no second "::".  Accept the entry as a class.
                             */
                            if (position > entry.size())
                            {
                                d_classNames.push_front(entry);
                                determined_entry = true;
                            }
                            else
                            {
                                /*
                                 * There is a second "::".  See if it is followed by a
                                 * "*".  If so, parse off the "::*" and accept entry as
                                 * a class.  If not, let it be determined below to be a
                                 * class::method entry.
                                 */
                                string_length = entry.size();
                                substring = entry.substr(string_length - 1, string_length);
                                if (substring == "*")
                                {
                                    entry = entry.substr(0, string_length - 3);
                                    d_classNames.push_front(entry);
                                    determined_entry = true;
                                }
                            }
                        }

                        if (!determined_entry)
                        {
                            /*
                             * The entry has not been identified as either a package or
                             * a class.  It must be a class::method combination. There
                             * are three options for entering class::method combinations:
                             *  1) Package::Foo::method
                             *  2) *::Foo::method
                             *  3) Foo::method
                             * We only want to maintain "Foo::method" in the package
                             * list.  Check first if there are two "::" in the entry.
                             * If there are, parse of whatever is in front of the
                             * first "::".  If not, just use the entry as is.
                             */
                            position = entry.find("::");
                            if (position < entry.size())
                            {
                                substring = entry.substr(position + 2);
                                position = substring.find("::");
                                if (position < substring.size())
                                {
                                    /*
                                     * There *are* two "::" so entry must contain a
                                     * package.  Parse it off.
                                     */
                                    position = entry.find("::");
                                    entry = entry.substr(position + 2);
                                }
                            }
                            d_classMethodNames.push_front(entry);
                        }
                    } 
                }
            }
        } 
    }
    
    void TimerManager::Quicksort(const std::vector<double>& a, int index[], int lo, int hi)
    {
        if (hi <= lo)
            return;

        /*
         * Put a[i] into position for i between lo and hi
         * (i.e. pivot point)
         */
        int i = lo - 1;
        int j = hi;
        double v = a[index[hi]];
        for ( ; ; )
        {
            while (a[index[++i]] > v)
                NULL_STATEMENT;
            while (v > a[index[--j]])
            {
                if (j == lo) break;
            }
            if (i >= j) break;

            // exchange i, j indices
            int temp = index[i];
            index[i] = index[j];
            index[j] = temp;
        }
        // exchange i, hi indices
        int temp = index[i];
        index[i] = index[hi];
        index[hi] = temp;

        Quicksort(a, index, lo, i - 1);
        Quicksort(a, index, i + 1, hi);
    }

    int TimerManager::ComputePercentageInt(const double frac, const double tot)
    {
        /*
         *  Put a cap on the percentage at 1000.  If tot = 0, this if
         *  test should catch it.
         */
        int perc = 0;
        if (tot > 0.1 * frac)
        {
            perc = int(frac / tot * 100.);
        }
        else
        {
            perc = 1000;
        }
        return perc;
    }

    double TimerManager::ComputePercentageDouble(const double frac, const double tot)
    {
        /*
         *  Put a cap on the percentage at 1000.  If tot = 0, this if
         *  test should catch it.
         */
        double perc = 0;
        if (tot > 0.1 * frac)
        {
            perc = frac / tot * 100.;
        }
        else
        {
            perc = 1000;
        }
        return perc;
    }

    void TimerManager::ComputeOverheadConstants()
    {
        if (d_timerActiveAccessTime < 0.0)
        {
            ClearArrays();
            d_timerActiveAccessTime  = ComputeOverheadConstantActiveOrInactive(false);

            ClearArrays();
            d_timerInactiveAccessTime  = ComputeOverheadConstantActiveOrInactive(true);

            ClearArrays();
        }
    }
   
    double TimerManager::ComputeOverheadConstantActiveOrInactive(bool active)
    {
        std::string outer_name("TimerManger::Outer");
        Timer* outer_timer(TimerManager::GetInstance()->GetTimer(outer_name, true));

        std::string inner_name("TimerMangerInner");
        Timer* inner_timer(TimerManager::GetInstance()->GetTimer(inner_name, active));

        const int ntest = 1000;
        for (int i = 0; i < ntest; ++i)
        {
            outer_timer->Start();
            inner_timer->Start();
            inner_timer->Stop();
            outer_timer->Stop();
        }

        return (outer_timer->GetTotalWallclockTime() - inner_timer->GetTotalWallclockTime())/(static_cast<double>(ntest));
    }

    void TimerManager::ClearArrays()
    {
        /*
         * Create a timer that measures overall solution time.  If the
         * application uses Tau, this timer will effectively measure
         * uninstrumented parts of the library.  Hence, use a different name
         * for the different cases to avoid confusion in the Tau analysis tool.
         */
        //d_mainTimer->Reset(new Timer("TOTAL RUN TIME"));
        d_mainTimer = new Timer("TOTAL RUN TIME");

        d_timers.clear();
        d_inactiveTimers.clear();

        d_exclusiveTimerStack.clear();
    }
    
    void TimerManager::FinalizeCallback()
    {
        if (d_instance)
        {
            delete d_instance;
            d_instance = NULL;
        }
    }
} // end namespace ARIES

