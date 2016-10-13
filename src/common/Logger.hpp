/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Utility class for logging.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    05-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_LOGGER_HPP
#define ARIES_LOGGER_HPP

/* Aries includes */
#include "StartupShutDownManager.hpp"

/* C++ includes */
#include <string>
#include <iostream>

/*
 *================================================================================
 *    Forward declarations
 *================================================================================
 */
#define ARIES_DEBUG(X)                          \
    do {                                        \
        std::ostringstream ariesos;             \
        ariesos << X << std::ends;              \
        ARIES::Logger::GetInstance()->DoLog(    \
            ARIES::Logger::Logger_Debug,        \
            ariesos.str(), __FILE__, __LINE__); \
    } while (0)

#define ARIES_INFO(X)                           \
    do {                                        \
        std::ostringstream ariesos;             \
        ariesos << X << std::ends;              \
        ARIES::Logger::GetInstance()->DoLog(    \
            ARIES::Logger::Logger_Info,         \
            ariesos.str(), __FILE__, __LINE__); \
    } while (0)

#define ARIES_WARNING(X)                        \
    do {                                        \
        std::ostringstream ariesos;             \
        ariesos << X << std::ends;              \
        ARIES::Logger::GetInstance()->DoLog(    \
            ARIES::Logger::Logger_Warning,      \
            ariesos.str(), __FILE__, __LINE__); \
    } while (0)

#define ARIES_ERROR(X)                          \
    do {                                        \
        std::ostringstream ariesos;             \
        ariesos << X << std::ends;              \
        ARIES::Logger::GetInstance()->DoLog(    \
            ARIES::Logger::Logger_Error,        \
            ariesos.str(), __FILE__, __LINE__); \
    } while (0)

#define ARIES_FATAL(X)                          \
    do {                                        \
        std::ostringstream ariesos;             \
        ariesos << X << std::ends;              \
        ARIES::Logger::GetInstance()->DoLog(    \
            ARIES::Logger::Logger_Fatal,        \
            ariesos.str(), __FILE__, __LINE__); \
    } while (0)

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    class Logger
    {
    public:
        typedef enum
        {
            Logger_None,
            Logger_Debug,
            Logger_Info,
            Logger_Warning,
            Logger_Error,
            Logger_Fatal
        } LoggerType;

        static Logger* GetInstance();

        void Startup(bool val_scrActive, bool val_logActive, const std::string& val_filename);
        void Shutdown();

        void SuspendLogging();
        void ResumeLogging();

        void DoLog(const LoggerType val_loggerType, const std::string& message, const std::string& filename, const int line);

    private:
        Logger();
        Logger(const Logger& other);
        ~Logger();
        Logger& operator = (const Logger& rhs);

        /*!
         * Frees instance of the singleton logger. NOTE: should be called by StartupShutdownManager only.
         */
        static void FinalizeCallback();

        static Logger* d_instance;          // instance of singleton
        static StartupShutdownManager::Handler d_finalizeHandler;
        int d_rank;                  // processor rank in MPI group
        std::ofstream* d_filestream;  // NULL or log filestream
       
        std::string d_prefix;
        
        bool d_isLogActive;
        bool d_isScrActive;
    };
}

#endif










