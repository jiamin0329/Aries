/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Parallel I/O classes pout, perr, and plog and control class
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    05-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "Logger.hpp"

/* Aries includes */
#include "AriesMPI.hpp"
#include "Utilities.hpp"
/* c++ includes */
#include <string>

/*
 *================================================================================
 *    Data definitions
 *================================================================================
 */

/*
 *================================================================================
 *    Static functions
 *================================================================================
 */

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    Logger* Logger::d_instance = NULL;
    StartupShutdownManager::Handler Logger::d_finalizeHandler(0, 0, 0, Logger::FinalizeCallback, StartupShutdownManager::priorityLogger);

    Logger* Logger::GetInstance()
    {
        if (d_instance == NULL)
        {
            d_instance = new Logger();
        }

        return d_instance;
    }

    void Logger::Startup(bool val_isScrActive, bool val_isLogActive, const std::string& val_filename)
    {
        d_isScrActive = val_isScrActive;
        d_isLogActive = val_isLogActive;

        if (d_isLogActive && !val_filename.empty())
        {
            std::string fullFilename;
            if (d_rank == 0)
                fullFilename = val_filename;
            else
                fullFilename = val_filename + "." + Utilities::ProcessorToString(d_rank);
        
            d_filestream = new std::ofstream(fullFilename.c_str());
        
            if (!(*d_filestream))
            {
                std::cout << "Cannot open log file!" << std::endl;
                AriesMPI::Abort();
            }
        }
    }

    void Logger::Shutdown()
    {
        d_isScrActive = false;
        d_isLogActive = false;
        if (*d_filestream)
            d_filestream->close();
    }

    Logger::Logger():
            d_rank(-1),
            d_isScrActive(false),
            d_isLogActive(false)
    {
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        mpi.Comm_rank(&d_rank);
        d_filestream = 0;
        d_prefix = "Processor= " + Utilities::ProcessorToString(d_rank) + ":";
    }
    
    Logger::~Logger()
    {
        if(d_filestream != NULL)
        {
            delete d_filestream;
            d_filestream = NULL;
        }
    }
 
    void Logger::FinalizeCallback()
    {
        if (d_instance)
        {
            delete d_instance;
            d_instance = NULL;
        }
    }

    void Logger::SuspendLogging()
    {
        d_isLogActive = false;
    }

    void Logger::ResumeLogging()
    {
        if (*d_filestream)
        {
            d_isLogActive = true;
        }
    }

    /*!
     * Logs an abort message with file & location
     */
    void Logger::DoLog(const Logger::LoggerType val_loggerType, const std::string& message, const std::string& filename, const int line)
    {
        std::string messageType = "INFO MESSAGE: ";

        switch(val_loggerType)
        {
            case Logger_Debug:
                messageType = "INFO MESSAGE: ";
                break;
            case Logger_Info:
                messageType = "DEBUG MESSAGE: ";
                break;
            case Logger_Warning:
                messageType = "WARNING MESSAGE: ";
                break;
            case Logger_Error:
                messageType = "ERROR MESSAGE: ";
                break;
            case Logger_Fatal:
                messageType = "FATAL ERROR MESSAGE: ";
                break;
            default:
                messageType = "INFO MESSAGE: ";
                break;
        }

        if (d_rank > 0)
            messageType = d_prefix + messageType;
            
        if (d_isScrActive)
        {
            std::cout << "***** in file: " << filename << std::endl;
            std::cout << "***** at line: " << line << std::endl;
            std::cout << messageType << std::endl << message.c_str() << std::endl;
            std::cout << std::flush;
        }

        if (d_isLogActive && *d_filestream)
        {
            *d_filestream << "***** in file: " << filename << std::endl;
            *d_filestream << "***** at line: " << line << std::endl;
            *d_filestream << messageType << std::endl << message.c_str() << std::endl;
            *d_filestream << std::flush;
        }

        if (val_loggerType == Logger_Fatal)
            AriesMPI::Abort();
    }

    
}
