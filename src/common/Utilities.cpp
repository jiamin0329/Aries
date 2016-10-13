/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Utility functions for error reporting, file manipulation, etc.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "Utilities.hpp"

#include "AriesMPI.hpp"
#include "Logger.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <iomanip>

namespace ARIES
{
    /*
     * Routine to recursively construct directories based on a relative path name.
     */
    void Utilities::RecursiveMkdir(const std::string& path, /*mode_t mode,*/ bool only_node_zero_creates)
    {
#ifdef _MSC_VER
        const char seperator = '/';
#define mkdir(path, mode) mkdir(path)
#else
        const char seperator = '/';
#endif
        int mode = 0; // temp setting
        
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        if ((!only_node_zero_creates) || (mpi.GetRank() == 0))
        {
            int length = static_cast<int>(path.length());
            char* path_buf = new char[length + 1];
            sprintf(path_buf, "%s", path.c_str());
            struct stat status;
            int pos = length - 1;

            /* find part of path that has not yet been created */
            while ((stat(path_buf, &status) != 0) && (pos >= 0))
            {
                /* slide backwards in string until next slash found */
                bool slash_found = false;
                while ((!slash_found) && (pos >= 0))
                {
                    if (path_buf[pos] == seperator)
                    {
                        slash_found = true;
                        if (pos >= 0) path_buf[pos] = '\0';
                    }
                    else
                        --pos;
                }
            }

            /*
             * if there is a part of the path that already exists make sure
             * it is really a directory
             */
            if (pos >= 0)
            {
                if (!S_ISDIR(status.st_mode))
                {
                    ARIES_ERROR("Error in Utilities::recursiveMkdir...\n"
                               << "    Cannot create directories in path = " << path
                               << "\n    because some intermediate item in path exists and"
                               << "is NOT a directory" << std::endl);
                }
            }

            /* make all directories that do not already exist */

            /*
             * if (pos < 0), then there is no part of the path that
             * already exists.  Need to make the first part of the
             * path before sliding along path_buf.
             */
            if (pos < 0)
            {
                if (mkdir(path_buf, mode) != 0)
                {
                    ARIES_ERROR("Error in Utilities::recursiveMkdir...\n"
                               << "    Cannot create directory  = "
                               << path_buf << std::endl);
                }
                pos = 0;
            }

            /* make rest of directories */
            do
            {
                /* slide forward in string until next '\0' found */
                bool null_found = false;
                while ((!null_found) && (pos < length))
                {
                    if (path_buf[pos] == '\0')
                    {
                        null_found = true;
                        path_buf[pos] = seperator;
                    }
                    ++pos;
                }

                /* make directory if not at end of path */
                if (pos < length)
                {
                    if (mkdir(path_buf, mode) != 0)
                    {
                        ARIES_ERROR("Error in Utilities::recursiveMkdir...\n"
                                   << "    Cannot create directory  = "
                                   << path_buf << std::endl);
                    }
                }
            } while (pos < length);

            delete[] path_buf;
        }

        /*
         * Make sure all processors wait until node zero creates
         * the directory structure.
         */
        if (only_node_zero_creates)
        {
            AriesMPI::GetAriesWorld().Barrier();
        }
    }

    /*
     * Routine to convert an integer to a string.
     */
    std::string Utilities::IntToString(int num, int min_width)
    {
        int tmp_width = (min_width > 0 ? min_width : 1);
        std::ostringstream os;
        if (num < 0) {
            os << '-' << std::setw(tmp_width - 1) << std::setfill('0') << -num;
        } else {
            os << std::setw(tmp_width) << std::setfill('0') << num;
        }
        os << std::flush;

        return os.str();  //returns the string form of the stringstream object
    }

    /*
     * Routine to convert a size_t to a string.
     */
    std::string Utilities::SizetToString(size_t num, int min_width)
    {
        int tmp_width = (min_width > 0 ? min_width : 1);
        std::ostringstream os;
        os << std::setw(tmp_width) << std::setfill('0') << num;
        os << std::flush;

        return os.str();  //returns the string form of the stringstream object
    }

    /*
     * Routine that calls abort and prints calling location to error stream.
     */
    void Utilities::Abort(const std::string& message,const std::string& filename, const int line)
    {
        Logger::GetInstance()->DoLog(Logger::Logger_Fatal, message, filename, line);
    }
}


















