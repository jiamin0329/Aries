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

#ifndef ARIES_UTILITES_HPP
#define ARIES_UTILITES_HPP

//#include "Logger.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

/*!
 * A statement that does nothing, for insure++ make it something
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if (0) int nullstatement = 0
#else
#define NULL_STATEMENT
#endif

/*!
 * A null use of a variable, use to avoid GNU compiler
 * warnings about unused variables.
 */
#define NULL_USE(variable)                                  \
    do {                                                    \
        if (0) { char* temp = (char *)&variable; ++temp; }  \
    } while (0)

/*!
 * Throw an error assertion from within any C++ source code if the
 * given expression is not true.  This is a parallel-friendly version
 * of assert.
 * The file and line number of the abort are also printed.
 */
#ifdef DEBUG_CHECK_ASSERTIONS

#define ARIES_ASSERT(EXP)                                               \
    do                                                                  \
    {                                                                   \
        if (!(EXP))                                                     \
        {                                                               \
            std::ostringstream ariesos;                                 \
            ariesos << "Failed assertion: " << # EXP << std::ends;      \
            ARIES::Utilities::Abort(ariesos.str(), __FILE__, __LINE__); \
        }                                                               \
    } while (0)
#else

/*
 * No assertion checking
 */
#define ARIES_ASSERT(EXP)

#endif

namespace ARIES
{
    /*!
     * Utilities is a Singleton class containing basic routines for error
     * reporting, file manipulations, etc.
     */
    struct Utilities
    {
        /*!
         * Create the directory specified by the path string.  Permissions are set
         * by default to rwx by user.  The intermediate directories in the
         * path are created if they do not already exist.  When
         * only_node_zero_creates is true, only node zero creates the
         * directories.  Otherwise, all nodes create the directories.
         */
        static void RecursiveMkdir(const std::string& path, /*mode_t mode = (S_IRUSR | S_IWUSR | S_IXUSR),*/ bool only_node_zero_creates = true);

        /*!
         * Rename a file from old file name to new file name.
         *
         * @pre !old_filename.empty()
         * @pre !new_filename.empty()
         */
        static void RenameFile(const std::string& old_filename, const std::string& new_filename)
        {
            ARIES_ASSERT(!old_filename.empty());
            ARIES_ASSERT(!new_filename.empty());
            rename(old_filename.c_str(), new_filename.c_str());
        };

        /*!
         * Convert an integer to a string.
         *
         * The returned string is padded with zeros as needed so that it
         * contains at least the number of characters indicated by the
         * minimum width argument.  When the number is positive, the
         * string is padded on the left. When the number is negative,
         * the '-' sign appears first, followed by the integer value
         * padded on the left with zeros.  For example, the statement
         * intToString(12, 5) returns "00012" and the statement
         * intToString(-12, 5) returns "-0012".
         */
        static std::string IntToString(int num, int min_width = 1);

        /*!
         * Convert a size_t to a string.
         *
         * The returned string is padded with zeros as needed so that it
         * contains at least the number of characters indicated by the
         * minimum width argument.  When the number is positive, the
         * string is padded on the left. When the number is negative,
         * the '-' sign appears first, followed by the integer value
         * padded on the left with zeros.  For example, the statement
         * intToString(12, 5) returns "00012" and the statement
         * intToString(-12, 5) returns "-0012".
         */
        static std::string SizetToString(size_t num, int min_width = 1);

        /*!
         * Convert common integer values to strings.
         *
         * These are simply wrappers around intToString that ensure the
         * same width is uniformally used when converting to string
         * representations.
         */
        static std::string NodeToString(int num) { return IntToString(num, d_nodeWidth); };
        static std::string ProcessorToString(int num) { return IntToString(num, d_processorWidth); };

        /*!
         * Aborts the run after printing an error message with file and
         * linenumber information.
         */
        static void Abort(const std::string& message, const std::string& filename, const int line);

    private:
        /*
         * Sizes for converting integers to fixed width strings
         * for things like filenames etc.
         */
        static const int d_nodeWidth = 7;
        static const int d_processorWidth = 7;
    };
} // end namespace ARIES

#endif



