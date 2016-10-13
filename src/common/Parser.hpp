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
 *    04-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_PARSER_HPP
#define ARIES_PARSER_HPP

#include <string>

using namespace std;

namespace ARIES
{

    class Parser
    {
    public:
        Parser(string val_configFile);
        ~Parser();


        bool DoParse();

        
        inline void StringToUpperCase(string & str)
        {
            std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        }


        inline string StringToUpperCase(const string & str)
        {
            string upp_str(str);
            std::transform(upp_str.begin(), upp_str.end(), upp_str.begin(), ::toupper);
            return upp_str;
        }


        
    private:
        string d_configFile;

    };
}

#endif
