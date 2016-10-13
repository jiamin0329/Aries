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

#ifndef ARIES_PARAM_HPP
#define ARIES_PARAM_HPP

#include <string>
#include <vector>
#include <map>

namespace ARIES
{
    class Param
    {
    public:
        Param();
        ~Param();

        static Param* GetParamObject();




    private:
        static Param* d_paramObject;
        
        map<int, vector<int>    > d_paramIntTable;
        map<int, vector<double> > d_paramDoubleTable;
        map<int, vector<string> > d_paramStringTable;
        map<int, vector<bool>   > d_paramBoolTable;
    };
}

#endif

