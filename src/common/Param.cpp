
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

#include "Param.hpp"

using namespace std;

namespace ARIES
{

    Param()
    {
        
    }

    ~Param()
    {
        
    }

    Param* Param::GetParamObject()
    {
        if(!d_paramObject)
        {
            d_paramObject = new Param();
        }

        return d_paramObject;
    }

    int Param::AskIntValue(int val_paramIndex)
    {
        int value = 0;
        ParamIntTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamIntTable.end() )
            value = paramIter->second[0];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }

    int Param::AskIntDefaultValue(int val_paramIndex)
    {
        int value = 0;
        ParamIntTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamIntTable.end() )
            value = paramIter->second[1];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }

    
    double Param::AskDoubleValue(int val_paramIndex)
    {
        int value = 0;
        ParamDoubleTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamDoubleTable.end() )
            value = paramIter->second[0];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }
    
    double Param::AskDoubleDefaultValue(int val_paramIndex)
    {
        int value = 0;
        ParamDoubleTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamDoubleTable.end() )
            value = paramIter->second[1];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }

    string Param::AskStringValue(int val_paramIndex)
    {
        int value = 0;
        ParamStringTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamStringTable.end() )
            value = paramIter->second[0];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }

    string Param::AskStringDefaultValue(int val_paramIndex)
    {
        int value = 0;
        ParamStringTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamStringTable.end() )
            value = paramIter->second[1];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }

    bool Param::AskBoolValue(int val_paramIndex)
    {
        int value = 0;
        ParamBoolTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamBoolTable.end() )
            value = paramIter->second[0];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }

    bool Param::AskBoolDefaultValue(int val_paramIndex)
    {
        int value = 0;
        ParamBoolTable::iterator paramIter;
        paramIter.find(val_paramIndex);
        if ( paramIter != ParamBoolTable.end() )
            value = paramIter->second[1];
        else
            cout<< paramIter->first << " not found! " <<endl;
        
        return value;
    }
    
}
