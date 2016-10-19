/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Aries param indices
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    19-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "ParamIndices.hpp"

#include <vector>

namespace 
{
    static ParamIndexTable table[] =
    {
        {PARAM_aries_component, NULL},
        {PARAM_cfd_solver_type, NULL}
    };

    static DefaultValueBool defBoolValueTable[] =
    {
        {PARAM_timemgr_print_exclusive,     true},
        {PARAM_timemgr_print_total,         true},
        {PARAM_timemgr_print_processor,     true},
        {PARAM_timemgr_print_max,           true},
        {PARAM_timemgr_print_summed,        true},
        {PARAM_timemgr_print_user,          true},
        {PARAM_timemgr_print_sys,           true},
        {PARAM_timemgr_print_wall,          true},
        {PARAM_timemgr_print_percentage,    true},
        {PARAM_timemgr_print_concurrent,    true},
        {PARAM_timemgr_print_tiemroverhead, true},
        {PARAM_timemgr_print_threshold,     true}
    };
    
   
}

namespace ARIES
{

    void InsertParam()
    {
        
        
        for (int paramIndex = PARAM_paramindices_start; paramIndex < PARAM_paramindices_end; ++paramIndex)
        {
            table[paramIndex].ParamIndexInsertCallback;
        }


    }










    
}





