/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Accesses system times.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "Clock.hpp"
/* c++ includes */
#include <cstdlib>

namespace ARIES
{
    struct tms Clock::d_tmsBuffer;
    clock_t Clock::d_nullClock_t;
}






