/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Main classes for defining the dual grid
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    24-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "DualGrid.hpp"

namespace ARIES
{
    unsigned short DualGrid::d_nDim = 0;
        
    DualGrid::DualGrid(unsigned short val_nDim) 
    { 
        d_nDim = val_nDim; 
    }

    DualGrid::~DualGrid() 
    {}
}
