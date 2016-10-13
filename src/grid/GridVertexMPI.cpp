/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for vertex grid
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridVertexMPI.hpp"

#include <iostream>

namespace ARIES
{
    GridVertexMPI::GridVertexMPI(unsigned long val_point, unsigned short val_nDim) : Grid() 
    {
        d_numNode = 1;
        d_numFace = 0;
        d_maxNodeFace = 0;
        d_numNeighborElement = 0;
        d_VTKType = 1;
        d_rotationType = 0;            // By default, no rotation in the solution

        d_numDim = val_nDim;
        for (unsigned short iDim = 0; iDim < d_numDim; iDim++)
            d_coordCG.push_back(0.0);

        d_node.push_back(val_point);
    }

    GridVertexMPI::~GridVertexMPI() 
    {
    }

    void GridVertexMPI::ChangeOrientation()
    { 
        std::cout << "Not defined orientation change" << std::endl; 
    }
}










