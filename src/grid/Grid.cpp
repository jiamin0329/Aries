/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for primal grid definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    02-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GRID.hpp"

#include <iostream>

namespace ARIES
{
    Grid::Grid(void) 
    {
        d_node.clear();
        d_neighborElement.clear();
        d_coordCG.clear();
        d_coordFaceElemCG.clear();
    }

    Grid::~Grid() 
    {
    }

    void Grid::SetCG(double **val_coord)
    {
        unsigned short iDim, iNode, nodeFace, iFace;
        
        for (iDim = 0; iDim < d_numDim; iDim++)
        {
            d_coordCG[iDim] = 0.0;
            for (iNode = 0; iNode < GetNumNode(); iNode++)
                d_coordCG[iDim] += val_coord[iNode][iDim] / double(GetNumNode());
        }

        for (iFace = 0; iFace < GetNumFace(); iFace++)
        {
            for (iDim = 0; iDim < d_numDim; iDim++) 
            {
                d_coordFaceElemCG[iFace][iDim] = 0.0;
                for (iNode = 0; iNode < GetNumNodeFace(iFace); iNode++)
                {
                    nodeFace = GetFace(iFace, iNode);
                    d_coordFaceElemCG[iFace][iDim] += val_coord[nodeFace][iDim] / double(GetNumNodeFace(iFace));
                }
            }
        }
    }

    void Grid::GetAllNeighborElement() 
    {
        std::cout << "( ";
        for (unsigned short iFace = 0; iFace < GetNumFace(); iFace++)
        {
            std::cout << GetNeighborElement(iFace) << ", ";
        }
        std::cout << ")" << std::endl;
    }
    
}
