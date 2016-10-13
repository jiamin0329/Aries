/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for line grid definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    02-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridLine.hpp"

namespace ARIES
{
    GridLine::GridLine(unsigned long val_point_0, unsigned long val_point_1, unsigned short val_nDim) : Grid()
    {
        d_face[0][0] = 0;
        d_face[0][1] = 1;
        
        d_neighborNode[0][0] = 1;
        d_neighborNode[1][0] = 0;

        d_numNodeFace[0] = 2;
        
        d_numNeighborNode[0] = 1;
        d_numNeighborNode[1] = 1;
        d_numFace = 1;
        d_numNode = 2;
        d_numNeighborElement = 1;
        d_VTKType = 3;
        d_maxNodeFace = 2;

        d_numDim = val_nDim;
        for (unsigned short iDim = 0; iDim < d_numDim; iDim++)
            d_coordCG.push_back(0.0);
        for (unsigned short iFace = 0; iFace < d_numFace; iFace++)
        {
            for (unsigned short iDim = 0; iDim < d_numDim; iDim++)
                d_coordFaceElemCG[iFace][iDim] = 0.0;
        }

        d_node.push_back(val_point_0);
        d_node.push_back(val_point_1);
    }

    GridLine::~GridLine()
    {
    }

    void GridLine::ChangeOrientation(void)
    {
        unsigned long iPoint, jPoint;

        iPoint = d_node[0];
        jPoint = d_node[1];
        d_node[0] = jPoint;
        d_node[1] = iPoint;
    }
}
