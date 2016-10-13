/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for triangle grid definiton
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    02-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridTriangle.hpp"

namespace ARIES
{
    GridTriangle::GridTriangle(unsigned long val_point_0, unsigned long val_point_1, unsigned long val_point_2, unsigned short val_nDim) : Grid()
    {
        d_face[0][0] = 0;
        d_face[0][1] = 1;
        d_face[1][0] = 1;
        d_face[1][1] = 2;
        d_face[2][0] = 0;
        d_face[2][1] = 2;

        d_neighborNode[0][0] = 1;
        d_neighborNode[0][1] = 2;
        d_neighborNode[1][0] = 2;
        d_neighborNode[1][1] = 0;
        d_neighborNode[2][0] = 0;
        d_neighborNode[2][1] = 1;
        
        d_numNodeFace[0] = 2;
        d_numNodeFace[0] = 2;
        d_numNodeFace[0] = 2;
        
        d_numNeighborNode[0] = 2;
        d_numNeighborNode[1] = 2;
        d_numNeighborNode[2] = 2;

        d_numFace = 3;
        d_numNode = 3;
        d_numNeighborElement = 3;
        d_VTKType = 5;
        d_maxNodeFace = 2;
        
        d_numDim = val_nDim;
        for (unsigned short iDim = 0; iDim < d_numDim; iDim++)
            d_coordCG.push_back(0.0);
        
        for (unsigned short iFace = 0; iFace < d_numFace; iFace++)
        {
            vector<double> coord;
            for (unsigned short iDim = 0; iDim < d_numDim; iDim++)
                coord.push_back(0.0);
            d_coordFaceElemCG.push_back(coord);
        }

        d_node.push_back(val_point_0);
        d_node.push_back(val_point_1);
        d_node.push_back(val_point_2);
 
        d_numNeighborElement = d_numFace;
        for (unsigned short iNeighborElement = 0; iNeighborElement < d_numNeighborElement; iNeighborElement++)
        {
            d_neighborElement.push_back(-1);
        }
    }

    GridTriangle::~GridTriangle()
    {
    }

    void GridTriangle::ChangeOrientation()
    {
        unsigned long iPoint, Point_2;
        iPoint =  d_node[0];
        Point_2 = d_node[2];
        d_node[0] = Point_2;
        d_node[2] = iPoint;
    }
    
}

