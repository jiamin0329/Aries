/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for quadrilateral
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    03-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridQuad.hpp"

namespace ARIES
{
    unsigned short GridQuad::d_face[4][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
    unsigned short GridQuad::d_neighborNode[4][2] = { { 1, 3 }, { 2, 0 }, { 3, 1 }, { 0, 2 } };
    unsigned short GridQuad::d_numNodeFace[4] = { 2, 2, 2, 2 };
    unsigned short GridQuad::d_numNeighborNode[4] = { 2, 2, 2, 2 };
    unsigned short GridQuad::d_numFace = 4;
    unsigned short GridQuad::d_numNode = 4;
    unsigned short GridQuad::d_numNeighborElement = 4;
    unsigned short GridQuad::d_VTKType = 9;
    unsigned short GridQuad::d_maxNodeFace = 2;
    
    GridQuad::GridQuad(unsigned long val_point_0, unsigned long val_point_1, unsigned long val_point_2, unsigned long val_point_3, unsigned short val_nDim) : Grid()
    {
        d_numDim = val_nDim;
        for (unsigned short iDim = 0; iDim <  d_numDim; iDim++)
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
        d_node.push_back(val_point_3);

        d_numNeighborElement = d_numFace;
        for (unsigned short iNeighborElement = 0; iNeighborElement < d_numNeighborElement; iNeighborElement++)
        {
            d_neighborElement.push_back(-1);
        }

    }

    GridQuad::~GridQuad()
    {
    }
    
    void GridQuad::ChangeOrientation()
    {
        unsigned long jPoint, Point_3;
        jPoint  = d_node[1];
        Point_3 = d_node[3];
        d_node[1] = Point_3;
        d_node[3] = jPoint;
    }
    
}


