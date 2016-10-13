/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for tetrahedron grid definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    03-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridTetra.hpp"

namespace ARIES
{

    unsigned short GridTetra::d_face[4][3] = { { 0, 2, 1 }, { 0, 1, 3 }, { 0, 3, 2 }, { 1, 2, 3 } };
    unsigned short GridTetra::d_neighborNode[4][3] = { { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } };
    unsigned short GridTetra::d_numNodeFace[4] = { 3, 3, 3, 3 };
    unsigned short GridTetra::d_numNeighborNode[4] = { 3, 3, 3, 3 };
    unsigned short GridTetra::d_numFace = 4;
    unsigned short GridTetra::d_numNode = 4;
    unsigned short GridTetra::d_numNeighborElement = 4;
    unsigned short GridTetra::d_VTKType = 10;
    unsigned short GridTetra::d_maxNodeFace = 3;

    GridTetra::GridTetra(unsigned long val_point_0, unsigned long val_point_1,
                         unsigned long val_point_2, unsigned long val_point_3): Grid()
    {
        d_numDim = 3;
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
        d_node.push_back(val_point_3);

        d_numNeighborElement = d_numFace;
        for (unsigned short iNeighborElement = 0; iNeighborElement < d_numNeighborElement; iNeighborElement++)
        {
            d_neighborElement.push_back(-1);
        }
    }

    GridTetra::~GridTetra()
    {
    }

    void GridTetra::ChangeOrientation()
    {
        unsigned long jPoint, Point_1;
        jPoint  = d_node[0];
        Point_1 = d_node[1];
        d_node[0] = Point_1;
        d_node[1] = jPoint;
    }
    
}
