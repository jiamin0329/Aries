/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for prism grid definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    03-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridPrism.hpp"

namespace ARIES
{
    unsigned short GridPrism::d_face[5][4] = { { 3, 4, 1, 0 }, { 5, 2, 1, 4 }, { 2, 5, 3, 0 }, { 0, 1, 2, 2 }, { 5, 4, 3, 3 } };
    unsigned short GridPrism::d_neighborNode[6][3] = { { 1, 2, 3 }, { 0, 2, 4 }, { 1, 0, 5 }, { 0, 4, 5 }, { 3, 5, 1 }, { 4, 3, 2 } };
    unsigned short GridPrism::d_numNodeFace[5] = { 4, 4, 4, 3, 3 };
    unsigned short GridPrism::d_numNeighborNode[6] = { 3, 3, 3, 3, 3, 3 };
    unsigned short GridPrism::d_numFace = 5;
    unsigned short GridPrism::d_numNode = 6;
    unsigned short GridPrism::d_numNeighborElement= 5;
    unsigned short GridPrism::d_VTKType = 13;
    unsigned short GridPrism::d_maxNodeFace = 4;

    GridPrism::GridPrism(unsigned long val_point_0, unsigned long val_point_1,
                         unsigned long val_point_2, unsigned long val_point_3,
                         unsigned long val_point_4, unsigned long val_point_5) : Grid()
    {
        unsigned short iDim, iFace, iNeighbor_Elements;

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
        d_node.push_back(val_point_4);
        d_node.push_back(val_point_5);
       
        for (unsigned short iNeighborElement = 0; iNeighborElement < d_numNeighborElement; iNeighborElement++)
        {
            d_neighborElement.push_back(-1);
        }

    }

    GridPrism::~GridPrism()
    {
    }

    void GridPrism::ChangeOrientation(void)
    {
        unsigned long Point_0, Point_1, Point_3, Point_4;
        Point_0 = d_node[0];
        Point_1 = d_node[1];
        Point_3 = d_node[3];
        Point_4 = d_node[4];

        d_node[0] = Point_1;
        d_node[1] = Point_0;
        d_node[3] = Point_4;
        d_node[4] = Point_3;
    }
}
