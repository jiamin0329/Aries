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

#include "GridHexa.hpp"

namespace ARIES
{

    unsigned short GridHexa::d_face[6][4] = { { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 }, { 0, 3, 2, 1 }, { 4, 5, 6, 7 } };
    unsigned short GridHexa::d_neighborNode[8][3] = { { 1, 3, 4 }, { 0, 2, 5 }, { 1, 3, 6 }, { 0, 2, 7 }, { 0, 5, 7 }, { 4, 6, 1 }, { 2, 5, 7 }, { 4, 3, 6 } };
    unsigned short GridHexa::d_numNodeFace[6] = { 4, 4, 4, 4, 4, 4 };
    unsigned short GridHexa::d_numNeighborNode[8] = { 3, 3, 3, 3, 3, 3, 3, 3 };
    unsigned short GridHexa::d_numFace = 6;
    unsigned short GridHexa::d_numNode = 8;
    unsigned short GridHexa::d_numNeighborElement = 6;
    unsigned short GridHexa::d_VTKType = 12;
    unsigned short GridHexa::d_maxNodeFace = 4;

    GridHexa::GridHexa(unsigned long val_point_0, unsigned long val_point_1,
                       unsigned long val_point_2, unsigned long val_point_3,
                       unsigned long val_point_4, unsigned long val_point_5,
                       unsigned long val_point_6, unsigned long val_point_7) : Grid() 
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
        d_node.push_back(val_point_4);
        d_node.push_back(val_point_5);
        d_node.push_back(val_point_6);
        d_node.push_back(val_point_7);

        d_numNeighborElement = d_numFace;
        for (unsigned short iNeighborElement = 0; iNeighborElement < d_numNeighborElement; iNeighborElement++) 
        {
            d_neighborElement.push_back(-1);
        }
    }

    GridHexa::~GridHexa() 
    {
    }

    void GridHexa::ChangeOrientation() 
    {
        unsigned long Point_0, Point_1, Point_2, Point_3, Point_4, Point_5, Point_6, Point_7;
        Point_0 = d_node[0];
        Point_1 = d_node[1];
        Point_2 = d_node[2];
        Point_3 = d_node[3];
        Point_4 = d_node[4];
        Point_5 = d_node[5];
        Point_6 = d_node[6];
        Point_7 = d_node[7];

        d_node[0] = Point_7;
        d_node[1] = Point_4;
        d_node[2] = Point_5;
        d_node[3] = Point_6;
        d_node[4] = Point_3;
        d_node[5] = Point_0;
        d_node[6] = Point_1;
        d_node[7] = Point_2;
    }
}
