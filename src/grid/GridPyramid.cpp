/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for pyramid grid definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GridPyramid.hpp"

#include <iostream>

namespace ARIES
{

        unsigned short GridPyramid::d_face[5][4] = { { 0, 3, 2, 1 }, { 4, 3, 0, 0 }, { 4, 0, 1, 1 }, { 2, 4, 1, 1 }, { 3, 4, 2, 2 } };
        unsigned short GridPyramid::d_neighborNode[5][4] = { { 1, 3, 4, 4 }, { 0, 2, 4, 4 }, { 1, 3, 4, 4 }, { 2, 0, 4, 4 }, { 0, 1, 2, 3 } };
        unsigned short GridPyramid::d_numNodeFace[5] = { 4, 3, 3, 3, 3 };
        unsigned short GridPyramid::d_numNeighborNode[5] = { 3, 3, 3, 3, 4 };
        unsigned short GridPyramid::d_numFace = 5;
        unsigned short GridPyramid::d_numNode = 5;
        unsigned short GridPyramid::d_numNeighborElement = 5;
        unsigned short GridPyramid::d_VTKType = 14;
        unsigned short GridPyramid::d_maxNodeFace = 4;

    GridPyramid::GridPyramid(unsigned long val_point_0, unsigned long val_point_1,
                             unsigned long val_point_2, unsigned long val_point_3,
                             unsigned long val_point_4)
            :Grid()
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

        for (unsigned short iNeighborElement = 0; iNeighborElement < d_numNeighborElement; iNeighborElement++)
        {
            d_neighborElement.push_back(-1);
        }

    }

    GridPyramid::~GridPyramid()
    {
    }

    void GridPyramid::ChangeOrientation()
    {
        std::cout << "Not defined orientation change" << std::endl;
    }
    
}
