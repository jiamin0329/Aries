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

#ifndef ARIES_GRIDLINE_HPP
#define ARIES_GRIDLINE_HPP

#include "Grid.hpp"

namespace ARIES
{

    class GridLine : public Grid
    {

    public:
        GridLine(unsigned long val_point_0, unsigned long val_point_1, unsigned short val_nDim);
        ~GridLine();

        unsigned long GetNode(unsigned short val_node) { return d_node[val_node]; };
        void SetNode(unsigned short val_node, unsigned long val_point) { d_node[val_node] = val_point; };

        unsigned short GetNumFace() { return d_numFace; };
        unsigned short GetNumNode() { return d_numNode; };

        unsigned short GetVTKType() { return d_VTKType; };

        void ChangeOrientation();

        unsigned short GetNumNeighborElement() { return d_numNeighborElement; };
        unsigned short GetNumNeighborNode(unsigned short val_node) { return d_numNeighborNode[val_node]; };

        unsigned short GetFace(unsigned short val_face, unsigned short val_index) { return d_face[val_face][val_index]; };
        unsigned short GetNeighborNode(unsigned short val_node, unsigned short val_index) { return d_neighborNode[val_node][val_index]; };
       
        unsigned short GetNumNodeFace(unsigned short val_face) { return d_numNodeFace[val_face]; };

        unsigned short GetMaxNodeFace() { return d_maxNodeFace; };

        unsigned long GetDomainElement() { return d_domainElement; };
        void SetDomainElement(unsigned long val_domainElement) { d_domainElement = val_domainElement; };
      
    private:
        unsigned short d_face[1][2];			/*!< \brief Matrix to store the local nodes of all the faces. */
        unsigned short d_neighborNode[2][1];    /*!< \brief Neighbor to a nodes in the element. */
        unsigned short d_numNodeFace[1];		/*!< \brief Number of nodes of each face of the element. */
            
        unsigned short d_numFace;				/*!< \brief Number of faces of the element. */
        unsigned short d_numNode;				/*!< \brief Number of nodes of the element. */
        unsigned short d_VTKType;				/*!< \brief Type of element using VTK nomenclature. */
        unsigned short d_maxNodeFace;			/*!< \brief Maximum number of nodes for a face. */
        unsigned short d_numNeighborNode[2];	/*!< \brief Number of Neighbor to a nodes in the element. */
        unsigned short d_numNeighborElement;	/*!< \brief Number of Neighbor_Elements. */
    };
}

#endif
