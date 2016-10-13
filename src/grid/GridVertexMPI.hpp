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

#ifndef ARIES_GRIDVERTEXMPI_HPP
#define ARIES_GRIDVERTEXMPI_HPP

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include <iostream>
#include <vector>
#include <cstdlib>

#include "Grid.hpp"

namespace ARIES
{
    class GridVertexMPI : public Grid
    {
    public:
        GridVertexMPI(unsigned long val_point, unsigned short val_nDim);
        ~GridVertexMPI();


        unsigned long GetNode(unsigned short val_node) { return d_node[val_node]; };
        void SetNode(unsigned short val_node, unsigned long val_point) { d_node[val_node] = val_point; };

        unsigned short GetNumFace() { return d_numFace; };
        unsigned short GetNumNode() { return d_numNode; };

        unsigned short GetVTKType() { return d_VTKType; };

        unsigned short GetRotationType() { return d_rotationType; };
        void SetRotationType(unsigned short val_rotationType) { d_rotationType = val_rotationType; };

        void ChangeOrientation();

        unsigned short GetNumNeighborElement() { return 0; };
        unsigned short GetNumNeighborNode(unsigned short val_node) { return 0; };
        unsigned short GetNumNodeFace(unsigned short val_face) { return 0; };
        unsigned short GetMaxNodeFace() { return 0; };
        unsigned short GetFace(unsigned short val_face, unsigned short val_index) { return 0; };
        unsigned short GetNeighborNode (unsigned short val_node, unsigned short val_index) { return 0; };
        
    private:
        unsigned short d_numNode;				/*!< \brief Number of nodes of the element. */
        unsigned short d_numFace;				/*!< \brief Number of faces of the element. */
        unsigned short d_maxNodeFace;			/*!< \brief Maximum number of nodes for a face. */
        unsigned short d_numNeighborElement;	/*!< \brief Number of Neighbor_Elements. */
        unsigned short d_VTKType;				/*!< \brief Type of element using VTK nomenclature. */
        unsigned short d_rotationType;			/*!< \brief Definition of the rotation, traslation of the solution at the vertex. */
    };
    
}

#endif

