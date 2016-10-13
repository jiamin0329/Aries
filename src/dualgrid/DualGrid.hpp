/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for controlling the dual volume definition
 *    The dual volume is compose by three main elements: points, edges, and vertices.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    24-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_DUALGRID_HPP
#define ARIES_DUALGRID_HPP

#ifdef ARIES_HAVE_MPI
#include "mpi.h"
#endif

namespace ARIES
{
    class DualGrid
    {
    public:
        DualGrid(unsigned short val_nDim);
        ~DualGrid();

        virtual unsigned short GetNumNode() = 0;
        virtual unsigned long GetNode(unsigned short val_node) = 0;
        
        virtual double GetCoord(unsigned short val_iDim) = 0;
        virtual void SetCoord(double *val_coord) = 0;

        virtual void SetNodesCoord(double *val_coord_Edge_CG, double *val_coord_Elem_CG) = 0;
        virtual void SetNodesCoord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG) = 0;

        virtual double GetNormal(unsigned short val_iDim) = 0;
        virtual void SetNormal(double *val_normal) = 0;
        virtual void AddNormal(double *val_normal) = 0;
        virtual void SetNormalZero() = 0;

    protected:
        static unsigned short d_nDim; /*!< \brief Number of dimensions of the problem. */
    };
}

#endif
