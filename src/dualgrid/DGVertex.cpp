/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Vertex class definition for dual grid
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "DGVertex.hpp"

namespace ARIES
{
    DGVertex::DGVertex(unsigned long val_point, unsigned short val_nDim) : DualGrid(val_nDim)
    {
        d_node = val_point;
        
        d_normal[0] = 0.0;
        d_normal[1] = 0.0;
        d_normal[2] = 0.0;

        d_varCoord[0] = 0.0;
        d_varCoord[1] = 0.0;
        d_varCoord[2] = 0.0;
    }

    DGVertex::~DGVertex() 
    {
    }

    void DGVertex::SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG) 
    {
        double vec_a[3] = { 0.0, 0.0, 0.0 };
        double vec_b[3] = { 0.0, 0.0, 0.0 };
        double dimNormal[3] = { 0.0, 0.0, 0.0 };
        unsigned short iDim;

        for (iDim = 0; iDim < d_nDim; iDim++) 
        {
            vec_a[iDim] =     val_coord_Elem_CG[iDim] - val_coord_Edge_CG[iDim];
            vec_b[iDim] = val_coord_FaceElem_CG[iDim] - val_coord_Edge_CG[iDim];
        }

        dimNormal[0] =  0.5*(vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]);
        dimNormal[1] = -0.5*(vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
        dimNormal[2] =  0.5*(vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

        d_normal[0] += dimNormal[0];
        d_normal[1] += dimNormal[1];
        d_normal[2] += dimNormal[2];
    }

    void DGVertex::SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_Elem_CG) 
    {
        double dimNormal[2];

        dimNormal[0] =   val_coord_Elem_CG[1] - val_coord_Edge_CG[1];
        dimNormal[1] = -(val_coord_Elem_CG[0] - val_coord_Edge_CG[0]);

        d_normal[0] += dimNormal[0];
        d_normal[1] += dimNormal[1];
    }

    void DGVertex::SetNormal(double *val_face_normal)
    {
        d_normal[0] = val_face_normal[0];
        d_normal[1] = val_face_normal[1];
        if (d_nDim == 3)
            d_normal[2] = val_face_normal[2];
    }
    
    void DGVertex::AddNormal(double *val_face_normal) 
    {
        d_normal[0] += val_face_normal[0];
        d_normal[1] += val_face_normal[1];
        if (d_nDim == 3)
            d_normal[2] += val_face_normal[2];
    }

    void DGVertex::SetNormalZero()
    {
        d_normal[0] = 0.0;
        d_normal[1] = 0.0;
        if (d_nDim == 3)
            d_normal[2] = 0.0;
    }

    void DGVertex::SetVarCoord(double *val_varCoord)
    {
        d_varCoord[0] = val_varCoord[0];
        d_varCoord[1] = val_varCoord[1];
        if (d_nDim == 3)
            d_varCoord[2] = val_varCoord[2];
    }

    void DGVertex::AddVarCoord(double *val_varCoord)
    {
        d_varCoord[0] += val_varCoord[0];
        d_varCoord[1] += val_varCoord[1];
        if (d_nDim == 3)
            d_varCoord[2] += val_varCoord[2];
    }

    void DGVertex::SetCoord(double *val_coord)
    {
        d_cartCoord[0] = val_coord[0];
        d_cartCoord[1] = val_coord[1];
        if (d_nDim == 3)
            d_cartCoord[2] = val_coord[2];
    }

    void DGVertex::AddCoord(double *val_coord)
    {
        d_cartCoord[0] += val_coord[0];
        d_cartCoord[1] += val_coord[1];
        if (d_nDim == 3)
            d_cartCoord[2] += val_coord[2];
    }
    
    void DGVertex::SetDonorPoint(long val_periodiGRID_DGPoint, long val_processor)
    {
        d_periodicGRID_DGPoint[0] = val_periodiGRID_DGPoint;
        d_periodicGRID_DGPoint[1] = val_processor;
    }
}













