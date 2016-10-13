/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for edge definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "DGEdge.hpp"

#include <cmath>

namespace ARIES
{
    DGEdge::DGEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim) : DualGrid(val_nDim) 
    {
        d_node[0] = val_iPoint;
        d_node[1] = val_jPoint;
        
        d_normal[0] = 0.0;
        d_normal[1] = 0.0;
        d_normal[2] = 0.0;
        
        d_coordCG[0] = 0.0;
        d_coordCG[1] = 0.0;
        d_coordCG[2] = 0.0;
    }

    DGEdge::~DGEdge() 
    {
    }

    void DGEdge::SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_Elem_CG)
    {
        double dimNormal[2];

        dimNormal[0] =   val_coord_Elem_CG[1] - val_coord_Edge_CG[1];
        dimNormal[1] = -(val_coord_Elem_CG[0] - val_coord_Edge_CG[0]);

        d_normal[0] += dimNormal[0];
        d_normal[1] += dimNormal[1];
    }

    void DGEdge::SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG) 
    {
        unsigned short iDim;
        double vec_a[3] = { 0.0, 0.0, 0.0 };
        double vec_b[3] = { 0.0, 0.0, 0.0 };
        double dimNormal[3];

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

    void DGEdge::SetNormal(double *val_normal)
    {
        d_normal[0] = val_normal[0];
        d_normal[1] = val_normal[1];
        if (d_nDim == 3)
            d_normal[2] = val_normal[2];
    }

    void DGEdge::AddNormal(double *val_normal)
    {
        d_normal[0] += val_normal[0];
        d_normal[1] += val_normal[1];
        if (d_nDim == 3)
            d_normal[2] += val_normal[2];
    }

    void DGEdge::SetNormalZero()
    {   
        d_normal[0] = 0.0;
        d_normal[1] = 0.0;
        d_normal[2] = 0.0;
    }

    double DGEdge::GetVolume(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG, double *val_coord_Point) 
    {
        unsigned short iDim;
        double vec_a[3] = { 0.0, 0.0, 0.0 };
        double vec_b[3] = { 0.0, 0.0, 0.0 };
        double vec_c[3] = { 0.0, 0.0, 0.0 };
        double vec_d[3] = { 0.0, 0.0, 0.0 };
        double localVolume;

        for (iDim = 0; iDim < d_nDim; iDim++) 
        {
            vec_a[iDim] =     val_coord_Edge_CG[iDim] - val_coord_Point[iDim];
            vec_b[iDim] = val_coord_FaceElem_CG[iDim] - val_coord_Point[iDim];
            vec_c[iDim] =     val_coord_Elem_CG[iDim] - val_coord_Point[iDim];
        }

        vec_d[0] =   vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
        vec_d[1] = -(vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
        vec_d[2] =   vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

        localVolume = fabs(vec_c[0] * vec_d[0] + vec_c[1] * vec_d[1] + vec_c[2] * vec_d[2]) / 6.0;

        return localVolume;
    }

    double DGEdge::GetVolume(double *val_coord_Edge_CG, double *val_coord_Elem_CG, double *val_coord_Point) 
    {
        unsigned short iDim;
        double vec_a[2] = { 0.0, 0.0 };
        double vec_b[2] = { 0.0, 0.0 };
        double localVolume;

        for (iDim = 0; iDim < d_nDim; iDim++) 
        {
            vec_a[iDim] = val_coord_Elem_CG[iDim] - val_coord_Point[iDim];
            vec_b[iDim] = val_coord_Edge_CG[iDim] - val_coord_Point[iDim];
        }

        localVolume = 0.5*fabs(vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

        return localVolume;
    }

    void DGEdge::SetCG(double **val_coord) 
    {
        unsigned short iDim, iNode;

        for (iDim = 0; iDim < d_nDim; iDim++) 
        {
            for (iNode = 0; iNode < 2; iNode++)
                d_coordCG[iDim] += val_coord[iNode][iDim] / 2.0;
        }
    }
}










