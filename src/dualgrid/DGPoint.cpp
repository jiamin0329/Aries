/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    25-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "DGPoint.hpp"

namespace ARIES
{
    DGPoint::DGPoint(unsigned short val_nDim, unsigned long val_globalIndex, IProcData* procData) : DualGrid(val_nDim) 
    {
        d_numElem = 0;
        d_numPoint = 0;
         
        d_elem.clear(); 
        d_point.clear(); 
        d_edge.clear();

        /*
         *  Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume
         */
        if (procData->GetUnsteadyType() == STEADY)
            d_volume.assign(1, 0.0);
        else
            d_volume.assign(3, 0.0);

        d_coord.assign(d_nDim, 0.0);
        
        /*--- Indicator if the control volume has been agglomerated ---*/
        d_isAgglomerate = false;

        /*--- Flip the normal orientation ---*/
        d_isFlipOri = false;

        /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
        d_isMove = true;

        /*--- Identify boundaries, physical boundaries (not send-receive
          condition), detect if an element belong to the domain or it must
          be computed with other processor  ---*/
        d_isDomain = true;
        d_isBoundary = false;
        d_isPhyBoundary = false;
        d_isSolidBoundary = false;

        /*--- Set the global index in the parallel simulation ---*/
        d_globalIndex = val_globalIndex;

        /*--- Set the color for mesh partitioning ---*/
        d_color = 0;

        /*--- For smoothing the numerical grid coordinates ---*/
        if (procData->GetSmoothNumGrid()) 
        {
           d_coord_old.assign(d_nDim, 0.0);
           d_coord_sum.assign(d_nDim, 0.0);
        }

        /*--- Storage of grid velocities for dynamic meshes ---*/
        if (procData->GetGridMovement()) 
        {
            d_gridVel.assign(d_nDim, 0.0);
            
            /*--- Gradient of the grid velocity ---*/
            //d_gridVelGrad.assign(d_nDim*d_nDim, 0.0);
            d_gridVelGrad.clear();
            d_gridVelGrad.resize(d_nDim);
            for (unsigned short idim = 0; idim < d_nDim; idim++)
                d_gridVelGrad[idim].assign(d_nDim, 0.0);
            
 
            /*--- Structures for storing old node coordinates for computing grid
              velocities via finite differencing with dynamically deforming meshes. ---*/
            if (procData->GetUnsteadyType() != STEADY) 
            {
                d_coord_p1.assign(d_nDim, 0.0);
                d_coord_n.assign(d_nDim, 0.0);
                d_coord_n1.assign(d_nDim, 0.0);;
            }
        }

        /*--- Intialize the value of the curvature ---*/
        d_curvature = 0.0;

    }

    DGPoint::DGPoint(double val_coord_0, double val_coord_1, unsigned long val_globalIndex, IProcData* procData): DualGrid(2) 
    {
        d_numElem = 0;
        d_numPoint = 0;
         
        d_elem.clear(); 
        d_point.clear(); 
        d_edge.clear();

        /*
         *  Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume
         */
        if (procData->GetUnsteadyType() == STEADY)
            d_volume.assign(1, 0.0);
        else
            d_volume.assign(3, 0.0);
        
        d_coord.assign(2, 0.0);
        d_coord[0] = val_coord_0;
        d_coord[1] = val_coord_1;

        /*--- Indicator if the control volume has been agglomerated ---*/
        d_isAgglomerate = false;

        /*--- Flip the normal orientation ---*/
        d_isFlipOri = false;

        /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
        d_isMove = true;

        /*--- Identify boundaries, physical boundaries (not send-receive
          condition), detect if an element belong to the domain or it must
          be computed with other processor  ---*/
        d_isDomain = true;
        d_isBoundary = false;
        d_isPhyBoundary = false;
        d_isSolidBoundary = false;

        /*--- Set the color for mesh partitioning ---*/
        d_color = 0;

        /*--- Set the global index in the parallel simulation ---*/
        d_globalIndex = val_globalIndex;

        /*--- For smoothing the numerical grid coordinates ---*/
        if (procData->GetSmoothNumGrid()) 
        {
           d_coord_old.assign(d_nDim, 0.0);
           d_coord_sum.assign(d_nDim, 0.0);
        }

        /*--- Storage of grid velocities for dynamic meshes ---*/
        if (procData->GetGridMovement()) 
        {
            d_gridVel.assign(d_nDim, 0.0);
            
            /*--- Gradient of the grid velocity ---*/
            //d_gridVelGrad.assign(d_nDim*d_nDim, 0.0);
            d_gridVelGrad.clear();
            d_gridVelGrad.resize(d_nDim);
            for (unsigned short idim = 0; idim < d_nDim; idim++)
                d_gridVelGrad[idim].assign(d_nDim, 0.0);
            
 
            /*--- Structures for storing old node coordinates for computing grid
              velocities via finite differencing with dynamically deforming meshes. ---*/
            if (procData->GetUnsteadyType() != STEADY) 
            {
                d_coord_p1.assign(d_nDim, 0.0);
                d_coord_n.assign(d_nDim, 0.0);
                d_coord_n1.assign(d_nDim, 0.0);;
            }
        }

        /*--- Intialize the value of the curvature ---*/
        d_curvature = 0.0;
    }

    DGPoint::DGPoint(double val_coord_0, double val_coord_1, double val_coord_2, unsigned long val_globalIndex, IProcData* procData): DualGrid(3) 
    {
        d_numElem = 0;
        d_numPoint = 0;
         
        d_elem.clear(); 
        d_point.clear(); 
        d_edge.clear();

        /*
         *  Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume
         */
        if (procData->GetUnsteadyType() == STEADY)
            d_volume.assign(1, 0.0);
        else
            d_volume.assign(3, 0.0);
        
        d_coord.assign(3, 0.0);
        d_coord[0] = val_coord_0;
        d_coord[1] = val_coord_1;
        d_coord[2] = val_coord_2;

        /*--- Indicator if the control volume has been agglomerated ---*/
        d_isAgglomerate = false;

        /*--- Flip the normal orientation ---*/
        d_isFlipOri = false;

        /*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
        d_isMove = true;

        /*--- Identify boundaries, physical boundaries (not send-receive
          condition), detect if an element belong to the domain or it must
          be computed with other processor  ---*/
        d_isDomain = true;
        d_isBoundary = false;
        d_isPhyBoundary = false;
        d_isSolidBoundary = false;

        /*--- Set the color for mesh partitioning ---*/
        d_color = 0;

        /*--- Set the global index in the parallel simulation ---*/
        d_globalIndex = val_globalIndex;

        /*--- For smoothing the numerical grid coordinates ---*/
        if (procData->GetSmoothNumGrid()) 
        {
            d_coord_old.assign(d_nDim, 0.0);
            d_coord_sum.assign(d_nDim, 0.0);
        }

        /*--- Storage of grid velocities for dynamic meshes ---*/
        if (procData->GetGridMovement()) 
        {
            d_gridVel.assign(d_nDim, 0.0);
            
            /*--- Gradient of the grid velocity ---*/
            //d_gridVelGrad.assign(d_nDim*d_nDim, 0.0);
            d_gridVelGrad.clear();
            d_gridVelGrad.resize(d_nDim);
            for (unsigned short idim = 0; idim < d_nDim; idim++)
                d_gridVelGrad[idim].assign(d_nDim, 0.0);
            
 
            /*--- Structures for storing old node coordinates for computing grid
              velocities via finite differencing with dynamically deforming meshes. ---*/
            if (procData->GetUnsteadyType() != STEADY) 
            {
                d_coord_p1.assign(d_nDim, 0.0);
                d_coord_n.assign(d_nDim, 0.0);
                d_coord_n1.assign(d_nDim, 0.0);;
            }
        }

        /*--- Intialize the value of the curvature ---*/
        d_curvature = 0.0;
    }

    DGPoint::~DGPoint() 
    {
        d_elem.~vector();
        d_point.~vector();
        d_edge.~vector();
        
        d_childrenCV.~vector();

        d_volume.~vector();
        d_vertex.~vector();
        
        d_coord.~vector();
        d_coord_old.~vector();
        d_coord_sum.~vector();
        d_coord_n.~vector();
        d_coord_n1.~vector();
        d_coord_p1.~vector();
                                         
        d_gridVel.~vector();
        d_gridVelGrad.~vector();
    }

    void DGPoint::SetCoord(double *val_coord)
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord[iDim] = val_coord[iDim];
    }

    void DGPoint::SetElem(unsigned long val_elem) 
    { 
        d_elem.push_back(val_elem); 
        d_numElem = d_elem.size(); 
    }

    void DGPoint::ResetElem()
    {
        d_elem.clear();
        d_numElem = 0;
    }
    
    void DGPoint::SetPoint(unsigned long val_point) 
    {
        unsigned short iPoint;
        bool new_point;

        /*--- Look for the point in the list ---*/
        new_point = true;
        for (iPoint = 0; iPoint < d_numPoint; iPoint++)
            if (d_point[iPoint] == val_point) 
            {
                new_point = false;
                break;
            }

        /*--- Store the point structure and dimensionalizate edge structure ---*/
        if (new_point) 
        {
            d_point.push_back(val_point);
            d_edge.push_back(-1);
            d_numPoint = d_point.size();
        }
    }

    void DGPoint::ResetPoint()
    {
        d_point.clear();
        d_edge.clear();
        d_numPoint = 0;
    }

    long DGPoint::GetVertex(unsigned short val_iVertex)
    {
        if (d_isBoundary)
            return d_vertex[val_iVertex];
        else
            return -1;
    }

    void DGPoint::SetVertex(long val_vertex, unsigned short val_iVertex)
    {
        if (d_isBoundary)
            d_vertex[val_iVertex] = val_vertex;
    }
    
    void DGPoint::SetBoundary(unsigned short val_numMarker) 
    {
        if (!d_isBoundary) 
        {
            d_vertex.assign (val_numMarker, 0.0);
            for (unsigned short iMarker = 0; iMarker < val_numMarker; iMarker++)
                d_vertex[iMarker] = -1;
        }
        d_isBoundary = true;
    }

    void DGPoint::ResetBoundary()
    {
        d_vertex.clear();
        d_isBoundary = false;
    }
    
    void DGPoint::SetChildrenCV(unsigned short val_iChildrenCV, unsigned long val_childrenCV)
    {
        if (d_childrenCV.size() <= val_iChildrenCV)
            d_childrenCV.resize(val_iChildrenCV + 1);
        d_childrenCV[val_iChildrenCV] = val_childrenCV;
    }

    void DGPoint::SetCoord_n()
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord_n[iDim] = d_coord[iDim];
    }

    void DGPoint::SetCoord_n1()
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord_n1[iDim] = d_coord_n[iDim];
    }

    void DGPoint::SetCoord_p1(double *val_coord)
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord_p1[iDim] = val_coord[iDim];
    }
    
    void DGPoint::AddCoordSum(double *val_coord_sum)
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord_sum[iDim] += val_coord_sum[iDim];
    }

    void DGPoint::SetCoordSumZero()
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord_sum[iDim] = 0.0;
    }
    
    void DGPoint::SetCoordOld(double *val_coord_old)
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_coord_old[iDim] = val_coord_old[iDim];
    }

    void DGPoint::SetGridVel(double *val_gridVel)
    {
        for (unsigned short iDim = 0; iDim < d_nDim; iDim++)
            d_gridVel[iDim] = val_gridVel[iDim];
    }
       

}

















