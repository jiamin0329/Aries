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

#ifndef ARIES_DGEDGE_HPP
#define ARIES_DGEDGE_HPP

#include "DualGrid.hpp"

#include <vector>

using namespace std;

namespace ARIES
{
    class DGEdge: public DualGrid
    {
    public:
        DGEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim);
        ~DGEdge(void);

        unsigned short GetNumNode() { return 2; };
        unsigned long GetNode(unsigned short val_node) { return d_node[val_node]; };

        double GetCoord(unsigned short val_iDim) { return 0.0; };
        void SetCoord(double *val_coord) {};
        
        void SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_Elem_CG);
        void SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG);
        
        double GetNormal(unsigned short val_iDim) { return d_normal[val_iDim]; };
        void SetNormal(double *val_normal);
        void AddNormal(double *val_normal);
        void SetNormalZero();

        double GetVolume(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG, double *val_coord_Point);
        double GetVolume(double *val_coord_Edge_CG, double *val_coord_Elem_CG, double *val_coord_Point);

        double GetCG(unsigned short val_iDim) { return d_coordCG[val_iDim]; };
        void SetCG(double **val_coord);
        
    private:
        double d_coordCG[3];			    /*!< \brief Center-of-gravity of the element. */
        double d_normal[3];				    /*!< \brief Normal al elemento y coordenadas de su centro de gravedad. */
        unsigned long d_node[2];		    /*!< \brief Vector to store the global nodes of an element. */
    };
}

#endif
















