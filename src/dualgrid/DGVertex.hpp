/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for vertex definition
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_DGVERTEX_HPP
#define ARIES_DGVERTEX_HPP

#include "DualGrid.hpp"

#include <vector>

using namespace std;

namespace ARIES
{
    class DGVertex: public DualGrid
    {
    public:
        DGVertex(unsigned long val_point, unsigned short val_nDim);
        ~DGVertex();

        unsigned short GetNumNode() { return 1; };
        unsigned long GetNode(unsigned short val_node) { return d_node; };

        void SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_Elem_CG);
        void SetNodeCoord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG);

        double GetNormal(unsigned short val_iDim) { return d_normal[val_iDim]; };
        void SetNormal(double *val_face_normal);
        void AddNormal(double *val_face_normal);
        void SetNormalZero();

        double GetAuxVar() { return d_auxVar; };
        void SetAuxVar(double val_auxVar) { d_auxVar = val_auxVar; };
        void AddAuxVar(double val_auxVar) { d_auxVar += val_auxVar; };
       
        double GetVarCoord(unsigned short val_iDim) { return d_varCoord[val_iDim]; };
        void SetVarCoord(double *val_varCoord);
        void AddVarCoord(double *val_varCoord);
        
        double GetCoord(unsigned short val_iDim) { return d_cartCoord[val_iDim]; };
        void SetCoord(double *val_coord);
        void AddCoord(double *val_coord);
        
        short GetRotationType() { return d_rotationType; };
        void SetRotationType(short val_rotationType) { d_rotationType = val_rotationType; };

        long GetDonorPoint() { return d_periodicGRID_DGPoint[0]; };
        long GetDonorProcessor() { return d_periodicGRID_DGPoint[1]; };
        long *GetPeriodiGRID_DGPointDomain(){ return d_periodicGRID_DGPoint; };
        void SetDonorPoint(long val_periodiGRID_DGPoint, long val_processor);

        long GetDonorElem() { return d_donorElem; };
        void SetDonorElem(long val_donorElem) { d_donorElem = val_donorElem; };
        
        double GetBasisFunction(unsigned short val_node) { return d_basisFunction[val_node]; }
        void SetBasisFunction(unsigned short val_node, double val_basis) { d_basisFunction[val_node] = val_basis; }

        unsigned long GetNormalNeighbor() { return d_normalNeighbor; };
        void SetNormalNeighbor(unsigned long val_normalNeighbor) { d_normalNeighbor = val_normalNeighbor; }
       
    private:
        unsigned long d_node;	                          /*!< \brief Vector to store the global nodes of an element. */
        double d_normal[3];			                      /*!< \brief Normal coordinates of the element and its center of gravity. */
        double d_auxVar;			                      /*!< \brief Auxiliar variable defined only on the surface. */
        double d_cartCoord[3];		                      /*!< \brief Vertex cartesians coordinates. */
        double d_varCoord[3];		                      /*!< \brief Used for storing the coordinate variation due to a surface modification. */
        long d_periodicGRID_DGPoint[2];			          /*!< \brief Store the periodic point of a boundary (iProcessor, iPoint) */
        short d_rotationType;			                  /*!< \brief Type of rotation associated with the vertex (MPI and periodic) */
        unsigned long d_normalNeighbor;                   /*!< \brief Index of the closest neighbor. */
        unsigned long d_donorElem;                        /*!< \brief Store the donor element for interpolation across zones. */
        double d_basisFunction[3];                        /*!< \brief Basis function values for interpolation across zones. */
    };
}

#endif





