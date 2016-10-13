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


#ifndef ARIES_DGPOINT_HPP
#define ARIES_DGPOINT_HPP

#include "DualGrid.hpp"

#include "IProcData.hpp"

#include <vector>


using namespace std;


namespace ARIES
{
    class DGPoint: public DualGrid
    {
    public:
        DGPoint(unsigned short val_nDim, unsigned long val_globalIndex, IProcData* procData);
        DGPoint(double val_coord_0, double val_coord_1, unsigned long val_globalindex, IProcData* procData);
        DGPoint(double val_coord_0, double val_coord_1, double val_coord_2, unsigned long val_globalindex, IProcData* procData);
        ~DGPoint();

        unsigned short GetNumNode() { return 0; };

        void SetNodesCoord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG, double *val_coord_Elem_CG) {};
        void SetNodesCoord(double *val_coord_Edge_CG, double *val_coord_Elem_CG) {};
     
        double GetNormal(unsigned short val_iDim) { return 0.0; };
        void SetNormal(double *val_normal) {};
        void AddNormal(double *val_normal) {};
        void SetNormalZero() {};
        
        bool GetDomain() { return d_isDomain; };
        void SetDomain(bool val_isDomain) { d_isDomain = val_isDomain; };
        
        double GetWallDist() { return d_wallDist; };
        void SetWallDist(double val_wallDist) {d_wallDist = val_wallDist; };

        double GetSharpEdgeDist() { return d_sharpEdgeDist; };
        void SetSharpEdgeDist(double val_distance) { d_sharpEdgeDist = val_distance; };
        
        double GetCurvature() { return d_curvature; };
        void SetCurvature(double val_curvature) { d_curvature = val_curvature; };

        unsigned short GetNumPoint() { return d_numPoint; };
        void SetNumPoint(unsigned short val_nPoint) { d_numPoint = val_nPoint; };

        unsigned short GetNumElem() { return d_numElem; };
        void SetNumElem(unsigned short val_numElem) { d_numElem = val_numElem; };

        bool GetFlipOrientation() { return d_isFlipOri; };
        void SetFlipOrientation(bool val_flipOri) { d_isFlipOri = val_flipOri; };
        
        double GetCoord(unsigned short val_dim) { return d_coord[val_dim]; };
        void SetCoord(double *val_coord);
        void SetCoord(unsigned short val_dim, double val_coord) { d_coord[val_dim]  = val_coord; };
        void AddCoord(unsigned short val_dim, double val_coord) { d_coord[val_dim] += val_coord; };
        
        unsigned long GetElem(unsigned short val_elem) { return d_elem[val_elem]; };
        void SetElem(unsigned long val_elem);
        void ResetElem();

        unsigned long GetPoint(unsigned short val_point) { return d_point[val_point]; };
        void SetPoint(unsigned long val_point);
        void ResetPoint();

        long GetEdge(unsigned short val_iEdge) { return d_edge[val_iEdge]; };
        void SetEdge(long val_edge, unsigned short val_iEdge) { d_edge[val_iEdge] = val_edge; };
        
        long GetVertex(unsigned short val_iVertex);
        void SetVertex(long val_vertex, unsigned short val_iVertex);

        double GetVolume() { return d_volume[0]; };
        void SetVolume(double val_volume) { d_volume[0]  = val_volume; };
        void AddVolume(double val_volume) { d_volume[0] += val_volume; };
        
        double GetVolume_n() { return d_volume[1]; };
        void SetVolume_n() { d_volume[1] = d_volume[0]; };
       
        double GetVolume_nM1() { return d_volume[2]; };
        void SetVolume_nM1() { d_volume[2] = d_volume[1]; };
        
        bool GetMove() { return d_isMove; };
        void SetMove(bool val_isMove) { d_isMove = val_isMove; };
        
        bool GetBoundary() { return d_isBoundary; };
        void SetBoundary(bool val_isBoundary) { d_isBoundary = val_isBoundary; };
        void SetBoundary(unsigned short val_nmarker);
        void ResetBoundary();

        bool GetPhysicalBoundary() { return d_isPhyBoundary; };
        void SetPhysicalBoundary(bool val_isPhyBoundary) { d_isPhyBoundary = val_isPhyBoundary; };

        bool GetSolidBoundary() { return d_isSolidBoundary; };
        void SetSolidBoundary(bool val_isSolidBoundary) { d_isSolidBoundary = val_isSolidBoundary; };
        
        unsigned short GetColor() { return d_color; };
        void SetColor(unsigned short val_color) { d_color = val_color; };
        
        unsigned short GetNumNeighbor() { return d_numNeighbor; };
        void SetNumNeighbor(unsigned short val_numNeighbor) { d_numNeighbor = val_numNeighbor; };
        
        unsigned long GetGlobalIndex() { return d_globalIndex; };
        void SetGlobalIndex(unsigned long val_globalIndex) { d_globalIndex = val_globalIndex; };
     
        double GetCoord_n(unsigned short val_iDim) { return d_coord_n[val_iDim]; };
        void SetCoord_n();
        
        double GetCoord_n1(unsigned short val_iDim) { return d_coord_n1[val_iDim]; };
        void SetCoord_n1();
        
        double GetCoord_p1(unsigned short val_iDim) { return d_coord_p1[val_iDim]; };
        void SetCoord_p1(double *val_coord);

        unsigned long GetParentCV() { return d_parentCV; };
        void SetParentCV(unsigned long val_parentCV) { d_parentCV = val_parentCV; d_isAgglomerate = true; };

        unsigned long GetChildrenCV(unsigned short val_iChildrenCV) { return d_childrenCV[val_iChildrenCV]; };
        void SetChildrenCV(unsigned short val_iChildrenCV, unsigned long val_childrenCV);

        bool GetAgglomerate() { return d_isAgglomerate; };
        bool GetAgglomerateIndirect() { return d_isAgglomerateIndirect; };
        void SetAgglomerateIndirect(bool val_isAgglomerateIndirect) { d_isAgglomerateIndirect = val_isAgglomerateIndirect; };

        unsigned short GetNumChildrenCV() { return d_numChildrenCV; };
        void SetNumChildrenCV(unsigned short val_numChildrenCV) { d_numChildrenCV = val_numChildrenCV; };

        double GetCoordSum(unsigned short val_iDim) { return d_coord_sum[val_iDim]; };
        void AddCoordSum(double *val_coord_sum);
        void SetCoordSumZero();

        double GetCoordOld(unsigned short val_iDim) { return d_coord_old[val_iDim]; };
        void SetCoordOld(double *val_coord_old);

        double GetGridVel(unsigned short val_iDim) { return d_gridVel[val_iDim]; };
        void SetGridVel(unsigned short val_iDim, double val_gridVel) { d_gridVel[val_iDim] = val_gridVel;};
        void SetGridVel(double *val_gridvel);
        
        double GetGridVelGrad(unsigned short val_iDim, unsigned short val_jDim) { return d_gridVelGrad[val_iDim][val_jDim]; };
        void SetGridVelGrad(unsigned short val_iDim, unsigned short val_jDim, double val_value) { d_gridVelGrad[val_iDim][val_jDim] = val_value; };
        
    private:
        unsigned short d_numElem;	                    /*!< \brief Number of elements that set up the control volume. */
        unsigned short d_numPoint;                      /*!< \brief Number of points that set up the control volume  */
                                                        
        vector<unsigned long> d_elem;		            /*!< \brief Elements that set up a control volume around a node. */
        vector<unsigned long> d_point;	                /*!< \brief Points surrounding the central node of the control volume. */
        vector<unsigned long> d_edge;		            /*!< \brief Edges that set up a control volume. */
                                                        
        vector<double> d_volume;	                    /*!< \brief Volume or Area of the control volume in 3D and 2D. */
                                                        
        bool d_isDomain;		                        /*!< \brief To see if a point must be computed or belong to another boundary */
        bool d_isBoundary;                              /*!< \brief To see if a point belong to the boundary (including MPI). */
        bool d_isPhyBoundary;			                /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
        bool d_isSolidBoundary;			                /*!< \brief To see if a point belong to the physical boundary (without includin MPI). */
                                                        
        vector<unsigned long> d_vertex;                 /*!< \brief Index of the vertex that correspond which the control volume (we need one for each marker in the same node). */
        
        vector<double> d_coord;	                        /*!< \brief vector with the coordinates of the node. */
        vector<double> d_coord_old;		                /*!< \brief Old coordinates vector for geometry smoothing. */
        vector<double> d_coord_sum;		                /*!< \brief Sum of coordinates vector for geometry smoothing. */
        vector<double> d_coord_n;		                /*!< \brief Coordinates at time n for use with dynamic meshes. */
        vector<double> d_coord_n1;		                /*!< \brief Coordinates at time n-1 for use with dynamic meshes. */
        vector<double> d_coord_p1;		                /*!< \brief Coordinates at time n+1 for use with dynamic meshes. */
                                                        
        vector<double> d_gridVel;	                    /*!< \brief Velocity of the grid for dynamic mesh cases. */
        vector<vector<double> > d_gridVelGrad;          /*!< \brief Gradient of the grid velocity for dynamic meshes. */
        
        unsigned long d_parentCV;			            /*!< \brief Index of the parent control volume in the agglomeration process. */
        unsigned short d_numChildrenCV;		            /*!< \brief Number of children in the agglomeration process. */
        vector<unsigned long> d_childrenCV;		        /*!< \brief Index of the children control volumes in the agglomeration process. */
        
        bool d_isAgglomerateIndirect;					/*!< \brief This flag indicates if the indirect points can be agglomerated. */
        bool d_isAgglomerate;					        /*!< \brief This flag indicates if the element has been agglomerated. */
        bool d_isMove;					                /*!< \brief This flag indicates if the point is going to move in the grid deformation process. */
        
        unsigned short d_color;	                        /*!< \brief Color of the point in the partitioning strategy. */
        
        double d_wallDist;	                            /*!< \brief Distance to the nearest wall. */
        double d_sharpEdgeDist;	                        /*!< \brief Distance to a sharp edge. */
        double d_curvature;	                            /*!< \brief Value of the surface curvature. */
        unsigned long d_globalIndex;	                /*!< \brief Global index in the parallel simulation. */
        unsigned short d_numNeighbor;	                /*!< \brief . */
        bool d_isFlipOri;	                            /*!< \brief Flip the orientation of the normal. */
    };
}

#endif
