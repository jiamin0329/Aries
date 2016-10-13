/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for mesh manager.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_MESHDATA_HPP
#define ARIES_MESHDATA_HPP

#include "const_def.h"

#include "IProcData.hpp"
#include "DGPoint.hpp"
#include "DGVertex.hpp"
#include "DGEdge.hpp"

#include <vector>
#include <string>

using namespace std;

namespace ARIES
{

    class MeshData
    {
    public:
        MeshData();
        virtual ~MeshData();

        unsigned short GetNumDim() { return d_numDim; };
        void SetNumDim(unsigned short val_numDim) { d_numDim = val_numDim; };

        unsigned short GetNumZone() { return d_numZone; };
        void SetNumZone(unsigned short val_numZone) { d_numZone = val_numZone; };
        
        unsigned short GetNumMarker() { return d_numMarker; };
        void SetNumMarker(unsigned short val_numMarker) { d_numMarker = val_numMarker; };

        unsigned long GetMaxGlobalPoint() { return d_maxGlobalPoint; };
        void SetMaxGlobalPoint(unsigned long val_maxGlobalPoint) { d_maxGlobalPoint = val_maxGlobalPoint; }
        
        unsigned long GetNumPoint() { return d_numPoint; };
        void SetNumPoint(unsigned long val_numPoint) { d_numPoint = val_numPoint; };
        
        unsigned long GetNumPointDomain() { return d_numPointDomain; };
        void SetNumPointDomain(unsigned long val_numPointDomain) { d_numPointDomain = val_numPointDomain; };

        unsigned long GetNumPointGhost() { return d_numPointGhost; };
        void SetNumPointGhost(unsigned long val_numPointGhost) { d_numPointGhost = val_numPointGhost; };
        
        unsigned long GetNumPointGlobal() { return d_numPointGlobal; };
        void SetNumPointGlobal(unsigned long val_numPointGlobal) { d_numPointGlobal = val_numPointGlobal; };
        
        unsigned long GetNumPointDomainGlobal() { return d_numPointDomainGlobal; };
        void SetNumPointDomainGlobal(unsigned long val_numPointDomainGlobal) { d_numPointDomainGlobal = val_numPointDomainGlobal; };

        unsigned long GetNumElem() { return d_numElem; };
        void SetNumElem(unsigned long val_numElem) { d_numElem = val_numElem; };
        
        unsigned long GetNumElemGlobal() { return d_numElemGlobal; };
        void SetNumElemGlobal(unsigned long val_numElemGlobal) { d_numElemGlobal = val_numElemGlobal; };
        
        unsigned long GetNumEdge() { return d_numEdge; };
        void SetNumEdge();
        
        unsigned long GetNumFace() { return d_numFace; };
        void SetNumFace();

        unsigned long GetNumVertexBound(unsigned short val_iMarker) { return d_numVertexBound[val_iMarker]; };
        void SetNumVertexBound(unsigned short val_iMarker, unsigned long val_numVertexBound) { d_numVertexBound[val_iMarker] = val_numVertexBound; };
        
        string GetMarkerTag(unsigned short val_iMarker) { return d_tagToMarker[val_iMarker]; };
        void SetMarkerTag(unsigned short val_iMarker, string val_tag) { d_tagToMarker[val_iMarker] = val_tag; };
        

        virtual unsigned long GetNumElemLine() { return d_numElemLine; };
        virtual unsigned long GetNumElemTria() { return d_numElemTria; };
        virtual unsigned long GetNumElemQuad() { return d_numElemQuad; };
        virtual unsigned long GetNumElemTetr() { return d_numElemTetr; };
        virtual unsigned long GetNumElemHexa() { return d_numElemHexa; };
        virtual unsigned long GetNumElemPris() { return d_numElemPris; };
        virtual unsigned long GetNumElemPyra() { return d_numElemPyra; };

        virtual unsigned long GetNumElemLineGlobal() { return d_numElemLineGlobal; };
        virtual unsigned long GetNumElemTriaGlobal() { return d_numElemTriaGlobal; };
        virtual unsigned long GetNumElemQuadGlobal() { return d_numElemQuadGlobal; };
        virtual unsigned long GetNumElemTetrGlobal() { return d_numElemTetrGlobal; };
        virtual unsigned long GetNumElemHexaGlobal() { return d_numElemHexaGlobal; };
        virtual unsigned long GetNumElemPrisGlobal() { return d_numElemPrisGlobal; };
        virtual unsigned long GetNumElemPyraGlobal() { return d_numElemPyraGlobal; };
        
        unsigned long GetNumElemBound(unsigned short val_marker) { return d_numElemBound[val_marker]; };
        void SetNumElemBound(unsigned short val_marker, unsigned long val_numElemBound);
        
        virtual vector<double> GetGeometryPlanes() { return d_xCoordList; };
        virtual void SetGeometryPlanes(IProcData* procData);

        virtual vector<vector<double> > GetXCoord() { return d_XcoordPlane; };
        virtual vector<vector<double> > GetYCoord() { return d_YcoordPlane; };
        virtual vector<vector<double> > GetZCoord() { return d_ZcoordPlane; };

        virtual vector<vector<unsigned long> > GetPlanePoints() { return d_planePoint; };

    private:
        unsigned short d_numDim;	                                    /*!< \brief Number of dimension of the problem. */
        unsigned short d_numZone;			                            /*!< \brief Number of zones in the problem. */
        unsigned short d_numMarker;	                                    /*!< \brief Number of different markers of the mesh. */
        unsigned long d_maxGlobalPoint;                                 /*!< \brief Greater global point in the domain local structure. */
                                                                        
        unsigned long d_numPoint;                                       /*!< \brief Number of points of the mesh. */
        unsigned long d_numPointDomain;		                            /*!< \brief Number of real points of the mesh. */
        unsigned long d_numPointGhost;		                            /*!< \brief Number of ghost points of the mesh. */
                                                                        
        unsigned long d_numPointGlobal;                                 /*!< \brief Total number of nodes in a simulation across all processors (including halos). */
        unsigned long d_numPointDomainGlobal;                           /*!< \brief Total number of nodes in a simulation across all processors (excluding halos). */
                                                                        
        unsigned long d_numElem;			                            /*!< \brief Number of elements of the mesh. */
        unsigned long d_numElemGlobal;	                                /*!< \brief Total number of elements in a simulation across all processors (all types). */
                                                                        
        unsigned long d_numEdge;			                            /*!< \brief Number of edges of the mesh. */
        unsigned long d_numFace;			                            /*!< \brief Number of faces of the mesh. */
                                                                        
        unsigned long d_numElemLine;                                    /*!< \brief Number of edges in the mesh. */
        unsigned long d_numElemLineGlobal;                              /*!< \brief Total number of edges in the mesh across all processors. */
        unsigned long d_numElemTria;                                    /*!< \brief Number of triangles in the mesh. */
        unsigned long d_numElemTriaGlobal;                              /*!< \brief Total number of triangles in the mesh across all processors. */
        unsigned long d_numElemQuad;                                    /*!< \brief Number of quadrangles in the mesh. */
        unsigned long d_numElemQuadGlobal;                              /*!< \brief Total number of quadrangles in the mesh across all processors. */
        unsigned long d_numElemTetr;                                    /*!< \brief Number of tetrahedra in the mesh. */
        unsigned long d_numElemTetrGlobal;                              /*!< \brief Total number of tetrahedra in the mesh across all processors. */
        unsigned long d_numElemHexa;                                    /*!< \brief Number of hexahedra in the mesh. */
        unsigned long d_numElemHexaGlobal;                              /*!< \brief Total number of hexahedra in the mesh across all processors. */
        unsigned long d_numElemPris;                                    /*!< \brief Number of prisms in the mesh. */
        unsigned long d_numElemPrisGlobal;                              /*!< \brief Total number of prisms in the mesh across all processors. */
        unsigned long d_numElemPyra;                                    /*!< \brief Number of pyramids in the mesh. */
        unsigned long d_numElemPyraGlobal;                              /*!< \brief Total number of pyramids in the mesh across all processors. */
                                                                        
        unsigned long d_numElemEdgeBound;                               /*!< \brief Number of edges on the mesh boundaries. */
        unsigned long d_numElemEdgeGlobalBound;                         /*!< \brief Total number of edges on the mesh boundaries across all processors. */
        unsigned long d_numElemTriaBound;                               /*!< \brief Number of triangles on the mesh boundaries. */
        unsigned long d_numElemTriaGlobalBound;                         /*!< \brief Total number of triangles on the mesh boundaries across all processors. */
        unsigned long d_numElemQuadBound;                               /*!< \brief Number of quads on the mesh boundaries. */
        unsigned long d_numElemQuadGlobalBound;                         /*!< \brief Total number of quads on the mesh boundaries across all processors. */
                                                                        
        vector<unsigned long> d_numElemBound;		                    /*!< \brief Number of elements of the boundary. */
                                                                        
        vector<string> d_tagToMarker;	                                /*!< \brief If you know the index of the boundary (depend of the
                                                                          grid definition), it gives you the maker (where the boundary
                                                                          is stored from 0 to boundaries). */
                                                                        
        //vector<Grid > d_elem;	                                        /*!< \brief Element std::vector (primal grid information). */
        //vector<Grid > d_face;			                                /*!< \brief Face std::vector (primal grid information). */
        //vector<Grid > d_bound;	                                        /*!< \brief Boundary std::vector (primal grid information). */
                                                                        
        vector<DGPoint  > d_node;			                            /*!< \brief Node std::vector (dual grid information). */
        vector<DGEdge   > d_edge;			                            /*!< \brief Edge std::vector (dual grid information). */
        vector<DGVertex > d_vertex;		                                /*!< \brief Boundary Vertex std::vector (dual grid information). */
                                                                        
        vector<unsigned long > d_numVertexBound;	                            /*!< \brief Number of vertex for each marker. */
        unsigned short d_numCommLevel;		                            /*!< \brief Number of non-blocking communication levels. */
        
        vector<unsigned long> d_periodicPoint[MAX_NUMBER_PERIODIC][2];  /*!< \brief PeriodiGRID::GRID_DGPoint[Periodic bc] and return the point that must be sent [0], and the image point in the periodic bc[1]. */
        vector<unsigned long> d_periodicElem[MAX_NUMBER_PERIODIC];      /*!< \brief PeriodicElem[Periodic bc] and return the elements that must be sent. */

        vector<short> d_markerAllSendRecv;

        /*
         *  Create std::vectors and distribute the values among the different planes queues
         */
        vector<vector<double> > d_XcoordPlane; /*!< \brief std::vector containing x coordinates of new points appearing on a single plane */
        vector<vector<double> > d_YcoordPlane; /*!< \brief std::vector containing y coordinates of  new points appearing on a single plane */
        vector<vector<double> > d_ZcoordPlane; 	/*!< \brief std::vector containing z coordinates of  new points appearing on a single plane */
        vector<vector<double> > d_faceAreaPlane; /*!< \brief std::vector containing area/volume associated with  new points appearing on a single plane */
        vector<vector<unsigned long> > d_planePoint; /*!< \brief std::vector containing points appearing on a single plane */

        vector<double> d_xCoordList;	/*!< \brief std::vector containing points appearing on a single plane */
        //vector<Grid > d_newBound;            /*!< \brief Boundary std::vector for new periodic elements (primal grid information). */
        vector<unsigned long> d_numNewElemBound;			/*!< \brief Number of new periodic elements of the boundary. */
        
        /*
         *  Parmetis variables
         */
        vector<unsigned long>  d_adjacency;
        vector<unsigned long>  d_xadj;
        unsigned long d_local_node;
        unsigned long d_local_elem;
        unsigned long d_xadj_size;
        unsigned long d_adjacency_size;
        vector<unsigned long> d_starting_node;
        vector<unsigned long> d_ending_node;
        vector<unsigned long> d_npoint_procs;
        unsigned long d_no_of_local_elements;
        vector<long> d_global_to_local_elem; 
    };
}

#endif
 








