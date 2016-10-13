












#include "../MACRO.hpp"
#include "GEOM_Geometry.hpp"

namespace ARIES
{
    namespace GEOM
    {
        GEOM_Geometry::GEOM_Geometry(void)
        {
            nEdge = 0;
            nPoint = 0;
            nElem = 0;

            nElem_Bound = NULL;
            Tag_to_Marker = NULL;
            elem = NULL;
            face = NULL;
            bound = NULL;
            node = NULL;
            edge = NULL;
            vertex = NULL;
            nVertex = NULL;
            newBound = NULL;
            nNewElem_Bound = NULL;
            Marker_All_SendRecv = NULL;

            //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].clear();
            //	PeriodicElem[MAX_NUMBER_PERIODIC].clear();
            //	XCoordList.clear();

            //	Xcoord_plane.clear();
            //	Ycoord_plane.clear();
            //	Zcoord_plane.clear();
            //	FaceArea_plane.clear();
            //	Plane_points.clear();
        }

        GEOM_Geometry::~GEOM_Geometry(void)
        {
            unsigned long iElem, iElem_Bound, iFace, iVertex, iEdge;
            unsigned short iMarker;

            if (elem != NULL)
            {
                for (iElem = 0; iElem < nElem; iElem++)
                    if (elem[iElem] != NULL) delete elem[iElem];
                delete[] elem;
            }

            if (bound != NULL)
            {
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                {
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
                    }
                }
                delete[] bound;
            }

            if (face != NULL)
            {
                for (iFace = 0; iFace < nFace; iFace++)
                    if (face[iFace] != NULL) delete face[iFace];
                delete[] face;
            }

            //  if (node != NULL) {
            //    for (iPoint = 0; iPoint < nPoint; iPoint ++)
            //      if (node[iPoint] != NULL) delete node[iPoint];
            //    delete[] node;
            //  }

            if (edge != NULL)
            {
                for (iEdge = 0; iEdge < nEdge; iEdge++)
                    if (edge[iEdge] != NULL) delete edge[iEdge];
                delete[] edge;
            }

            if (vertex != NULL)
            {
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                {
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        if (vertex[iMarker][iVertex] != NULL) delete vertex[iMarker][iVertex];
                    }
                }
                delete[] vertex;
            }

            if (newBound != NULL)
            {
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                {
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        if (newBound[iMarker][iElem_Bound] != NULL) delete newBound[iMarker][iElem_Bound];
                    }
                }
                delete[] newBound;
            }

            if (nElem_Bound != NULL) delete[] nElem_Bound;
            if (nVertex != NULL) delete[] nVertex;
            if (nNewElem_Bound != NULL) delete[] nNewElem_Bound;
            if (Marker_All_SendRecv != NULL) delete[] Marker_All_SendRecv;
            if (Tag_to_Marker != NULL) delete[] Tag_to_Marker;

            //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].~vector();
            //	PeriodicElem[MAX_NUMBER_PERIODIC].~vector();
            //	XCoordList.~vector();

            //	Xcoord_plane.~vector()
            //	Ycoord_plane.~vector()
            //	Zcoord_plane.~vector()
            //	FaceArea_plane.~vector()
            //	Plane_points.~vector()
        }



        long GEOM_Geometry::FindEdge(unsigned long first_point, unsigned long second_point)
        {
            unsigned long iPoint = 0;
            unsigned short iNode;
            for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++)
            {
                iPoint = node[first_point]->GetPoint(iNode);
                if (iPoint == second_point) break;
            }

            if (iPoint == second_point)
                return node[first_point]->GetEdge(iNode);
            else
            {
                std::cout << "\n\n   !!! Error !!!\n" << std::endl;
                std::cout << "Can't find the edge that connects " << first_point << " and " << second_point << "." << std::endl;
                exit(EXIT_FAILURE);
                return -1;
            }
        }

        bool GEOM_Geometry::CheckEdge(unsigned long first_point, unsigned long second_point)
        {
            unsigned long iPoint = 0;
            unsigned short iNode;
            for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++)
            {
                iPoint = node[first_point]->GetPoint(iNode);
                if (iPoint == second_point) break;
            }

            if (iPoint == second_point)
                return true;
            else
                return false;
        }

        void GEOM_Geometry::SetEdges(void)
        {
            unsigned long iPoint, jPoint;
            long iEdge;
            unsigned short jNode, iNode;
            long TestEdge = 0;

            nEdge = 0;
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++)
                {
                    jPoint = node[iPoint]->GetPoint(iNode);
                    for (jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
                        if (node[jPoint]->GetPoint(jNode) == iPoint)
                        {
                            TestEdge = node[jPoint]->GetEdge(jNode);
                            break;
                        }
                    if (TestEdge == -1)
                    {
                        node[iPoint]->SetEdge(nEdge, iNode);
                        node[jPoint]->SetEdge(nEdge, jNode);
                        nEdge++;
                    }
                }

            edge = new GRID::GRID_DGEdge*[nEdge];

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++)
                {
                    jPoint = node[iPoint]->GetPoint(iNode);
                    iEdge = FindEdge(iPoint, jPoint);
                    if (iPoint < jPoint) edge[iEdge] = new GRID::GRID_DGEdge(iPoint, jPoint, nDim);
                }
        }

        void GEOM_Geometry::SetFaces(void)
        {
            //	unsigned long iPoint, jPoint, iFace;
            //	unsigned short jNode, iNode;
            //	long TestFace = 0;
            //
            //	nFace = 0;
            //	for (iPoint = 0; iPoint < nPoint; iPoint++)
            //		for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
            //			jPoint = node[iPoint]->GetPoint(iNode);
            //			for (jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
            //				if (node[jPoint]->GetPoint(jNode) == iPoint) {
            //					TestFace = node[jPoint]->GetFace(jNode);
            //					break;
            //				}
            //			if (TestFace == -1) {
            //				node[iPoint]->SetFace(nFace, iNode);
            //				node[jPoint]->SetFace(nFace, jNode);
            //				nFace++;
            //			}
            //		}
            //
            //	face = new CFace*[nFace];
            //
            //	for (iPoint = 0; iPoint < nPoint; iPoint++)
            //		for (iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
            //			jPoint = node[iPoint]->GetPoint(iNode);
            //			iFace = FindFace(iPoint, jPoint);
            //			if (iPoint < jPoint) face[iFace] = new CFace(iPoint, jPoint, nDim);
            //		}
        }

        void GEOM_Geometry::TestGeometry(void)
        {
            std::ofstream para_file;

            para_file.open("test_geometry.dat", std::ios::out);

            double *Normal = new double[nDim];

            for (unsigned long iEdge = 0; iEdge < nEdge; iEdge++)
            {
                para_file << "Edge index: " << iEdge << std::endl;
                para_file << "   Point index: " << edge[iEdge]->GetNode(0) << "\t" << edge[iEdge]->GetNode(1) << std::endl;
                edge[iEdge]->GetNormal(Normal);
                para_file << "      Face normal : ";
                for (unsigned short iDim = 0; iDim < nDim; iDim++)
                    para_file << Normal[iDim] << "\t";
                para_file << std::endl;
            }

            para_file << std::endl;
            para_file << std::endl;
            para_file << std::endl;
            para_file << std::endl;

            for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++)
            {
                para_file << "Marker index: " << iMarker << std::endl;
                for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                {
                    para_file << "   Vertex index: " << iVertex << std::endl;
                    para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << std::endl;
                    para_file << "      Point coordinates : ";
                    for (unsigned short iDim = 0; iDim < nDim; iDim++)
                    {
                        para_file << node[vertex[iMarker][iVertex]->GetNode()]->GetCoord(iDim) << "\t";
                    }
                    para_file << std::endl;
                    vertex[iMarker][iVertex]->GetNormal(Normal);
                    para_file << "         Face normal : ";
                    for (unsigned short iDim = 0; iDim < nDim; iDim++)
                        para_file << Normal[iDim] << "\t";
                    para_file << std::endl;
                }
            }
        }
        
}
