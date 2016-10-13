




#include "GEOM_GeometryPhysical.hpp"

#include "../Grid/GRID_DualGrid.hpp"

#include "../Grid/GRID_VertexMPI.hpp"
#include "../Grid/GRID_Line.hpp"
#include "../Grid/GRID_Triangle.hpp"
#include "../Grid/GRID_Rectangle.hpp"
#include "../Grid/GRID_Tetrahedron.hpp"
#include "../Grid/GRID_Hexahedron.hpp"
#include "../Grid/GRID_Prism.hpp"
#include "../Grid/GRID_Pyramid.hpp"

namespace ARIES
{
    namespace GEOM
    {
        GEOM_GeometryPhysical::GEOM_GeometryPhysical() : GEOM_Geometry()
        {
            Global_to_Local_Point = NULL;
            Local_to_Global_Point = NULL;
            Local_to_Global_Marker = NULL;
            Global_to_Local_Marker = NULL;
        }

        GEOM_GeometryPhysical::GEOM_GeometryPhysical(TBOX::TBOX_Config *config, unsigned short val_iZone, unsigned short val_nZone) : GEOM_Geometry()
        {
            Global_to_Local_Point = NULL;
            Local_to_Global_Point = NULL;
            Local_to_Global_Marker = NULL;
            Global_to_Local_Marker = NULL;

            std::string text_line, Marker_Tag;
            std::ifstream mesh_file;
            unsigned short iDim, iMarker, iNodes;
            unsigned long iPoint, LocaNodes, iElem_Bound;
            double *NewCoord;
            nZone = val_nZone;
            std::ofstream boundary_file;
            std::string Grid_Marker;

            int rank = TBOX::MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            std::string val_mesh_filename = config->GetMesh_FileName();
            unsigned short val_format = config->GetMesh_FileFormat();

            /*--- Initialize counters for local/global points & elements ---*/
            if (rank == TBOX::MASTER_NODE)
                std::cout << std::endl << "---------------------- Read Grid File Information -----------------------" << std::endl;

            switch (val_format)
            {
            case TBOX::SU2:
                Read_SU2_Format_Parallel(config, val_mesh_filename, val_iZone, val_nZone);
                LocaNodes = local_node;
                break;
            case TBOX::CGNS:
                Read_CGNS_Format_Parallel(config, val_mesh_filename, val_iZone, val_nZone);
                LocaNodes = local_node;
                break;
            default:
                if (rank == TBOX::MASTER_NODE)
                    std::cout << "Unrecognized mesh format specified!" << std::endl;
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD, 1);
                MPI_Finalize();
#endif
                break;
            }

            /*--- Loop over the points element to re-scale the mesh, and plot it (only SU2_CFD) ---*/
            if (config->GetKind_SU2() == TBOX::SU2_CFD)
            {
                NewCoord = new double[nDim];

                /*--- The US system uses feet, but SU2 assumes that the grid is in inches ---*/
                if (config->GetSystemMeasurements() == TBOX::US)
                {
                    for (iPoint = 0; iPoint < LocaNodes; iPoint++)
                    {
                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            NewCoord[iDim] = node[iPoint]->GetCoord(iDim) / 12.0;
                        }
                        node[iPoint]->SetCoord(NewCoord);
                    }
                }

                delete[] NewCoord;
            }

            /*--- If SU2_DEF then write a file with the boundary information ---*/
            if ((config->GetKind_SU2() == TBOX::SU2_DEF) && (rank == TBOX::MASTER_NODE))
            {
                /*--- Open .su2 grid file ---*/
                boundary_file.open("boundary.su2", std::ios::out);

                /*--- Loop through and write the boundary info ---*/
                boundary_file << "NMARK= " << nMarker << std::endl;

                for (iMarker = 0; iMarker < nMarker; iMarker++)
                {
                    Grid_Marker = config->GetMarker_All_TagBound(iMarker);
                    boundary_file << "MARKER_TAG= " << Grid_Marker << std::endl;
                    boundary_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] << std::endl;

                    if (nDim == 2)
                    {
                        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                        {
                            boundary_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t";
                            for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
                                boundary_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t";
                            boundary_file << iElem_Bound << std::endl;
                        }
                    }

                    if (nDim == 3)
                    {
                        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                        {
                            boundary_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t";
                            for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
                                boundary_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t";
                            boundary_file << iElem_Bound << std::endl;
                        }
                    }
                }

                boundary_file.close();
            }
        }

        GEOM_GeometryPhysical::GEOM_GeometryPhysical(GEOM_Geometry *geometry, TBOX::TBOX_Config *config)
        {
            Global_to_Local_Point = NULL;
            Local_to_Global_Point = NULL;
            Local_to_Global_Marker = NULL;
            Global_to_Local_Marker = NULL;

            unsigned long iter, iPoint, jPoint, iElem, jElem, iVertex;
            unsigned long nElemTotal = 0, nPointTotal = 0, nPointDomainTotal = 0, nPointGhost = 0, nPointPeriodic = 0, nElemTriangle = 0, nElemRectangle = 0, nElemTetrahedron = 0, nElemHexahedron = 0, nElemPrism = 0, nElemPyramid = 0;
            unsigned long iElemTotal, iPointTotal, iPointGhost, iPointDomain, iPointPeriodic, iElemTriangle, iElemRectangle, iElemTetrahedron, iElemHexahedron, iElemPrism, iElemPyramid;
            unsigned long nBoundLineTotal = 0, iBoundLineTotal;
            unsigned long nBoundTriangleTotal = 0, iBoundTriangleTotal;
            unsigned long nBoundRectangleTotal = 0, iBoundRectangleTotal;
            unsigned long ReceptorColor = 0, DonorColor = 0, Transformation;
            unsigned long *nElem_Color = NULL, **Elem_Color = NULL, Max_nElem_Color = 0;
            unsigned long nTotalSendDomain_Periodic = 0, iTotalSendDomain_Periodic = 0, nTotalReceivedDomain_Periodic = 0, iTotalReceivedDomain_Periodic = 0, *nSendDomain_Periodic = NULL, *nReceivedDomain_Periodic = NULL;
            unsigned long Buffer_Send_nPointTotal = 0, Buffer_Send_nPointDomainTotal = 0, Buffer_Send_nPointGhost = 0, Buffer_Send_nPointPeriodic = 0;
            unsigned long Buffer_Send_nElemTotal, Buffer_Send_nElemTriangle = 0, Buffer_Send_nElemRectangle = 0, Buffer_Send_nElemTetrahedron = 0, Buffer_Send_nElemHexahedron = 0, Buffer_Send_nElemPrism = 0, Buffer_Send_nElemPyramid = 0;
            unsigned long Buffer_Send_nTotalSendDomain_Periodic = 0, Buffer_Send_nTotalReceivedDomain_Periodic = 0, *Buffer_Send_nSendDomain_Periodic = NULL, *Buffer_Send_nReceivedDomain_Periodic = NULL;
            unsigned long Buffer_Send_nBoundLineTotal = 0, Buffer_Send_nBoundTriangleTotal = 0, Buffer_Send_nBoundRectangleTotal = 0;
            unsigned long iVertexDomain, iBoundLine, iBoundTriangle, iBoundRectangle;

            /*--- Need to double-check these shorts in case we go to nprocs > ~32,000 ---*/
            unsigned short iNode, iDim, iMarker, jMarker, nMarkerDomain = 0, iMarkerDomain;
            unsigned short nDomain = 0, iDomain, jDomain, nPeriodic = 0, iPeriodic, overhead = 4, Buffer_Send_nMarkerDomain = 0, Buffer_Send_nDim = 0, Buffer_Send_nZone = 0, Buffer_Send_nPeriodic = 0;

            bool *MarkerIn = NULL, **VertexIn = NULL, CheckDomain;
            long vnodes_local[8], *Global2Local_Point = NULL;
            std::vector<long> DomainList;
            short *Marker_All_SendRecv_Copy = NULL;
            std::string *Marker_All_TagBound_Copy = NULL;
            unsigned short nMarker_Max = config->GetnMarker_Max();

            int rank = TBOX::MASTER_NODE;
            int size = TBOX::SINGLE_NODE;

            /*--- Some dynamic arrays so we're not allocating too much on the stack ---*/
            unsigned long *nVertexDomain = new unsigned long[nMarker_Max];
            unsigned long *nBoundLine = new unsigned long[nMarker_Max];
            unsigned long *nBoundTriangle = new unsigned long[nMarker_Max];
            unsigned long *nBoundRectangle = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nVertexDomain = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nBoundLine = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nBoundTriangle = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nBoundRectangle = new unsigned long[nMarker_Max];
            short *Buffer_Send_Marker_All_SendRecv = new short[nMarker_Max];
            char *Marker_All_TagBound = new char[nMarker_Max*TBOX::MAX_STRING_SIZE];
            char *Buffer_Send_Marker_All_TagBound = new char[nMarker_Max*TBOX::MAX_STRING_SIZE];

#ifdef HAVE_MPI
            /*--- MPI initialization ---*/
            MPI_Status status;
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            /*--- MPI status and request arrays for non-blocking communications ---*/
            MPI_Status send_stat[31], recv_stat[31];
            MPI_Request send_req[31], recv_req[31];
#endif

            /*--- Define buffer std::vector interior domain ---*/
            double        *Buffer_Send_Coord = NULL, *Buffer_Receive_Coord = NULL;
            unsigned long *Buffer_Send_Color = NULL, *Buffer_Receive_Color = NULL;
            unsigned long *Buffer_Send_GlobalPointIndex = NULL, *Buffer_Receive_GlobalPointIndex = NULL;
            unsigned long *Buffer_Send_Triangle = NULL, *Buffer_Receive_Triangle = NULL;
            unsigned long *Buffer_Send_Rectangle = NULL, *Buffer_Receive_Rectangle = NULL;
            unsigned long *Buffer_Send_Tetrahedron = NULL, *Buffer_Receive_Tetrahedron = NULL;
            unsigned long *Buffer_Send_Hexahedron = NULL, *Buffer_Receive_Hexahedron = NULL;
            unsigned long *Buffer_Send_Prism = NULL, *Buffer_Receive_Prism = NULL;
            unsigned long *Buffer_Send_Pyramid = NULL, *Buffer_Receive_Pyramid = NULL;

            /*--- Define buffer std::vector boundary ---*/
            unsigned long *Buffer_Send_BoundLine = NULL, *Buffer_Receive_BoundLine = NULL;
            unsigned long *Buffer_Send_BoundTriangle = NULL, *Buffer_Receive_BoundTriangle = NULL;
            unsigned long *Buffer_Send_BoundRectangle = NULL, *Buffer_Receive_BoundRectangle = NULL;
            unsigned long *Buffer_Send_Local2Global_Marker = NULL, *Buffer_Receive_Local2Global_Marker = NULL;

            /*--- Define buffer std::vector periodic boundary conditions ---*/
            double *Buffer_Send_Center = NULL, *Buffer_Receive_Center = NULL;
            double *Buffer_Send_Rotation = NULL, *Buffer_Receive_Rotation = NULL;
            double *Buffer_Send_Translate = NULL, *Buffer_Receive_Translate = NULL;

            /*--- Define buffer std::vector periodic boundary conditions ---*/
            unsigned long *Buffer_Send_SendDomain_Periodic = NULL, *Buffer_Receive_SendDomain_Periodic = NULL;
            unsigned long *Buffer_Send_SendDomain_PeriodicTrans = NULL, *Buffer_Receive_SendDomain_PeriodicTrans = NULL;
            unsigned long *Buffer_Send_SendDomain_PeriodicReceptor = NULL, *Buffer_Receive_SendDomain_PeriodicReceptor = NULL;
            unsigned long *Buffer_Send_ReceivedDomain_Periodic = NULL, *Buffer_Receive_ReceivedDomain_Periodic = NULL;
            unsigned long *Buffer_Send_ReceivedDomain_PeriodicTrans = NULL, *Buffer_Receive_ReceivedDomain_PeriodicTrans = NULL;
            unsigned long *Buffer_Send_ReceivedDomain_PeriodicDonor = NULL, *Buffer_Receive_ReceivedDomain_PeriodicDonor = NULL;


            /*--- Basic dimensionalization ---*/
            nDomain = size;
            Marker_All_SendRecv = new short[nMarker_Max];
            nSendDomain_Periodic = new unsigned long[nDomain];
            nReceivedDomain_Periodic = new unsigned long[nDomain];

            /*--- Auxiliar std::vector based on the original geometry ---*/

            if (rank == TBOX::MASTER_NODE)
            {
                MarkerIn = new bool[geometry->GetnMarker()];

                VertexIn = new bool*[geometry->GetnMarker()];
                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    VertexIn[iMarker] = new bool[geometry->GetnElem_Bound(iMarker)];

                Global2Local_Point = new long[geometry->GetnPoint()];

                Buffer_Send_nDim = geometry->GetnDim();
                Buffer_Send_nZone = geometry->GetnZone();

                Buffer_Send_nPeriodic = config->GetnPeriodicIndex();
                Buffer_Send_Center = new double[Buffer_Send_nPeriodic * 3];
                Buffer_Send_Rotation = new double[Buffer_Send_nPeriodic * 3];
                Buffer_Send_Translate = new double[Buffer_Send_nPeriodic * 3];

                Buffer_Send_nSendDomain_Periodic = new unsigned long[nDomain];
                Buffer_Send_nReceivedDomain_Periodic = new unsigned long[nDomain];

                /*--- Divide the elements in color list to speed up the grid partitioning ---*/
                nElem_Color = new unsigned long[nDomain];
                for (iDomain = 0; iDomain < nDomain; iDomain++) nElem_Color[iDomain] = 0;

                for (iElem = 0; iElem < geometry->GetnElem(); iElem++)
                {
                    DomainList.clear();
                    for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = geometry->elem[iElem]->GetNode(iNode);
                        iDomain = geometry->node[iPoint]->GetColor();

                        CheckDomain = true;
                        for (jDomain = 0; jDomain < DomainList.size(); jDomain++)
                        {
                            if (DomainList[jDomain] == iDomain)
                            {
                                CheckDomain = false;
                                break;
                            }
                        }

                        /*--- If the element is not in the list, then add it ---*/
                        if (CheckDomain)
                        {
                            DomainList.push_back(iDomain);
                            nElem_Color[iDomain]++;
                        }
                    }
                }

                /*--- Find the maximum number of elements per color to allocate the list ---*/

                Max_nElem_Color = 0;
                for (iDomain = 0; iDomain < nDomain; iDomain++)
                {
                    if (nElem_Color[iDomain] > Max_nElem_Color) Max_nElem_Color = nElem_Color[iDomain];
                }

                /*--- Allocate the element color array ---*/
                Elem_Color = new unsigned long*[nDomain];
                for (iDomain = 0; iDomain < nDomain; iDomain++)
                {
                    Elem_Color[iDomain] = new unsigned long[Max_nElem_Color];
                    nElem_Color[iDomain] = 0;
                }

                /*--- Create the lement list based on the color ---*/
                for (iElem = 0; iElem < geometry->GetnElem(); iElem++)
                {
                    DomainList.clear();
                    for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = geometry->elem[iElem]->GetNode(iNode);
                        iDomain = geometry->node[iPoint]->GetColor();

                        /*--- Check if the element has been already added to the color ---*/
                        CheckDomain = true;
                        for (jDomain = 0; jDomain < DomainList.size(); jDomain++)
                        {
                            if (DomainList[jDomain] == iDomain)
                            {
                                CheckDomain = false;
                                break;
                            }
                        }

                        if (CheckDomain)
                        {
                            DomainList.push_back(iDomain);
                            Elem_Color[iDomain][nElem_Color[iDomain]] = iElem;
                            nElem_Color[iDomain]++;
                        }
                    }
                }


                /*--- Create a local copy of config->GetMarker_All_SendRecv and
                config->GetMarker_All_TagBound in the master node ---*/
                Marker_All_SendRecv_Copy = new short[geometry->GetnMarker()];
                Marker_All_TagBound_Copy = new std::string[geometry->GetnMarker()];

                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                {
                    Marker_All_SendRecv_Copy[iMarker] = config->GetMarker_All_SendRecv(iMarker);
                    Marker_All_TagBound_Copy[iMarker] = config->GetMarker_All_TagBound(iMarker);
                }
            }

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                if (rank == TBOX::MASTER_NODE)
                {
                    /*--- Interior dimensionalization. Loop over the original grid to perform the
                    dimensionalizaton of the domain variables ---*/
                    Buffer_Send_nElemTotal = 0; Buffer_Send_nPointTotal = 0; Buffer_Send_nPointGhost = 0; Buffer_Send_nPointDomainTotal = 0; Buffer_Send_nPointPeriodic = 0;
                    Buffer_Send_nElemTriangle = 0; Buffer_Send_nElemRectangle = 0; Buffer_Send_nElemTetrahedron = 0; Buffer_Send_nElemHexahedron = 0; Buffer_Send_nElemPrism = 0; Buffer_Send_nElemPyramid = 0;

                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Global2Local_Point[iPoint] = -1;

                    for (jElem = 0; jElem < nElem_Color[iDomain]; jElem++)
                    {

                        iElem = Elem_Color[iDomain][jElem];

                        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                        {
                            iPoint = geometry->elem[iElem]->GetNode(iNode);
                            if (Global2Local_Point[iPoint] == -1)
                            {
                                Global2Local_Point[iPoint] = 1;
                                Buffer_Send_nPointTotal++;
                                if (geometry->node[iPoint]->GetColor() != iDomain) Buffer_Send_nPointGhost++;
                                else
                                {
                                    if (iPoint > geometry->GetnPointDomain() - 1)
                                    {
                                        Buffer_Send_nPointGhost++;
                                        Buffer_Send_nPointPeriodic++;
                                    }
                                    else Buffer_Send_nPointDomainTotal++;
                                }
                            }
                        }

                        switch (geometry->elem[iElem]->GetVTK_Type())
                        {
                        case TBOX::TRIANGLE: Buffer_Send_nElemTriangle++; break;
                        case TBOX::RECTANGLE: Buffer_Send_nElemRectangle++; break;
                        case TBOX::TETRAHEDRON: Buffer_Send_nElemTetrahedron++; break;
                        case TBOX::HEXAHEDRON: Buffer_Send_nElemHexahedron++; break;
                        case TBOX::PRISM: Buffer_Send_nElemPrism++; break;
                        case TBOX::PYRAMID: Buffer_Send_nElemPyramid++; break;
                        }
                        Buffer_Send_nElemTotal++;
                    }

                    /*--- Boundary dimensionalization. Dimensionalization with physical boundaries, compute Buffer_Send_nMarkerDomain,
                    Buffer_Send_nVertexDomain[nMarkerDomain] ---*/

                    Buffer_Send_nMarkerDomain = 0; Buffer_Send_nBoundLineTotal = 0; Buffer_Send_nBoundTriangleTotal = 0; Buffer_Send_nBoundRectangleTotal = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        Buffer_Send_nVertexDomain[iMarker] = 0;
                        Buffer_Send_nBoundLine[iMarker] = 0;
                        Buffer_Send_nBoundTriangle[iMarker] = 0;
                        Buffer_Send_nBoundRectangle[iMarker] = 0;

                        Buffer_Send_Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
                        sprintf(&Buffer_Send_Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE], "%s", Marker_All_TagBound_Copy[iMarker].c_str());
                    }

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
                        {
                            MarkerIn[iMarker] = false; Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain] = 0;

                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                VertexIn[iMarker][iVertex] = false;
                                for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++)
                                {
                                    iPoint = geometry->bound[iMarker][iVertex]->GetNode(iNode);
                                    if (geometry->node[iPoint]->GetColor() == iDomain) VertexIn[iMarker][iVertex] = true;
                                }

                                if (VertexIn[iMarker][iVertex])
                                {
                                    switch (geometry->bound[iMarker][iVertex]->GetVTK_Type())
                                    {
                                    case TBOX::LINE: Buffer_Send_nBoundLine[Buffer_Send_nMarkerDomain]++; Buffer_Send_nBoundLineTotal++; break;
                                    case TBOX::TRIANGLE: Buffer_Send_nBoundTriangle[Buffer_Send_nMarkerDomain]++; Buffer_Send_nBoundTriangleTotal++; break;
                                    case TBOX::RECTANGLE: Buffer_Send_nBoundRectangle[Buffer_Send_nMarkerDomain]++; Buffer_Send_nBoundRectangleTotal++; break;
                                    }

                                    Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain] ++;
                                    MarkerIn[iMarker] = true;
                                }
                            }

                            if (MarkerIn[iMarker])
                            {
                                Buffer_Send_nMarkerDomain++;
                            }

                        }
                    }

                    /*--- Copy periodic information from the config file ---*/
                    for (iPeriodic = 0; iPeriodic < Buffer_Send_nPeriodic; iPeriodic++)
                    {
                        for (iDim = 0; iDim < 3; iDim++)
                        {
                            Buffer_Send_Center[iDim + iPeriodic * 3] = config->GetPeriodicCenter(iPeriodic)[iDim];
                            Buffer_Send_Rotation[iDim + iPeriodic * 3] = config->GetPeriodicRotation(iPeriodic)[iDim];
                            Buffer_Send_Translate[iDim + iPeriodic * 3] = config->GetPeriodicTranslate(iPeriodic)[iDim];
                        }
                    }

                    /*--- Dimensionalization of the periodic auxiliar std::vectors ---*/
                    for (jDomain = 0; jDomain < nDomain; jDomain++)
                    {
                        Buffer_Send_nSendDomain_Periodic[jDomain] = 0;
                        Buffer_Send_nReceivedDomain_Periodic[jDomain] = 0;
                    }
                    Buffer_Send_nTotalSendDomain_Periodic = 0;
                    Buffer_Send_nTotalReceivedDomain_Periodic = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE)
                        {
                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
                                if (iDomain == geometry->node[iPoint]->GetColor())
                                {
                                    if (config->GetMarker_All_SendRecv(iMarker) > 0)
                                    {
                                        /*--- Identify the color of the receptor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                ReceptorColor = geometry->node[jPoint]->GetColor();
                                            }
                                        }

                                        Buffer_Send_nSendDomain_Periodic[ReceptorColor]++;
                                        Buffer_Send_nTotalSendDomain_Periodic++;

                                    }
                                    if (config->GetMarker_All_SendRecv(iMarker) < 0)
                                    {
                                        /*--- Identify the color of the donor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                DonorColor = geometry->node[jPoint]->GetColor();
                                            }
                                        }

                                        Buffer_Send_nReceivedDomain_Periodic[DonorColor]++;
                                        Buffer_Send_nTotalReceivedDomain_Periodic++;
                                    }
                                }
                            }
                        }
                    }

                    /*--- Allocate the buffer std::vectors in the appropiate domain (master, iDomain) ---*/
                    Buffer_Send_Coord = new double[Buffer_Send_nPointTotal*Buffer_Send_nDim];
                    Buffer_Send_Color = new unsigned long[Buffer_Send_nPointTotal];
                    Buffer_Send_GlobalPointIndex = new unsigned long[Buffer_Send_nPointTotal];
                    Buffer_Send_Triangle = new unsigned long[Buffer_Send_nElemTriangle * 3];
                    Buffer_Send_Rectangle = new unsigned long[Buffer_Send_nElemRectangle * 4];
                    Buffer_Send_Tetrahedron = new unsigned long[Buffer_Send_nElemTetrahedron * 4];
                    Buffer_Send_Hexahedron = new unsigned long[Buffer_Send_nElemHexahedron * 8];
                    Buffer_Send_Prism = new unsigned long[Buffer_Send_nElemPrism * 6];
                    Buffer_Send_Pyramid = new unsigned long[Buffer_Send_nElemPyramid * 5];

                    Buffer_Send_BoundLine = new unsigned long[Buffer_Send_nBoundLineTotal * 2];
                    Buffer_Send_BoundTriangle = new unsigned long[Buffer_Send_nBoundTriangleTotal * 3];
                    Buffer_Send_BoundRectangle = new unsigned long[Buffer_Send_nBoundRectangleTotal * 4];
                    Buffer_Send_Local2Global_Marker = new unsigned long[Buffer_Send_nMarkerDomain];

                    Buffer_Send_SendDomain_Periodic = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
                    Buffer_Send_SendDomain_PeriodicTrans = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
                    Buffer_Send_SendDomain_PeriodicReceptor = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
                    Buffer_Send_ReceivedDomain_Periodic = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
                    Buffer_Send_ReceivedDomain_PeriodicTrans = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
                    Buffer_Send_ReceivedDomain_PeriodicDonor = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];

                    if (iDomain != TBOX::MASTER_NODE)
                    {
#ifdef HAVE_MPI

                        MPI_Isend(&Buffer_Send_nDim, 1, MPI_UNSIGNED_SHORT, iDomain, 0, MPI_COMM_WORLD, &send_req[0]);
                        MPI_Isend(&Buffer_Send_nZone, 1, MPI_UNSIGNED_SHORT, iDomain, 1, MPI_COMM_WORLD, &send_req[1]);
                        MPI_Isend(&Buffer_Send_nPointTotal, 1, MPI_UNSIGNED_LONG, iDomain, 2, MPI_COMM_WORLD, &send_req[2]);
                        MPI_Isend(&Buffer_Send_nPointDomainTotal, 1, MPI_UNSIGNED_LONG, iDomain, 3, MPI_COMM_WORLD, &send_req[3]);
                        MPI_Isend(&Buffer_Send_nPointGhost, 1, MPI_UNSIGNED_LONG, iDomain, 4, MPI_COMM_WORLD, &send_req[4]);
                        MPI_Isend(&Buffer_Send_nPointPeriodic, 1, MPI_UNSIGNED_LONG, iDomain, 5, MPI_COMM_WORLD, &send_req[5]);
                        MPI_Isend(&Buffer_Send_nElemTotal, 1, MPI_UNSIGNED_LONG, iDomain, 6, MPI_COMM_WORLD, &send_req[6]);
                        MPI_Isend(&Buffer_Send_nElemTriangle, 1, MPI_UNSIGNED_LONG, iDomain, 7, MPI_COMM_WORLD, &send_req[7]);
                        MPI_Isend(&Buffer_Send_nElemRectangle, 1, MPI_UNSIGNED_LONG, iDomain, 8, MPI_COMM_WORLD, &send_req[8]);
                        MPI_Isend(&Buffer_Send_nElemTetrahedron, 1, MPI_UNSIGNED_LONG, iDomain, 9, MPI_COMM_WORLD, &send_req[9]);
                        MPI_Isend(&Buffer_Send_nElemHexahedron, 1, MPI_UNSIGNED_LONG, iDomain, 10, MPI_COMM_WORLD, &send_req[10]);
                        MPI_Isend(&Buffer_Send_nElemPrism, 1, MPI_UNSIGNED_LONG, iDomain, 11, MPI_COMM_WORLD, &send_req[11]);
                        MPI_Isend(&Buffer_Send_nElemPyramid, 1, MPI_UNSIGNED_LONG, iDomain, 12, MPI_COMM_WORLD, &send_req[12]);

                        MPI_Isend(&Buffer_Send_nBoundLineTotal, 1, MPI_UNSIGNED_LONG, iDomain, 13, MPI_COMM_WORLD, &send_req[13]);
                        MPI_Isend(&Buffer_Send_nBoundTriangleTotal, 1, MPI_UNSIGNED_LONG, iDomain, 14, MPI_COMM_WORLD, &send_req[14]);
                        MPI_Isend(&Buffer_Send_nBoundRectangleTotal, 1, MPI_UNSIGNED_LONG, iDomain, 15, MPI_COMM_WORLD, &send_req[15]);
                        MPI_Isend(&Buffer_Send_nMarkerDomain, 1, MPI_UNSIGNED_SHORT, iDomain, 16, MPI_COMM_WORLD, &send_req[16]);
                        MPI_Isend(Buffer_Send_nVertexDomain, nMarker_Max, MPI_UNSIGNED_LONG, iDomain, 17, MPI_COMM_WORLD, &send_req[17]);
                        MPI_Isend(Buffer_Send_nBoundLine, nMarker_Max, MPI_UNSIGNED_LONG, iDomain, 18, MPI_COMM_WORLD, &send_req[18]);
                        MPI_Isend(Buffer_Send_nBoundTriangle, nMarker_Max, MPI_UNSIGNED_LONG, iDomain, 19, MPI_COMM_WORLD, &send_req[19]);
                        MPI_Isend(Buffer_Send_nBoundRectangle, nMarker_Max, MPI_UNSIGNED_LONG, iDomain, 20, MPI_COMM_WORLD, &send_req[20]);
                        MPI_Isend(Buffer_Send_Marker_All_SendRecv, nMarker_Max, MPI_SHORT, iDomain, 21, MPI_COMM_WORLD, &send_req[21]);
                        MPI_Isend(Buffer_Send_Marker_All_TagBound, nMarker_Max*MAX_STRING_SIZE, MPI_CHAR, iDomain, 22, MPI_COMM_WORLD, &send_req[22]);

                        MPI_Isend(&Buffer_Send_nPeriodic, 1, MPI_UNSIGNED_SHORT, iDomain, 23, MPI_COMM_WORLD, &send_req[23]);
                        MPI_Isend(Buffer_Send_Center, nPeriodic * 3, MPI_DOUBLE, iDomain, 24, MPI_COMM_WORLD, &send_req[24]);
                        MPI_Isend(Buffer_Send_Rotation, nPeriodic * 3, MPI_DOUBLE, iDomain, 25, MPI_COMM_WORLD, &send_req[25]);
                        MPI_Isend(Buffer_Send_Translate, nPeriodic * 3, MPI_DOUBLE, iDomain, 26, MPI_COMM_WORLD, &send_req[26]);

                        MPI_Isend(&Buffer_Send_nTotalSendDomain_Periodic, 1, MPI_UNSIGNED_LONG, iDomain, 27, MPI_COMM_WORLD, &send_req[27]);
                        MPI_Isend(&Buffer_Send_nTotalReceivedDomain_Periodic, 1, MPI_UNSIGNED_LONG, iDomain, 28, MPI_COMM_WORLD, &send_req[28]);
                        MPI_Isend(Buffer_Send_nSendDomain_Periodic, nDomain, MPI_UNSIGNED_LONG, iDomain, 29, MPI_COMM_WORLD, &send_req[29]);
                        MPI_Isend(Buffer_Send_nReceivedDomain_Periodic, nDomain, MPI_UNSIGNED_LONG, iDomain, 30, MPI_COMM_WORLD, &send_req[30]);

                        /*--- Wait for this set of non-blocking comm. to complete ---*/
                        MPI_Waitall(31, send_req, send_stat);
#endif
                    }
                    else
                    {
                        /*--- We are the master node, so simply copy values into place ---*/
                        nDim = Buffer_Send_nDim;
                        nZone = Buffer_Send_nZone;

                        nPeriodic = Buffer_Send_nPeriodic;

                        nPointTotal = Buffer_Send_nPointTotal;
                        nPointDomainTotal = Buffer_Send_nPointDomainTotal;
                        nPointGhost = Buffer_Send_nPointGhost;
                        nPointPeriodic = Buffer_Send_nPointPeriodic;

                        nElemTotal = Buffer_Send_nElemTotal;
                        nElemTriangle = Buffer_Send_nElemTriangle;
                        nElemRectangle = Buffer_Send_nElemRectangle;
                        nElemTetrahedron = Buffer_Send_nElemTetrahedron;
                        nElemHexahedron = Buffer_Send_nElemHexahedron;
                        nElemPrism = Buffer_Send_nElemPrism;
                        nElemPyramid = Buffer_Send_nElemPyramid;

                        nelem_triangle = nElemTriangle;
                        nelem_quad = nElemRectangle;
                        nelem_tetra = nElemTetrahedron;
                        nelem_hexa = nElemHexahedron;
                        nelem_prism = nElemPrism;
                        nelem_pyramid = nElemPyramid;

                        nBoundLineTotal = Buffer_Send_nBoundLineTotal;
                        nBoundTriangleTotal = Buffer_Send_nBoundTriangleTotal;
                        nBoundRectangleTotal = Buffer_Send_nBoundRectangleTotal;
                        nMarkerDomain = Buffer_Send_nMarkerDomain;

                        for (iMarker = 0; iMarker < nMarker_Max; iMarker++)
                        {
                            nVertexDomain[iMarker] = Buffer_Send_nVertexDomain[iMarker];
                            nBoundLine[iMarker] = Buffer_Send_nBoundLine[iMarker];
                            nBoundTriangle[iMarker] = Buffer_Send_nBoundTriangle[iMarker];
                            nBoundRectangle[iMarker] = Buffer_Send_nBoundRectangle[iMarker];
                            Marker_All_SendRecv[iMarker] = Buffer_Send_Marker_All_SendRecv[iMarker];
                            for (iter = 0; iter < TBOX::MAX_STRING_SIZE; iter++)
                                Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE + iter] = Buffer_Send_Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE + iter];
                        }

                        Buffer_Receive_Center = new double[nPeriodic * 3];
                        Buffer_Receive_Rotation = new double[nPeriodic * 3];
                        Buffer_Receive_Translate = new double[nPeriodic * 3];

                        for (iter = 0; iter < nPeriodic * 3; iter++)
                        {
                            Buffer_Receive_Center[iter] = Buffer_Send_Center[iter];
                            Buffer_Receive_Rotation[iter] = Buffer_Send_Rotation[iter];
                            Buffer_Receive_Translate[iter] = Buffer_Send_Translate[iter];
                        }

                        nTotalSendDomain_Periodic = Buffer_Send_nTotalSendDomain_Periodic;
                        nTotalReceivedDomain_Periodic = Buffer_Send_nTotalReceivedDomain_Periodic;

                        for (iter = 0; iter < nDomain; iter++)
                        {
                            nSendDomain_Periodic[iter] = Buffer_Send_nSendDomain_Periodic[iter];
                            nReceivedDomain_Periodic[iter] = Buffer_Send_nReceivedDomain_Periodic[iter];
                        }

                    }
                }

                /*--- Receive the size of buffers---*/

                if (rank == iDomain)
                {
                    /*--- Receive the size of buffers---*/
                    if (rank != TBOX::MASTER_NODE)
                    {
#ifdef HAVE_MPI
                        MPI_Irecv(&nDim, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, 0, MPI_COMM_WORLD, &recv_req[0]);
                        MPI_Irecv(&nZone, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, 1, MPI_COMM_WORLD, &recv_req[1]);
                        MPI_Irecv(&nPointTotal, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 2, MPI_COMM_WORLD, &recv_req[2]);
                        MPI_Irecv(&nPointDomainTotal, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 3, MPI_COMM_WORLD, &recv_req[3]);
                        MPI_Irecv(&nPointGhost, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 4, MPI_COMM_WORLD, &recv_req[4]);
                        MPI_Irecv(&nPointPeriodic, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 5, MPI_COMM_WORLD, &recv_req[5]);
                        MPI_Irecv(&nElemTotal, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 6, MPI_COMM_WORLD, &recv_req[6]);
                        MPI_Irecv(&nElemTriangle, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 7, MPI_COMM_WORLD, &recv_req[7]);
                        MPI_Irecv(&nElemRectangle, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 8, MPI_COMM_WORLD, &recv_req[8]);
                        MPI_Irecv(&nElemTetrahedron, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 9, MPI_COMM_WORLD, &recv_req[9]);
                        MPI_Irecv(&nElemHexahedron, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 10, MPI_COMM_WORLD, &recv_req[10]);
                        MPI_Irecv(&nElemPrism, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 11, MPI_COMM_WORLD, &recv_req[11]);
                        MPI_Irecv(&nElemPyramid, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 12, MPI_COMM_WORLD, &recv_req[12]);

                        MPI_Irecv(&nBoundLineTotal, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 13, MPI_COMM_WORLD, &recv_req[13]);
                        MPI_Irecv(&nBoundTriangleTotal, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 14, MPI_COMM_WORLD, &recv_req[14]);
                        MPI_Irecv(&nBoundRectangleTotal, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 15, MPI_COMM_WORLD, &recv_req[15]);
                        MPI_Irecv(&nMarkerDomain, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, 16, MPI_COMM_WORLD, &recv_req[16]);
                        MPI_Irecv(nVertexDomain, nMarker_Max, MPI_UNSIGNED_LONG, MASTER_NODE, 17, MPI_COMM_WORLD, &recv_req[17]);
                        MPI_Irecv(nBoundLine, nMarker_Max, MPI_UNSIGNED_LONG, MASTER_NODE, 18, MPI_COMM_WORLD, &recv_req[18]);
                        MPI_Irecv(nBoundTriangle, nMarker_Max, MPI_UNSIGNED_LONG, MASTER_NODE, 19, MPI_COMM_WORLD, &recv_req[19]);
                        MPI_Irecv(nBoundRectangle, nMarker_Max, MPI_UNSIGNED_LONG, MASTER_NODE, 20, MPI_COMM_WORLD, &recv_req[20]);
                        MPI_Irecv(Marker_All_SendRecv, nMarker_Max, MPI_SHORT, MASTER_NODE, 21, MPI_COMM_WORLD, &recv_req[21]);
                        MPI_Irecv(Marker_All_TagBound, nMarker_Max*MAX_STRING_SIZE, MPI_CHAR, MASTER_NODE, 22, MPI_COMM_WORLD, &recv_req[22]);
                        MPI_Irecv(&nPeriodic, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, 23, MPI_COMM_WORLD, &recv_req[23]);

                        /*--- Wait for the this set of non-blocking recv's to complete ---*/

                        MPI_Waitall(24, recv_req, recv_stat);
#endif

                        /*--- Update the number of elements (local) ---*/
                        nelem_triangle = nElemTriangle;
                        nelem_quad = nElemRectangle;
                        nelem_tetra = nElemTetrahedron;
                        nelem_hexa = nElemHexahedron;
                        nelem_prism = nElemPrism;
                        nelem_pyramid = nElemPyramid;

                        /*--- Marker_All_TagBound and Marker_All_SendRecv, set the same values in the config files of all the files ---*/
                        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                        {
                            config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
                            config->SetMarker_All_TagBound(iMarker, std::string(&Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE]));
                        }

                        /*--- Periodic boundary conditions, set the values in the config files of all the files ---*/
                        Buffer_Receive_Center = new double[nPeriodic * 3];
                        Buffer_Receive_Rotation = new double[nPeriodic * 3];
                        Buffer_Receive_Translate = new double[nPeriodic * 3];

#ifdef HAVE_MPI

                        MPI_Irecv(Buffer_Receive_Center, nPeriodic * 3, MPI_DOUBLE, MASTER_NODE, 24, MPI_COMM_WORLD, &recv_req[0]);
                        MPI_Irecv(Buffer_Receive_Rotation, nPeriodic * 3, MPI_DOUBLE, MASTER_NODE, 25, MPI_COMM_WORLD, &recv_req[1]);
                        MPI_Irecv(Buffer_Receive_Translate, nPeriodic * 3, MPI_DOUBLE, MASTER_NODE, 26, MPI_COMM_WORLD, &recv_req[2]);

                        MPI_Irecv(&nTotalSendDomain_Periodic, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 27, MPI_COMM_WORLD, &recv_req[3]);
                        MPI_Irecv(&nTotalReceivedDomain_Periodic, 1, MPI_UNSIGNED_LONG, MASTER_NODE, 28, MPI_COMM_WORLD, &recv_req[4]);
                        MPI_Irecv(nSendDomain_Periodic, nDomain, MPI_UNSIGNED_LONG, MASTER_NODE, 29, MPI_COMM_WORLD, &recv_req[5]);
                        MPI_Irecv(nReceivedDomain_Periodic, nDomain, MPI_UNSIGNED_LONG, MASTER_NODE, 30, MPI_COMM_WORLD, &recv_req[6]);

                        /*--- Wait for this set of non-blocking comm. to complete ---*/
                        MPI_Waitall(7, recv_req, recv_stat);
#endif

                        config->SetnPeriodicIndex(nPeriodic);

                        for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++)
                        {
                            double* center = new double[3];       // Do not deallocate the memory
                            double* rotation = new double[3];    // Do not deallocate the memory
                            double* translate = new double[3];    // Do not deallocate the memory

                            for (iDim = 0; iDim < 3; iDim++)
                            {
                                center[iDim] = Buffer_Receive_Center[iDim + iPeriodic * 3];
                                rotation[iDim] = Buffer_Receive_Rotation[iDim + iPeriodic * 3];
                                translate[iDim] = Buffer_Receive_Translate[iDim + iPeriodic * 3];
                            }
                            config->SetPeriodicCenter(iPeriodic, center);
                            config->SetPeriodicRotation(iPeriodic, rotation);
                            config->SetPeriodicTranslate(iPeriodic, translate);
                        }
                    }

                    delete[] Buffer_Receive_Center;
                    delete[] Buffer_Receive_Rotation;
                    delete[] Buffer_Receive_Translate;

                    /*--- Allocate the receive buffer std::vector ---*/
                    Buffer_Receive_Coord = new double[nPointTotal*nDim];
                    Buffer_Receive_Color = new unsigned long[nPointTotal];
                    Buffer_Receive_GlobalPointIndex = new unsigned long[nPointTotal];
                    Buffer_Receive_Triangle = new unsigned long[nElemTriangle * 3];
                    Buffer_Receive_Rectangle = new unsigned long[nElemRectangle * 4];
                    Buffer_Receive_Tetrahedron = new unsigned long[nElemTetrahedron * 4];
                    Buffer_Receive_Hexahedron = new unsigned long[nElemHexahedron * 8];
                    Buffer_Receive_Prism = new unsigned long[nElemPrism * 6];
                    Buffer_Receive_Pyramid = new unsigned long[nElemPyramid * 5];
                    Buffer_Receive_BoundLine = new unsigned long[nBoundLineTotal * 2];
                    Buffer_Receive_BoundTriangle = new unsigned long[nBoundTriangleTotal * 3];
                    Buffer_Receive_BoundRectangle = new unsigned long[nBoundRectangleTotal * 4];
                    Buffer_Receive_Local2Global_Marker = new unsigned long[nMarkerDomain];

                    Buffer_Receive_SendDomain_Periodic = new unsigned long[nTotalSendDomain_Periodic];
                    Buffer_Receive_SendDomain_PeriodicTrans = new unsigned long[nTotalSendDomain_Periodic];
                    Buffer_Receive_SendDomain_PeriodicReceptor = new unsigned long[nTotalSendDomain_Periodic];
                    Buffer_Receive_ReceivedDomain_Periodic = new unsigned long[nTotalReceivedDomain_Periodic];
                    Buffer_Receive_ReceivedDomain_PeriodicTrans = new unsigned long[nTotalReceivedDomain_Periodic];
                    Buffer_Receive_ReceivedDomain_PeriodicDonor = new unsigned long[nTotalReceivedDomain_Periodic];
                }

                /*--- Set the value of the Send buffers ---*/
                if (rank == TBOX::MASTER_NODE)
                {
                    /*--- Set the value of the interior geometry ---*/
                    iElemTotal = 0; iPointDomain = 0; iPointPeriodic = Buffer_Send_nPointDomainTotal; iPointGhost = Buffer_Send_nPointDomainTotal + Buffer_Send_nPointPeriodic;
                    iElemTriangle = 0; iElemRectangle = 0; iElemTetrahedron = 0; iElemHexahedron = 0; iElemPrism = 0; iElemPyramid = 0;

                    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                        Global2Local_Point[iPoint] = -1;

                    for (jElem = 0; jElem < nElem_Color[iDomain]; jElem++)
                    {
                        iElem = Elem_Color[iDomain][jElem];

                        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                        {
                            iPoint = geometry->elem[iElem]->GetNode(iNode);
                            if (Global2Local_Point[iPoint] == -1)
                            {
                                if (geometry->node[iPoint]->GetColor() == iDomain)
                                {
                                    if (iPoint > geometry->GetnPointDomain() - 1)
                                        iPointTotal = iPointPeriodic;
                                    else
                                        iPointTotal = iPointDomain;
                                }
                                else
                                    iPointTotal = iPointGhost;

                                Global2Local_Point[iPoint] = iPointTotal;
                                Buffer_Send_Color[iPointTotal] = geometry->node[iPoint]->GetColor();
                                Buffer_Send_GlobalPointIndex[iPointTotal] = iPoint;
                                for (iDim = 0; iDim < Buffer_Send_nDim; iDim++)
                                    Buffer_Send_Coord[Buffer_Send_nDim*iPointTotal + iDim] = geometry->node[iPoint]->GetCoord(iDim);

                                if (geometry->node[iPoint]->GetColor() == iDomain)
                                {
                                    if (iPoint > geometry->GetnPointDomain() - 1) iPointPeriodic++;
                                    else iPointDomain++;
                                }
                                else iPointGhost++;
                            }

                            vnodes_local[iNode] = Global2Local_Point[iPoint];
                        }

                        switch (geometry->elem[iElem]->GetVTK_Type())
                        {
                        case TBOX::TRIANGLE:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Triangle[3 * iElemTriangle + iNode] = vnodes_local[iNode];
                            iElemTriangle++; break;
                        case TBOX::RECTANGLE:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Rectangle[4 * iElemRectangle + iNode] = vnodes_local[iNode];
                            iElemRectangle++; break;
                        case TBOX::TETRAHEDRON:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Tetrahedron[4 * iElemTetrahedron + iNode] = vnodes_local[iNode];
                            iElemTetrahedron++; break;
                        case TBOX::HEXAHEDRON:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Hexahedron[8 * iElemHexahedron + iNode] = vnodes_local[iNode];
                            iElemHexahedron++; break;
                        case TBOX::PRISM:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Prism[6 * iElemPrism + iNode] = vnodes_local[iNode];
                            iElemPrism++; break;
                        case TBOX::PYRAMID:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Pyramid[5 * iElemPyramid + iNode] = vnodes_local[iNode];
                            iElemPyramid++; break;
                        }

                        iElemTotal++;
                    }

                    /*--- Set the value of the boundary geometry ---*/
                    iMarkerDomain = 0;
                    iBoundLineTotal = 0; iBoundTriangleTotal = 0; iBoundRectangleTotal = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if ((config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE) && (MarkerIn[iMarker]))
                        {
                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                if (VertexIn[iMarker][iVertex])
                                {
                                    for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++)
                                    {
                                        vnodes_local[iNode] = Global2Local_Point[geometry->bound[iMarker][iVertex]->GetNode(iNode)];
                                    }

                                    switch (geometry->bound[iMarker][iVertex]->GetVTK_Type())
                                    {
                                    case TBOX::LINE:
                                        Buffer_Send_BoundLine[2 * iBoundLineTotal + 0] = vnodes_local[0];
                                        Buffer_Send_BoundLine[2 * iBoundLineTotal + 1] = vnodes_local[1];
                                        iBoundLineTotal++;
                                        break;
                                    case TBOX::TRIANGLE:
                                        Buffer_Send_BoundTriangle[3 * iBoundTriangleTotal + 0] = vnodes_local[0];
                                        Buffer_Send_BoundTriangle[3 * iBoundTriangleTotal + 1] = vnodes_local[1];
                                        Buffer_Send_BoundTriangle[3 * iBoundTriangleTotal + 2] = vnodes_local[2];
                                        iBoundTriangleTotal++;
                                        break;
                                    case TBOX::RECTANGLE:
                                        Buffer_Send_BoundRectangle[4 * iBoundRectangleTotal + 0] = vnodes_local[0];
                                        Buffer_Send_BoundRectangle[4 * iBoundRectangleTotal + 1] = vnodes_local[1];
                                        Buffer_Send_BoundRectangle[4 * iBoundRectangleTotal + 2] = vnodes_local[2];
                                        Buffer_Send_BoundRectangle[4 * iBoundRectangleTotal + 3] = vnodes_local[3];
                                        iBoundRectangleTotal++;
                                        break;
                                    }
                                }
                            }

                            Buffer_Send_Local2Global_Marker[iMarkerDomain] = iMarker;
                            iMarkerDomain++;
                        }
                    }

                    /*--- Evaluate the number of already existing periodic boundary conditions ---*/
                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE)
                        {
                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
                                Transformation = geometry->bound[iMarker][iVertex]->GetRotation_Type();

                                if (iDomain == geometry->node[iPoint]->GetColor())
                                {
                                    /*--- If the information is going to be sended, find the
                                    domain of the receptor ---*/

                                    if (config->GetMarker_All_SendRecv(iMarker) > 0)
                                    {
                                        /*--- Identify the color of the receptor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                ReceptorColor = geometry->node[jPoint]->GetColor();
                                            }
                                        }

                                        /*--- For each color of the receptor we will han an extra marker (+) ---*/
                                        Buffer_Send_SendDomain_Periodic[iTotalSendDomain_Periodic] = Global2Local_Point[iPoint];
                                        Buffer_Send_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic] = Transformation;
                                        Buffer_Send_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] = ReceptorColor;

                                        iTotalSendDomain_Periodic++;
                                    }

                                    /*--- If the information is goint to be received, find the domain if the donor ---*/
                                    if (config->GetMarker_All_SendRecv(iMarker) < 0)
                                    {

                                        /*--- Identify the color of the donor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                DonorColor = geometry->node[jPoint]->GetColor();
                                            }
                                        }

                                        /*--- For each color of the donor we will han an extra marker (-) ---*/
                                        Buffer_Send_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic] = Global2Local_Point[iPoint];
                                        Buffer_Send_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic] = Transformation;
                                        Buffer_Send_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] = DonorColor;

                                        iTotalReceivedDomain_Periodic++;
                                    }
                                }
                            }
                        }
                    }

                    /*--- Send the buffers with the geometrical information ---*/
                    if (iDomain != TBOX::MASTER_NODE)
                    {

#ifdef HAVE_MPI
                        MPI_Isend(Buffer_Send_Coord, Buffer_Send_nPointTotal*Buffer_Send_nDim, MPI_DOUBLE, iDomain, 0, MPI_COMM_WORLD, &send_req[0]);
                        MPI_Isend(Buffer_Send_GlobalPointIndex, Buffer_Send_nPointTotal, MPI_UNSIGNED_LONG, iDomain, 1, MPI_COMM_WORLD, &send_req[1]);
                        MPI_Isend(Buffer_Send_Color, Buffer_Send_nPointTotal, MPI_UNSIGNED_LONG, iDomain, 2, MPI_COMM_WORLD, &send_req[2]);
                        MPI_Isend(Buffer_Send_Triangle, Buffer_Send_nElemTriangle * 3, MPI_UNSIGNED_LONG, iDomain, 3, MPI_COMM_WORLD, &send_req[3]);
                        MPI_Isend(Buffer_Send_Rectangle, Buffer_Send_nElemRectangle * 4, MPI_UNSIGNED_LONG, iDomain, 4, MPI_COMM_WORLD, &send_req[4]);
                        MPI_Isend(Buffer_Send_Tetrahedron, Buffer_Send_nElemTetrahedron * 4, MPI_UNSIGNED_LONG, iDomain, 5, MPI_COMM_WORLD, &send_req[5]);
                        MPI_Isend(Buffer_Send_Hexahedron, Buffer_Send_nElemHexahedron * 8, MPI_UNSIGNED_LONG, iDomain, 6, MPI_COMM_WORLD, &send_req[6]);
                        MPI_Isend(Buffer_Send_Prism, Buffer_Send_nElemPrism * 6, MPI_UNSIGNED_LONG, iDomain, 7, MPI_COMM_WORLD, &send_req[7]);
                        MPI_Isend(Buffer_Send_Pyramid, Buffer_Send_nElemPyramid * 5, MPI_UNSIGNED_LONG, iDomain, 8, MPI_COMM_WORLD, &send_req[8]);
                        MPI_Isend(Buffer_Send_BoundLine, Buffer_Send_nBoundLineTotal * 2, MPI_UNSIGNED_LONG, iDomain, 9, MPI_COMM_WORLD, &send_req[9]);
                        MPI_Isend(Buffer_Send_BoundTriangle, Buffer_Send_nBoundTriangleTotal * 3, MPI_UNSIGNED_LONG, iDomain, 10, MPI_COMM_WORLD, &send_req[10]);
                        MPI_Isend(Buffer_Send_BoundRectangle, Buffer_Send_nBoundRectangleTotal * 4, MPI_UNSIGNED_LONG, iDomain, 11, MPI_COMM_WORLD, &send_req[11]);
                        MPI_Isend(Buffer_Send_Local2Global_Marker, Buffer_Send_nMarkerDomain, MPI_UNSIGNED_LONG, iDomain, 12, MPI_COMM_WORLD, &send_req[12]);

                        MPI_Isend(Buffer_Send_SendDomain_Periodic, Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 13, MPI_COMM_WORLD, &send_req[13]);
                        MPI_Isend(Buffer_Send_SendDomain_PeriodicTrans, Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 14, MPI_COMM_WORLD, &send_req[14]);
                        MPI_Isend(Buffer_Send_SendDomain_PeriodicReceptor, Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 15, MPI_COMM_WORLD, &send_req[15]);
                        MPI_Isend(Buffer_Send_ReceivedDomain_Periodic, Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 16, MPI_COMM_WORLD, &send_req[16]);
                        MPI_Isend(Buffer_Send_ReceivedDomain_PeriodicTrans, Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 17, MPI_COMM_WORLD, &send_req[17]);
                        MPI_Isend(Buffer_Send_ReceivedDomain_PeriodicDonor, Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 18, MPI_COMM_WORLD, &send_req[18]);

                        /*--- Wait for this set of non-blocking comm. to complete ---*/

                        MPI_Waitall(19, send_req, send_stat);
#endif

                    }
                    else
                    {

                        for (iter = 0; iter < Buffer_Send_nPointTotal*Buffer_Send_nDim; iter++)
                            Buffer_Receive_Coord[iter] = Buffer_Send_Coord[iter];

                        for (iter = 0; iter < Buffer_Send_nPointTotal; iter++)
                        {
                            Buffer_Receive_GlobalPointIndex[iter] = Buffer_Send_GlobalPointIndex[iter];
                            Buffer_Receive_Color[iter] = Buffer_Send_Color[iter];
                        }

                        for (iter = 0; iter < Buffer_Send_nElemTriangle * 3; iter++)
                            Buffer_Receive_Triangle[iter] = Buffer_Send_Triangle[iter];

                        for (iter = 0; iter < Buffer_Send_nElemRectangle * 4; iter++)
                            Buffer_Receive_Rectangle[iter] = Buffer_Send_Rectangle[iter];

                        for (iter = 0; iter < Buffer_Send_nElemTetrahedron * 4; iter++)
                            Buffer_Receive_Tetrahedron[iter] = Buffer_Send_Tetrahedron[iter];

                        for (iter = 0; iter < Buffer_Send_nElemHexahedron * 8; iter++)
                            Buffer_Receive_Hexahedron[iter] = Buffer_Send_Hexahedron[iter];

                        for (iter = 0; iter < Buffer_Send_nElemPrism * 6; iter++)
                            Buffer_Receive_Prism[iter] = Buffer_Send_Prism[iter];

                        for (iter = 0; iter < Buffer_Send_nElemPyramid * 5; iter++)
                            Buffer_Receive_Pyramid[iter] = Buffer_Send_Pyramid[iter];

                        for (iter = 0; iter < Buffer_Send_nBoundLineTotal * 2; iter++)
                            Buffer_Receive_BoundLine[iter] = Buffer_Send_BoundLine[iter];

                        for (iter = 0; iter < Buffer_Send_nBoundTriangleTotal * 3; iter++)
                            Buffer_Receive_BoundTriangle[iter] = Buffer_Send_BoundTriangle[iter];

                        for (iter = 0; iter < Buffer_Send_nBoundRectangleTotal * 4; iter++)
                            Buffer_Receive_BoundRectangle[iter] = Buffer_Send_BoundRectangle[iter];

                        for (iter = 0; iter < Buffer_Send_nMarkerDomain; iter++)
                            Buffer_Receive_Local2Global_Marker[iter] = Buffer_Send_Local2Global_Marker[iter];

                        for (iter = 0; iter < Buffer_Send_nTotalSendDomain_Periodic; iter++)
                        {
                            Buffer_Receive_SendDomain_Periodic[iter] = Buffer_Send_SendDomain_Periodic[iter];
                            Buffer_Receive_SendDomain_PeriodicTrans[iter] = Buffer_Send_SendDomain_PeriodicTrans[iter];
                            Buffer_Receive_SendDomain_PeriodicReceptor[iter] = Buffer_Send_SendDomain_PeriodicReceptor[iter];
                        }

                        for (iter = 0; iter < Buffer_Send_nTotalReceivedDomain_Periodic; iter++)
                        {
                            Buffer_Receive_ReceivedDomain_Periodic[iter] = Buffer_Send_ReceivedDomain_Periodic[iter];
                            Buffer_Receive_ReceivedDomain_PeriodicTrans[iter] = Buffer_Send_ReceivedDomain_PeriodicTrans[iter];
                            Buffer_Receive_ReceivedDomain_PeriodicDonor[iter] = Buffer_Send_ReceivedDomain_PeriodicDonor[iter];
                        }

                    }

                    delete[] Buffer_Send_Coord;
                    delete[] Buffer_Send_GlobalPointIndex;
                    delete[] Buffer_Send_Color;
                    delete[] Buffer_Send_Triangle;
                    delete[] Buffer_Send_Rectangle;
                    delete[] Buffer_Send_Tetrahedron;
                    delete[] Buffer_Send_Hexahedron;
                    delete[] Buffer_Send_Prism;
                    delete[] Buffer_Send_Pyramid;
                    delete[] Buffer_Send_BoundLine;
                    delete[] Buffer_Send_BoundTriangle;
                    delete[] Buffer_Send_BoundRectangle;
                    delete[] Buffer_Send_Local2Global_Marker;

                    delete[] Buffer_Send_SendDomain_Periodic;
                    delete[] Buffer_Send_SendDomain_PeriodicTrans;
                    delete[] Buffer_Send_SendDomain_PeriodicReceptor;
                    delete[] Buffer_Send_ReceivedDomain_Periodic;
                    delete[] Buffer_Send_ReceivedDomain_PeriodicTrans;
                    delete[] Buffer_Send_ReceivedDomain_PeriodicDonor;
                }

                if (rank == iDomain)
                {
                    if (rank != TBOX::MASTER_NODE)
                    {
                        /*--- Receive the buffers with the geometrical information ---*/
#ifdef HAVE_MPI
                        MPI_Irecv(Buffer_Receive_Coord, nPointTotal*nDim, MPI_DOUBLE, MASTER_NODE, 0, MPI_COMM_WORLD, &recv_req[0]);
                        MPI_Irecv(Buffer_Receive_GlobalPointIndex, nPointTotal, MPI_UNSIGNED_LONG, MASTER_NODE, 1, MPI_COMM_WORLD, &recv_req[1]);
                        MPI_Irecv(Buffer_Receive_Color, nPointTotal, MPI_UNSIGNED_LONG, MASTER_NODE, 2, MPI_COMM_WORLD, &recv_req[2]);
                        MPI_Irecv(Buffer_Receive_Triangle, nElemTriangle * 3, MPI_UNSIGNED_LONG, MASTER_NODE, 3, MPI_COMM_WORLD, &recv_req[3]);
                        MPI_Irecv(Buffer_Receive_Rectangle, nElemRectangle * 4, MPI_UNSIGNED_LONG, MASTER_NODE, 4, MPI_COMM_WORLD, &recv_req[4]);
                        MPI_Irecv(Buffer_Receive_Tetrahedron, nElemTetrahedron * 4, MPI_UNSIGNED_LONG, MASTER_NODE, 5, MPI_COMM_WORLD, &recv_req[5]);
                        MPI_Irecv(Buffer_Receive_Hexahedron, nElemHexahedron * 8, MPI_UNSIGNED_LONG, MASTER_NODE, 6, MPI_COMM_WORLD, &recv_req[6]);
                        MPI_Irecv(Buffer_Receive_Prism, nElemPrism * 6, MPI_UNSIGNED_LONG, MASTER_NODE, 7, MPI_COMM_WORLD, &recv_req[7]);
                        MPI_Irecv(Buffer_Receive_Pyramid, nElemPyramid * 5, MPI_UNSIGNED_LONG, MASTER_NODE, 8, MPI_COMM_WORLD, &recv_req[8]);
                        MPI_Irecv(Buffer_Receive_BoundLine, nBoundLineTotal * 2, MPI_UNSIGNED_LONG, MASTER_NODE, 9, MPI_COMM_WORLD, &recv_req[9]);
                        MPI_Irecv(Buffer_Receive_BoundTriangle, nBoundTriangleTotal * 3, MPI_UNSIGNED_LONG, MASTER_NODE, 10, MPI_COMM_WORLD, &recv_req[10]);
                        MPI_Irecv(Buffer_Receive_BoundRectangle, nBoundRectangleTotal * 4, MPI_UNSIGNED_LONG, MASTER_NODE, 11, MPI_COMM_WORLD, &recv_req[11]);
                        MPI_Irecv(Buffer_Receive_Local2Global_Marker, nMarkerDomain, MPI_UNSIGNED_LONG, MASTER_NODE, 12, MPI_COMM_WORLD, &recv_req[12]);

                        MPI_Irecv(Buffer_Receive_SendDomain_Periodic, nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 13, MPI_COMM_WORLD, &recv_req[13]);
                        MPI_Irecv(Buffer_Receive_SendDomain_PeriodicTrans, nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 14, MPI_COMM_WORLD, &recv_req[14]);
                        MPI_Irecv(Buffer_Receive_SendDomain_PeriodicReceptor, nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 15, MPI_COMM_WORLD, &recv_req[15]);
                        MPI_Irecv(Buffer_Receive_ReceivedDomain_Periodic, nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 16, MPI_COMM_WORLD, &recv_req[16]);
                        MPI_Irecv(Buffer_Receive_ReceivedDomain_PeriodicTrans, nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 17, MPI_COMM_WORLD, &recv_req[17]);
                        MPI_Irecv(Buffer_Receive_ReceivedDomain_PeriodicDonor, nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 18, MPI_COMM_WORLD, &recv_req[18]);

                        /*--- Wait for this set of non-blocking recv's to complete ---*/

                        MPI_Waitall(19, recv_req, recv_stat);
#endif
                    }

                    /*--- Create the domain structures for the points ---*/
                    nPoint = nPointTotal;
                    nPointDomain = nPointDomainTotal;
                    node = new GRID::GRID_DGPoint*[nPoint];
                    Local_to_Global_Point = new long[nPoint];

                    for (iPoint = 0; iPoint < nPoint; iPoint++)
                    {
                        Local_to_Global_Point[iPoint] = Buffer_Receive_GlobalPointIndex[iPoint];
                        if (nDim == 2)
                            node[iPoint] = new GRID::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0], Buffer_Receive_Coord[iPoint*nDim + 1], Local_to_Global_Point[iPoint], config);
                        if (nDim == 3)
                            node[iPoint] = new GRID::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0], Buffer_Receive_Coord[iPoint*nDim + 1], Buffer_Receive_Coord[iPoint*nDim + 2], Local_to_Global_Point[iPoint], config);
                        node[iPoint]->SetColor(Buffer_Receive_Color[iPoint]);
                    }

                    delete[] Buffer_Receive_Coord;
                    delete[] Buffer_Receive_GlobalPointIndex;
                    delete[] Buffer_Receive_Color;

                    /*--- Create the domain structures for the elements ---*/
                    nElem = nElemTotal; iElem = 0;
                    elem = new GRID::GRID_Primal*[nElem];

                    for (iElemTriangle = 0; iElemTriangle < nElemTriangle; iElemTriangle++)
                    {
                        elem[iElem] = new GRID::GRID_Triangle(Buffer_Receive_Triangle[iElemTriangle * 3 + 0], Buffer_Receive_Triangle[iElemTriangle * 3 + 1], Buffer_Receive_Triangle[iElemTriangle * 3 + 2], 2);
                        iElem++;
                    }
                    for (iElemRectangle = 0; iElemRectangle < nElemRectangle; iElemRectangle++)
                    {
                        elem[iElem] = new GRID::GRID_Rectangle(Buffer_Receive_Rectangle[iElemRectangle * 4 + 0], Buffer_Receive_Rectangle[iElemRectangle * 4 + 1], Buffer_Receive_Rectangle[iElemRectangle * 4 + 2], Buffer_Receive_Rectangle[iElemRectangle * 4 + 3], 2);
                        iElem++;
                    }
                    for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron; iElemTetrahedron++)
                    {
                        elem[iElem] = new GRID::GRID_Tetrahedron(Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 0], Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 1], Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 2], Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 3]);
                        iElem++;
                    }
                    for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron; iElemHexahedron++)
                    {
                        elem[iElem] = new GRID::GRID_Hexahedron(Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 0], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 1], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 2], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 3], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 4], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 5], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 6], Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 7]);
                        iElem++;
                    }
                    for (iElemPrism = 0; iElemPrism < nElemPrism; iElemPrism++)
                    {
                        elem[iElem] = new GRID::GRID_Prism(Buffer_Receive_Prism[iElemPrism * 6 + 0], Buffer_Receive_Prism[iElemPrism * 6 + 1], Buffer_Receive_Prism[iElemPrism * 6 + 2], Buffer_Receive_Prism[iElemPrism * 6 + 3], Buffer_Receive_Prism[iElemPrism * 6 + 4], Buffer_Receive_Prism[iElemPrism * 6 + 5]);
                        iElem++;
                    }
                    for (iElemPyramid = 0; iElemPyramid < nElemPyramid; iElemPyramid++)
                    {
                        elem[iElem] = new GRID::GRID_Pyramid(Buffer_Receive_Pyramid[iElemPyramid * 5 + 0], Buffer_Receive_Pyramid[iElemPyramid * 5 + 1], Buffer_Receive_Pyramid[iElemPyramid * 5 + 2], Buffer_Receive_Pyramid[iElemPyramid * 5 + 3], Buffer_Receive_Pyramid[iElemPyramid * 5 + 4]);
                        iElem++;
                    }

                    delete[] Buffer_Receive_Triangle;
                    delete[] Buffer_Receive_Rectangle;
                    delete[] Buffer_Receive_Tetrahedron;
                    delete[] Buffer_Receive_Hexahedron;
                    delete[] Buffer_Receive_Prism;
                    delete[] Buffer_Receive_Pyramid;

                    /*--- Create the domain structures for the boundaries ---*/

                    nMarker = nMarkerDomain;

                    nElem_Bound = new unsigned long[nMarker_Max];
                    Local_to_Global_Marker = new unsigned short[nMarker_Max];
                    Tag_to_Marker = new std::string[nMarker_Max];
                    std::string *TagBound_Copy = new std::string[nMarker_Max];
                    short *SendRecv_Copy = new short[nMarker_Max];

                    for (iMarker = 0; iMarker < nMarker; iMarker++) nElem_Bound[iMarker] = nVertexDomain[iMarker];

                    bound = new GRID::GRID_Primal**[nMarker + (overhead*nDomain)];
                    for (iMarker = 0; iMarker < nMarker; iMarker++) bound[iMarker] = new GRID::GRID_Primal*[nElem_Bound[iMarker]];

                    iBoundLineTotal = 0; iBoundTriangleTotal = 0; iBoundRectangleTotal = 0;

                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                    {

                        iVertexDomain = 0;

                        for (iBoundLine = 0; iBoundLine < nBoundLine[iMarker]; iBoundLine++)
                        {
                            bound[iMarker][iVertexDomain] = new GRID::GRID_Line(Buffer_Receive_BoundLine[iBoundLineTotal * 2 + 0],
                                Buffer_Receive_BoundLine[iBoundLineTotal * 2 + 1], 2);
                            iVertexDomain++; iBoundLineTotal++;
                        }
                        for (iBoundTriangle = 0; iBoundTriangle < nBoundTriangle[iMarker]; iBoundTriangle++)
                        {
                            bound[iMarker][iVertexDomain] = new GRID::GRID_Triangle(Buffer_Receive_BoundTriangle[iBoundTriangleTotal * 3 + 0],
                                Buffer_Receive_BoundTriangle[iBoundTriangleTotal * 3 + 1],
                                Buffer_Receive_BoundTriangle[iBoundTriangleTotal * 3 + 2], 3);
                            iVertexDomain++; iBoundTriangleTotal++;
                        }
                        for (iBoundRectangle = 0; iBoundRectangle < nBoundRectangle[iMarker]; iBoundRectangle++)
                        {
                            bound[iMarker][iVertexDomain] = new GRID::GRID_Rectangle(Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 0],
                                Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 1],
                                Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 2],
                                Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 3], 3);
                            iVertexDomain++; iBoundRectangleTotal++;
                        }

                        Local_to_Global_Marker[iMarker] = Buffer_Receive_Local2Global_Marker[iMarker];

                        /*--- Now each domain has the right information ---*/
                        std::string Grid_Marker = config->GetMarker_All_TagBound(Local_to_Global_Marker[iMarker]);
                        short SendRecv = config->GetMarker_All_SendRecv(Local_to_Global_Marker[iMarker]);
                        TagBound_Copy[iMarker] = Grid_Marker;
                        SendRecv_Copy[iMarker] = SendRecv;
                    }

                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                    {
                        config->SetMarker_All_TagBound(iMarker, TagBound_Copy[iMarker]);
                        config->SetMarker_All_SendRecv(iMarker, SendRecv_Copy[iMarker]);
                    }

                    /*--- Add the new periodic markers to the domain ---*/
                    iTotalSendDomain_Periodic = 0;
                    iTotalReceivedDomain_Periodic = 0;

                    for (jDomain = 0; jDomain < nDomain; jDomain++)
                    {
                        if (nSendDomain_Periodic[jDomain] != 0)
                        {
                            nVertexDomain[nMarker] = 0;
                            bound[nMarker] = new GRID::GRID_Primal*[nSendDomain_Periodic[jDomain]];

                            iVertex = 0;
                            for (iTotalSendDomain_Periodic = 0; iTotalSendDomain_Periodic < nTotalSendDomain_Periodic; iTotalSendDomain_Periodic++)
                            {
                                if (Buffer_Receive_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] == jDomain)
                                {
                                    bound[nMarker][iVertex] = new GRID::GRID_VertexMPI(Buffer_Receive_SendDomain_Periodic[iTotalSendDomain_Periodic], nDim);
                                    bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic]);
                                    nVertexDomain[nMarker]++; iVertex++;
                                }
                            }

                            Marker_All_SendRecv[nMarker] = jDomain + 1;
                            nElem_Bound[nMarker] = nVertexDomain[nMarker];
                            nMarker++;
                        }

                        if (nReceivedDomain_Periodic[jDomain] != 0)
                        {
                            nVertexDomain[nMarker] = 0;
                            bound[nMarker] = new GRID::GRID_Primal*[nReceivedDomain_Periodic[jDomain]];

                            iVertex = 0;
                            for (iTotalReceivedDomain_Periodic = 0; iTotalReceivedDomain_Periodic < nTotalReceivedDomain_Periodic; iTotalReceivedDomain_Periodic++)
                            {
                                if (Buffer_Receive_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] == jDomain)
                                {
                                    bound[nMarker][iVertex] = new GRID::GRID_VertexMPI(Buffer_Receive_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic], nDim);
                                    bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic]);
                                    nVertexDomain[nMarker]++; iVertex++;
                                }
                            }

                            Marker_All_SendRecv[nMarker] = -(jDomain + 1);
                            nElem_Bound[nMarker] = nVertexDomain[nMarker];
                            nMarker++;
                        }

                    }

                    delete[] TagBound_Copy;
                    delete[] SendRecv_Copy;

                    delete[] Buffer_Receive_BoundLine;
                    delete[] Buffer_Receive_BoundTriangle;
                    delete[] Buffer_Receive_BoundRectangle;
                    delete[] Buffer_Receive_Local2Global_Marker;

                    delete[] Buffer_Receive_SendDomain_Periodic;
                    delete[] Buffer_Receive_SendDomain_PeriodicTrans;
                    delete[] Buffer_Receive_SendDomain_PeriodicReceptor;
                    delete[] Buffer_Receive_ReceivedDomain_Periodic;
                    delete[] Buffer_Receive_ReceivedDomain_PeriodicTrans;
                    delete[] Buffer_Receive_ReceivedDomain_PeriodicDonor;
                }

            }


            /*--- Set the value of Marker_All_SendRecv and Marker_All_TagBound in the config structure ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
            }

            /*--- Set the value of Global_nPoint and Global_nPointDomain ---*/
            unsigned long Local_nPoint = nPoint;
            unsigned long Local_nPointDomain = nPointDomain;

#ifdef HAVE_MPI

            MPI_Allreduce(&Local_nPoint, &Global_nPoint, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

#else

            Global_nPoint = Local_nPoint;
            Global_nPointDomain = Local_nPointDomain;

#endif

            if (rank == TBOX::MASTER_NODE)
            {
                delete[] MarkerIn;
                delete[] nElem_Color;

                delete[] Buffer_Send_Center;
                delete[] Buffer_Send_Rotation;
                delete[] Buffer_Send_Translate;

                delete[] Buffer_Send_nSendDomain_Periodic;
                delete[] Buffer_Send_nReceivedDomain_Periodic;

                delete[] Marker_All_SendRecv_Copy;
                delete[] Marker_All_TagBound_Copy;

                for (iDomain = 0; iDomain < nDomain; iDomain++)
                {
                    delete[] Elem_Color[iDomain];
                }
                delete[] Elem_Color;
                delete[] Global2Local_Point;

                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    delete VertexIn[iMarker];
                delete[] VertexIn;

            }

            delete[] nSendDomain_Periodic;
            delete[] nReceivedDomain_Periodic;

            delete[] nVertexDomain;
            delete[] nBoundLine;
            delete[] nBoundTriangle;
            delete[] nBoundRectangle;
            delete[] Buffer_Send_nVertexDomain;
            delete[] Buffer_Send_nBoundLine;
            delete[] Buffer_Send_nBoundTriangle;
            delete[] Buffer_Send_nBoundRectangle;
            delete[] Buffer_Send_Marker_All_SendRecv;
            delete[] Marker_All_TagBound;
            delete[] Buffer_Send_Marker_All_TagBound;
        }

        GEOM_GeometryPhysical::GEOM_GeometryPhysical(GEOM_Geometry *geometry, TBOX::TBOX_Config *config, int option)
        {
            Global_to_Local_Point = NULL;
            Local_to_Global_Point = NULL;
            Local_to_Global_Marker = NULL;
            Global_to_Local_Marker = NULL;

            unsigned long iter, iPoint, jPoint, iElem, iVertex;
            unsigned long nElemTotal = 0, nPointTotal = 0, nPointDomainTotal = 0, nPointGhost = 0, nPointPeriodic = 0, nElemTriangle = 0, nElemRectangle = 0, nElemTetrahedron = 0, nElemHexahedron = 0, nElemPrism = 0, nElemPyramid = 0;
            unsigned long iElemTotal, iPointTotal, iPointGhost, iPointDomain, iPointPeriodic, iElemTriangle, iElemRectangle, iElemTetrahedron, iElemHexahedron, iElemPrism, iElemPyramid, iPointCurrent;
            unsigned long nBoundLineTotal = 0, iBoundLineTotal;
            unsigned long nBoundTriangleTotal = 0, iBoundTriangleTotal;
            unsigned long nBoundRectangleTotal = 0, iBoundRectangleTotal;
            unsigned long ReceptorColor = 0, DonorColor = 0, Transformation;
            unsigned long nTotalSendDomain_Periodic = 0, iTotalSendDomain_Periodic, nTotalReceivedDomain_Periodic = 0, iTotalReceivedDomain_Periodic, *nSendDomain_Periodic = NULL, *nReceivedDomain_Periodic = NULL;
            unsigned long Buffer_Send_nPointTotal = 0, Buffer_Send_nPointDomainTotal = 0, Buffer_Send_nPointGhost = 0, Buffer_Send_nPointPeriodic = 0;
            unsigned long Buffer_Send_nElemTotal, Buffer_Send_nElemTriangle = 0, Buffer_Send_nElemRectangle = 0, Buffer_Send_nElemTetrahedron = 0, Buffer_Send_nElemHexahedron = 0, Buffer_Send_nElemPrism = 0, Buffer_Send_nElemPyramid = 0;
            unsigned long Buffer_Send_nTotalSendDomain_Periodic = 0, Buffer_Send_nTotalReceivedDomain_Periodic = 0, *Buffer_Send_nSendDomain_Periodic = NULL, *Buffer_Send_nReceivedDomain_Periodic = NULL;
            unsigned long Buffer_Send_nBoundLineTotal = 0, Buffer_Send_nBoundTriangleTotal = 0, Buffer_Send_nBoundRectangleTotal = 0;
            unsigned long iVertexDomain, iBoundLine, iBoundTriangle, iBoundRectangle;

            /*--- Need to double-check these shorts in case we go to nprocs > ~32,000 ---*/
            unsigned long iNode, iDim, iMarker, jMarker, nMarkerDomain = 0, iMarkerDomain;
            unsigned long nDomain = 0, iDomain, jDomain, nPeriodic = 0, iPeriodic, overhead = 4, Buffer_Send_nMarkerDomain = 0, Buffer_Send_nDim = 0, Buffer_Send_nZone = 0, Buffer_Send_nPeriodic = 0;

            bool *MarkerIn = NULL, **VertexIn = NULL, *PointIn = NULL, *ElemIn = NULL;
            long vnodes_local[8];

            std::vector<long> DomainList;
            short *Marker_All_SendRecv_Copy = NULL;
            std::string *Marker_All_TagBound_Copy = NULL;

            int rank = TBOX::MASTER_NODE;
            int size = TBOX::SINGLE_NODE;
            unsigned short nMarker_Max = config->GetnMarker_Max();

            /*--- Some dynamic arrays so we're not allocating too much on the stack ---*/
            unsigned long *nVertexDomain = new unsigned long[nMarker_Max];
            unsigned long *nBoundLine = new unsigned long[nMarker_Max];
            unsigned long *nBoundTriangle = new unsigned long[nMarker_Max];
            unsigned long *nBoundRectangle = new unsigned long[nMarker_Max];

            unsigned long *Buffer_Send_nVertexDomain = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nBoundLine = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nBoundTriangle = new unsigned long[nMarker_Max];
            unsigned long *Buffer_Send_nBoundRectangle = new unsigned long[nMarker_Max];

            short *Buffer_Send_Marker_All_SendRecv = new short[nMarker_Max];

            char *Marker_All_TagBound = new char[nMarker_Max*TBOX::MAX_STRING_SIZE];
            char *Buffer_Send_Marker_All_TagBound = new char[nMarker_Max*TBOX::MAX_STRING_SIZE];

#ifdef HAVE_MPI

            /*--- MPI initialization ---*/
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            /*--- MPI status and request arrays for non-blocking communications ---*/
            MPI_Status status, status2;
            unsigned long source;
            int recv_count = 0;

            int offset = 17;
            MPI_Status *send_stat = new MPI_Status[offset + size];
            MPI_Status *recv_stat = new MPI_Status[offset + size];

            MPI_Request *send_req = new MPI_Request[offset + size];
            MPI_Request *recv_req = new MPI_Request[offset + size];

#endif

            if (rank == TBOX::MASTER_NODE && size > TBOX::SINGLE_NODE)
                std::cout << "Communicating partition data and creating halo layers." << std::endl;

            /*--- Define buffer std::vector interior domain ---*/
            double        *Buffer_Send_Coord = NULL;
            unsigned long *Buffer_Send_Color = NULL;
            unsigned long *Buffer_Send_GlobalPointIndex = NULL;
            unsigned long *Buffer_Send_Triangle = NULL;
            unsigned long *Buffer_Send_Rectangle = NULL;
            unsigned long *Buffer_Send_Tetrahedron = NULL;
            unsigned long *Buffer_Send_Hexahedron = NULL;
            unsigned long *Buffer_Send_Prism = NULL;
            unsigned long *Buffer_Send_Pyramid = NULL;
            unsigned long *Buffer_Send_GlobElem = NULL;

            /*--- Define buffer std::vector boundary ---*/
            unsigned long *Buffer_Send_BoundLine = NULL, *Buffer_Receive_BoundLine = NULL;
            unsigned long *Buffer_Send_BoundTriangle = NULL, *Buffer_Receive_BoundTriangle = NULL;
            unsigned long *Buffer_Send_BoundRectangle = NULL, *Buffer_Receive_BoundRectangle = NULL;
            unsigned long *Buffer_Send_Local2Global_Marker = NULL, *Buffer_Receive_Local2Global_Marker = NULL;

            /*--- Define buffer std::vector periodic boundary conditions ---*/
            double *Buffer_Send_Center = NULL, *Buffer_Receive_Center = NULL;
            double *Buffer_Send_Rotation = NULL, *Buffer_Receive_Rotation = NULL;
            double *Buffer_Send_Translate = NULL, *Buffer_Receive_Translate = NULL;

            /*--- Define buffer std::vector periodic boundary conditions ---*/
            unsigned long *Buffer_Send_SendDomain_Periodic = NULL, *Buffer_Receive_SendDomain_Periodic = NULL;
            unsigned long *Buffer_Send_SendDomain_PeriodicTrans = NULL, *Buffer_Receive_SendDomain_PeriodicTrans = NULL;
            unsigned long *Buffer_Send_SendDomain_PeriodicReceptor = NULL, *Buffer_Receive_SendDomain_PeriodicReceptor = NULL;
            unsigned long *Buffer_Send_ReceivedDomain_Periodic = NULL, *Buffer_Receive_ReceivedDomain_Periodic = NULL;
            unsigned long *Buffer_Send_ReceivedDomain_PeriodicTrans = NULL, *Buffer_Receive_ReceivedDomain_PeriodicTrans = NULL;
            unsigned long *Buffer_Send_ReceivedDomain_PeriodicDonor = NULL, *Buffer_Receive_ReceivedDomain_PeriodicDonor = NULL;

            /*--- Variables below are needed specifically for the ParMETIS version ---*/
            unsigned long *Global_to_local_Point_recv;
            unsigned long *local_colour_values;
            unsigned long *local_colour_temp;
            unsigned long *Local_to_global_elem;

            unsigned short *nDim_s = new unsigned short[size];
            unsigned short *nDim_r = new unsigned short[size];
            unsigned short *nZone_s = new unsigned short[size];
            unsigned short *nZone_r = new unsigned short[size];

            unsigned long *nPointTotal_s = new unsigned long[size];
            unsigned long *nPointDomainTotal_s = new unsigned long[size];
            unsigned long *nPointGhost_s = new unsigned long[size];
            unsigned long *nPointPeriodic_s = new unsigned long[size];
            unsigned long *nElemTotal_s = new unsigned long[size];
            unsigned long *nElemTriangle_s = new unsigned long[size];
            unsigned long *nElemRectangle_s = new unsigned long[size];
            unsigned long *nElemTetrahedron_s = new unsigned long[size];
            unsigned long *nElemHexahedron_s = new unsigned long[size];
            unsigned long *nElemPrism_s = new unsigned long[size];
            unsigned long *nElemPyramid_s = new unsigned long[size];

            unsigned long *nPointTotal_r = new unsigned long[size];
            unsigned long *nPointDomainTotal_r = new unsigned long[size];
            unsigned long *nPointGhost_r = new unsigned long[size];
            unsigned long *nPointPeriodic_r = new unsigned long[size];
            unsigned long *nElemTotal_r = new unsigned long[size];
            unsigned long *nElemTriangle_r = new unsigned long[size];
            unsigned long *nElemRectangle_r = new unsigned long[size];
            unsigned long *nElemTetrahedron_r = new unsigned long[size];
            unsigned long *nElemHexahedron_r = new unsigned long[size];
            unsigned long *nElemPrism_r = new unsigned long[size];
            unsigned long *nElemPyramid_r = new unsigned long[size];

            unsigned long nPointTotal_r_tot = 0;
            unsigned long nPointDomainTotal_r_tot = 0;
            unsigned long nPointGhost_r_tot = 0;
            unsigned long nPointPeriodic_r_tot = 0;
            unsigned long nElemTotal_r_tot = 0;
            unsigned long nElemTriangle_r_tot = 0;
            unsigned long nElemRectangle_r_tot = 0;
            unsigned long nElemTetrahedron_r_tot = 0;
            unsigned long nElemHexahedron_r_tot = 0;
            unsigned long nElemPrism_r_tot = 0;
            unsigned long nElemPyramid_r_tot = 0;

            unsigned long Buffer_Size_Coord = 0;
            unsigned long Buffer_Size_Color = 0;
            unsigned long Buffer_Size_GlobalPointIndex = 0;
            unsigned long Buffer_Size_Triangle = 0;
            unsigned long Buffer_Size_Rectangle = 0;
            unsigned long Buffer_Size_Tetrahedron = 0;
            unsigned long Buffer_Size_Hexahedron = 0;
            unsigned long Buffer_Size_Prism = 0;
            unsigned long Buffer_Size_Pyramid = 0;
            unsigned long Buffer_Size_GlobElem = 0;

            unsigned long ElemTotal_Counter = 0;
            unsigned long PointTotal_Counter = 0;
            unsigned long PointDomain_Counter = 0;

            /*--- WARNING: check the next two counters ---*/
            unsigned long PointPeriodic_Counter = 0;
            unsigned long PointGhost_Counter = 0;
            unsigned long ElemTriangle_Counter = 0;
            unsigned long ElemRectangle_Counter = 0;
            unsigned long ElemTetrahedron_Counter = 0;
            unsigned long ElemHexahedron_Counter = 0;
            unsigned long ElemPrism_Counter = 0;
            unsigned long ElemPyramid_Counter = 0;

            unsigned long *Local_to_global_Triangle;
            unsigned long *Local_to_global_Rectangle;
            unsigned long *Local_to_global_Tetrahedron;
            unsigned long *Local_to_global_Hexahedron;
            unsigned long *Local_to_global_Prism;
            unsigned long *Local_to_global_Pyramid;

            bool *Triangle_presence;
            bool *Rectangle_presence;
            bool *Tetrahedron_presence;
            bool *Hexahedron_presence;
            bool *Prism_presence;
            bool *Pyramid_presence;
            bool *Element_presence;

            Element_presence = new bool[geometry->GetnElem()];
            Triangle_presence = new bool[geometry->GetnElem()];
            Rectangle_presence = new bool[geometry->GetnElem()];
            Tetrahedron_presence = new bool[geometry->GetnElem()];
            Hexahedron_presence = new bool[geometry->GetnElem()];
            Prism_presence = new bool[geometry->GetnElem()];
            Pyramid_presence = new bool[geometry->GetnElem()];

            for (unsigned long i = 0; i < geometry->GetnElem(); i++)
            {
                Element_presence[i] = false;
                Triangle_presence[i] = false;
                Rectangle_presence[i] = false;
                Tetrahedron_presence[i] = false;
                Hexahedron_presence[i] = false;
                Prism_presence[i] = false;
                Pyramid_presence[i] = false;
            }

            double *Buffer_Receive_Coord_loc = NULL;

            unsigned long *Buffer_Receive_Color_loc = NULL;
            unsigned long *Buffer_Receive_GlobalPointIndex_loc = NULL;
            unsigned long *Buffer_Receive_Triangle_loc = NULL;
            unsigned long *Buffer_Receive_Rectangle_loc = NULL;
            unsigned long *Buffer_Receive_Tetrahedron_loc = NULL;
            unsigned long *Buffer_Receive_Hexahedron_loc = NULL;
            unsigned long *Buffer_Receive_Prism_loc = NULL;
            unsigned long *Buffer_Receive_Pyramid_loc = NULL;

            unsigned long *Buffer_Receive_GlobElem_loc = NULL;
            unsigned long *Buffer_Receive_Triangle_presence_loc = NULL;
            unsigned long *Buffer_Receive_Rectangle_presence_loc = NULL;
            unsigned long *Buffer_Receive_Tetrahedron_presence_loc = NULL;
            unsigned long *Buffer_Receive_Hexahedron_presence_loc = NULL;
            unsigned long *Buffer_Receive_Prism_presence_loc = NULL;
            unsigned long *Buffer_Receive_Pyramid_presence_loc = NULL;

            /*--- Allocate the memory that we only need if we have MPI support ---*/

#ifdef HAVE_MPI

            unsigned long temp_element_count = 0;

            double        *Buffer_Receive_Coord = NULL;
            unsigned long *Buffer_Receive_Color = NULL;
            unsigned long *Buffer_Receive_GlobalPointIndex = NULL;
            unsigned long *Buffer_Receive_Triangle = NULL;
            unsigned long *Buffer_Receive_Rectangle = NULL;
            unsigned long *Buffer_Receive_Tetrahedron = NULL;
            unsigned long *Buffer_Receive_Hexahedron = NULL;
            unsigned long *Buffer_Receive_Prism = NULL;
            unsigned long *Buffer_Receive_Pyramid = NULL;
            unsigned long *Buffer_Receive_GlobElem = NULL;

            unsigned long **Buffer_Receive_Triangle_presence = new unsigned long*[size];
            unsigned long **Buffer_Receive_Rectangle_presence = new unsigned long*[size];
            unsigned long **Buffer_Receive_Tetrahedron_presence = new unsigned long*[size];
            unsigned long **Buffer_Receive_Hexahedron_presence = new unsigned long*[size];
            unsigned long **Buffer_Receive_Prism_presence = new unsigned long*[size];
            unsigned long **Buffer_Receive_Pyramid_presence = new unsigned long*[size];

#endif

            /*--- Basic dimensionalization ---*/
            nDomain = size;

            Marker_All_SendRecv = new short[nMarker_Max];
            nSendDomain_Periodic = new unsigned long[nDomain];
            nReceivedDomain_Periodic = new unsigned long[nDomain];

            /*--- Auxiliar std::vector based on the original geometry ---*/
            ElemIn = new bool[geometry->no_of_local_elements];
            PointIn = new bool[geometry->GetnPoint()];


            Buffer_Send_nDim = geometry->GetnDim();
            Buffer_Send_nZone = geometry->GetnZone();

            // DOUBLE CHECK THESE, SINCE WE DO THIS AGAIN AT BOTTOM WITH THE MASTER
            //  MarkerIn = new bool [geometry->GetnMarker()];
            //  VertexIn = new bool* [geometry->GetnMarker()];
            //  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
            //    VertexIn[iMarker] = new bool [geometry->GetnElem_Bound(iMarker)];
            //
            //  
            //  Buffer_Send_nPeriodic = config->GetnPeriodicIndex();
            //  Buffer_Send_Center    = new double[Buffer_Send_nPeriodic*3];
            //  Buffer_Send_Rotation  = new double[Buffer_Send_nPeriodic*3];
            //  Buffer_Send_Translate = new double[Buffer_Send_nPeriodic*3];
            //  
            //  Buffer_Send_nSendDomain_Periodic = new unsigned long [nDomain];
            //  Buffer_Send_nReceivedDomain_Periodic = new unsigned long [nDomain];

            /*--- Divide the elements in color list to speed up the grid partitioning ---*/

            Local_to_global_elem = new unsigned long[geometry->no_of_local_elements];
            for (unsigned long i = 0; i < geometry->GetnElem(); i++)
            {
                if (geometry->Global_to_local_elem[i] != -1)
                {
                    Local_to_global_elem[geometry->Global_to_local_elem[i]] = i;
                }
            }

            Global_to_local_Point_recv = new unsigned long[geometry->GetnPoint()];
            for (unsigned long i = 0; i < geometry->GetnPoint(); i++)
            {
                Global_to_local_Point_recv[i] = -1;
            }

            //  unsigned long *Global_to_Local_Point_loc;
            //  Global_to_Local_Point_loc = new unsigned long[geometry->GetnPoint()];
            //  for (iPoint=0; iPoint<geometry->GetnPoint(); iPoint++) {
            //    Global_to_Local_Point_loc[iPoint]=-1;
            //  }

            local_colour_values = new unsigned long[geometry->GetnPoint()];
            local_colour_temp = new unsigned long[geometry->ending_node[rank] - geometry->starting_node[rank]];

            for (unsigned long i = 0; i < geometry->ending_node[rank] - geometry->starting_node[rank]; i++)
            {
                local_colour_temp[i] = geometry->node[i]->GetColor();
                local_colour_values[geometry->starting_node[rank] + i] = local_colour_temp[i];
            }

            /*--- Communicate the grid coloring to all partitions. This information
            will be repeatedly used throughout the organization of the partitions
            and sorting out their ghost points/elements. ---*/

#ifdef HAVE_MPI

            int comm_counter = 0;
            for (iDomain = 0; iDomain < size; iDomain++) {
                if (iDomain != rank) {
                    MPI_Isend(local_colour_temp, geometry->ending_node[rank] - geometry->starting_node[rank],
                        MPI_UNSIGNED_LONG, iDomain, iDomain, MPI_COMM_WORLD, &send_req[comm_counter]);
                    comm_counter++;
                }
            }

            unsigned long*recv_buffer;
            for (iDomain = 0; iDomain < size - 1; iDomain++) {
                MPI_Probe(MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status2);
                source = status2.MPI_SOURCE;
                MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                MPI_Recv(&local_colour_values[geometry->starting_node[source]], recv_count,
                    MPI_UNSIGNED_LONG, source, rank, MPI_COMM_WORLD, &status2);
            }

            /*--- Wait for the sends to complete (will be true since we're using
            blocking recv's above. ---*/

            MPI_Waitall(size - 1, send_req, send_stat);

#endif

            /*--- Free temporary buffer for communicating colors. ---*/
            delete[] local_colour_temp;

#ifdef HAVE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            //cout << " ==== Rank " << rank << " starting first send " << endl;

            /*--- This loop gets the array sizes of points, elements, etc. for each
            rank to send to each other rank. ---*/

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                /*--- Interior dimensionalization. Loop over the original grid to
                perform the dimensionalizaton of the domain variables ---*/

                Buffer_Send_nElemTotal = 0;
                Buffer_Send_nPointTotal = 0;
                Buffer_Send_nPointGhost = 0;
                Buffer_Send_nPointDomainTotal = 0;
                Buffer_Send_nPointPeriodic = 0;
                Buffer_Send_nElemTriangle = 0;
                Buffer_Send_nElemRectangle = 0;
                Buffer_Send_nElemTetrahedron = 0;
                Buffer_Send_nElemHexahedron = 0;
                Buffer_Send_nElemPrism = 0;
                Buffer_Send_nElemPyramid = 0;

                /*--- Initialize the global to local mapping ---*/
                for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                {
                    PointIn[iPoint] = false;
                }

                /*--- Loop over all of the local elements and count the number of each
                type of point and element that needs to be sent. ---*/

                for (iElem = 0; iElem < geometry->no_of_local_elements; iElem++)
                {
                    /*--- Check if the element belongs to the domain ---*/
                    ElemIn[iElem] = false;
                    for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = geometry->elem[iElem]->GetNode(iNode);
                        if (local_colour_values[iPoint] == iDomain)
                        {
                            ElemIn[iElem] = true; break;
                        }
                    }

                    /*--- If this element is needed by iDomain, get information
                    about the number of points and element type. ---*/
                    if (ElemIn[iElem])
                    {
                        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                        {
                            iPoint = geometry->elem[iElem]->GetNode(iNode);

                            /*--- If we haven't already found this point... ---*/

                            if (PointIn[iPoint] == false)
                            {
                                /*--- Mark point as found and collect information ---*/
                                PointIn[iPoint] = true;

                                if ((iPoint >= geometry->starting_node[rank]) &&
                                    (iPoint < geometry->ending_node[rank]))
                                {
                                    Buffer_Send_nPointTotal++;

                                    /*--- Increment our counters ---*/
                                    if (local_colour_values[iPoint] == iDomain)
                                    {
                                        if (iPoint > geometry->GetnPointDomain() - 1)
                                            Buffer_Send_nPointPeriodic++;
                                        else
                                            Buffer_Send_nPointDomainTotal++;
                                    }
                                    else
                                        Buffer_Send_nPointGhost++;
                                }
                            }
                        }

                        /*--- Increment the counter for the current type of element ---*/
                        switch (geometry->elem[iElem]->GetVTK_Type())
                        {
                        case TBOX::TRIANGLE:    Buffer_Send_nElemTriangle++;    break;
                        case TBOX::RECTANGLE:   Buffer_Send_nElemRectangle++;   break;
                        case TBOX::TETRAHEDRON: Buffer_Send_nElemTetrahedron++; break;
                        case TBOX::HEXAHEDRON:  Buffer_Send_nElemHexahedron++;  break;
                        case TBOX::PRISM:       Buffer_Send_nElemPrism++;       break;
                        case TBOX::PYRAMID:     Buffer_Send_nElemPyramid++;     break;
                        }

                        /*--- Increment the total number of elements for iDomain ---*/
                        Buffer_Send_nElemTotal++;
                    }
                }

                /*--- Store the counts on a partition by partition basis. ---*/
                nDim_s[iDomain] = geometry->GetnDim();
                nZone_s[iDomain] = Buffer_Send_nZone;
                nPointTotal_s[iDomain] = Buffer_Send_nPointTotal;
                nPointDomainTotal_s[iDomain] = Buffer_Send_nPointDomainTotal;
                nPointGhost_s[iDomain] = Buffer_Send_nPointGhost;
                nPointPeriodic_s[iDomain] = Buffer_Send_nPointPeriodic;
                nElemTotal_s[iDomain] = Buffer_Send_nElemTotal;
                nElemTriangle_s[iDomain] = Buffer_Send_nElemTriangle;
                nElemRectangle_s[iDomain] = Buffer_Send_nElemRectangle;
                nElemTetrahedron_s[iDomain] = Buffer_Send_nElemTetrahedron;
                nElemHexahedron_s[iDomain] = Buffer_Send_nElemHexahedron;
                nElemPrism_s[iDomain] = Buffer_Send_nElemPrism;
                nElemPyramid_s[iDomain] = Buffer_Send_nElemPyramid;

                /*--- Total counts for allocating send buffers below ---*/
                Buffer_Size_Coord += nPointTotal_s[iDomain] * nDim_s[iDomain];
                Buffer_Size_Color += nPointTotal_s[iDomain];
                Buffer_Size_GlobalPointIndex += nPointTotal_s[iDomain];
                Buffer_Size_Triangle += nElemTriangle_s[iDomain];
                Buffer_Size_Rectangle += nElemRectangle_s[iDomain];
                Buffer_Size_Tetrahedron += nElemTetrahedron_s[iDomain];
                Buffer_Size_Hexahedron += nElemHexahedron_s[iDomain];
                Buffer_Size_Prism += nElemPrism_s[iDomain];
                Buffer_Size_Pyramid += nElemPyramid_s[iDomain];
                Buffer_Size_GlobElem += nElemTotal_s[iDomain];
            }

            /*--- Allocate the buffer std::vectors in the appropiate domain (master, iDomain) ---*/
            Buffer_Send_Coord = new double[Buffer_Size_Coord];
            Buffer_Send_Color = new unsigned long[Buffer_Size_Color];
            Buffer_Send_GlobalPointIndex = new unsigned long[Buffer_Size_GlobalPointIndex];
            Buffer_Send_Triangle = new unsigned long[Buffer_Size_Triangle*TBOX::N_POINTS_TRIANGLE];
            Buffer_Send_Rectangle = new unsigned long[Buffer_Size_Rectangle*TBOX::N_POINTS_QUADRILATERAL];
            Buffer_Send_Tetrahedron = new unsigned long[Buffer_Size_Tetrahedron*TBOX::N_POINTS_TETRAHEDRON];
            Buffer_Send_Hexahedron = new unsigned long[Buffer_Size_Hexahedron*TBOX::N_POINTS_HEXAHEDRON];
            Buffer_Send_Prism = new unsigned long[Buffer_Size_Prism*TBOX::N_POINTS_PRISM];
            Buffer_Send_Pyramid = new unsigned long[Buffer_Size_Pyramid*TBOX::N_POINTS_PYRAMID];
            Buffer_Send_GlobElem = new unsigned long[Buffer_Size_GlobElem];

            Local_to_global_Triangle = new unsigned long[Buffer_Size_Triangle];
            Local_to_global_Rectangle = new unsigned long[Buffer_Size_Rectangle];
            Local_to_global_Tetrahedron = new unsigned long[Buffer_Size_Tetrahedron];
            Local_to_global_Hexahedron = new unsigned long[Buffer_Size_Hexahedron];
            Local_to_global_Prism = new unsigned long[Buffer_Size_Prism];
            Local_to_global_Pyramid = new unsigned long[Buffer_Size_Pyramid];

            /*--- Initialize the counters for the larger send buffers (by domain) ---*/
            ElemTotal_Counter = 0;
            PointTotal_Counter = 0;
            PointDomain_Counter = 0;
            /*--- WARNING: check the next two counters ---*/
            PointPeriodic_Counter = 0;
            PointGhost_Counter = 0;
            ElemTriangle_Counter = 0;
            ElemRectangle_Counter = 0;
            ElemTetrahedron_Counter = 0;
            ElemHexahedron_Counter = 0;
            ElemPrism_Counter = 0;
            ElemPyramid_Counter = 0;

            /*--- Now that we know the sizes of the point, elem, etc. arrays, we can
            allocate and send the information in large chunks to all processors. ---*/

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {

                /*--- A rank does not communicate with itself through MPI ---*/
                if (rank != iDomain)
                {
#ifdef HAVE_MPI

                    /*--- Communicate the counts to iDomain with non-blocking sends ---*/

                    MPI_Isend(&nDim_s[iDomain], 1, MPI_UNSIGNED_SHORT, iDomain,
                        iDomain * 13 + 0, MPI_COMM_WORLD, &send_req[0]);

                    MPI_Isend(&nZone_s[iDomain], 1, MPI_UNSIGNED_SHORT, iDomain,
                        iDomain * 13 + 1, MPI_COMM_WORLD, &send_req[1]);

                    MPI_Isend(&nPointTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 2, MPI_COMM_WORLD, &send_req[2]);

                    MPI_Isend(&nPointDomainTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 3, MPI_COMM_WORLD, &send_req[3]);

                    MPI_Isend(&nPointGhost_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 4, MPI_COMM_WORLD, &send_req[4]);

                    MPI_Isend(&nPointPeriodic_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 5, MPI_COMM_WORLD, &send_req[5]);

                    MPI_Isend(&nElemTotal_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 6, MPI_COMM_WORLD, &send_req[6]);

                    MPI_Isend(&nElemTriangle_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 7, MPI_COMM_WORLD, &send_req[7]);

                    MPI_Isend(&nElemRectangle_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 8, MPI_COMM_WORLD, &send_req[8]);

                    MPI_Isend(&nElemTetrahedron_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 9, MPI_COMM_WORLD, &send_req[9]);

                    MPI_Isend(&nElemHexahedron_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 10, MPI_COMM_WORLD, &send_req[10]);

                    MPI_Isend(&nElemPrism_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 11, MPI_COMM_WORLD, &send_req[11]);

                    MPI_Isend(&nElemPyramid_s[iDomain], 1, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 13 + 12, MPI_COMM_WORLD, &send_req[12]);

#endif

                }
                else
                {
                    /*--- If iDomain = rank, we simply copy values into place in memory ---*/
                    nDim = nDim_s[iDomain];
                    nZone = nZone_s[iDomain];
                    nPointTotal = nPointTotal_s[iDomain];
                    nPointDomainTotal = nPointDomainTotal_s[iDomain];
                    nPointGhost = nPointGhost_s[iDomain];
                    nPointPeriodic = nPointPeriodic_s[iDomain];
                    nElemTotal = nElemTotal_s[iDomain];
                    nElemTriangle = nElemTriangle_s[iDomain];
                    nElemRectangle = nElemRectangle_s[iDomain];
                    nElemTetrahedron = nElemTetrahedron_s[iDomain];
                    nElemHexahedron = nElemHexahedron_s[iDomain];
                    nElemPrism = nElemPrism_s[iDomain];
                    nElemPyramid = nElemPyramid_s[iDomain];

                    nDim_r[iDomain] = nDim_s[iDomain];
                    nZone_r[iDomain] = nZone_s[iDomain];
                    nPointTotal_r[iDomain] = nPointTotal_s[iDomain];
                    nPointDomainTotal_r[iDomain] = nPointDomainTotal_s[iDomain];
                    nPointPeriodic_r[iDomain] = nPointPeriodic_s[iDomain];
                    nElemTotal_r[iDomain] = nElemTotal_s[iDomain];
                    nElemTriangle_r[iDomain] = nElemTriangle_s[iDomain];
                    nElemRectangle_r[iDomain] = nElemRectangle_s[iDomain];
                    nElemTetrahedron_r[iDomain] = nElemTetrahedron_s[iDomain];
                    nElemHexahedron_r[iDomain] = nElemHexahedron_s[iDomain];
                    nElemPrism_r[iDomain] = nElemPrism_s[iDomain];
                    nElemPyramid_r[iDomain] = nElemPyramid_s[iDomain];

                    nPointTotal_r_tot += nPointTotal_r[iDomain];
                    nPointDomainTotal_r_tot += nPointDomainTotal_r[iDomain];
                    nPointGhost_r_tot += nPointGhost_r[iDomain];
                    nPointPeriodic_r_tot += nPointPeriodic_r[iDomain];
                    nElemTotal_r_tot += nElemTotal_r[iDomain];
                    nElemTriangle_r_tot += nElemTriangle_r[iDomain];
                    nElemRectangle_r_tot += nElemRectangle_r[iDomain];
                    nElemTetrahedron_r_tot += nElemTetrahedron_r[iDomain];
                    nElemHexahedron_r_tot += nElemHexahedron_r[iDomain];
                    nElemPrism_r_tot += nElemPrism_r[iDomain];
                    nElemPyramid_r_tot += nElemPyramid_r[iDomain];

                }

                /*--- Receive the counts. All processors are sending their counters to
                iDomain up above, so only iDomain needs to perform the recv here from
                all other ranks. ---*/
                if (rank == iDomain)
                {
                    for (jDomain = 0; jDomain < size; jDomain++)
                    {
                        /*--- A rank does not communicate with itself through MPI ---*/
                        if (rank != jDomain)
                        {

#ifdef HAVE_MPI

                            /*--- Recv the data by probing for the current sender, jDomain,
                            first and then receiving the values from it. ---*/

                            MPI_Probe(jDomain, 13 * rank + 0, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nDim_r[jDomain], 1, MPI_UNSIGNED_SHORT, jDomain,
                                rank * 13 + 0, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 1, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nZone_r[jDomain], 1, MPI_UNSIGNED_SHORT, jDomain,
                                rank * 13 + 1, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 2, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nPointTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 2, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 3, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nPointDomainTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 3, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 4, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nPointGhost_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 4, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 5, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nPointPeriodic_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 5, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 6, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemTotal_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 6, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 7, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemTriangle_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 7, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 8, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemRectangle_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 8, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 9, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemTetrahedron_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 9, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 10, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemHexahedron_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 10, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 11, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemPrism_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 11, MPI_COMM_WORLD, &status2);

                            MPI_Probe(jDomain, 13 * rank + 12, MPI_COMM_WORLD, &status2);
                            MPI_Recv(&nElemPyramid_r[jDomain], 1, MPI_UNSIGNED_LONG, jDomain,
                                rank * 13 + 12, MPI_COMM_WORLD, &status2);

#endif

                            /*--- These are the cumulative totals that we will recv below. ----*/

                            nPointTotal_r_tot += nPointTotal_r[jDomain];
                            nPointDomainTotal_r_tot += nPointDomainTotal_r[jDomain];
                            nPointGhost_r_tot += nPointGhost_r[jDomain];
                            nPointPeriodic_r_tot += nPointPeriodic_r[jDomain];
                            nElemTotal_r_tot += nElemTotal_r[jDomain];
                            nElemTriangle_r_tot += nElemTriangle_r[jDomain];
                            nElemRectangle_r_tot += nElemRectangle_r[jDomain];
                            nElemTetrahedron_r_tot += nElemTetrahedron_r[jDomain];
                            nElemHexahedron_r_tot += nElemHexahedron_r[jDomain];
                            nElemPrism_r_tot += nElemPrism_r[jDomain];
                            nElemPyramid_r_tot += nElemPyramid_r[jDomain];
                        }
                    }
                }
            }

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                /*--- Wait for the non-blocking sends to complete. ---*/

#ifdef HAVE_MPI
                if (rank != iDomain) MPI_Waitall(13, send_req, send_stat);
                MPI_Barrier(MPI_COMM_WORLD);
#endif
                //cout << " ==== Rank " << rank << " finished sending counts " << endl;
            }

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                /*--- Above was number of elements to send and receive, and here is where
                we send/recv the actual elements. Here you're sending global index values,
                which are later changed to local. ---*/

                /*--- Set the value of the interior geometry. Initialize counters. ---*/

                iElemTotal = 0;
                iPointTotal = 0;
                iPointDomain = 0;
                iPointPeriodic = nPointDomainTotal_s[iDomain];
                iPointGhost = nPointDomainTotal_s[iDomain] + nPointPeriodic_s[iDomain];
                iElemTriangle = 0;
                iElemRectangle = 0;
                iElemTetrahedron = 0;
                iElemHexahedron = 0;
                iElemPrism = 0;
                iElemPyramid = 0;

                /*--- Initialize the global to local mapping ---*/
                for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
                    PointIn[iPoint] = false;

                /*--- Load up the actual elements into the buffers for sending. ---*/
                for (iElem = 0; iElem < geometry->no_of_local_elements; iElem++)
                {

                    /*--- Check if the element belongs to the domain ---*/
                    ElemIn[iElem] = false;
                    for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = geometry->elem[iElem]->GetNode(iNode);
                        if (local_colour_values[iPoint] == iDomain)
                        {
                            ElemIn[iElem] = true; break;
                        }
                    }

                    /*--- If this element should be sent ---*/

                    if (ElemIn[iElem])
                    {

                        /*--- We need to send this element, so add it to the send buffer. The
                        local to global mapping has already been done as a class data member. ---*/

                        Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal] = Local_to_global_elem[iElem];

                        /*--- Loop through the nodes of the current element ---*/

                        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                        {

                            /*--- Get the global index for this node in the element ---*/
                            iPoint = geometry->elem[iElem]->GetNode(iNode);

                            /*--- Store the connectivity for this element for each node ---*/
                            vnodes_local[iNode] = iPoint;

                            /*--- Check if this point has been found previously ---*/

                            if (PointIn[iPoint] == false)
                            {

                                /*--- Check if this node lives on the current rank based on the
                                initial linear partitioning. We are only ever sending nodes that
                                we own in the linear partitioning (no duplicate nodes are sent) ---*/

                                if ((iPoint >= geometry->starting_node[rank]) &&
                                    (iPoint < geometry->ending_node[rank]))
                                {

                                    /*--- Decide whether this is an interior, periodic, or ghost node ---*/

                                    if (local_colour_values[iPoint] == iDomain)
                                    {

                                        /*--- If iDomain owns the point, it must be either an interior
                                        node (iPoint < nPointDomain) or a periodic node. ---*/

                                        if (iPoint > geometry->GetnPointDomain() - 1)
                                            iPointCurrent = iPointPeriodic;
                                        else
                                            iPointCurrent = iPointDomain;

                                    }
                                    else
                                    {

                                        /*--- Otherwise, it must be a ghost point for iDomain ---*/
                                        iPointCurrent = iPointGhost;

                                    }

                                    /*--- Setting global to local, the color, and index. ---*/

                                    PointIn[iPoint] = true;

                                    Buffer_Send_Color[PointTotal_Counter + iPointCurrent] = local_colour_values[iPoint];
                                    Buffer_Send_GlobalPointIndex[PointTotal_Counter + iPointCurrent] = iPoint;

                                    /*--- Get the coordinates for this point ---*/

                                    for (iDim = 0; iDim < nDim_s[iDomain]; iDim++)
                                    {

                                        /*--- iPoint is the global index, but we store everything local
                                        to this rank. So we need to subtract the starting index. All
                                        ranks re-index their points from zero. ---*/
                                        Buffer_Send_Coord[nDim_s[iDomain] * (PointTotal_Counter + iPointCurrent) + iDim] = geometry->node[iPoint - geometry->starting_node[rank]]->GetCoord(iDim);
                                    }

                                    /*--- Increment our counters ---*/
                                    if (local_colour_values[iPoint] == iDomain)
                                    {
                                        if (iPoint > geometry->GetnPointDomain() - 1)
                                            iPointPeriodic++;
                                        else
                                            iPointDomain++;
                                    }
                                    else
                                        iPointGhost++;

                                    /*--- Increment the total number of points we're sending ---*/
                                    iPointTotal++;
                                }
                            }
                        }

                        /*--- Load the connectivity for the current element into the send buffer.
                        Also store the local to global mapping for the elements.
                        Note that we are using the vnode_local array we filled above to store
                        the connectivity. Loop through each element type. ---*/

                        switch (geometry->elem[iElem]->GetVTK_Type())
                        {
                        case TBOX::TRIANGLE:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Triangle[3 * (ElemTriangle_Counter + iElemTriangle) + iNode] = vnodes_local[iNode];
                            Local_to_global_Triangle[ElemTriangle_Counter + iElemTriangle] = Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal];
                            iElemTriangle++; break;
                        case TBOX::RECTANGLE:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Rectangle[4 * (ElemRectangle_Counter + iElemRectangle) + iNode] = vnodes_local[iNode];
                            Local_to_global_Rectangle[ElemRectangle_Counter + iElemRectangle] = Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal];
                            iElemRectangle++; break;
                        case TBOX::TETRAHEDRON:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Tetrahedron[4 * (ElemTetrahedron_Counter + iElemTetrahedron) + iNode] = vnodes_local[iNode];
                            Local_to_global_Tetrahedron[ElemTetrahedron_Counter + iElemTetrahedron] = Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal];
                            iElemTetrahedron++; break;
                        case TBOX::HEXAHEDRON:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Hexahedron[8 * (ElemHexahedron_Counter + iElemHexahedron) + iNode] = vnodes_local[iNode];
                            Local_to_global_Hexahedron[ElemHexahedron_Counter + iElemHexahedron] = Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal];
                            iElemHexahedron++; break;
                        case TBOX::PRISM:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Prism[6 * (ElemPrism_Counter + iElemPrism) + iNode] = vnodes_local[iNode];
                            Local_to_global_Prism[ElemPrism_Counter + iElemPrism] = Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal];
                            iElemPrism++; break;
                        case TBOX::PYRAMID:
                            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
                                Buffer_Send_Pyramid[5 * (ElemPyramid_Counter + iElemPyramid) + iNode] = vnodes_local[iNode];
                            Local_to_global_Pyramid[ElemPyramid_Counter + iElemPyramid] = Buffer_Send_GlobElem[ElemTotal_Counter + iElemTotal];
                            iElemPyramid++; break;
                        }

                        /*--- Regardless of the type, increment the total count ---*/
                        iElemTotal++;

                    }
                }

                /*--- Send the buffers with the geometrical information ---*/
                if (iDomain != rank)
                {

#ifdef HAVE_MPI

                    /*--- Communicate the coordinates, global index, colors, and element
                    date to iDomain with non-blocking sends. ---*/

                    MPI_Isend(&Buffer_Send_Coord[PointTotal_Counter*nDim_s[iDomain]],
                        nPointTotal_s[iDomain] * nDim_s[iDomain], MPI_DOUBLE, iDomain,
                        iDomain * 16 + 0, MPI_COMM_WORLD, &send_req[0]);

                    MPI_Isend(&Buffer_Send_GlobalPointIndex[PointTotal_Counter],
                        nPointTotal_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 1, MPI_COMM_WORLD, &send_req[1]);

                    MPI_Isend(&Buffer_Send_Color[PointTotal_Counter],
                        nPointTotal_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 2, MPI_COMM_WORLD, &send_req[2]);

                    MPI_Isend(&Buffer_Send_Triangle[ElemTriangle_Counter * 3],
                        nElemTriangle_s[iDomain] * 3, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 3, MPI_COMM_WORLD, &send_req[3]);

                    MPI_Isend(&Buffer_Send_Rectangle[ElemRectangle_Counter * 4],
                        nElemRectangle_s[iDomain] * 4, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 4, MPI_COMM_WORLD, &send_req[4]);

                    MPI_Isend(&Buffer_Send_Tetrahedron[ElemTetrahedron_Counter * 4],
                        nElemTetrahedron_s[iDomain] * 4, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 5, MPI_COMM_WORLD, &send_req[5]);

                    MPI_Isend(&Buffer_Send_Hexahedron[ElemHexahedron_Counter * 8],
                        nElemHexahedron_s[iDomain] * 8, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 6, MPI_COMM_WORLD, &send_req[6]);

                    MPI_Isend(&Buffer_Send_Prism[ElemPrism_Counter * 6],
                        nElemPrism_s[iDomain] * 6, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 7, MPI_COMM_WORLD, &send_req[7]);

                    MPI_Isend(&Buffer_Send_Pyramid[ElemPyramid_Counter * 5],
                        nElemPyramid_s[iDomain] * 5, MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 8, MPI_COMM_WORLD, &send_req[8]);

                    MPI_Isend(&Buffer_Send_GlobElem[ElemTotal_Counter],
                        nElemTotal_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 9, MPI_COMM_WORLD, &send_req[9]);

                    MPI_Isend(&Local_to_global_Triangle[ElemTriangle_Counter],
                        nElemTriangle_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 10, MPI_COMM_WORLD, &send_req[10]);

                    MPI_Isend(&Local_to_global_Rectangle[ElemRectangle_Counter],
                        nElemRectangle_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 11, MPI_COMM_WORLD, &send_req[11]);

                    MPI_Isend(&Local_to_global_Tetrahedron[ElemTetrahedron_Counter],
                        nElemTetrahedron_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 12, MPI_COMM_WORLD, &send_req[12]);

                    MPI_Isend(&Local_to_global_Hexahedron[ElemHexahedron_Counter],
                        nElemHexahedron_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 13, MPI_COMM_WORLD, &send_req[13]);

                    MPI_Isend(&Local_to_global_Prism[ElemPrism_Counter],
                        nElemPrism_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 14, MPI_COMM_WORLD, &send_req[14]);

                    MPI_Isend(&Local_to_global_Pyramid[ElemPyramid_Counter],
                        nElemPyramid_s[iDomain], MPI_UNSIGNED_LONG, iDomain,
                        iDomain * 16 + 15, MPI_COMM_WORLD, &send_req[15]);

#endif

                }
                else
                {

                    /*--- Allocate local memory for the local recv of the elements ---*/

                    Buffer_Receive_Coord_loc = new double[nPointTotal_s[iDomain] * nDim_s[iDomain]];

                    Buffer_Receive_GlobalPointIndex_loc = new unsigned long[nPointTotal_s[iDomain]];
                    Buffer_Receive_Color_loc = new unsigned long[nPointTotal_s[iDomain]];
                    Buffer_Receive_Triangle_loc = new unsigned long[nElemTriangle_s[iDomain] * TBOX::N_POINTS_TRIANGLE];
                    Buffer_Receive_Rectangle_loc = new unsigned long[nElemRectangle_s[iDomain] * TBOX::N_POINTS_QUADRILATERAL];
                    Buffer_Receive_Tetrahedron_loc = new unsigned long[nElemTetrahedron_s[iDomain] * TBOX::N_POINTS_TETRAHEDRON];
                    Buffer_Receive_Hexahedron_loc = new unsigned long[nElemHexahedron_s[iDomain] * TBOX::N_POINTS_HEXAHEDRON];
                    Buffer_Receive_Prism_loc = new unsigned long[nElemPrism_s[iDomain] * TBOX::N_POINTS_PRISM];
                    Buffer_Receive_Pyramid_loc = new unsigned long[nElemPyramid_s[iDomain] * TBOX::N_POINTS_PYRAMID];
                    Buffer_Receive_GlobElem_loc = new unsigned long[nElemTotal_s[iDomain]];

                    Buffer_Receive_Triangle_presence_loc = new unsigned long[nElemTriangle_s[iDomain]];
                    Buffer_Receive_Rectangle_presence_loc = new unsigned long[nElemRectangle_s[iDomain]];
                    Buffer_Receive_Tetrahedron_presence_loc = new unsigned long[nElemTetrahedron_s[iDomain]];
                    Buffer_Receive_Hexahedron_presence_loc = new unsigned long[nElemHexahedron_s[iDomain]];
                    Buffer_Receive_Prism_presence_loc = new unsigned long[nElemPrism_s[iDomain]];
                    Buffer_Receive_Pyramid_presence_loc = new unsigned long[nElemPyramid_s[iDomain]];

                    for (iter = 0; iter < nPointTotal_s[iDomain] * nDim_s[iDomain]; iter++)
                        Buffer_Receive_Coord_loc[iter] = Buffer_Send_Coord[PointTotal_Counter*nDim_s[iDomain] + iter];

                    for (iter = 0; iter < nPointTotal_s[iDomain]; iter++)
                    {
                        Buffer_Receive_GlobalPointIndex_loc[iter] = Buffer_Send_GlobalPointIndex[PointTotal_Counter + iter];
                        Buffer_Receive_Color_loc[iter] = Buffer_Send_Color[PointTotal_Counter + iter];
                    }

                    for (iter = 0; iter < nElemTriangle_s[iDomain] * TBOX::N_POINTS_TRIANGLE; iter++)
                        Buffer_Receive_Triangle_loc[iter] = Buffer_Send_Triangle[ElemTriangle_Counter*TBOX::N_POINTS_TRIANGLE + iter];

                    for (iter = 0; iter < nElemRectangle_s[iDomain] * TBOX::N_POINTS_QUADRILATERAL; iter++)
                        Buffer_Receive_Rectangle_loc[iter] = Buffer_Send_Rectangle[ElemRectangle_Counter*TBOX::N_POINTS_QUADRILATERAL + iter];

                    for (iter = 0; iter < nElemTetrahedron_s[iDomain] * TBOX::N_POINTS_TETRAHEDRON; iter++)
                        Buffer_Receive_Tetrahedron_loc[iter] = Buffer_Send_Tetrahedron[ElemTetrahedron_Counter*TBOX::N_POINTS_TETRAHEDRON + iter];

                    for (iter = 0; iter < nElemHexahedron_s[iDomain] * TBOX::N_POINTS_HEXAHEDRON; iter++)
                        Buffer_Receive_Hexahedron_loc[iter] = Buffer_Send_Hexahedron[ElemHexahedron_Counter*TBOX::N_POINTS_HEXAHEDRON + iter];

                    for (iter = 0; iter < nElemPrism_s[iDomain] * TBOX::N_POINTS_PRISM; iter++)
                        Buffer_Receive_Prism_loc[iter] = Buffer_Send_Prism[ElemPrism_Counter*TBOX::N_POINTS_PRISM + iter];

                    for (iter = 0; iter < nElemPyramid_s[iDomain] * TBOX::N_POINTS_PYRAMID; iter++)
                        Buffer_Receive_Pyramid_loc[iter] = Buffer_Send_Pyramid[ElemPyramid_Counter*TBOX::N_POINTS_PYRAMID + iter];

                    for (unsigned long i = 0; i < nElemTotal_s[iDomain]; i++)
                    {
                        Buffer_Receive_GlobElem_loc[i] = Buffer_Send_GlobElem[ElemTotal_Counter + i];
                    }

                    for (unsigned long i = 0; i < nElemTriangle_s[iDomain]; i++)
                    {
                        Buffer_Receive_Triangle_presence_loc[i] = Local_to_global_Triangle[ElemTriangle_Counter + i];
                    }

                    for (unsigned long i = 0; i < nElemRectangle_s[iDomain]; i++)
                    {
                        Buffer_Receive_Rectangle_presence_loc[i] = Local_to_global_Rectangle[ElemRectangle_Counter + i];
                    }

                    for (unsigned long i = 0; i < nElemTetrahedron_s[iDomain]; i++)
                    {
                        Buffer_Receive_Tetrahedron_presence_loc[i] = Local_to_global_Tetrahedron[ElemTetrahedron_Counter + i];
                    }

                    for (unsigned long i = 0; i < nElemHexahedron_s[iDomain]; i++)
                    {
                        Buffer_Receive_Hexahedron_presence_loc[i] = Local_to_global_Hexahedron[ElemHexahedron_Counter + i];
                    }

                    for (unsigned long i = 0; i < nElemPrism_s[iDomain]; i++)
                    {
                        Buffer_Receive_Prism_presence_loc[i] = Local_to_global_Prism[ElemPrism_Counter + i];
                    }

                    for (unsigned long i = 0; i < nElemPyramid_s[iDomain]; i++)
                    {
                        Buffer_Receive_Pyramid_presence_loc[i] = Local_to_global_Pyramid[ElemPyramid_Counter + i];
                    }
                }

                /*--- Increment the counters for the send buffers (iDomain loop) ---*/

                ElemTotal_Counter += iElemTotal;
                PointTotal_Counter += iPointTotal;
                PointDomain_Counter += iPointDomain;
                /*--- WARNING: check the next two counters ---*/
                PointPeriodic_Counter += iPointPeriodic;
                PointGhost_Counter += iPointGhost;
                ElemTriangle_Counter += iElemTriangle;
                ElemRectangle_Counter += iElemRectangle;
                ElemTetrahedron_Counter += iElemTetrahedron;
                ElemHexahedron_Counter += iElemHexahedron;
                ElemPrism_Counter += iElemPrism;
                ElemPyramid_Counter += iElemPyramid;
            }

#ifdef HAVE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            //cout << " ==== Rank " << rank << " sent all point elem data " << endl;

            /*--- The next section begins the recv of all data for the interior
            points/elements in the mesh. First, create the domain structures for
            the points on this rank ---*/
            nPoint = nPointTotal_r_tot; iPoint = 0;
            nPointDomain = nPointDomainTotal_r_tot;
            node = new GRID::GRID_DGPoint*[nPoint];
            Local_to_Global_Point = new long[nPoint];

            /*--- Array initialization ---*/
            for (iPoint = 0; iPoint < nPointTotal_r_tot; iPoint++)
            {
                Local_to_Global_Point[iPoint] = -1;
            }

            /*--- Initialize some counters ---*/

            unsigned long temp_node_count = 0;
            unsigned long temp_node_count_periodic = nPointDomainTotal_r_tot;
            unsigned long temp_node_count_ghost = nPointDomainTotal_r_tot + nPointPeriodic_r_tot;

            /*--- First, we recv all of the point data ---*/
            for (iDomain = 0; iDomain < size; iDomain++)
            {

                if (rank != iDomain)
                {
#ifdef HAVE_MPI

                    /*--- Allocate the receive buffer std::vector. Send the colors so that we
                    know whether what we recv is an owned or halo node. ---*/

                    Buffer_Receive_Coord = new double[nPointTotal_r[iDomain] * nDim_r[iDomain]];
                    Buffer_Receive_Color = new unsigned long[nPointTotal_r[iDomain]];
                    Buffer_Receive_GlobalPointIndex = new unsigned long[nPointTotal_r[iDomain]];

                    /*--- Receive the buffers with the coords, global index, and colors ---*/

                    MPI_Probe(iDomain, rank * 16 + 0, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_DOUBLE, &recv_count);
                    MPI_Recv(Buffer_Receive_Coord, recv_count, MPI_DOUBLE,
                        source, rank * 16 + 0, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 1, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_GlobalPointIndex, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 1, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 2, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Color, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 2, MPI_COMM_WORLD, &status2);

                    /*--- Wait for the three recv above to complete ---*/

                    //if (rank != iDomain)  MPI_Waitall(3, send_req, send_stat);

                    /*--- Loop over all of the points that we have recv'd and store the
                    coords, global index, and colors ---*/

                    unsigned long index = 0;
                    for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++) {

                        /*--- If this rank owns the current point ---*/

                        if (Buffer_Receive_Color[iPoint] == rank) {

                            /*--- If iDomain owns the point, it must be either an interior
                            node (iPoint < nPointDomain) or a periodic node. ---*/

                            if (Buffer_Receive_GlobalPointIndex[iPoint] > geometry->GetnPointDomain() - 1) {

                                /*--- Set the starting point for the local index of the recv points.
                                The temp_node_count increments for the interior nodes, between 0 up
                                to nPointDomain-1. ---*/
                                index = temp_node_count_periodic;

                                /*--- Get the global index ---*/
                                Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex[iPoint];

                                /*--- Allocating the Point object ---*/
                                if (nDim == 2) node[index] = new GRID_DGPoint::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0],
                                    Buffer_Receive_Coord[iPoint*nDim + 1],
                                    Local_to_Global_Point[index], config);
                                if (nDim == 3) node[index] = new GRID_DGPoint::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0],
                                    Buffer_Receive_Coord[iPoint*nDim + 1],
                                    Buffer_Receive_Coord[iPoint*nDim + 2],
                                    Local_to_Global_Point[index], config);

                                /*--- Set the color ---*/
                                node[index]->SetColor(Buffer_Receive_Color[iPoint]);

                                /*--- Increment the interior node counter ---*/
                                temp_node_count_periodic++;


                            }

                            else {


                                /*--- Set the starting point for the local index of the recv points.
                                The temp_node_count increments for the interior nodes, between 0 up
                                to nPointDomain-1. ---*/
                                index = temp_node_count;

                                /*--- Get the global index ---*/
                                Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex[iPoint];

                                /*--- Allocating the Point object ---*/
                                if (nDim == 2) node[index] = new GRID_DGPoint::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0],
                                    Buffer_Receive_Coord[iPoint*nDim + 1],
                                    Local_to_Global_Point[index], config);
                                if (nDim == 3) node[index] = new GRID_DGPoint::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0],
                                    Buffer_Receive_Coord[iPoint*nDim + 1],
                                    Buffer_Receive_Coord[iPoint*nDim + 2],
                                    Local_to_Global_Point[index], config);

                                /*--- Set the color ---*/
                                node[index]->SetColor(Buffer_Receive_Color[iPoint]);

                                /*--- Increment the interior node counter ---*/
                                temp_node_count++;
                            }
                        }
                        else {

                            /*--- Set the starting point for the local index of the recv points.
                            The temp_node_count_domain increments for the ghost nodes, between
                            nPointDomain up to nPoint. ---*/

                            index = temp_node_count_ghost;

                            /*--- Get the global index ---*/
                            Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex[iPoint];

                            /*--- Allocating the Point object ---*/
                            if (nDim == 2) node[index] = new GRID_DGPoint::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0],
                                Buffer_Receive_Coord[iPoint*nDim + 1],
                                Local_to_Global_Point[index], config);
                            if (nDim == 3) node[index] = new GRID_DGPoint::GRID_DGPoint(Buffer_Receive_Coord[iPoint*nDim + 0],
                                Buffer_Receive_Coord[iPoint*nDim + 1],
                                Buffer_Receive_Coord[iPoint*nDim + 2],
                                Local_to_Global_Point[index], config);

                            /*--- Set the color ---*/
                            node[index]->SetColor(Buffer_Receive_Color[iPoint]);

                            /*--- Increment the ghost node counter ---*/
                            temp_node_count_ghost++;

                        }
                    }

                    /*--- Delete memory for recv the point stuff ---*/
                    delete[] Buffer_Receive_Coord;
                    delete[] Buffer_Receive_Color;
                    delete[] Buffer_Receive_GlobalPointIndex;

#endif

                }
                else
                {

                    /*--- Recv the point data from ourselves (same procedure as above) ---*/

                    unsigned long index = 0;
                    for (iPoint = 0; iPoint < nPointTotal_r[iDomain]; iPoint++)
                    {

                        if (Buffer_Receive_Color_loc[iPoint] == rank)
                        {

                            /*--- If iDomain owns the point, it must be either an interior
                            node (iPoint < nPointDomain) or a periodic node. ---*/

                            if (Buffer_Receive_GlobalPointIndex_loc[iPoint] > geometry->GetnPointDomain() - 1)
                            {

                                index = temp_node_count_periodic;

                                Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex_loc[iPoint];
                                if (nDim == 2) node[index] = new GRID::GRID_DGPoint(Buffer_Receive_Coord_loc[iPoint*nDim + 0],
                                    Buffer_Receive_Coord_loc[iPoint*nDim + 1],
                                    Local_to_Global_Point[index], config);
                                if (nDim == 3) node[index] = new GRID::GRID_DGPoint(Buffer_Receive_Coord_loc[iPoint*nDim + 0],
                                    Buffer_Receive_Coord_loc[iPoint*nDim + 1],
                                    Buffer_Receive_Coord_loc[iPoint*nDim + 2],
                                    Local_to_Global_Point[index], config);
                                node[index]->SetColor(Buffer_Receive_Color_loc[iPoint]);
                                temp_node_count_periodic++;
                            }
                            else
                            {
                                index = temp_node_count;
                                Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex_loc[iPoint];
                                if (nDim == 2) node[index] = new GRID::GRID_DGPoint(Buffer_Receive_Coord_loc[iPoint*nDim + 0],
                                    Buffer_Receive_Coord_loc[iPoint*nDim + 1],
                                    Local_to_Global_Point[index], config);
                                if (nDim == 3) node[index] = new GRID::GRID_DGPoint(Buffer_Receive_Coord_loc[iPoint*nDim + 0],
                                    Buffer_Receive_Coord_loc[iPoint*nDim + 1],
                                    Buffer_Receive_Coord_loc[iPoint*nDim + 2],
                                    Local_to_Global_Point[index], config);
                                node[index]->SetColor(Buffer_Receive_Color_loc[iPoint]);
                                temp_node_count++;
                            }
                        }
                        else
                        {

                            index = temp_node_count_ghost;
                            Local_to_Global_Point[index] = Buffer_Receive_GlobalPointIndex_loc[iPoint];
                            if (nDim == 2) node[index] = new GRID::GRID_DGPoint(Buffer_Receive_Coord_loc[iPoint*nDim + 0],
                                Buffer_Receive_Coord_loc[iPoint*nDim + 1],
                                Local_to_Global_Point[index], config);
                            if (nDim == 3) node[index] = new GRID::GRID_DGPoint(Buffer_Receive_Coord_loc[iPoint*nDim + 0],
                                Buffer_Receive_Coord_loc[iPoint*nDim + 1],
                                Buffer_Receive_Coord_loc[iPoint*nDim + 2],
                                Local_to_Global_Point[index], config);
                            node[index]->SetColor(Buffer_Receive_Color_loc[iPoint]);
                            temp_node_count_ghost++;
                        }
                    }

                    delete[] Buffer_Receive_Coord_loc;
                    delete[] Buffer_Receive_Color_loc;
                    delete[] Buffer_Receive_GlobalPointIndex_loc;

                }
            }

            /*--- Get the global to local mapping --- */
            for (iPoint = 0; iPoint < nPointTotal_r_tot; iPoint++)
            {
                Global_to_local_Point_recv[Local_to_Global_Point[iPoint]] = iPoint;
            }


            //cout << " ==== Rank " << rank << " recv of point data finished" << endl;
#ifdef HAVE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            /*--- Recv all of the element data. First decide which elements we need to own on each proc ---*/
            iElem = 0;
            for (iDomain = 0; iDomain < size; iDomain++)
            {
                if (rank != iDomain)
                {
#ifdef HAVE_MPI

                    /*--- Allocate memory for the element recv ---*/

                    Buffer_Receive_Triangle_presence[iDomain] = new unsigned long[nElemTriangle_r[iDomain]];
                    Buffer_Receive_Rectangle_presence[iDomain] = new unsigned long[nElemRectangle_r[iDomain]];
                    Buffer_Receive_Tetrahedron_presence[iDomain] = new unsigned long[nElemTetrahedron_r[iDomain]];
                    Buffer_Receive_Hexahedron_presence[iDomain] = new unsigned long[nElemHexahedron_r[iDomain]];
                    Buffer_Receive_Prism_presence[iDomain] = new unsigned long[nElemPrism_r[iDomain]];
                    Buffer_Receive_Pyramid_presence[iDomain] = new unsigned long[nElemPyramid_r[iDomain]];

                    /*--- Recv the element data ---*/

                    MPI_Probe(iDomain, rank * 16 + 10, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(&Buffer_Receive_Triangle_presence[iDomain][0],
                        recv_count, MPI_UNSIGNED_LONG, source,
                        rank * 16 + 10, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 11, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(&Buffer_Receive_Rectangle_presence[iDomain][0],
                        recv_count, MPI_UNSIGNED_LONG, source,
                        rank * 16 + 11, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 12, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(&Buffer_Receive_Tetrahedron_presence[iDomain][0],
                        recv_count, MPI_UNSIGNED_LONG, source,
                        rank * 16 + 12, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 13, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(&Buffer_Receive_Hexahedron_presence[iDomain][0],
                        recv_count, MPI_UNSIGNED_LONG, source,
                        rank * 16 + 13, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 14, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(&Buffer_Receive_Prism_presence[iDomain][0],
                        recv_count, MPI_UNSIGNED_LONG, source,
                        rank * 16 + 14, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 15, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(&Buffer_Receive_Pyramid_presence[iDomain][0],
                        recv_count, MPI_UNSIGNED_LONG, source,
                        rank * 16 + 15, MPI_COMM_WORLD, &status2);

                    /*--- Wait to complete the above sends ---*/

                    //if (rank!=iDomain)  MPI_Waitall(6, &send_req[10], &send_stat[10]);

                    /*--- Allocating the elements after the recv ---*/

                    for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++) {
                        if (Triangle_presence[Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]] == false) {
                            Triangle_presence[Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]] = true;
                            iElem++;
                        }
                    }

                    for (iElemRectangle = 0; iElemRectangle < nElemRectangle_r[iDomain]; iElemRectangle++) {
                        if (Rectangle_presence[Buffer_Receive_Rectangle_presence[iDomain][iElemRectangle]] == false) {
                            Rectangle_presence[Buffer_Receive_Rectangle_presence[iDomain][iElemRectangle]] = true;
                            iElem++;
                        }
                    }

                    for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++) {
                        if (Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]] == false) {
                            Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]] = true;
                            iElem++;
                        }
                    }

                    for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++) {
                        if (Hexahedron_presence[Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]] == false) {
                            Hexahedron_presence[Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]] = true;
                            iElem++;
                        }
                    }

                    for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++) {
                        if (Prism_presence[Buffer_Receive_Prism_presence[iDomain][iElemPrism]] == false) {
                            Prism_presence[Buffer_Receive_Prism_presence[iDomain][iElemPrism]] = true;
                            iElem++;
                        }
                    }

                    for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++) {
                        if (Pyramid_presence[Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]] == false) {
                            Pyramid_presence[Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]] = true;
                            iElem++;
                        }
                    }

#endif

                }
                else
                {
                    /*--- Store the element data from our own local rank info ---*/
                    for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++)
                    {
                        if (Triangle_presence[Buffer_Receive_Triangle_presence_loc[iElemTriangle]] == false)
                        {
                            Triangle_presence[Buffer_Receive_Triangle_presence_loc[iElemTriangle]] = true;
                            iElem++;
                        }
                    }

                    for (iElemRectangle = 0; iElemRectangle < nElemRectangle_r[iDomain]; iElemRectangle++)
                    {
                        if (Rectangle_presence[Buffer_Receive_Rectangle_presence_loc[iElemRectangle]] == false)
                        {
                            Rectangle_presence[Buffer_Receive_Rectangle_presence_loc[iElemRectangle]] = true;
                            iElem++;
                        }
                    }

                    for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++)
                    {
                        if (Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]] == false)
                        {
                            Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]] = true;
                            iElem++;
                        }
                    }

                    for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++)
                    {
                        if (Hexahedron_presence[Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]] == false)
                        {
                            Hexahedron_presence[Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]] = true;
                            iElem++;
                        }
                    }

                    for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++)
                    {
                        if (Prism_presence[Buffer_Receive_Prism_presence_loc[iElemPrism]] == false)
                        {
                            Prism_presence[Buffer_Receive_Prism_presence_loc[iElemPrism]] = true;
                            iElem++;
                        }
                    }

                    for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++)
                    {
                        if (Pyramid_presence[Buffer_Receive_Pyramid_presence_loc[iElemPyramid]] == false)
                        {
                            Pyramid_presence[Buffer_Receive_Pyramid_presence_loc[iElemPyramid]] = true;
                            iElem++;
                        }
                    }
                }
            }

#ifdef HAVE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            /*--- iElem now contains the number of elements that this processor needs in
            total. Now we can complete the recv of the element connectivity and only
            store the elements that we need on this particular rank. Initialize space
            for the elements on this rank. ---*/

            nElem = iElem; iElem = 0;
            elem = new GRID::GRID_Primal*[nElem];
            unsigned long iElemTria = 0;
            unsigned long iElemRect = 0;
            unsigned long iElemTetr = 0;
            unsigned long iElemHexa = 0;
            unsigned long iElemPris = 0;
            unsigned long iElemPyra = 0;

            /*--- Reset presence before storing elems now that we know nElem ---*/

            for (unsigned long i = 0; i < geometry->GetnElem(); i++)
            {
                Element_presence[i] = false;
                Triangle_presence[i] = false;
                Rectangle_presence[i] = false;
                Tetrahedron_presence[i] = false;
                Hexahedron_presence[i] = false;
                Prism_presence[i] = false;
                Pyramid_presence[i] = false;
            }

            /*--- Now recv all of the element connectivity data ---*/
            for (iDomain = 0; iDomain < size; iDomain++)
            {
                if (rank != iDomain)
                {
#ifdef HAVE_MPI

                    /*--- Allocate memory for the element recv ---*/

                    Buffer_Receive_Triangle = new unsigned long[nElemTriangle_r[iDomain] * N_POINTS_TRIANGLE];
                    Buffer_Receive_Rectangle = new unsigned long[nElemRectangle_r[iDomain] * N_POINTS_QUADRILATERAL];
                    Buffer_Receive_Tetrahedron = new unsigned long[nElemTetrahedron_r[iDomain] * N_POINTS_TETRAHEDRON];
                    Buffer_Receive_Hexahedron = new unsigned long[nElemHexahedron_r[iDomain] * N_POINTS_HEXAHEDRON];
                    Buffer_Receive_Prism = new unsigned long[nElemPrism_r[iDomain] * N_POINTS_PRISM];
                    Buffer_Receive_Pyramid = new unsigned long[nElemPyramid_r[iDomain] * N_POINTS_PYRAMID];
                    Buffer_Receive_GlobElem = new unsigned long[nElemTotal_r[iDomain]];

                    /*--- Recv the element data ---*/

                    MPI_Probe(iDomain, rank * 16 + 3, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Triangle, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 3, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 4, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Rectangle, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 4, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 5, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Tetrahedron, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 5, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 6, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Hexahedron, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 6, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 7, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Prism, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 7, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 8, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_Pyramid, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 8, MPI_COMM_WORLD, &status2);

                    MPI_Probe(iDomain, rank * 16 + 9, MPI_COMM_WORLD, &status2);
                    source = status2.MPI_SOURCE;
                    MPI_Get_count(&status2, MPI_UNSIGNED_LONG, &recv_count);
                    MPI_Recv(Buffer_Receive_GlobElem, recv_count, MPI_UNSIGNED_LONG,
                        source, rank * 16 + 9, MPI_COMM_WORLD, &status2);

                    /*--- Wait to complete the above sends ---*/

                    //if (rank!=iDomain)  MPI_Waitall(7, &send_req[3], &send_stat[3]);
                    //cout << " ==== Rank " << rank << " recv from " << iDomain << " would be waiting here... " << endl;

                    /*--- Allocating the elements after the recv. Note that here we are
                    reusing the presence arrays to make sure that we find the exact same
                    set of elements that were counted above to get nElem. ---*/

                    for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++) {
                        if (Triangle_presence[Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]] == false) {
                            Triangle_presence[Buffer_Receive_Triangle_presence[iDomain][iElemTriangle]] = true;
                            elem[iElem] = new GRID::GRID_Triangle(Global_to_local_Point_recv[Buffer_Receive_Triangle[iElemTriangle * 3 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Triangle[iElemTriangle * 3 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Triangle[iElemTriangle * 3 + 2]], 2);
                            iElem++; iElemTria++;
                        }
                    }

                    for (iElemRectangle = 0; iElemRectangle < nElemRectangle_r[iDomain]; iElemRectangle++) {
                        if (Rectangle_presence[Buffer_Receive_Rectangle_presence[iDomain][iElemRectangle]] == false) {
                            Rectangle_presence[Buffer_Receive_Rectangle_presence[iDomain][iElemRectangle]] = true;
                            elem[iElem] = new GRID::GRID_Rectangle(Global_to_local_Point_recv[Buffer_Receive_Rectangle[iElemRectangle * 4 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Rectangle[iElemRectangle * 4 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Rectangle[iElemRectangle * 4 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Rectangle[iElemRectangle * 4 + 3]], 2);
                            iElem++; iElemRect++;
                        }
                    }

                    for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++) {
                        if (Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]] == false) {
                            Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence[iDomain][iElemTetrahedron]] = true;
                            elem[iElem] = new GRID::GRID_Tetrahedron(Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Tetrahedron[iElemTetrahedron * 4 + 3]]);
                            iElem++; iElemTetr++;
                        }
                    }

                    for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++) {
                        if (Hexahedron_presence[Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]] == false) {
                            Hexahedron_presence[Buffer_Receive_Hexahedron_presence[iDomain][iElemHexahedron]] = true;
                            elem[iElem] = new GRID::GRID_Hexahedron(Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 3]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 4]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 5]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 6]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron[iElemHexahedron * 8 + 7]]);
                            iElem++; iElemHexa++;
                        }
                    }

                    for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++) {
                        if (Prism_presence[Buffer_Receive_Prism_presence[iDomain][iElemPrism]] == false) {
                            Prism_presence[Buffer_Receive_Prism_presence[iDomain][iElemPrism]] = true;
                            elem[iElem] = new GRID::GRID_Prism(Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism * 6 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism * 6 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism * 6 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism * 6 + 3]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism * 6 + 4]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism[iElemPrism * 6 + 5]]);
                            iElem++; iElemPris++;
                        }
                    }

                    for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++) {
                        if (Pyramid_presence[Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]] == false) {
                            Pyramid_presence[Buffer_Receive_Pyramid_presence[iDomain][iElemPyramid]] = true;
                            elem[iElem] = new GRID::GRID_Pyramid(Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid * 5 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid * 5 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid * 5 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid * 5 + 3]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid[iElemPyramid * 5 + 4]]);
                            iElem++; iElemPyra++;
                        }
                    }

                    /*--- Free memory for the element data --*/

                    delete[] Buffer_Receive_Triangle;
                    delete[] Buffer_Receive_Rectangle;
                    delete[] Buffer_Receive_Tetrahedron;
                    delete[] Buffer_Receive_Hexahedron;
                    delete[] Buffer_Receive_Prism;
                    delete[] Buffer_Receive_Pyramid;

                    delete[] Buffer_Receive_Triangle_presence[iDomain];
                    delete[] Buffer_Receive_Rectangle_presence[iDomain];
                    delete[] Buffer_Receive_Tetrahedron_presence[iDomain];
                    delete[] Buffer_Receive_Hexahedron_presence[iDomain];
                    delete[] Buffer_Receive_Prism_presence[iDomain];
                    delete[] Buffer_Receive_Pyramid_presence[iDomain];

#endif

                }
                else
                {
                    /*--- Store the element data from our local rank ---*/
                    for (iElemTriangle = 0; iElemTriangle < nElemTriangle_r[iDomain]; iElemTriangle++)
                    {
                        if (Triangle_presence[Buffer_Receive_Triangle_presence_loc[iElemTriangle]] == false)
                        {
                            Triangle_presence[Buffer_Receive_Triangle_presence_loc[iElemTriangle]] = true;
                            elem[iElem] = new GRID::GRID_Triangle(Global_to_local_Point_recv[Buffer_Receive_Triangle_loc[iElemTriangle * 3 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Triangle_loc[iElemTriangle * 3 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Triangle_loc[iElemTriangle * 3 + 2]], 2);
                            iElem++; iElemTria++;
                        }
                    }

                    for (iElemRectangle = 0; iElemRectangle < nElemRectangle_r[iDomain]; iElemRectangle++)
                    {
                        if (Rectangle_presence[Buffer_Receive_Rectangle_presence_loc[iElemRectangle]] == false)
                        {
                            Rectangle_presence[Buffer_Receive_Rectangle_presence_loc[iElemRectangle]] = true;
                            elem[iElem] = new GRID::GRID_Rectangle(Global_to_local_Point_recv[Buffer_Receive_Rectangle_loc[iElemRectangle * 4 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Rectangle_loc[iElemRectangle * 4 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Rectangle_loc[iElemRectangle * 4 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Rectangle_loc[iElemRectangle * 4 + 3]], 2);
                            iElem++; iElemRect++;
                        }
                    }

                    for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron_r[iDomain]; iElemTetrahedron++)
                    {
                        if (Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]] == false)
                        {
                            Tetrahedron_presence[Buffer_Receive_Tetrahedron_presence_loc[iElemTetrahedron]] = true;
                            elem[iElem] = new GRID::GRID_Tetrahedron(Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron * 4 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron * 4 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron * 4 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Tetrahedron_loc[iElemTetrahedron * 4 + 3]]);
                            iElem++; iElemTetr++;
                        }
                    }

                    for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron_r[iDomain]; iElemHexahedron++)
                    {
                        if (Hexahedron_presence[Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]] == false)
                        {
                            Hexahedron_presence[Buffer_Receive_Hexahedron_presence_loc[iElemHexahedron]] = true;
                            elem[iElem] = new GRID::GRID_Hexahedron(Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 3]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 4]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 5]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 6]],
                                Global_to_local_Point_recv[Buffer_Receive_Hexahedron_loc[iElemHexahedron * 8 + 7]]);
                            iElem++; iElemHexa++;
                        }
                    }

                    for (iElemPrism = 0; iElemPrism < nElemPrism_r[iDomain]; iElemPrism++)
                    {
                        if (Prism_presence[Buffer_Receive_Prism_presence_loc[iElemPrism]] == false)
                        {
                            Prism_presence[Buffer_Receive_Prism_presence_loc[iElemPrism]] = true;
                            elem[iElem] = new GRID::GRID_Prism(Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism * 6 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism * 6 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism * 6 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism * 6 + 3]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism * 6 + 4]],
                                Global_to_local_Point_recv[Buffer_Receive_Prism_loc[iElemPrism * 6 + 5]]);
                            iElem++; iElemPris++;
                        }
                    }

                    for (iElemPyramid = 0; iElemPyramid < nElemPyramid_r[iDomain]; iElemPyramid++)
                    {
                        if (Pyramid_presence[Buffer_Receive_Pyramid_presence_loc[iElemPyramid]] == false)
                        {
                            Pyramid_presence[Buffer_Receive_Pyramid_presence_loc[iElemPyramid]] = true;
                            elem[iElem] = new GRID::GRID_Pyramid(Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid * 5 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid * 5 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid * 5 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid * 5 + 3]],
                                Global_to_local_Point_recv[Buffer_Receive_Pyramid_loc[iElemPyramid * 5 + 4]]);
                            iElem++; iElemPyra++;
                        }
                    }

                    /*--- Free memory for element data ---*/

                    delete[] Buffer_Receive_Triangle_loc;
                    delete[] Buffer_Receive_Rectangle_loc;
                    delete[] Buffer_Receive_Tetrahedron_loc;
                    delete[] Buffer_Receive_Hexahedron_loc;
                    delete[] Buffer_Receive_Prism_loc;
                    delete[] Buffer_Receive_Pyramid_loc;

                    delete[] Buffer_Receive_Triangle_presence_loc;
                    delete[] Buffer_Receive_Rectangle_presence_loc;
                    delete[] Buffer_Receive_Tetrahedron_presence_loc;
                    delete[] Buffer_Receive_Hexahedron_presence_loc;
                    delete[] Buffer_Receive_Prism_presence_loc;
                    delete[] Buffer_Receive_Pyramid_presence_loc;
                }
            }

#ifdef HAVE_MPI
            for (iDomain = 0; iDomain < size; iDomain++) {
                if (rank != iDomain) MPI_Waitall(16, send_req, send_stat);
            }
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            /*--- Free all of the memory used for communicating points and elements ---*/
            delete[] Buffer_Send_Coord;
            delete[] Buffer_Send_GlobalPointIndex;
            delete[] Buffer_Send_Color;
            delete[] Buffer_Send_Triangle;
            delete[] Buffer_Send_Rectangle;
            delete[] Buffer_Send_Tetrahedron;
            delete[] Buffer_Send_Hexahedron;
            delete[] Buffer_Send_Prism;
            delete[] Buffer_Send_Pyramid;
            delete[] Buffer_Send_BoundLine;
            delete[] Buffer_Send_BoundTriangle;
            delete[] Buffer_Send_BoundRectangle;
            delete[] Buffer_Send_Local2Global_Marker;

            delete[] Buffer_Send_SendDomain_Periodic;
            delete[] Buffer_Send_SendDomain_PeriodicTrans;
            delete[] Buffer_Send_SendDomain_PeriodicReceptor;
            delete[] Buffer_Send_ReceivedDomain_Periodic;
            delete[] Buffer_Send_ReceivedDomain_PeriodicTrans;
            delete[] Buffer_Send_ReceivedDomain_PeriodicDonor;

            delete[] Local_to_global_Triangle;
            delete[] Local_to_global_Rectangle;
            delete[] Local_to_global_Tetrahedron;
            delete[] Local_to_global_Hexahedron;
            delete[] Local_to_global_Prism;
            delete[] Local_to_global_Pyramid;


            /*--- Communicate the number of each element type to all processors. These
            values are important for merging and writing output later. ---*/

#ifdef HAVE_MPI
            unsigned long Local_nElem = nElem;
            MPI_Allreduce(&Local_nElem, &Global_nElem, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
            Global_nElem = nElem;
#endif

            if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
                std::cout << Global_nElem << " interior elements including halo cells. " << std::endl;

            /*--- Store total number of each element type after incrementing the
            counters in the recv loop above (to make sure there aren't repeats). ---*/

            nelem_triangle = iElemTria;
            nelem_quad = iElemRect;
            nelem_tetra = iElemTetr;
            nelem_hexa = iElemHexa;
            nelem_prism = iElemPris;
            nelem_pyramid = iElemPyra;

#ifdef HAVE_MPI
            unsigned long Local_nElemTri = nelem_triangle;
            unsigned long Local_nElemQuad = nelem_quad;
            unsigned long Local_nElemTet = nelem_tetra;
            unsigned long Local_nElemHex = nelem_hexa;
            unsigned long Local_nElemPrism = nelem_prism;
            unsigned long Local_nElemPyramid = nelem_pyramid;
            MPI_Allreduce(&Local_nElemTri, &Global_nelem_triangle, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nElemQuad, &Global_nelem_quad, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nElemTet, &Global_nelem_tetra, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nElemHex, &Global_nelem_hexa, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nElemPrism, &Global_nelem_prism, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nElemPyramid, &Global_nelem_pyramid, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
            Global_nelem_triangle = nelem_triangle;
            Global_nelem_quad = nelem_quad;
            Global_nelem_tetra = nelem_tetra;
            Global_nelem_hexa = nelem_hexa;
            Global_nelem_prism = nelem_prism;
            Global_nelem_pyramid = nelem_pyramid;
#endif

            /*--- Print information about the elements to the console ---*/

            if (rank == TBOX::MASTER_NODE)
            {
                if (Global_nelem_triangle > 0)  std::cout << Global_nelem_triangle << " triangles." << std::endl;
                if (Global_nelem_quad > 0)      std::cout << Global_nelem_quad << " quadrilaterals." << std::endl;
                if (Global_nelem_tetra > 0)     std::cout << Global_nelem_tetra << " tetrahedra." << std::endl;
                if (Global_nelem_hexa > 0)      std::cout << Global_nelem_hexa << " hexahedra." << std::endl;
                if (Global_nelem_prism > 0)     std::cout << Global_nelem_prism << " prisms." << std::endl;
                if (Global_nelem_pyramid > 0)   std::cout << Global_nelem_pyramid << " pyramids." << std::endl;
            }

            delete[] Triangle_presence;
            delete[] Rectangle_presence;
            delete[] Tetrahedron_presence;
            delete[] Hexahedron_presence;
            delete[] Prism_presence;
            delete[] Pyramid_presence;

            /*--- Now partition the boundary elements on the markers. Note that, for
            now, we are still performing the boundary partitioning using the master
            node alone. The boundaries should make up a much smaller portion of the
            mesh, so this is ok for now, but we will transition to a parallel version
            of this soon that follows the same procedure above for the interior. ---*/

            if (rank == TBOX::MASTER_NODE)
            {
                /*--- Create auxiliary std::vectors based on the original geometry ---*/

                MarkerIn = new bool[geometry->GetnMarker()];
                VertexIn = new bool*[geometry->GetnMarker()];

                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    VertexIn[iMarker] = new bool[geometry->GetnElem_Bound(iMarker)];

                Buffer_Send_nDim = geometry->GetnDim();
                Buffer_Send_nZone = geometry->GetnZone();
                Buffer_Send_nPeriodic = config->GetnPeriodicIndex();
                Buffer_Send_Center = new double[Buffer_Send_nPeriodic * 3];
                Buffer_Send_Rotation = new double[Buffer_Send_nPeriodic * 3];
                Buffer_Send_Translate = new double[Buffer_Send_nPeriodic * 3];

                Buffer_Send_nSendDomain_Periodic = new unsigned long[nDomain];
                Buffer_Send_nReceivedDomain_Periodic = new unsigned long[nDomain];

                /*--- Create a local copy of config->GetMarker_All_SendRecv and
                config->GetMarker_All_TagBound in the master node ---*/

                Marker_All_SendRecv_Copy = new short[geometry->GetnMarker()];
                Marker_All_TagBound_Copy = new std::string[geometry->GetnMarker()];

                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                {
                    Marker_All_SendRecv_Copy[iMarker] = config->GetMarker_All_SendRecv(iMarker);
                    Marker_All_TagBound_Copy[iMarker] = config->GetMarker_All_TagBound(iMarker);
                }

            }

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {

                if (rank == TBOX::MASTER_NODE)
                {

                    /*--- Interior dimensionalization. Loop over the original grid
                    to perform the dimensionalizaton of the domain variables ---*/
                    Buffer_Send_nElemTotal = 0;
                    Buffer_Send_nPointTotal = 0;
                    Buffer_Send_nPointGhost = 0;
                    Buffer_Send_nPointDomainTotal = 0;
                    Buffer_Send_nPointPeriodic = 0;
                    Buffer_Send_nElemTriangle = 0;
                    Buffer_Send_nElemRectangle = 0;
                    Buffer_Send_nElemTetrahedron = 0;
                    Buffer_Send_nElemHexahedron = 0;
                    Buffer_Send_nElemPrism = 0;
                    Buffer_Send_nElemPyramid = 0;

                    /*--- Boundary dimensionalization. Dimensionalization with physical
                    boundaries, compute Buffer_Send_nMarkerDomain,
                    Buffer_Send_nVertexDomain[nMarkerDomain] ---*/

                    Buffer_Send_nMarkerDomain = 0;
                    Buffer_Send_nBoundLineTotal = 0;
                    Buffer_Send_nBoundTriangleTotal = 0;
                    Buffer_Send_nBoundRectangleTotal = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        Buffer_Send_nVertexDomain[iMarker] = 0;
                        Buffer_Send_nBoundLine[iMarker] = 0;
                        Buffer_Send_nBoundTriangle[iMarker] = 0;
                        Buffer_Send_nBoundRectangle[iMarker] = 0;
                        Buffer_Send_Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
                        sprintf(&Buffer_Send_Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE], "%s",
                            Marker_All_TagBound_Copy[iMarker].c_str());
                    }

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
                        {

                            MarkerIn[iMarker] = false;
                            Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain] = 0;

                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                VertexIn[iMarker][iVertex] = false;
                                for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++)
                                {
                                    iPoint = geometry->bound[iMarker][iVertex]->GetNode(iNode);
                                    if (local_colour_values[iPoint] == iDomain) VertexIn[iMarker][iVertex] = true;
                                }

                                /*--- If this vertex should be sent, increment the element type ---*/
                                if (VertexIn[iMarker][iVertex])
                                {
                                    switch (geometry->bound[iMarker][iVertex]->GetVTK_Type())
                                    {
                                    case TBOX::LINE:
                                        Buffer_Send_nBoundLine[Buffer_Send_nMarkerDomain]++;
                                        Buffer_Send_nBoundLineTotal++;
                                        break;
                                    case TBOX::TRIANGLE:
                                        Buffer_Send_nBoundTriangle[Buffer_Send_nMarkerDomain]++;
                                        Buffer_Send_nBoundTriangleTotal++;
                                        break;
                                    case TBOX::RECTANGLE:
                                        Buffer_Send_nBoundRectangle[Buffer_Send_nMarkerDomain]++;
                                        Buffer_Send_nBoundRectangleTotal++;
                                        break;
                                    }

                                    /*--- Increment the total number of vertices to be sent ---*/
                                    Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain]++;
                                    MarkerIn[iMarker] = true;
                                }
                            }

                            /*--- Increment the number of markers to be sent ---*/
                            if (MarkerIn[iMarker])
                            {
                                Buffer_Send_nMarkerDomain++;
                            }
                        }
                    }

                    /*--- Copy periodic information from the config file ---*/
                    for (iPeriodic = 0; iPeriodic < Buffer_Send_nPeriodic; iPeriodic++)
                    {
                        for (iDim = 0; iDim < 3; iDim++)
                        {
                            Buffer_Send_Center[iDim + iPeriodic * 3] = config->GetPeriodicCenter(iPeriodic)[iDim];
                            Buffer_Send_Rotation[iDim + iPeriodic * 3] = config->GetPeriodicRotation(iPeriodic)[iDim];
                            Buffer_Send_Translate[iDim + iPeriodic * 3] = config->GetPeriodicTranslate(iPeriodic)[iDim];
                        }
                    }

                    /*--- Dimensionalization of the periodic auxiliary std::vectors ---*/

                    for (jDomain = 0; jDomain < nDomain; jDomain++)
                    {
                        Buffer_Send_nSendDomain_Periodic[jDomain] = 0;
                        Buffer_Send_nReceivedDomain_Periodic[jDomain] = 0;
                    }
                    Buffer_Send_nTotalSendDomain_Periodic = 0;
                    Buffer_Send_nTotalReceivedDomain_Periodic = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE)
                        {
                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
                                if (iDomain == local_colour_values[iPoint])
                                {
                                    if (config->GetMarker_All_SendRecv(iMarker) > 0)
                                    {

                                        /*--- Identify the color of the receptor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                ReceptorColor = local_colour_values[jPoint];
                                            }
                                        }

                                        Buffer_Send_nSendDomain_Periodic[ReceptorColor]++;
                                        Buffer_Send_nTotalSendDomain_Periodic++;

                                    }
                                    if (config->GetMarker_All_SendRecv(iMarker) < 0)
                                    {

                                        /*--- Identify the color of the donor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                DonorColor = local_colour_values[jPoint];
                                            }
                                        }

                                        Buffer_Send_nReceivedDomain_Periodic[DonorColor]++;
                                        Buffer_Send_nTotalReceivedDomain_Periodic++;
                                    }
                                }
                            }
                        }
                    }

                    /*--- Allocate the buffer std::vectors in the appropiate domain (master, iDomain) ---*/

                    Buffer_Send_BoundLine = new unsigned long[Buffer_Send_nBoundLineTotal*TBOX::N_POINTS_LINE];
                    Buffer_Send_BoundTriangle = new unsigned long[Buffer_Send_nBoundTriangleTotal*TBOX::N_POINTS_TRIANGLE];
                    Buffer_Send_BoundRectangle = new unsigned long[Buffer_Send_nBoundRectangleTotal*TBOX::N_POINTS_QUADRILATERAL];
                    Buffer_Send_Local2Global_Marker = new unsigned long[Buffer_Send_nMarkerDomain];

                    Buffer_Send_SendDomain_Periodic = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
                    Buffer_Send_SendDomain_PeriodicTrans = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
                    Buffer_Send_SendDomain_PeriodicReceptor = new unsigned long[Buffer_Send_nTotalSendDomain_Periodic];
                    Buffer_Send_ReceivedDomain_Periodic = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
                    Buffer_Send_ReceivedDomain_PeriodicTrans = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];
                    Buffer_Send_ReceivedDomain_PeriodicDonor = new unsigned long[Buffer_Send_nTotalReceivedDomain_Periodic];

                    if (iDomain != TBOX::MASTER_NODE)
                    {
                        //cout << " Rank " << rank << " iDomain " << iDomain << endl;

#ifdef HAVE_MPI

                        MPI_Isend(&Buffer_Send_nBoundLineTotal, 1,
                            MPI_UNSIGNED_LONG, iDomain,
                            0, MPI_COMM_WORLD, &send_req[0]);

                        MPI_Isend(&Buffer_Send_nBoundTriangleTotal, 1,
                            MPI_UNSIGNED_LONG, iDomain,
                            1, MPI_COMM_WORLD, &send_req[1]);

                        MPI_Isend(&Buffer_Send_nBoundRectangleTotal, 1,
                            MPI_UNSIGNED_LONG, iDomain,
                            2, MPI_COMM_WORLD, &send_req[2]);

                        MPI_Isend(&Buffer_Send_nMarkerDomain, 1,
                            MPI_UNSIGNED_SHORT, iDomain,
                            3, MPI_COMM_WORLD, &send_req[3]);

                        MPI_Isend(Buffer_Send_nVertexDomain,
                            nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                            4, MPI_COMM_WORLD, &send_req[4]);

                        MPI_Isend(Buffer_Send_nBoundLine,
                            nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                            5, MPI_COMM_WORLD, &send_req[5]);

                        MPI_Isend(Buffer_Send_nBoundTriangle,
                            nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                            6, MPI_COMM_WORLD, &send_req[6]);

                        MPI_Isend(Buffer_Send_nBoundRectangle,
                            nMarker_Max, MPI_UNSIGNED_LONG, iDomain,
                            7, MPI_COMM_WORLD, &send_req[7]);

                        MPI_Isend(Buffer_Send_Marker_All_SendRecv,
                            nMarker_Max, MPI_SHORT, iDomain,
                            8, MPI_COMM_WORLD, &send_req[8]);

                        MPI_Isend(Buffer_Send_Marker_All_TagBound,
                            nMarker_Max*MAX_STRING_SIZE, MPI_CHAR, iDomain,
                            9, MPI_COMM_WORLD, &send_req[9]);

                        MPI_Isend(&Buffer_Send_nPeriodic,
                            1, MPI_UNSIGNED_SHORT, iDomain,
                            10, MPI_COMM_WORLD, &send_req[10]);

                        MPI_Isend(Buffer_Send_Center,
                            nPeriodic * 3, MPI_DOUBLE, iDomain,
                            11, MPI_COMM_WORLD, &send_req[11]);

                        MPI_Isend(Buffer_Send_Rotation,
                            nPeriodic * 3, MPI_DOUBLE, iDomain,
                            12, MPI_COMM_WORLD, &send_req[12]);

                        MPI_Isend(Buffer_Send_Translate,
                            nPeriodic * 3, MPI_DOUBLE, iDomain,
                            13, MPI_COMM_WORLD, &send_req[13]);

                        MPI_Isend(&Buffer_Send_nTotalSendDomain_Periodic,
                            1, MPI_UNSIGNED_LONG, iDomain,
                            14, MPI_COMM_WORLD, &send_req[14]);

                        MPI_Isend(&Buffer_Send_nTotalReceivedDomain_Periodic,
                            1, MPI_UNSIGNED_LONG, iDomain,
                            15, MPI_COMM_WORLD, &send_req[15]);

                        MPI_Isend(Buffer_Send_nSendDomain_Periodic,
                            nDomain, MPI_UNSIGNED_LONG, iDomain,
                            16, MPI_COMM_WORLD, &send_req[16]);

                        MPI_Isend(Buffer_Send_nReceivedDomain_Periodic,
                            nDomain, MPI_UNSIGNED_LONG, iDomain,
                            17, MPI_COMM_WORLD, &send_req[17]);

                        /*--- Wait for this set of non-blocking comm. to complete ---*/


                        MPI_Waitall(18, send_req, send_stat);
                        //cout << " Rank " << rank << " iDomain " << iDomain << " just waited for first sends" << endl;

#endif

                    }
                    else
                    {

                        /*--- We are the master node, so simply copy values into place ---*/
                        nDim = Buffer_Send_nDim;
                        nZone = Buffer_Send_nZone;

                        nPeriodic = Buffer_Send_nPeriodic;
                        nPointGhost = Buffer_Send_nPointGhost;
                        nPointPeriodic = Buffer_Send_nPointPeriodic;

                        nBoundLineTotal = Buffer_Send_nBoundLineTotal;
                        nBoundTriangleTotal = Buffer_Send_nBoundTriangleTotal;
                        nBoundRectangleTotal = Buffer_Send_nBoundRectangleTotal;
                        nMarkerDomain = Buffer_Send_nMarkerDomain;

                        for (iMarker = 0; iMarker < nMarker_Max; iMarker++)
                        {
                            nVertexDomain[iMarker] = Buffer_Send_nVertexDomain[iMarker];
                            nBoundLine[iMarker] = Buffer_Send_nBoundLine[iMarker];
                            nBoundTriangle[iMarker] = Buffer_Send_nBoundTriangle[iMarker];
                            nBoundRectangle[iMarker] = Buffer_Send_nBoundRectangle[iMarker];
                            Marker_All_SendRecv[iMarker] = Buffer_Send_Marker_All_SendRecv[iMarker];
                            for (iter = 0; iter < TBOX::MAX_STRING_SIZE; iter++)
                                Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE + iter] = Buffer_Send_Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE + iter];
                        }

                        Buffer_Receive_Center = new double[nPeriodic * 3];
                        Buffer_Receive_Rotation = new double[nPeriodic * 3];
                        Buffer_Receive_Translate = new double[nPeriodic * 3];

                        for (iter = 0; iter < nPeriodic * 3; iter++)
                        {
                            Buffer_Receive_Center[iter] = Buffer_Send_Center[iter];
                            Buffer_Receive_Rotation[iter] = Buffer_Send_Rotation[iter];
                            Buffer_Receive_Translate[iter] = Buffer_Send_Translate[iter];
                        }

                        nTotalSendDomain_Periodic = Buffer_Send_nTotalSendDomain_Periodic;
                        nTotalReceivedDomain_Periodic = Buffer_Send_nTotalReceivedDomain_Periodic;

                        for (iter = 0; iter < nDomain; iter++)
                        {
                            nSendDomain_Periodic[iter] = Buffer_Send_nSendDomain_Periodic[iter];
                            nReceivedDomain_Periodic[iter] = Buffer_Send_nReceivedDomain_Periodic[iter];
                        }

                    }
                }

                /*--- Each rank now begins to receive information from the master ---*/

                if (rank == iDomain)
                {
                    /*--- First, receive the size of buffers before receiving the data ---*/
                    if (rank != TBOX::MASTER_NODE)
                    {

#ifdef HAVE_MPI

                        MPI_Probe(MASTER_NODE, 0, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(&nBoundLineTotal, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 0, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 1, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(&nBoundTriangleTotal, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 1, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 2, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(&nBoundRectangleTotal, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 2, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 3, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_SHORT, &recv_count);
                        MPI_Recv(&nMarkerDomain, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 3, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 4, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(nVertexDomain, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 4, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 5, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(nBoundLine, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 5, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 6, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(nBoundTriangle, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 6, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 7, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(nBoundRectangle, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 7, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 8, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_SHORT, &recv_count);
                        MPI_Recv(Marker_All_SendRecv, recv_count, MPI_SHORT,
                            MASTER_NODE, 8, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 9, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &recv_count);
                        MPI_Recv(Marker_All_TagBound, recv_count, MPI_CHAR,
                            MASTER_NODE, 9, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 10, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_SHORT, &recv_count);
                        MPI_Recv(&nPeriodic, recv_count, MPI_UNSIGNED_SHORT,
                            MASTER_NODE, 10, MPI_COMM_WORLD, &status);

#endif

                        /*--- Marker_All_TagBound and Marker_All_SendRecv, set the same
                        values in the config files of all the files ---*/

                        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                        {
                            config->SetMarker_All_SendRecv(iMarker,
                                Marker_All_SendRecv[iMarker]);
                            config->SetMarker_All_TagBound(iMarker,
                                std::string(&Marker_All_TagBound[iMarker*TBOX::MAX_STRING_SIZE]));
                        }


                        /*--- Periodic boundary conditions ---*/

                        Buffer_Receive_Center = new double[nPeriodic * 3];
                        Buffer_Receive_Rotation = new double[nPeriodic * 3];
                        Buffer_Receive_Translate = new double[nPeriodic * 3];

#ifdef HAVE_MPI

                        MPI_Probe(MASTER_NODE, 11, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
                        MPI_Recv(Buffer_Receive_Center, recv_count, MPI_DOUBLE,
                            MASTER_NODE, 11, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 12, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
                        MPI_Recv(Buffer_Receive_Rotation, recv_count, MPI_DOUBLE,
                            MASTER_NODE, 12, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 13, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
                        MPI_Recv(Buffer_Receive_Translate, recv_count, MPI_DOUBLE,
                            MASTER_NODE, 13, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 14, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(&nTotalSendDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 14, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 15, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(&nTotalReceivedDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 15, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 16, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(nSendDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 16, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 17, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(nReceivedDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 17, MPI_COMM_WORLD, &status);

#endif

                        config->SetnPeriodicIndex(nPeriodic);

                        for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++)
                        {

                            double* center = new double[3];       // Do not deallocate the memory
                            double* rotation = new double[3];    // Do not deallocate the memory
                            double* translate = new double[3];    // Do not deallocate the memory

                            for (iDim = 0; iDim < 3; iDim++)
                            {
                                center[iDim] = Buffer_Receive_Center[iDim + iPeriodic * 3];
                                rotation[iDim] = Buffer_Receive_Rotation[iDim + iPeriodic * 3];
                                translate[iDim] = Buffer_Receive_Translate[iDim + iPeriodic * 3];
                            }
                            config->SetPeriodicCenter(iPeriodic, center);
                            config->SetPeriodicRotation(iPeriodic, rotation);
                            config->SetPeriodicTranslate(iPeriodic, translate);
                        }

                    }

                    delete[] Buffer_Receive_Center;
                    delete[] Buffer_Receive_Rotation;
                    delete[] Buffer_Receive_Translate;

                    /*--- Allocate the receive buffer std::vector ---*/

                    Buffer_Receive_BoundLine = new unsigned long[nBoundLineTotal * 2];
                    Buffer_Receive_BoundTriangle = new unsigned long[nBoundTriangleTotal * 3];
                    Buffer_Receive_BoundRectangle = new unsigned long[nBoundRectangleTotal * 4];
                    Buffer_Receive_Local2Global_Marker = new unsigned long[nMarkerDomain];

                    Buffer_Receive_SendDomain_Periodic = new unsigned long[nTotalSendDomain_Periodic];
                    Buffer_Receive_SendDomain_PeriodicTrans = new unsigned long[nTotalSendDomain_Periodic];
                    Buffer_Receive_SendDomain_PeriodicReceptor = new unsigned long[nTotalSendDomain_Periodic];
                    Buffer_Receive_ReceivedDomain_Periodic = new unsigned long[nTotalReceivedDomain_Periodic];
                    Buffer_Receive_ReceivedDomain_PeriodicTrans = new unsigned long[nTotalReceivedDomain_Periodic];
                    Buffer_Receive_ReceivedDomain_PeriodicDonor = new unsigned long[nTotalReceivedDomain_Periodic];

                }


                //cout << " &&&& Rank " << rank << " about to start bound elems " << endl;

                /*--- Set the value of the Send buffers ---*/

                if (rank == TBOX::MASTER_NODE)
                {

                    /*--- Set the value of the boundary geometry ---*/

                    iMarkerDomain = 0;
                    iBoundLineTotal = 0; iBoundTriangleTotal = 0; iBoundRectangleTotal = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if ((config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE) && (MarkerIn[iMarker]))
                        {
                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {

                                if (VertexIn[iMarker][iVertex])
                                {

                                    /*--- Send global index here and then convert to local on the recv ---*/

                                    for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++)
                                    {
                                        vnodes_local[iNode] = geometry->bound[iMarker][iVertex]->GetNode(iNode);
                                    }

                                    switch (geometry->bound[iMarker][iVertex]->GetVTK_Type())
                                    {
                                    case TBOX::LINE:
                                        Buffer_Send_BoundLine[TBOX::N_POINTS_LINE*iBoundLineTotal + 0] = vnodes_local[0];
                                        Buffer_Send_BoundLine[TBOX::N_POINTS_LINE*iBoundLineTotal + 1] = vnodes_local[1];
                                        iBoundLineTotal++;
                                        break;
                                    case TBOX::TRIANGLE:
                                        Buffer_Send_BoundTriangle[TBOX::N_POINTS_TRIANGLE*iBoundTriangleTotal + 0] = vnodes_local[0];
                                        Buffer_Send_BoundTriangle[TBOX::N_POINTS_TRIANGLE*iBoundTriangleTotal + 1] = vnodes_local[1];
                                        Buffer_Send_BoundTriangle[TBOX::N_POINTS_TRIANGLE*iBoundTriangleTotal + 2] = vnodes_local[2];
                                        iBoundTriangleTotal++;
                                        break;
                                    case TBOX::RECTANGLE:
                                        Buffer_Send_BoundRectangle[TBOX::N_POINTS_QUADRILATERAL*iBoundRectangleTotal + 0] = vnodes_local[0];
                                        Buffer_Send_BoundRectangle[TBOX::N_POINTS_QUADRILATERAL*iBoundRectangleTotal + 1] = vnodes_local[1];
                                        Buffer_Send_BoundRectangle[TBOX::N_POINTS_QUADRILATERAL*iBoundRectangleTotal + 2] = vnodes_local[2];
                                        Buffer_Send_BoundRectangle[TBOX::N_POINTS_QUADRILATERAL*iBoundRectangleTotal + 3] = vnodes_local[3];
                                        iBoundRectangleTotal++;
                                        break;
                                    }
                                }
                            }

                            Buffer_Send_Local2Global_Marker[iMarkerDomain] = iMarker;
                            iMarkerDomain++;

                        }
                    }

                    /*--- Evaluate the number of already existing periodic boundary conditions ---*/

                    iTotalSendDomain_Periodic = 0;
                    iTotalReceivedDomain_Periodic = 0;

                    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE)
                        {
                            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++)
                            {
                                iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
                                Transformation = geometry->bound[iMarker][iVertex]->GetRotation_Type();

                                if (iDomain == local_colour_values[iPoint])
                                {

                                    /*--- If the information is going to be sended, find the
                                    domain of the receptor ---*/

                                    if (config->GetMarker_All_SendRecv(iMarker) > 0)
                                    {

                                        /*--- Identify the color of the receptor ---*/

                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                ReceptorColor = local_colour_values[jPoint];
                                            }
                                        }

                                        /*--- For each color of the receptor we will han an extra marker (+) ---*/

                                        Buffer_Send_SendDomain_Periodic[iTotalSendDomain_Periodic] = iPoint;
                                        Buffer_Send_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic] = Transformation;
                                        Buffer_Send_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] = ReceptorColor;

                                        iTotalSendDomain_Periodic++;

                                    }

                                    /*--- If the information is goint to be received, find the domain if the donor ---*/

                                    if (config->GetMarker_All_SendRecv(iMarker) < 0)
                                    {

                                        /*--- Identify the color of the donor ---*/
                                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++)
                                        {
                                            if ((config->GetMarker_All_KindBC(jMarker) == TBOX::SEND_RECEIVE) &&
                                                (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker)))
                                            {
                                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                                DonorColor = local_colour_values[jPoint];
                                            }
                                        }

                                        /*--- For each color of the donor we will han an extra marker (-) ---*/

                                        Buffer_Send_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic] = iPoint;
                                        Buffer_Send_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic] = Transformation;
                                        Buffer_Send_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] = DonorColor;

                                        iTotalReceivedDomain_Periodic++;

                                    }
                                }
                            }
                        }
                    }

                    /*--- Send the buffers with the geometrical information ---*/

                    if (iDomain != TBOX::MASTER_NODE)
                    {

#ifdef HAVE_MPI

                        MPI_Isend(Buffer_Send_BoundLine,
                            Buffer_Send_nBoundLineTotal*N_POINTS_LINE, MPI_UNSIGNED_LONG, iDomain,
                            0, MPI_COMM_WORLD, &send_req[0]);

                        MPI_Isend(Buffer_Send_BoundTriangle,
                            Buffer_Send_nBoundTriangleTotal*N_POINTS_TRIANGLE, MPI_UNSIGNED_LONG, iDomain,
                            1, MPI_COMM_WORLD, &send_req[1]);

                        MPI_Isend(Buffer_Send_BoundRectangle,
                            Buffer_Send_nBoundRectangleTotal*N_POINTS_QUADRILATERAL, MPI_UNSIGNED_LONG, iDomain,
                            2, MPI_COMM_WORLD, &send_req[2]);

                        MPI_Isend(Buffer_Send_Local2Global_Marker,
                            Buffer_Send_nMarkerDomain, MPI_UNSIGNED_LONG, iDomain,
                            3, MPI_COMM_WORLD, &send_req[3]);

                        MPI_Isend(Buffer_Send_SendDomain_Periodic,
                            Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                            4, MPI_COMM_WORLD, &send_req[4]);

                        MPI_Isend(Buffer_Send_SendDomain_PeriodicTrans,
                            Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                            5, MPI_COMM_WORLD, &send_req[5]);

                        MPI_Isend(Buffer_Send_SendDomain_PeriodicReceptor,
                            Buffer_Send_nTotalSendDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                            6, MPI_COMM_WORLD, &send_req[6]);

                        MPI_Isend(Buffer_Send_ReceivedDomain_Periodic,
                            Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                            7, MPI_COMM_WORLD, &send_req[7]);

                        MPI_Isend(Buffer_Send_ReceivedDomain_PeriodicTrans,
                            Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                            8, MPI_COMM_WORLD, &send_req[8]);

                        MPI_Isend(Buffer_Send_ReceivedDomain_PeriodicDonor,
                            Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain,
                            9, MPI_COMM_WORLD, &send_req[9]);

                        /*--- Wait for this set of non-blocking comm. to complete ---*/

                        MPI_Waitall(10, send_req, send_stat);

#endif

                    }
                    else
                    {

                        /*--- Copy the data directly from our own rank ---*/
                        for (iter = 0; iter < Buffer_Send_nBoundLineTotal*TBOX::N_POINTS_LINE; iter++)
                            Buffer_Receive_BoundLine[iter] = Buffer_Send_BoundLine[iter];

                        for (iter = 0; iter < Buffer_Send_nBoundTriangleTotal*TBOX::N_POINTS_TRIANGLE; iter++)
                            Buffer_Receive_BoundTriangle[iter] = Buffer_Send_BoundTriangle[iter];

                        for (iter = 0; iter < Buffer_Send_nBoundRectangleTotal*TBOX::N_POINTS_QUADRILATERAL; iter++)
                            Buffer_Receive_BoundRectangle[iter] = Buffer_Send_BoundRectangle[iter];

                        for (iter = 0; iter < Buffer_Send_nMarkerDomain; iter++)
                            Buffer_Receive_Local2Global_Marker[iter] = Buffer_Send_Local2Global_Marker[iter];

                        for (iter = 0; iter < Buffer_Send_nTotalSendDomain_Periodic; iter++)
                        {
                            Buffer_Receive_SendDomain_Periodic[iter] = Buffer_Send_SendDomain_Periodic[iter];
                            Buffer_Receive_SendDomain_PeriodicTrans[iter] = Buffer_Send_SendDomain_PeriodicTrans[iter];
                            Buffer_Receive_SendDomain_PeriodicReceptor[iter] = Buffer_Send_SendDomain_PeriodicReceptor[iter];
                        }

                        for (iter = 0; iter < Buffer_Send_nTotalReceivedDomain_Periodic; iter++)
                        {
                            Buffer_Receive_ReceivedDomain_Periodic[iter] = Buffer_Send_ReceivedDomain_Periodic[iter];
                            Buffer_Receive_ReceivedDomain_PeriodicTrans[iter] = Buffer_Send_ReceivedDomain_PeriodicTrans[iter];
                            Buffer_Receive_ReceivedDomain_PeriodicDonor[iter] = Buffer_Send_ReceivedDomain_PeriodicDonor[iter];
                        }

                    }

                    delete[] Buffer_Send_BoundLine;
                    delete[] Buffer_Send_BoundTriangle;
                    delete[] Buffer_Send_BoundRectangle;
                    delete[] Buffer_Send_Local2Global_Marker;

                    delete[] Buffer_Send_SendDomain_Periodic;
                    delete[] Buffer_Send_SendDomain_PeriodicTrans;
                    delete[] Buffer_Send_SendDomain_PeriodicReceptor;
                    delete[] Buffer_Send_ReceivedDomain_Periodic;
                    delete[] Buffer_Send_ReceivedDomain_PeriodicTrans;
                    delete[] Buffer_Send_ReceivedDomain_PeriodicDonor;

                }

                //cout << " Rank " << rank << " about to recv of bound elems " << endl;

                if (rank == iDomain)
                {

                    if (rank != TBOX::MASTER_NODE)
                    {

                        /*--- Receive the buffers with the geometrical information ---*/
#ifdef HAVE_MPI

                        MPI_Probe(MASTER_NODE, 0, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_BoundLine, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 0, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 1, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_BoundTriangle, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 1, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 2, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_BoundRectangle, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 2, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 3, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_Local2Global_Marker, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 3, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 4, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_SendDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 4, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 5, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_SendDomain_PeriodicTrans, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 5, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 6, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_SendDomain_PeriodicReceptor, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 6, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 7, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_ReceivedDomain_Periodic, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 7, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 8, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_ReceivedDomain_PeriodicTrans, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 8, MPI_COMM_WORLD, &status);

                        MPI_Probe(MASTER_NODE, 9, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &recv_count);
                        MPI_Recv(Buffer_Receive_ReceivedDomain_PeriodicDonor, recv_count, MPI_UNSIGNED_LONG,
                            MASTER_NODE, 9, MPI_COMM_WORLD, &status);

#endif

                    }

                    /*--- Create the domain structures for the boundaries ---*/

                    nMarker = nMarkerDomain;
                    nElem_Bound = new unsigned long[nMarker_Max];
                    Local_to_Global_Marker = new unsigned short[nMarker_Max];
                    Tag_to_Marker = new std::string[nMarker_Max];
                    std::string *TagBound_Copy = new std::string[nMarker_Max];
                    short *SendRecv_Copy = new short[nMarker_Max];

                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                        nElem_Bound[iMarker] = nVertexDomain[iMarker];

                    bound = new GRID::GRID_Primal**[nMarker + (overhead*nDomain)];
                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                        bound[iMarker] = new GRID::GRID_Primal*[nElem_Bound[iMarker]];

                    /*--- Initialize boundary element counters ---*/
                    iBoundLineTotal = 0;
                    iBoundTriangleTotal = 0;
                    iBoundRectangleTotal = 0;

                    /*--- Store the boundary element connectivity. Note here that we have
                    communicated the global index values for the elements, so we need to
                    convert this to the local index when instantiating the element. ---*/

                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                    {

                        iVertexDomain = 0;

                        for (iBoundLine = 0; iBoundLine < nBoundLine[iMarker]; iBoundLine++)
                        {
                            bound[iMarker][iVertexDomain] = new GRID::GRID_Line(Global_to_local_Point_recv[Buffer_Receive_BoundLine[iBoundLineTotal * 2 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_BoundLine[iBoundLineTotal * 2 + 1]], 2);
                            iVertexDomain++; iBoundLineTotal++;
                        }
                        for (iBoundTriangle = 0; iBoundTriangle < nBoundTriangle[iMarker]; iBoundTriangle++)
                        {
                            bound[iMarker][iVertexDomain] = new GRID::GRID_Triangle(Global_to_local_Point_recv[Buffer_Receive_BoundTriangle[iBoundTriangleTotal * 3 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_BoundTriangle[iBoundTriangleTotal * 3 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_BoundTriangle[iBoundTriangleTotal * 3 + 2]], 3);
                            iVertexDomain++; iBoundTriangleTotal++;
                        }
                        for (iBoundRectangle = 0; iBoundRectangle < nBoundRectangle[iMarker]; iBoundRectangle++)
                        {
                            bound[iMarker][iVertexDomain] = new GRID::GRID_Rectangle(Global_to_local_Point_recv[Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 0]],
                                Global_to_local_Point_recv[Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 1]],
                                Global_to_local_Point_recv[Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 2]],
                                Global_to_local_Point_recv[Buffer_Receive_BoundRectangle[iBoundRectangleTotal * 4 + 3]], 3);
                            iVertexDomain++; iBoundRectangleTotal++;
                        }

                        Local_to_Global_Marker[iMarker] = Buffer_Receive_Local2Global_Marker[iMarker];

                        /*--- Now each domain has the right information ---*/

                        std::string Grid_Marker = config->GetMarker_All_TagBound(Local_to_Global_Marker[iMarker]);
                        short SendRecv = config->GetMarker_All_SendRecv(Local_to_Global_Marker[iMarker]);
                        TagBound_Copy[iMarker] = Grid_Marker;
                        SendRecv_Copy[iMarker] = SendRecv;

                    }

                    /*--- Store total number of each boundary element type ---*/

                    nelem_edge_bound = iBoundLineTotal;
                    nelem_triangle_bound = iBoundTriangleTotal;
                    nelem_quad_bound = iBoundRectangleTotal;

                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                    {
                        config->SetMarker_All_TagBound(iMarker, TagBound_Copy[iMarker]);
                        config->SetMarker_All_SendRecv(iMarker, SendRecv_Copy[iMarker]);
                    }

                    /*--- Add the new periodic markers to the domain ---*/

                    iTotalSendDomain_Periodic = 0;
                    iTotalReceivedDomain_Periodic = 0;

                    for (jDomain = 0; jDomain < nDomain; jDomain++)
                    {

                        if (nSendDomain_Periodic[jDomain] != 0)
                        {
                            nVertexDomain[nMarker] = 0;
                            bound[nMarker] = new GRID::GRID_Primal*[nSendDomain_Periodic[jDomain]];

                            iVertex = 0;
                            for (iTotalSendDomain_Periodic = 0; iTotalSendDomain_Periodic < nTotalSendDomain_Periodic; iTotalSendDomain_Periodic++)
                            {
                                if (Buffer_Receive_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] == jDomain)
                                {
                                    bound[nMarker][iVertex] = new GRID::GRID_VertexMPI(Global_to_local_Point_recv[Buffer_Receive_SendDomain_Periodic[iTotalSendDomain_Periodic]], nDim);
                                    bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic]);
                                    nVertexDomain[nMarker]++; iVertex++;
                                }
                            }

                            Marker_All_SendRecv[nMarker] = jDomain + 1;
                            nElem_Bound[nMarker] = nVertexDomain[nMarker];
                            nMarker++;
                        }

                        if (nReceivedDomain_Periodic[jDomain] != 0)
                        {
                            nVertexDomain[nMarker] = 0;
                            bound[nMarker] = new GRID::GRID_Primal*[nReceivedDomain_Periodic[jDomain]];

                            iVertex = 0;
                            for (iTotalReceivedDomain_Periodic = 0; iTotalReceivedDomain_Periodic < nTotalReceivedDomain_Periodic; iTotalReceivedDomain_Periodic++)
                            {
                                if (Buffer_Receive_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] == jDomain)
                                {
                                    bound[nMarker][iVertex] = new GRID::GRID_VertexMPI(Global_to_local_Point_recv[Buffer_Receive_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic]], nDim);
                                    bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic]);
                                    nVertexDomain[nMarker]++; iVertex++;
                                }
                            }

                            Marker_All_SendRecv[nMarker] = -1 * (jDomain + 1);
                            nElem_Bound[nMarker] = nVertexDomain[nMarker];
                            nMarker++;
                        }
                    }

                    delete[] TagBound_Copy;
                    delete[] SendRecv_Copy;

                    delete[] Buffer_Receive_BoundLine;
                    delete[] Buffer_Receive_BoundTriangle;
                    delete[] Buffer_Receive_BoundRectangle;
                    delete[] Buffer_Receive_Local2Global_Marker;

                    delete[] Buffer_Receive_SendDomain_Periodic;
                    delete[] Buffer_Receive_SendDomain_PeriodicTrans;
                    delete[] Buffer_Receive_SendDomain_PeriodicReceptor;
                    delete[] Buffer_Receive_ReceivedDomain_Periodic;
                    delete[] Buffer_Receive_ReceivedDomain_PeriodicTrans;
                    delete[] Buffer_Receive_ReceivedDomain_PeriodicDonor;
                }
            }

            /*--- The MASTER should wait for the sends above to complete ---*/

#ifdef HAVE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            /*--- Set the value of Marker_All_SendRecv and Marker_All_TagBound in the config structure ---*/

            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
            }

            /*--- Set the value of Global_nPoint and Global_nPointDomain ---*/

            unsigned long Local_nPoint = nPoint;
            unsigned long Local_nPointDomain = nPointDomain;

#ifdef HAVE_MPI
            MPI_Allreduce(&Local_nPoint, &Global_nPoint, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
            Global_nPoint = Local_nPoint;
            Global_nPointDomain = Local_nPointDomain;
#endif

            if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
                std::cout << Global_nPoint << " vertices including ghost points. " << std::endl;

            /*--- Release all of the temporary memory ---*/

            delete[] nDim_s;
            delete[] nDim_r;

            delete[] nPointTotal_s;
            delete[] nPointDomainTotal_s;
            delete[] nPointGhost_s;
            delete[] nPointPeriodic_s;
            delete[] nElemTotal_s;
            delete[] nElemTriangle_s;
            delete[] nElemRectangle_s;
            delete[] nElemTetrahedron_s;
            delete[] nElemHexahedron_s;
            delete[] nElemPrism_s;
            delete[] nElemPyramid_s;
            delete[] nZone_s;

            delete[] nPointTotal_r;
            delete[] nPointDomainTotal_r;
            delete[] nPointGhost_r;
            delete[] nPointPeriodic_r;
            delete[] nElemTotal_r;
            delete[] nElemTriangle_r;
            delete[] nElemRectangle_r;
            delete[] nElemTetrahedron_r;
            delete[] nElemHexahedron_r;
            delete[] nElemPrism_r;
            delete[] nElemPyramid_r;
            delete[] nZone_r;

            if (rank == TBOX::MASTER_NODE)
            {
                delete[] MarkerIn;
                delete[] Buffer_Send_Center;
                delete[] Buffer_Send_Rotation;
                delete[] Buffer_Send_Translate;
                delete[] Buffer_Send_nSendDomain_Periodic;
                delete[] Buffer_Send_nReceivedDomain_Periodic;
                delete[] Marker_All_SendRecv_Copy;
                delete[] Marker_All_TagBound_Copy;
                delete[] PointIn;
                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    delete VertexIn[iMarker];
                delete[] VertexIn;
            }

            delete[] Marker_All_TagBound;
            delete[] Buffer_Send_Marker_All_TagBound;

            delete[] nSendDomain_Periodic;
            delete[] nReceivedDomain_Periodic;
            delete[] nVertexDomain;
            delete[] nBoundLine;
            delete[] nBoundTriangle;
            delete[] nBoundRectangle;
            delete[] Buffer_Send_nVertexDomain;
            delete[] Buffer_Send_nBoundLine;
            delete[] Buffer_Send_nBoundTriangle;
            delete[] Buffer_Send_nBoundRectangle;
            delete[] Buffer_Send_Marker_All_SendRecv;

#ifdef HAVE_MPI
            delete[] send_stat;
            delete[] recv_stat;
            delete[] send_req;
            delete[] recv_req;
#endif
        }

        GEOM_GeometryPhysical::~GEOM_GeometryPhysical(void)
        {
            if (Global_to_Local_Point != NULL)
                delete[] Global_to_Local_Point;
            if (Local_to_Global_Point != NULL)
                delete[] Local_to_Global_Point;
            if (Global_to_Local_Marker != NULL)
                delete[] Global_to_Local_Marker;
            if (Local_to_Global_Marker != NULL)
                delete[] Local_to_Global_Marker;
        }

        void GEOM_GeometryPhysical::SetSendReceive(TBOX::TBOX_Config *config)
        {
            unsigned short Counter_Send, Counter_Receive, iMarkerSend, iMarkerReceive;
            unsigned long iVertex, LocalNode;
            unsigned short nMarker_Max = config->GetnMarker_Max();
            unsigned long  iPoint, jPoint, iElem;
            unsigned short nDomain, iNode, iDomain, jDomain, jNode;
            std::vector<unsigned long>::iterator it;

            std::vector<std::vector<unsigned long> > SendTransfLocal;	/*!< \brief std::vector to store the type of transformation for this send point. */
            std::vector<std::vector<unsigned long> > ReceivedTransfLocal;	/*!< \brief std::vector to store the type of transformation for this received point. */
            std::vector<std::vector<unsigned long> > SendDomainLocal; /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */
            std::vector<std::vector<unsigned long> > ReceivedDomainLocal; /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */

            int rank = TBOX::MASTER_NODE;
            int size = TBOX::SINGLE_NODE;

            unsigned long *nVertexDomain = new unsigned long[nMarker_Max];
#ifdef HAVE_MPI
            /*--- MPI initialization ---*/
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            if (rank == TBOX::MASTER_NODE && size > TBOX::SINGLE_NODE)
                std::cout << "Establishing MPI communication patterns." << std::endl;

            nDomain = size;

            SendTransfLocal.resize(nDomain);
            ReceivedTransfLocal.resize(nDomain);
            SendDomainLocal.resize(nDomain);
            ReceivedDomainLocal.resize(nDomain);

            /*--- Loop over the all the points of the element
            to find the points with different colours, and create the send/received list ---*/
            for (iElem = 0; iElem < nElem; iElem++)
            {
                for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++)
                {
                    iPoint = elem[iElem]->GetNode(iNode);
                    iDomain = node[iPoint]->GetColor();

                    if (iDomain == rank)
                    {
                        for (jNode = 0; jNode < elem[iElem]->GetnNodes(); jNode++)
                        {
                            jPoint = elem[iElem]->GetNode(jNode);
                            jDomain = node[jPoint]->GetColor();

                            /*--- If different color and connected by an edge, then we add them to the list ---*/
                            if (iDomain != jDomain)
                            {
                                /*--- We send from iDomain to jDomain the value of iPoint, we save the
                                global value becuase we need to sort the lists ---*/
                                SendDomainLocal[jDomain].push_back(Local_to_Global_Point[iPoint]);
                                /*--- We send from jDomain to iDomain the value of jPoint, we save the
                                global value becuase we need to sort the lists ---*/
                                ReceivedDomainLocal[jDomain].push_back(Local_to_Global_Point[jPoint]);
                            }
                        }
                    }
                }
            }

            /*--- Sort the points that must be sended and delete repeated points, note
            that the sorting should be done with the global point (not the local) ---*/
            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                sort(SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
                it = unique(SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
                SendDomainLocal[iDomain].resize(it - SendDomainLocal[iDomain].begin());
            }

            /*--- Sort the points that must be received and delete repeated points, note
            that the sorting should be done with the global point (not the local) ---*/
            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                sort(ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
                it = unique(ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
                ReceivedDomainLocal[iDomain].resize(it - ReceivedDomainLocal[iDomain].begin());
            }

            /*--- Create Global to Local Point array, note that the array is smaller (Max_GlobalPoint) than the total
            number of points in the simulation  ---*/
            Max_GlobalPoint = 0;
            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                if (Local_to_Global_Point[iPoint] > Max_GlobalPoint)
                    Max_GlobalPoint = Local_to_Global_Point[iPoint];
            }
            Global_to_Local_Point = new long[Max_GlobalPoint + 1]; // +1 to include the bigger point.

            /*--- Initialization of the array with -1 this is important for the FFD ---*/
            for (iPoint = 0; iPoint < Max_GlobalPoint + 1; iPoint++)
                Global_to_Local_Point[iPoint] = -1;

            /*--- Set the value of some of the points ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                Global_to_Local_Point[Local_to_Global_Point[iPoint]] = iPoint;

            /*--- Add the new MPI send receive boundaries, reset the transformation, and save the local value ---*/
            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                if (SendDomainLocal[iDomain].size() != 0)
                {
                    nVertexDomain[nMarker] = SendDomainLocal[iDomain].size();
                    for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++)
                    {
                        SendDomainLocal[iDomain][iVertex] = Global_to_Local_Point[SendDomainLocal[iDomain][iVertex]];
                        SendTransfLocal[iDomain].push_back(0);
                    }
                    nElem_Bound[nMarker] = nVertexDomain[nMarker];
                    bound[nMarker] = new GRID::GRID_Primal*[nElem_Bound[nMarker]];
                    nMarker++;
                }
            }

            /*--- Add the new MPI receive boundaries, reset the transformation, and save the local value ---*/
            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                if (ReceivedDomainLocal[iDomain].size() != 0)
                {
                    nVertexDomain[nMarker] = ReceivedDomainLocal[iDomain].size();
                    for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++)
                    {
                        ReceivedDomainLocal[iDomain][iVertex] = Global_to_Local_Point[ReceivedDomainLocal[iDomain][iVertex]];
                        ReceivedTransfLocal[iDomain].push_back(0);
                    }
                    nElem_Bound[nMarker] = nVertexDomain[nMarker];
                    bound[nMarker] = new GRID::GRID_Primal*[nElem_Bound[nMarker]];
                    nMarker++;
                }
            }

            /*--- First compute the Send/Receive boundaries ---*/
            Counter_Send = 0; 	Counter_Receive = 0;
            for (iDomain = 0; iDomain < nDomain; iDomain++)
                if (SendDomainLocal[iDomain].size() != 0) Counter_Send++;

            for (iDomain = 0; iDomain < nDomain; iDomain++)
                if (ReceivedDomainLocal[iDomain].size() != 0) Counter_Receive++;

            iMarkerSend = nMarker - Counter_Send - Counter_Receive;
            iMarkerReceive = nMarker - Counter_Receive;

            /*--- First we do the send ---*/
            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                if (SendDomainLocal[iDomain].size() != 0)
                {
                    for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerSend); iVertex++)
                    {
                        LocalNode = SendDomainLocal[iDomain][iVertex];
                        bound[iMarkerSend][iVertex] = new GRID::GRID_VertexMPI(LocalNode, nDim);
                        bound[iMarkerSend][iVertex]->SetRotation_Type(SendTransfLocal[iDomain][iVertex]);
                    }
                    Marker_All_SendRecv[iMarkerSend] = iDomain + 1;
                    iMarkerSend++;
                }
            }

            /*--- Second we do the receive ---*/
            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                if (ReceivedDomainLocal[iDomain].size() != 0)
                {
                    for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerReceive); iVertex++)
                    {
                        LocalNode = ReceivedDomainLocal[iDomain][iVertex];
                        bound[iMarkerReceive][iVertex] = new GRID::GRID_VertexMPI(LocalNode, nDim);
                        bound[iMarkerReceive][iVertex]->SetRotation_Type(ReceivedTransfLocal[iDomain][iVertex]);
                    }
                    Marker_All_SendRecv[iMarkerReceive] = -(iDomain + 1);
                    iMarkerReceive++;
                }
            }
            delete[] nVertexDomain;
        }

        void GEOM_GeometryPhysical::SetBoundaries(TBOX::TBOX_Config *config)
        {
            unsigned long iElem_Bound, TotalElem, *nElem_Bound_Copy, iVertex_;
            std::string Grid_Marker;
            unsigned short iDomain, nDomain, iMarkersDomain, iLoop, *DomainCount, nMarker_Physical, Duplicate_SendReceive, *DomainSendCount, **DomainSendMarkers, *DomainReceiveCount, **DomainReceiveMarkers, nMarker_SendRecv, iMarker, iMarker_;
            GRID::GRID_Primal*** bound_Copy;
            short *Marker_All_SendRecv_Copy;
            bool CheckStart;

            int size = TBOX::SINGLE_NODE;

#ifdef HAVE_MPI
            /*--- MPI initialization ---*/
            MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

            nDomain = size + 1;
            /*--- Count the number of physical markers
            in the boundaries ---*/
            nMarker_Physical = 0;
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                if (bound[iMarker][0]->GetVTK_Type() != TBOX::VERTEX)
                {
                    nMarker_Physical++;
                }
            }

            /*--- Identify if there are markers that send/received with the same domain,
            they should be together---*/

            Duplicate_SendReceive = 0;
            for (iLoop = 0; iLoop < 2; iLoop++)
            {
                DomainCount = new unsigned short[nDomain];

                for (iDomain = 0; iDomain < nDomain; iDomain++)
                    DomainCount[iDomain] = 0;

                if (iLoop == 0)
                {
                    for (iDomain = 0; iDomain < nDomain; iDomain++)
                        for (iMarker = 0; iMarker < nMarker; iMarker++)
                        {
                            nMarker_Physical++;
                            if (bound[iMarker][0]->GetVTK_Type() == TBOX::VERTEX)
                                if (Marker_All_SendRecv[iMarker] == iDomain) DomainCount[iDomain]++;
                        }
                }
                else
                {
                    for (iDomain = 0; iDomain < nDomain; iDomain++)
                        for (iMarker = 0; iMarker < nMarker; iMarker++)
                            if (bound[iMarker][0]->GetVTK_Type() == TBOX::VERTEX)
                                if (Marker_All_SendRecv[iMarker] == -iDomain) DomainCount[iDomain]++;
                }

                for (iDomain = 0; iDomain < nDomain; iDomain++)
                    if (DomainCount[iDomain] > 1) Duplicate_SendReceive++;

                delete[] DomainCount;

            }

            DomainSendCount = new unsigned short[nDomain];
            DomainSendMarkers = new unsigned short *[nDomain];
            DomainReceiveCount = new unsigned short[nDomain];
            DomainReceiveMarkers = new unsigned short *[nDomain];

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                DomainSendCount[iDomain] = 0;
                DomainSendMarkers[iDomain] = new unsigned short[nMarker];

                DomainReceiveCount[iDomain] = 0;
                DomainReceiveMarkers[iDomain] = new unsigned short[nMarker];
            }

            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                {
                    if (bound[iMarker][0]->GetVTK_Type() == TBOX::VERTEX)
                    {
                        if (Marker_All_SendRecv[iMarker] == iDomain)
                        {
                            DomainSendMarkers[iDomain][DomainSendCount[iDomain]] = iMarker;
                            DomainSendCount[iDomain]++;
                        }
                        if (Marker_All_SendRecv[iMarker] == -iDomain)
                        {
                            DomainReceiveMarkers[iDomain][DomainReceiveCount[iDomain]] = iMarker;
                            DomainReceiveCount[iDomain]++;
                        }
                    }
                }
            }

            /*--- Create an structure to store the Send/Receive
            boundaries, because they require some reorganization ---*/

            nMarker_SendRecv = nMarker - nMarker_Physical - Duplicate_SendReceive;
            bound_Copy = new GRID::GRID_Primal**[nMarker_Physical + nMarker_SendRecv];
            nElem_Bound_Copy = new unsigned long[nMarker_Physical + nMarker_SendRecv];
            Marker_All_SendRecv_Copy = new short[nMarker_Physical + nMarker_SendRecv];
            iMarker_ = nMarker_Physical;
            iVertex_ = 0;
            CheckStart = false;

            /*--- Copy and allocate the physical markers in the data structure ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                if (bound[iMarker][0]->GetVTK_Type() != TBOX::VERTEX)
                {
                    nElem_Bound_Copy[iMarker] = nElem_Bound[iMarker];
                    bound_Copy[iMarker] = new GRID::GRID_Primal*[nElem_Bound[iMarker]];

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TBOX::LINE)
                            bound_Copy[iMarker][iElem_Bound] = new GRID::GRID_Line(bound[iMarker][iElem_Bound]->GetNode(0),
                            bound[iMarker][iElem_Bound]->GetNode(1), 2);
                        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TBOX::TRIANGLE)

                            bound_Copy[iMarker][iElem_Bound] = new GRID::GRID_Triangle(bound[iMarker][iElem_Bound]->GetNode(0),
                            bound[iMarker][iElem_Bound]->GetNode(1),
                            bound[iMarker][iElem_Bound]->GetNode(2), 3);
                        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TBOX::RECTANGLE)
                            bound_Copy[iMarker][iElem_Bound] = new GRID::GRID_Rectangle(bound[iMarker][iElem_Bound]->GetNode(0),
                            bound[iMarker][iElem_Bound]->GetNode(1),
                            bound[iMarker][iElem_Bound]->GetNode(2),
                            bound[iMarker][iElem_Bound]->GetNode(3), 3);
                    }
                }
            }


            for (iDomain = 0; iDomain < nDomain; iDomain++)
            {

                /*--- Compute the total number of elements (adding all the
                boundaries with the same Send/Receive ---*/

                if (DomainSendCount[iDomain] != 0)
                {
                    TotalElem = 0;
                    for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++)
                    {
                        iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
                        TotalElem += nElem_Bound[iMarker];
                    }
                    if (CheckStart) iMarker_++;
                    CheckStart = true;
                    iVertex_ = 0;
                    nElem_Bound_Copy[iMarker_] = TotalElem;
                    bound_Copy[iMarker_] = new GRID::GRID_Primal*[TotalElem];
                }

                for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++)
                {
                    iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
                    Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        bound_Copy[iMarker_][iVertex_] = new GRID::GRID_VertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
                        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
                        iVertex_++;
                    }

                }

                /*--- Compute the total number of elements (adding all the
                boundaries with the same Send/Receive ---*/

                if (DomainReceiveCount[iDomain] != 0)
                {
                    TotalElem = 0;
                    for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++)
                    {
                        iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
                        TotalElem += nElem_Bound[iMarker];
                    }
                    if (CheckStart) iMarker_++;
                    CheckStart = true;
                    iVertex_ = 0;
                    nElem_Bound_Copy[iMarker_] = TotalElem;
                    bound_Copy[iMarker_] = new GRID::GRID_Primal*[TotalElem];

                }

                for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++)
                {
                    iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
                    Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        bound_Copy[iMarker_][iVertex_] = new GRID::GRID_VertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
                        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
                        iVertex_++;
                    }

                }

            }

            delete[] DomainSendCount;
            for (iDomain = 0; iDomain < nDomain; iDomain++)
                delete DomainSendMarkers[iDomain];
            delete[] DomainSendMarkers;

            delete[] DomainReceiveCount;
            for (iDomain = 0; iDomain < nDomain; iDomain++)
                delete DomainReceiveMarkers[iDomain];
            delete[] DomainReceiveMarkers;

            /*--- Deallocate the bound variables ---*/

            for (iMarker = 0; iMarker < nMarker; iMarker++)
                delete bound[iMarker];
            delete bound;

            /*--- Allocate the new bound variables, and set the number of markers ---*/

            bound = bound_Copy;
            nMarker = nMarker_Physical + nMarker_SendRecv;

            config->SetnMarker_All(nMarker);

            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                nElem_Bound[iMarker] = nElem_Bound_Copy[iMarker];
            }
            for (iMarker = nMarker_Physical; iMarker < nMarker; iMarker++)
            {
                Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
                config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
                config->SetMarker_All_TagBound(iMarker, "SEND_RECEIVE");
            }

            /*--- Update config information storing the boundary information in the right place ---*/

            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {

                std::string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

                if (Marker_Tag != "SEND_RECEIVE")
                {
                    /*--- Update config information storing the boundary information in the right place ---*/
                    Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
                    config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
                    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
                    config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
                    config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
                    config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
                    config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
                    config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
                    config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
                    config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
                }

                /*--- Send-Receive boundaries definition ---*/
                else
                {
                    config->SetMarker_All_KindBC(iMarker, TBOX::SEND_RECEIVE);
                    config->SetMarker_All_Monitoring(iMarker, TBOX::NO);
                    config->SetMarker_All_GeoEval(iMarker, TBOX::NO);
                    config->SetMarker_All_Designing(iMarker, TBOX::NO);
                    config->SetMarker_All_Plotting(iMarker, TBOX::NO);
                    config->SetMarker_All_DV(iMarker, TBOX::NO);
                    config->SetMarker_All_Moving(iMarker, TBOX::NO);
                    config->SetMarker_All_PerBound(iMarker, TBOX::NO);
                    config->SetMarker_All_Out_1D(iMarker, TBOX::NO);

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        if (config->GetMarker_All_SendRecv(iMarker) < 0)
                            node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
                    }
                }

                /*--- Loop over the surface element to set the boundaries ---*/
                unsigned long Point_Surface, iElem_Surface;
                unsigned short iNode_Surface;

                for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
                {
                    for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++)
                    {
                        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
                        node[Point_Surface]->SetBoundary(nMarker);
                        if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE &&
                            config->GetMarker_All_KindBC(iMarker) != TBOX::INTERFACE_BOUNDARY &&
                            config->GetMarker_All_KindBC(iMarker) != TBOX::NEARFIELD_BOUNDARY &&
                            config->GetMarker_All_KindBC(iMarker) != TBOX::PERIODIC_BOUNDARY)
                            node[Point_Surface]->SetPhysicalBoundary(true);

                        if (config->GetMarker_All_KindBC(iMarker) == TBOX::EULER_WALL &&
                            config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX &&
                            config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL)
                            node[Point_Surface]->SetSolidBoundary(true);
                    }
                }

            }

        }

        void GEOM_GeometryPhysical::Check_IntElem_Orientation(TBOX::TBOX_Config *config)
        {
            unsigned long Point_1, Point_2, Point_3, Point_4, Point_5, Point_6, iElem;
            double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
                *Coord_5, *Coord_6, a[3] = { 0.0, 0.0, 0.0 }, b[3] = { 0.0, 0.0, 0.0 }, c[3] = { 0.0, 0.0, 0.0 }, n[3] = { 0.0, 0.0, 0.0 }, test;
            unsigned short iDim;

            /*--- Loop over all the elements ---*/
            for (iElem = 0; iElem < nElem; iElem++)
            {
                /*--- 2D grid, triangle case ---*/
                if (elem[iElem]->GetVTK_Type() == TBOX::TRIANGLE)
                {
                    Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++) {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                    }
                    test = a[0] * b[1] - b[0] * a[1];

                    if (test < 0.0) elem[iElem]->Change_Orientation();
                }

                /*--- 2D grid, rectangle case ---*/

                if (elem[iElem]->GetVTK_Type() == TBOX::RECTANGLE)
                {
                    Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                    }
                    test_1 = a[0] * b[1] - b[0] * a[1];

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_3[iDim] - Coord_2[iDim]);
                        b[iDim] = 0.5*(Coord_4[iDim] - Coord_2[iDim]);
                    }
                    test_2 = a[0] * b[1] - b[0] * a[1];

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_4[iDim] - Coord_3[iDim]);
                        b[iDim] = 0.5*(Coord_1[iDim] - Coord_3[iDim]);
                    }
                    test_3 = a[0] * b[1] - b[0] * a[1];

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_1[iDim] - Coord_4[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_4[iDim]);
                    }
                    test_4 = a[0] * b[1] - b[0] * a[1];

                    if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
                        elem[iElem]->Change_Orientation();
                }

                /*--- 3D grid, tetrahedron case ---*/
                if (elem[iElem]->GetVTK_Type() == TBOX::TETRAHEDRON)
                {

                    Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];
                    if (test < 0.0) elem[iElem]->Change_Orientation();
                }

                /*--- 3D grid, prism case ---*/
                if (elem[iElem]->GetVTK_Type() == TBOX::PRISM)
                {
                    Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
                    Point_5 = elem[iElem]->GetNode(4); Coord_5 = node[Point_5]->GetCoord();
                    Point_6 = elem[iElem]->GetNode(5); Coord_6 = node[Point_6]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        c[iDim] = (Coord_4[iDim] - Coord_1[iDim]) +
                            (Coord_5[iDim] - Coord_2[iDim]) +
                            (Coord_6[iDim] - Coord_3[iDim]);
                    }

                    /*--- The normal std::vector should point to the interior of the element ---*/

                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_5[iDim] - Coord_4[iDim]);
                        b[iDim] = 0.5*(Coord_6[iDim] - Coord_4[iDim]);
                        c[iDim] = (Coord_1[iDim] - Coord_4[iDim]) +
                            (Coord_2[iDim] - Coord_5[iDim]) +
                            (Coord_3[iDim] - Coord_6[iDim]);
                    }

                    /*--- The normal std::vector should point to the interior of the element ---*/

                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    if ((test_1 < 0.0) || (test_2 < 0.0))
                        elem[iElem]->Change_Orientation();

                }

                if (elem[iElem]->GetVTK_Type() == TBOX::HEXAHEDRON)
                {

                    Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(5); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(7); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    Point_1 = elem[iElem]->GetNode(1); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(2); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(3); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(6); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_3 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    Point_1 = elem[iElem]->GetNode(3); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(0); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(1); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_4 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    if ((test_1 < 0.0) || (test_2 < 0.0) || (test_3 < 0.0)
                        || (test_4 < 0.0)) elem[iElem]->Change_Orientation();

                }

                if (elem[iElem]->GetVTK_Type() == TBOX::PYRAMID)
                {

                    Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
                    Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();

                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                        b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                    }
                    n[0] = a[1] * b[2] - b[1] * a[2];
                    n[1] = -(a[0] * b[2] - b[0] * a[2]);
                    n[2] = a[0] * b[1] - b[0] * a[1];

                    test_2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                    if ((test_1 < 0.0) || (test_2 < 0.0))
                        elem[iElem]->Change_Orientation();
                }
            }
        }

        void GEOM_GeometryPhysical::Check_BoundElem_Orientation(TBOX::TBOX_Config *config)
        {
            unsigned long Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_4_Surface,
                iElem_Domain, Point_Domain = 0, Point_Surface, iElem_Surface;
            double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
                *Coord_5, a[3] = { 0.0, 0.0, 0.0 }, b[3] = { 0.0, 0.0, 0.0 }, c[3] = { 0.0, 0.0, 0.0 }, n[3] = { 0.0, 0.0, 0.0 }, test;
            unsigned short iDim, iMarker, iNode_Domain, iNode_Surface;
            bool find;

            /*--- Surface elements ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
                {
                    iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
                    for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++)
                    {
                        Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
                        find = false;
                        for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++)
                        {
                            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
                            if (Point_Surface == Point_Domain)
                            {
                                find = true;
                                break;
                            }
                        }
                        if (!find) break;
                    }

                    /*--- 2D grid, line case ---*/

                    if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TBOX::LINE)
                    {
                        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
                        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
                        Coord_3 = node[Point_Domain]->GetCoord();

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                            b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                        }
                        test = a[0] * b[1] - b[0] * a[1];

                        if (test < 0.0)
                        {
                            bound[iMarker][iElem_Surface]->Change_Orientation();
                            node[Point_1_Surface]->SetFlip_Orientation();
                            node[Point_2_Surface]->SetFlip_Orientation();
                        }

                    }

                    /*--- 3D grid, triangle case ---*/
                    if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TBOX::TRIANGLE)
                    {
                        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
                        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
                        Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
                        Coord_4 = node[Point_Domain]->GetCoord();

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                            b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                            c[iDim] = Coord_4[iDim] - Coord_1[iDim];
                        }
                        n[0] = a[1] * b[2] - b[1] * a[2];
                        n[1] = -(a[0] * b[2] - b[0] * a[2]);
                        n[2] = a[0] * b[1] - b[0] * a[1];

                        test = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];
                        if (test < 0.0)
                        {
                            bound[iMarker][iElem_Surface]->Change_Orientation();
                            node[Point_1_Surface]->SetFlip_Orientation();
                            node[Point_2_Surface]->SetFlip_Orientation();
                            node[Point_3_Surface]->SetFlip_Orientation();
                        }
                    }

                    if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TBOX::RECTANGLE)
                    {
                        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
                        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
                        Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
                        Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3); Coord_4 = node[Point_4_Surface]->GetCoord();
                        Coord_5 = node[Point_Domain]->GetCoord();

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            a[iDim] = 0.5*(Coord_2[iDim] - Coord_1[iDim]);
                            b[iDim] = 0.5*(Coord_3[iDim] - Coord_1[iDim]);
                            c[iDim] = Coord_5[iDim] - Coord_1[iDim];
                        }
                        n[0] = a[1] * b[2] - b[1] * a[2];
                        n[1] = -(a[0] * b[2] - b[0] * a[2]);
                        n[2] = a[0] * b[1] - b[0] * a[1];
                        test_1 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            a[iDim] = 0.5*(Coord_3[iDim] - Coord_2[iDim]);
                            b[iDim] = 0.5*(Coord_4[iDim] - Coord_2[iDim]);
                            c[iDim] = Coord_5[iDim] - Coord_2[iDim];
                        }
                        n[0] = a[1] * b[2] - b[1] * a[2];
                        n[1] = -(a[0] * b[2] - b[0] * a[2]);
                        n[2] = a[0] * b[1] - b[0] * a[1];
                        test_2 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            a[iDim] = 0.5*(Coord_4[iDim] - Coord_3[iDim]);
                            b[iDim] = 0.5*(Coord_1[iDim] - Coord_3[iDim]);
                            c[iDim] = Coord_5[iDim] - Coord_3[iDim];
                        }
                        n[0] = a[1] * b[2] - b[1] * a[2];
                        n[1] = -(a[0] * b[2] - b[0] * a[2]);
                        n[2] = a[0] * b[1] - b[0] * a[1];
                        test_3 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            a[iDim] = 0.5*(Coord_1[iDim] - Coord_4[iDim]);
                            b[iDim] = 0.5*(Coord_3[iDim] - Coord_4[iDim]);
                            c[iDim] = Coord_5[iDim] - Coord_4[iDim];
                        }
                        n[0] = a[1] * b[2] - b[1] * a[2];
                        n[1] = -(a[0] * b[2] - b[0] * a[2]);
                        n[2] = a[0] * b[1] - b[0] * a[1];
                        test_4 = n[0] * c[0] + n[1] * c[1] + n[2] * c[2];

                        if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
                        {
                            bound[iMarker][iElem_Surface]->Change_Orientation();
                            node[Point_1_Surface]->SetFlip_Orientation();
                            node[Point_2_Surface]->SetFlip_Orientation();
                            node[Point_3_Surface]->SetFlip_Orientation();
                            node[Point_4_Surface]->SetFlip_Orientation();
                        }
                    }
                }
            }
        }

        void GEOM_GeometryPhysical::ComputeWall_Distance(TBOX::TBOX_Config *config) 
        {
            double *coord, dist2, dist;
            unsigned short iDim, iMarker;
            unsigned long iPoint, iVertex, nVertex_SolidWall;

            int rank = TBOX::MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
            if (rank == TBOX::MASTER_NODE)
                std::cout << "Computing wall distances." << std::endl;

#ifndef HAVE_MPI

            /*--- Compute the total number of nodes on no-slip boundaries ---*/
            nVertex_SolidWall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if ((config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_NONCATALYTIC))
                    nVertex_SolidWall += GetnVertex(iMarker);

            /*--- Allocate an array to hold boundary node coordinates ---*/
            double **Coord_bound;
            Coord_bound = new double*[nVertex_SolidWall];
            for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
                Coord_bound[iVertex] = new double[nDim];

            /*--- Retrieve and store the coordinates of the no-slip boundary nodes ---*/
            nVertex_SolidWall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                if ((config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_NONCATALYTIC))
                    for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) 
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Coord_bound[nVertex_SolidWall][iDim] = node[iPoint]->GetCoord(iDim);
                        nVertex_SolidWall++;
                    }
            }

            /*--- Loop over all interior mesh nodes and compute the distances to each
            of the no-slip boundary nodes. Store the minimum distance to the wall for
            each interior mesh node. ---*/

            if (nVertex_SolidWall != 0) 
            {
                for (iPoint = 0; iPoint < GetnPoint(); iPoint++) 
                {
                    coord = node[iPoint]->GetCoord();
                    dist = 1E20;
                    for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++) 
                    {
                        dist2 = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++)
                            dist2 += (coord[iDim] - Coord_bound[iVertex][iDim])
                            *(coord[iDim] - Coord_bound[iVertex][iDim]);
                        if (dist2 < dist) dist = dist2;
                    }
                    node[iPoint]->SetWall_Distance(sqrt(dist));
                }
            }
            else 
            {
                for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
                    node[iPoint]->SetWall_Distance(0.0);
            }

            /*--- Deallocate the std::vector of boundary coordinates. ---*/

            for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
                delete[] Coord_bound[iVertex];
            delete[] Coord_bound;


#else

            /*--- Variables and buffers needed for MPI ---*/

            int iProcessor, nProcessor;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            unsigned long nLocalVertex_NS = 0, nGlobalVertex_NS = 0, MaxLocalVertex_NS = 0;
            unsigned long *Buffer_Send_nVertex = new unsigned long[1];
            unsigned long *Buffer_Receive_nVertex = new unsigned long[nProcessor];

            /*--- Count the total number of nodes on no-slip boundaries within the
            local partition. ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
                    (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC))
                    nLocalVertex_NS += GetnVertex(iMarker);

            /*--- Communicate to all processors the total number of no-slip boundary
            nodes, the maximum number of no-slip boundary nodes on any single single
            partition, and the number of no-slip nodes on each partition. ---*/

            Buffer_Send_nVertex[0] = nLocalVertex_NS;
            MPI_Allreduce(&nLocalVertex_NS, &nGlobalVertex_NS, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&nLocalVertex_NS, &MaxLocalVertex_NS, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

            /*--- Create and initialize to zero some buffers to hold the coordinates
            of the boundary nodes that are communicated from each partition (all-to-all). ---*/

            double *Buffer_Send_Coord = new double[MaxLocalVertex_NS*nDim];
            double *Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex_NS*nDim];
            unsigned long nBuffer = MaxLocalVertex_NS*nDim;

            for (iVertex = 0; iVertex < MaxLocalVertex_NS; iVertex++)
                for (iDim = 0; iDim < nDim; iDim++)
                    Buffer_Send_Coord[iVertex*nDim + iDim] = 0.0;

            /*--- Retrieve and store the coordinates of the no-slip boundary nodes on
            the local partition and broadcast them to all partitions. ---*/

            nVertex_SolidWall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
                    (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC))
                    for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_Coord[nVertex_SolidWall*nDim + iDim] = node[iPoint]->GetCoord(iDim);
                        nVertex_SolidWall++;
                    }

            MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);

            /*--- Loop over all interior mesh nodes on the local partition and compute
            the distances to each of the no-slip boundary nodes in the entire mesh.
            Store the minimum distance to the wall for each interior mesh node. ---*/

            nVertex_SolidWall = 0;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                nVertex_SolidWall += Buffer_Receive_nVertex[iProcessor];
            }

            if (nVertex_SolidWall != 0) {
                for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
                    coord = node[iPoint]->GetCoord();
                    dist = 1E20;
                    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
                            dist2 = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++)
                                dist2 += (coord[iDim] - Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NS + iVertex)*nDim + iDim])*
                                (coord[iDim] - Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NS + iVertex)*nDim + iDim]);
                            if (dist2 < dist) dist = dist2;
                        }
                    node[iPoint]->SetWall_Distance(sqrt(dist));
                }
            }
            else {
                for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
                    node[iPoint]->SetWall_Distance(0.0);
            }

            /*--- Deallocate the buffers needed for the MPI communication. ---*/

            delete[] Buffer_Send_Coord;
            delete[] Buffer_Receive_Coord;
            delete[] Buffer_Send_nVertex;
            delete[] Buffer_Receive_nVertex;

#endif

        }

        void GEOM_GeometryPhysical::SetPositive_ZArea(TBOX::TBOX_Config *config)
        {
            unsigned short iMarker, Boundary, Monitoring;
            unsigned long iVertex, iPoint;
            double *Normal, PositiveZArea;
            int rank = TBOX::MASTER_NODE;

#ifndef HAVE_MPI

            PositiveZArea = 0.0;
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                Boundary = config->GetMarker_All_KindBC(iMarker);
                Monitoring = config->GetMarker_All_Monitoring(iMarker);

                if (((Boundary == TBOX::EULER_WALL) ||
                    (Boundary == TBOX::HEAT_FLUX) ||
                    (Boundary == TBOX::HEAT_FLUX_CATALYTIC) ||
                    (Boundary == TBOX::HEAT_FLUX_NONCATALYTIC) ||
                    (Boundary == TBOX::ISOTHERMAL) ||
                    (Boundary == TBOX::ISOTHERMAL_CATALYTIC) ||
                    (Boundary == TBOX::ISOTHERMAL_NONCATALYTIC) ||
                    (Boundary == TBOX::LOAD_BOUNDARY) ||
                    (Boundary == TBOX::DISPLACEMENT_BOUNDARY)) && (Monitoring == TBOX::YES))
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        if (node[iPoint]->GetDomain())
                        {
                            Normal = vertex[iMarker][iVertex]->GetNormal();
                            if (Normal[nDim - 1] < 0) PositiveZArea -= Normal[nDim - 1];
                        }
                    }
            }

#else

            double TotalPositiveZArea;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            PositiveZArea = 0.0;
            for (iMarker = 0; iMarker < nMarker; iMarker++) {
                Boundary = config->GetMarker_All_KindBC(iMarker);
                Monitoring = config->GetMarker_All_Monitoring(iMarker);

                if (((Boundary == EULER_WALL) ||
                    (Boundary == HEAT_FLUX) ||
                    (Boundary == HEAT_FLUX_CATALYTIC) ||
                    (Boundary == HEAT_FLUX_NONCATALYTIC) ||
                    (Boundary == ISOTHERMAL) ||
                    (Boundary == ISOTHERMAL_CATALYTIC) ||
                    (Boundary == ISOTHERMAL_NONCATALYTIC) ||
                    (Boundary == LOAD_BOUNDARY) ||
                    (Boundary == DISPLACEMENT_BOUNDARY)) && (Monitoring == YES))
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        if (node[iPoint]->GetDomain()) {
                            Normal = vertex[iMarker][iVertex]->GetNormal();
                            if (Normal[nDim - 1] < 0) PositiveZArea -= Normal[nDim - 1];
                        }
                    }
            }
            MPI_Reduce(&PositiveZArea, &TotalPositiveZArea, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == MASTER_NODE) PositiveZArea = TotalPositiveZArea;
            MPI_Bcast(&PositiveZArea, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

#endif

            if (config->GetRefAreaCoeff() == 0.0)
                config->SetRefAreaCoeff(PositiveZArea);

            if (rank == TBOX::MASTER_NODE)
            {
                if (nDim == 2) 
                    std::cout << "Area projection in the y-plane = " << PositiveZArea << "." << std::endl;
                else 
                    std::cout << "Area projection in the z-plane = " << PositiveZArea << "." << std::endl;
            }
        }

        void GEOM_GeometryPhysical::SetPoint_Connectivity(void) 
        {
            unsigned short Node_Neighbor, iNode, iNeighbor;
            unsigned long jElem, Point_Neighbor, iPoint, iElem;

            /*--- Loop over all the elements ---*/
            for (iElem = 0; iElem < nElem; iElem++)
            {
                /*--- Loop over all the nodes of an element ---*/
                for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++)
                {
                    iPoint = elem[iElem]->GetNode(iNode);

                    /*--- Store the element into the point ---*/
                    node[iPoint]->SetElem(iElem);
                }
            }

            /*--- Loop over all the points ---*/

            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                /*--- Loop over all elements shared by the point ---*/
                for (iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++)
                {

                    jElem = node[iPoint]->GetElem(iElem);

                    /*--- If we find the point iPoint in the surronding element ---*/

                    for (iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)
                    {
                        if (elem[jElem]->GetNode(iNode) == iPoint)
                        {
                            /*--- Localize the local index of the neighbor of iPoint in the element ---*/
                            for (iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++)
                            {
                                Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode, iNeighbor);
                                Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);

                                /*--- Store the point into the point ---*/
                                node[iPoint]->SetPoint(Point_Neighbor);
                            }
                        }
                    }
                }
            }

            /*--- Set the number of neighbors variable, this is
            important for JST and multigrid in parallel ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());

        }

        void GEOM_GeometryPhysical::SetRCM_Ordering(TBOX::TBOX_Config *config)
        {
            unsigned long iPoint, AdjPoint, AuxPoint, AddPoint, iElem, iNode, jNode;
            std::vector<unsigned long> Queue, AuxQueue, Result;
            unsigned short Degree, MinDegree, iDim, iMarker;
            bool *inQueue;

            inQueue = new bool[nPoint];

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                inQueue[iPoint] = false;

            /*--- Select the node with the lowest degree in the grid. ---*/

            MinDegree = node[0]->GetnNeighbor(); AddPoint = 0;
            for (iPoint = 1; iPoint < nPointDomain; iPoint++)
            {
                Degree = node[iPoint]->GetnPoint();
                if (Degree < MinDegree)
                {
                    MinDegree = Degree;
                    AddPoint = iPoint;
                }
            }

            /*--- Add the node in the first free position. ---*/
            Result.push_back(AddPoint); inQueue[AddPoint] = true;

            /*--- Loop until reorganize all the nodes ---*/
            do {

                /*--- Add to the queue all the nodes adjacent in the increasing
                order of their degree, checking if the element is already
                in the Queue. ---*/
                AuxQueue.clear();
                for (iNode = 0; iNode < node[AddPoint]->GetnPoint(); iNode++)
                {
                    AdjPoint = node[AddPoint]->GetPoint(iNode);
                    if ((!inQueue[AdjPoint]) && (AdjPoint < nPointDomain))
                    {
                        AuxQueue.push_back(AdjPoint);
                    }
                }

                if (AuxQueue.size() != 0)
                {
                    /*--- Sort the auxiliar queue based on the number of neighbors ---*/
                    for (iNode = 0; iNode < AuxQueue.size(); iNode++)
                    {
                        for (jNode = 0; jNode < AuxQueue.size() - 1 - iNode; jNode++)
                        {
                            if (node[AuxQueue[jNode]]->GetnPoint() > node[AuxQueue[jNode + 1]]->GetnPoint())
                            {
                                AuxPoint = AuxQueue[jNode];
                                AuxQueue[jNode] = AuxQueue[jNode + 1];
                                AuxQueue[jNode + 1] = AuxPoint;
                            }
                        }
                    }

                    Queue.insert(Queue.end(), AuxQueue.begin(), AuxQueue.end());
                    for (iNode = 0; iNode < AuxQueue.size(); iNode++)
                    {
                        inQueue[AuxQueue[iNode]] = true;
                    }

                }

                /*--- Extract the first node from the queue and add it in the first free
                position. ---*/
                if (Queue.size() != 0)
                {
                    AddPoint = Queue[0];
                    Result.push_back(Queue[0]);
                    Queue.erase(Queue.begin(), Queue.begin() + 1);
                }

                /*--- Add to the queue all the nodes adjacent in the increasing
                order of their degree, checking if the element is already
                in the Queue. ---*/
            } while (Queue.size() != 0);

            /*--- Check that all the points have been added ---*/
            for (iPoint = 0; iPoint < nPointDomain; iPoint++)
            {
                if (inQueue[iPoint] == false)
                    Result.push_back(iPoint);
            }

            delete[] inQueue;

            reverse(Result.begin(), Result.end());

            /*--- Add the MPI points ---*/
            for (iPoint = nPointDomain; iPoint < nPoint; iPoint++)
            {
                Result.push_back(iPoint);
            }

            /*--- Reset old data structures ---*/

            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                node[iPoint]->ResetElem();
                node[iPoint]->ResetPoint();
                node[iPoint]->ResetBoundary();
                node[iPoint]->SetPhysicalBoundary(false);
                node[iPoint]->SetSolidBoundary(false);
                node[iPoint]->SetDomain(true);
            }

            /*--- Set the new coordinates ---*/

            double **AuxCoord;
            unsigned long *AuxGlobalIndex;

            AuxGlobalIndex = new unsigned long[nPoint];
            AuxCoord = new double*[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                AuxCoord[iPoint] = new double[nDim];

            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                AuxGlobalIndex[iPoint] = node[iPoint]->GetGlobalIndex();
                for (iDim = 0; iDim < nDim; iDim++)
                {
                    AuxCoord[iPoint][iDim] = node[iPoint]->GetCoord(iDim);
                }
            }

            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                node[iPoint]->SetGlobalIndex(AuxGlobalIndex[Result[iPoint]]);
                for (iDim = 0; iDim < nDim; iDim++)
                    node[iPoint]->SetCoord(iDim, AuxCoord[Result[iPoint]][iDim]);
            }

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                delete[] AuxCoord[iPoint];
            delete[] AuxCoord;
            delete[] AuxGlobalIndex;

            /*--- Set the new conectivities ---*/

            unsigned long *InvResult;
            InvResult = new unsigned long[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                InvResult[Result[iPoint]] = iPoint;

            for (iElem = 0; iElem < nElem; iElem++)
            {
                for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++)
                {
                    iPoint = elem[iElem]->GetNode(iNode);
                    elem[iElem]->SetNode(iNode, InvResult[iPoint]);
                }
            }

            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                {

                    std::string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    if (Marker_Tag == "SEND_RECEIVE")
                    {
                        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                        {
                            if (config->GetMarker_All_SendRecv(iMarker) < 0)
                                node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
                        }
                    }

                    for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = bound[iMarker][iElem]->GetNode(iNode);
                        bound[iMarker][iElem]->SetNode(iNode, InvResult[iPoint]);
                        node[InvResult[iPoint]]->SetBoundary(nMarker);
                        if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE &&
                            config->GetMarker_All_KindBC(iMarker) != TBOX::INTERFACE_BOUNDARY &&
                            config->GetMarker_All_KindBC(iMarker) != TBOX::NEARFIELD_BOUNDARY &&
                            config->GetMarker_All_KindBC(iMarker) != TBOX::PERIODIC_BOUNDARY)
                            node[InvResult[iPoint]]->SetPhysicalBoundary(true);

                        if (config->GetMarker_All_KindBC(iMarker) == TBOX::EULER_WALL ||
                            config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX ||
                            config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL)
                            node[InvResult[iPoint]]->SetSolidBoundary(true);
                    }
                }
            }

            delete[] InvResult;
        }

        void GEOM_GeometryPhysical::SetElement_Connectivity(void)
        {
            unsigned short first_elem_face, second_elem_face, iFace, iNode, jElem;
            unsigned long face_point, Test_Elem, iElem;

            /*--- Loop over all the elements, faces and nodes ---*/

            for (iElem = 0; iElem < nElem; iElem++)
            {
                for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
                {
                    for (iNode = 0; iNode < elem[iElem]->GetnNodesFace(iFace); iNode++)
                    {
                        face_point = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iNode));
                        /*--- Loop over all elements sharing the face point ---*/
                        for (jElem = 0; jElem < node[face_point]->GetnElem(); jElem++)
                        {
                            Test_Elem = node[face_point]->GetElem(jElem);

                            /*--- If it is a new element in this face ---*/
                            if ((elem[iElem]->GetNeighbor_Elements(iFace) == -1) && (iElem < Test_Elem) &&
                                (FindFace(iElem, Test_Elem, first_elem_face, second_elem_face)))
                            {

                                /*--- Localice which faces are sharing both elements ---*/
                                elem[iElem]->SetNeighbor_Elements(Test_Elem, first_elem_face);

                                /*--- Store the element for both elements ---*/
                                elem[Test_Elem]->SetNeighbor_Elements(iElem, second_elem_face);
                            }
                        }
                    }
                }
            }
        }

        void GEOM_GeometryPhysical::SetBoundVolume(void)
        {
            unsigned short cont, iMarker, iElem, iNode_Domain, iNode_Surface;
            unsigned long Point_Domain, Point_Surface, Point, iElem_Surface, iElem_Domain;
            bool CheckVol;

            for (iMarker = 0; iMarker < nMarker; iMarker++)
                for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
                {

                    /*--- Choose and arbitrary point from the surface --*/
                    Point = bound[iMarker][iElem_Surface]->GetNode(0);
                    CheckVol = false;

                    for (iElem = 0; iElem < node[Point]->GetnElem(); iElem++)
                    {
                        /*--- Look for elements surronding that point --*/
                        cont = 0; iElem_Domain = node[Point]->GetElem(iElem);
                        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++)
                        {
                            Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
                            for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++)
                            {
                                Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
                                if (Point_Surface == Point_Domain) cont++;
                                if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
                            }
                            if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
                        }

                        if (cont == bound[iMarker][iElem_Surface]->GetnNodes())
                        {
                            bound[iMarker][iElem_Surface]->SetDomainElement(iElem_Domain);
                            CheckVol = true;
                            break;
                        }
                    }
                    if (!CheckVol)
                    {
                        std::cout << "The surface element (" << iMarker << ", " << iElem_Surface << ") doesn't have an associated volume element." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
        }

        void GEOM_GeometryPhysical::SetVertex(TBOX::TBOX_Config *config)
        {
            unsigned long  iPoint, iVertex, iElem;
            unsigned short iMarker, iNode;

            /*--- Initialize the Vertex std::vector for each node of the grid ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                    node[iPoint]->SetVertex(-1, iMarker);

            /*--- Create and compute the std::vector with the number of vertex per marker ---*/
            nVertex = new unsigned long[nMarker];
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {

                /*--- Initialize the number of Bound Vertex for each Marker ---*/
                nVertex[iMarker] = 0;
                for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                    for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = bound[iMarker][iElem]->GetNode(iNode);

                        /*--- Set the vertex in the node information ---*/

                        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE))
                        {
                            node[iPoint]->SetVertex(nVertex[iMarker], iMarker);
                            nVertex[iMarker]++;
                        }
                    }
            }

            /*--- Initialize the Vertex std::vector for each node, the previous result is deleted ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                    node[iPoint]->SetVertex(-1, iMarker);

            /*--- Create the bound vertex structure, note that the order
            is the same as in the input file, this is important for Send/Receive part ---*/
            vertex = new GRID::GRID_DGVertex**[nMarker];
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                vertex[iMarker] = new GRID::GRID_DGVertex*[nVertex[iMarker]];
                nVertex[iMarker] = 0;

                /*--- Initialize the number of Bound Vertex for each Marker ---*/

                for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                    for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = bound[iMarker][iElem]->GetNode(iNode);

                        /*--- Set the vertex in the node information ---*/

                        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE))
                        {
                            iVertex = nVertex[iMarker];
                            vertex[iMarker][iVertex] = new GRID::GRID_DGVertex(iPoint, nDim);

                            if (config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE)
                            {
                                vertex[iMarker][iVertex]->SetRotation_Type(bound[iMarker][iElem]->GetRotation_Type());
                            }
                            node[iPoint]->SetVertex(nVertex[iMarker], iMarker);
                            nVertex[iMarker]++;
                        }
                    }
            }
        }

        void GEOM_GeometryPhysical::SetCG(void)
        {
            unsigned short nNode, iDim, iMarker, iNode;
            unsigned long elem_poin, edge_poin, iElem, iEdge;
            double **Coord;

            /*--- Compute the center of gravity for elements ---*/
            for (iElem = 0; iElem < nElem; iElem++)
            {
                nNode = elem[iElem]->GetnNodes();
                Coord = new double*[nNode];

                /*--- Store the coordinates for all the element nodes ---*/
                for (iNode = 0; iNode < nNode; iNode++)
                {
                    elem_poin = elem[iElem]->GetNode(iNode);
                    Coord[iNode] = new double[nDim];
                    for (iDim = 0; iDim < nDim; iDim++)
                        Coord[iNode][iDim] = node[elem_poin]->GetCoord(iDim);
                }

                /*--- Compute the element CG coordinates ---*/
                elem[iElem]->SetCG(Coord);

                for (iNode = 0; iNode < nNode; iNode++)
                    if (Coord[iNode] != NULL) delete[] Coord[iNode];
                if (Coord != NULL) delete[] Coord;
            }

            /*--- Center of gravity for face elements ---*/

            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                {
                    nNode = bound[iMarker][iElem]->GetnNodes();
                    Coord = new double*[nNode];

                    /*--- Store the coordinates for all the element nodes ---*/
                    for (iNode = 0; iNode < nNode; iNode++) {
                        elem_poin = bound[iMarker][iElem]->GetNode(iNode);
                        Coord[iNode] = new double[nDim];
                        for (iDim = 0; iDim < nDim; iDim++)
                            Coord[iNode][iDim] = node[elem_poin]->GetCoord(iDim);
                    }
                    /*--- Compute the element CG coordinates ---*/
                    bound[iMarker][iElem]->SetCG(Coord);
                    for (iNode = 0; iNode < nNode; iNode++)
                        if (Coord[iNode] != NULL) delete[] Coord[iNode];
                    if (Coord != NULL) delete[] Coord;
                }
            }

            /*--- Center of gravity for edges ---*/
            for (iEdge = 0; iEdge < nEdge; iEdge++)
            {
                nNode = edge[iEdge]->GetnNodes();
                Coord = new double*[nNode];

                /*--- Store the coordinates for all the element nodes ---*/
                for (iNode = 0; iNode < nNode; iNode++)
                {
                    edge_poin = edge[iEdge]->GetNode(iNode);
                    Coord[iNode] = new double[nDim];
                    for (iDim = 0; iDim < nDim; iDim++)
                        Coord[iNode][iDim] = node[edge_poin]->GetCoord(iDim);
                }

                /*--- Compute the edge CG coordinates ---*/
                edge[iEdge]->SetCG(Coord);
                for (iNode = 0; iNode < nNode; iNode++)
                    if (Coord[iNode] != NULL) delete[] Coord[iNode];
                if (Coord != NULL) delete[] Coord;
            }
        }

        void GEOM_GeometryPhysical::SetBoundControlVolume(TBOX::TBOX_Config *config, unsigned short action)
        {
            unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
            unsigned long Neighbor_Point, iVertex, iPoint, iElem;
            long iEdge;
            double Area, *NormalFace = NULL;

            /*--- Update values of faces of the edge ---*/

            if (action != TBOX::ALLOCATE)
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                        vertex[iMarker][iVertex]->SetZeroValues();

            double *Coord_Edge_CG = new double[nDim];
            double *Coord_Elem_CG = new double[nDim];
            double *Coord_Vertex  = new double[nDim];

            /*--- Loop over all the markers ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
                /*--- Loop over all the boundary elements ---*/
                for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                    /*--- Loop over all the nodes of the boundary ---*/
                    for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = bound[iMarker][iElem]->GetNode(iNode);
                        iVertex = node[iPoint]->GetVertex(iMarker);

                        /*--- Loop over the neighbor nodes, there is a face for each one ---*/
                        for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++)
                        {
                            Neighbor_Node = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                            Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);

                            /*--- Shared edge by the Neighbor Point and the point ---*/

                            iEdge = FindEdge(iPoint, Neighbor_Point);
                            for (iDim = 0; iDim < nDim; iDim++)
                            {
                                Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
                                Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
                                Coord_Vertex[iDim] = node[iPoint]->GetCoord(iDim);
                            }
                            switch (nDim)
                            {
                            case 2:
                                /*--- Store the 2D face ---*/
                                if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex);
                                if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG);
                                break;
                            case 3:
                                /*--- Store the 3D face ---*/
                                if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex);
                                if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex);
                                break;
                            }
                        }
                    }

            delete[] Coord_Edge_CG;
            delete[] Coord_Elem_CG;
            delete[] Coord_Vertex;

            /*--- Check if there is a normal with null area ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                {
                    NormalFace = vertex[iMarker][iVertex]->GetNormal();
                    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim] * NormalFace[iDim];
                    Area = sqrt(Area);
                    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = TBOX::EPS*TBOX::EPS;
                }
        }

        void GEOM_GeometryPhysical::MatchInterface(TBOX::TBOX_Config *config) 
        {
            double epsilon = 1.5e-1;

            unsigned short nMarker_InterfaceBound = config->GetnMarker_InterfaceBound();

            if (nMarker_InterfaceBound != 0) 
            {
#ifndef HAVE_MPI

                unsigned short iMarker, jMarker;
                unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
                double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;

                std::cout << "Set Interface boundary conditions." << std::endl;

                maxdist = 0.0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_KindBC(iMarker) == TBOX::INTERFACE_BOUNDARY) 
                    {
                        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                        {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            Coord_i = node[iPoint]->GetCoord();

                            mindist = 1E6;
                            for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
                                if ((config->GetMarker_All_KindBC(jMarker) == TBOX::INTERFACE_BOUNDARY) && (iMarker != jMarker))
                                    for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++)
                                    {
                                        jPoint = vertex[jMarker][jVertex]->GetNode();
                                        Coord_j = node[jPoint]->GetCoord();
                                        if (nDim == 2) dist = sqrt(pow(Coord_j[0] - Coord_i[0], 2.0) + pow(Coord_j[1] - Coord_i[1], 2.0));
                                        if (nDim == 3) dist = sqrt(pow(Coord_j[0] - Coord_i[0], 2.0) + pow(Coord_j[1] - Coord_i[1], 2.0) + pow(Coord_j[2] - Coord_i[2], 2.0));
                                        if (dist < mindist) 
                                        { 
                                            mindist = dist; 
                                            pPoint = jPoint; 
                                        }
                                    }
                            maxdist = std::max(maxdist, mindist);
                            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, TBOX::MASTER_NODE);

                            if (mindist > epsilon) 
                            {
                                std::cout.precision(10);
                                std::cout << std::endl;
                                std::cout << "   Bad match for point " << iPoint << ".\tNearest";
                                std::cout << " donor distance: " << std::scientific << mindist << ".";
                                vertex[iMarker][iVertex]->SetDonorPoint(iPoint, TBOX::MASTER_NODE);
                                maxdist = std::min(maxdist, 0.0);
                            }

                        }
                        std::cout << "The max distance between points is: " << maxdist << "." << std::endl;
                    }

#else
                unsigned short iMarker, iDim;
                unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
                double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
                int iProcessor, pProcessor = 0;
                unsigned long nLocalVertex_Interface = 0, nGlobalVertex_Interface = 0, MaxLocalVertex_Interface = 0;
                int rank, nProcessor;

                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

                unsigned long *Buffer_Send_nVertex = new unsigned long[1];
                unsigned long *Buffer_Receive_nVertex = new unsigned long[nProcessor];

                if (rank == MASTER_NODE) cout << "Set Interface boundary conditions (if any)." << endl;

                /*--- Compute the number of vertex that have interfase boundary condition
                without including the ghost nodes ---*/

                nLocalVertex_Interface = 0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
                        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            if (node[iPoint]->GetDomain()) nLocalVertex_Interface++;
                        }

                Buffer_Send_nVertex[0] = nLocalVertex_Interface;

                /*--- Send Interface vertex information --*/

                MPI_Allreduce(&nLocalVertex_Interface, &nGlobalVertex_Interface, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&nLocalVertex_Interface, &MaxLocalVertex_Interface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

                double *Buffer_Send_Coord = new double[MaxLocalVertex_Interface*nDim];
                unsigned long *Buffer_Send_Point = new unsigned long[MaxLocalVertex_Interface];

                double *Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex_Interface*nDim];
                unsigned long *Buffer_Receive_Point = new unsigned long[nProcessor*MaxLocalVertex_Interface];

                unsigned long nBuffer_Coord = MaxLocalVertex_Interface*nDim;
                unsigned long nBuffer_Point = MaxLocalVertex_Interface;

                for (iVertex = 0; iVertex < MaxLocalVertex_Interface; iVertex++) {
                    Buffer_Send_Point[iVertex] = 0;
                    for (iDim = 0; iDim < nDim; iDim++)
                        Buffer_Send_Coord[iVertex*nDim + iDim] = 0.0;
                }

                /*--- Copy coordinates and point to the auxiliar std::vector --*/

                nLocalVertex_Interface = 0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
                        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            if (node[iPoint]->GetDomain()) {
                                Buffer_Send_Point[nLocalVertex_Interface] = iPoint;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    Buffer_Send_Coord[nLocalVertex_Interface*nDim + iDim] = node[iPoint]->GetCoord(iDim);
                                nLocalVertex_Interface++;
                            }
                        }

                MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

                /*--- Compute the closest point to a Near-Field boundary point ---*/

                maxdist_local = 0.0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
                        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            if (node[iPoint]->GetDomain()) {

                                /*--- Coordinates of the boundary point ---*/

                                Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;

                                /*--- Loop over all the boundaries to find the pair ---*/
                                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                                    for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                                        jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Interface + jVertex];

                                        /*--- Compute the distance ---*/

                                        dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                                            Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Interface + jVertex)*nDim + iDim];
                                            dist += pow(Coord_j[iDim] - Coord_i[iDim], 2.0);
                                        } dist = sqrt(dist);

                                        if (((dist < mindist) && (iProcessor != rank)) ||
                                            ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                                            mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                                        }
                                    }

                                /*--- Store the value of the pair ---*/

                                maxdist_local = max(maxdist_local, mindist);
                                vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);

                                if (mindist > epsilon) {
                                    cout.precision(10);
                                    cout << endl;
                                    cout << "   Bad match for point " << iPoint << ".\tNearest";
                                    cout << " donor distance: " << scientific << mindist << ".";
                                    vertex[iMarker][iVertex]->SetDonorPoint(iPoint, pProcessor);
                                    maxdist_local = min(maxdist_local, 0.0);
                                }

                            }
                        }
                    }
                }

                MPI_Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);

                if (rank == MASTER_NODE) cout << "The max distance between points is: " << maxdist_global << "." << endl;

                delete[] Buffer_Send_Coord;
                delete[] Buffer_Send_Point;

                delete[] Buffer_Receive_Coord;
                delete[] Buffer_Receive_Point;

                delete[] Buffer_Send_nVertex;
                delete[] Buffer_Receive_nVertex;

#endif
            }
        }

        void GEOM_GeometryPhysical::MatchNearField(TBOX::TBOX_Config *config) 
        {
            double epsilon = 1e-1;

            unsigned short nMarker_NearfieldBound = config->GetnMarker_NearFieldBound();

            if (nMarker_NearfieldBound != 0) 
            {
#ifndef HAVE_MPI
                unsigned short iMarker, jMarker;
                unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
                double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;

                std::cout << "Set Near-Field boundary conditions. " << std::endl;

                maxdist = 0.0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_KindBC(iMarker) == TBOX::NEARFIELD_BOUNDARY) 
                    {
                        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                        {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            Coord_i = node[iPoint]->GetCoord();

                            mindist = 1e10;
                            for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
                                if ((config->GetMarker_All_KindBC(jMarker) == TBOX::NEARFIELD_BOUNDARY) && (iMarker != jMarker))
                                    for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) 
                                    {
                                        jPoint = vertex[jMarker][jVertex]->GetNode();
                                        Coord_j = node[jPoint]->GetCoord();
                                        if (nDim == 2) dist = sqrt(pow(Coord_j[0] - Coord_i[0], 2.0) + pow(Coord_j[1] - Coord_i[1], 2.0));
                                        if (nDim == 3) dist = sqrt(pow(Coord_j[0] - Coord_i[0], 2.0) + pow(Coord_j[1] - Coord_i[1], 2.0) + pow(Coord_j[2] - Coord_i[2], 2.0));
                                        if (dist < mindist) 
                                        { 
                                            mindist = dist; 
                                            pPoint = jPoint; 
                                        }
                                    }
                            maxdist = std::max(maxdist, mindist);
                            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, TBOX::MASTER_NODE);

                            if (mindist > epsilon) 
                            {
                                std::cout.precision(10);
                                std::cout << std::endl;
                                std::cout << "   Bad match for point " << iPoint << ".\tNearest";
                                std::cout << " donor distance: " << std::scientific << mindist << ".";
                                vertex[iMarker][iVertex]->SetDonorPoint(iPoint, TBOX::MASTER_NODE);
                                maxdist = std::min(maxdist, 0.0);
                            }
                        }
                        std::cout << "The max distance between points is: " << maxdist << "." << std::endl;
                    }

#else

                unsigned short iMarker, iDim;
                unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
                double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
                int iProcessor, pProcessor = 0;
                unsigned long nLocalVertex_NearField = 0, nGlobalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
                int rank, nProcessor;

                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

                unsigned long *Buffer_Send_nVertex = new unsigned long[1];
                unsigned long *Buffer_Receive_nVertex = new unsigned long[nProcessor];

                if (rank == MASTER_NODE) cout << "Set Near-Field boundary conditions." << endl;

                /*--- Compute the number of vertex that have nearfield boundary condition
                without including the ghost nodes ---*/

                nLocalVertex_NearField = 0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
                        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            if (node[iPoint]->GetDomain()) nLocalVertex_NearField++;
                        }

                Buffer_Send_nVertex[0] = nLocalVertex_NearField;

                /*--- Send Near-Field vertex information --*/

                MPI_Allreduce(&nLocalVertex_NearField, &nGlobalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

                double *Buffer_Send_Coord = new double[MaxLocalVertex_NearField*nDim];
                unsigned long *Buffer_Send_Point = new unsigned long[MaxLocalVertex_NearField];

                double *Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex_NearField*nDim];
                unsigned long *Buffer_Receive_Point = new unsigned long[nProcessor*MaxLocalVertex_NearField];

                unsigned long nBuffer_Coord = MaxLocalVertex_NearField*nDim;
                unsigned long nBuffer_Point = MaxLocalVertex_NearField;

                for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
                    Buffer_Send_Point[iVertex] = 0;
                    for (iDim = 0; iDim < nDim; iDim++)
                        Buffer_Send_Coord[iVertex*nDim + iDim] = 0.0;
                }

                /*--- Copy coordinates and point to the auxiliar std::vector --*/

                nLocalVertex_NearField = 0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
                        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            if (node[iPoint]->GetDomain()) {
                                Buffer_Send_Point[nLocalVertex_NearField] = iPoint;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    Buffer_Send_Coord[nLocalVertex_NearField*nDim + iDim] = node[iPoint]->GetCoord(iDim);
                                nLocalVertex_NearField++;
                            }
                        }

                MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

                /*--- Compute the closest point to a Near-Field boundary point ---*/

                maxdist_local = 0.0;
                maxdist_global = 0.0;
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {

                        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                            iPoint = vertex[iMarker][iVertex]->GetNode();
                            if (node[iPoint]->GetDomain()) {

                                /*--- Coordinates of the boundary point ---*/

                                Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;

                                /*--- Loop over all the boundaries to find the pair ---*/

                                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                                    for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                                        jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_NearField + jVertex];

                                        /*--- Compute the distance ---*/

                                        dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                                            Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NearField + jVertex)*nDim + iDim];
                                            dist += pow(Coord_j[iDim] - Coord_i[iDim], 2.0);
                                        } dist = sqrt(dist);

                                        if (((dist < mindist) && (iProcessor != rank)) ||
                                            ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                                            mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                                        }
                                    }

                                /*--- Store the value of the pair ---*/

                                maxdist_local = max(maxdist_local, mindist);
                                vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);

                                if (mindist > epsilon) {
                                    cout.precision(10);
                                    cout << endl;
                                    cout << "   Bad match for point " << iPoint << ".\tNearest";
                                    cout << " donor distance: " << scientific << mindist << ".";
                                    vertex[iMarker][iVertex]->SetDonorPoint(iPoint, pProcessor);
                                    maxdist_local = min(maxdist_local, 0.0);
                                }

                            }
                        }

                    }
                }

                MPI_Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);

                if (rank == MASTER_NODE) cout << "The max distance between points is: " << maxdist_global << "." << endl;

                delete[] Buffer_Send_Coord;
                delete[] Buffer_Send_Point;

                delete[] Buffer_Receive_Coord;
                delete[] Buffer_Receive_Point;

                delete[] Buffer_Send_nVertex;
                delete[] Buffer_Receive_nVertex;

#endif
            }
        }

        void GEOM_GeometryPhysical::MatchActuator_Disk(TBOX::TBOX_Config *config) 
        {
            double epsilon = 1e-1;
            unsigned short iMarker, iDim;
            unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
            double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local = 0.0, maxdist_global = 0.0;
            int iProcessor, pProcessor = 0;
            unsigned long nLocalVertex_ActDisk = 0, MaxLocalVertex_ActDisk = 0;
            int rank, nProcessor;
            unsigned short Beneficiary = 0, Donor = 0, iBC;

            unsigned short nMarker_ActDisk_Inlet = config->GetnMarker_ActDisk_Inlet();

            if (nMarker_ActDisk_Inlet != 0) 
            {
                for (iBC = 0; iBC < 2; iBC++) 
                {
                    if (iBC == 0) 
                    { 
                        Beneficiary = TBOX::ACTDISK_INLET; 
                        Donor = TBOX::ACTDISK_OUTLET; 
                    }
                    if (iBC == 1) 
                    { 
                        Beneficiary = TBOX::ACTDISK_OUTLET; 
                        Donor = TBOX::ACTDISK_INLET; 
                }

#ifndef HAVE_MPI
                    rank = TBOX::MASTER_NODE;
                    nProcessor = TBOX::SINGLE_NODE;
#else
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

                    unsigned long *Buffer_Send_nVertex = new unsigned long[1];
                    unsigned long *Buffer_Receive_nVertex = new unsigned long[nProcessor];

                    if ((iBC == 0) && (rank == TBOX::MASTER_NODE)) std::cout << "Set Actuator Disk inlet boundary conditions." << std::endl;
                    if ((iBC == 1) && (rank == TBOX::MASTER_NODE)) std::cout << "Set Actuator Disk outlet boundary conditions." << std::endl;

                    /*--- Compute the number of vertex that have an actuator disk outlet boundary condition
                    without including the ghost nodes ---*/
                    nLocalVertex_ActDisk = 0;
                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == Donor) 
                        {
                            for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) 
                            {
                                iPoint = vertex[iMarker][iVertex]->GetNode();
                                if (node[iPoint]->GetDomain()) nLocalVertex_ActDisk++;
                            }
                        }
                    }

                    Buffer_Send_nVertex[0] = nLocalVertex_ActDisk;

                    /*--- Send actuator disk vertex information --*/
#ifndef HAVE_MPI
                    MaxLocalVertex_ActDisk = nLocalVertex_ActDisk;
                    Buffer_Receive_nVertex[0] = Buffer_Send_nVertex[0];
#else
                    MPI_Allreduce(&nLocalVertex_ActDisk, &MaxLocalVertex_ActDisk, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
                    MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif

                    /*--- Array dimentionalization --*/

                    double *Buffer_Send_Coord = new double[MaxLocalVertex_ActDisk*nDim];
                    unsigned long *Buffer_Send_Point = new unsigned long[MaxLocalVertex_ActDisk];

                    double *Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex_ActDisk*nDim];
                    unsigned long *Buffer_Receive_Point = new unsigned long[nProcessor*MaxLocalVertex_ActDisk];

                    unsigned long nBuffer_Coord = MaxLocalVertex_ActDisk*nDim;
                    unsigned long nBuffer_Point = MaxLocalVertex_ActDisk;

                    for (iVertex = 0; iVertex < MaxLocalVertex_ActDisk; iVertex++) 
                    {
                        Buffer_Send_Point[iVertex] = 0;
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_Coord[iVertex*nDim + iDim] = 0.0;
                    }

                    /*--- Copy coordinates and point to the auxiliar std::vector --*/
                    nLocalVertex_ActDisk = 0;
                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == Donor) 
                        {
                            for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++)
                            {
                                iPoint = vertex[iMarker][iVertex]->GetNode();
                                if (node[iPoint]->GetDomain()) 
                                {
                                    Buffer_Send_Point[nLocalVertex_ActDisk] = iPoint;
                                    for (iDim = 0; iDim < nDim; iDim++)
                                        Buffer_Send_Coord[nLocalVertex_ActDisk*nDim + iDim] = node[iPoint]->GetCoord(iDim);
                                    nLocalVertex_ActDisk++;
                                }
                            }
                        }
                    }

#ifndef HAVE_MPI
                    for (unsigned long iBuffer_Coord = 0; iBuffer_Coord < nBuffer_Coord; iBuffer_Coord++)
                        Buffer_Receive_Coord[iBuffer_Coord] = Buffer_Send_Coord[iBuffer_Coord];
                    for (unsigned long iBuffer_Point = 0; iBuffer_Point < nBuffer_Point; iBuffer_Point++)
                        Buffer_Receive_Point[iBuffer_Point] = Buffer_Send_Point[iBuffer_Point];
#else
                    MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
                    MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif

                    /*--- Compute the closest point to an actuator disk inlet point ---*/

                    maxdist_local = 0.0;

                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
                    {
                        if (config->GetMarker_All_KindBC(iMarker) == Beneficiary)
                        {
                            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                            {
                                iPoint = vertex[iMarker][iVertex]->GetNode();
                                if (node[iPoint]->GetDomain()) 
                                {
                                    /*--- Coordinates of the boundary point ---*/
                                    Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;

                                    /*--- Loop over all the boundaries to find the pair ---*/
                                    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) 
                                    {
                                        for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) 
                                        {
                                            jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_ActDisk + jVertex];

                                            /*--- Compute the distance ---*/
                                            dist = 0.0;
                                            for (iDim = 0; iDim < nDim; iDim++) 
                                            {
                                                Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_ActDisk + jVertex)*nDim + iDim];
                                                dist += pow(Coord_j[iDim] - Coord_i[iDim], 2.0);
                                            }
                                            dist = sqrt(dist);

                                            if (dist < mindist) 
                                            {
                                                mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                                                if (dist == 0.0) break;
                                            }
                                        }
                                    }

                                    /*--- Store the value of the pair ---*/

                                    maxdist_local = std::max(maxdist_local, mindist);
                                    vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);

                                    if (mindist > epsilon) 
                                    {
                                        std::cout.precision(10);
                                        std::cout << std::endl;
                                        std::cout << "   Bad match for point " << iPoint << ".\tNearest";
                                        std::cout << " donor distance: " << std::scientific << mindist << ".";
                                        vertex[iMarker][iVertex]->SetDonorPoint(iPoint, pProcessor);
                                        maxdist_local = std::min(maxdist_local, 0.0);
                                    }
                                }
                            }
                        }
                    }

#ifndef HAVE_MPI
                    maxdist_global = maxdist_local;
#else
                    MPI_Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
#endif

                    if (rank == TBOX::MASTER_NODE) std::cout << "The max distance between points is: " << maxdist_global << "." << std::endl;

                    delete[] Buffer_Send_Coord;
                    delete[] Buffer_Send_Point;

                    delete[] Buffer_Receive_Coord;
                    delete[] Buffer_Receive_Point;

                    delete[] Buffer_Send_nVertex;
                    delete[] Buffer_Receive_nVertex;
                }
            }
        }

        void GEOM_GeometryPhysical::MatchZone(TBOX::TBOX_Config *config, GEOM_Geometry *geometry_donor, TBOX::TBOX_Config *config_donor, unsigned short val_iZone, unsigned short val_nZone) 
        {

#ifndef HAVE_MPI

            unsigned short iMarker, jMarker;
            unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
            double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;

            if (val_iZone == TBOX::ZONE_0)
                std::cout << "Set zone boundary conditions (if any)." << std::endl;

            maxdist = 0.0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    Coord_i = node[iPoint]->GetCoord();

                    mindist = 1E6;
                    for (jMarker = 0; jMarker < config_donor->GetnMarker_All(); jMarker++)
                        for (jVertex = 0; jVertex < geometry_donor->GetnVertex(jMarker); jVertex++) 
                        {
                            jPoint = geometry_donor->vertex[jMarker][jVertex]->GetNode();
                            Coord_j = geometry_donor->node[jPoint]->GetCoord();
                            if (nDim == 2) dist = sqrt(pow(Coord_j[0] - Coord_i[0], 2.0) + pow(Coord_j[1] - Coord_i[1], 2.0));
                            if (nDim == 3) dist = sqrt(pow(Coord_j[0] - Coord_i[0], 2.0) + pow(Coord_j[1] - Coord_i[1], 2.0) + pow(Coord_j[2] - Coord_i[2], 2.0));
                            if (dist < mindist) 
                            { 
                                mindist = dist; 
                                pPoint = jPoint; 
                            }
                        }

                    maxdist = std::max(maxdist, mindist);
                    vertex[iMarker][iVertex]->SetDonorPoint(pPoint, TBOX::MASTER_NODE);
                }
            }

#else

            unsigned short iMarker, iDim;
            unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
            double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist;
            int iProcessor, pProcessor = 0;
            unsigned long nLocalVertex_Zone = 0, nGlobalVertex_Zone = 0, MaxLocalVertex_Zone = 0;
            int rank, nProcessor;

            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            unsigned long *Buffer_Send_nVertex = new unsigned long[1];
            unsigned long *Buffer_Receive_nVertex = new unsigned long[nProcessor];

            if (val_iZone == ZONE_0) cout << "Set zone boundary conditions (if any)." << endl;

            nLocalVertex_Zone = 0;
            for (iMarker = 0; iMarker < config_donor->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < geometry_donor->GetnVertex(iMarker); iVertex++) {
                    iPoint = geometry_donor->vertex[iMarker][iVertex]->GetNode();
                    if (geometry_donor->node[iPoint]->GetDomain()) nLocalVertex_Zone++;
                }

            Buffer_Send_nVertex[0] = nLocalVertex_Zone;

            /*--- Send Interface vertex information --*/

            MPI_Allreduce(&nLocalVertex_Zone, &nGlobalVertex_Zone, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&nLocalVertex_Zone, &MaxLocalVertex_Zone, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

            double *Buffer_Send_Coord = new double[MaxLocalVertex_Zone*nDim];
            unsigned long *Buffer_Send_Point = new unsigned long[MaxLocalVertex_Zone];

            double *Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex_Zone*nDim];
            unsigned long *Buffer_Receive_Point = new unsigned long[nProcessor*MaxLocalVertex_Zone];

            unsigned long nBuffer_Coord = MaxLocalVertex_Zone*nDim;
            unsigned long nBuffer_Point = MaxLocalVertex_Zone;

            for (iVertex = 0; iVertex < MaxLocalVertex_Zone; iVertex++) {
                Buffer_Send_Point[iVertex] = 0;
                for (iDim = 0; iDim < nDim; iDim++)
                    Buffer_Send_Coord[iVertex*nDim + iDim] = 0.0;
            }

            /*--- Copy coordinates and point to the auxiliar std::vector --*/
            nLocalVertex_Zone = 0;
            for (iMarker = 0; iMarker < config_donor->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < geometry_donor->GetnVertex(iMarker); iVertex++) {
                    iPoint = geometry_donor->vertex[iMarker][iVertex]->GetNode();
                    if (geometry_donor->node[iPoint]->GetDomain()) {
                        Buffer_Send_Point[nLocalVertex_Zone] = iPoint;
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_Coord[nLocalVertex_Zone*nDim + iDim] = geometry_donor->node[iPoint]->GetCoord(iDim);
                        nLocalVertex_Zone++;
                    }
                }

            MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

            /*--- Compute the closest point to a Near-Field boundary point ---*/
            maxdist = 0.0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                    iPoint = vertex[iMarker][iVertex]->GetNode();

                    if (node[iPoint]->GetDomain()) {

                        /*--- Coordinates of the boundary point ---*/
                        Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;

                        /*--- Loop over all the boundaries to find the pair ---*/
                        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                            for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                                jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Zone + jVertex];

                                /*--- Compute the distance ---*/
                                dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                                    Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Zone + jVertex)*nDim + iDim];
                                    dist += pow(Coord_j[iDim] - Coord_i[iDim], 2.0);
                                } dist = sqrt(dist);

                                if (((dist < mindist) && (iProcessor != rank)) ||
                                    ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                                    mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                                }
                            }

                        /*--- Store the value of the pair ---*/
                        maxdist = max(maxdist, mindist);
                        vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);


                    }
                }
            }

            delete[] Buffer_Send_Coord;
            delete[] Buffer_Send_Point;

            delete[] Buffer_Receive_Coord;
            delete[] Buffer_Receive_Point;

            delete[] Buffer_Send_nVertex;
            delete[] Buffer_Receive_nVertex;

#endif
        }

        void GEOM_GeometryPhysical::SetControlVolume(TBOX::TBOX_Config *config, unsigned short action)
        {
            unsigned long face_iPoint = 0, face_jPoint = 0, iPoint, iElem;
            long iEdge;
            unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
            double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint, *Coord_FacejPoint, Area,
                Volume, DomainVolume, my_DomainVolume, *NormalFace = NULL;
            bool change_face_orientation;
            int rank;

#ifndef HAVE_MPI
            rank = TBOX::MASTER_NODE;
#else
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Update values of faces of the edge ---*/
            if (action != TBOX::ALLOCATE)
            {
                for (iEdge = 0; iEdge < nEdge; iEdge++)
                    edge[iEdge]->SetZeroValues();
                for (iPoint = 0; iPoint < nPoint; iPoint++)
                    node[iPoint]->SetVolume(0.0);
            }

            Coord_Edge_CG = new double[nDim];
            Coord_FaceElem_CG = new double[nDim];
            Coord_Elem_CG = new double[nDim];
            Coord_FaceiPoint = new double[nDim];
            Coord_FacejPoint = new double[nDim];

            my_DomainVolume = 0.0;
            for (iElem = 0; iElem < nElem; iElem++)
                for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
                {

                    /*--- In 2D all the faces have only one edge ---*/
                    if (nDim == 2) nEdgesFace = 1;
                    /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
                    if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);

                    /*-- Loop over the edges of a face ---*/
                    for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++)
                    {

                        /*--- In 2D only one edge (two points) per edge ---*/
                        if (nDim == 2)
                        {
                            face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, 0));
                            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, 1));
                        }

                        /*--- In 3D there are several edges in each face ---*/
                        if (nDim == 3)
                        {
                            face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
                            if (iEdgesFace != nEdgesFace - 1)
                                face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace + 1));
                            else
                                face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, 0));
                        }

                        /*--- We define a direction (from the smalest index to the greatest) --*/
                        change_face_orientation = false;
                        if (face_iPoint > face_jPoint) change_face_orientation = true;
                        iEdge = FindEdge(face_iPoint, face_jPoint);

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
                            Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
                            Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
                            Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
                            Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
                        }

                        switch (nDim)
                        {
                        case 2:
                            /*--- Two dimensional problem ---*/
                            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG);
                            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG);
                            Area = edge[iEdge]->GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_Elem_CG);
                            node[face_iPoint]->AddVolume(Area); my_DomainVolume += Area;
                            Area = edge[iEdge]->GetVolume(Coord_FacejPoint, Coord_Edge_CG, Coord_Elem_CG);
                            node[face_jPoint]->AddVolume(Area); my_DomainVolume += Area;
                            break;
                        case 3:
                            /*--- Three dimensional problem ---*/
                            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_FaceElem_CG, Coord_Edge_CG, Coord_Elem_CG);
                            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
                            Volume = edge[iEdge]->GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
                            node[face_iPoint]->AddVolume(Volume); my_DomainVolume += Volume;
                            Volume = edge[iEdge]->GetVolume(Coord_FacejPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
                            node[face_jPoint]->AddVolume(Volume); my_DomainVolume += Volume;
                            break;
                        }
                    }
                }

            /*--- Check if there is a normal with null area ---*/
            for (iEdge = 0; iEdge < nEdge; iEdge++)
            {
                NormalFace = edge[iEdge]->GetNormal();
                Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim] * NormalFace[iDim];
                Area = sqrt(Area);
                if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = TBOX::EPS*TBOX::EPS;
            }


#ifdef HAVE_MPI
            MPI_Allreduce(&my_DomainVolume, &DomainVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
            DomainVolume = my_DomainVolume;
#endif

            if ((rank == TBOX::MASTER_NODE) && (action == TBOX::ALLOCATE))
            {
                if (nDim == 2) std::cout << "Area of the computational grid: " << DomainVolume << "." << std::endl;
                if (nDim == 3) std::cout << "Volume of the computational grid: " << DomainVolume << "." << std::endl;
            }

            config->SetDomainVolume(DomainVolume);

            delete[] Coord_Edge_CG;
            delete[] Coord_FaceElem_CG;
            delete[] Coord_Elem_CG;
            delete[] Coord_FaceiPoint;
            delete[] Coord_FacejPoint;
        }

        void GEOM_GeometryPhysical::VisualizeControlVolume(TBOX::TBOX_Config *config, unsigned short action)
        {
            /*--- This routine is only meant for visualization in serial currently ---*/
#ifndef HAVE_MPI

            unsigned long face_iPoint = 0, face_jPoint = 0, iElem, iPoint_Viz;
            long iEdge;
            unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
            double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint,
                *Coord_FacejPoint;
            int counter = 0;
            char cstr[TBOX::MAX_STRING_SIZE], buffer[50];
            std::ofstream Tecplot_File;
            std::string mesh_filename;
            std::vector<double> X, Y, Z, X_n, Y_n, Z_n;
            double r1[3], r2[3], CrossProduct[3];

            /*--- Access the point number for control volume we want to vizualize ---*/

            iPoint_Viz = config->GetVisualize_CV();

            /*--- Allocate some structures for building the dual CVs ---*/

            Coord_Edge_CG = new double[nDim];
            Coord_FaceElem_CG = new double[nDim];
            Coord_Elem_CG = new double[nDim];
            Coord_FaceiPoint = new double[nDim];
            Coord_FacejPoint = new double[nDim];

            /*--- Loop over each face of each element ---*/

            CrossProduct[0] = 0.0; CrossProduct[1] = 0.0; CrossProduct[2] = 0.0;

            for (iElem = 0; iElem < nElem; iElem++)
            {
                for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
                {
                    /*--- In 2D all the faces have only one edge ---*/
                    if (nDim == 2) nEdgesFace = 1;
                    /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
                    if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);

                    /*-- Loop over the edges of a face ---*/
                    for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++)
                    {

                        /*--- In 2D only one edge (two points) per edge ---*/
                        if (nDim == 2)
                        {
                            face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, 0));
                            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, 1));
                        }

                        /*--- In 3D there are several edges in each face ---*/
                        if (nDim == 3)
                        {
                            face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
                            if (iEdgesFace != nEdgesFace - 1)
                                face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace + 1));
                            else
                                face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, 0));
                        }

                        /*--- We define a direction (from the smallest index to the greatest) --*/
                        iEdge = FindEdge(face_iPoint, face_jPoint);

                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
                            Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
                            Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
                            Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
                            Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
                        }

                        /*--- Print out the coordinates for a set of triangles making
                        up a single dual control volume for visualization. ---*/

                        if (face_iPoint == iPoint_Viz || face_jPoint == iPoint_Viz)
                        {

                            if (nDim == 2)
                            {
                                X.push_back(Coord_Elem_CG[0]); X.push_back(Coord_Edge_CG[0]);
                                Y.push_back(Coord_Elem_CG[1]); Y.push_back(Coord_Edge_CG[1]);
                            }
                            else if (nDim == 3)
                            {
                                X.push_back(Coord_FaceElem_CG[0]); X.push_back(Coord_Edge_CG[0]); X.push_back(Coord_Elem_CG[0]);
                                Y.push_back(Coord_FaceElem_CG[1]); Y.push_back(Coord_Edge_CG[1]); Y.push_back(Coord_Elem_CG[1]);
                                Z.push_back(Coord_FaceElem_CG[2]); Z.push_back(Coord_Edge_CG[2]); Z.push_back(Coord_Elem_CG[2]);

                                for (iDim = 0; iDim < nDim; iDim++)
                                {
                                    r1[iDim] = Coord_FaceElem_CG[iDim] - Coord_Elem_CG[iDim];
                                    r2[iDim] = Coord_Edge_CG[iDim] - Coord_Elem_CG[iDim];
                                }
                                CrossProduct[0] += 0.5*(r1[1] * r2[2] - r1[2] * r2[1]);
                                CrossProduct[1] += 0.5*(r1[2] * r2[0] - r1[0] * r2[2]);
                                CrossProduct[2] += 0.5*(r1[0] * r2[1] - r1[1] * r2[0]);
                            }
                            counter++;
                        }
                    }
                }
            }

            /*--- Write a Tecplot file to visualize the CV ---*/

            strcpy(cstr, "dual_cv");
            sprintf(buffer, "_%d.dat", int(iPoint_Viz));
            strcat(cstr, buffer);

            Tecplot_File.open(cstr, std::ios::out);
            Tecplot_File << "TITLE= \"Visualization of the control volume\"" << std::endl;

            if (nDim == 2)
            {
                Tecplot_File << "VARIABLES = \"x\",\"y\" " << std::endl;
                Tecplot_File << "ZONE NODES= " << counter * 2 << ", ELEMENTS= ";
                Tecplot_File << counter << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << std::endl;
            }
            if (nDim == 3)
            {
                Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << std::endl;
                Tecplot_File << "ZONE NODES= " << counter * 3 << ", ELEMENTS= ";
                Tecplot_File << counter << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
            }

            /*--- Write coordinates for the nodes in the order that they were found
            for each of the edges/triangles making up a dual control volume. ---*/
            for (std::vector<double>::size_type i = 0; i != X.size(); i++)
            {
                Tecplot_File << X[i] << "\t" << Y[i];
                if (nDim == 3) Tecplot_File << "\t" << Z[i];
                Tecplot_File << "\n";
            }

            /*--- Create a new connectivity table in the order the faces were found ---*/

            int j;
            for (int i = 0; i < counter; i++)
            {
                if (nDim == 2)
                {
                    j = i * 2;
                    Tecplot_File << j + 1 << "\t" << j + 2 << "\t" << j + 2 << "\t" << j + 2 << std::endl;
                }
                if (nDim == 3)
                {
                    j = i * 3;
                    Tecplot_File << j + 1 << "\t" << j + 2 << "\t" << j + 3 << "\t" << j + 3 << "\t";
                    Tecplot_File << j + 3 << "\t" << j + 3 << "\t" << j + 3 << "\t" << j + 3 << std::endl;
                }
            }

            Tecplot_File.close();
            X.clear();
            Y.clear();
            Z.clear();

            delete[] Coord_Edge_CG;
            delete[] Coord_FaceElem_CG;
            delete[] Coord_Elem_CG;
            delete[] Coord_FaceiPoint;
            delete[] Coord_FacejPoint;
#endif
        }

        void GEOM_GeometryPhysical::SetMeshFile(TBOX::TBOX_Config *config, std::string val_mesh_out_filename) 
        {
            unsigned long iElem, iPoint, iElem_Bound;
            unsigned short iMarker, iNodes, iDim;
            unsigned short iPeriodic, nPeriodic = 0;
            std::ofstream output_file;
            std::string Grid_Marker;
            char *cstr;
            double *center, *angles, *transl;

            cstr = new char[val_mesh_out_filename.size() + 1];
            strcpy(cstr, val_mesh_out_filename.c_str());

            /*--- Open .su2 grid file ---*/

            output_file.precision(15);
            output_file.open(cstr, std::ios::out);

            /*--- Write dimension, number of elements and number of points ---*/

            output_file << "NDIME= " << nDim << std::endl;
            output_file << "NELEM= " << nElem << std::endl;
            for (iElem = 0; iElem < nElem; iElem++) 
            {
                output_file << elem[iElem]->GetVTK_Type();
                for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
                    output_file << "\t" << elem[iElem]->GetNode(iNodes);
                output_file << "\t" << iElem << std::endl;
            }

            /*--- Write the node coordinates ---*/

            output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << std::endl;
            output_file.precision(15);
            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                for (iDim = 0; iDim < nDim; iDim++)
                    output_file << std::scientific << "\t" << node[iPoint]->GetCoord(iDim);
#ifndef HAVE_MPI
                output_file << "\t" << iPoint << std::endl;
#else
                output_file << "\t" << iPoint << "\t" << node[iPoint]->GetGlobalIndex() << endl;
#endif

            }

            /*--- Loop through and write the boundary info ---*/
            output_file << "NMARK= " << nMarker << std::endl;
            for (iMarker = 0; iMarker < nMarker; iMarker++) 
            {
                /*--- Ignore SEND_RECEIVE for the moment ---*/
                if (bound[iMarker][0]->GetVTK_Type() != TBOX::VERTEX) 
                {
                    Grid_Marker = config->GetMarker_All_TagBound(iMarker);
                    output_file << "MARKER_TAG= " << Grid_Marker << std::endl;
                    output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] << std::endl;

                    if (nDim == 2) 
                    {
                        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) 
                        {
                            output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t";
                            for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
                                output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t";
                            output_file << iElem_Bound << std::endl;
                        }
                    }

                    if (nDim == 3) 
                    {
                        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) 
                        {
                            output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t";
                            for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
                                output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t";
                            output_file << iElem_Bound << std::endl;
                        }
                    }

                }
                else if (bound[iMarker][0]->GetVTK_Type() == TBOX::VERTEX) 
                {
                    output_file << "MARKER_TAG= SEND_RECEIVE" << std::endl;
                    output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] << std::endl;
                    if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << std::endl;
                    if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << std::endl;

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) 
                    {
                        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
                            bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
                            bound[iMarker][iElem_Bound]->GetRotation_Type() << std::endl;
                    }

                }
            }

            /*--- Get the total number of periodic transformations ---*/

            nPeriodic = config->GetnPeriodicIndex();
            output_file << "NPERIODIC= " << nPeriodic << std::endl;

            /*--- From iPeriodic obtain the iMarker ---*/
            for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) 
            {
                /*--- Retrieve the supplied periodic information. ---*/
                center = config->GetPeriodicCenter(iPeriodic);
                angles = config->GetPeriodicRotation(iPeriodic);
                transl = config->GetPeriodicTranslate(iPeriodic);

                output_file << "PERIODIC_INDEX= " << iPeriodic << std::endl;
                output_file << center[0] << "\t" << center[1] << "\t" << center[2] << std::endl;
                output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << std::endl;
                output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << std::endl;
            }

            output_file.close();
        }

        void GEOM_GeometryPhysical::SetCoord_Smoothing(unsigned short val_nSmooth, double val_smooth_coeff, TBOX::TBOX_Config *config)
        {
            unsigned short iSmooth, nneigh, iMarker;
            double *Coord_Old, *Coord_Sum, *Coord, *Coord_i, *Coord_j, Position_Plane = 0.0;
            unsigned long iEdge, iPoint, jPoint, iVertex;
            double eps = 1E-6;
            bool NearField = false;

            Coord = new double[nDim];

            for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
            {
                double *Coord = node[iPoint]->GetCoord();
                node[iPoint]->SetCoord_Old(Coord);
            }

            /*--- Jacobi iterations ---*/
            for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++)
            {

                for (iPoint = 0; iPoint < nPoint; iPoint++)
                    node[iPoint]->SetCoord_SumZero();

                /*--- Loop over Interior edges ---*/
                for (iEdge = 0; iEdge < nEdge; iEdge++)
                {
                    iPoint = edge[iEdge]->GetNode(0);
                    Coord_i = node[iPoint]->GetCoord();

                    jPoint = edge[iEdge]->GetNode(1);
                    Coord_j = node[jPoint]->GetCoord();

                    /*--- Accumulate nearest neighbor Coord to Res_sum for each variable ---*/
                    node[iPoint]->AddCoord_Sum(Coord_j);
                    node[jPoint]->AddCoord_Sum(Coord_i);
                }

                /*--- Loop over all mesh points (Update Coords with averaged sum) ---*/
                for (iPoint = 0; iPoint < nPoint; iPoint++)
                {
                    nneigh = node[iPoint]->GetnPoint();
                    Coord_Sum = node[iPoint]->GetCoord_Sum();
                    Coord_Old = node[iPoint]->GetCoord_Old();

                    if (nDim == 2)
                    {
                        Coord[0] = (Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) / (1.0 + val_smooth_coeff*double(nneigh));
                        Coord[1] = (Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) / (1.0 + val_smooth_coeff*double(nneigh));
                        if ((NearField) && ((Coord_Old[1] > Position_Plane - eps) && (Coord_Old[1] < Position_Plane + eps)))
                            Coord[1] = Coord_Old[1];
                    }

                    if (nDim == 3)
                    {
                        Coord[0] = (Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) / (1.0 + val_smooth_coeff*double(nneigh));
                        Coord[1] = (Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) / (1.0 + val_smooth_coeff*double(nneigh));
                        Coord[2] = (Coord_Old[2] + val_smooth_coeff*Coord_Sum[2]) / (1.0 + val_smooth_coeff*double(nneigh));
                        if ((NearField) && ((Coord_Old[2] > Position_Plane - eps) && (Coord_Old[2] < Position_Plane + eps)))
                            Coord[2] = Coord_Old[2];
                    }

                    node[iPoint]->SetCoord(Coord);
                }

                /*--- Copy boundary values ---*/
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        Coord_Old = node[iPoint]->GetCoord_Old();
                        node[iPoint]->SetCoord(Coord_Old);
                    }
            }

            delete[] Coord;
        }

        bool GEOM_GeometryPhysical::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem, unsigned short &face_second_elem) 
        {
            /*--- Find repeated nodes between two elements to identify the common face ---*/
            unsigned long iPoint = 0, jPoint = 0;
            unsigned short face_node, iFace, iNode, jNode, nNodesFace;
            std::vector<unsigned long> CommonPoints, PointFaceFirst, PointFaceSecond;
            std::vector<unsigned long>::iterator IterPoint;
            std::pair<std::vector <unsigned long>::iterator, std::vector <unsigned long>::iterator> mypair1, mypair2;
            bool face_first_found = false, face_second_found = false;

            if (first_elem == second_elem)
                return 0;

            for (iNode = 0; iNode < elem[first_elem]->GetnNodes(); iNode++) 
            {
                iPoint = elem[first_elem]->GetNode(iNode);
                for (jNode = 0; jNode < elem[second_elem]->GetnNodes(); jNode++) 
                {
                    jPoint = elem[second_elem]->GetNode(jNode);
                    if (iPoint == jPoint) 
                    {
                        CommonPoints.push_back(iPoint);
                        break;
                    }
                }
            }

            /*--- Sort point in face and check that the list is unique ---*/
            sort(CommonPoints.begin(), CommonPoints.end());
            IterPoint = unique(CommonPoints.begin(), CommonPoints.end());
            CommonPoints.resize(distance(CommonPoints.begin(), IterPoint));

            if (CommonPoints.size() <= 1) 
                return false;

            /*--- Search the secuence in the first element ---*/
            for (iFace = 0; iFace < elem[first_elem]->GetnFaces(); iFace++)
            {
                nNodesFace = elem[first_elem]->GetnNodesFace(iFace);
                for (iNode = 0; iNode < nNodesFace; iNode++) 
                {
                    face_node = elem[first_elem]->GetFaces(iFace, iNode);
                    PointFaceFirst.push_back(elem[first_elem]->GetNode(face_node));
                }

                /*--- Sort face_poin to perform comparison ---*/
                sort(PointFaceFirst.begin(), PointFaceFirst.end());

                /*--- List comparison ---*/
                mypair1 = std::mismatch(PointFaceFirst.begin(), PointFaceFirst.end(), CommonPoints.begin());
                if (mypair1.first == PointFaceFirst.end()) 
                {
                    face_first_elem = iFace;
                    face_first_found = true;
                    break;
                }

                PointFaceFirst.erase(PointFaceFirst.begin(), PointFaceFirst.end());
            }

            /*--- Search the secuence in the second element ---*/
            for (iFace = 0; iFace < elem[second_elem]->GetnFaces(); iFace++) 
            {
                nNodesFace = elem[second_elem]->GetnNodesFace(iFace);
                for (iNode = 0; iNode < nNodesFace; iNode++) 
                {
                    face_node = elem[second_elem]->GetFaces(iFace, iNode);
                    PointFaceSecond.push_back(elem[second_elem]->GetNode(face_node));
                }

                /*--- Sort face_poin to perform comparison ---*/
                sort(PointFaceSecond.begin(), PointFaceSecond.end());

                /*--- List comparison ---*/
                mypair2 = std::mismatch(PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin());
                if (mypair2.first == PointFaceSecond.end()) 
                //if (equal(PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin()))
                {
                    face_second_elem = iFace;
                    face_second_found = true;
                    break;
                }

                PointFaceSecond.erase(PointFaceSecond.begin(), PointFaceSecond.end());
            }

            if (face_first_found && face_second_found) 
                return true;
            else 
                return false;
        }

        void GEOM_GeometryPhysical::SetTecPlot(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file) 
        {
            unsigned long iElem, iPoint;
            unsigned short iDim;
            std::ofstream Tecplot_File;

            /*--- Open the tecplot file and write the header ---*/

            if (new_file) 
            {
                Tecplot_File.open(mesh_filename, std::ios::out);
                Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << std::endl;
                if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << std::endl;
                if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << std::endl;
            }
            else Tecplot_File.open(mesh_filename, std::ios::out | std::ios::app);

            Tecplot_File << "ZONE T= ";
            if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
            else Tecplot_File << "\"Deformed grid\", C=RED, ";
            Tecplot_File << "NODES= " << nPoint << ", ELEMENTS= " << nElem << ", DATAPACKING= POINT";
            if (nDim == 2) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL" << std::endl;
            if (nDim == 3) Tecplot_File << ", ZONETYPE= FEBRICK" << std::endl;

            /*--- Adding coordinates ---*/

            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                for (iDim = 0; iDim < nDim; iDim++)
                    Tecplot_File << std::scientific << node[iPoint]->GetCoord(iDim) << "\t";
                Tecplot_File << "\n";
            }

            /*--- Adding conectivity ---*/

            for (iElem = 0; iElem < nElem; iElem++) 
            {
                if (elem[iElem]->GetVTK_Type() == TBOX::TRIANGLE) 
                {
                    Tecplot_File <<
                        elem[iElem]->GetNode(0) + 1 << " " << elem[iElem]->GetNode(1) + 1 << " " <<
                        elem[iElem]->GetNode(2) + 1 << " " << elem[iElem]->GetNode(2) + 1 << std::endl;
                }
                if (elem[iElem]->GetVTK_Type() == TBOX::RECTANGLE) 
                {
                    Tecplot_File <<
                        elem[iElem]->GetNode(0) + 1 << " " << elem[iElem]->GetNode(1) + 1 << " " <<
                        elem[iElem]->GetNode(2) + 1 << " " << elem[iElem]->GetNode(3) + 1 << std::endl;
                }
                if (elem[iElem]->GetVTK_Type() == TBOX::TETRAHEDRON) 
                {
                    Tecplot_File <<
                        elem[iElem]->GetNode(0) + 1 << " " << elem[iElem]->GetNode(1) + 1 << " " <<
                        elem[iElem]->GetNode(2) + 1 << " " << elem[iElem]->GetNode(2) + 1 << " " <<
                        elem[iElem]->GetNode(3) + 1 << " " << elem[iElem]->GetNode(3) + 1 << " " <<
                        elem[iElem]->GetNode(3) + 1 << " " << elem[iElem]->GetNode(3) + 1 << std::endl;
                }
                if (elem[iElem]->GetVTK_Type() == TBOX::HEXAHEDRON) 
                {
                    Tecplot_File <<
                        elem[iElem]->GetNode(0) + 1 << " " << elem[iElem]->GetNode(1) + 1 << " " <<
                        elem[iElem]->GetNode(2) + 1 << " " << elem[iElem]->GetNode(3) + 1 << " " <<
                        elem[iElem]->GetNode(4) + 1 << " " << elem[iElem]->GetNode(5) + 1 << " " <<
                        elem[iElem]->GetNode(6) + 1 << " " << elem[iElem]->GetNode(7) + 1 << std::endl;
                }
                if (elem[iElem]->GetVTK_Type() == TBOX::PYRAMID) 
                {
                    Tecplot_File <<
                        elem[iElem]->GetNode(0) + 1 << " " << elem[iElem]->GetNode(1) + 1 << " " <<
                        elem[iElem]->GetNode(2) + 1 << " " << elem[iElem]->GetNode(3) + 1 << " " <<
                        elem[iElem]->GetNode(4) + 1 << " " << elem[iElem]->GetNode(4) + 1 << " " <<
                        elem[iElem]->GetNode(4) + 1 << " " << elem[iElem]->GetNode(4) + 1 << std::endl;
                }
                if (elem[iElem]->GetVTK_Type() == TBOX::PRISM) 
                {
                    Tecplot_File <<
                        elem[iElem]->GetNode(0) + 1 << " " << elem[iElem]->GetNode(1) + 1 << " " <<
                        elem[iElem]->GetNode(1) + 1 << " " << elem[iElem]->GetNode(2) + 1 << " " <<
                        elem[iElem]->GetNode(3) + 1 << " " << elem[iElem]->GetNode(4) + 1 << " " <<
                        elem[iElem]->GetNode(4) + 1 << " " << elem[iElem]->GetNode(5) + 1 << std::endl;
                }
            }

            Tecplot_File.close();
        }
        
        void GEOM_GeometryPhysical::SetBoundTecPlot(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file, TBOX::TBOX_Config *config) 
        {
            std::ofstream Tecplot_File;
            unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
            unsigned short Coord_i, iMarker;

            /*--- It is important to do a renumbering to don't add points
            that do not belong to the surfaces ---*/
            PointSurface = new unsigned long[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                if (node[iPoint]->GetBoundary()) 
                {
                    PointSurface[iPoint] = nPointSurface;
                    nPointSurface++;
                }

            /*--- Compute the total number of elements ---*/

            Total_nElem_Bound = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                if (config->GetMarker_All_Plotting(iMarker) == TBOX::YES) 
                {
                    Total_nElem_Bound += nElem_Bound[iMarker];
                }
            }

            /*--- Open the tecplot file and write the header ---*/

            if (new_file) 
            {
                Tecplot_File.open(mesh_filename, std::ios::out);
                Tecplot_File << "TITLE= \"Visualization of the surface grid\"" << std::endl;
                if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << std::endl;
                if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << std::endl;
            }
            else Tecplot_File.open(mesh_filename, std::ios::out | std::ios::app);

            if (Total_nElem_Bound != 0) 
            {
                /*--- Write the header of the file ---*/
                Tecplot_File << "ZONE T= ";
                if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
                else Tecplot_File << "\"Deformed grid\", C=RED, ";
                Tecplot_File << "NODES= " << nPointSurface << ", ELEMENTS= " << Total_nElem_Bound << ", DATAPACKING= POINT";
                if (nDim == 2) Tecplot_File << ", ZONETYPE= FELINESEG" << std::endl;
                if (nDim == 3) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL" << std::endl;

                /*--- Only write the coordiantes of the points that are on the surfaces ---*/

                if (nDim == 3) 
                {
                    for (iPoint = 0; iPoint < nPoint; iPoint++)
                        if (node[iPoint]->GetBoundary())
                        {
                            for (Coord_i = 0; Coord_i < nDim - 1; Coord_i++)
                                Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
                            Tecplot_File << node[iPoint]->GetCoord(nDim - 1) << "\n";
                        }
                }
                else 
                {
                    for (iPoint = 0; iPoint < nPoint; iPoint++)
                        if (node[iPoint]->GetBoundary()) 
                        {
                            for (Coord_i = 0; Coord_i < nDim; Coord_i++)
                                Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
                            Tecplot_File << "\n";
                        }
                }

                /*--- Write the cells using the new numbering ---*/

                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    if (config->GetMarker_All_Plotting(iMarker) == TBOX::YES)
                        for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) 
                        {
                            if (nDim == 2) {
                                Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)] + 1 << " "
                                    << PointSurface[bound[iMarker][iElem]->GetNode(1)] + 1 << std::endl;
                            }
                            if (nDim == 3) {
                                if (bound[iMarker][iElem]->GetnNodes() == 3) 
                                {
                                    Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)] + 1 << " "
                                        << PointSurface[bound[iMarker][iElem]->GetNode(1)] + 1 << " "
                                        << PointSurface[bound[iMarker][iElem]->GetNode(2)] + 1 << " "
                                        << PointSurface[bound[iMarker][iElem]->GetNode(2)] + 1 << std::endl;
                                }
                                if (bound[iMarker][iElem]->GetnNodes() == 4) 
                                {
                                    Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)] + 1 << " "
                                        << PointSurface[bound[iMarker][iElem]->GetNode(1)] + 1 << " "
                                        << PointSurface[bound[iMarker][iElem]->GetNode(2)] + 1 << " "
                                        << PointSurface[bound[iMarker][iElem]->GetNode(3)] + 1 << std::endl;
                                }
                            }
                        }
            }
            else 
            {
                /*--- No elements in the surface ---*/
                if (nDim == 2) 
                {
                    Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG" << std::endl;
                    Tecplot_File << "0.0 0.0" << std::endl;
                    Tecplot_File << "1 1" << std::endl;
                }
                if (nDim == 3) 
                {
                    Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << std::endl;
                    Tecplot_File << "0.0 0.0 0.0" << std::endl;
                    Tecplot_File << "1 1 1 1" << std::endl;
                }
            }

            /*--- Dealocate memory and close the file ---*/
            delete[] PointSurface;
            Tecplot_File.close();
        }

        void GEOM_GeometryPhysical::SetBoundSTL(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file, TBOX::TBOX_Config *config)
        {
            std::ofstream STL_File;
            unsigned long this_node, iNode, nNode, iElem;
            unsigned short iDim, iMarker;
            double p[3] = { 0.0, 0.0, 0.0 }, u[3] = { 0.0, 0.0, 0.0 }, v[3] = { 0.0, 0.0, 0.0 }, n[3] = { 0.0, 0.0, 0.0 }, a;

            /*---	STL format:
            solid NAME
            ...
            facet normal 0.00 0.00 1.00
            outer loop
            vertex  2.00  2.00  0.00
            vertex -1.00  1.00  0.00
            vertex  0.00 -1.00  0.00
            endloop
            endfacet
            ...
            end solid
            --- */

            /*--- Open the STL file ---*/

            if (new_file) STL_File.open(mesh_filename, std::ios::out);
            else STL_File.open(mesh_filename, std::ios::out | std::ios::app);

            /*--- Write the header of the file ---*/
            STL_File << "solid surface_mesh" << std::endl;

            /*--- Write facets of surface markers ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_Plotting(iMarker) == TBOX::YES)
                    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                    {
                        /*--- number of nodes for this elemnt ---*/
                        nNode = bound[iMarker][iElem]->GetnNodes();

                        /*--- Calculate Normal std::vector ---*/
                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            p[0] = node[bound[iMarker][iElem]->GetNode(0)]->GetCoord(iDim);
                            p[1] = node[bound[iMarker][iElem]->GetNode(1)]->GetCoord(iDim);
                            p[2] = node[bound[iMarker][iElem]->GetNode(nNode - 1)]->GetCoord(iDim);
                            u[iDim] = p[1] - p[0];
                            v[iDim] = p[2] - p[0];
                        }

                        n[0] = u[1] * v[2] - u[2] * v[1];
                        n[1] = u[2] * v[0] - u[0] * v[2];
                        n[2] = u[0] * v[1] - u[1] * v[0];
                        a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

                        /*--- Print normal std::vector ---*/
                        STL_File << "  facet normal ";
                        for (iDim = 0; iDim < nDim; iDim++)
                        {
                            STL_File << n[iDim] / a << " ";
                        }
                        STL_File << std::endl;

                        /*--- STL Facet Loop --*/

                        STL_File << "    outer loop" << std::endl;

                        /*--- Print Nodes for Facet ---*/
                        for (iNode = 0; iNode < nNode; iNode++)
                        {
                            this_node = bound[iMarker][iElem]->GetNode(iNode);
                            STL_File << "      vertex ";
                            for (iDim = 0; iDim < nDim; iDim++)
                                STL_File << node[this_node]->GetCoord(iDim) << " ";
                            if (nDim == 2)
                                STL_File << 0.0 << " ";
                            STL_File << std::endl;
                        }
                        STL_File << "    endloop" << std::endl;
                        STL_File << "  endfacet" << std::endl;
                    }

            /*--- Done with Surface Mesh ---*/
            STL_File << "endsolid" << std::endl;

            /*--- Close the file ---*/
            STL_File.close();
        }

        void GEOM_GeometryPhysical::SetColorGrid(TBOX::TBOX_Config *config)
        {
#ifdef HAVE_MPI

#ifdef HAVE_METIS
            unsigned long iPoint, iElem, iElem_Triangle, iElem_Tetrahedron, nElem_Triangle,
                nElem_Tetrahedron;
            unsigned short iNode;
            idx_t ne = 0, nn, *elmnts = NULL, etype, *epart = NULL, *npart = NULL, numflag, nparts, edgecut, *eptr;
            int rank, size;

            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            if (size != SINGLE_ZONE)
                cout << endl << "---------------------------- Grid partitioning --------------------------" << endl;

            unsigned short nDomain = size;

            nElem_Triangle = 0;
            nElem_Tetrahedron = 0;
            for (iElem = 0; iElem < GetnElem(); iElem++) {
                if (elem[iElem]->GetVTK_Type() == TRIANGLE) nElem_Triangle = nElem_Triangle + 1;
                if (elem[iElem]->GetVTK_Type() == RECTANGLE) nElem_Triangle = nElem_Triangle + 2;
                if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) nElem_Tetrahedron = nElem_Tetrahedron + 1;
                if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) nElem_Tetrahedron = nElem_Tetrahedron + 5;
                if (elem[iElem]->GetVTK_Type() == PYRAMID) nElem_Tetrahedron = nElem_Tetrahedron + 2;
                if (elem[iElem]->GetVTK_Type() == PRISM) nElem_Tetrahedron = nElem_Tetrahedron + 3;
            }

            if (GetnDim() == 2) {
                ne = nElem_Triangle;
                elmnts = new idx_t[ne * 3];
                etype = 1;
            }
            if (GetnDim() == 3) {
                ne = nElem_Tetrahedron;
                elmnts = new idx_t[ne * 4];
                etype = 2;
            }

            nn = nPoint;
            numflag = 0;
            nparts = nDomain;
            epart = new idx_t[ne];
            npart = new idx_t[nn];
            eptr = new idx_t[ne + 1];

            /*--- Initialize the color std::vector ---*/

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                node[iPoint]->SetColor(0);

            if (nparts > 1) {

                iElem_Triangle = 0; iElem_Tetrahedron = 0;
                for (iElem = 0; iElem < GetnElem(); iElem++) {
                    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
                        elmnts[3 * iElem_Triangle + 0] = elem[iElem]->GetNode(0);
                        elmnts[3 * iElem_Triangle + 1] = elem[iElem]->GetNode(1);
                        elmnts[3 * iElem_Triangle + 2] = elem[iElem]->GetNode(2);
                        eptr[iElem_Triangle] = 3 * iElem_Triangle;
                        iElem_Triangle++;
                    }
                    if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
                        elmnts[3 * iElem_Triangle + 0] = elem[iElem]->GetNode(0);
                        elmnts[3 * iElem_Triangle + 1] = elem[iElem]->GetNode(1);
                        elmnts[3 * iElem_Triangle + 2] = elem[iElem]->GetNode(2);
                        eptr[iElem_Triangle] = 3 * iElem_Triangle;
                        iElem_Triangle++;
                        elmnts[3 * iElem_Triangle + 0] = elem[iElem]->GetNode(0);
                        elmnts[3 * iElem_Triangle + 1] = elem[iElem]->GetNode(2);
                        elmnts[3 * iElem_Triangle + 2] = elem[iElem]->GetNode(3);
                        eptr[iElem_Triangle] = 3 * iElem_Triangle;
                        iElem_Triangle++;
                    }
                    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(1);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(3);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                    }
                    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(1);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(5);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(3);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(7);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(5);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(7);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(4);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(7);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(5);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(6);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(7);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(5);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                    }
                    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(1);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(4);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(3);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(4);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                    }
                    if (elem[iElem]->GetVTK_Type() == PRISM) {
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(1);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(4);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(2);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(0);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(2);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(3);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(4);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                        elmnts[4 * iElem_Tetrahedron + 0] = elem[iElem]->GetNode(3);
                        elmnts[4 * iElem_Tetrahedron + 1] = elem[iElem]->GetNode(4);
                        elmnts[4 * iElem_Tetrahedron + 2] = elem[iElem]->GetNode(5);
                        elmnts[4 * iElem_Tetrahedron + 3] = elem[iElem]->GetNode(2);
                        eptr[iElem_Tetrahedron] = 4 * iElem_Tetrahedron;
                        iElem_Tetrahedron++;
                    }
                }

                /*--- Add final value to element pointer array ---*/

                if (GetnDim() == 2) eptr[ne] = 3 * ne;
                else eptr[ne] = 4 * ne;

                METIS_PartMeshNodal(&ne, &nn, eptr, elmnts, NULL, NULL, &nparts, NULL, NULL, &edgecut, epart, npart);

                cout << "Finished partitioning using METIS. (" << edgecut << " edge cuts)." << endl;

                for (iPoint = 0; iPoint < nPoint; iPoint++)
                    node[iPoint]->SetColor(npart[iPoint]);
            }

            delete[] epart;
            delete[] npart;
            delete[] elmnts;
            delete[] eptr;
#endif

#endif
        }

        void GEOM_GeometryPhysical::SetColorGrid_Parallel(TBOX::TBOX_Config *config) 
        {
            /*--- Initialize the color std::vector ---*/
            for (unsigned long iPoint = 0; iPoint < local_node; iPoint++) node[iPoint]->SetColor(0);
            /*--- This routine should only ever be called if we have parallel support
            with MPI and have the ParMETIS library compiled and linked. ---*/

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS

            unsigned long iPoint;
            int rank, size;
            MPI_Comm comm = MPI_COMM_WORLD;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            /*--- Only call ParMETIS if we have more than one rank to avoid errors ---*/

            if (size > SINGLE_NODE) {

                /*--- Create some structures that ParMETIS needs for partitioning. ---*/

                idx_t numflag, nparts, edgecut, vwgt, adjwgt, wgtflag, ncon, ncommonnodes;
                idx_t *vtxdist = new idx_t[size + 1];
                idx_t *xadj_l = new idx_t[xadj_size];
                idx_t *adjacency_l = new idx_t[adjacency_size];
                idx_t *elmwgt = new idx_t[local_node];
                idx_t *part = new idx_t[local_node];

                real_t ubvec;
                real_t *tpwgts = new real_t[size];

                /*--- Some recommended defaults for the various ParMETIS options. ---*/

                wgtflag = 0;
                numflag = 0;
                ncon = 1;
                ubvec = 1.05;
                nparts = (idx_t)size;
                idx_t options[METIS_NOPTIONS];
                METIS_SetDefaultOptions(options);
                options[1] = 0;

                /*--- Fill the necessary ParMETIS data arrays ---*/

                for (int i = 0; i < size; i++) {
                    tpwgts[i] = 1.0 / ((real_t)size);
                }

                vtxdist[0] = 0;
                for (int i = 0; i < size; i++) {
                    vtxdist[i + 1] = (idx_t)ending_node[i];
                }

                for (int i = 0; i < xadj_size; i++) {
                    xadj_l[i] = (idx_t)xadj[i];
                }

                for (int i = 0; i < adjacency_size; i++) {
                    adjacency_l[i] = (idx_t)adjacency[i];
                }

                /*--- Calling ParMETIS ---*/
                if (rank == MASTER_NODE) cout << "Calling ParMETIS..." << endl;
                ParMETIS_V3_PartKway(vtxdist, xadj_l, adjacency_l, NULL, NULL, &wgtflag,
                    &numflag, &ncon, &nparts, tpwgts, &ubvec, options,
                    &edgecut, part, &comm);
                if (rank == MASTER_NODE) {
                    cout << "Finished partitioning using ParMETIS (";
                    cout << edgecut << " edge cuts)." << endl;
                }

                /*--- Store the results of the partitioning (note that this is local
                since each processor is calling ParMETIS in parallel and storing the
                results for its initial piece of the grid. ---*/

                for (iPoint = 0; iPoint < local_node; iPoint++) {
                    node[iPoint]->SetColor(part[iPoint]);
                }

                /*--- Free all memory needed for the ParMETIS structures ---*/

                delete[] vtxdist;
                delete[] xadj_l;
                delete[] adjacency_l;
                delete[] elmwgt;
                delete[] part;
                delete[] tpwgts;

            }

            /*--- Delete the memory from the geometry class that carried the
            adjacency structure. ---*/

            delete[] xadj;
            delete[] adjacency;

#endif
#endif
        }

        void GEOM_GeometryPhysical::GetQualityStatistics(double *statistics) 
        {
            unsigned long jPoint, Point_2, Point_3, iElem;
            double *Coord_j, *Coord_2, *Coord_3;
            unsigned short iDim;

            statistics[0] = 1e06;
            statistics[1] = 0;

            /*--- Loop interior edges ---*/
            for (iElem = 0; iElem < this->GetnElem(); iElem++) 
            {
                if ((this->GetnDim() == 2) && (elem[iElem]->GetVTK_Type() == TBOX::TRIANGLE)) 
                {
                    jPoint = elem[iElem]->GetNode(0); Coord_j = node[jPoint]->GetCoord();
                    Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
                    Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();

                    /*--- Compute sides of the triangle ---*/
                    double a = 0, b = 0, c = 0;
                    for (iDim = 0; iDim < nDim; iDim++) 
                    {
                        a += (Coord_2[iDim] - Coord_j[iDim])*(Coord_2[iDim] - Coord_j[iDim]);
                        b += (Coord_3[iDim] - Coord_j[iDim])*(Coord_3[iDim] - Coord_j[iDim]);
                        c += (Coord_3[iDim] - Coord_2[iDim])*(Coord_3[iDim] - Coord_2[iDim]);
                    }
                    a = sqrt(a); b = sqrt(b); c = sqrt(c);

                    /*--- Compute semiperimeter (s) and area ---*/
                    double s = 0.5*(a + b + c);
                    double Area = sqrt(s*(s - a)*(s - b)*(s - c));

                    /*--- Compute radius of the circumcircle (R) and of the incircle (r) ---*/
                    double R = (a*b*c) / (4.0*Area);
                    double r = Area / s;
                    double roR = r / R;

                    /*--- Update statistics ---*/
                    if (roR < statistics[0])
                        statistics[0] = roR;
                    statistics[1] += roR;
                }
            }
            statistics[1] /= this->GetnElem();
        }

        void GEOM_GeometryPhysical::SetRotationalVelocity(TBOX::TBOX_Config *config) 
        {
            unsigned long iPoint;
            double RotVel[3], Distance[3], *Coord, Center[3], Omega[3], L_Ref;

            int rank = TBOX::MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Center of rotation & angular velocity std::vector from config ---*/
            Center[0] = config->GetMotion_Origin_X(TBOX::ZONE_0);
            Center[1] = config->GetMotion_Origin_Y(TBOX::ZONE_0);
            Center[2] = config->GetMotion_Origin_Z(TBOX::ZONE_0);
            Omega[0] = config->GetRotation_Rate_X(TBOX::ZONE_0) / config->GetOmega_Ref();
            Omega[1] = config->GetRotation_Rate_Y(TBOX::ZONE_0) / config->GetOmega_Ref();
            Omega[2] = config->GetRotation_Rate_Z(TBOX::ZONE_0) / config->GetOmega_Ref();
            L_Ref = config->GetLength_Ref();

            /*--- Print some information to the console ---*/

            if (rank == TBOX::MASTER_NODE)
            {
                std::cout << " Rotational origin (x, y, z): ( " << Center[0] << ", " << Center[1];
                std::cout << ", " << Center[2] << " )" << std::endl;
                std::cout << " Angular velocity about x, y, z axes: ( " << Omega[0] << ", ";
                std::cout << Omega[1] << ", " << Omega[2] << " ) rad/s" << std::endl;
            }

            /*--- Loop over all nodes and set the rotational velocity ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                /*--- Get the coordinates of the current node ---*/
                Coord = node[iPoint]->GetCoord();

                /*--- Calculate the non-dim. distance from the rotation center ---*/
                Distance[0] = (Coord[0] - Center[0]) / L_Ref;
                Distance[1] = (Coord[1] - Center[1]) / L_Ref;
                Distance[2] = (Coord[2] - Center[2]) / L_Ref;

                /*--- Calculate the angular velocity as omega X r ---*/
                RotVel[0] = Omega[1] * (Distance[2]) - Omega[2] * (Distance[1]);
                RotVel[1] = Omega[2] * (Distance[0]) - Omega[0] * (Distance[2]);
                RotVel[2] = Omega[0] * (Distance[1]) - Omega[1] * (Distance[0]);

                /*--- Store the grid velocity at this node ---*/
                node[iPoint]->SetGridVel(RotVel);
            }
        }

        void GEOM_GeometryPhysical::SetGridVelocity(TBOX::TBOX_Config *config, unsigned long iter)
        {
            /*--- Local variables ---*/
            double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
            double TimeStep, GridVel = 0.0;
            unsigned long iPoint;
            unsigned short iDim;

            /*--- Compute the velocity of each node in the volume mesh ---*/
            for (iPoint = 0; iPoint < GetnPoint(); iPoint++) 
            {
                /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
                Coord_nM1 = node[iPoint]->GetCoord_n1();
                Coord_n = node[iPoint]->GetCoord_n();
                Coord_nP1 = node[iPoint]->GetCoord();

                /*--- Unsteady time step ---*/
                TimeStep = config->GetDelta_UnstTimeND();

                /*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
                for (iDim = 0; iDim < nDim; iDim++) 
                {
                    if (config->GetUnsteady_Simulation() == TBOX::DT_STEPPING_1ST)
                        GridVel = (Coord_nP1[iDim] - Coord_n[iDim]) / TimeStep;
                    if (config->GetUnsteady_Simulation() == TBOX::DT_STEPPING_2ND)
                        GridVel = (3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                        + 1.0*Coord_nM1[iDim]) / (2.0*TimeStep);

                    /*--- Store grid velocity for this point ---*/
                    node[iPoint]->SetGridVel(iDim, GridVel);
                }
            }
        }

        void GEOM_GeometryPhysical::Set_MPI_Coord(TBOX::TBOX_Config *config)
        {
            unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
            unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_vector, nBufferR_vector;
            double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_Coord = NULL, *Buffer_Send_Coord = NULL, *Coord = NULL, *newCoord = NULL;

            newCoord = new double[nDim];

#ifdef HAVE_MPI
            int send_to, receive_from;
            MPI_Status status;
#endif

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            {

                if ((config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE) &&
                    (config->GetMarker_All_SendRecv(iMarker) > 0))
                {

                    MarkerS = iMarker;  MarkerR = iMarker + 1;

#ifdef HAVE_MPI
                    send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
                    receive_from = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
#endif

                    nVertexS = nVertex[MarkerS];
                    nVertexR = nVertex[MarkerR];
                    nBufferS_vector = nVertexS*nDim;
                    nBufferR_vector = nVertexR*nDim;

                    /*--- Allocate Receive and send buffers  ---*/

                    Buffer_Receive_Coord = new double[nBufferR_vector];
                    Buffer_Send_Coord = new double[nBufferS_vector];

                    /*--- Copy the coordinates that should be sended ---*/

                    for (iVertex = 0; iVertex < nVertexS; iVertex++)
                    {
                        iPoint = vertex[MarkerS][iVertex]->GetNode();
                        Coord = node[iPoint]->GetCoord();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_Coord[iDim*nVertexS + iVertex] = Coord[iDim];
                    }

#ifdef HAVE_MPI
                    /*--- Send/Receive information using Sendrecv ---*/
                    MPI_Sendrecv(Buffer_Send_Coord, nBufferS_std::vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_Coord, nBufferR_std::vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else

                    /*--- Receive information without MPI ---*/
                    for (iVertex = 0; iVertex < nVertexR; iVertex++)
                    {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Receive_Coord[iDim*nVertexR + iVertex] = Buffer_Send_Coord[iDim*nVertexR + iVertex];
                    }

#endif

                    /*--- Deallocate send buffer ---*/

                    delete[] Buffer_Send_Coord;

                    /*--- Do the coordinate transformation ---*/

                    for (iVertex = 0; iVertex < nVertexR; iVertex++)
                    {

                        /*--- Find point and its type of transformation ---*/

                        iPoint = vertex[MarkerR][iVertex]->GetNode();
                        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();

                        /*--- Retrieve the supplied periodic information. ---*/

                        angles = config->GetPeriodicRotation(iPeriodic_Index);

                        /*--- Store angles separately for clarity. ---*/

                        theta = angles[0];   phi = angles[1];     psi = angles[2];
                        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
                        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

                        /*--- Compute the rotation matrix. Note that the implicit
                        ordering is rotation about the x-axis, y-axis,
                        then z-axis. Note that this is the transpose of the matrix
                        used during the preprocessing stage. ---*/

                        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

                        /*--- Copy coordinates before performing transformation. ---*/

                        for (iDim = 0; iDim < nDim; iDim++)
                            newCoord[iDim] = Buffer_Receive_Coord[iDim*nVertexR + iVertex];

                        /*--- Rotate the coordinates. ---*/

                        if (nDim == 2)
                        {
                            newCoord[0] = (rotMatrix[0][0] * Buffer_Receive_Coord[0 * nVertexR + iVertex] +
                                rotMatrix[0][1] * Buffer_Receive_Coord[1 * nVertexR + iVertex]);
                            newCoord[1] = (rotMatrix[1][0] * Buffer_Receive_Coord[0 * nVertexR + iVertex] +
                                rotMatrix[1][1] * Buffer_Receive_Coord[1 * nVertexR + iVertex]);
                        }
                        else
                        {
                            newCoord[0] = (rotMatrix[0][0] * Buffer_Receive_Coord[0 * nVertexR + iVertex] +
                                rotMatrix[0][1] * Buffer_Receive_Coord[1 * nVertexR + iVertex] +
                                rotMatrix[0][2] * Buffer_Receive_Coord[2 * nVertexR + iVertex]);
                            newCoord[1] = (rotMatrix[1][0] * Buffer_Receive_Coord[0 * nVertexR + iVertex] +
                                rotMatrix[1][1] * Buffer_Receive_Coord[1 * nVertexR + iVertex] +
                                rotMatrix[1][2] * Buffer_Receive_Coord[2 * nVertexR + iVertex]);
                            newCoord[2] = (rotMatrix[2][0] * Buffer_Receive_Coord[0 * nVertexR + iVertex] +
                                rotMatrix[2][1] * Buffer_Receive_Coord[1 * nVertexR + iVertex] +
                                rotMatrix[2][2] * Buffer_Receive_Coord[2 * nVertexR + iVertex]);
                        }

                        /*--- Copy transformed coordinates back into buffer. ---*/
                        for (iDim = 0; iDim < nDim; iDim++)
                            node[iPoint]->SetCoord(iDim, newCoord[iDim]);

                    }

                    /*--- Deallocate receive buffer. ---*/
                    delete[] Buffer_Receive_Coord;
                }
            }
            delete[] newCoord;
        }

        void GEOM_GeometryPhysical::Set_MPI_GridVel(TBOX::TBOX_Config *config) 
        {
            unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
            unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_vector, nBufferR_vector;
            double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_GridVel = NULL, *Buffer_Send_GridVel = NULL, *GridVel = NULL, *newGridVel = NULL;

            newGridVel = new double[nDim];

#ifdef HAVE_MPI
            int send_to, receive_from;
            MPI_Status status;
#endif

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            {
                if ((config->GetMarker_All_KindBC(iMarker) == TBOX::SEND_RECEIVE) &&
                    (config->GetMarker_All_SendRecv(iMarker) > 0)) 
                {
                    MarkerS = iMarker;  MarkerR = iMarker + 1;

#ifdef HAVE_MPI
                    send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
                    receive_from = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
#endif

                    nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
                    nBufferS_vector = nVertexS*nDim;        
                    nBufferR_vector = nVertexR*nDim;

                    /*--- Allocate Receive and send buffers  ---*/

                    Buffer_Receive_GridVel = new double[nBufferR_vector];
                    Buffer_Send_GridVel = new double[nBufferS_vector];

                    /*--- Copy the grid velocity that should be sended ---*/

                    for (iVertex = 0; iVertex < nVertexS; iVertex++) 
                    {
                        iPoint = vertex[MarkerS][iVertex]->GetNode();
                        GridVel = node[iPoint]->GetGridVel();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_GridVel[iDim*nVertexS + iVertex] = GridVel[iDim];
                    }

#ifdef HAVE_MPI
                    /*--- Send/Receive information using Sendrecv ---*/
                    MPI_Sendrecv(Buffer_Send_GridVel, nBufferS_std::vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_GridVel, nBufferR_std::vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else

                    /*--- Receive information without MPI ---*/
                    for (iVertex = 0; iVertex < nVertexR; iVertex++) 
                    {
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Receive_GridVel[iDim*nVertexR + iVertex] = Buffer_Send_GridVel[iDim*nVertexR + iVertex];
                    }

#endif

                    /*--- Deallocate send buffer ---*/
                    delete[] Buffer_Send_GridVel;

                    /*--- Do the coordinate transformation ---*/
                    for (iVertex = 0; iVertex < nVertexR; iVertex++) 
                    {

                        /*--- Find point and its type of transformation ---*/
                        iPoint = vertex[MarkerR][iVertex]->GetNode();
                        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();

                        /*--- Retrieve the supplied periodic information. ---*/

                        angles = config->GetPeriodicRotation(iPeriodic_Index);

                        /*--- Store angles separately for clarity. ---*/
                        theta = angles[0];   phi = angles[1];     psi = angles[2];
                        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
                        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

                        /*--- Compute the rotation matrix. Note that the implicit
                        ordering is rotation about the x-axis, y-axis,
                        then z-axis. Note that this is the transpose of the matrix
                        used during the preprocessing stage. ---*/

                        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

                        /*--- Copy grid velocity before performing transformation. ---*/

                        for (iDim = 0; iDim < nDim; iDim++)
                            newGridVel[iDim] = Buffer_Receive_GridVel[iDim*nVertexR + iVertex];

                        if (nDim == 2) 
                        {
                            newGridVel[0] = (rotMatrix[0][0] * Buffer_Receive_GridVel[0 * nVertexR + iVertex] +
                                rotMatrix[0][1] * Buffer_Receive_GridVel[1 * nVertexR + iVertex]);
                            newGridVel[1] = (rotMatrix[1][0] * Buffer_Receive_GridVel[0 * nVertexR + iVertex] +
                                rotMatrix[1][1] * Buffer_Receive_GridVel[1 * nVertexR + iVertex]);
                        }
                        else 
                        {
                            newGridVel[0] = (rotMatrix[0][0] * Buffer_Receive_GridVel[0 * nVertexR + iVertex] +
                                rotMatrix[0][1] * Buffer_Receive_GridVel[1 * nVertexR + iVertex] +
                                rotMatrix[0][2] * Buffer_Receive_GridVel[2 * nVertexR + iVertex]);
                            newGridVel[1] = (rotMatrix[1][0] * Buffer_Receive_GridVel[0 * nVertexR + iVertex] +
                                rotMatrix[1][1] * Buffer_Receive_GridVel[1 * nVertexR + iVertex] +
                                rotMatrix[1][2] * Buffer_Receive_GridVel[2 * nVertexR + iVertex]);
                            newGridVel[2] = (rotMatrix[2][0] * Buffer_Receive_GridVel[0 * nVertexR + iVertex] +
                                rotMatrix[2][1] * Buffer_Receive_GridVel[1 * nVertexR + iVertex] +
                                rotMatrix[2][2] * Buffer_Receive_GridVel[2 * nVertexR + iVertex]);
                        }

                        /*--- Copy transformed grid velocity back into buffer. ---*/
                        for (iDim = 0; iDim < nDim; iDim++)
                            node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);

                    }

                    /*--- Deallocate receive buffer ---*/

                    delete[] Buffer_Receive_GridVel;
                }
            }
            delete[] newGridVel;
        }

        void GEOM_GeometryPhysical::SetPeriodicBoundary(TBOX::TBOX_Config *config) 
        {
            unsigned short iMarker, jMarker, kMarker = 0, iPeriodic, iDim, nPeriodic = 0, VTK_Type;
            unsigned long iNode, iIndex, iVertex, iPoint, iElem, kElem;
            unsigned long jElem, kPoint = 0, jVertex = 0, jPoint = 0, pPoint = 0, nPointPeriodic, newNodes[4] = { 0, 0, 0, 0 };
            std::vector<unsigned long>::iterator IterElem, IterPoint[TBOX::MAX_NUMBER_PERIODIC][2];
            double *center, *angles, rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
                translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
                dx, dy, dz, rotCoord[3], epsilon = 1e-10, mindist = 1e6, *Coord_i, *Coord_j, dist = 0.0;
            bool isBadMatch = false;

            /*--- Check this dimensionalization ---*/
            std::vector<unsigned long> OldBoundaryElems[100];
            std::vector<unsigned long>::iterator IterNewElem[100];

            /*--- It only create the mirror structure for the second boundary ---*/
            bool CreateMirror[10];
            CreateMirror[1] = false;
            CreateMirror[2] = true;

            /*--- Send an initial message to the console. ---*/
            std::cout << "Setting the periodic boundary conditions." << std::endl;

            /*--- Loop through each marker to find any periodic boundaries. ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY) 
                {
                    /*--- Evaluate the number of periodic boundary conditions defined
                    in the geometry file ---*/
                    nPeriodic++;

                    /*--- Get marker index of the periodic donor boundary. ---*/
                    jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_TagBound(iMarker));

                    /*--- Write some info to the console. ---*/
                    std::cout << "Checking " << config->GetMarker_All_TagBound(iMarker);
                    std::cout << " boundary against periodic donor, " << config->GetMarker_All_TagBound(jMarker) << ". ";

                    /*--- Retrieve the supplied periodic information. ---*/
                    center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
                    angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
                    trans = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));

                    /*--- Store (center+trans) as it is constant and will be added on. ---*/
                    translation[0] = center[0] + trans[0];
                    translation[1] = center[1] + trans[1];
                    translation[2] = center[2] + trans[2];

                    /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
                    theta = angles[0];
                    phi = angles[1];
                    psi = angles[2];

                    cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
                    sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);

                    /*--- Compute the rotation matrix. Note that the implicit
                    ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
                    rotMatrix[0][0] = cosPhi*cosPsi;
                    rotMatrix[1][0] = cosPhi*sinPsi;
                    rotMatrix[2][0] = -sinPhi;

                    rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
                    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
                    rotMatrix[2][1] = sinTheta*cosPhi;

                    rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                    rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                    rotMatrix[2][2] = cosTheta*cosPhi;

                    /*--- Loop through all vertices and find/set the periodic point. ---*/
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                    {

                        /*--- Retrieve node information for this boundary point. ---*/
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        Coord_i = node[iPoint]->GetCoord();

                        /*--- Get the position std::vector from rot center to point. ---*/
                        dx = Coord_i[0] - center[0];
                        dy = Coord_i[1] - center[1];
                        if (nDim == 3) 
                        {
                            dz = Coord_i[2] - center[2];
                        }
                        else {
                            dz = 0.0;
                        }

                        /*--- Compute transformed point coordinates. ---*/
                        rotCoord[0] = rotMatrix[0][0] * dx
                            + rotMatrix[0][1] * dy
                            + rotMatrix[0][2] * dz + translation[0];

                        rotCoord[1] = rotMatrix[1][0] * dx
                            + rotMatrix[1][1] * dy
                            + rotMatrix[1][2] * dz + translation[1];

                        rotCoord[2] = rotMatrix[2][0] * dx
                            + rotMatrix[2][1] * dy
                            + rotMatrix[2][2] * dz + translation[2];

                        /*--- Perform a search to find the closest donor point. ---*/
                        mindist = 1e10;
                        for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) 
                        {
                            /*--- Retrieve information for this jPoint. ---*/
                            jPoint = vertex[jMarker][jVertex]->GetNode();
                            Coord_j = node[jPoint]->GetCoord();

                            /*--- Check the distance between the computed periodic
                            location and this jPoint. ---*/
                            dist = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++) 
                            {
                                dist += (Coord_j[iDim] - rotCoord[iDim])*(Coord_j[iDim] - rotCoord[iDim]);
                            }
                            dist = sqrt(dist);

                            /*---  Store vertex information if this is the closest
                            point found thus far. ---*/
                            if (dist < mindist) 
                            { 
                                mindist = dist; 
                                pPoint = jPoint; }
                        }

                        /*--- Set the periodic point for this iPoint. ---*/
                        vertex[iMarker][iVertex]->SetDonorPoint(pPoint, TBOX::MASTER_NODE);

                        /*--- Print warning if the nearest point was not within
                        the specified tolerance. Computation will continue. ---*/
                        if (mindist > epsilon) 
                        {
                            isBadMatch = true;
                            std::cout.precision(10);
                            std::cout << std::endl;
                            std::cout << "   Bad match for point " << iPoint << ".\tNearest";
                            std::cout << " donor distance: " << std::scientific << mindist << ".";
                        }
                    }

                    /*--- Print final warning when finding bad matches. ---*/
                    if (isBadMatch) 
                    {
                        std::cout << std::endl;
                        std::cout << "\n !!! Warning !!!" << std::endl;
                        std::cout << "Bad matches found. Computation will continue, but be cautious.\n";
                    }
                    std::cout << std::endl;
                    isBadMatch = false;
                }

            /*--- Create a std::vector to identify the points that belong to each periodic boundary condition ---*/
            bool *PeriodicBC = new bool[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++) PeriodicBC[iPoint] = false;

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY)
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        PeriodicBC[iPoint] = true;
                    }

            /*--- Determine the new points that must be added to each periodic boundary,
            note that only one of the boundaries require the extra data ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY) 
                {
                    iPeriodic = config->GetMarker_All_PerBound(iMarker);
                    /*--- An integer identify the periodic boundary condition --*/
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                    {
                        /*--- iPoint is the original point on the surface and jPoint is the
                        equivalent point in the other periodic surface ---*/
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        jPoint = vertex[iMarker][iVertex]->GetDonorPoint();

                        /*--- First the case in which it is necessary to create a mirror set of elements ---*/
                        if (CreateMirror[iPeriodic]) 
                        {
                            /*--- Now we must determine the neighbor points (including indirect ones) to the periodic points
                            and store all the information (in this list we do not include the points
                            that already belong to the periodic boundary), we also add the elements that
                            share a point with the periodic boundary condition ---*/
                            for (iIndex = 0; iIndex < node[jPoint]->GetnElem(); iIndex++) 
                            {
                                iElem = node[jPoint]->GetElem(iIndex);
                                PeriodicElem[iPeriodic].push_back(iElem);
                                for (unsigned short iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) 
                                {
                                    kPoint = elem[iElem]->GetNode(iNode);
                                    if (!PeriodicBC[kPoint]) PeriodicPoint[iPeriodic][0].push_back(kPoint);
                                }
                            }
                        }
                        /*--- Second the case where no new element is added, neither points ---*/
                        else 
                        {
                            PeriodicPoint[iPeriodic][0].push_back(jPoint);
                            PeriodicPoint[iPeriodic][1].push_back(iPoint);
                        }
                    }
                }
            }

            /*--- Sort the points that must be sended and delete repeated points ---*/
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (CreateMirror[iPeriodic]) 
                {
                    sort(PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
                    IterPoint[iPeriodic][0] = unique(PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
                    PeriodicPoint[iPeriodic][0].resize(IterPoint[iPeriodic][0] - PeriodicPoint[iPeriodic][0].begin());
                }
            }

            /*--- Create a list of the points that receive the values (only the new points) ---*/
            nPointPeriodic = nPoint;
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (CreateMirror[iPeriodic]) 
                {
                    for (iPoint = 0; iPoint < PeriodicPoint[iPeriodic][0].size(); iPoint++) 
                    {
                        PeriodicPoint[iPeriodic][1].push_back(nPointPeriodic);
                        nPointPeriodic++;
                    }
                }
            }

            /*--- Sort the elements that must be replicated in the periodic boundary
            and delete the repeated elements ---*/
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (CreateMirror[iPeriodic]) 
                {
                    sort(PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
                    IterElem = unique(PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
                    PeriodicElem[iPeriodic].resize(IterElem - PeriodicElem[iPeriodic].begin());
                }
            }

            /*--- Check all SEND points to see if they also lie on another boundary. ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                {
                    /*--- iPoint is a node that lies on the current marker. ---*/
                    iPoint = vertex[iMarker][iVertex]->GetNode();

                    /*--- Search through SEND points to check for iPoint. ---*/
                    for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
                    {
                        if (CreateMirror[iPeriodic]) 
                        {
                            /*--- jPoint is the SEND point. ---*/
                            for (iElem = 0; iElem < PeriodicPoint[iPeriodic][0].size(); iElem++)
                            {
                                jPoint = PeriodicPoint[iPeriodic][0][iElem];

                                /*--- If the two match, then jPoint lies on this boundary.
                                However, we are concerned with the new points, so we
                                will store kPoint instead. ---*/
                                if (iPoint == jPoint) 
                                {
                                    //              kPoint = PeriodicPoint[iPeriodic][1][iElem];

                                    /*--- We also want the type of boundary element that this point
                                    was within, so that we know what type of element to add
                                    built from the new points. ---*/
                                    bool isJPoint, isPeriodic;
                                    for (jElem = 0; jElem < nElem_Bound[iMarker]; jElem++) 
                                    {
                                        isJPoint = false; isPeriodic = false;
                                        for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) 
                                        {
                                            if (bound[iMarker][jElem]->GetNode(iNode) == jPoint) isJPoint = true;
                                            if (PeriodicBC[bound[iMarker][jElem]->GetNode(iNode)]) isPeriodic = true;
                                        }

                                        /*--- If both points were found, store this element. ---*/
                                        if (isJPoint && isPeriodic) 
                                        {
                                            OldBoundaryElems[iMarker].push_back(jElem);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /*--- Sort the elements that must be added and remove duplicates. ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                sort(OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
                IterNewElem[iMarker] = unique(OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
                OldBoundaryElems[iMarker].resize(IterNewElem[iMarker] - OldBoundaryElems[iMarker].begin());
            }

            /*--- Create the new boundary elements. Points making up these new
            elements must either be SEND/RECEIVE or periodic points. ---*/
            nNewElem_Bound = new unsigned long[nMarker];
            newBound = new GRID::GRID_Primal**[nMarker];
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {

                nNewElem_Bound[iMarker] = OldBoundaryElems[iMarker].size();
                newBound[iMarker] = new GRID::GRID_Primal*[nNewElem_Bound[iMarker]];

                /*--- Loop through all new elements to be added. ---*/
                for (iElem = 0; iElem < nNewElem_Bound[iMarker]; iElem++) 
                {
                    jElem = OldBoundaryElems[iMarker][iElem];

                    /*--- Loop through all nodes of this element. ---*/
                    for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) 
                    {
                        pPoint = bound[iMarker][jElem]->GetNode(iNode);

                        /*--- Check if this node is a send point. If so, the corresponding
                        receive point will be used in making the new boundary element. ---*/
                        for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
                        {
                            for (kElem = 0; kElem < PeriodicPoint[iPeriodic][0].size(); kElem++) 
                            {
                                if (pPoint == PeriodicPoint[iPeriodic][0][kElem]) newNodes[iNode] = PeriodicPoint[iPeriodic][1][kElem];
                            }
                        }

                        /*--- Check if this node is a periodic point. If so, the corresponding
                        periodic point will be used in making the new boundary element. ---*/
                        if (PeriodicBC[pPoint]) 
                        {
                            /*--- Find the corresponding periodic point. ---*/
                            for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++) 
                            {
                                if (config->GetMarker_All_KindBC(jMarker) == TBOX::PERIODIC_BOUNDARY) 
                                {
                                    for (iVertex = 0; iVertex < nVertex[jMarker]; iVertex++) 
                                    {
                                        if (pPoint == vertex[jMarker][iVertex]->GetNode()) 
                                        { 
                                            kMarker = jMarker; 
                                            jVertex = iVertex; 
                                        }
                                    }
                                }
                            }
                            newNodes[iNode] = vertex[kMarker][jVertex]->GetDonorPoint();
                        }
                    }

                    /*--- Now instantiate the new element. ---*/
                    VTK_Type = bound[iMarker][jElem]->GetVTK_Type();
                    switch (VTK_Type) 
                    {
                    case TBOX::LINE:
                        newBound[iMarker][iElem] = new GRID::GRID_Line(newNodes[0], newNodes[1], 2);
                        break;
                    case TBOX::TRIANGLE:
                        newBound[iMarker][iElem] = new GRID::GRID_Triangle(newNodes[0], newNodes[1], newNodes[2], 3);
                        break;
                    case TBOX::RECTANGLE:
                        newBound[iMarker][iElem] = new GRID::GRID_Rectangle(newNodes[0], newNodes[1], newNodes[2], newNodes[3], 3);
                        break;
                    }
                }
            }

            delete[] PeriodicBC;
        }

        void GEOM_GeometryPhysical::FindNormal_Neighbor(TBOX::TBOX_Config *config) 
        {
            double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord, *Normal;
            unsigned long Point_Normal, jPoint;
            unsigned short iNeigh, iMarker, iDim;
            unsigned long iPoint, iVertex;

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) 
            {
                if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE &&
                    config->GetMarker_All_KindBC(iMarker) != TBOX::INTERFACE_BOUNDARY &&
                    config->GetMarker_All_KindBC(iMarker) != TBOX::NEARFIELD_BOUNDARY) 
                {

                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        Normal = vertex[iMarker][iVertex]->GetNormal();

                        /*--- Compute closest normal neighbor, note that the normal are oriented inwards ---*/
                        Point_Normal = 0; cos_max = -1.0;
                        for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++)
                        {
                            jPoint = node[iPoint]->GetPoint(iNeigh);
                            scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++) 
                            {
                                diff_coord = node[jPoint]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                                scalar_prod += diff_coord*Normal[iDim];
                                norm_vect += diff_coord*diff_coord;
                                norm_Normal += Normal[iDim] * Normal[iDim];
                            }
                            norm_vect = sqrt(norm_vect);
                            norm_Normal = sqrt(norm_Normal);
                            cos_alpha = scalar_prod / (norm_vect*norm_Normal);

                            /*--- Get maximum cosine ---*/
                            if (cos_alpha >= cos_max) 
                            {
                                Point_Normal = jPoint;
                                cos_max = cos_alpha;
                            }
                        }
                        vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
                    }
                }
            }
        }

        void GEOM_GeometryPhysical::SetGeometryPlanes(TBOX::TBOX_Config *config) 
        {

            bool loop_on;
            unsigned short iMarker = 0;
            double auxXCoord, auxYCoord, auxZCoord, *Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
            unsigned long jVertex, iVertex, ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;

            /*--- Compute the total number of points on the near-field ---*/
            nVertex_Wall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if ((config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::EULER_WALL))
                    nVertex_Wall += nVertex[iMarker];

            /*--- Create an array with all the coordinates, points, pressures, face area,
            equivalent area, and nearfield weight ---*/
            Xcoord = new double[nVertex_Wall];
            Ycoord = new double[nVertex_Wall];
            if (nDim == 3)	
                Zcoord = new double[nVertex_Wall];
            FaceArea = new double[nVertex_Wall];

            /*--- Copy the boundary information to an array ---*/
            iVertex_Wall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if ((config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::HEAT_FLUX_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_CATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::ISOTHERMAL_NONCATALYTIC) ||
                    (config->GetMarker_All_KindBC(iMarker) == TBOX::EULER_WALL))
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
                        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
                        if (nDim == 3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
                        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
                        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim - 1]);
                        iVertex_Wall++;
                    }


            //std::vector<double> XCoordList;
            std::vector<double>::iterator IterXCoordList;

            for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
                XCoordList.push_back(Xcoord[iVertex]);

            sort(XCoordList.begin(), XCoordList.end());
            IterXCoordList = unique(XCoordList.begin(), XCoordList.end());
            XCoordList.resize(IterXCoordList - XCoordList.begin());

            /*--- Create std::vectors and distribute the values among the different PhiAngle queues ---*/
            Xcoord_plane.resize(XCoordList.size());
            Ycoord_plane.resize(XCoordList.size());
            if (nDim == 3) Zcoord_plane.resize(XCoordList.size());
            FaceArea_plane.resize(XCoordList.size());
            Plane_points.resize(XCoordList.size());


            double dist_ratio;
            unsigned long iCoord;

            /*--- Distribute the values among the different PhiAngles ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                if (node[iPoint]->GetDomain()) 
                {
                    loop_on = true;
                    for (ixCoord = 0; ixCoord < XCoordList.size() - 1 && loop_on; ixCoord++) 
                    {
                        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord]) / (XCoordList[ixCoord + 1] - XCoordList[ixCoord]);
                        if (dist_ratio >= 0 && dist_ratio <= 1.0) 
                        {
                            if (dist_ratio <= 0.5) iCoord = ixCoord;
                            else iCoord = ixCoord + 1;
                            Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0));
                            Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1));
                            if (nDim == 3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2));
                            FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
                            Plane_points[iCoord].push_back(iPoint);
                            loop_on = false;
                        }
                    }
                }
            }

            unsigned long auxPoint;
            /*--- Order the arrays in ascending values of y ---*/
            for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
                for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
                    for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
                        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex + 1])
                        {
                            auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex + 1]; Xcoord_plane[ixCoord][jVertex + 1] = auxXCoord;
                            auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex + 1]; Ycoord_plane[ixCoord][jVertex + 1] = auxYCoord;
                            auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex + 1]; Plane_points[ixCoord][jVertex + 1] = auxPoint;
                            if (nDim == 3) 
                            {
                                auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex + 1]; Zcoord_plane[ixCoord][jVertex + 1] = auxZCoord;
                            }
                            auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex + 1]; FaceArea_plane[ixCoord][jVertex + 1] = auxArea;
                        }

            /*--- Delete structures ---*/
            delete[] Xcoord; delete[] Ycoord;
            if (nDim == 3) delete[] Zcoord;
            delete[] FaceArea;
        }

        void GEOM_GeometryPhysical::SetBoundSensitivity(TBOX::TBOX_Config *config) 
        {
            unsigned short iMarker, icommas;
            unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
            double Sensitivity;
            bool *PointInDomain;

#ifdef HAVE_MPI
            int rank = MASTER_NODE;
            int size = SINGLE_NODE;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

            nPointLocal = nPoint;
#ifdef HAVE_MPI
            MPI_Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
            nPointGlobal = nPointLocal;
#endif

            Point2Vertex = new unsigned long[nPointGlobal][2];
            PointInDomain = new bool[nPointGlobal];

            for (iPoint = 0; iPoint < nPointGlobal; iPoint++)
                PointInDomain[iPoint] = false;

            for (iMarker = 0; iMarker < nMarker; iMarker++)
                if (config->GetMarker_All_DV(iMarker) == TBOX::YES)
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) 
                    {

                        /*--- The sensitivity file uses the global numbering ---*/
#ifndef HAVE_MPI
                        iPoint = vertex[iMarker][iVertex]->GetNode();
#else
                        iPoint = node[vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
                        if (vertex[iMarker][iVertex]->GetNode() < GetnPointDomain())
                        {
                            Point2Vertex[iPoint][0] = iMarker;
                            Point2Vertex[iPoint][1] = iVertex;
                            PointInDomain[iPoint] = true;
                            vertex[iMarker][iVertex]->SetAuxVar(0.0);
                        }
                    }

            /*--- Time-average any unsteady surface sensitivities ---*/

            unsigned long iExtIter, nExtIter;
            double delta_T, total_T;
            if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) 
            {
                nExtIter = config->GetUnst_AdjointIter();
                delta_T = config->GetDelta_UnstTimeND();
                total_T = (double)nExtIter*delta_T;
            }
            else if (config->GetUnsteady_Simulation() == TBOX::TIME_SPECTRAL) 
            {

                /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
                double period = config->GetTimeSpectral_Period();
                nExtIter = config->GetnTimeInstances();
                delta_T = period / (double)nExtIter;
                total_T = period;
            }
            else 
            {
                nExtIter = 1;
                delta_T = 1.0;
                total_T = 1.0;
            }

            for (iExtIter = 0; iExtIter < nExtIter; iExtIter++) 
            {
                /*--- Prepare to read surface sensitivity files (CSV) ---*/
                std::string text_line;
                std::ifstream Surface_file;
                char buffer[50];
                char cstr[TBOX::MAX_STRING_SIZE];
                std::string surfadj_filename = config->GetSurfAdjCoeff_FileName();
                strcpy(cstr, surfadj_filename.c_str());

                /*--- Write file name with extension if unsteady or steady ---*/

                if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
                    (config->GetUnsteady_Simulation() == TBOX::TIME_SPECTRAL)) 
                {
                    if ((int(iExtIter) >= 0) && (int(iExtIter) < 10))    sprintf(buffer, "_0000%d.csv", int(iExtIter));
                    if ((int(iExtIter) >= 10) && (int(iExtIter) < 100))   sprintf(buffer, "_000%d.csv", int(iExtIter));
                    if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000))  sprintf(buffer, "_00%d.csv", int(iExtIter));
                    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf(buffer, "_0%d.csv", int(iExtIter));
                    if (int(iExtIter) >= 10000) sprintf(buffer, "_%d.csv", int(iExtIter));
                }
                else
                    sprintf(buffer, ".csv");

                strcat(cstr, buffer);

                /*--- Read the sensitivity file ---*/

                std::string::size_type position;

                Surface_file.open(cstr, std::ios::in);
                getline(Surface_file, text_line);

                while (getline(Surface_file, text_line)) 
                {
                    for (icommas = 0; icommas < 50; icommas++) 
                    {
                        position = text_line.find(",", 0);
                        if (position != std::string::npos) text_line.erase(position, 1);
                    }
                    std::stringstream  point_line(text_line);
                    point_line >> iPoint >> Sensitivity;

                    if (PointInDomain[iPoint]) 
                    {
                        /*--- Find the vertex for the Point and Marker ---*/
                        iMarker = Point2Vertex[iPoint][0];
                        iVertex = Point2Vertex[iPoint][1];

                        /*--- Increment the auxiliary variable with the contribution of
                        this unsteady timestep. For steady problems, this reduces to
                        a single sensitivity value multiplied by 1.0. ---*/
                        vertex[iMarker][iVertex]->AddAuxVar(Sensitivity*(delta_T / total_T));
                    }
                }
                Surface_file.close();
            }
            delete[] Point2Vertex;
        }

 
    }
}
