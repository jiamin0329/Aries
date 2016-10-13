


#include "GEOM_GeometryPeriodic.hpp"

namespace ARIES
{
    namespace GEOM
    {

        GEOM_GeometryPeriodic::GEOM_GeometryPeriodic(GEOM_Geometry *geometry, TBOX::TBOX_Config *config)
        {
            unsigned long nElem_new, nPoint_new, jPoint, iPoint, iElem, jElem, iVertex,
                nelem_triangle = 0, nelem_quad = 0, nelem_tetra = 0, nelem_hexa = 0, nelem_prism = 0,
                nelem_pyramid = 0, iIndex, newElementsBound = 0;
            unsigned short  iMarker, nPeriodic = 0, iPeriodic;
            double *center, *angles, rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
                translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
                dx, dy, dz, rotCoord[3], *Coord_i;
            unsigned short nMarker_Max = config->GetnMarker_Max();

            /*--- It only create the mirror structure for the second boundary ---*/
            bool CreateMirror[10];
            CreateMirror[1] = false;
            CreateMirror[2] = true;

            /*--- Compute the number of periodic bc on the geometry ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY)
                    nPeriodic++;

            /*--- Write the number of dimensions of the problem ---*/
            nDim = geometry->GetnDim();

            /*--- Copy the new boundary element information from the geometry class.
            Be careful, as these are pointers to vectors/objects. ---*/
            nNewElem_BoundPer = geometry->nNewElem_Bound;
            newBoundPer = geometry->newBound;

            /*--- Count the number of new boundary elements. ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                newElementsBound += nNewElem_BoundPer[iMarker];

            /*--- Loop over the original grid to perform the dimensionalizaton of the new vectors ---*/
            nElem_new = 0; nPoint_new = 0;
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
            {
                if (CreateMirror[iPeriodic]) 
                {
                    nElem_new += geometry->PeriodicElem[iPeriodic].size();
                    nPoint_new += geometry->PeriodicPoint[iPeriodic][0].size();
                }
            }

            std::cout << "Number of new points: " << nPoint_new << "." << std::endl;
            std::cout << "Number of new interior elements: " << nElem_new << "." << std::endl;
            std::cout << "Number of new boundary elements added to preexisting markers: " << newElementsBound << "." << std::endl;

            /*--- Create a copy of the original grid ---*/
            elem = new GRID::GRID_Primal*[geometry->GetnElem() + nElem_new];
            for (iElem = 0; iElem < geometry->GetnElem(); iElem++) 
            {
                switch (geometry->elem[iElem]->GetVTK_Type()) 
                {
                case TBOX::TRIANGLE:
                    elem[iElem] = new GRID::GRID_Triangle(geometry->elem[iElem]->GetNode(0),
                        geometry->elem[iElem]->GetNode(1),
                        geometry->elem[iElem]->GetNode(2), 2);
                    nelem_triangle++;
                    break;

                case TBOX::RECTANGLE:
                    elem[iElem] = new GRID::GRID_Rectangle(geometry->elem[iElem]->GetNode(0),
                        geometry->elem[iElem]->GetNode(1),
                        geometry->elem[iElem]->GetNode(2),
                        geometry->elem[iElem]->GetNode(3), 2);
                    nelem_quad++;
                    break;

                case TBOX::TETRAHEDRON:
                    elem[iElem] = new GRID::GRID_Tetrahedron(geometry->elem[iElem]->GetNode(0),
                        geometry->elem[iElem]->GetNode(1),
                        geometry->elem[iElem]->GetNode(2),
                        geometry->elem[iElem]->GetNode(3));
                    nelem_tetra++;
                    break;

                case TBOX::HEXAHEDRON:
                    elem[iElem] = new GRID::GRID_Hexahedron(geometry->elem[iElem]->GetNode(0),
                        geometry->elem[iElem]->GetNode(1),
                        geometry->elem[iElem]->GetNode(2),
                        geometry->elem[iElem]->GetNode(3),
                        geometry->elem[iElem]->GetNode(4),
                        geometry->elem[iElem]->GetNode(5),
                        geometry->elem[iElem]->GetNode(6),
                        geometry->elem[iElem]->GetNode(7));
                    nelem_hexa++;
                    break;

                case TBOX::PRISM:
                    elem[iElem] = new GRID::GRID_Prism(geometry->elem[iElem]->GetNode(0),
                        geometry->elem[iElem]->GetNode(1),
                        geometry->elem[iElem]->GetNode(2),
                        geometry->elem[iElem]->GetNode(3),
                        geometry->elem[iElem]->GetNode(4),
                        geometry->elem[iElem]->GetNode(5));
                    nelem_prism++;
                    break;

                case TBOX::PYRAMID:
                    elem[iElem] = new GRID::GRID_Pyramid(geometry->elem[iElem]->GetNode(0),
                        geometry->elem[iElem]->GetNode(1),
                        geometry->elem[iElem]->GetNode(2),
                        geometry->elem[iElem]->GetNode(3),
                        geometry->elem[iElem]->GetNode(4));
                    nelem_pyramid++;
                    break;
                }
            }

            /*--- Create a list with all the points and the new index ---*/
            unsigned long *Index = new unsigned long[geometry->GetnPoint()];
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Index[iPoint] = 0;

            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (CreateMirror[iPeriodic]) 
                {
                    for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) 
                    {
                        iPoint = geometry->PeriodicPoint[iPeriodic][0][iIndex];
                        Index[iPoint] = geometry->PeriodicPoint[iPeriodic][1][iIndex];
                    }
                }
            }

            for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY)
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) 
                    {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
                        Index[iPoint] = jPoint;
                    }

            /*--- Add the new elements due to the periodic boundary condtion ---*/
            iElem = geometry->GetnElem();

            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (CreateMirror[iPeriodic]) 
                {
                    for (iIndex = 0; iIndex < geometry->PeriodicElem[iPeriodic].size(); iIndex++) 
                    {
                        jElem = geometry->PeriodicElem[iPeriodic][iIndex];

                        switch (geometry->elem[jElem]->GetVTK_Type()) 
                        {
                        case TBOX::TRIANGLE:
                            elem[iElem] = new GRID::GRID_Triangle(Index[geometry->elem[jElem]->GetNode(0)],
                                Index[geometry->elem[jElem]->GetNode(1)],
                                Index[geometry->elem[jElem]->GetNode(2)], 2);
                            iElem++; nelem_triangle++;
                            break;

                        case TBOX::RECTANGLE:
                            elem[iElem] = new GRID::GRID_Rectangle(Index[geometry->elem[jElem]->GetNode(0)],
                                Index[geometry->elem[jElem]->GetNode(1)],
                                Index[geometry->elem[jElem]->GetNode(2)],
                                Index[geometry->elem[jElem]->GetNode(3)], 2);
                            iElem++; nelem_quad++;
                            break;

                        case TBOX::TETRAHEDRON:
                            elem[iElem] = new GRID::GRID_Tetrahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                Index[geometry->elem[jElem]->GetNode(1)],
                                Index[geometry->elem[jElem]->GetNode(2)],
                                Index[geometry->elem[jElem]->GetNode(3)]);
                            iElem++; nelem_tetra++;
                            break;

                        case TBOX::HEXAHEDRON:
                            elem[iElem] = new GRID::GRID_Hexahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                Index[geometry->elem[jElem]->GetNode(1)],
                                Index[geometry->elem[jElem]->GetNode(2)],
                                Index[geometry->elem[jElem]->GetNode(3)],
                                Index[geometry->elem[jElem]->GetNode(4)],
                                Index[geometry->elem[jElem]->GetNode(5)],
                                Index[geometry->elem[jElem]->GetNode(6)],
                                Index[geometry->elem[jElem]->GetNode(7)]);
                            iElem++; nelem_hexa++;
                            break;

                        case TBOX::PRISM:
                            elem[iElem] = new GRID::GRID_Prism(Index[geometry->elem[jElem]->GetNode(0)],
                                Index[geometry->elem[jElem]->GetNode(1)],
                                Index[geometry->elem[jElem]->GetNode(2)],
                                Index[geometry->elem[jElem]->GetNode(3)],
                                Index[geometry->elem[jElem]->GetNode(4)],
                                Index[geometry->elem[jElem]->GetNode(5)]);
                            iElem++; nelem_prism++;
                            break;

                        case TBOX::PYRAMID:
                            elem[iElem] = new GRID::GRID_Pyramid(Index[geometry->elem[jElem]->GetNode(0)],
                                Index[geometry->elem[jElem]->GetNode(1)],
                                Index[geometry->elem[jElem]->GetNode(2)],
                                Index[geometry->elem[jElem]->GetNode(3)],
                                Index[geometry->elem[jElem]->GetNode(4)]);
                            iElem++; nelem_pyramid++;
                            break;

                        }
                    }
                }
            }

            nElem = geometry->GetnElem() + nElem_new;

            /*--- Add the old points ---*/
            node = new GRID::GRID_DGPoint*[geometry->GetnPoint() + nPoint_new];
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
            {
                if (geometry->GetnDim() == 2)
                    node[iPoint] = new GRID::GRID_DGPoint(geometry->node[iPoint]->GetCoord(0),
                    geometry->node[iPoint]->GetCoord(1), iPoint, config);
                if (geometry->GetnDim() == 3)
                    node[iPoint] = new GRID::GRID_DGPoint(geometry->node[iPoint]->GetCoord(0),
                    geometry->node[iPoint]->GetCoord(1),
                    geometry->node[iPoint]->GetCoord(2), iPoint, config);
            }

            /*--- Add the new points due to the periodic boundary condtion (only in the mirror part) ---*/
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (CreateMirror[iPeriodic]) 
                {
                    for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) 
                    {
                        /*--- From iPeriodic obtain the iMarker ---*/
                        for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                            if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;

                        /*--- Retrieve the supplied periodic information. ---*/
                        center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
                        angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
                        trans = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));

                        /*--- Store center - trans as it is constant and will be added on.
                        Note the subtraction, as this is the inverse translation. ---*/
                        translation[0] = center[0] - trans[0];
                        translation[1] = center[1] - trans[1];
                        translation[2] = center[2] - trans[2];

                        /*--- Store angles separately for clarity. Compute sines/cosines.
                        Note the negative sign, as this is the inverse rotation. ---*/
                        theta = -angles[0];
                        phi = -angles[1];
                        psi = -angles[2];

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

                        /*--- Retrieve node information for this boundary point. ---*/
                        iPoint = geometry->PeriodicPoint[iPeriodic][0][iIndex];
                        jPoint = geometry->PeriodicPoint[iPeriodic][1][iIndex];
                        Coord_i = geometry->node[iPoint]->GetCoord();

                        /*--- Get the position vector from rot center to point. ---*/
                        dx = Coord_i[0] - center[0];
                        dy = Coord_i[1] - center[1];
                        if (nDim == 3) 
                        {
                            dz = Coord_i[2] - center[2];
                        }
                        else 
                        {
                            dz = 0.0;
                        }

                        /*--- Compute transformed point coordinates. ---*/
                        rotCoord[0] = rotMatrix[0][0] * dx + rotMatrix[0][1] * dy + rotMatrix[0][2] * dz + translation[0];
                        rotCoord[1] = rotMatrix[1][0] * dx + rotMatrix[1][1] * dy + rotMatrix[1][2] * dz + translation[1];
                        rotCoord[2] = rotMatrix[2][0] * dx + rotMatrix[2][1] * dy + rotMatrix[2][2] * dz + translation[2];

                        /*--- Save the new points with the new coordinates. ---*/
                        if (geometry->GetnDim() == 2)
                            node[jPoint] = new GRID::GRID_DGPoint(rotCoord[0], rotCoord[1], jPoint, config);
                        if (geometry->GetnDim() == 3)
                            node[jPoint] = new GRID::GRID_DGPoint(rotCoord[0], rotCoord[1], rotCoord[2], jPoint, config);

                    }
                }
            }

            nPoint = geometry->GetnPoint() + nPoint_new;

            /*--- Add the old boundary, reserving space for two new bc (send/recive periodic bc) ---*/
            nMarker = geometry->GetnMarker() + 2;
            nElem_Bound = new unsigned long[nMarker];
            bound = new GRID::GRID_Primal**[nMarker];
            Tag_to_Marker = new std::string[nMarker_Max];
            config->SetnMarker_All(nMarker);

            /*--- Copy the olf boundary ---*/
            for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) 
            {

                bound[iMarker] = new GRID::GRID_Primal*[geometry->GetnElem_Bound(iMarker)];

                for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) 
                {
                    if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TBOX::LINE)
                        bound[iMarker][iVertex] = new GRID::GRID_Line(geometry->bound[iMarker][iVertex]->GetNode(0),
                        geometry->bound[iMarker][iVertex]->GetNode(1), 2);
                    if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TBOX::TRIANGLE)
                        bound[iMarker][iVertex] = new GRID::GRID_Triangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                        geometry->bound[iMarker][iVertex]->GetNode(1),
                        geometry->bound[iMarker][iVertex]->GetNode(2), 3);
                    if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TBOX::RECTANGLE)
                        bound[iMarker][iVertex] = new GRID::GRID_Rectangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                        geometry->bound[iMarker][iVertex]->GetNode(1),
                        geometry->bound[iMarker][iVertex]->GetNode(2),
                        geometry->bound[iMarker][iVertex]->GetNode(3), 3);
                }

                nElem_Bound[iMarker] = geometry->GetnElem_Bound(iMarker);
                Tag_to_Marker[iMarker] = geometry->GetMarker_Tag(iMarker);

            }

            delete[] Index;

        }

        GEOM_GeometryPeriodic::~GEOM_GeometryPeriodic(void) 
        {
            unsigned long iElem_Bound;
            unsigned short iMarker;

            for (iMarker = 0; iMarker < nMarker; iMarker++) 
            {
                for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) 
                {
                    if (newBoundPer[iMarker][iElem_Bound] != NULL) delete newBoundPer[iMarker][iElem_Bound];
                }
            }
            if (newBoundPer != NULL) delete[] newBoundPer;

            if (nNewElem_BoundPer != NULL) delete[] nNewElem_BoundPer;
        }

        void GEOM_GeometryPeriodic::SetPeriodicBoundary(GEOM_Geometry *geometry, TBOX::TBOX_Config *config)
        {
            unsigned short iMarker, iPeriodic, nPeriodic = 0, iMarkerSend, iMarkerReceive;
            unsigned long iVertex, Counter_Send = 0, Counter_Receive = 0, iIndex;

            /*--- Compute the number of periodic bc on the geometry ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY)
                    nPeriodic++;

            /*--- First compute the Send/Receive boundaries, count the number of points ---*/
            Counter_Send = 0; 	Counter_Receive = 0;
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
                    Counter_Send += geometry->PeriodicPoint[iPeriodic][0].size();
                if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
                    Counter_Receive += geometry->PeriodicPoint[iPeriodic][1].size();
            }

            /*--- Adimensionalization of the new boundaries ---*/
            iMarkerSend = nMarker - 2; iMarkerReceive = nMarker - 1;
            config->SetMarker_All_SendRecv(iMarkerSend, 1);
            config->SetMarker_All_SendRecv(iMarkerReceive, -1);
            nElem_Bound[iMarkerSend] = Counter_Send;
            nElem_Bound[iMarkerReceive] = Counter_Receive;
            bound[iMarkerSend] = new GRID::GRID_Primal*[Counter_Send];
            bound[iMarkerReceive] = new GRID::GRID_Primal*[Counter_Receive];

            /*--- First we do the send ---*/
            iVertex = 0;
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
                if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
                    for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) 
                    {
                        bound[iMarkerSend][iVertex] = new GRID::GRID_VertexMPI(geometry->PeriodicPoint[iPeriodic][0][iIndex], nDim);
                        bound[iMarkerSend][iVertex]->SetRotation_Type(iPeriodic);
                        iVertex++;
                    }

            /*--- Second we do the receive ---*/
            iVertex = 0;
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
                if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
                    for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][1].size(); iIndex++) 
                    {
                        bound[iMarkerReceive][iVertex] = new GRID::GRID_VertexMPI(geometry->PeriodicPoint[iPeriodic][1][iIndex], nDim);
                        bound[iMarkerReceive][iVertex]->SetRotation_Type(iPeriodic);
                        iVertex++;
                    }

        }

        void GEOM_GeometryPeriodic::SetMeshFile(GEOM_Geometry *geometry, TBOX::TBOX_Config *config, std::string val_mesh_out_filename) 
        {
            unsigned long iElem, iPoint, iElem_Bound, GhostPoints;
            unsigned short iMarker, iNodes, iDim;
            unsigned short iMarkerReceive, iPeriodic, nPeriodic = 0;
            std::ofstream output_file;
            std::string Grid_Marker;
            char *cstr;
            double *center, *angles, *transl;

            cstr = new char[val_mesh_out_filename.size() + 1];
            strcpy(cstr, val_mesh_out_filename.c_str());

            /*--- Open .su2 grid file ---*/
            output_file.precision(15);
            output_file.open(cstr, std::ios::out);

            /*--- Ghost points, look at the nodes in the send receive ---*/
            iMarkerReceive = nMarker - 1;
            GhostPoints = nElem_Bound[iMarkerReceive];

            /*--- Change the numbering to guarantee that the all the receive
            points are at the end of the file ---*/
            unsigned long OldnPoint = geometry->GetnPoint();
            //unsigned long NewSort[nPoint];
            unsigned long *NewSort = new unsigned long[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                NewSort[iPoint] = iPoint;
            }

            unsigned long Index = OldnPoint - 1;
            for (iMarker = 0; iMarker < nMarker; iMarker++) 
            {
                if (bound[iMarker][0]->GetVTK_Type() == TBOX::VERTEX) 
                {
                    if (config->GetMarker_All_SendRecv(iMarker) < 0) 
                    {
                        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) 
                        {
                            if (bound[iMarker][iElem_Bound]->GetNode(0) < geometry->GetnPoint()) 
                            {
                                NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] = Index;
                                NewSort[Index] = bound[iMarker][iElem_Bound]->GetNode(0);
                                Index--;
                            }
                        }
                    }
                }
            }


            /*--- Write dimension, number of elements and number of points ---*/
            output_file << "NDIME= " << nDim << std::endl;
            output_file << "NELEM= " << nElem << std::endl;
            for (iElem = 0; iElem < nElem; iElem++) 
            {
                output_file << elem[iElem]->GetVTK_Type();
                for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
                    output_file << "\t" << NewSort[elem[iElem]->GetNode(iNodes)];
                output_file << "\t" << iElem << std::endl;
            }

            output_file << "NPOIN= " << nPoint << "\t" << nPoint - GhostPoints << std::endl;
            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                for (iDim = 0; iDim < nDim; iDim++)
                    output_file << std::scientific << "\t" << node[NewSort[iPoint]]->GetCoord(iDim);
                output_file << "\t" << iPoint << std::endl;
            }

            output_file << "NMARK= " << nMarker << std::endl;
            for (iMarker = 0; iMarker < nMarker; iMarker++) 
            {
                if (bound[iMarker][0]->GetVTK_Type() != TBOX::VERTEX) 
                {
                    Grid_Marker = config->GetMarker_All_TagBound(iMarker);
                    output_file << "MARKER_TAG= " << Grid_Marker << std::endl;
                    output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] + nNewElem_BoundPer[iMarker] << std::endl;

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) 
                    {
                        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t";
                        for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes() - 1; iNodes++)
                            output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t";
                        iNodes = bound[iMarker][iElem_Bound]->GetnNodes() - 1;
                        output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << std::endl;
                    }

                    /*--- Write any new elements at the end of the list. ---*/
                    if (nNewElem_BoundPer[iMarker] > 0) 
                    {
                        for (iElem_Bound = 0; iElem_Bound < nNewElem_BoundPer[iMarker]; iElem_Bound++) 
                        {
                            output_file << newBoundPer[iMarker][iElem_Bound]->GetVTK_Type() << "\t";
                            for (iNodes = 0; iNodes < newBoundPer[iMarker][iElem_Bound]->GetnNodes() - 1; iNodes++)
                                output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t";
                            iNodes = newBoundPer[iMarker][iElem_Bound]->GetnNodes() - 1;
                            output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << std::endl;
                        }
                    }
                }

                if (bound[iMarker][0]->GetVTK_Type() == TBOX::VERTEX) 
                {
                    output_file << "MARKER_TAG= SEND_RECEIVE" << std::endl;
                    output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] << std::endl;
                    if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << std::endl;
                    if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << std::endl;

                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
                            NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] << "\t" <<
                            bound[iMarker][iElem_Bound]->GetRotation_Type() << std::endl;
                    }
                }
            }

            /*--- Compute the number of periodic bc on the geometry ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_KindBC(iMarker) == TBOX::PERIODIC_BOUNDARY)
                    nPeriodic++;

            output_file << "NPERIODIC= " << nPeriodic + 1 << std::endl;

            /*--- Periodic 0 correspond with no movement of the surface ---*/
            output_file << "PERIODIC_INDEX= 0" << std::endl;
            output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << std::endl;
            output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << std::endl;
            output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << std::endl;

            /*--- From iPeriodic obtain the iMarker ---*/
            for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) 
            {
                for (iMarker = 0; iMarker < nMarker; iMarker++)
                    if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;

                /*--- Retrieve the supplied periodic information. ---*/
                center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
                angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
                transl = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));

                output_file << "PERIODIC_INDEX= " << iPeriodic << std::endl;
                output_file << center[0] << "\t" << center[1] << "\t" << center[2] << std::endl;
                output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << std::endl;
                output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << std::endl;

            }

            delete []NewSort;
            output_file.close();
        }

        void GEOM_GeometryPeriodic::SetTecPlot(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file) 
        {
            unsigned long iElem, iPoint;
            unsigned short iDim;
            std::ofstream Tecplot_File;

            Tecplot_File.open(mesh_filename, std::ios::out);
            Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << std::endl;

            if (nDim == 2) 
            {
                Tecplot_File << "VARIABLES = \"x\",\"y\" " << std::endl;
                Tecplot_File << "ZONE NODES= " << nPoint << ", ELEMENTS= " << nElem << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << std::endl;
            }
            if (nDim == 3) 
            {
                Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << std::endl;
                Tecplot_File << "ZONE NODES= " << nPoint << ", ELEMENTS= " << nElem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
            }

            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                for (iDim = 0; iDim < nDim; iDim++)
                    Tecplot_File << std::scientific << node[iPoint]->GetCoord(iDim) << "\t";
                Tecplot_File << "\n";
            }

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
    }
}