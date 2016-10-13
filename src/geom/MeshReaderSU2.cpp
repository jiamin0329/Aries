/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for mesh data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "MeshReaderSU2.hpp"

namespace ARIES
{
  
    bool MeshReaderSU2::ReadMesh(IProcData* procData)
    {
        std::string text_line, Marker_Tag;
        std::ifstream mesh_file;
        unsigned short nMarker_Max = config->GetnMarker_Max();
        unsigned long VTK_Type, iMarker, iChar;
        unsigned long iCount = 0;
        unsigned long iElem_Bound = 0, iPoint = 0, ielem_div = 0, ielem = 0;
        unsigned long vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4];
        unsigned long vnodes_tetra[4], vnodes_hexa[8], vnodes_prism[6], vnodes_pyramid[5], dummyLong, GlobalIndex;
        char cstr[200];
        double Coord_2D[2], Coord_3D[3];
        std::string::size_type position;
        int rank = TBOX::MASTER_NODE, size = TBOX::SINGLE_NODE;
        bool domain_flag = false;
        bool found_transform = false;
        nZone = val_nZone;

        /*--- Initialize some additional counters for the parallel partitioning ---*/
        unsigned long total_pt_accounted = 0;
        unsigned long rem_points = 0;
        unsigned long element_count = 0;
        unsigned long node_count = 0;
        unsigned long loc_element_count = 0;
        bool elem_reqd = false;

        /*--- Initialize counters for local/global points & elements ---*/
#ifdef HAVE_MPI
        unsigned long LocalIndex;
        unsigned long Local_nPoint, Local_nPointDomain;
        unsigned long Local_nElem;
        unsigned long Local_nElemTri, Local_nElemQuad, Local_nElemTet;
        unsigned long Local_nElemHex, Local_nElemPrism, Local_nElemPyramid;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
        Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
        nelem_edge = 0; Global_nelem_edge = 0;
        nelem_triangle = 0; Global_nelem_triangle = 0;
        nelem_quad = 0; Global_nelem_quad = 0;
        nelem_tetra = 0; Global_nelem_tetra = 0;
        nelem_hexa = 0; Global_nelem_hexa = 0;
        nelem_prism = 0; Global_nelem_prism = 0;
        nelem_pyramid = 0; Global_nelem_pyramid = 0;

        /*--- Allocate memory for the linear partition of the mesh. These
          arrays are the size of the number of ranks. ---*/

        starting_node = new unsigned long[size];
        ending_node = new unsigned long[size];
        npoint_procs = new unsigned long[size];

        /*--- Open grid file ---*/

        strcpy(cstr, val_mesh_filename.c_str());
        mesh_file.open(cstr, std::ios::in);

        /*--- Check the grid ---*/
        if (mesh_file.fail())
        {
            std::cout << "There is no mesh file (GEOM_GeometryPhysical)!! " << cstr << std::endl;
#ifndef HAVE_MPI
            exit(EXIT_FAILURE);
#else
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
#endif
        }

        /*--- Read grid file with format SU2 ---*/
        while (getline(mesh_file, text_line))
        {
            /*--- Read the dimension of the problem ---*/
            position = text_line.find("NDIME=", 0);
            if (position != std::string::npos)
            {
                if (domain_flag == false)
                {
                    text_line.erase(0, 6); nDim = atoi(text_line.c_str());
                    if (rank == TBOX::MASTER_NODE)
                    {
                        if (nDim == 2) std::cout << "Two dimensional problem." << std::endl;
                        if (nDim == 3) std::cout << "Three dimensional problem." << std::endl;
                    }
                    domain_flag = true;
                }
                else
                {
                    break;
                }
            }

            /*--- Read number of points ---*/
            position = text_line.find("NPOIN=", 0);
            if (position != std::string::npos)
            {
                text_line.erase(0, 6);

                /*--- Check for ghost points. ---*/
                std::stringstream test_line(text_line);
                while (test_line >> dummyLong)
                    iCount++;

                /*--- Now read and store the number of points and possible ghost points. ---*/

                std::stringstream  stream_line(text_line);
                if (iCount == 2)
                {
                    stream_line >> nPoint;
                    stream_line >> nPointDomain;

                    /*--- Set some important point information for parallel simulations. ---*/
                    Global_nPoint = nPoint;
                    Global_nPointDomain = nPointDomain;
                    if (rank == TBOX::MASTER_NODE && size > TBOX::SINGLE_NODE)
                    {
                        std::cout << Global_nPointDomain << " points and " << Global_nPoint - Global_nPointDomain;
                        std::cout << " ghost points before parallel partitioning." << std::endl;
                    }
                    else if (rank == TBOX::MASTER_NODE)
                    {
                        std::cout << Global_nPointDomain << " points and " << Global_nPoint - Global_nPointDomain;
                        std::cout << " ghost points." << std::endl;
                    }

                }
                else if (iCount == 1)
                {
                    stream_line >> nPoint;
                    nPointDomain = nPoint;
                    Global_nPointDomain = nPoint;
                    Global_nPoint = nPoint;
                    if (rank == TBOX::MASTER_NODE && size > TBOX::SINGLE_NODE)
                    {
                        std::cout << nPoint << " points before parallel partitioning." << std::endl;
                    }
                    else if (rank == TBOX::MASTER_NODE)
                    {
                        std::cout << nPoint << " points." << std::endl;
                    }
                }
                else
                {
                    std::cout << "NPOIN improperly specified!!" << std::endl;
#ifndef HAVE_MPI
                    exit(EXIT_FAILURE);
#else
                    MPI_Abort(MPI_COMM_WORLD, 1);
                    MPI_Finalize();
#endif
                }

                if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
                    std::cout << "Performing linear partitioning of the grid nodes." << std::endl;

                /*--- Compute the number of points that will be on each processor.
                  This is a linear partitioning with the addition of a simple load
                  balancing for any remainder points. ---*/

                total_pt_accounted = 0;
                for (unsigned long i = 0; i < size; i++)
                {
                    npoint_procs[i] = nPoint / size;
                    total_pt_accounted = total_pt_accounted + npoint_procs[i];
                }

                /*--- Get the number of remainder points after the even division ---*/
                rem_points = nPoint - total_pt_accounted;
                for (unsigned long i = 0; i < rem_points; i++)
                {
                    npoint_procs[i]++;
                }

                /*--- Store the local number of nodes and the beginning/end index ---*/
                local_node = npoint_procs[rank];
                starting_node[0] = 0;
                ending_node[0] = starting_node[0] + npoint_procs[0];
                for (unsigned long i = 1; i < size; i++)
                {
                    starting_node[i] = ending_node[i - 1];
                    ending_node[i] = starting_node[i] + npoint_procs[i];
                }

                /*--- Here we check if a point in the mesh file lies in the domain
                  and if so then store it on the local processor. We only create enough
                  space in the node container for the local nodes at this point. ---*/

                node = new GRID::GRID_DGPoint*[local_node];
                iPoint = 0; node_count = 0;
                while (node_count < nPoint)
                {
                    getline(mesh_file, text_line);
                    std::istringstream point_line(text_line);

                    /*--- We only read information for this node if it is owned by this
                      rank based upon our initial linear partitioning. ---*/

                    if ((node_count >= starting_node[rank]) && (node_count < ending_node[rank]))
                    {
                        switch (nDim)
                        {
                            case 2:
                                GlobalIndex = node_count;
#ifndef HAVE_MPI
                                point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
#else
                                if (size > SINGLE_NODE) { point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; point_line >> LocalIndex; point_line >> GlobalIndex; }
                                else { point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; LocalIndex = iPoint; GlobalIndex = node_count; }
#endif
                                node[iPoint] = new GRID::GRID_DGPoint(Coord_2D[0], Coord_2D[1], GlobalIndex, config);
                                iPoint++; break;
                            case 3:
                                GlobalIndex = node_count;
#ifndef HAVE_MPI
                                point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
#else
                                if (size > SINGLE_NODE) { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; point_line >> LocalIndex; point_line >> GlobalIndex; }
                                else { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; LocalIndex = iPoint; GlobalIndex = node_count; }
#endif
                                node[iPoint] = new GRID::GRID_DGPoint(Coord_3D[0], Coord_3D[1], Coord_3D[2], GlobalIndex, config);
                                iPoint++; break;
                        }
                    }
                    node_count++;
                }
            }
        }

        mesh_file.close();
        strcpy(cstr, val_mesh_filename.c_str());

        /*--- Initialize some arrays for the adjacency information (ParMETIS). ---*/

        //  unsigned long *adj_counter = new unsigned long[local_node];
        //  unsigned long **adjacent_elem = new unsigned long*[local_node];
        adj_counter = new unsigned long[local_node];
        adjacent_elem = new unsigned long*[local_node];

        for (iPoint = 0; iPoint < local_node; iPoint++)
        {
            //adjacent_elem[iPoint] = new unsigned long[2000]; modified by jiamin, 12-Nov-2015
            adjacent_elem[iPoint] = new unsigned long[1000];
            adj_counter[iPoint] = 0;
        }

        mesh_file.open(cstr, std::ios::in);
        while (getline(mesh_file, text_line))
        {
            /*--- Read the information about inner elements ---*/
            position = text_line.find("NELEM=", 0);
            if (position != std::string::npos)
            {
                text_line.erase(0, 6); nElem = atoi(text_line.c_str());

                /*--- Store total number of elements in the original mesh ---*/
                Global_nElem = nElem;
                if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
                    std::cout << Global_nElem << " interior elements before parallel partitioning." << std::endl;

                /*--- Allocate space for elements ---*/

                elem = new GRID::GRID_Primal*[nElem];
                for (int iElem = 0; iElem < nElem; iElem++) elem[iElem] = NULL;


                /*--- Set up the global to local element mapping. ---*/
                Global_to_local_elem = new long[nElem];
                for (unsigned long i = 0; i < nElem; i++)
                {
                    Global_to_local_elem[i] = -1;
                }

                if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
                    std::cout << "Distributing elements across all ranks." << std::endl;

                /*--- Loop over all the volumetric elements and store any element that
                  contains at least one of an owned node for this rank (i.e., there will
                  be element redundancy, since multiple ranks will store the same elems
                  on the boundaries of the initial linear partitioning. ---*/

                // TO DO: remove redundant edges (quads have extra diagonals for instance)
                element_count = 0; elem_reqd = false; loc_element_count = 0; ielem_div = 0;
                while (ielem_div < nElem)
                {
                    getline(mesh_file, text_line);
                    std::istringstream elem_line(text_line);

                    elem_line >> VTK_Type;
                    elem_reqd = false;

                    /*--- Decide whether this rank needs each element. If so, build the
                      adjacency arrays needed by ParMETIS and store the element connectivity.
                      Note that every proc starts it's node indexing from zero. ---*/

                    switch (VTK_Type)
                    {

                        case TBOX::TRIANGLE:
                            elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
                            for (unsigned long i = 0; i < TBOX::N_POINTS_TRIANGLE; i++)
                            {
                                if ((vnodes_triangle[i] >= starting_node[rank]) && (vnodes_triangle[i] < ending_node[rank]))
                                {
                                    elem_reqd = true;
                                    for (unsigned long j = 0; j < TBOX::N_POINTS_TRIANGLE; j++)
                                    {
                                        if (i != j)
                                        {
                                            adjacent_elem[vnodes_triangle[i] - starting_node[rank]][adj_counter[vnodes_triangle[i] - starting_node[rank]]] = vnodes_triangle[j];
                                            adj_counter[vnodes_triangle[i] - starting_node[rank]]++;
                                        }
                                    }
                                }
                            }
                            if (elem_reqd)
                            {
                                Global_to_local_elem[element_count] = loc_element_count;
                                elem[loc_element_count] = new GRID::GRID_Triangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2], 2);
                                nelem_triangle++; loc_element_count++;
                            }
                            break;

                        case TBOX::RECTANGLE:
                            elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
                            for (unsigned long i = 0; i < TBOX::N_POINTS_QUADRILATERAL; i++)
                            {
                                if ((vnodes_quad[i] >= starting_node[rank]) && (vnodes_quad[i] < ending_node[rank]))
                                {
                                    elem_reqd = true;

                                    for (unsigned long j = 0; j < TBOX::N_POINTS_QUADRILATERAL; j++)
                                    {
                                        if (i != j)
                                        {
                                            adjacent_elem[vnodes_quad[i] - starting_node[rank]][adj_counter[vnodes_quad[i] - starting_node[rank]]] = vnodes_quad[j];
                                            adj_counter[vnodes_quad[i] - starting_node[rank]]++;
                                        }
                                    }
                                }
                            }
                            if (elem_reqd)
                            {
                                Global_to_local_elem[element_count] = loc_element_count;
                                elem[loc_element_count] = new GRID::GRID_Rectangle(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3], 2);
                                loc_element_count++; nelem_quad++;
                            }
                            break;

                        case TBOX::TETRAHEDRON:
                            elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
                            for (unsigned long i = 0; i < TBOX::N_POINTS_TETRAHEDRON; i++)
                            {
                                if ((vnodes_tetra[i] >= starting_node[rank]) && (vnodes_tetra[i] < ending_node[rank]))
                                {
                                    elem_reqd = true;
                                    for (unsigned long j = 0; j < TBOX::N_POINTS_TETRAHEDRON; j++)
                                    {
                                        if (i != j)
                                        {
                                            adjacent_elem[vnodes_tetra[i] - starting_node[rank]][adj_counter[vnodes_tetra[i] - starting_node[rank]]] = vnodes_tetra[j];
                                            adj_counter[vnodes_tetra[i] - starting_node[rank]]++;
                                        }
                                    }
                                }
                            }
                            if (elem_reqd)
                            {
                                Global_to_local_elem[element_count] = loc_element_count;
                                elem[loc_element_count] = new GRID::GRID_Tetrahedron(vnodes_tetra[0], vnodes_tetra[1], vnodes_tetra[2], vnodes_tetra[3]);
                                loc_element_count++; nelem_tetra++;
                            }
                            break;

                        case TBOX::HEXAHEDRON:

                            elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
                            elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
                            elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
                            for (unsigned long i = 0; i < TBOX::N_POINTS_HEXAHEDRON; i++)
                            {
                                if ((vnodes_hexa[i] >= starting_node[rank]) && (vnodes_hexa[i] < ending_node[rank]))
                                {
                                    elem_reqd = true;
                                    for (unsigned long j = 0; j < TBOX::N_POINTS_HEXAHEDRON; j++)
                                    {
                                        if (i != j)
                                        {
                                            adjacent_elem[vnodes_hexa[i] - starting_node[rank]][adj_counter[vnodes_hexa[i] - starting_node[rank]]] = vnodes_hexa[j];
                                            adj_counter[vnodes_hexa[i] - starting_node[rank]]++;
                                        }
                                    }
                                }
                            }
                            if (elem_reqd)
                            {
                                Global_to_local_elem[element_count] = loc_element_count;
                                elem[loc_element_count] = new GRID::GRID_Hexahedron(vnodes_hexa[0], vnodes_hexa[1], vnodes_hexa[2], vnodes_hexa[3],
                                                                                    vnodes_hexa[4], vnodes_hexa[5], vnodes_hexa[6], vnodes_hexa[7]);
                                loc_element_count++; nelem_hexa++;
                            }
                            break;

                        case TBOX::PRISM:

                            elem_line >> vnodes_prism[0]; elem_line >> vnodes_prism[1]; elem_line >> vnodes_prism[2];
                            elem_line >> vnodes_prism[3]; elem_line >> vnodes_prism[4]; elem_line >> vnodes_prism[5];
                            for (unsigned long i = 0; i < TBOX::N_POINTS_PRISM; i++) {
                                if ((vnodes_prism[i] >= starting_node[rank]) && (vnodes_prism[i] < ending_node[rank]))
                                {
                                    elem_reqd = true;
                                    for (unsigned long j = 0; j < TBOX::N_POINTS_PRISM; j++)
                                    {
                                        if (i != j)
                                        {
                                            adjacent_elem[vnodes_prism[i] - starting_node[rank]][adj_counter[vnodes_prism[i] - starting_node[rank]]] = vnodes_prism[j];
                                            adj_counter[vnodes_prism[i] - starting_node[rank]]++;
                                        }
                                    }
                                }
                            }
                            if (elem_reqd)
                            {
                                Global_to_local_elem[element_count] = loc_element_count;
                                elem[loc_element_count] = new GRID::GRID_Prism(vnodes_prism[0], vnodes_prism[1], vnodes_prism[2], vnodes_prism[3], vnodes_prism[4], vnodes_prism[5]);
                                loc_element_count++; nelem_prism++;
                            }
                            break;

                        case TBOX::PYRAMID:
                            elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
                            elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
                            for (unsigned long i = 0; i < TBOX::N_POINTS_PYRAMID; i++)
                            {
                                if ((vnodes_pyramid[i] >= starting_node[rank]) && (vnodes_pyramid[i] < ending_node[rank]))
                                {
                                    elem_reqd = true;
                                    for (unsigned long j = 0; j < TBOX::N_POINTS_PYRAMID; j++)
                                    {
                                        if (i != j)
                                        {
                                            adjacent_elem[vnodes_pyramid[i] - starting_node[rank]][adj_counter[vnodes_pyramid[i] - starting_node[rank]]] = vnodes_pyramid[j];
                                            adj_counter[vnodes_pyramid[i] - starting_node[rank]]++;
                                        }
                                    }
                                }
                            }
                            if (elem_reqd)
                            {
                                Global_to_local_elem[element_count] = loc_element_count;
                                elem[loc_element_count] = new GRID::GRID_Pyramid(vnodes_pyramid[0], vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                                loc_element_count++; nelem_pyramid++;
                            }
                            break;
                    }
                    ielem_div++;
                    element_count++;
                }
            }
        }

        mesh_file.close();


        if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
            std::cout << "Calling the partitioning functions." << std::endl;

        /*--- Store the number of local elements on each rank after determining
          which elements must be kept in the loop above. ---*/

        no_of_local_elements = loc_element_count;

        /*--- Post process the adjacency information in order to get it into the
          proper format before sending the data to ParMETIS. We need to remove
          repeats and adjust the size of the array for each local node. ---*/
        if ((rank == TBOX::MASTER_NODE) && (size > TBOX::SINGLE_NODE))
            std::cout << "Building the graph adjacency structure." << std::endl;

        unsigned long loc_adjc_size = 0;
        std::vector<unsigned long> adjac_vec;
        unsigned long adj_elem_size;
        std::vector<unsigned long>::iterator it;
        local_elem = loc_element_count;

        xadj = new unsigned long[npoint_procs[rank] + 1];
        xadj[0] = 0;
        std::vector<unsigned long> temp_adjacency;
        unsigned long local_count = 0;

        for (unsigned long i = 0; i < local_node; i++)
        {
            for (unsigned long j = 0; j < adj_counter[i]; j++)
            {
                temp_adjacency.push_back(adjacent_elem[i][j]);
            }

            sort(temp_adjacency.begin(), temp_adjacency.end());
            it = unique(temp_adjacency.begin(), temp_adjacency.end());
            loc_adjc_size = it - temp_adjacency.begin();

            temp_adjacency.resize(loc_adjc_size);
            xadj[local_count + 1] = xadj[local_count] + loc_adjc_size;
            local_count++;

            for (unsigned long j = 0; j < loc_adjc_size; j++)
            {
                adjac_vec.push_back(temp_adjacency[j]);
            }
            temp_adjacency.clear();
        }

        /*--- Now that we know the size, create the final adjacency array ---*/
        adj_elem_size = xadj[npoint_procs[rank]];
        adjacency = new unsigned long[adj_elem_size];
        copy(adjac_vec.begin(), adjac_vec.end(), adjacency);

        xadj_size = npoint_procs[rank] + 1;
        adjacency_size = adj_elem_size;

        /*--- Free temporary memory used to build the adjacency. ---*/
        adjac_vec.clear();
        delete[] adj_counter;
        for (iPoint = 0; iPoint < local_node; iPoint++)
        {
            delete[] adjacent_elem[iPoint];
        }
        delete[] adjacent_elem;

        /*--- For now, the boundary marker information is still read by the
          master node alone (and eventually distributed by the master as well).
          In the future, this component will also be performed in parallel. ---*/

        mesh_file.open(cstr, std::ios::in);

        if (rank == TBOX::MASTER_NODE)
        {

            while (getline(mesh_file, text_line))
            {
                /*--- Read number of markers ---*/
                position = text_line.find("NMARK=", 0);
                if (position != std::string::npos)
                {
                    text_line.erase(0, 6); nMarker = atoi(text_line.c_str());
                    if (rank == TBOX::MASTER_NODE)
                        std::cout << nMarker << " surface markers." << std::endl;
                    config->SetnMarker_All(nMarker);
                    bound = new GRID::GRID_Primal**[nMarker];
                    nElem_Bound = new unsigned long[nMarker];
                    Tag_to_Marker = new std::string[nMarker_Max];

                    for (iMarker = 0; iMarker < nMarker; iMarker++)
                    {
                        getline(mesh_file, text_line);
                        text_line.erase(0, 11);
                        std::string::size_type position;
                        for (iChar = 0; iChar < 20; iChar++)
                        {
                            position = text_line.find(" ", 0);
                            if (position != std::string::npos) text_line.erase(position, 1);
                            position = text_line.find("\r", 0);
                            if (position != std::string::npos) text_line.erase(position, 1);
                            position = text_line.find("\n", 0);
                            if (position != std::string::npos) text_line.erase(position, 1);
                        }
                        Marker_Tag = text_line.c_str();

                        /*--- Physical boundaries definition ---*/
                        if (Marker_Tag != "SEND_RECEIVE")
                        {
                            getline(mesh_file, text_line);
                            text_line.erase(0, 13); nElem_Bound[iMarker] = atoi(text_line.c_str());
                            if (rank == TBOX::MASTER_NODE)
                                std::cout << nElem_Bound[iMarker] << " boundary elements in index " << iMarker << " (Marker = " << Marker_Tag << ")." << std::endl;

                            /*--- Allocate space for elements ---*/
                            bound[iMarker] = new GRID::GRID_Primal*[nElem_Bound[iMarker]];

                            nelem_edge_bound = 0; nelem_triangle_bound = 0; nelem_quad_bound = 0; ielem = 0;
                            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                            {
                                getline(mesh_file, text_line);
                                std::istringstream bound_line(text_line);
                                bound_line >> VTK_Type;
                                switch (VTK_Type)
                                {
                                    case TBOX::LINE:

                                        if (nDim == 3)
                                        {
                                            std::cout << "Please remove line boundary conditions from the mesh file!" << std::endl;
#ifndef HAVE_MPI
                                            exit(EXIT_FAILURE);
#else
                                            MPI_Abort(MPI_COMM_WORLD, 1);
                                            MPI_Finalize();
#endif
                                        }

                                        bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                                        bound[iMarker][ielem] = new GRID::GRID_Line(vnodes_edge[0], vnodes_edge[1], 2);
                                        ielem++; nelem_edge_bound++; break;

                                    case TBOX::TRIANGLE:
                                        bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                                        bound[iMarker][ielem] = new GRID::GRID_Triangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2], 3);
                                        ielem++; nelem_triangle_bound++; break;

                                    case TBOX::RECTANGLE:

                                        bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];

                                        bound[iMarker][ielem] = new GRID::GRID_Rectangle(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3], 3);
                                        ielem++; nelem_quad_bound++;

                                        break;
                                }
                            }

                            /*--- Update config information storing the boundary information in the right place ---*/
                            Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
                            config->SetMarker_All_TagBound(iMarker, Marker_Tag);
                            config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
                            config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
                            config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
                            config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
                            config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
                            config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
                            config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
                            config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
                            config->SetMarker_All_SendRecv(iMarker, TBOX::NONE);
                            config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
                        }

                        /*--- Send-Receive boundaries definition ---*/
                        else
                        {
                            unsigned long nelem_vertex = 0, vnodes_vertex;
                            unsigned short transform;
                            getline(mesh_file, text_line);
                            text_line.erase(0, 13); nElem_Bound[iMarker] = atoi(text_line.c_str());
                            bound[iMarker] = new GRID::GRID_Primal*[nElem_Bound[iMarker]];

                            nelem_vertex = 0; ielem = 0;
                            getline(mesh_file, text_line); text_line.erase(0, 8);
                            config->SetMarker_All_KindBC(iMarker, TBOX::SEND_RECEIVE);
                            config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));

                            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                            {
                                getline(mesh_file, text_line);
                                std::istringstream bound_line(text_line);
                                bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;

                                bound[iMarker][ielem] = new GRID::GRID_VertexMPI(vnodes_vertex, nDim);
                                bound[iMarker][ielem]->SetRotation_Type(transform);
                                ielem++; nelem_vertex++;
                            }
                        }
                    }
                }

                /*--- Read periodic transformation info (center, rotation, translation) ---*/
                position = text_line.find("NPERIODIC=", 0);
                if (position != std::string::npos)
                {
                    unsigned short nPeriodic, iPeriodic, iIndex;

                    /*--- Set bool signifying that periodic transormations were found ---*/
                    found_transform = true;

                    /*--- Read and store the number of transformations. ---*/
                    text_line.erase(0, 10); nPeriodic = atoi(text_line.c_str());
                    if (rank == TBOX::MASTER_NODE)
                    {
                        if (nPeriodic - 1 != 0)
                            std::cout << nPeriodic - 1 << " periodic transformations." << std::endl;
                    }
                    config->SetnPeriodicIndex(nPeriodic);

                    /*--- Store center, rotation, & translation in that order for each. ---*/
                    for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++)
                    {
                        getline(mesh_file, text_line);
                        position = text_line.find("PERIODIC_INDEX=", 0);
                        if (position != std::string::npos)
                        {
                            text_line.erase(0, 15); iIndex = atoi(text_line.c_str());
                            if (iIndex != iPeriodic)
                            {
                                std::cout << "PERIODIC_INDEX out of order in SU2 file!!" << std::endl;
#ifndef HAVE_MPI
                                exit(EXIT_FAILURE);
#else
                                MPI_Abort(MPI_COMM_WORLD, 1);
                                MPI_Finalize();
#endif
                            }
                        }
                        double* center = new double[3];
                        double* rotation = new double[3];
                        double* translate = new double[3];
                        getline(mesh_file, text_line);
                        std::istringstream cent(text_line);
                        cent >> center[0]; cent >> center[1]; cent >> center[2];
                        config->SetPeriodicCenter(iPeriodic, center);
                        getline(mesh_file, text_line);
                        std::istringstream rot(text_line);
                        rot >> rotation[0]; rot >> rotation[1]; rot >> rotation[2];
                        config->SetPeriodicRotation(iPeriodic, rotation);
                        getline(mesh_file, text_line);
                        std::istringstream tran(text_line);
                        tran >> translate[0]; tran >> translate[1]; tran >> translate[2];
                        config->SetPeriodicTranslate(iPeriodic, translate);
                    }
                }
            }

            /*--- If no periodic transormations were found, store default zeros ---*/
            if (!found_transform)
            {
                unsigned short nPeriodic = 1, iPeriodic = 0;
                config->SetnPeriodicIndex(nPeriodic);
                double* center = new double[3];
                double* rotation = new double[3];
                double* translate = new double[3];
                for (unsigned short iDim = 0; iDim < 3; iDim++)
                {
                    center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
                }
                config->SetPeriodicCenter(iPeriodic, center);
                config->SetPeriodicRotation(iPeriodic, rotation);
                config->SetPeriodicTranslate(iPeriodic, translate);
            }
        }

        /*--- Close the input file ---*/
        mesh_file.close();
    }

        
    
}
