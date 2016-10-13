/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for cgns mesh data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "MeshReaderCGNS.hpp"

namespace ARIES
{
    MeshReaderCGNS::MeshReaderCGNS()
    {

    }

    MeshReaderCGNS::~MeshReaderCGNS()
    {
        
    }



    bool MeshReaderCGNS::ReadMesh(IProcData* procData)
    {
        /*--- Original CGNS reader implementation by Thomas D. Economon,
          Francisco Palacios. Improvements for mixed-element meshes generated
          by ICEM added by Martin Spel (3D) & Shlomy Shitrit (2D), April 2014.
          Parallel version by Thomas D. Economon, February 2015. ---*/
#ifdef ARIES_HAVE_CGNS
        /*--- Local variables and initialization ---*/
        std::string text_line, Marker_Tag;
        std::ifstream mesh_file;
        unsigned short VTK_Type = 0, iMarker = 0;
        unsigned short nMarker_Max = config->GetnMarker_Max();
        unsigned long iPoint = 0, iProcessor = 0, ielem = 0, GlobalIndex = 0;
        unsigned long globalOffset = 0;
        int rank = TBOX::MASTER_NODE, size = TBOX::SINGLE_NODE;
        nZone = val_nZone;

        /*--- Local variables needed when calling the CGNS mid-level API. ---*/
        unsigned long vnodes_cgns[8];
        double Coord_cgns[3];
        int fn, nbases = 0, nzones = 0, ngrids = 0, ncoords = 0, nsections = 0;
        int *vertices = NULL, *cells = NULL, nMarkers = 0, *boundVerts = NULL, npe;
        int interiorElems = 0, totalVerts = 0;
        int cell_dim = 0, phys_dim = 0, nbndry, parent_flag, file_type;
        char basename[TBOX::CGNS_STRING_SIZE], zonename[TBOX::CGNS_STRING_SIZE];
        char coordname[CGNS_STRING_SIZE];
        cgsize_t* cgsize; cgsize = new cgsize_t[3];
        ZoneType_t zonetype;
        DataType_t datatype;
        double** coordArray = NULL;
        double*** gridCoords = NULL;
        ElementType_t elemType;
            cgsize_t range_min, range_max, startE, endE;
            range_min = 1;
            std::string currentElem;
            int** elemTypeVTK = NULL;
            int** elemIndex = NULL;
            int** elemBegin = NULL;
            int** elemEnd = NULL;
            int** nElems = NULL;
            int indexMax, elemMax; indexMax = elemMax = 0;
            cgsize_t**** connElems = NULL;
            cgsize_t* connElemCGNS = NULL;
            cgsize_t* connElemTemp = NULL;
            cgsize_t ElementDataSize = 0;
            cgsize_t* parentData = NULL;
            int** dataSize = NULL;
            bool** isInternal = NULL;
            char*** sectionNames = NULL;

            /*--- Initialize counters for local/global points & elements ---*/

#ifdef ARIES_HAVE_MPI
            unsigned long Local_nElem;
            unsigned long Local_nElemTri, Local_nElemQuad, Local_nElemTet;
            unsigned long Local_nElemHex, Local_nElemPrism, Local_nElemPyramid;

            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            MPI_Request *send_req, *recv_req;
            MPI_Status  status;
            int ind;
#endif

            /*--- Initialize counters for local/global points & elements ---*/

            Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
            nelem_edge = 0; Global_nelem_edge = 0;
            nelem_triangle = 0; Global_nelem_triangle = 0;
            nelem_quad = 0; Global_nelem_quad = 0;
            nelem_tetra = 0; Global_nelem_tetra = 0;
            nelem_hexa = 0; Global_nelem_hexa = 0;
            nelem_prism = 0; Global_nelem_prism = 0;
            nelem_pyramid = 0; Global_nelem_pyramid = 0;

            /*--- Initialize some additional counters for the parallel partitioning ---*/

            unsigned long total_pt_accounted = 0;
            unsigned long rem_points = 0;
            unsigned long element_count = 0;
            unsigned long loc_element_count = 0;
            unsigned long element_remainder = 0;
            unsigned long total_elems = 0;

            /*--- Allocate memory for the linear partitioning of the mesh. These
            arrays are the size of the number of ranks. ---*/

            starting_node = new unsigned long[size];
            ending_node = new unsigned long[size];
            npoint_procs = new unsigned long[size];

            unsigned long *nPoint_Linear = new unsigned long[size + 1];
            unsigned long *nElem_Linear = new unsigned long[size];

            unsigned long *elemB = new unsigned long[size];
            unsigned long *elemE = new unsigned long[size];

            unsigned long *elemGlobalID = NULL;

            unsigned short *nPoinPerElem = NULL;
            unsigned short *elemTypes = NULL;

            bool *isMixed = NULL;

            unsigned short connSize = 10;

            /*--- Check whether the supplied file is truly a CGNS file. ---*/
            if (cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK) {
                if (rank == MASTER_NODE) {
                    printf("\n\n   !!! Error !!!\n");
                    printf(" %s is not a CGNS file.\n", val_mesh_filename.c_str());
                    printf(" Now exiting...\n\n");
                }
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD, 1);
                MPI_Finalize();
#endif
            }

            /*--- Open the CGNS file for reading. The value of fn returned
            is the specific index number for this file and will be
            repeatedly used in the function calls. ---*/

            if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn)) cg_error_exit();
            if (rank == MASTER_NODE) {
                cout << "Reading the CGNS file: ";
                cout << val_mesh_filename.c_str() << "." << endl;
            }

            /*--- Get the number of databases. This is the highest node
            in the CGNS heirarchy. ---*/

            if (cg_nbases(fn, &nbases)) cg_error_exit();
            if (rank == MASTER_NODE)
                cout << "CGNS file contains " << nbases << " database(s)." << endl;

            /*--- Check if there is more than one database. Throw an
            error if there is because this reader can currently
            only handle one database. ---*/

            if (nbases > 1) {
                if (rank == MASTER_NODE) {
                    printf("\n\n   !!! Error !!!\n");
                    printf("CGNS reader currently incapable of handling more than 1 database.");
                    printf("Now exiting...\n\n");
                }
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD, 1);
                MPI_Finalize();
#endif
            }

            /*--- Read the databases. Note that the CGNS indexing starts at 1. ---*/

            for (int i = 1; i <= nbases; i++) {

                if (cg_base_read(fn, i, basename, &cell_dim, &phys_dim)) cg_error_exit();

                /*--- Get the number of zones for this base. ---*/

                if (cg_nzones(fn, i, &nzones)) cg_error_exit();
                if (rank == MASTER_NODE) {
                    cout << "Database " << i << ", " << basename << ": " << nzones;
                    cout << " zone(s), cell dimension of " << cell_dim << ", physical ";
                    cout << "dimension of " << phys_dim << "." << endl;
                }

                /*--- Check if there is more than one zone. Throw an
                error if there is, because this reader can currently
                only handle one zone. This could be extended in the future. ---*/

                if (nzones > 1) {
                    if (rank == MASTER_NODE) {
                        printf("\n\n   !!! Error !!!\n");
                        printf("CGNS reader currently incapable of handling more than 1 zone.");
                        printf("Now exiting...\n\n");
                    }
#ifndef HAVE_MPI
                    exit(EXIT_FAILURE);
#else
                    MPI_Abort(MPI_COMM_WORLD, 1);
                    MPI_Finalize();
#endif
                }

                /*--- Initialize some data structures for  all zones. ---*/

                vertices = new int[nzones];
                cells = new int[nzones];
                boundVerts = new int[nzones];
                coordArray = new double*[nzones];
                gridCoords = new double**[nzones];
                elemTypeVTK = new int*[nzones];
                elemIndex = new int*[nzones];
                elemBegin = new int*[nzones];
                elemEnd = new int*[nzones];
                nElems = new int*[nzones];
                dataSize = new int*[nzones];
                isInternal = new bool*[nzones];
                nMarkers = 0;
                sectionNames = new char**[nzones];
                connElems = new cgsize_t***[nzones];

                /*--- Loop over all zones in this base. Again, indexing starts at 1. ---*/

                for (int j = 1; j <= nzones; j++) {

                    /*--- Read the basic information for this zone, including
                    the name and the number of vertices, cells, and
                    boundary cells which are stored in the cgsize variable. ---*/

                    if (cg_zone_read(fn, i, j, zonename, cgsize)) cg_error_exit();

                    /*--- Rename the zone size information for clarity.
                    NOTE: The number of cells here may be only the number of
                    interior elements or it may be the total. This needs to
                    be counted explicitly later. ---*/

                    vertices[j - 1] = cgsize[0];
                    cells[j - 1] = cgsize[1];
                    boundVerts[j - 1] = cgsize[2];

                    /*--- Increment the total number of vertices from all zones. ---*/

                    nPoint = vertices[j - 1];
                    nPointDomain = vertices[j - 1];

                    Global_nPoint = vertices[j - 1];
                    Global_nPointDomain = vertices[j - 1];

                    totalVerts += vertices[j - 1];

                    /*--- Print some information about the current zone. ---*/

                    if (cg_zone_type(fn, i, j, &zonetype)) cg_error_exit();
                    if (rank == MASTER_NODE) {
                        cout << "Zone " << j << ", " << zonename << ": " << vertices[j - 1];
                        cout << " vertices, " << cells[j - 1] << " cells, " << boundVerts[j - 1];
                        cout << " boundary vertices." << endl;
                    }

                    /*--- Retrieve the number of grids in this zone. For now, we know
                    this is one, but to be more general, this will need to check and
                    allow for a loop over all grids. ---*/

                    if (cg_ngrids(fn, i, j, &ngrids)) cg_error_exit();
                    if (ngrids > 1) {
                        if (rank == MASTER_NODE) {
                            printf("\n\n   !!! Error !!!\n");
                            printf("CGNS reader currently handles only 1 grid per zone.");
                            printf("Now exiting...\n\n");
                        }
#ifndef HAVE_MPI
                        exit(EXIT_FAILURE);
#else
                        MPI_Abort(MPI_COMM_WORLD, 1);
                        MPI_Finalize();
#endif
                    }

                    /*--- Check the number of coordinate arrays stored in this zone.
                    Should be 2 for 2-D grids and 3 for 3-D grids. ---*/

                    if (cg_ncoords(fn, i, j, &ncoords)) cg_error_exit();
                    if (rank == MASTER_NODE) {
                        cout << "Reading grid coordinates." << endl;
                        cout << "Number of coordinate dimensions is " << ncoords << "." << endl;
                    }

                    /*--- Compute the number of points that will be on each processor.
                    This is a linear partitioning with the addition of a simple load
                    balancing for any remainder points. ---*/

                    total_pt_accounted = 0;
                    for (int ii = 0; ii < size; ii++) {
                        npoint_procs[ii] = vertices[j - 1] / size;
                        total_pt_accounted = total_pt_accounted + npoint_procs[ii];
                    }

                    /*--- Get the number of remainder points after the even division ---*/

                    rem_points = vertices[j - 1] - total_pt_accounted;
                    for (int ii = 0; ii < rem_points; ii++) {
                        npoint_procs[ii]++;
                    }

                    /*--- Store the local number of nodes and the beginning/end index ---*/

                    local_node = npoint_procs[rank];
                    starting_node[0] = 0;
                    ending_node[0] = starting_node[0] + npoint_procs[0];
                    nPoint_Linear[0] = 0;
                    for (int ii = 1; ii < size; ii++) {
                        starting_node[ii] = ending_node[ii - 1];
                        ending_node[ii] = starting_node[ii] + npoint_procs[ii];
                        nPoint_Linear[ii] = nPoint_Linear[ii - 1] + npoint_procs[ii - 1];
                    }
                    nPoint_Linear[size] = vertices[j - 1];

                    /*--- Set the value of range_max to the total number of nodes in
                    the unstructured mesh. Also allocate memory for the temporary array
                    that will hold the grid coordinates as they are extracted. Note the
                    +1 for CGNS convention. ---*/

                    range_min = (cgsize_t)starting_node[rank] + 1;
                    range_max = (cgsize_t)ending_node[rank];
                    coordArray[j - 1] = new double[local_node];

                    /*--- Allocate memory for the 2-D array that will store the x, y,
                    & z (if required) coordinates for writing into the SU2 mesh. ---*/

                    gridCoords[j - 1] = new double*[ncoords];
                    for (int ii = 0; ii < ncoords; ii++) {
                        *(gridCoords[j - 1] + ii) = new double[local_node];
                    }

                    /*--- Loop over each set of coordinates. Note again
                    that the indexing starts at 1. ---*/

                    for (int k = 1; k <= ncoords; k++) {

                        /*--- Read the coordinate info. This will retrieve the
                        data type (either RealSingle or RealDouble) as
                        well as the coordname which will specifiy the
                        type of data that it is based in the SIDS convention.
                        This might be "CoordinateX," for instance. ---*/

                        if (cg_coord_info(fn, i, j, k, &datatype, coordname))
                            cg_error_exit();
                        if (rank == MASTER_NODE) {
                            cout << "Loading " << coordname;
                            cout << " values into linear partitions." << endl;
                        }

                        /*--- Always retrieve the grid coords in double precision. ---*/

                        if (datatype != RealDouble) {
                            printf("\n\n   !!! Error !!!\n");
                            printf(" CGNS coordinates are not double precision.\n");
                            printf(" Now exiting...\n\n");
#ifndef HAVE_MPI
                            exit(EXIT_FAILURE);
#else
                            MPI_Abort(MPI_COMM_WORLD, 1);
                            MPI_Finalize();
#endif
                        }
                        if (cg_coord_read(fn, i, j, coordname, datatype, &range_min,
                            &range_max, coordArray[j - 1])) cg_error_exit();

                        /*--- Copy these coords into the array for storage until
                        writing the SU2 mesh. ---*/

                        for (int m = 0; m < local_node; m++) {
                            gridCoords[j - 1][k - 1][m] = coordArray[j - 1][m];
                        }

                    }

                    /*--- Begin section for retrieving the connectivity info. ---*/

                    if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
                        cout << "Distributing connectivity across all ranks." << endl;

                    /*--- First check the number of sections. ---*/

                    if (cg_nsections(fn, i, j, &nsections)) cg_error_exit();
                    if (rank == MASTER_NODE) {
                        cout << "Number of connectivity sections is ";
                        cout << nsections << "." << endl;
                    }

                    /*--- Allocate several data structures to hold the various
                    pieces of information describing each section. It is
                    stored in this manner so that it can be written to
                    SU2 memory later. ---*/

                    elemTypeVTK[j - 1] = new int[nsections];
                    elemIndex[j - 1] = new int[nsections];
                    elemBegin[j - 1] = new int[nsections];
                    elemEnd[j - 1] = new int[nsections];
                    nElems[j - 1] = new int[nsections];
                    dataSize[j - 1] = new int[nsections];
                    isInternal[j - 1] = new bool[nsections];

                    sectionNames[j - 1] = new char*[nsections];
                    for (int ii = 0; ii < nsections; ii++) {
                        sectionNames[j - 1][ii] = new char[CGNS_STRING_SIZE];
                    }

                    connElems[j - 1] = new cgsize_t**[nsections];

                    /*--- Loop over each section. This will include the main
                    connectivity information for the grid cells, as well
                    as any boundaries which were labeled before export. ---*/

                    for (int s = 1; s <= nsections; s++) {

                        /*--- Read the connectivity details for this section.
                        Store the total number of elements in this section
                        to be used later for memory allocation. ---*/

                        if (cg_section_read(fn, i, j, s, sectionNames[j - 1][s - 1],
                            &elemType, &startE, &endE, &nbndry,
                            &parent_flag)) cg_error_exit();

                        /*--- Store the beginning and ending index for this section. ---*/

                        elemBegin[j - 1][s - 1] = (int)startE;
                        elemEnd[j - 1][s - 1] = (int)endE;

                        /*--- Compute element linear partitioning ---*/

                        element_count = (int)(endE - startE + 1);
                        total_elems = 0;
                        for (int ii = 0; ii < size; ii++) {
                            nElem_Linear[ii] = element_count / size;
                            total_elems += nElem_Linear[ii];
                        }

                        /*--- Get the number of remainder elements after even division ---*/

                        element_remainder = element_count - total_elems;
                        for (int ii = 0; ii < element_remainder; ii++) {
                            nElem_Linear[ii]++;
                        }

                        /*--- Store the number of elements that this rank is responsible for
                        in the current section. ---*/

                        nElems[j - 1][s - 1] = (int)nElem_Linear[rank];

                        /*--- Get starting and end element index for my rank. ---*/

                        elemB[0] = startE;
                        elemE[0] = startE + nElem_Linear[0] - 1;
                        for (unsigned long ii = 1; ii < size; ii++) {
                            elemB[ii] = elemE[ii - 1] + 1;
                            elemE[ii] = elemB[ii] + nElem_Linear[ii] - 1;
                        }

                        /*--- Allocate some memory for the handling the connectivity
                        and auxiliary data that we are need to communicate. ---*/

                        connElemCGNS = new cgsize_t[nElems[j - 1][s - 1] * connSize];
                        nPoinPerElem = new unsigned short[nElems[j - 1][s - 1]];
                        elemGlobalID = new unsigned long[nElems[j - 1][s - 1]];
                        elemTypes = new unsigned short[nElems[j - 1][s - 1]];

                        isMixed = new bool[nElems[j - 1][s - 1]];
                        for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) isMixed[ii] = false;

                        /*--- Retrieve the connectivity information and store. Note that
                        we are only accessing our rank's piece of the data here in the
                        partial read function in the CGNS API. ---*/

                        if (cg_elements_partial_read(fn, i, j, s, (cgsize_t)elemB[rank],
                            (cgsize_t)elemE[rank], connElemCGNS,
                            parentData) != CG_OK) cg_error_exit();

                        /*--- Find the number of nodes required to represent
                        this type of element. ---*/

                        ElementType_t elmt_type;
                        if (cg_npe(elemType, &npe)) cg_error_exit();

                        /*--- Loop through all of the elements in this section to get more
                        information and to decide whether it has internal elements. ---*/

                        int counter = 0;
                        for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) {

                            /*--- If we have a mixed element section, we need to check the elem
                            type one by one. Set the flag to true if mixed. ---*/

                            if (elemType == MIXED) {
                                elmt_type = ElementType_t(connElemCGNS[counter]);
                                cg_npe(elmt_type, &npe);
                                counter++; for (int jj = 0; jj < npe; jj++) counter++;
                                isMixed[ii] = true;
                            }
                            else {
                                elmt_type = elemType;
                            }

                            /*--- Store the number of verts per elem for the current elem. ---*/

                            nPoinPerElem[ii] = npe;

                            /*--- Store the global ID for this element. Note the -1 to move
                            from CGNS convention to SU2 convention. We also subtract off
                            an additional offset in case we have found boundary sections
                            prior to this one, in order to keep the internal element global
                            IDs indexed starting from zero. ---*/

                            elemGlobalID[ii] = elemB[rank] + ii - 1 - globalOffset;

                            /*--- Need to check the element type and correctly specify the
                            VTK identifier for that element. SU2 recognizes elements by
                            their VTK number. ---*/

                            switch (elmt_type) {
                            case NODE:
                                currentElem = "Vertex";
                                elemTypes[ii] = 1;
                                break;
                            case BAR_2:
                                currentElem = "Line";
                                elemTypes[ii] = 3;
                                break;
                            case BAR_3:
                                currentElem = "Line";
                                elemTypes[ii] = 3;
                                break;
                            case TRI_3:
                                currentElem = "Triangle";
                                elemTypes[ii] = 5;
                                break;
                            case QUAD_4:
                                currentElem = "Quadrilateral";
                                elemTypes[ii] = 9;
                                break;
                            case TETRA_4:
                                currentElem = "Tetrahedron";
                                elemTypes[ii] = 10;
                                break;
                            case HEXA_8:
                                currentElem = "Hexahedron";
                                elemTypes[ii] = 12;
                                break;
                            case PENTA_6:
                                currentElem = "Prism";
                                elemTypes[ii] = 13;
                                break;
                            case PYRA_5:
                                currentElem = "Pyramid";
                                elemTypes[ii] = 14;
                                break;
                            case HEXA_20:
                                if (rank == MASTER_NODE) {
                                    printf("\n\n   !!! Error !!!\n");
                                    printf(" HEXA-20 element type not supported\n");
                                    printf(" Section %d, npe=%d\n", s, npe);
                                    printf(" startE %d, endE %d\n", startE, endE);
                                    printf(" Now exiting...\n\n");
                                }
#ifndef HAVE_MPI
                                exit(EXIT_FAILURE);
#else
                                MPI_Abort(MPI_COMM_WORLD, 1);
                                MPI_Finalize();
#endif
                                break;
                            default:
                                if (rank == MASTER_NODE) {
                                    printf("\n\n   !!! Error !!!\n");
                                    printf(" Unknown elem: (type %d, npe=%d)\n", elemType, npe);
                                    printf(" Section %d\n", s);
                                    printf(" startE %d, endE %d\n", startE, endE);
                                    printf(" Now exiting...\n\n");
                                }
#ifndef HAVE_MPI
                                exit(EXIT_FAILURE);
#else
                                MPI_Abort(MPI_COMM_WORLD, 1);
                                MPI_Finalize();
#endif
                                break;
                            }

                            /*--- Check if the elements in this section are part
                            of the internal domain or are part of the boundary
                            surfaces. This will be used to separate the
                            internal connectivity from the boundary connectivity.
                            We will check for quad and tri elements for 3-D meshes
                            because these will be the boundaries. Similarly, line
                            elements will be boundaries to 2-D problems. ---*/

                            if (cell_dim == 2) {

                                /*--- In 2-D check for line elements, VTK type 3. ---*/

                                if (elemTypes[ii] == 3) {
                                    isInternal[j - 1][s - 1] = false;
                                }
                                else {
                                    isInternal[j - 1][s - 1] = true;
                                    interiorElems++;
                                }

                            }
                            else if (cell_dim == 3) {

                                /*--- In 3-D check for tri/quad elements, VTK types 5 or 9. ---*/

                                switch (elemTypes[ii]) {
                                case 5:
                                case 9:
                                    isInternal[j - 1][s - 1] = false;
                                    break;
                                default:
                                    isInternal[j - 1][s - 1] = true;
                                    interiorElems++;
                                    break;
                                }

                            }
                        }

                        /*--- Print some information to the console. ---*/

                        if (rank == MASTER_NODE) {
                            for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++)
                                if (isMixed[ii]) { currentElem = "Mixed"; break; }
                            cout << "Loading section " << sectionNames[j - 1][s - 1];
                            cout << " of element type " << currentElem << "." << endl;
                        }

                        /*--- If we have found that this is a boundary section (we assume
                        that internal cells and boundary cells do not exist in the same
                        section together), the master node read the boundary section.
                        Otherwise, we have all ranks read and communicate the internals. ---*/

                        if (!isInternal[j - 1][s - 1]) {

                            /*--- Master node should read this entire marker section. Free
                            the memory for the conn. from the CGNS file since we are going
                            to read the section again with the master. ---*/

                            delete[] connElemCGNS;
                            delete[] nPoinPerElem;
                            delete[] elemTypes;
                            delete[] elemGlobalID;
                            delete[] isMixed;

                            /*--- Since we found an internal section, we should adjust the
                            element global ID offset by the total size of the section. ---*/

                            globalOffset += element_count;

                            if (rank == MASTER_NODE) {

                                /*--- First increment the markers ---*/

                                nMarkers++;

                                /*--- Read the section info again ---*/

                                if (cg_section_read(fn, i, j, s, sectionNames[j - 1][s - 1],
                                    &elemType, &startE, &endE, &nbndry,
                                    &parent_flag)) cg_error_exit();

                                /*--- Store the number of elems (all on the master). ---*/

                                nElems[j - 1][s - 1] = (int)(endE - startE + 1);

                                /*--- Read and store the total amount of data that will be
                                listed when reading this section. ---*/

                                if (cg_ElementDataSize(fn, i, j, s, &ElementDataSize))
                                    cg_error_exit();
                                dataSize[j - 1][s - 1] = ElementDataSize;

                                /*--- Find the number of nodes required to represent
                                this type of element. ---*/

                                if (cg_npe(elemType, &npe)) cg_error_exit();
                                elemIndex[j - 1][s - 1] = npe;

                                /*--- Need to check the element type and correctly
                                specify the VTK identifier for that element.
                                SU2 recognizes elements by their VTK number. ---*/

                                switch (elemType) {
                                case NODE:
                                    elemTypeVTK[j - 1][s - 1] = 1;
                                    break;
                                case BAR_2:
                                    elemTypeVTK[j - 1][s - 1] = 3;
                                    break;
                                case BAR_3:
                                    elemTypeVTK[j - 1][s - 1] = 3;
                                    break;
                                case TRI_3:
                                    elemTypeVTK[j - 1][s - 1] = 5;
                                    break;
                                case QUAD_4:
                                    elemTypeVTK[j - 1][s - 1] = 9;
                                    break;
                                case TETRA_4:
                                    elemTypeVTK[j - 1][s - 1] = 10;
                                    break;
                                case HEXA_8:
                                    elemTypeVTK[j - 1][s - 1] = 12;
                                    break;
                                case PENTA_6:
                                    elemTypeVTK[j - 1][s - 1] = 13;
                                    break;
                                case PYRA_5:
                                    elemTypeVTK[j - 1][s - 1] = 14;
                                    break;
                                case HEXA_20:
                                    printf("\n\n   !!! Error !!!\n");
                                    printf(" HEXA-20 element type not supported\n");
                                    printf(" Section %d, npe=%d\n", s, npe);
                                    printf(" startE %d, endE %d\n", startE, endE);
                                    printf(" Now exiting...\n\n");
#ifndef HAVE_MPI
                                    exit(EXIT_FAILURE);
#else
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                    MPI_Finalize();
#endif
                                    break;
                                case MIXED:
                                    currentElem = "Mixed";
                                    elemTypeVTK[j - 1][s - 1] = -1;
                                    break;
                                default:
                                    printf("\n\n   !!! Error !!!\n");
                                    printf(" Unknown elem: (type %d, npe=%d)\n", elemType, npe);
                                    printf(" Section %d\n", s);
                                    printf(" startE %d, endE %d\n", startE, endE);
                                    printf(" Now exiting...\n\n");
#ifndef HAVE_MPI
                                    exit(EXIT_FAILURE);
#else
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                    MPI_Finalize();
#endif
                                    break;
                                }

                                /*--- In case of mixed data type, allocate place for 8 nodes
                                maximum (hex), plus element type. ---*/

                                if (elemTypeVTK[j - 1][s - 1] == -1) elemIndex[j - 1][s - 1] = 9;

                                /*--- Allocate memory for accessing the connectivity and to
                                store it in the proper data structure for post-processing. ---*/

                                connElemTemp = new cgsize_t[dataSize[j - 1][s - 1]];

                                connElems[j - 1][s - 1] = new cgsize_t*[elemIndex[j - 1][s - 1]];
                                for (int jj = 0; jj < elemIndex[j - 1][s - 1]; jj++) {
                                    connElems[j - 1][s - 1][jj] = new cgsize_t[nElems[j - 1][s - 1]];
                                }

                                /*--- Retrieve the connectivity information and store. ---*/

                                if (cg_elements_read(fn, i, j, s, connElemTemp, parentData))
                                    cg_error_exit();

                                /*--- Copy these values into the larger array for
                                storage until writing the SU2 file. ---*/

                                if (elemTypeVTK[j - 1][s - 1] == -1) {
                                    int counter = 0;
                                    for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) {
                                        ElementType_t elmt_type = ElementType_t(connElemTemp[counter]);
                                        cg_npe(elmt_type, &npe);
                                        counter++;
                                        connElems[j - 1][s - 1][0][ii] = elmt_type;
                                        for (int jj = 0; jj < npe; jj++) {
                                            connElems[j - 1][s - 1][jj + 1][ii] = connElemTemp[counter] - 1;
                                            counter++;
                                        }
                                    }
                                }
                                else {
                                    int counter = 0;
                                    for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) {
                                        for (int jj = 0; jj < elemIndex[j - 1][s - 1]; jj++) {
                                            connElems[j - 1][s - 1][jj][ii] = connElemTemp[counter] - 1;
                                            counter++;
                                        }
                                    }
                                }
                                delete[] connElemTemp;

                            } // end master

                        }
                        else {

                            /*--- These are internal elems. Allocate memory on each proc. ---*/

                            connElemTemp = new cgsize_t[nElems[j - 1][s - 1] * connSize];

                            /*--- Copy these values into the larger array for
                            storage until writing the SU2 file. ---*/

                            int counterTemp = 0, counterCGNS = 0;
                            for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) {

                                /*--- Store the conn in chunks of connSize for simplicity. ---*/

                                counterTemp = ii*connSize;

                                /*--- Store the connectivity values. Note we subtract one from
                                the CGNS 1-based convention. We may also need to remove the first
                                entry is this is a mixed element section. ---*/

                                if (isMixed[ii]) counterCGNS++;
                                for (int jj = 0; jj < nPoinPerElem[ii]; jj++) {
                                    connElemTemp[counterTemp] = connElemCGNS[counterCGNS + jj] - 1;
                                    counterTemp++;
                                }
                                counterCGNS += nPoinPerElem[ii];

                            }

                            /*--- Free the memory for the conn. from the CGNS file. ---*/

                            delete[] connElemCGNS;
                            delete[] isMixed;

                            /*--- We now have the connectivity stored in linearly partitioned
                            chunks. We need to loop through and decide how many elements we
                            must send to each rank in order to have all elements that
                            surround a particular "owned" node on each rank (i.e., elements
                            will appear on multiple ranks). First, initialize a counter
                            and flag. ---*/

                            int *nElem_Send = new int[size + 1]; nElem_Send[0] = 0;
                            int *nElem_Recv = new int[size + 1]; nElem_Recv[0] = 0;
                            int *nElem_Flag = new int[size];

                            for (int ii = 0; ii < size; ii++) {
                                nElem_Send[ii] = 0;
                                nElem_Recv[ii] = 0;
                                nElem_Flag[ii] = -1;
                            }
                            nElem_Send[size] = 0; nElem_Recv[size] = 0;

                            for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) {
                                for (int jj = 0; jj < nPoinPerElem[ii]; jj++) {

                                    /*--- Get the index of the current point. ---*/

                                    iPoint = connElemTemp[ii*connSize + jj];

                                    /*--- Search for the processor that owns this point ---*/

                                    iProcessor = iPoint / npoint_procs[0];
                                    if (iProcessor >= size) iProcessor = size - 1;
                                    if (iPoint >= nPoint_Linear[iProcessor])
                                        while (iPoint >= nPoint_Linear[iProcessor + 1]) iProcessor++;
                                    else
                                        while (iPoint < nPoint_Linear[iProcessor])   iProcessor--;

                                    /*--- If we have not visted this element yet, increment our
                                    number of elements that must be sent to a particular proc. ---*/

                                    if (nElem_Flag[iProcessor] != ii) {
                                        nElem_Flag[iProcessor] = ii;
                                        nElem_Send[iProcessor + 1]++;
                                    }

                                }
                            }

                            /*--- Communicate the number of cells to be sent/recv'd amongst
                            all processors. After this communication, each proc knows how
                            many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
                            MPI_Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                                &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
                            nElem_Recv[1] = nElem_Send[1];
#endif

                            /*--- Prepare to send connectivities. First check how many
                            messages we will be sending and receiving. Here we also put
                            the counters into cumulative storage format to make the
                            communications simpler. ---*/

                            int nSends = 0, nRecvs = 0;
                            for (int ii = 0; ii < size; ii++) nElem_Flag[ii] = -1;

                            for (int ii = 0; ii < size; ii++) {

                                if ((ii != rank) && (nElem_Send[ii + 1] > 0)) nSends++;
                                if ((ii != rank) && (nElem_Recv[ii + 1] > 0)) nRecvs++;

                                nElem_Send[ii + 1] += nElem_Send[ii];
                                nElem_Recv[ii + 1] += nElem_Recv[ii];
                            }

                            /*--- Allocate memory to hold the connectivity that we are
                            sending. Note that we are also sending the VTK element type
                            in the first position and also the global ID. We have assumed
                            a constant message size of a hex element + 2 extra vals. ---*/

                            unsigned long *connSend = NULL;
                            connSend = new unsigned long[connSize*nElem_Send[size]];
                            for (int ii = 0; ii < connSize*nElem_Send[size]; ii++)
                                connSend[ii] = 0;

                            /*--- Create an index variable to keep track of our index
                            position as we load up the send buffer. ---*/

                            unsigned long *index = new unsigned long[size];
                            for (int ii = 0; ii < size; ii++) index[ii] = connSize*nElem_Send[ii];

                            /*--- Loop through our elements and load the elems and their
                            additional data that we will send to the other procs. ---*/

                            for (int ii = 0; ii < nElems[j - 1][s - 1]; ii++) {
                                for (int jj = 0; jj < nPoinPerElem[ii]; jj++) {

                                    /*--- Get the index of the current point. ---*/

                                    iPoint = connElemTemp[ii*connSize + jj];

                                    /*--- Search for the processor that owns this point ---*/

                                    iProcessor = iPoint / npoint_procs[0];
                                    if (iProcessor >= size) iProcessor = size - 1;
                                    if (iPoint >= nPoint_Linear[iProcessor])
                                        while (iPoint >= nPoint_Linear[iProcessor + 1]) iProcessor++;
                                    else
                                        while (iPoint < nPoint_Linear[iProcessor])   iProcessor--;

                                    /*--- Load connectivity into the buffer for sending ---*/

                                    if (nElem_Flag[iProcessor] != ii) {

                                        nElem_Flag[iProcessor] = ii;
                                        unsigned long nn = index[iProcessor];

                                        /*--- Load the VTK type first into the conn array,
                                        then the connectivity vals, and last, the global ID. ---*/

                                        connSend[nn] = elemTypes[ii]; nn++;
                                        for (int kk = 0; kk < nPoinPerElem[ii]; kk++) {
                                            connSend[nn] = connElemTemp[ii*connSize + kk]; nn++;
                                        }
                                        connSend[nn] = (cgsize_t)elemGlobalID[ii];

                                        /*--- Increment the index by the message length ---*/

                                        index[iProcessor] += connSize;

                                    }
                                }
                            }

                            /*--- Free memory after loading up the send buffer. ---*/

                            delete[] connElemTemp;
                            delete[] elemTypes;
                            delete[] nPoinPerElem;
                            delete[] elemGlobalID;
                            delete[] index;

                            /*--- Allocate the memory that we need for receiving the conn
                            values and then cue up the non-blocking receives. Note that
                            we do not include our own rank in the communications. We will
                            directly copy our own data later. ---*/

                            unsigned long *connRecv = NULL;
                            connRecv = new unsigned long[connSize*nElem_Recv[size]];
                            for (int ii = 0; ii < connSize*nElem_Recv[size]; ii++)
                                connRecv[ii] = 0;

#ifdef HAVE_MPI
                            send_req = new MPI_Request[nSends];
                            recv_req = new MPI_Request[nRecvs];
                            unsigned long iMessage = 0;
                            for (int ii = 0; ii < size; ii++) {
                                if ((ii != rank) && (nElem_Recv[ii + 1] > nElem_Recv[ii])) {
                                    int ll = connSize*nElem_Recv[ii];
                                    int kk = nElem_Recv[ii + 1] - nElem_Recv[ii];
                                    int count = connSize*kk;
                                    int source = ii;
                                    int tag = ii + 1;
                                    MPI_Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                                        MPI_COMM_WORLD, &(recv_req[iMessage]));
                                    iMessage++;
                                }
                            }

                            /*--- Launch the non-blocking sends of the connectivity. ---*/

                            iMessage = 0;
                            for (int ii = 0; ii < size; ii++) {
                                if ((ii != rank) && (nElem_Send[ii + 1] > nElem_Send[ii])) {
                                    int ll = connSize*nElem_Send[ii];
                                    int kk = nElem_Send[ii + 1] - nElem_Send[ii];
                                    int count = connSize*kk;
                                    int dest = ii;
                                    int tag = rank + 1;
                                    MPI_Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                                        MPI_COMM_WORLD, &(send_req[iMessage]));
                                    iMessage++;
                                }
                            }
#endif

                            /*--- Copy my own rank's data into the recv buffer directly. ---*/

                            int mm = connSize*nElem_Recv[rank];
                            int ll = connSize*nElem_Send[rank];
                            int kk = connSize*nElem_Send[rank + 1];

                            for (int nn = ll; nn < kk; nn++, mm++) connRecv[mm] = connSend[nn];

                            /*--- Wait for the non-blocking sends and recvs to complete ---*/

#ifdef HAVE_MPI
                            int number = nSends;
                            for (int ii = 0; ii < nSends; ii++)
                                MPI_Waitany(number, send_req, &ind, &status);

                            number = nRecvs;
                            for (int ii = 0; ii < nRecvs; ii++)
                                MPI_Waitany(number, recv_req, &ind, &status);

                            delete[] send_req;
                            delete[] recv_req;
#endif

                            /*--- Store the connectivity for this rank in the proper data
                            structure before post-processing below. First, allocate the
                            appropriate amount of memory for this section. ---*/

                            connElems[j - 1][s - 1] = new cgsize_t*[connSize];
                            for (int jj = 0; jj < connSize; jj++) {
                                connElems[j - 1][s - 1][jj] = new cgsize_t[nElem_Recv[size]];
                            }

                            for (int ii = 0; ii < nElem_Recv[size]; ii++) {
                                for (int jj = 0; jj < connSize; jj++) {
                                    connElems[j - 1][s - 1][jj][ii] = (cgsize_t)connRecv[ii*connSize + jj];
                                }
                            }

                            /*--- Store the total number of elements I now have for
                            the current section after completing the communications. ---*/

                            nElems[j - 1][s - 1] = nElem_Recv[size];

                            /*--- Free temporary memory from communications ---*/

                            delete[] connSend;
                            delete[] connRecv;
                            delete[] nElem_Recv;
                            delete[] nElem_Send;
                            delete[] nElem_Flag;

                        }

                    } // end section

                } // end zone

            } // end database

            /*--- Close the CGNS file. ---*/

            if (cg_close(fn)) cg_error_exit();
            if (rank == MASTER_NODE)
                cout << "Successfully closed the CGNS file." << endl;

            /*--- Load the data from the CGNS file into SU2 memory. ---*/

            if (rank == MASTER_NODE)
                cout << endl << "Loading CGNS data into SU2 data structures." << endl;

            /*--- Read the dimension of the problem ---*/

            nDim = cell_dim;
            if (rank == MASTER_NODE) {
                if (nDim == 2) cout << "Two dimensional problem." << endl;
                if (nDim == 3) cout << "Three dimensional problem." << endl;
            }

            /*--- Initialize some arrays for the adjacency information (ParMETIS). ---*/

            unsigned long *adj_counter = new unsigned long[local_node];
            unsigned long **adjacent_elem = new unsigned long*[local_node];
            for (iPoint = 0; iPoint < local_node; iPoint++) {
                adjacent_elem[iPoint] = new unsigned long[2000];
                adj_counter[iPoint] = 0;
            }

            /*--- Store the total number of interior elements (global). ---*/

#ifdef HAVE_MPI
            Local_nElem = interiorElems;
            MPI_Allreduce(&Local_nElem, &Global_nElem, 1, MPI_UNSIGNED_LONG,
                MPI_SUM, MPI_COMM_WORLD);
            nElem = Global_nElem;
#else
            Global_nElem = interiorElems;
            nElem = Global_nElem;
#endif

            if ((rank == MASTER_NODE) && (size > SINGLE_NODE)) {
                cout << nElem << " interior elements before linear partitioning." << endl;
            }
            else if (rank == MASTER_NODE) {
                cout << nElem << " interior elements." << endl;
            }

            /*--- Set up the global to local element mapping. ---*/

            Global_to_local_elem = new long[nElem];
            for (unsigned long i = 0; i < nElem; i++) {
                Global_to_local_elem[i] = -1;
            }

            /*--- Allocate space for elements. We allocate enough for all interior
            elements globally, but we will only instantiate our local set. ---*/

            elem = new GRID::GRID_Primal*[nElem];
            for (int iElem = 0; iElem < nElem; iElem++) elem[iElem] = NULL;
            loc_element_count = 0; ielem = 0;
            unsigned long global_id = 0;

            /*--- Loop over all the internal, local volumetric elements. ---*/

            for (int k = 0; k < nzones; k++) {
                for (int s = 0; s < nsections; s++) {
                    if (isInternal[k][s]) {
                        for (int i = 0; i < nElems[k][s]; i++) {

                            /*--- Get the VTK type for this element. This is stored in the
                            first entry of the connectivity structure. ---*/

                            VTK_Type = connElems[k][s][0][i];

                            /*--- Instantiate this element and build adjacency structure. ---*/

                            switch (VTK_Type) {
                            case TRIANGLE:
                                for (int j = 0; j < N_POINTS_TRIANGLE; j++) {
                                    vnodes_cgns[j] = connElems[k][s][j + 1][i];
                                }
                                global_id = connElems[k][s][N_POINTS_TRIANGLE + 1][i];
                                for (unsigned short ii = 0; ii < N_POINTS_TRIANGLE; ii++) {
                                    if ((vnodes_cgns[ii] >= starting_node[rank]) && (vnodes_cgns[ii] < ending_node[rank])) {
                                        for (unsigned short j = 0; j < N_POINTS_TRIANGLE; j++) {
                                            if (ii != j) {
                                                adjacent_elem[vnodes_cgns[ii] - starting_node[rank]][adj_counter[vnodes_cgns[ii] - starting_node[rank]]] = vnodes_cgns[j];
                                                adj_counter[vnodes_cgns[ii] - starting_node[rank]]++;
                                            }
                                        }
                                    }
                                }
                                Global_to_local_elem[global_id] = ielem;
                                elem[ielem] = new GRID::GRID_Triangle(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], nDim);
                                ielem++; nelem_triangle++; break;
                            case RECTANGLE:
                                for (int j = 0; j < N_POINTS_QUADRILATERAL; j++) {
                                    vnodes_cgns[j] = connElems[k][s][j + 1][i];
                                }
                                global_id = connElems[k][s][N_POINTS_QUADRILATERAL + 1][i];
                                for (unsigned short ii = 0; ii < N_POINTS_QUADRILATERAL; ii++) {
                                    if ((vnodes_cgns[ii] >= starting_node[rank]) && (vnodes_cgns[ii] < ending_node[rank])) {
                                        for (unsigned short j = 0; j < N_POINTS_QUADRILATERAL; j++) {
                                            if (ii != j) {
                                                adjacent_elem[vnodes_cgns[ii] - starting_node[rank]][adj_counter[vnodes_cgns[ii] - starting_node[rank]]] = vnodes_cgns[j];
                                                adj_counter[vnodes_cgns[ii] - starting_node[rank]]++;
                                            }
                                        }
                                    }
                                }
                                Global_to_local_elem[global_id] = ielem;
                                elem[ielem] = new GRID::GRID_Rectangle(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], nDim);
                                ielem++; nelem_quad++;
                                break;
                            case TETRAHEDRON:
                                for (int j = 0; j < N_POINTS_TETRAHEDRON; j++) {
                                    vnodes_cgns[j] = connElems[k][s][j + 1][i];
                                }
                                global_id = connElems[k][s][N_POINTS_TETRAHEDRON + 1][i];
                                for (unsigned short ii = 0; ii < N_POINTS_TETRAHEDRON; ii++) {
                                    if ((vnodes_cgns[ii] >= starting_node[rank]) && (vnodes_cgns[ii] < ending_node[rank])) {
                                        for (unsigned short j = 0; j < N_POINTS_TETRAHEDRON; j++) {
                                            if (ii != j) {
                                                adjacent_elem[vnodes_cgns[ii] - starting_node[rank]][adj_counter[vnodes_cgns[ii] - starting_node[rank]]] = vnodes_cgns[j];
                                                adj_counter[vnodes_cgns[ii] - starting_node[rank]]++;
                                            }
                                        }
                                    }
                                }
                                Global_to_local_elem[global_id] = ielem;
                                elem[ielem] = new GRID::GRID_Tetrahedron(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3]);
                                ielem++; nelem_tetra++; break;
                            case HEXAHEDRON:
                                for (int j = 0; j < N_POINTS_HEXAHEDRON; j++) {
                                    vnodes_cgns[j] = connElems[k][s][j + 1][i];
                                }
                                global_id = connElems[k][s][N_POINTS_HEXAHEDRON + 1][i];
                                for (unsigned short ii = 0; ii < N_POINTS_HEXAHEDRON; ii++) {
                                    if ((vnodes_cgns[ii] >= starting_node[rank]) && (vnodes_cgns[ii] < ending_node[rank])) {
                                        for (unsigned short j = 0; j < N_POINTS_HEXAHEDRON; j++) {
                                            if (ii != j) {
                                                adjacent_elem[vnodes_cgns[ii] - starting_node[rank]][adj_counter[vnodes_cgns[ii] - starting_node[rank]]] = vnodes_cgns[j];
                                                adj_counter[vnodes_cgns[ii] - starting_node[rank]]++;
                                            }
                                        }
                                    }
                                }
                                Global_to_local_elem[global_id] = ielem;
                                elem[ielem] = new GRID::GRID_Hexahedron(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], vnodes_cgns[4], vnodes_cgns[5], vnodes_cgns[6], vnodes_cgns[7]);
                                ielem++; nelem_hexa++;
                                break;
                            case PRISM:
                                for (int j = 0; j < N_POINTS_PRISM; j++) {
                                    vnodes_cgns[j] = connElems[k][s][j + 1][i];
                                }
                                global_id = connElems[k][s][N_POINTS_PRISM + 1][i];
                                for (unsigned short ii = 0; ii < N_POINTS_PRISM; ii++) {
                                    if ((vnodes_cgns[ii] >= starting_node[rank]) && (vnodes_cgns[ii] < ending_node[rank])) {
                                        for (unsigned short j = 0; j < N_POINTS_PRISM; j++) {
                                            if (ii != j) {
                                                adjacent_elem[vnodes_cgns[ii] - starting_node[rank]][adj_counter[vnodes_cgns[ii] - starting_node[rank]]] = vnodes_cgns[j];
                                                adj_counter[vnodes_cgns[ii] - starting_node[rank]]++;
                                            }
                                        }
                                    }
                                }
                                Global_to_local_elem[global_id] = ielem;
                                elem[ielem] = new GRID::GRID_Prism(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], vnodes_cgns[4], vnodes_cgns[5]);
                                ielem++; nelem_prism++;
                                break;
                            case PYRAMID:
                                for (int j = 0; j < N_POINTS_PYRAMID; j++) {
                                    vnodes_cgns[j] = connElems[k][s][j + 1][i];
                                }
                                global_id = connElems[k][s][N_POINTS_PYRAMID + 1][i];
                                for (unsigned short ii = 0; ii < N_POINTS_PYRAMID; ii++) {
                                    if ((vnodes_cgns[ii] >= starting_node[rank]) && (vnodes_cgns[ii] < ending_node[rank])) {
                                        for (unsigned short j = 0; j < N_POINTS_PYRAMID; j++) {
                                            if (ii != j) {
                                                adjacent_elem[vnodes_cgns[ii] - starting_node[rank]][adj_counter[vnodes_cgns[ii] - starting_node[rank]]] = vnodes_cgns[j];
                                                adj_counter[vnodes_cgns[ii] - starting_node[rank]]++;
                                            }
                                        }
                                    }
                                }
                                Global_to_local_elem[global_id] = ielem;
                                elem[ielem] = new GRID::GRID_Pyramid(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], vnodes_cgns[4]);
                                ielem++; nelem_pyramid++;
                                break;
                            default:
                                if (rank == MASTER_NODE)
                                    cout << "Element type not suppported!" << endl;
#ifndef HAVE_MPI
                                exit(EXIT_FAILURE);
#else
                                MPI_Abort(MPI_COMM_WORLD, 1);
                                MPI_Finalize();
#endif
                                break;
                            }
                        }
                    }
                }
            }

#ifdef HAVE_MPI
            Local_nElemTri = nelem_triangle;
            Local_nElemQuad = nelem_quad;
            Local_nElemTet = nelem_tetra;
            Local_nElemHex = nelem_hexa;
            Local_nElemPrism = nelem_prism;
            Local_nElemPyramid = nelem_pyramid;
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

            /*--- Store the number of local elements on each rank after determining
            which elements must be kept in the loop above. ---*/

            no_of_local_elements = ielem;

            /*--- Post process the adjacency information in order to get it into the
            proper format before sending the data to ParMETIS. We need to remove
            repeats and adjust the size of the array for each local node. ---*/

            if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
                cout << "Building the graph adjacency structure." << endl;

            unsigned long loc_adjc_size = 0;
            std::vector<unsigned long> adjac_vec;
            unsigned long adj_elem_size;
            std::vector<unsigned long>::iterator it;
            local_elem = ielem;

            xadj = new unsigned long[npoint_procs[rank] + 1];
            xadj[0] = 0;
            std::vector<unsigned long> temp_adjacency;
            unsigned long local_count = 0;

            for (unsigned long i = 0; i < local_node; i++) {

                for (unsigned long j = 0; j < adj_counter[i]; j++) {
                    temp_adjacency.push_back(adjacent_elem[i][j]);
                }

                sort(temp_adjacency.begin(), temp_adjacency.end());
                it = unique(temp_adjacency.begin(), temp_adjacency.end());
                loc_adjc_size = it - temp_adjacency.begin();

                temp_adjacency.resize(loc_adjc_size);
                xadj[local_count + 1] = xadj[local_count] + loc_adjc_size;
                local_count++;

                for (unsigned long j = 0; j < loc_adjc_size; j++) {
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
            for (iPoint = 0; iPoint < local_node; iPoint++) {
                delete[] adjacent_elem[iPoint];
            }
            delete[] adjacent_elem;

            /*--- Store the nodal coordinates from the linear partitioning. ---*/

            if ((rank == MASTER_NODE) && (size > SINGLE_NODE)) {
                cout << nPoint << " grid points before linear partitioning." << endl;
            }
            else if (rank == MASTER_NODE) {
                cout << nPoint << " grid points." << endl;
            }

            iPoint = 0;
            node = new GRID_DGPoint::GRID_DGPoint*[local_node];
            GlobalIndex = starting_node[rank];
            for (int k = 0; k < nzones; k++) {
                for (int i = 0; i < local_node; i++) {
                    for (int j = 0; j < cell_dim; j++) Coord_cgns[j] = gridCoords[k][j][i];
                    switch (nDim) {
                    case 2:
                        node[iPoint] = new GRID_DGPoint::GRID_DGPoint(Coord_cgns[0], Coord_cgns[1], GlobalIndex, config);
                        iPoint++; break;
                    case 3:
                        node[iPoint] = new GRID_DGPoint::GRID_DGPoint(Coord_cgns[0], Coord_cgns[1], Coord_cgns[2], GlobalIndex, config);
                        iPoint++; break;
                    }
                    GlobalIndex++;
                }
            }

            /*--- For now, the master node takes care of all markers. ---*/

            if (rank == MASTER_NODE) {

                /*--- Read number of markers ---*/

                nMarker = nMarkers;
                cout << nMarker << " surface markers." << endl;
                config->SetnMarker_All(nMarker);
                bound = new GRID::GRID_Primal**[nMarker];
                nElem_Bound = new unsigned long[nMarker];
                Tag_to_Marker = new std::string[nMarker_Max];

                iMarker = 0;
                for (int k = 0; k < nzones; k++) {
                    for (int s = 0; s < nsections; s++) {
                        if (!isInternal[k][s]) {

                            /*--- Initialize some counter variables ---*/

                            nelem_edge_bound = 0; nelem_triangle_bound = 0;
                            nelem_quad_bound = 0; ielem = 0;

                            Marker_Tag = sectionNames[k][s];

                            /*--- Remove whitespaces from the marker names ---*/
                            Marker_Tag.erase(remove(Marker_Tag.begin(), Marker_Tag.end(), ' '), Marker_Tag.end());

                            if (Marker_Tag != "SEND_RECEIVE") {
                                nElem_Bound[iMarker] = nElems[k][s];
                                if (rank == MASTER_NODE) {
                                    cout << nElem_Bound[iMarker] << " boundary elements in index ";
                                    cout << iMarker << " (Marker = " << Marker_Tag << ")." << endl;
                                }
                                bound[iMarker] = new GRID::GRID_Primal*[nElem_Bound[iMarker]];

                                for (int i = 0; i < nElems[k][s]; i++) {

                                    /*--- Get the VTK type for this element. Check for mixed
                                    elements. ---*/

                                    if (elemTypeVTK[k][s] == -1) {

                                        /*--- Mixed-element support. Check the elem type. ---*/

                                        ElementType_t elmt_type = ElementType_t(connElems[k][s][0][i]);
                                        cg_npe(elmt_type, &npe);

                                        switch (elmt_type) {
                                        case NODE:    VTK_Type = 1;  break;
                                        case BAR_2:   VTK_Type = 3;  break;
                                        case BAR_3:   VTK_Type = 3;  break;
                                        case TRI_3:   VTK_Type = 5;  break;
                                        case QUAD_4:  VTK_Type = 9;  break;
                                        case TETRA_4: VTK_Type = 10; break;
                                        case HEXA_8:  VTK_Type = 12; break;
                                        case PENTA_6: VTK_Type = 13; break;
                                        case PYRA_5:  VTK_Type = 14; break;
                                        default:
                                            cout << "Kind of element not suppported!" << endl;
#ifndef HAVE_MPI
                                            exit(EXIT_FAILURE);
#else
                                            MPI_Abort(MPI_COMM_WORLD, 1);
                                            MPI_Finalize();
#endif
                                            break;
                                        }

                                        /*--- Transfer the nodes for this element. ---*/

                                        for (int j = 1; j < npe + 1; j++) {
                                            vnodes_cgns[j - 1] = connElems[k][s][j][i];
                                        }

                                    }
                                    else {

                                        /*--- Not a mixed section. We know the element type. ---*/

                                        VTK_Type = elemTypeVTK[k][s];

                                        /*--- Transfer the nodes for this element. ---*/

                                        for (int j = 0; j < elemIndex[k][s]; j++) {
                                            vnodes_cgns[j] = connElems[k][s][j][i];
                                        }

                                    }

                                    /*--- Instantiate the boundary elements. ---*/

                                    switch (VTK_Type) {
                                    case LINE:
                                        if (nDim == 3) {
                                            if (rank == MASTER_NODE) {
                                                printf("\n\n   !!! Error !!!\n");
                                                printf("Remove line boundary elems from the mesh.");
                                            }
#ifndef HAVE_MPI
                                            exit(EXIT_FAILURE);
#else
                                            MPI_Abort(MPI_COMM_WORLD, 1);
                                            MPI_Finalize();
#endif
                                        }
                                        bound[iMarker][ielem] = new GRID::GRID_Line(vnodes_cgns[0], vnodes_cgns[1], 2);
                                        ielem++; nelem_edge_bound++; break;
                                    case TRIANGLE:
                                        bound[iMarker][ielem] = new GRID::GRID_Triangle(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], 3);
                                        ielem++; nelem_triangle_bound++; break;
                                    case RECTANGLE:
                                        bound[iMarker][ielem] = new GRID::GRID_Rectangle(vnodes_cgns[0], vnodes_cgns[1], vnodes_cgns[2], vnodes_cgns[3], 3);
                                        ielem++; nelem_quad_bound++; break;
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
                                config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
                                config->SetMarker_All_SendRecv(iMarker, NONE);

                            }
                            iMarker++;
                        }
                    }
                }

                /*--- Periodic transormations is not implement, store default zeros ---*/
                unsigned short nPeriodic = 1, iPeriodic = 0;
                config->SetnPeriodicIndex(nPeriodic);
                double* center = new double[3];
                double* rotation = new double[3];
                double* translate = new double[3];
                for (unsigned short iDim = 0; iDim < 3; iDim++) {
                    center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
                }
                config->SetPeriodicCenter(iPeriodic, center);
                config->SetPeriodicRotation(iPeriodic, rotation);
                config->SetPeriodicTranslate(iPeriodic, translate);

            }

            /*--- Deallocate temporary memory. ---*/

            delete[] vertices;
            delete[] cells;
            delete[] boundVerts;

            for (int j = 0; j < nzones; j++) {
                delete[] coordArray[j];
                delete[] elemTypeVTK[j];
                delete[] elemIndex[j];
                delete[] nElems[j];
                delete[] dataSize[j];
                delete[] isInternal[j];
                delete[] sectionNames[j];
            }
            delete[] coordArray;
            delete[] elemTypeVTK;
            delete[] elemIndex;
            delete[] nElems;
            delete[] dataSize;
            delete[] isInternal;
            delete[] sectionNames;

            for (int j = 0; j < nzones; j++) {
                for (int i = 0; i < ncoords; i++) {
                    delete[] gridCoords[j][i];
                }
                delete[] gridCoords[j];
            }
            delete[] gridCoords;

            //  for ( int kk = 0; kk < nzones; kk++) {
            //    for (int ii = 0; ii < nsections; ii++) {
            //      for (int jj = 0; jj < indexMax; jj++) {
            //        delete[] connElems[kk][ii][jj];
            //      }
            //      delete connElems[kk][ii];
            //    }
            //    delete connElems[kk];
            //  }
            //  delete[] connElems;

            delete[] nPoint_Linear;
            delete[] nElem_Linear;

            delete[] elemB;
            delete[] elemE;

#else
            std::cout << "SU2 built without CGNS support!!" << std::endl;
            std::cout << "To use CGNS, remove the -DNO_CGNS directive ";
            std::cout << "from the makefile and supply the correct path";
            std::cout << " to the CGNS library." << std::endl;
            exit(EXIT_FAILURE);
#endif

        }


    
}










