
#include "GRID_SurfaceMovement.hpp"
using namespace std;
namespace ARIES
{
    namespace GRID
    {
        GRID_SurfaceMovement::GRID_SurfaceMovement(void) : GRID_Gridmovement() {
            nFFDBox = 0;
            nLevel = 0;
            FFDBoxDefinition = false;
        }

        GRID_SurfaceMovement::~GRID_SurfaceMovement(void) {}

        void GRID_SurfaceMovement::SetSurface_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) 
        {

            unsigned short iFFDBox, iDV, iLevel, iChild, iParent, jFFDBox;
            int rank = TBOX::MASTER_NODE;
            std::string FFDBoxTag;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Setting the Free Form Deformation ---*/

            if (config->GetDesign_Variable(0) == TBOX::FFD_SETTING) {

                /*--- Definition of the FFD deformation class ---*/

                FFDBox = new GRID_FreeFormDefBox*[TBOX::MAX_NUMBER_FFD];

                /*--- Read the FFD information from the config file ---*/

                ReadFFDInfo(geometry, config, FFDBox);

                /*--- If there is a FFDBox in the input file ---*/

                if (nFFDBox != 0) {

                    /*--- If the FFDBox was not defined in the input file ---*/

                    if ((rank == TBOX::MASTER_NODE) && (GetnFFDBox() != 0))
                        cout << endl << "----------------- FFD technique (cartesian -> parametric) ---------------" << endl;

                    /*--- Create a unitary FFDBox as baseline for other FFDBoxes shapes ---*/

                    GRID_FreeFormDefBox FFDBox_unitary(1, 1, 1);
                    FFDBox_unitary.SetUnitCornerPoints();

                    /*--- Compute the control points of the unitary box, in this case the degree is 1 and the order is 2 ---*/

                    FFDBox_unitary.SetControlPoints_Parallelepiped();

                    for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {

                        /*--- Compute the support control points for the final FFD using the unitary box ---*/

                        FFDBox_unitary.SetSupportCP(FFDBox[iFFDBox]);

                        /*--- Compute control points in the support box ---*/

                        FFDBox_unitary.SetSupportCPChange(FFDBox[iFFDBox]);

                        /*--- Compute the parametric coordinates, it also find the points in
                        the FFDBox using the parametrics coordinates ---*/

                        SetParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

                        /*--- Output original FFD FFDBox ---*/

                        if (rank == TBOX::MASTER_NODE) {
                            cout << "Writing a Tecplot file of the FFD boxes." << endl;
                            FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, true);
                        }

                    }

                }

                else {

                    cout << "There are not FFD boxes in the mesh file!!" << endl;
                    exit(EXIT_FAILURE);

                }

            }

            /*--- Free Form deformation based ---*/

            if ((config->GetDesign_Variable(0) == TBOX::FFD_CONTROL_POINT_2D) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_CAMBER_2D) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_THICKNESS_2D) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_CONTROL_POINT) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_DIHEDRAL_ANGLE) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_TWIST_ANGLE) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_ROTATION) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_CONTROL_SURFACE) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_CAMBER) ||
                (config->GetDesign_Variable(0) == TBOX::FFD_THICKNESS)) {

                /*--- Definition of the FFD deformation class ---*/

                FFDBox = new GRID_FreeFormDefBox*[TBOX::MAX_NUMBER_FFD];

                /*--- Read the FFD information from the grid file ---*/

                ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());

                /*--- If there is a FFDBox in the input file ---*/

                if (nFFDBox != 0) {

                    /*--- If the FFDBox was not defined in the input file ---*/

                    if (!GetFFDBoxDefinition()) {

                        cout << endl << "There is not FFD box definition in the mesh file," << endl;
                        cout << "run DV_KIND=FFD_SETTING first !!" << endl;
                        exit(EXIT_FAILURE);

                    }

                    /*--- Output original FFD FFDBox ---*/

                    if (rank == TBOX::MASTER_NODE) {
                        cout << "Writing a Tecplot file of the FFD boxes." << endl;
                        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                            FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, true);
                        }
                    }

                    /*--- Apply the deformation to the orifinal FFD box ---*/

                    if ((rank == TBOX::MASTER_NODE) && (GetnFFDBox() != 0))
                        cout << endl << "----------------- FFD technique (parametric -> cartesian) ---------------" << endl;

                    /*--- Loop over all the FFD boxes levels ---*/

                    for (iLevel = 0; iLevel < GetnLevel(); iLevel++) {

                        /*--- Loop over all FFD FFDBoxes ---*/

                        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {

                            /*--- Check the level of the FFD box ---*/

                            if (FFDBox[iFFDBox]->GetLevel() == iLevel) {


                                /*--- Compute intersections of the FFD box with the surface to eliminate design
                                variables and satisfy surface continuity ---*/

                                if (rank == TBOX::MASTER_NODE)
                                    cout << "Checking FFD box intersections with the solid surfaces." << endl;

                                CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);

                                /*--- Compute the parametric coordinates of the child box
                                control points (using the parent FFDBox)  ---*/

                                for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
                                    FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
                                    for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                                        if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
                                    SetParametricCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
                                }

                                /*--- Update the parametric coordinates if it is a child FFDBox ---*/

                                if (iLevel > 0) UpdateParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

                                /*--- Apply the design variables to the control point position ---*/

                                for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                                    switch (config->GetDesign_Variable(iDV)) {
                                    case TBOX::FFD_CONTROL_POINT_2D: SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_CAMBER_2D:        SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_THICKNESS_2D:     SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_CONTROL_POINT:    SetFFDCPChange(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_DIHEDRAL_ANGLE:   SetFFDDihedralAngle(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_TWIST_ANGLE:      SetFFDTwistAngle(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_ROTATION:         SetFFDRotation(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_CONTROL_SURFACE:  SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_CAMBER:           SetFFDCamber(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    case TBOX::FFD_THICKNESS:        SetFFDThickness(geometry, config, FFDBox[iFFDBox], iDV, false); break;
                                    }
                                }

                                /*--- Recompute cartesian coordinates using the new control point location ---*/

                                SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

                                /*--- Reparametrization of the parent FFD box ---*/

                                for (iParent = 0; iParent < FFDBox[iFFDBox]->GetnParentFFDBox(); iParent++) {
                                    FFDBoxTag = FFDBox[iFFDBox]->GetParentFFDBoxTag(iParent);
                                    for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                                        if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
                                    UpdateParametricCoord(geometry, config, FFDBox[jFFDBox], jFFDBox);
                                }

                                /*--- Compute the new location of the control points of the child boxes
                                (using the parent FFDBox) ---*/

                                for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
                                    FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
                                    for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                                        if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
                                    GetCartesianCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
                                }
                            }
                        }

                        /*--- Output the deformed FFD Boxes ---*/

                        if (rank == TBOX::MASTER_NODE) {
                            cout << "Writing a Tecplot file of the FFD boxes." << endl;
                            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                                FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, false);
                            }
                        }

                    }
                }

                else {

                    cout << "There are not FFD boxes in the mesh file!!" << endl;
                    exit(EXIT_FAILURE);

                }

            }

            /*--- External surface file based ---*/

            else if (config->GetDesign_Variable(0) == TBOX::SURFACE_FILE) {

                /*--- Check whether a surface file exists for input ---*/
                ofstream Surface_File;
                string filename = config->GetMotion_FileName();
                Surface_File.open(filename.c_str(), ios::in);

                /*--- A surface file does not exist, so write a new one for the
                markers that are specified as part of the motion. ---*/
                if (Surface_File.fail()) {

                    if (rank == TBOX::MASTER_NODE)
                        cout << "No surface file found. Writing a new file: " << filename << "." << endl;

                    Surface_File.open(filename.c_str(), ios::out);
                    Surface_File.precision(15);
                    unsigned long iMarker, jPoint, GlobalIndex, iVertex; double *Coords;
                    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                        if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                                jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                                GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
                                Coords = geometry->node[jPoint]->GetCoord();
                                Surface_File << GlobalIndex << "\t" << Coords[0] << "\t" << Coords[1];
                                if (geometry->GetnDim() == 2) Surface_File << endl;
                                else Surface_File << "\t" << Coords[2] << endl;
                            }
                        }
                    }
                    Surface_File.close();

                    /*--- A surface file exists, so read in the coordinates ---*/

                }

                else {
                    Surface_File.close();
                    if (rank == TBOX::MASTER_NODE) cout << "Updating the surface coordinates from the input file." << endl;
                    SetExternal_Deformation(geometry, config, TBOX::ZONE_0, 0);
                }

            }

            /*--- General 2D airfoil deformations ---*/

            else if ((config->GetDesign_Variable(0) == TBOX::ROTATION) ||
                (config->GetDesign_Variable(0) == TBOX::TRANSLATION) ||
                (config->GetDesign_Variable(0) == TBOX::SCALE) ||
                (config->GetDesign_Variable(0) == TBOX::HICKS_HENNE)) {

                /*--- Apply rotation, displacement and stretching design variables (this
                should be done before the bump function design variables) ---*/

                for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                    switch (config->GetDesign_Variable(iDV)) {
                    case TBOX::SCALE:  SetScale(geometry, config, iDV, false); break;
                    case TBOX::TRANSLATION:  SetTranslation(geometry, config, iDV, false); break;
                    case TBOX::ROTATION:     SetRotation(geometry, config, iDV, false); break;
                    }
                }

                /*--- Apply the design variables to the control point position ---*/

                for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                    switch (config->GetDesign_Variable(iDV)) {
                        case TBOX::HICKS_HENNE:  SetHicksHenne(geometry, config, iDV, false); break;
                    }
                }

            }

            /*--- NACA_4Digits design variable ---*/

            else if (config->GetDesign_Variable(0) == TBOX::NACA_4DIGITS) { SetNACA_4Digits(geometry, config); }

            /*--- Parabolic airfoil design variable ---*/

            else if (config->GetDesign_Variable(0) == TBOX::PARABOLIC) { SetParabolic(geometry, config); }

            /*--- Airfoil from file design variable ---*/

            else if (config->GetDesign_Variable(0) == TBOX::AIRFOIL) { SetAirfoil(geometry, config); }

            /*--- FFD setting ---*/

            else if (config->GetDesign_Variable(0) == TBOX::FFD_SETTING) {
                if (rank == TBOX::MASTER_NODE)
                    cout << "No surface deformation (setting FFD)." << endl;
            }

            /*--- Design variable not implement ---*/

            else {
                if (rank == TBOX::MASTER_NODE)
                    cout << "Design Variable not implement yet" << endl;
            }

        }

        void GRID_SurfaceMovement::CopyBoundary(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            unsigned short iMarker;
            unsigned long iVertex, iPoint;
            double *Coord;

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                    Coord = geometry->node[iPoint]->GetCoord();
                    geometry->vertex[iMarker][iVertex]->SetCoord(Coord);
                }
            }

        }

        void GRID_SurfaceMovement::SetParametricCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox) {

            unsigned short iMarker, iDim, iOrder, jOrder, kOrder, lOrder, mOrder, nOrder;
            unsigned long iVertex, iPoint, TotalVertex = 0;
            double *CartCoordNew, *ParamCoord, CartCoord[3], ParamCoordGuess[3], MaxDiff, my_MaxDiff = 0.0, Diff, *Coord;
            int rank;
            unsigned short nDim = geometry->GetnDim();

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Change order and control points reduce the
            complexity of the point inversion (this only works with boxes,
            and we maintain an internal copy) ---*/

            for (iOrder = 0; iOrder < 2; iOrder++) {
                for (jOrder = 0; jOrder < 2; jOrder++) {
                    for (kOrder = 0; kOrder < 2; kOrder++) {

                        lOrder = 0; mOrder = 0; nOrder = 0;
                        if (iOrder == 1) { lOrder = FFDBox->GetlOrder() - 1; }
                        if (jOrder == 1) { mOrder = FFDBox->GetmOrder() - 1; }
                        if (kOrder == 1) { nOrder = FFDBox->GetnOrder() - 1; }

                        Coord = FFDBox->GetCoordControlPoints(lOrder, mOrder, nOrder);

                        FFDBox->SetCoordControlPoints(Coord, iOrder, jOrder, kOrder);

                    }
                }
            }

            FFDBox->SetlOrder(2); FFDBox->SetmOrder(2); FFDBox->SetnOrder(2);
            FFDBox->SetnControlPoints();

            /*--- Point inversion algorithm with a basic box ---*/

            ParamCoordGuess[0] = 0.5; ParamCoordGuess[1] = 0.5; ParamCoordGuess[2] = 0.5;
            CartCoord[0] = 0.0; CartCoord[1] = 0.0; CartCoord[2] = 0.0;

            /*--- Count the number of vertices ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_DV(iMarker) == TBOX::YES)
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
                        TotalVertex++;

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

                if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {

                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                        /*--- Get the cartesian coordinates ---*/

                        for (iDim = 0; iDim < nDim; iDim++)
                            CartCoord[iDim] = geometry->vertex[iMarker][iVertex]->GetCoord(iDim);

                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

                        /*--- If the point is inside the FFD, compute the value of the parametric coordinate ---*/

                        if (FFDBox->GetPointFFD(geometry, config, iPoint)) {

                            /*--- Find the parametric coordinate ---*/

                            ParamCoord = FFDBox->GetParametricCoord_Iterative(iPoint, CartCoord, ParamCoordGuess, config);

                            /*--- If the parametric coordinates are in (0,1) the point belongs to the FFDBox ---*/

                            if (((ParamCoord[0] >= -TBOX::EPS) && (ParamCoord[0] <= 1.0 + TBOX::EPS)) &&
                                ((ParamCoord[1] >= -TBOX::EPS) && (ParamCoord[1] <= 1.0 + TBOX::EPS)) &&
                                ((ParamCoord[2] >= -TBOX::EPS) && (ParamCoord[2] <= 1.0 + TBOX::EPS))) {

                                /*--- Set the value of the parametric coordinate ---*/

                                FFDBox->Set_MarkerIndex(iMarker);
                                FFDBox->Set_VertexIndex(iVertex);
                                FFDBox->Set_PointIndex(iPoint);
                                FFDBox->Set_ParametricCoord(ParamCoord);
                                FFDBox->Set_CartesianCoord(CartCoord);

                                /*--- Compute the cartesian coordinates using the parametric coordinates
                                to check that everithing is right ---*/

                                CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);

                                /*--- Compute max difference between original value and the recomputed value ---*/

                                Diff = 0.0;
                                for (iDim = 0; iDim < nDim; iDim++)
                                    Diff += (CartCoordNew[iDim] - CartCoord[iDim])*(CartCoordNew[iDim] - CartCoord[iDim]);
                                Diff = sqrt(Diff);
                                my_MaxDiff = max(my_MaxDiff, Diff);

                                ParamCoordGuess[0] = ParamCoord[0]; ParamCoordGuess[1] = ParamCoord[1]; ParamCoordGuess[2] = ParamCoord[2];

                            }
                            else {
                                cout << "Please check this point: (" << ParamCoord[0] << " " << ParamCoord[1] << " " << ParamCoord[2] << ") <-> ("
                                    << CartCoord[0] << " " << CartCoord[1] << " " << CartCoord[2] << ")." << endl;
                            }

                        }
                    }
                }
            }

#ifdef HAVE_MPI
            MPI_Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
            MaxDiff = my_MaxDiff;
#endif

            if (rank == TBOX::MASTER_NODE)
                cout << "Compute parametric coord      | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff << "." << endl;


            /*--- After the point inversion, copy the original information back ---*/

            FFDBox->SetOriginalControlPoints();

        }

        void GRID_SurfaceMovement::SetParametricCoordCP(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBoxParent, GRID_FreeFormDefBox *FFDBoxChild) {
            unsigned short iOrder, jOrder, kOrder;
            double *CartCoord, *ParamCoord, ParamCoordGuess[3];
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
                for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
                    for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
                        CartCoord = FFDBoxChild->GetCoordControlPoints(iOrder, jOrder, kOrder);
                        ParamCoord = FFDBoxParent->GetParametricCoord_Iterative(0, CartCoord, ParamCoordGuess, config);
                        FFDBoxChild->SetParCoordControlPoints(ParamCoord, iOrder, jOrder, kOrder);
                    }

            if (rank == TBOX::MASTER_NODE)
                cout << "Compute parametric coord (CP) | FFD parent box: " << FFDBoxParent->GetTag() << ". FFD child box: " << FFDBoxChild->GetTag() << "." << endl;


        }

        void GRID_SurfaceMovement::GetCartesianCoordCP(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBoxParent, GRID_FreeFormDefBox *FFDBoxChild) {
            unsigned short iOrder, jOrder, kOrder, iDim;
            double *CartCoord, *ParamCoord;
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
                for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
                    for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
                        ParamCoord = FFDBoxChild->GetParCoordControlPoints(iOrder, jOrder, kOrder);

                        /*--- Clip the value of the parametric coordinates (just in case)  ---*/
                        for (iDim = 0; iDim < 3; iDim++) {
                            if (ParamCoord[iDim] >= 1.0) ParamCoord[iDim] = 1.0;
                            if (ParamCoord[iDim] <= 0.0) ParamCoord[iDim] = 0.0;
                        }

                        CartCoord = FFDBoxParent->EvalCartesianCoord(ParamCoord);
                        FFDBoxChild->SetCoordControlPoints(CartCoord, iOrder, jOrder, kOrder);
                        FFDBoxChild->SetCoordControlPoints_Copy(CartCoord, iOrder, jOrder, kOrder);

                    }

            if (rank == TBOX::MASTER_NODE)
                cout << "Update cartesian coord (CP)   | FFD parent box: " << FFDBoxParent->GetTag() << ". FFD child box: " << FFDBoxChild->GetTag() << "." << endl;

        }

        void GRID_SurfaceMovement::CheckFFDIntersections(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox) {

            double *Coord_0, *Coord_1;
            unsigned short iMarker, iNode, jNode, lDegree, mDegree, nDegree;
            unsigned long iElem, iPoint, jPoint;

            unsigned short Kind_SU2 = config->GetKind_SU2();

            int rank = TBOX::MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            lDegree = FFDBox->GetlOrder() - 1;
            mDegree = FFDBox->GetmOrder() - 1;
            nDegree = FFDBox->GetnOrder() - 1;

            /*--- Check intersection with plane i=0 ---*/

            double *IPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
            double *IPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, 0, nDegree);
            double *IPlane_Coord_2_A = FFDBox->GetCoordControlPoints(0, mDegree, 0);

            double *IPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
            double *IPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(0, mDegree, 0);
            double *IPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, 0, nDegree);

            /*--- Check intersection with plane i=lDegree ---*/

            double *IPlane_Coord_0_B = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
            double *IPlane_Coord_1_B = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
            double *IPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);

            double *IPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
            double *IPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
            double *IPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);

            /*--- Check intersection with plane j=0 ---*/

            double *JPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
            double *JPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, 0, nDegree);
            double *JPlane_Coord_2_A = FFDBox->GetCoordControlPoints(lDegree, 0, 0);

            double *JPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
            double *JPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
            double *JPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, 0, nDegree);

            /*--- Check intersection with plane j=mDegree ---*/

            double *JPlane_Coord_0_B = FFDBox->GetCoordControlPoints(0, mDegree, 0);
            double *JPlane_Coord_1_B = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
            double *JPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);

            double *JPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
            double *JPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
            double *JPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);

            /*--- Check intersection with plane k=0 ---*/

            double *KPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
            double *KPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, mDegree, 0);
            double *KPlane_Coord_2_A = FFDBox->GetCoordControlPoints(lDegree, 0, 0);

            double *KPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
            double *KPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
            double *KPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, mDegree, 0);

            /*--- Check intersection with plane k=nDegree ---*/

            double *KPlane_Coord_0_B = FFDBox->GetCoordControlPoints(0, 0, nDegree);
            double *KPlane_Coord_1_B = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
            double *KPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);

            double *KPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
            double *KPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
            double *KPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);

            /*--- Loop over all the grid triangles ---*/

            bool IPlane_Intersect_A = false, IPlane_Intersect_B = false;
            bool JPlane_Intersect_A = false, JPlane_Intersect_B = false;
            bool KPlane_Intersect_A = false, KPlane_Intersect_B = false;

            /*--- Only the markers in the moving list ---*/

            for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
                if (((config->GetMarker_All_Moving(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_DEF)) ||
                    ((config->GetMarker_All_DV(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_GEO)) ||
                    ((config->GetMarker_All_DV(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_DOT))) {
                    for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
                        for (iNode = 0; iNode < geometry->bound[iMarker][iElem]->GetnNodes(); iNode++) {
                            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
                            for (jNode = 0; jNode < geometry->bound[iMarker][iElem]->GetnNodes(); jNode++) {
                                jPoint = geometry->bound[iMarker][iElem]->GetNode(jNode);

                                if (jPoint > iPoint) {

                                    Coord_0 = geometry->node[iPoint]->GetCoord();
                                    Coord_1 = geometry->node[jPoint]->GetCoord();

                                    if (!IPlane_Intersect_A) {
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_A, IPlane_Coord_1_A, IPlane_Coord_2_A)) { IPlane_Intersect_A = true; }
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_A_, IPlane_Coord_1_A_, IPlane_Coord_2_A_)) { IPlane_Intersect_A = true; }
                                    }

                                    if (!IPlane_Intersect_B) {
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_B, IPlane_Coord_1_B, IPlane_Coord_2_B)) { IPlane_Intersect_B = true; }
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_B_, IPlane_Coord_1_B_, IPlane_Coord_2_B_)) { IPlane_Intersect_B = true; }
                                    }

                                    if (!JPlane_Intersect_A) {
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_A, JPlane_Coord_1_A, JPlane_Coord_2_A)) { JPlane_Intersect_A = true; }
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_A_, JPlane_Coord_1_A_, JPlane_Coord_2_A_)) { JPlane_Intersect_A = true; }
                                    }

                                    if (!JPlane_Intersect_B) {
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_1_B, JPlane_Coord_2_B)) { JPlane_Intersect_B = true; }
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B_, JPlane_Coord_1_B_, JPlane_Coord_2_B_)) { JPlane_Intersect_B = true; }
                                    }

                                    if (!KPlane_Intersect_A) {
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_A, KPlane_Coord_1_A, KPlane_Coord_2_A)) { KPlane_Intersect_A = true; }
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_A_, KPlane_Coord_1_A_, KPlane_Coord_2_A_)) { KPlane_Intersect_A = true; }
                                    }

                                    if (!KPlane_Intersect_B) {
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_B, KPlane_Coord_1_B, KPlane_Coord_2_B)) { KPlane_Intersect_B = true; }
                                        if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_B_, KPlane_Coord_1_B_, KPlane_Coord_2_B_)) { KPlane_Intersect_B = true; }
                                    }

                                }
                            }
                        }
                    }
                }
            }

            /*--- Comunicate the planes that interesect the surface ---*/

            unsigned short MyCode[6] = { 0, 0, 0, 0, 0, 0 }, Code[6] = { 0, 0, 0, 0, 0, 0 };

            if (IPlane_Intersect_A) MyCode[0] = 1;
            if (IPlane_Intersect_B) MyCode[1] = 1;
            if (JPlane_Intersect_A) MyCode[2] = 1;
            if (JPlane_Intersect_B) MyCode[3] = 1;
            if (KPlane_Intersect_A) MyCode[4] = 1;
            if (KPlane_Intersect_B) MyCode[5] = 1;

#ifdef HAVE_MPI

            /*--- Add MPI_Allreduce information using all the nodes ---*/

            MPI_Allreduce(&MyCode, &Code, 6, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

#else

            Code[0] = MyCode[0]; Code[1] = MyCode[1]; Code[2] = MyCode[2];
            Code[3] = MyCode[3]; Code[4] = MyCode[4]; Code[5] = MyCode[5];

#endif

            if (Code[0] != 0) IPlane_Intersect_A = true; else IPlane_Intersect_A = false;
            if (Code[1] != 0) IPlane_Intersect_B = true; else IPlane_Intersect_B = false;
            if (Code[2] != 0) JPlane_Intersect_A = true; else JPlane_Intersect_A = false;
            if (Code[3] != 0) JPlane_Intersect_B = true; else JPlane_Intersect_B = false;
            if (Code[4] != 0) KPlane_Intersect_A = true; else KPlane_Intersect_A = false;
            if (Code[5] != 0) KPlane_Intersect_B = true; else KPlane_Intersect_B = false;

            /*--- Screen output ---*/

            if (rank == TBOX::MASTER_NODE) {

                if (IPlane_Intersect_A || IPlane_Intersect_B ||
                    JPlane_Intersect_A || JPlane_Intersect_B ||
                    KPlane_Intersect_A || KPlane_Intersect_B) {
                    cout << "The FFD planes ";
                    if (IPlane_Intersect_A) cout << "i=0, ";
                    if (IPlane_Intersect_B) cout << "i=" << lDegree << ", ";
                    if (JPlane_Intersect_A) cout << "j=0, ";
                    if (JPlane_Intersect_B) cout << "j=" << mDegree << ", ";
                    if (KPlane_Intersect_A) cout << "k=0, ";
                    if (KPlane_Intersect_B) cout << "k=" << nDegree << ", ";
                    cout << "intersect solid surfaces." << endl;
                }

            }

            /*--- Fix the FFD planes based on the intersections with solid surfaces,
            and the continuity level, check that we have enough degree for the continuity
            that we are looking for ---*/

            if (IPlane_Intersect_A) { FFDBox->Set_Fix_IPlane(0); FFDBox->Set_Fix_IPlane(1); }
            if (IPlane_Intersect_B) { FFDBox->Set_Fix_IPlane(lDegree); FFDBox->Set_Fix_IPlane(lDegree - 1); }

            if (JPlane_Intersect_A) { FFDBox->Set_Fix_JPlane(0); FFDBox->Set_Fix_JPlane(1); }
            if (JPlane_Intersect_B) { FFDBox->Set_Fix_JPlane(mDegree); FFDBox->Set_Fix_JPlane(mDegree - 1); }

            if (KPlane_Intersect_A) { FFDBox->Set_Fix_KPlane(0); FFDBox->Set_Fix_KPlane(1); }
            if (KPlane_Intersect_B) { FFDBox->Set_Fix_KPlane(nDegree); FFDBox->Set_Fix_KPlane(nDegree - 1); }

            if (config->GetFFD_Continuity() == TBOX::DERIVATIVE_2ND) {

                if ((IPlane_Intersect_A) && (lDegree > 1)) { FFDBox->Set_Fix_IPlane(2); }
                if ((IPlane_Intersect_B) && (lDegree > 1)) { FFDBox->Set_Fix_IPlane(lDegree - 2); }

                if ((JPlane_Intersect_A) && (mDegree > 1)) { FFDBox->Set_Fix_JPlane(2); }
                if ((JPlane_Intersect_B) && (mDegree > 1)) { FFDBox->Set_Fix_JPlane(mDegree - 2); }

                if ((KPlane_Intersect_A) && (nDegree > 1)) { FFDBox->Set_Fix_KPlane(2); }
                if ((KPlane_Intersect_B) && (nDegree > 1)) { FFDBox->Set_Fix_KPlane(nDegree - 2); }

            }

        }

        void GRID_SurfaceMovement::UpdateParametricCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox) {
            unsigned short iMarker, iDim;
            unsigned long iVertex, iPoint, iSurfacePoints;
            double CartCoord[3], *CartCoordNew, *CartCoordOld, *ParamCoord, *var_coord, ParamCoordGuess[3], MaxDiff,
                my_MaxDiff = 0.0, Diff;
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Recompute the parametric coordinates ---*/

            for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {

                /*--- Get the marker of the surface point ---*/

                iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);

                if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {

                    /*--- Get the vertex of the surface point ---*/

                    iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
                    iPoint = FFDBox->Get_PointIndex(iSurfacePoints);

                    /*--- Get the parametric and cartesians coordinates of the
                    surface point (they don't mach) ---*/

                    ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);

                    /*--- Compute and set the cartesian coord using the variation computed
                    with the previous deformation ---*/

                    var_coord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
                    CartCoordOld = geometry->node[iPoint]->GetCoord();
                    for (iDim = 0; iDim < 3; iDim++)
                        CartCoord[iDim] = CartCoordOld[iDim] + var_coord[iDim];
                    FFDBox->Set_CartesianCoord(CartCoord, iSurfacePoints);

                    /*--- Find the parametric coordinate using as ParamCoordGuess the previous value ---*/

                    ParamCoordGuess[0] = ParamCoord[0]; ParamCoordGuess[1] = ParamCoord[1]; ParamCoordGuess[2] = ParamCoord[2];
                    ParamCoord = FFDBox->GetParametricCoord_Iterative(iPoint, CartCoord, ParamCoordGuess, config);

                    /*--- Set the new value of the parametric coordinates ---*/

                    FFDBox->Set_ParametricCoord(ParamCoord, iSurfacePoints);

                    /*--- Compute the cartesian coordinates using the parametric coordinates
                    to check that everithing is right ---*/

                    CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);

                    /*--- Compute max difference between original value and the recomputed value ---*/

                    Diff = 0.0;
                    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                        Diff += (CartCoordNew[iDim] - CartCoord[iDim])*(CartCoordNew[iDim] - CartCoord[iDim]);
                    Diff = sqrt(Diff);
                    my_MaxDiff = max(my_MaxDiff, Diff);

                }
            }

#ifdef HAVE_MPI
            MPI_Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
            MaxDiff = my_MaxDiff;
#endif

            if (rank == TBOX::MASTER_NODE)
                cout << "Update parametric coord       | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff << "." << endl;

        }

        void GRID_SurfaceMovement::SetCartesianCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox, unsigned short iFFDBox) {

            double *CartCoordNew, Diff, my_MaxDiff = 0.0, MaxDiff,
                *ParamCoord, VarCoord[3] = { 0.0, 0.0, 0.0 }, CartCoordOld[3] = { 0.0, 0.0, 0.0 };
            unsigned short iMarker, iDim;
            unsigned long iVertex, iPoint, iSurfacePoints;
            int rank;

            unsigned short nDim = geometry->GetnDim();

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Recompute the cartesians coordinates ---*/

            for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {

                /*--- Get the marker of the surface point ---*/

                iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);

                if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {

                    /*--- Get the vertex of the surface point ---*/

                    iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
                    iPoint = FFDBox->Get_PointIndex(iSurfacePoints);

                    /*--- Set to zero the variation of the coordinates ---*/

                    geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                    /*--- Get the parametric coordinate of the surface point ---*/

                    ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);

                    /*--- Compute the new cartesian coordinate, and set the value in
                    the FFDBox structure ---*/

                    CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);
                    FFDBox->Set_CartesianCoord(CartCoordNew, iSurfacePoints);

                    /*--- Get the original cartesian coordinates of the surface point ---*/

                    for (iDim = 0; iDim < nDim; iDim++) {
                        CartCoordOld[iDim] = geometry->node[iPoint]->GetCoord(iDim);
                    }

                    /*--- Set the value of the variation of the coordinates ---*/

                    Diff = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                        VarCoord[iDim] = CartCoordNew[iDim] - CartCoordOld[iDim];
                        if (fabs(VarCoord[iDim]) <= TBOX::EPS) VarCoord[iDim] = 0.0;
                        Diff += (VarCoord[iDim] * VarCoord[iDim]);
                    }
                    Diff = sqrt(Diff);

                    my_MaxDiff = max(my_MaxDiff, Diff);

                    /*--- Set the variation of the coordinates ---*/

                    geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                }
            }

#ifdef HAVE_MPI
            MPI_Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
            MaxDiff = my_MaxDiff;
#endif

            if (rank == TBOX::MASTER_NODE)
                cout << "Update cartesian coord        | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff << "." << endl;

        }


        void GRID_SurfaceMovement::SetFFDCPChange_2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {

            double movement[3], Ampl;
            unsigned short index[3], i, j;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                /*--- Compute deformation ---*/

                Ampl = config->GetDV_Value(iDV);

                movement[0] = config->GetParamDV(iDV, 3)*Ampl;
                movement[1] = config->GetParamDV(iDV, 4)*Ampl;
                movement[2] = 0.0;

                index[0] = int(config->GetParamDV(iDV, 1));
                index[1] = int(config->GetParamDV(iDV, 2));

                /*--- Lower surface ---*/

                index[2] = 0;

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1)) {
                    for (j = 0; j < FFDBox->GetmOrder(); j++) {
                        index[1] = j;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        for (j = 0; j < FFDBox->GetmOrder(); j++) {
                            index[1] = j;
                            FFDBox->SetControlPoints(index, movement);
                        }
                    }
                }
                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1)) {

                    FFDBox->SetControlPoints(index, movement);
                }


                /*--- Upper surface ---*/

                index[2] = 1;

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1)) {
                    for (j = 0; j < FFDBox->GetmOrder(); j++) {
                        index[1] = j;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        for (j = 0; j < FFDBox->GetmOrder(); j++) {
                            index[1] = j;
                            FFDBox->SetControlPoints(index, movement);
                        }
                    }
                }
                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1)) {

                    FFDBox->SetControlPoints(index, movement);
                }
            }

        }

        void GRID_SurfaceMovement::SetFFDCPChange(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {

            double movement[3], Ampl;
            unsigned short index[3], i, j, k, iPlane;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                /*--- Compute deformation ---*/

                Ampl = config->GetDV_Value(iDV);

                movement[0] = config->GetParamDV(iDV, 4)*Ampl;
                movement[1] = config->GetParamDV(iDV, 5)*Ampl;
                movement[2] = config->GetParamDV(iDV, 6)*Ampl;

                index[0] = int(config->GetParamDV(iDV, 1));
                index[1] = int(config->GetParamDV(iDV, 2));
                index[2] = int(config->GetParamDV(iDV, 3));

                /*--- Check that it is possible to move the control point ---*/

                for (iPlane = 0; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
                    if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return;
                }

                for (iPlane = 0; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
                    if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return;
                }

                for (iPlane = 0; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
                    if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return;
                }

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1) &&
                    (int(config->GetParamDV(iDV, 3)) != -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1) &&
                    (int(config->GetParamDV(iDV, 3)) != -1)) {
                    for (j = 0; j < FFDBox->GetmOrder(); j++) {
                        index[1] = j;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1) &&
                    (int(config->GetParamDV(iDV, 3)) == -1)) {
                    for (k = 0; k < FFDBox->GetnOrder(); k++) {
                        index[2] = k;
                        FFDBox->SetControlPoints(index, movement);
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1) &&
                    (int(config->GetParamDV(iDV, 3)) != -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        for (j = 0; j < FFDBox->GetmOrder(); j++) {
                            index[1] = j;
                            FFDBox->SetControlPoints(index, movement);
                        }
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) == -1) &&
                    (int(config->GetParamDV(iDV, 3)) == -1)) {
                    for (j = 0; j < FFDBox->GetmOrder(); j++) {
                        index[1] = j;
                        for (k = 0; k < FFDBox->GetnOrder(); k++) {
                            index[2] = k;
                            FFDBox->SetControlPoints(index, movement);
                        }
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) == -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1) &&
                    (int(config->GetParamDV(iDV, 3)) == -1)) {
                    for (i = 0; i < FFDBox->GetlOrder(); i++) {
                        index[0] = i;
                        for (k = 0; k < FFDBox->GetnOrder(); k++) {
                            index[2] = k;
                            FFDBox->SetControlPoints(index, movement);
                        }
                    }
                }

                if ((int(config->GetParamDV(iDV, 1)) != -1) &&
                    (int(config->GetParamDV(iDV, 2)) != -1) &&
                    (int(config->GetParamDV(iDV, 3)) != -1)) {
                    FFDBox->SetControlPoints(index, movement);
                }

            }

        }

        void GRID_SurfaceMovement::SetFFDCamber_2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            double Ampl, movement[3];
            unsigned short index[3], kIndex;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                for (kIndex = 0; kIndex < 2; kIndex++) {

                    Ampl = config->GetDV_Value(iDV);

                    movement[0] = 0.0;
                    if (kIndex == 0) movement[1] = Ampl;
                    else movement[1] = Ampl;
                    movement[2] = 0.0;

                    index[0] = int(config->GetParamDV(iDV, 1)); index[1] = kIndex; index[2] = 0;
                    FFDBox->SetControlPoints(index, movement);

                    index[2] = 1;
                    FFDBox->SetControlPoints(index, movement);

                }

            }

        }

        void GRID_SurfaceMovement::SetFFDThickness_2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            double Ampl, movement[3];
            unsigned short index[3], kIndex;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                for (kIndex = 0; kIndex < 2; kIndex++) {

                    Ampl = config->GetDV_Value(iDV);

                    movement[0] = 0.0;
                    if (kIndex == 0) movement[1] = -Ampl;
                    else movement[1] = Ampl;
                    movement[2] = 0.0;

                    index[0] = int(config->GetParamDV(iDV, 1)); index[1] = kIndex; index[2] = 0;
                    FFDBox->SetControlPoints(index, movement);

                    index[2] = 1;
                    FFDBox->SetControlPoints(index, movement);

                }

            }

        }

        void GRID_SurfaceMovement::SetFFDCamber(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            double Ampl, movement[3];
            unsigned short index[3], kIndex;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                for (kIndex = 0; kIndex < 2; kIndex++) {

                    Ampl = config->GetDV_Value(iDV);

                    index[0] = int(config->GetParamDV(iDV, 1));
                    index[1] = int(config->GetParamDV(iDV, 2));
                    index[2] = kIndex;

                    movement[0] = 0.0; movement[1] = 0.0;
                    if (kIndex == 0) movement[2] = Ampl;
                    else movement[2] = Ampl;

                    FFDBox->SetControlPoints(index, movement);

                }

            }

        }

        void GRID_SurfaceMovement::SetFFDThickness(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            double Ampl, movement[3];
            unsigned short index[3], kIndex;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                for (kIndex = 0; kIndex < 2; kIndex++) {

                    Ampl = config->GetDV_Value(iDV);

                    index[0] = int(config->GetParamDV(iDV, 1));
                    index[1] = int(config->GetParamDV(iDV, 2));
                    index[2] = kIndex;

                    movement[0] = 0.0; movement[1] = 0.0;
                    if (kIndex == 0) movement[2] = -Ampl;
                    else movement[2] = Ampl;

                    FFDBox->SetControlPoints(index, movement);

                }

            }

        }


        void GRID_SurfaceMovement::SetFFDDihedralAngle(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            unsigned short iOrder, jOrder, kOrder, index[3];
            double movement[3], theta;
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                /*--- The angle of rotation. ---*/

                theta = config->GetDV_Value(iDV)*TBOX::PI_NUMBER / 180.0;

                /*--- Change the value of the control point if move is true ---*/
                for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
                    for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
                        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
                            index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
                            double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
                            movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1] * tan(theta);

                            FFDBox->SetControlPoints(index, movement);
                        }

            }

        }

        void GRID_SurfaceMovement::SetFFDTwistAngle(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            unsigned short iOrder, jOrder, kOrder;
            double  x, y, z, movement[3];
            unsigned short index[3];
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                /*--- xyz-coordinates of a point on the line of rotation. ---*/

                double a = config->GetParamDV(iDV, 1);
                double b = config->GetParamDV(iDV, 2);
                double c = config->GetParamDV(iDV, 3);

                /*--- xyz-coordinate of the line's direction vector. ---*/

                double u = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
                double v = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);
                double w = config->GetParamDV(iDV, 6) - config->GetParamDV(iDV, 3);

                /*--- The angle of rotation. ---*/

                double theta = config->GetDV_Value(iDV)*TBOX::PI_NUMBER / 180.0;

                /*--- An intermediate value used in computations. ---*/

                double u2 = u*u; double v2 = v*v; double w2 = w*w;
                double l2 = u2 + v2 + w2; double l = sqrt(l2);
                double cosT; double sinT;

                /*--- Change the value of the control point if move is true ---*/

                for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
                    for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
                        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
                            index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
                            double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
                            x = coord[0]; y = coord[1]; z = coord[2];

                            double factor = 0.0;
                            if (y < config->GetParamDV(iDV, 2))
                                factor = 0.0;
                            if ((y >= config->GetParamDV(iDV, 2)) && (y <= config->GetParamDV(iDV, 5)))
                                factor = (y - config->GetParamDV(iDV, 2)) / (config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2));
                            if (y > config->GetParamDV(iDV, 5))
                                factor = 1.0;

                            cosT = cos(theta*factor);
                            sinT = sin(theta*factor);

                            movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
                                + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
                                + l*(-c*v + b*w - w*y + v*z)*sinT;
                            movement[0] = movement[0] / l2 - x;

                            movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
                                + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
                                + l*(c*u - a*w + w*x - u*z)*sinT;
                            movement[1] = movement[1] / l2 - y;

                            movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
                                + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
                                + l*(-b*u + a*v - v*x + u*y)*sinT;
                            movement[2] = movement[2] / l2 - z;

                            FFDBox->SetControlPoints(index, movement);

                        }

            }

        }


        void GRID_SurfaceMovement::SetFFDRotation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            unsigned short iOrder, jOrder, kOrder;
            double  movement[3], x, y, z;
            unsigned short index[3];
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                /*--- xyz-coordinates of a point on the line of rotation. ---*/

                double a = config->GetParamDV(iDV, 1);
                double b = config->GetParamDV(iDV, 2);
                double c = config->GetParamDV(iDV, 3);

                /*--- xyz-coordinate of the line's direction vector. ---*/

                double u = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
                double v = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);
                double w = config->GetParamDV(iDV, 6) - config->GetParamDV(iDV, 3);

                /*--- The angle of rotation. ---*/

                double theta = config->GetDV_Value(iDV)*TBOX::PI_NUMBER / 180.0;

                /*--- An intermediate value used in computations. ---*/

                double u2 = u*u; double v2 = v*v; double w2 = w*w;
                double cosT = cos(theta); double sinT = sin(theta);
                double l2 = u2 + v2 + w2; double l = sqrt(l2);

                /*--- Change the value of the control point if move is true ---*/

                for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
                    for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
                        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
                            index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
                            double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
                            x = coord[0]; y = coord[1]; z = coord[2];
                            movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
                                + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
                                + l*(-c*v + b*w - w*y + v*z)*sinT;
                            movement[0] = movement[0] / l2 - x;

                            movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
                                + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
                                + l*(c*u - a*w + w*x - u*z)*sinT;
                            movement[1] = movement[1] / l2 - y;

                            movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
                                + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
                                + l*(-b*u + a*v - v*x + u*y)*sinT;
                            movement[2] = movement[2] / l2 - z;

                            FFDBox->SetControlPoints(index, movement);

                        }
            }

        }

        void GRID_SurfaceMovement::SetFFDControl_Surface(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox *FFDBox,
            unsigned short iDV, bool ResetDef) {
            unsigned short iOrder, jOrder, kOrder;
            double  movement[3], x, y, z;
            unsigned short index[3];
            string design_FFDBox;

            /*--- Set control points to its original value (even if the
            design variable is not in this box) ---*/

            if (ResetDef == true) FFDBox->SetOriginalControlPoints();

            design_FFDBox = config->GetFFDTag(iDV);

            if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {

                /*--- xyz-coordinates of a point on the line of rotation. ---*/

                double a = config->GetParamDV(iDV, 1);
                double b = config->GetParamDV(iDV, 2);
                double c = config->GetParamDV(iDV, 3);

                /*--- xyz-coordinate of the line's direction vector. ---*/

                double u = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
                double v = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);
                double w = config->GetParamDV(iDV, 6) - config->GetParamDV(iDV, 3);

                /*--- The angle of rotation. ---*/

                double theta = -config->GetDV_Value(iDV)*TBOX::PI_NUMBER / 180.0;

                /*--- An intermediate value used in computations. ---*/

                double u2 = u*u; double v2 = v*v; double w2 = w*w;
                double cosT = cos(theta); double sinT = sin(theta);
                double l2 = u2 + v2 + w2; double l = sqrt(l2);

                /*--- Change the value of the control point if move is true ---*/

                for (iOrder = 0; iOrder < FFDBox->GetlOrder() - 2; iOrder++)
                    for (jOrder = 2; jOrder < FFDBox->GetmOrder() - 2; jOrder++)
                        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
                            index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
                            double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
                            x = coord[0]; y = coord[1]; z = coord[2];
                            movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
                                + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
                                + l*(-c*v + b*w - w*y + v*z)*sinT;
                            movement[0] = movement[0] / l2 - x;

                            movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
                                + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
                                + l*(c*u - a*w + w*x - u*z)*sinT;
                            movement[1] = movement[1] / l2 - y;

                            movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
                                + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
                                + l*(-b*u + a*v - v*x + u*y)*sinT;
                            movement[2] = movement[2] / l2 - z;

                            FFDBox->SetControlPoints(index, movement);

                        }
            }

        }

        void GRID_SurfaceMovement::SetHicksHenne(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef) {
            unsigned long iVertex;
            unsigned short iMarker;
            double VarCoord[3] = { 0.0, 0.0, 0.0 }, VarCoord_[3] = { 0.0, 0.0, 0.0 }, *Coord_, *Normal_, ek, fk, BumpSize = 1.0,
                BumpLoc = 0.0, Coord[3] = { 0.0, 0.0, 0.0 }, Normal[3] = { 0.0, 0.0, 0.0 },
                xCoord, TPCoord[2] = { 0.0, 0.0 }, LPCoord[2] = { 0.0, 0.0 }, Distance, Chord, AoA, ValCos, ValSin;

            bool upper = true, double_surface = false;

            /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

            if ((iDV == 0) || (ResetDef == true)) {
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
                    }
            }

            /*--- Compute the angle of attack to apply the deformation ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                    Coord_ = boundary->vertex[iMarker][0]->GetCoord();
                    TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1];
                    for (iVertex = 1; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
                        if (Coord_[0] > TPCoord[0]) { TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1]; }
                    }
                }
            }

#ifdef HAVE_MPI

            unsigned long *Buffer_Send_nVertex, *Buffer_Receive_nVertex;
            int iProcessor, nProcessor;
            double *Buffer_Send_Coord, *Buffer_Receive_Coord;

            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            Buffer_Receive_Coord = new double[nProcessor * 2];
            Buffer_Send_Coord = new double[2];

            Buffer_Send_Coord[0] = TPCoord[0]; Buffer_Send_Coord[1] = TPCoord[1];

            MPI_Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);

            TPCoord[0] = Buffer_Receive_Coord[0]; TPCoord[1] = Buffer_Receive_Coord[1];
            for (iProcessor = 1; iProcessor < nProcessor; iProcessor++) {
                Coord[0] = Buffer_Receive_Coord[iProcessor * 2 + 0];
                Coord[1] = Buffer_Receive_Coord[iProcessor * 2 + 1];
                if (Coord[0] > TPCoord[0]) { TPCoord[0] = Coord[0]; TPCoord[1] = Coord[1]; }
            }

            delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;

#endif


            Chord = 0.0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
                        Distance = sqrt(pow(Coord_[0] - TPCoord[0], 2.0) + pow(Coord_[1] - TPCoord[1], 2.0));
                        if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord_[0]; LPCoord[1] = Coord_[1]; }
                    }
                }
            }

#ifdef HAVE_MPI

            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            Buffer_Receive_Coord = new double[nProcessor * 2];
            Buffer_Send_Coord = new double[2];

            Buffer_Send_Coord[0] = LPCoord[0]; Buffer_Send_Coord[1] = LPCoord[1];

            MPI_Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);

            Chord = 0.0;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                Coord[0] = Buffer_Receive_Coord[iProcessor * 2 + 0];
                Coord[1] = Buffer_Receive_Coord[iProcessor * 2 + 1];
                Distance = sqrt(pow(Coord[0] - TPCoord[0], 2.0) + pow(Coord[1] - TPCoord[1], 2.0));
                if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord[0]; LPCoord[1] = Coord[1]; }
            }

            delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;

#endif

            AoA = atan((LPCoord[1] - TPCoord[1]) / (TPCoord[0] - LPCoord[0])) * 180 / TBOX::PI_NUMBER;
            AoA = 0.0;

            /*--- Perform multiple airfoil deformation ---*/

            double Ampl = config->GetDV_Value(iDV);
            double xk = config->GetParamDV(iDV, 1);
            const double t2 = 3.0;

            if (config->GetParamDV(iDV, 0) == TBOX::NO) { upper = false; double_surface = true; }
            if (config->GetParamDV(iDV, 0) == TBOX::YES) { upper = true; double_surface = true; }

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

                    if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {

                        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
                        Normal_ = boundary->vertex[iMarker][iVertex]->GetNormal();

                        /*--- The Hicks Henne bump functions should be applied to a basic airfoil without AoA,
                        and unitary chord, a tranformation is required ---*/

                        ValCos = cos(AoA*TBOX::PI_NUMBER / 180.0);
                        ValSin = sin(AoA*TBOX::PI_NUMBER / 180.0);

                        Coord[0] = Coord_[0] * ValCos - Coord_[1] * ValSin;
                        Coord[0] = max(0.0, Coord[0]); // Coord x should be always positive
                        Coord[1] = Coord_[1] * ValCos + Coord_[0] * ValSin;

                        Normal[0] = Normal_[0] * ValCos - Normal_[1] * ValSin;
                        Normal[1] = Normal_[1] * ValCos + Normal_[0] * ValSin;

                        /*--- Bump computation ---*/

                        if (double_surface) {
                            ek = log10(0.5) / log10(xk);
                            fk = pow(sin(TBOX::PI_NUMBER * pow(Coord[0], ek)), t2);

                            /*--- Upper and lower surface ---*/

                            if ((upper) && (Normal[1] > 0)) { VarCoord[1] = Ampl*fk; }
                            if ((!upper) && (Normal[1] < 0)) { VarCoord[1] = -Ampl*fk; }

                        }
                        else {
                            xCoord = Coord[0] - BumpLoc;
                            ek = log10(0.5) / log10(xk / BumpSize);
                            fk = pow(sin(TBOX::PI_NUMBER * pow(xCoord / BumpSize, ek)), t2);

                            /*--- Only one surface ---*/

                            if ((xCoord <= 0.0) || (xCoord >= BumpSize)) VarCoord[1] = 0.0;
                            else VarCoord[1] = Ampl*fk;


                        }
                    }

                    /*--- Apply the transformation to the coordinate variation ---*/

                    ValCos = cos(-AoA*TBOX::PI_NUMBER / 180.0);
                    ValSin = sin(-AoA*TBOX::PI_NUMBER / 180.0);

                    VarCoord_[0] = VarCoord[0] * ValCos - VarCoord[1] * ValSin;
                    VarCoord_[1] = VarCoord[1] * ValCos + VarCoord[0] * ValSin;

                    boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord_);

                }
            }

        }

        void GRID_SurfaceMovement::SetRotation(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef) {
            unsigned long iVertex;
            unsigned short iMarker;
            double VarCoord[3], *Coord;
            double  movement[3], x, y, z;

            /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

            if ((iDV == 0) || (ResetDef == true)) {
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
                    }
            }

            /*--- xyz-coordinates of a point on the line of rotation. */

            double a = config->GetParamDV(iDV, 0);
            double b = config->GetParamDV(iDV, 1);
            double c = 0.0;
            if (boundary->GetnDim() == 3) c = config->GetParamDV(0, 2);

            /*--- xyz-coordinate of the line's direction vector. ---*/

            double u = config->GetParamDV(iDV, 3) - config->GetParamDV(iDV, 0);
            double v = config->GetParamDV(iDV, 4) - config->GetParamDV(iDV, 1);
            double w = 1.0;
            if (boundary->GetnDim() == 3) w = config->GetParamDV(iDV, 5) - config->GetParamDV(iDV, 2);

            /*--- The angle of rotation. ---*/

            double theta = config->GetDV_Value(iDV)*TBOX::PI_NUMBER / 180.0;

            /*--- An intermediate value used in computations. ---*/

            double u2 = u*u; double v2 = v*v; double w2 = w*w;
            double cosT = cos(theta); double sinT = sin(theta);
            double l2 = u2 + v2 + w2; double l = sqrt(l2);

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                    if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
                        x = Coord[0]; y = Coord[1]; z = Coord[2];

                        movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
                            + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
                            + l*(-c*v + b*w - w*y + v*z)*sinT;
                        movement[0] = movement[0] / l2 - x;

                        movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
                            + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
                            + l*(c*u - a*w + w*x - u*z)*sinT;
                        movement[1] = movement[1] / l2 - y;

                        movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
                            + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
                            + l*(-b*u + a*v - v*x + u*y)*sinT;
                        if (boundary->GetnDim() == 3) movement[2] = movement[2] / l2 - z;
                        else movement[2] = 0.0;

                        VarCoord[0] = movement[0];
                        VarCoord[1] = movement[1];
                        if (boundary->GetnDim() == 3) VarCoord[2] = movement[2];

                    }
                    boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
                }
        }

        void GRID_SurfaceMovement::SetTranslation(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef) {
            unsigned long iVertex;
            unsigned short iMarker;
            double VarCoord[3];
            double Ampl = config->GetDV_Value(iDV);

            /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

            if ((iDV == 0) || (ResetDef == true)) {
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
                    }
            }

            double xDispl = config->GetParamDV(iDV, 0);
            double yDispl = config->GetParamDV(iDV, 1);
            double zDispl = 0;
            if (boundary->GetnDim() == 3) zDispl = config->GetParamDV(iDV, 2);

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                    if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                        VarCoord[0] = Ampl*xDispl;
                        VarCoord[1] = Ampl*yDispl;
                        if (boundary->GetnDim() == 3) VarCoord[2] = Ampl*zDispl;
                    }
                    boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
                }

        }

        void GRID_SurfaceMovement::SetScale(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config, unsigned short iDV, bool ResetDef) {
            unsigned long iVertex;
            unsigned short iMarker;
            double VarCoord[3], x, y, z, *Coord;
            double Ampl = config->GetDV_Value(iDV);

            /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

            if ((iDV == 0) || (ResetDef == true)) {
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
                    }
            }

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                    if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
                        x = Coord[0]; y = Coord[1]; z = Coord[2];
                        VarCoord[0] = (Ampl - 1.0)*x;
                        VarCoord[1] = (Ampl - 1.0)*y;
                        if (boundary->GetnDim() == 3) VarCoord[2] = (Ampl - 1.0)*z;
                    }
                    boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
                }

        }

        void GRID_SurfaceMovement::Moving_Walls(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            unsigned short iZone, unsigned long iter) {

            int rank = TBOX::MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables ---*/
            unsigned short iMarker, jMarker, iDim, nDim = geometry->GetnDim();
            unsigned long iPoint, iVertex;
            double xDot[3] = { 0.0, 0.0, 0.0 }, *Coord, Center[3] = { 0.0, 0.0, 0.0 }, Omega[3] = { 0.0, 0.0, 0.0 }, r[3] = { 0.0, 0.0, 0.0 }, GridVel[3] = { 0.0, 0.0, 0.0 };
            double L_Ref = config->GetLength_Ref();
            double Omega_Ref = config->GetOmega_Ref();
            double Vel_Ref = config->GetVelocity_Ref();
            string Marker_Tag;

            /*--- Store grid velocity for each node on the moving surface(s).
            Sum and store the x, y, & z velocities due to translation and rotation. ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (config->GetMarker_All_Moving(iMarker) == TBOX::YES) {

                    /*--- Identify iMarker from the list of those under MARKER_MOVING ---*/

                    Marker_Tag = config->GetMarker_All_TagBound(iMarker);
                    jMarker = config->GetMarker_Moving(Marker_Tag);

                    /*--- Get prescribed wall speed from config for this marker ---*/

                    Center[0] = config->GetMotion_Origin_X(jMarker);
                    Center[1] = config->GetMotion_Origin_Y(jMarker);
                    Center[2] = config->GetMotion_Origin_Z(jMarker);
                    Omega[0] = config->GetRotation_Rate_X(jMarker) / Omega_Ref;
                    Omega[1] = config->GetRotation_Rate_Y(jMarker) / Omega_Ref;
                    Omega[2] = config->GetRotation_Rate_Z(jMarker) / Omega_Ref;
                    xDot[0] = config->GetTranslation_Rate_X(jMarker) / Vel_Ref;
                    xDot[1] = config->GetTranslation_Rate_Y(jMarker) / Vel_Ref;
                    xDot[2] = config->GetTranslation_Rate_Z(jMarker) / Vel_Ref;

                    if (rank == TBOX::MASTER_NODE && iter == 0) {
                        cout << " Storing grid velocity for marker: ";
                        cout << Marker_Tag << "." << endl;
                        cout << " Translational velocity: (" << xDot[0] << ", " << xDot[1];
                        cout << ", " << xDot[2] << ") m/s." << endl;
                        cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
                        cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
                        cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
                    }

                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                        /*--- Get the index and coordinates of the current point ---*/

                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        Coord = geometry->node[iPoint]->GetCoord();

                        /*--- Calculate non-dim. position from rotation center ---*/
                        for (iDim = 0; iDim < nDim; iDim++)
                            r[iDim] = (Coord[iDim] - Center[iDim]) / L_Ref;
                        if (nDim == 2) r[nDim] = 0.0;

                        /*--- Cross Product of angular velocity and distance from center to
                        get the rotational velocity. Note that we are adding on the velocity
                        due to pure translation as well. ---*/

                        GridVel[0] = xDot[0] + Omega[1] * r[2] - Omega[2] * r[1];
                        GridVel[1] = xDot[1] + Omega[2] * r[0] - Omega[0] * r[2];
                        GridVel[2] = xDot[2] + Omega[0] * r[1] - Omega[1] * r[0];

                        /*--- Store the moving wall velocity for this node ---*/

                        for (iDim = 0; iDim < nDim; iDim++)
                            geometry->node[iPoint]->SetGridVel(iDim, GridVel[iDim]);

                    }
                }
            }
        }

        void GRID_SurfaceMovement::Surface_Translating(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            unsigned long iter, unsigned short iZone) {

            double deltaT, time_new, time_old;
            double Center[3], VarCoord[3], xDot[3];
            unsigned short iMarker, jMarker, Moving;
            unsigned long iVertex;
            string Marker_Tag, Moving_Tag;
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Initialize the delta variation in coordinates ---*/
            VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

            /*--- Retrieve values from the config file ---*/

            deltaT = config->GetDelta_UnstTimeND();

            /*--- Compute delta time based on physical time step ---*/
            time_new = static_cast<double>(iter)*deltaT;
            if (iter == 0) {
                time_old = time_new;
            }
            else {
                time_old = static_cast<double>(iter - 1)*deltaT;
            }

            /*--- Store displacement of each node on the translating surface ---*/
            /*--- Loop over markers and find the particular marker(s) (surface) to translate ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Moving = config->GetMarker_All_Moving(iMarker);
                if (Moving == TBOX::YES) {
                    for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                        Moving_Tag = config->GetMarker_Moving(jMarker);
                        Marker_Tag = config->GetMarker_All_TagBound(iMarker);

                        if (Marker_Tag == Moving_Tag) {

                            /*--- Translation velocity from config. ---*/

                            xDot[0] = config->GetTranslation_Rate_X(jMarker);
                            xDot[1] = config->GetTranslation_Rate_Y(jMarker);
                            xDot[2] = config->GetTranslation_Rate_Z(jMarker);

                            /*--- Print some information to the console. Be verbose at the first
                            iteration only (mostly for debugging purposes). ---*/
                            // Note that the TBOX::MASTER_NODE might not contain all the markers being moved.

                            if (rank == TBOX::MASTER_NODE) {
                                cout << " Storing translating displacement for marker: ";
                                cout << Marker_Tag << "." << endl;
                                if (iter == 0) {
                                    cout << " Translational velocity: (" << xDot[0] << ", " << xDot[1];
                                    cout << ", " << xDot[2] << ") m/s." << endl;
                                }
                            }

                            /*--- Compute delta change in the position in the x, y, & z directions. ---*/

                            VarCoord[0] = xDot[0] * (time_new - time_old);
                            VarCoord[1] = xDot[1] * (time_new - time_old);
                            VarCoord[2] = xDot[2] * (time_new - time_old);

                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                                /*--- Set node displacement for volume deformation ---*/
                                geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                            }
                        }
                    }
                }
            }

            /*--- When updating the origins it is assumed that all markers have the
            same translational velocity, because we use the last VarCoord set ---*/

            /*--- Set the mesh motion center to the new location after
            incrementing the position with the translation. This new
            location will be used for subsequent mesh motion for the given marker.---*/

            for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                /*-- Check if we want to update the motion origin for the given marker ---*/

                if (config->GetMoveMotion_Origin(jMarker) == TBOX::YES) {
                    Center[0] = config->GetMotion_Origin_X(jMarker) + VarCoord[0];
                    Center[1] = config->GetMotion_Origin_Y(jMarker) + VarCoord[1];
                    Center[2] = config->GetMotion_Origin_Z(jMarker) + VarCoord[2];
                    config->SetMotion_Origin_X(jMarker, Center[0]);
                    config->SetMotion_Origin_Y(jMarker, Center[1]);
                    config->SetMotion_Origin_Z(jMarker, Center[2]);
                }
            }

            /*--- Set the moment computation center to the new location after
            incrementing the position with the translation. ---*/

            for (jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
                Center[0] = config->GetRefOriginMoment_X(jMarker) + VarCoord[0];
                Center[1] = config->GetRefOriginMoment_Y(jMarker) + VarCoord[1];
                Center[2] = config->GetRefOriginMoment_Z(jMarker) + VarCoord[2];
                config->SetRefOriginMoment_X(jMarker, Center[0]);
                config->SetRefOriginMoment_Y(jMarker, Center[1]);
                config->SetRefOriginMoment_Z(jMarker, Center[2]);
            }
        }

        void GRID_SurfaceMovement::Surface_Plunging(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            unsigned long iter, unsigned short iZone) {

            double deltaT, time_new, time_old, Lref;
            double Center[3], VarCoord[3], Omega[3], Ampl[3];
            double DEG2RAD = TBOX::PI_NUMBER / 180.0;
            unsigned short iMarker, jMarker, Moving;
            unsigned long iVertex;
            string Marker_Tag, Moving_Tag;
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Initialize the delta variation in coordinates ---*/
            VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

            /*--- Retrieve values from the config file ---*/

            deltaT = config->GetDelta_UnstTimeND();
            Lref = config->GetLength_Ref();

            /*--- Compute delta time based on physical time step ---*/
            time_new = static_cast<double>(iter)*deltaT;
            if (iter == 0) {
                time_old = time_new;
            }
            else {
                time_old = static_cast<double>(iter - 1)*deltaT;
            }

            /*--- Store displacement of each node on the plunging surface ---*/
            /*--- Loop over markers and find the particular marker(s) (surface) to plunge ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Moving = config->GetMarker_All_Moving(iMarker);
                if (Moving == TBOX::YES) {
                    for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                        Moving_Tag = config->GetMarker_Moving(jMarker);
                        Marker_Tag = config->GetMarker_All_TagBound(iMarker);

                        if (Marker_Tag == Moving_Tag) {

                            /*--- Plunging frequency and amplitude from config. ---*/

                            Omega[0] = config->GetPlunging_Omega_X(jMarker) / config->GetOmega_Ref();
                            Omega[1] = config->GetPlunging_Omega_Y(jMarker) / config->GetOmega_Ref();
                            Omega[2] = config->GetPlunging_Omega_Z(jMarker) / config->GetOmega_Ref();
                            Ampl[0] = config->GetPlunging_Ampl_X(jMarker) / Lref;
                            Ampl[1] = config->GetPlunging_Ampl_Y(jMarker) / Lref;
                            Ampl[2] = config->GetPlunging_Ampl_Z(jMarker) / Lref;

                            /*--- Print some information to the console. Be verbose at the first
                            iteration only (mostly for debugging purposes). ---*/
                            // Note that the TBOX::MASTER_NODE might not contain all the markers being moved.

                            if (rank == TBOX::MASTER_NODE) {
                                cout << " Storing plunging displacement for marker: ";
                                cout << Marker_Tag << "." << endl;
                                if (iter == 0) {
                                    cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
                                    cout << ", " << Omega[2] << ") rad/s." << endl;
                                    cout << " Plunging amplitude: (" << Ampl[0] / DEG2RAD;
                                    cout << ", " << Ampl[1] / DEG2RAD << ", " << Ampl[2] / DEG2RAD;
                                    cout << ") degrees." << endl;
                                }
                            }

                            /*--- Compute delta change in the position in the x, y, & z directions. ---*/

                            VarCoord[0] = -Ampl[0] * (sin(Omega[0] * time_new) - sin(Omega[0] * time_old));
                            VarCoord[1] = -Ampl[1] * (sin(Omega[1] * time_new) - sin(Omega[1] * time_old));
                            VarCoord[2] = -Ampl[2] * (sin(Omega[2] * time_new) - sin(Omega[2] * time_old));

                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                                /*--- Set node displacement for volume deformation ---*/
                                geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                            }
                        }
                    }
                }
            }

            /*--- When updating the origins it is assumed that all markers have the
            same plunging movement, because we use the last VarCoord set ---*/

            /*--- Set the mesh motion center to the new location after
            incrementing the position with the translation. This new
            location will be used for subsequent mesh motion for the given marker.---*/

            for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                /*-- Check if we want to update the motion origin for the given marker ---*/

                if (config->GetMoveMotion_Origin(jMarker) == TBOX::YES) {
                    Center[0] = config->GetMotion_Origin_X(jMarker) + VarCoord[0];
                    Center[1] = config->GetMotion_Origin_Y(jMarker) + VarCoord[1];
                    Center[2] = config->GetMotion_Origin_Z(jMarker) + VarCoord[2];
                    config->SetMotion_Origin_X(jMarker, Center[0]);
                    config->SetMotion_Origin_Y(jMarker, Center[1]);
                    config->SetMotion_Origin_Z(jMarker, Center[2]);
                }
            }

            /*--- Set the moment computation center to the new location after
            incrementing the position with the plunging. ---*/

            for (jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
                Center[0] = config->GetRefOriginMoment_X(jMarker) + VarCoord[0];
                Center[1] = config->GetRefOriginMoment_Y(jMarker) + VarCoord[1];
                Center[2] = config->GetRefOriginMoment_Z(jMarker) + VarCoord[2];
                config->SetRefOriginMoment_X(jMarker, Center[0]);
                config->SetRefOriginMoment_Y(jMarker, Center[1]);
                config->SetRefOriginMoment_Z(jMarker, Center[2]);
            }
        }

        void GRID_SurfaceMovement::Surface_Pitching(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            unsigned long iter, unsigned short iZone) {

            double deltaT, time_new, time_old, Lref, *Coord;
            double Center[3], VarCoord[3], Omega[3], Ampl[3], Phase[3], rotCoord[3], r[3];
            double rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
            double dtheta, dphi, dpsi, cosTheta, sinTheta;
            double cosPhi, sinPhi, cosPsi, sinPsi;
            double DEG2RAD = TBOX::PI_NUMBER / 180.0;
            unsigned short iMarker, jMarker, Moving, iDim, nDim = geometry->GetnDim();
            unsigned long iPoint, iVertex;
            string Marker_Tag, Moving_Tag;
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Initialize the delta variation in coordinates ---*/
            VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

            /*--- Retrieve values from the config file ---*/

            deltaT = config->GetDelta_UnstTimeND();
            Lref = config->GetLength_Ref();

            /*--- Compute delta time based on physical time step ---*/
            time_new = static_cast<double>(iter)*deltaT;
            if (iter == 0) {
                time_old = time_new;
            }
            else {
                time_old = static_cast<double>(iter - 1)*deltaT;
            }

            /*--- Store displacement of each node on the pitching surface ---*/
            /*--- Loop over markers and find the particular marker(s) (surface) to pitch ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Moving = config->GetMarker_All_Moving(iMarker);
                if (Moving == TBOX::YES) {
                    for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                        Moving_Tag = config->GetMarker_Moving(jMarker);
                        Marker_Tag = config->GetMarker_All_TagBound(iMarker);

                        if (Marker_Tag == Moving_Tag) {

                            /*--- Pitching origin, frequency, and amplitude from config. ---*/

                            Center[0] = config->GetMotion_Origin_X(jMarker);
                            Center[1] = config->GetMotion_Origin_Y(jMarker);
                            Center[2] = config->GetMotion_Origin_Z(jMarker);
                            Omega[0] = config->GetPitching_Omega_X(jMarker) / config->GetOmega_Ref();
                            Omega[1] = config->GetPitching_Omega_Y(jMarker) / config->GetOmega_Ref();
                            Omega[2] = config->GetPitching_Omega_Z(jMarker) / config->GetOmega_Ref();
                            Ampl[0] = config->GetPitching_Ampl_X(jMarker)*DEG2RAD;
                            Ampl[1] = config->GetPitching_Ampl_Y(jMarker)*DEG2RAD;
                            Ampl[2] = config->GetPitching_Ampl_Z(jMarker)*DEG2RAD;
                            Phase[0] = config->GetPitching_Phase_X(jMarker)*DEG2RAD;
                            Phase[1] = config->GetPitching_Phase_Y(jMarker)*DEG2RAD;
                            Phase[2] = config->GetPitching_Phase_Z(jMarker)*DEG2RAD;

                            /*--- Print some information to the console. Be verbose at the first
                            iteration only (mostly for debugging purposes). ---*/
                            // Note that the TBOX::MASTER_NODE might not contain all the markers being moved.

                            if (rank == TBOX::MASTER_NODE) {
                                cout << " Storing pitching displacement for marker: ";
                                cout << Marker_Tag << "." << endl;
                                if (iter == 0) {
                                    cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
                                    cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
                                    cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
                                    cout << " Pitching amplitude about origin: (" << Ampl[0] / DEG2RAD;
                                    cout << ", " << Ampl[1] / DEG2RAD << ", " << Ampl[2] / DEG2RAD;
                                    cout << ") degrees." << endl;
                                    cout << " Pitching phase lag about origin: (" << Phase[0] / DEG2RAD;
                                    cout << ", " << Phase[1] / DEG2RAD << ", " << Phase[2] / DEG2RAD;
                                    cout << ") degrees." << endl;
                                }
                            }

                            /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

                            dtheta = -Ampl[0] * (sin(Omega[0] * time_new + Phase[0])
                                - sin(Omega[0] * time_old + Phase[0]));
                            dphi = -Ampl[1] * (sin(Omega[1] * time_new + Phase[1])
                                - sin(Omega[1] * time_old + Phase[1]));
                            dpsi = -Ampl[2] * (sin(Omega[2] * time_new + Phase[2])
                                - sin(Omega[2] * time_old + Phase[2]));

                            /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

                            cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
                            sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);

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

                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                                /*--- Index and coordinates of the current point ---*/

                                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                                Coord = geometry->node[iPoint]->GetCoord();

                                /*--- Calculate non-dim. position from rotation center ---*/

                                for (iDim = 0; iDim < nDim; iDim++)
                                    r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;
                                if (nDim == 2) r[nDim] = 0.0;

                                /*--- Compute transformed point coordinates ---*/

                                rotCoord[0] = rotMatrix[0][0] * r[0]
                                    + rotMatrix[0][1] * r[1]
                                    + rotMatrix[0][2] * r[2] + Center[0];

                                rotCoord[1] = rotMatrix[1][0] * r[0]
                                    + rotMatrix[1][1] * r[1]
                                    + rotMatrix[1][2] * r[2] + Center[1];

                                rotCoord[2] = rotMatrix[2][0] * r[0]
                                    + rotMatrix[2][1] * r[1]
                                    + rotMatrix[2][2] * r[2] + Center[2];

                                /*--- Calculate delta change in the x, y, & z directions ---*/
                                for (iDim = 0; iDim < nDim; iDim++)
                                    VarCoord[iDim] = (rotCoord[iDim] - Coord[iDim]) / Lref;
                                if (nDim == 2) VarCoord[nDim] = 0.0;

                                /*--- Set node displacement for volume deformation ---*/
                                geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                            }
                        }
                    }
                }
            }
            /*--- For pitching we don't update the motion origin and moment reference origin. ---*/
        }

        void GRID_SurfaceMovement::Surface_Rotating(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            unsigned long iter, unsigned short iZone) {

            double deltaT, time_new, time_old, Lref, *Coord;
            double Center[3] = { 0.0, 0.0, 0.0 }, VarCoord[3] = { 0.0, 0.0, 0.0 }, Omega[3] = { 0.0, 0.0, 0.0 },
                rotCoord[3] = { 0.0, 0.0, 0.0 }, r[3] = { 0.0, 0.0, 0.0 }, Center_Aux[3] = { 0.0, 0.0, 0.0 };
            double rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
            double dtheta, dphi, dpsi, cosTheta, sinTheta;
            double cosPhi, sinPhi, cosPsi, sinPsi;
            unsigned short iMarker, jMarker, Moving, iDim, nDim = geometry->GetnDim();
            unsigned long iPoint, iVertex;
            string Marker_Tag, Moving_Tag;
            int rank;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Initialize the delta variation in coordinates ---*/
            VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

            /*--- Retrieve values from the config file ---*/

            deltaT = config->GetDelta_UnstTimeND();
            Lref = config->GetLength_Ref();

            /*--- Compute delta time based on physical time step ---*/
            time_new = static_cast<double>(iter)*deltaT;
            if (iter == 0) {
                time_old = time_new;
            }
            else {
                time_old = static_cast<double>(iter - 1)*deltaT;
            }

            /*--- Store displacement of each node on the rotating surface ---*/
            /*--- Loop over markers and find the particular marker(s) (surface) to rotate ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Moving = config->GetMarker_All_Moving(iMarker);
                if (Moving == TBOX::YES) {
                    for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                        Moving_Tag = config->GetMarker_Moving(jMarker);
                        Marker_Tag = config->GetMarker_All_TagBound(iMarker);

                        if (Marker_Tag == Moving_Tag) {

                            /*--- Rotation origin and angular velocity from config. ---*/

                            Center[0] = config->GetMotion_Origin_X(jMarker);
                            Center[1] = config->GetMotion_Origin_Y(jMarker);
                            Center[2] = config->GetMotion_Origin_Z(jMarker);
                            Omega[0] = config->GetRotation_Rate_X(jMarker) / config->GetOmega_Ref();
                            Omega[1] = config->GetRotation_Rate_Y(jMarker) / config->GetOmega_Ref();
                            Omega[2] = config->GetRotation_Rate_Z(jMarker) / config->GetOmega_Ref();

                            /*--- Print some information to the console. Be verbose at the first
                            iteration only (mostly for debugging purposes). ---*/
                            // Note that the TBOX::MASTER_NODE might not contain all the markers being moved.

                            if (rank == TBOX::MASTER_NODE) {
                                cout << " Storing rotating displacement for marker: ";
                                cout << Marker_Tag << "." << endl;
                                if (iter == 0) {
                                    cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
                                    cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
                                    cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
                                }
                            }

                            /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

                            dtheta = Omega[0] * (time_new - time_old);
                            dphi = Omega[1] * (time_new - time_old);
                            dpsi = Omega[2] * (time_new - time_old);

                            /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

                            cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
                            sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);

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

                            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

                                /*--- Index and coordinates of the current point ---*/

                                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                                Coord = geometry->node[iPoint]->GetCoord();

                                /*--- Calculate non-dim. position from rotation center ---*/

                                for (iDim = 0; iDim < nDim; iDim++)
                                    r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;
                                if (nDim == 2) r[nDim] = 0.0;

                                /*--- Compute transformed point coordinates ---*/

                                rotCoord[0] = rotMatrix[0][0] * r[0]
                                    + rotMatrix[0][1] * r[1]
                                    + rotMatrix[0][2] * r[2] + Center[0];

                                rotCoord[1] = rotMatrix[1][0] * r[0]
                                    + rotMatrix[1][1] * r[1]
                                    + rotMatrix[1][2] * r[2] + Center[1];

                                rotCoord[2] = rotMatrix[2][0] * r[0]
                                    + rotMatrix[2][1] * r[1]
                                    + rotMatrix[2][2] * r[2] + Center[2];

                                /*--- Calculate delta change in the x, y, & z directions ---*/
                                for (iDim = 0; iDim < nDim; iDim++)
                                    VarCoord[iDim] = (rotCoord[iDim] - Coord[iDim]) / Lref;
                                if (nDim == 2) VarCoord[nDim] = 0.0;

                                /*--- Set node displacement for volume deformation ---*/
                                geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                            }
                        }
                    }
                }
            }

            /*--- When updating the origins it is assumed that all markers have the
            same rotation movement, because we use the last markers rotation matrix and center ---*/

            /*--- Set the mesh motion center to the new location after
            incrementing the position with the rotation. This new
            location will be used for subsequent mesh motion for the given marker.---*/

            for (jMarker = 0; jMarker < config->GetnMarker_Moving(); jMarker++) {

                /*-- Check if we want to update the motion origin for the given marker ---*/

                if (config->GetMoveMotion_Origin(jMarker) == TBOX::YES) {

                    Center_Aux[0] = config->GetMotion_Origin_X(jMarker);
                    Center_Aux[1] = config->GetMotion_Origin_Y(jMarker);
                    Center_Aux[2] = config->GetMotion_Origin_Z(jMarker);

                    /*--- Calculate non-dim. position from rotation center ---*/

                    for (iDim = 0; iDim < nDim; iDim++)
                        r[iDim] = (Center_Aux[iDim] - Center[iDim]) / Lref;
                    if (nDim == 2) r[nDim] = 0.0;

                    /*--- Compute transformed point coordinates ---*/

                    rotCoord[0] = rotMatrix[0][0] * r[0]
                        + rotMatrix[0][1] * r[1]
                        + rotMatrix[0][2] * r[2] + Center[0];

                    rotCoord[1] = rotMatrix[1][0] * r[0]
                        + rotMatrix[1][1] * r[1]
                        + rotMatrix[1][2] * r[2] + Center[1];

                    rotCoord[2] = rotMatrix[2][0] * r[0]
                        + rotMatrix[2][1] * r[1]
                        + rotMatrix[2][2] * r[2] + Center[2];

                    /*--- Calculate delta change in the x, y, & z directions ---*/
                    for (iDim = 0; iDim < nDim; iDim++)
                        VarCoord[iDim] = (rotCoord[iDim] - Center_Aux[iDim]) / Lref;
                    if (nDim == 2) VarCoord[nDim] = 0.0;
                    config->SetMotion_Origin_X(jMarker, Center_Aux[0] + VarCoord[0]);
                    config->SetMotion_Origin_Y(jMarker, Center_Aux[1] + VarCoord[1]);
                    config->SetMotion_Origin_Z(jMarker, Center_Aux[2] + VarCoord[2]);
                }
            }

            /*--- Set the moment computation center to the new location after
            incrementing the position with the rotation. ---*/

            for (jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {

                Center_Aux[0] = config->GetRefOriginMoment_X(jMarker);
                Center_Aux[1] = config->GetRefOriginMoment_Y(jMarker);
                Center_Aux[2] = config->GetRefOriginMoment_Z(jMarker);

                /*--- Calculate non-dim. position from rotation center ---*/

                for (iDim = 0; iDim < nDim; iDim++)
                    r[iDim] = (Center_Aux[iDim] - Center[iDim]) / Lref;
                if (nDim == 2) r[nDim] = 0.0;

                /*--- Compute transformed point coordinates ---*/

                rotCoord[0] = rotMatrix[0][0] * r[0]
                    + rotMatrix[0][1] * r[1]
                    + rotMatrix[0][2] * r[2] + Center[0];

                rotCoord[1] = rotMatrix[1][0] * r[0]
                    + rotMatrix[1][1] * r[1]
                    + rotMatrix[1][2] * r[2] + Center[1];

                rotCoord[2] = rotMatrix[2][0] * r[0]
                    + rotMatrix[2][1] * r[1]
                    + rotMatrix[2][2] * r[2] + Center[2];

                /*--- Calculate delta change in the x, y, & z directions ---*/
                for (iDim = 0; iDim < nDim; iDim++)
                    VarCoord[iDim] = (rotCoord[iDim] - Center_Aux[iDim]) / Lref;
                if (nDim == 2) VarCoord[nDim] = 0.0;

                config->SetRefOriginMoment_X(jMarker, Center_Aux[0] + VarCoord[0]);
                config->SetRefOriginMoment_Y(jMarker, Center_Aux[1] + VarCoord[1]);
                config->SetRefOriginMoment_Z(jMarker, Center_Aux[2] + VarCoord[2]);
            }
        }

        void GRID_SurfaceMovement::AeroelasticDeform(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned long ExtIter, unsigned short iMarker, unsigned short iMarker_Monitoring, double displacements[4]) {

            /* The sign conventions of these are those of the Typical Section Wing Model, below the signs are corrected */
            double dh = -displacements[0];           // relative plunge
            double dalpha = -displacements[1];       // relative pitch
            double dh_x, dh_y;
            double Center[2];
            unsigned short iDim;
            double Lref = config->GetLength_Ref();
            double *Coord;
            unsigned long iPoint, iVertex;
            double x_new, y_new;
            double VarCoord[3];
            string Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);

            /*--- Calculate the plunge displacement for the Typical Section Wing Model taking into account rotation ---*/
            if (config->GetKind_GridMovement(TBOX::ZONE_0) == TBOX::AEROELASTIC_RIGID_MOTION) {
                double Omega, dt, psi;
                dt = config->GetDelta_UnstTimeND();
                Omega = (config->GetRotation_Rate_Z(TBOX::ZONE_0) / config->GetOmega_Ref());
                psi = Omega*(dt*ExtIter);

                /* --- Correct for the airfoil starting position (This is hardcoded in here) --- */
                if (Monitoring_Tag == "Airfoil1") {
                    psi = psi + 0.0;
                }
                else if (Monitoring_Tag == "Airfoil2") {
                    psi = psi + 2.0 / 3.0*TBOX::PI_NUMBER;
                }
                else if (Monitoring_Tag == "Airfoil3") {
                    psi = psi + 4.0 / 3.0*TBOX::PI_NUMBER;
                }
                else
                    cout << "WARNING: There is a marker that we are monitoring that doesn't match the values hardcoded above!" << endl;

                dh_x = -dh*sin(psi);
                dh_y = dh*cos(psi);

            }
            else {
                dh_x = 0;
                dh_y = dh;
            }

            /*--- Pitching origin from config. ---*/

            Center[0] = config->GetRefOriginMoment_X(iMarker_Monitoring);
            Center[1] = config->GetRefOriginMoment_Y(iMarker_Monitoring);

            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                /*--- Coordinates of the current point ---*/
                Coord = geometry->node[iPoint]->GetCoord();

                /*--- Calculate non-dim. position from rotation center ---*/
                double r[2] = { 0, 0 };
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
                    r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;

                /*--- Compute delta of transformed point coordinates ---*/
                // The deltas are needed for the FEA grid deformation Method.
                // rotation contribution - previous position + plunging contribution
                x_new = cos(dalpha)*r[0] - sin(dalpha)*r[1] - r[0] + dh_x;
                y_new = sin(dalpha)*r[0] + cos(dalpha)*r[1] - r[1] + dh_y;

                VarCoord[0] = x_new;
                VarCoord[1] = y_new;
                VarCoord[2] = 0.0;

                /*--- Store new delta node locations for the surface ---*/
                geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            }
            /*--- Set the elastic axis to the new location after incrementing the position with the plunge ---*/
            config->SetRefOriginMoment_X(iMarker_Monitoring, Center[0] + dh_x);
            config->SetRefOriginMoment_Y(iMarker_Monitoring, Center[1] + dh_y);


        }

        void GRID_SurfaceMovement::SetBoundary_Flutter3D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            GRID_FreeFormDefBox **FFDBox, unsigned long iter, unsigned short iZone) {

            double omega, deltaT;
            double alpha, alpha_new, alpha_old;
            double time_new, time_old;
            double Center[3], Omega[3], Ampl[3], Phase[3];
            double DEG2RAD = TBOX::PI_NUMBER / 180.0;
            int rank;
            bool adjoint = config->GetAdjoint();

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
            rank = TBOX::MASTER_NODE;
#endif

            /*--- Retrieve values from the config file ---*/

            deltaT = config->GetDelta_UnstTimeND();

            /*--- Pitching origin, frequency, and amplitude from config. ---*/

            Center[0] = config->GetMotion_Origin_X(iZone);
            Center[1] = config->GetMotion_Origin_Y(iZone);
            Center[2] = config->GetMotion_Origin_Z(iZone);
            Omega[0] = (config->GetPitching_Omega_X(iZone) / config->GetOmega_Ref());
            Omega[1] = (config->GetPitching_Omega_Y(iZone) / config->GetOmega_Ref());
            Omega[2] = (config->GetPitching_Omega_Z(iZone) / config->GetOmega_Ref());
            Ampl[0] = config->GetPitching_Ampl_X(iZone)*DEG2RAD;
            Ampl[1] = config->GetPitching_Ampl_Y(iZone)*DEG2RAD;
            Ampl[2] = config->GetPitching_Ampl_Z(iZone)*DEG2RAD;
            Phase[0] = config->GetPitching_Phase_X(iZone)*DEG2RAD;
            Phase[1] = config->GetPitching_Phase_Y(iZone)*DEG2RAD;
            Phase[2] = config->GetPitching_Phase_Z(iZone)*DEG2RAD;

            /*--- Compute delta time based on physical time step ---*/

            if (adjoint) {

                /*--- For the unsteady adjoint, we integrate backwards through
                physical time, so perform mesh motion in reverse. ---*/

                unsigned long nFlowIter = config->GetnExtIter();
                unsigned long directIter = nFlowIter - iter - 1;
                time_new = static_cast<double>(directIter)*deltaT;
                time_old = time_new;
                if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
            }
            else {

                /*--- Forward time for the direct problem ---*/

                time_new = static_cast<double>(iter)*deltaT;
                time_old = time_new;
                if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
            }

            /*--- Update the pitching angle at this time step. Flip sign for
            nose-up positive convention. ---*/

            omega = Omega[2];
            alpha_new = Ampl[2] * sin(omega*time_new);
            alpha_old = Ampl[2] * sin(omega*time_old);
            alpha = (1E-10 + (alpha_new - alpha_old))*(-TBOX::PI_NUMBER / 180.0);

            if (rank == TBOX::MASTER_NODE)
                cout << "New dihedral angle (alpha): " << alpha_new / DEG2RAD << " degrees." << endl;

            unsigned short iOrder, jOrder, kOrder;
            short iFFDBox;
            double movement[3] = { 0.0, 0.0, 0.0 };
            bool *move = new bool[nFFDBox];
            unsigned short *index = new unsigned short[3];

            move[0] = true; move[1] = true; move[2] = true;

            /*--- Change the value of the control point if move is true ---*/

            for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
                if (move[iFFDBox])
                    for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
                        for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
                            for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
                                index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
                                double *coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
                                movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1] * tan(alpha);
                                FFDBox[iFFDBox]->SetControlPoints(index, movement);
                            }

            /*--- Recompute cartesian coordinates using the new control points position ---*/

            for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
                SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);

        }

        void GRID_SurfaceMovement::SetExternal_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter) {

            int rank = TBOX::MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables ---*/

            unsigned short iDim, nDim;
            unsigned long iPoint = 0, flowIter = 0;
            unsigned long jPoint, GlobalIndex;
            double VarCoord[3], *Coord_Old = NULL, *Coord_New = NULL, Center[3] = { 0.0, 0.0, 0.0 };
            double Lref = config->GetLength_Ref();
            double NewCoord[3] = { 0.0, 0.0, 0.0 }, rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
            double r[3] = { 0.0, 0.0, 0.0 }, rotCoord[3] = { 0.0, 0.0, 0.0 };
            unsigned long iVertex;
            unsigned short iMarker;
            char buffer[50];
            string motion_filename, UnstExt, text_line;
            ifstream motion_file;
            bool unsteady = config->GetUnsteady_Simulation();
            bool adjoint = config->GetAdjoint();

            /*--- Load stuff from config ---*/

            nDim = geometry->GetnDim();
            motion_filename = config->GetMotion_FileName();

            /*--- Set the extension for the correct unsteady mesh motion file ---*/

            if (unsteady) {
                if (adjoint) {
                    /*--- For the unsteady adjoint, we integrate backwards through
                    physical time, so perform mesh motion in reverse. ---*/
                    unsigned long nFlowIter = config->GetnExtIter() - 1;
                    flowIter = nFlowIter - iter;
                    unsigned short lastindex = motion_filename.find_last_of(".");
                    motion_filename = motion_filename.substr(0, lastindex);
                    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf(buffer, "_0000%d.dat", int(flowIter));
                    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf(buffer, "_000%d.dat", int(flowIter));
                    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf(buffer, "_00%d.dat", int(flowIter));
                    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf(buffer, "_0%d.dat", int(flowIter));
                    if (int(flowIter) >= 10000) sprintf(buffer, "_%d.dat", int(flowIter));
                    UnstExt = string(buffer);
                    motion_filename.append(UnstExt);
                }
                else {
                    /*--- Forward time for the direct problem ---*/
                    flowIter = iter;
                    unsigned short lastindex = motion_filename.find_last_of(".");
                    motion_filename = motion_filename.substr(0, lastindex);
                    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf(buffer, "_0000%d.dat", int(flowIter));
                    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf(buffer, "_000%d.dat", int(flowIter));
                    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf(buffer, "_00%d.dat", int(flowIter));
                    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf(buffer, "_0%d.dat", int(flowIter));
                    if (int(flowIter) >= 10000) sprintf(buffer, "_%d.dat", int(flowIter));
                    UnstExt = string(buffer);
                    motion_filename.append(UnstExt);
                }

                if (rank == TBOX::MASTER_NODE)
                    cout << "Reading in the arbitrary mesh motion from direct iteration " << flowIter << "." << endl;
            }

            /*--- Open the motion file ---*/

            motion_file.open(motion_filename.data(), ios::in);
            /*--- Throw error if there is no file ---*/
            if (motion_file.fail()) {
                cout << "There is no mesh motion file!" << endl;
                exit(EXIT_FAILURE);
            }

            /*--- Read in and store the new mesh node locations ---*/

            while (getline(motion_file, text_line)) {
                istringstream point_line(text_line);
                if (nDim == 2) point_line >> iPoint >> NewCoord[0] >> NewCoord[1];
                if (nDim == 3) point_line >> iPoint >> NewCoord[0] >> NewCoord[1] >> NewCoord[2];
                for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                    if (config->GetMarker_All_Moving(iMarker) == TBOX::YES) {
                        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                            jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                            GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
                            if (GlobalIndex == iPoint) {
                                geometry->vertex[iMarker][iVertex]->SetVarCoord(NewCoord);
                                break;
                            }
                        }
                    }
                }
            }
            /*--- Close the restart file ---*/
            motion_file.close();

            /*--- If rotating as well, prepare the rotation matrix ---*/

            if (config->GetGrid_Movement() &&
                config->GetKind_GridMovement(iZone) == TBOX::EXTERNAL_ROTATION) {

                /*--- Variables needed only for rotation ---*/

                double Omega[3], dt;
                double dtheta, dphi, dpsi, cosTheta, sinTheta;
                double cosPhi, sinPhi, cosPsi, sinPsi;

                /*--- Center of rotation & angular velocity vector from config ---*/
                Center[0] = config->GetMotion_Origin_X(iZone);
                Center[1] = config->GetMotion_Origin_Y(iZone);
                Center[2] = config->GetMotion_Origin_Z(iZone);

                /*--- Angular velocity vector from config ---*/

                dt = static_cast<double>(iter)*config->GetDelta_UnstTimeND();
                Omega[0] = config->GetRotation_Rate_X(iZone);
                Omega[1] = config->GetRotation_Rate_Y(iZone);
                Omega[2] = config->GetRotation_Rate_Z(iZone);

                /*--- For the unsteady adjoint, use reverse time ---*/
                if (adjoint) {
                    /*--- Set the first adjoint mesh position to the final direct one ---*/
                    if (iter == 0) dt = ((double)config->GetnExtIter() - 1) * dt;
                    /*--- Reverse the rotation direction for the adjoint ---*/
                    else dt = -1.0*dt;
                }
                else {
                    /*--- No rotation at all for the first direct solution ---*/
                    if (iter == 0) dt = 0;
                }

                /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

                dtheta = Omega[0] * dt;
                dphi = Omega[1] * dt;
                dpsi = Omega[2] * dt;

                /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

                cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
                sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);

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

            }

            /*--- Loop through to find only moving surface markers ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (config->GetMarker_All_Moving(iMarker) == TBOX::YES) {

                    /*--- Loop over all surface points for this marker ---*/

                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

                        /*--- Get current and new coordinates from file ---*/

                        Coord_Old = geometry->node[iPoint]->GetCoord();
                        Coord_New = geometry->vertex[iMarker][iVertex]->GetVarCoord();

                        /*--- If we're also rotating, multiply each point by the
                        rotation matrix. It is assumed that the coordinates in
                        Coord_Old have already been rotated using SetRigid_Rotation(). ---*/

                        if (config->GetGrid_Movement() &&
                            config->GetKind_GridMovement(iZone) == TBOX::EXTERNAL_ROTATION) {

                            /*--- Calculate non-dim. position from rotation center ---*/

                            for (iDim = 0; iDim < nDim; iDim++)
                                r[iDim] = (Coord_New[iDim] - Center[iDim]) / Lref;
                            if (nDim == 2) r[nDim] = 0.0;

                            /*--- Compute transformed point coordinates ---*/

                            rotCoord[0] = rotMatrix[0][0] * r[0]
                                + rotMatrix[0][1] * r[1]
                                + rotMatrix[0][2] * r[2] + Center[0];

                            rotCoord[1] = rotMatrix[1][0] * r[0]
                                + rotMatrix[1][1] * r[1]
                                + rotMatrix[1][2] * r[2] + Center[1];

                            rotCoord[2] = rotMatrix[2][0] * r[0]
                                + rotMatrix[2][1] * r[1]
                                + rotMatrix[2][2] * r[2] + Center[2];

                            /*--- Copy rotated coords back to original array for consistency ---*/
                            for (iDim = 0; iDim < nDim; iDim++)
                                Coord_New[iDim] = rotCoord[iDim];
                        }

                        /*--- Calculate delta change in the x, y, & z directions ---*/
                        for (iDim = 0; iDim < nDim; iDim++)
                            VarCoord[iDim] = (Coord_New[iDim] - Coord_Old[iDim]) / Lref;
                        if (nDim == 2) VarCoord[nDim] = 0.0;

                        /*--- Set position changes to be applied by the spring analogy ---*/
                        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                    }
                }
            }
        }

        void GRID_SurfaceMovement::SetNACA_4Digits(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config) {
            unsigned long iVertex;
            unsigned short iMarker;
            double VarCoord[3], *Coord, *Normal, Ycurv, Yesp;

            if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get(); }

            double Ya = config->GetParamDV(0, 0) / 100.0; /*--- Maximum camber as a fraction of the chord
                                                          (100 m is the first of the four digits) ---*/
            double Xa = config->GetParamDV(0, 1) / 10.0; /*--- Location of maximum camber as a fraction of
                                                         the chord (10 p is the second digit in the NACA xxxx description) ---*/
            double t = config->GetParamDV(0, 2) / 100.0; /*--- Maximum thickness as a fraction of the
                                                         chord (so 100 t gives the last two digits in
                                                         the NACA 4-digit denomination) ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                    if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
                        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();

                        if (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0] - pow(Coord[0], 2.0))*(Ya / pow(Xa, 2.0));
                        else Ycurv = ((1.0 - 2.0*Xa) + 2.0*Xa*Coord[0] - pow(Coord[0], 2.0))*(Ya / pow((1.0 - Xa), 2.0));

                        Yesp = t*(1.4845*sqrt(Coord[0]) - 0.6300*Coord[0] - 1.7580*pow(Coord[0], 2.0) +
                            1.4215*pow(Coord[0], 3.0) - 0.518*pow(Coord[0], 4.0));

                        if (Normal[1] > 0) VarCoord[1] = (Ycurv + Yesp) - Coord[1];
                        if (Normal[1] < 0) VarCoord[1] = (Ycurv - Yesp) - Coord[1];

                    }
                    boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
                }
        }

        void GRID_SurfaceMovement::SetParabolic(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config) {
            unsigned long iVertex;
            unsigned short iMarker;
            double VarCoord[3], *Coord, *Normal;

            if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get(); }

            double c = config->GetParamDV(0, 0); /*--- Center of the parabola ---*/
            double t = config->GetParamDV(0, 1) / 100.0; /*--- Thickness of the parabola ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                    if (config->GetMarker_All_DV(iMarker) == TBOX::YES) {
                        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
                        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();

                        if (Normal[1] > 0) {
                            VarCoord[1] = t*(Coord[0] * Coord[0] - Coord[0]) / (2.0*(c*c - c)) - Coord[1];
                        }
                        if (Normal[1] < 0) {
                            VarCoord[1] = t*(Coord[0] - Coord[0] * Coord[0]) / (2.0*(c*c - c)) - Coord[1];
                        }
                    }
                    boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
                }
        }

        void GRID_SurfaceMovement::SetAirfoil(GEOM::GEOM_Geometry *boundary, TBOX::TBOX_Config *config) {
            unsigned long iVertex, n_Airfoil = 0;
            unsigned short iMarker, nUpper, nLower, iUpper, iLower, iVar, iDim;
            double *VarCoord, *Coord, NewYCoord, NewXCoord, *Coord_i, *Coord_ip1, yp1, ypn,
                Airfoil_Coord[2] = { 0.0, 0.0 }, factor, coeff = 10000, Upper, Lower, Arch = 0.0, TotalArch = 0.0,
                x_i, x_ip1, y_i, y_ip1, AirfoilScale;
            vector<double> Svalue, Xcoord, Ycoord, Xcoord2, Ycoord2, Xcoord_Aux, Ycoord_Aux;
            bool AddBegin = true, AddEnd = true;
            char AirfoilFile[256], AirfoilFormat[15], MeshOrientation[15], AirfoilClose[15];
            ifstream airfoil_file;
            string text_line;

            unsigned short nDim = boundary->GetnDim();

            VarCoord = new double[nDim];
            for (iDim = 0; iDim < nDim; iDim++)
                VarCoord[iDim] = 0.0;

            /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
            deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
            meshes after imposing design variable surface deformations (DV_MARKER). ---*/

            unsigned short Kind_SU2 = config->GetKind_SU2();

            /*--- Read the coordinates. Two main formats:
            - Selig are in an x, y format starting from trailing edge, along the upper surface to the leading
            edge and back around the lower surface to trailing edge.
            - Lednicer are upper surface points leading edge to trailing edge and then lower surface leading
            edge to trailing edge.
            ---*/

            /*--- Open the restart file, throw an error if this fails. ---*/

            cout << "Enter the name of file with the airfoil information: ";
            scanf("%s", AirfoilFile);
            airfoil_file.open(AirfoilFile, ios::in);
            if (airfoil_file.fail()) {
                cout << "There is no airfoil file!! " << endl;
                exit(EXIT_FAILURE);
            }
            cout << "Enter the format of the airfoil (Selig or Lednicer): ";
            scanf("%s", AirfoilFormat);

            cout << "Thickness scaling (1.0 means no scaling)?: ";
            scanf("%lf", &AirfoilScale);

            cout << "Close the airfoil (Yes or No)?: ";
            scanf("%s", AirfoilClose);

            cout << "Surface mesh orientation (clockwise, or anticlockwise): ";
            scanf("%s", MeshOrientation);

            /*--- The first line is the header ---*/

            getline(airfoil_file, text_line);
            cout << "File info: " << text_line << endl;

            if (strcmp(AirfoilFormat, "Selig") == 0) {

                while (getline(airfoil_file, text_line)) {
                    istringstream point_line(text_line);

                    /*--- Read the x & y coordinates from this line of the file (anticlockwise) ---*/

                    point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];

                    /*--- Close the arifoil ---*/

                    if (strcmp(AirfoilClose, "Yes") == 0)
                        factor = -atan(coeff*(Airfoil_Coord[0] - 1.0))*2.0 / TBOX::PI_NUMBER;
                    else factor = 1.0;

                    /*--- Store the coordinates in vectors ---*/

                    Xcoord.push_back(Airfoil_Coord[0]);
                    Ycoord.push_back(Airfoil_Coord[1] * factor*AirfoilScale);
                }

            }
            if (strcmp(AirfoilFormat, "Lednicer") == 0) {

                /*--- The second line is the number of points ---*/

                getline(airfoil_file, text_line);
                istringstream point_line(text_line);
                point_line >> Upper >> Lower;

                nUpper = int(Upper);
                nLower = int(Lower);

                Xcoord.resize(nUpper + nLower - 1);
                Ycoord.resize(nUpper + nLower - 1);

                /*--- White line ---*/

                getline(airfoil_file, text_line);

                for (iUpper = 0; iUpper < nUpper; iUpper++) {
                    getline(airfoil_file, text_line);
                    istringstream point_line(text_line);
                    point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
                    Xcoord[nUpper - iUpper - 1] = Airfoil_Coord[0];

                    if (strcmp(AirfoilClose, "Yes") == 0)
                        factor = -atan(coeff*(Airfoil_Coord[0] - 1.0))*2.0 / TBOX::PI_NUMBER;
                    else factor = 1.0;

                    Ycoord[nUpper - iUpper - 1] = Airfoil_Coord[1] * AirfoilScale*factor;
                }

                getline(airfoil_file, text_line);

                for (iLower = 0; iLower < nLower; iLower++) {
                    getline(airfoil_file, text_line);
                    istringstream point_line(text_line);
                    point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];

                    if (strcmp(AirfoilClose, "Yes") == 0)
                        factor = -atan(coeff*(Airfoil_Coord[0] - 1.0))*2.0 / TBOX::PI_NUMBER;
                    else factor = 1.0;

                    Xcoord[nUpper + iLower - 1] = Airfoil_Coord[0];
                    Ycoord[nUpper + iLower - 1] = Airfoil_Coord[1] * AirfoilScale*factor;
                }

            }

            /*--- Check the coordinate (1,0) at the beginning and end of the file ---*/

            if (Xcoord[0] == 1.0) AddBegin = false;
            if (Xcoord[Xcoord.size() - 1] == 1.0) AddEnd = false;

            if (AddBegin) { Xcoord.insert(Xcoord.begin(), 1.0);   Ycoord.insert(Ycoord.begin(), 0.0); }
            if (AddEnd) { Xcoord.push_back(1.0);                Ycoord.push_back(0.0); }

            /*--- Change the orientation (depend on the input file, and the mesh file) ---*/

            if (strcmp(MeshOrientation, "clockwise") == 0) {
                for (iVar = 0; iVar < Xcoord.size(); iVar++) {
                    Xcoord_Aux.push_back(Xcoord[iVar]);
                    Ycoord_Aux.push_back(Ycoord[iVar]);
                }

                for (iVar = 0; iVar < Xcoord.size(); iVar++) {
                    Xcoord[iVar] = Xcoord_Aux[Xcoord.size() - iVar - 1];
                    Ycoord[iVar] = Ycoord_Aux[Xcoord.size() - iVar - 1];
                }
            }

            /*--- Compute the total arch length ---*/

            Arch = 0.0; Svalue.push_back(Arch);

            for (iVar = 0; iVar < Xcoord.size() - 1; iVar++) {
                x_i = Xcoord[iVar];  x_ip1 = Xcoord[iVar + 1];
                y_i = Ycoord[iVar];  y_ip1 = Ycoord[iVar + 1];
                Arch += sqrt((x_ip1 - x_i)*(x_ip1 - x_i) + (y_ip1 - y_i)*(y_ip1 - y_i));
                Svalue.push_back(Arch);
            }
            x_i = Xcoord[Xcoord.size() - 1];  x_ip1 = Xcoord[0];
            y_i = Ycoord[Xcoord.size() - 1];  y_ip1 = Ycoord[0];
            Arch += sqrt((x_ip1 - x_i)*(x_ip1 - x_i) + (y_ip1 - y_i)*(y_ip1 - y_i));

            /*--- Non dimensionalization ---*/

            for (iVar = 0; iVar < Svalue.size(); iVar++) { Svalue[iVar] /= Arch; }

            /*--- Close the restart file ---*/

            airfoil_file.close();

            /*--- Create a spline for X and Y coordiantes using the arch length ---*/

            n_Airfoil = Svalue.size();
            yp1 = (Xcoord[1] - Xcoord[0]) / (Svalue[1] - Svalue[0]);
            ypn = (Xcoord[n_Airfoil - 1] - Xcoord[n_Airfoil - 2]) / (Svalue[n_Airfoil - 1] - Svalue[n_Airfoil - 2]);

            Xcoord2.resize(n_Airfoil + 1);
            boundary->SetSpline(Svalue, Xcoord, n_Airfoil, yp1, ypn, Xcoord2);

            n_Airfoil = Svalue.size();
            yp1 = (Ycoord[1] - Ycoord[0]) / (Svalue[1] - Svalue[0]);
            ypn = (Ycoord[n_Airfoil - 1] - Ycoord[n_Airfoil - 2]) / (Svalue[n_Airfoil - 1] - Svalue[n_Airfoil - 2]);

            Ycoord2.resize(n_Airfoil + 1);
            boundary->SetSpline(Svalue, Ycoord, n_Airfoil, yp1, ypn, Ycoord2);

            TotalArch = 0.0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (((config->GetMarker_All_Moving(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_DEF))) {
                    for (iVertex = 0; iVertex < boundary->nVertex[iMarker] - 1; iVertex++) {
                        Coord_i = boundary->vertex[iMarker][iVertex]->GetCoord();
                        Coord_ip1 = boundary->vertex[iMarker][iVertex + 1]->GetCoord();

                        x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
                        y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];

                        TotalArch += sqrt((x_ip1 - x_i)*(x_ip1 - x_i) + (y_ip1 - y_i)*(y_ip1 - y_i));
                    }
                    Coord_i = boundary->vertex[iMarker][boundary->nVertex[iMarker] - 1]->GetCoord();
                    Coord_ip1 = boundary->vertex[iMarker][0]->GetCoord();
                    x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
                    y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
                    TotalArch += sqrt((x_ip1 - x_i)*(x_ip1 - x_i) + (y_ip1 - y_i)*(y_ip1 - y_i));
                }
            }


            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                Arch = 0.0;
                for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
                    VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
                    if (((config->GetMarker_All_Moving(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_CFD)) ||
                        ((config->GetMarker_All_DV(iMarker) == TBOX::YES) && (Kind_SU2 == TBOX::SU2_DEF))) {
                        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();

                        if (iVertex == 0) Arch = 0.0;
                        else {
                            Coord_i = boundary->vertex[iMarker][iVertex - 1]->GetCoord();
                            Coord_ip1 = boundary->vertex[iMarker][iVertex]->GetCoord();
                            x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
                            y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
                            Arch += sqrt((x_ip1 - x_i)*(x_ip1 - x_i) + (y_ip1 - y_i)*(y_ip1 - y_i)) / TotalArch;
                        }

                        NewXCoord = boundary->GetSpline(Svalue, Xcoord, Xcoord2, n_Airfoil, Arch);
                        NewYCoord = boundary->GetSpline(Svalue, Ycoord, Ycoord2, n_Airfoil, Arch);

                        /*--- Store the delta change in the x & y coordinates ---*/

                        VarCoord[0] = NewXCoord - Coord[0];
                        VarCoord[1] = NewYCoord - Coord[1];
                    }

                    boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

                }
            }

            delete[] VarCoord;

        }

        void GRID_SurfaceMovement::ReadFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox **FFDBox, string val_mesh_filename) {

            string text_line, iTag;
            ifstream mesh_file;
            double coord[3];
            unsigned short degree[3], iFFDBox, iCornerPoints, iControlPoints, iMarker, iDegree, jDegree, kDegree,
                iChar, LevelFFDBox, nParentFFDBox, iParentFFDBox, nChildFFDBox, iChildFFDBox, nMarker, *nCornerPoints,
                *nControlPoints;
            unsigned long iSurfacePoints, iPoint, jPoint, iVertex, nVertex, nPoint, iElem = 0,
                nElem, my_nSurfPoints, nSurfPoints, *nSurfacePoints;

            unsigned short nDim = geometry->GetnDim();
            int rank = TBOX::MASTER_NODE;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            char *cstr = new char[val_mesh_filename.size() + 1];
            strcpy(cstr, val_mesh_filename.c_str());

            mesh_file.open(cstr, ios::in);
            if (mesh_file.fail()) {
                cout << "There is no geometry file (ReadFFDInfo)!!" << endl;
                exit(EXIT_FAILURE);
            }

            while (getline(mesh_file, text_line)) {

                /*--- Read the inner elements ---*/

                string::size_type position = text_line.find("NELEM=", 0);
                if (position != string::npos) {
                    text_line.erase(0, 6); nElem = atoi(text_line.c_str());
                    for (iElem = 0; iElem < nElem; iElem++) {
                        getline(mesh_file, text_line);
                    }
                }

                /*--- Read the inner points ---*/

                position = text_line.find("NPOIN=", 0);
                if (position != string::npos) {
                    text_line.erase(0, 6); nPoint = atoi(text_line.c_str());
                    for (iPoint = 0; iPoint < nPoint; iPoint++) {
                        getline(mesh_file, text_line);
                    }
                }

                /*--- Read the boundaries  ---*/

                position = text_line.find("NMARK=", 0);
                if (position != string::npos) {
                    text_line.erase(0, 6); nMarker = atoi(text_line.c_str());
                    for (iMarker = 0; iMarker < nMarker; iMarker++) {
                        getline(mesh_file, text_line);
                        getline(mesh_file, text_line);
                        text_line.erase(0, 13); nVertex = atoi(text_line.c_str());
                        for (iVertex = 0; iVertex < nVertex; iVertex++) {
                            getline(mesh_file, text_line);
                        }
                    }
                }

                /*--- Read the FFDBox information  ---*/

                position = text_line.find("FFD_NBOX=", 0);
                if (position != string::npos) {
                    text_line.erase(0, 9);
                    nFFDBox = atoi(text_line.c_str());

                    if (rank == TBOX::MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;

                    nCornerPoints = new unsigned short[nFFDBox];
                    nControlPoints = new unsigned short[nFFDBox];
                    nSurfacePoints = new unsigned long[nFFDBox];

                    getline(mesh_file, text_line);
                    text_line.erase(0, 11);
                    nLevel = atoi(text_line.c_str());

                    if (rank == TBOX::MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;

                    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {

                        /*--- Read the name of the FFD box ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 8);

                        /*--- Remove extra data from the FFDBox name ---*/

                        string::size_type position;
                        for (iChar = 0; iChar < 20; iChar++) {
                            position = text_line.find(" ", 0);
                            if (position != string::npos) text_line.erase(position, 1);
                            position = text_line.find("\r", 0);
                            if (position != string::npos) text_line.erase(position, 1);
                            position = text_line.find("\n", 0);
                            if (position != string::npos) text_line.erase(position, 1);
                        }

                        string TagFFDBox = text_line.c_str();

                        if (rank == TBOX::MASTER_NODE) cout << "FFD box tag: " << TagFFDBox << ". ";

                        /*--- Read the level of the FFD box ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 10);
                        LevelFFDBox = atoi(text_line.c_str());

                        if (rank == TBOX::MASTER_NODE) cout << "FFD box level: " << LevelFFDBox << ". ";

                        /*--- Read the degree of the FFD box ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 13); degree[0] = atoi(text_line.c_str());
                        getline(mesh_file, text_line);
                        text_line.erase(0, 13); degree[1] = atoi(text_line.c_str());

                        if (nDim == 2) {
                            degree[2] = 1;
                        }
                        else {
                            getline(mesh_file, text_line);
                            text_line.erase(0, 13); degree[2] = atoi(text_line.c_str());
                        }

                        if (rank == TBOX::MASTER_NODE) {
                            cout << "Degrees: " << degree[0] << ", " << degree[1];
                            if (nDim == 3) cout << ", " << degree[2];
                            cout << ". " << endl;
                        }

                        FFDBox[iFFDBox] = new GRID_FreeFormDefBox(int(degree[0]), int(degree[1]), int(degree[2]));
                        FFDBox[iFFDBox]->SetTag(TagFFDBox); FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

                        /*--- Read the number of parents boxes ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 12);
                        nParentFFDBox = atoi(text_line.c_str());
                        if (rank == TBOX::MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox << ". ";
                        for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
                            getline(mesh_file, text_line);

                            /*--- Remove extra data from the FFDBox name ---*/

                            string::size_type position;
                            for (iChar = 0; iChar < 20; iChar++) {
                                position = text_line.find(" ", 0);
                                if (position != string::npos) text_line.erase(position, 1);
                                position = text_line.find("\r", 0);
                                if (position != string::npos) text_line.erase(position, 1);
                                position = text_line.find("\n", 0);
                                if (position != string::npos) text_line.erase(position, 1);
                            }

                            string ParentFFDBox = text_line.c_str();
                            FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
                        }

                        /*--- Read the number of children boxes ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 13);
                        nChildFFDBox = atoi(text_line.c_str());
                        if (rank == TBOX::MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox << "." << endl;

                        for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
                            getline(mesh_file, text_line);

                            /*--- Remove extra data from the FFDBox name ---*/

                            string::size_type position;
                            for (iChar = 0; iChar < 20; iChar++) {
                                position = text_line.find(" ", 0);
                                if (position != string::npos) text_line.erase(position, 1);
                                position = text_line.find("\r", 0);
                                if (position != string::npos) text_line.erase(position, 1);
                                position = text_line.find("\n", 0);
                                if (position != string::npos) text_line.erase(position, 1);
                            }

                            string ChildFFDBox = text_line.c_str();
                            FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
                        }

                        /*--- Read the number of the corner points ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 18); nCornerPoints[iFFDBox] = atoi(text_line.c_str());
                        if (rank == TBOX::MASTER_NODE) cout << "Corner points: " << nCornerPoints[iFFDBox] << ". ";
                        if (nDim == 2) nCornerPoints[iFFDBox] = nCornerPoints[iFFDBox] * int(2);

                        /*--- Read the coordinates of the corner points ---*/

                        for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {

                            if (nDim == 2) {
                                if (iCornerPoints < nCornerPoints[iFFDBox] / int(2)) {
                                    getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
                                    FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; coord[2] = -0.5;
                                }
                                else {
                                    coord[0] = FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints - nCornerPoints[iFFDBox] / int(2));
                                    coord[1] = FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints - nCornerPoints[iFFDBox] / int(2));
                                    coord[2] = 0.5;
                                }
                            }
                            else {
                                getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
                                FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; FFDBox_line >> coord[2];
                            }

                            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, iCornerPoints);

                        }

                        /*--- Read the number of the control points ---*/

                        getline(mesh_file, text_line);
                        text_line.erase(0, 19); nControlPoints[iFFDBox] = atoi(text_line.c_str());

                        if (rank == TBOX::MASTER_NODE) cout << "Control points: " << nControlPoints[iFFDBox] << ". ";

                        /*--- Method to identify if there is a FFDBox definition ---*/

                        if (nControlPoints[iFFDBox] != 0) FFDBoxDefinition = true;

                        /*--- Read the coordinates of the control points ---*/

                        for (iControlPoints = 0; iControlPoints < nControlPoints[iFFDBox]; iControlPoints++) {
                            getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
                            FFDBox_line >> iDegree; FFDBox_line >> jDegree; FFDBox_line >> kDegree;
                            FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; FFDBox_line >> coord[2];
                            FFDBox[iFFDBox]->SetCoordControlPoints(coord, iDegree, jDegree, kDegree);
                            FFDBox[iFFDBox]->SetCoordControlPoints_Copy(coord, iDegree, jDegree, kDegree);
                        }

                        getline(mesh_file, text_line);
                        text_line.erase(0, 19); nSurfacePoints[iFFDBox] = atoi(text_line.c_str());

                        /*--- The surface points parametric coordinates, all the nodes read the FFD
                        information but they only store their part ---*/

                        my_nSurfPoints = 0;
                        for (iSurfacePoints = 0; iSurfacePoints < nSurfacePoints[iFFDBox]; iSurfacePoints++) {
                            getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
                            FFDBox_line >> iTag; FFDBox_line >> iPoint;

                            if (config->GetMarker_All_TagBound(iTag) != -1) {

                                iMarker = config->GetMarker_All_TagBound(iTag);
                                FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; FFDBox_line >> coord[2];

                                for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                                    jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                                    if (iPoint == geometry->node[jPoint]->GetGlobalIndex()) {
                                        FFDBox[iFFDBox]->Set_MarkerIndex(iMarker);
                                        FFDBox[iFFDBox]->Set_VertexIndex(iVertex);
                                        FFDBox[iFFDBox]->Set_PointIndex(jPoint);
                                        FFDBox[iFFDBox]->Set_ParametricCoord(coord);
                                        FFDBox[iFFDBox]->Set_CartesianCoord(geometry->node[jPoint]->GetCoord());
                                        my_nSurfPoints++;
                                    }
                                }

                            }

                        }

                        nSurfacePoints[iFFDBox] = my_nSurfPoints;

#ifdef HAVE_MPI
                        nSurfPoints = 0;
                        MPI_Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
                        if (rank == TBOX::MASTER_NODE) cout << "Surface points: " << nSurfPoints << "." << endl;
#else
                        nSurfPoints = my_nSurfPoints;
                        if (rank == TBOX::MASTER_NODE) cout << "Surface points: " << nSurfPoints << "." << endl;
#endif

                    }

                    delete[] nCornerPoints;
                    delete[] nControlPoints;
                    delete[] nSurfacePoints;
                }
            }
            mesh_file.close();

            if (nFFDBox == 0) {
                if (rank == TBOX::MASTER_NODE) cout << "There is no FFD box definition. Just in case, check the .su2 file" << endl;
            }

        }

        void GRID_SurfaceMovement::ReadFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, GRID_FreeFormDefBox **FFDBox) {

            string text_line, iTag;
            ifstream mesh_file;
            double coord[3];
            unsigned short degree[3], iFFDBox, iCornerPoints, LevelFFDBox, nParentFFDBox,
                iParentFFDBox, nChildFFDBox, iChildFFDBox, *nCornerPoints;

            unsigned short nDim = geometry->GetnDim();
            int rank = TBOX::MASTER_NODE;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


            /*--- Read the FFDBox information from the config file ---*/

            nFFDBox = config->GetnFFDBox();

            if (rank == TBOX::MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;

            nCornerPoints = new unsigned short[nFFDBox];

            nLevel = 1; // Nested FFD is not active

            if (rank == TBOX::MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;

            for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {

                /*--- Read the name of the FFD box ---*/

                string TagFFDBox = config->GetTagFFDBox(iFFDBox);

                if (rank == TBOX::MASTER_NODE) cout << "FFD box tag: " << TagFFDBox << ". ";

                /*--- Read the level of the FFD box ---*/

                LevelFFDBox = 0; // Nested FFD is not active

                if (rank == TBOX::MASTER_NODE) cout << "FFD box level: " << LevelFFDBox << ". ";

                /*--- Read the degree of the FFD box ---*/

                degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
                degree[1] = config->GetDegreeFFDBox(iFFDBox, 1);

                if (nDim == 2) { degree[2] = 1; }
                else { degree[2] = config->GetDegreeFFDBox(iFFDBox, 2); }

                if (rank == TBOX::MASTER_NODE) {
                    cout << "Degrees: " << degree[0] << ", " << degree[1];
                    if (nDim == 3) cout << ", " << degree[2];
                    cout << ". " << endl;
                }

                FFDBox[iFFDBox] = new GRID_FreeFormDefBox(int(degree[0]), int(degree[1]), int(degree[2]));
                FFDBox[iFFDBox]->SetTag(TagFFDBox); FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

                /*--- Read the number of parents boxes ---*/

                nParentFFDBox = 0; // Nested FFD is not active
                if (rank == TBOX::MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox << ". ";

                for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
                    string ParentFFDBox = "NONE"; // Nested FFD is not active
                    FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
                }

                /*--- Read the number of children boxes ---*/

                nChildFFDBox = 0; // Nested FFD is not active
                if (rank == TBOX::MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox << "." << endl;

                for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
                    string ChildFFDBox = "NONE"; // Nested FFD is not active
                    FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
                }

                /*--- Read the number of the corner points ---*/

                nCornerPoints[iFFDBox] = 8;

                /*--- Read the coordinates of the corner points ---*/

                for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {

                    if (nDim == 2) {
                        if (iCornerPoints < nCornerPoints[iFFDBox] / int(2)) {
                            coord[0] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3);
                            coord[1] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3 + 1);
                            coord[2] = -0.5;
                        }
                        else {
                            coord[0] = FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints - nCornerPoints[iFFDBox] / int(2));
                            coord[1] = FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints - nCornerPoints[iFFDBox] / int(2));
                            coord[2] = 0.5;
                        }
                    }
                    else {
                        coord[0] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3);
                        coord[1] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3 + 1);
                        coord[2] = config->GetCoordFFDBox(iFFDBox, iCornerPoints * 3 + 2);
                    }

                    FFDBox[iFFDBox]->SetCoordCornerPoints(coord, iCornerPoints);

                }

                /*--- Method to identify if there is a FFDBox definition ---*/

                FFDBoxDefinition = false;

            }

            delete[] nCornerPoints;

            if (nFFDBox == 0) {
                if (rank == TBOX::MASTER_NODE) cout << "There is no FFD box definition. Check the config file." << endl;
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD, 1);
                MPI_Finalize();
#endif
            }

        }

        void GRID_SurfaceMovement::MergeFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            /*--- Local variables needed on all processors ---*/

            unsigned long iPoint;
            unsigned short iFFDBox;

#ifndef HAVE_MPI

            /*--- In serial, the single process has access to all geometry, so simply
            load the coordinates into the data structure. ---*/

            /*--- Total number of points in each FFD box. ---*/

            for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {

                /*--- Loop over the mesh to collect the coords of the local points. ---*/

                for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {

                    /*--- Retrieve the current parametric coordinates at this node. ---*/

                    GlobalCoordX[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[0]);
                    GlobalCoordY[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[1]);
                    GlobalCoordZ[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[2]);
                    GlobalPoint[iFFDBox].push_back(FFDBox[iFFDBox]->Get_PointIndex(iPoint));

                    /*--- Marker of the boundary in the local domain. ---*/

                    unsigned short MarkerIndex = FFDBox[iFFDBox]->Get_MarkerIndex(iPoint);
                    string TagBound = config->GetMarker_All_TagBound(MarkerIndex);

                    /*--- Find the Marker of the boundary in the config file. ---*/

                    unsigned short MarkerIndex_CfgFile = config->GetMarker_CfgFile_TagBound(TagBound);
                    string TagBound_CfgFile = config->GetMarker_CfgFile_TagBound(MarkerIndex_CfgFile);

                    /*--- Set the value of the tag at this node. ---*/

                    GlobalTag[iFFDBox].push_back(TagBound_CfgFile);

                }

            }

#else

            /*--- MPI preprocessing ---*/

            int iProcessor, nProcessor, rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            /*--- Local variables needed for merging the geometry with MPI. ---*/

            unsigned long jPoint, iPointLocal;
            unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
            unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
            unsigned long nBuffer_Scalar = 0;

            if (rank == TBOX::MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];

            for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {

                nLocalPoint = 0;
                for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {

                    iPointLocal = FFDBox[iFFDBox]->Get_PointIndex(iPoint);

                    if (iPointLocal < geometry->GetnPointDomain()) {
                        nLocalPoint++;
                    }

                }
                Buffer_Send_nPoint[0] = nLocalPoint;

                /*--- Communicate the total number of nodes on this domain. ---*/

                MPI_Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG,
                    Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, TBOX::MASTER_NODE, MPI_COMM_WORLD);
                MPI_Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

                nBuffer_Scalar = MaxLocalPoint;

                /*--- Send and Recv buffers. ---*/

                double *Buffer_Send_X = new double[MaxLocalPoint];
                double *Buffer_Recv_X = NULL;

                double *Buffer_Send_Y = new double[MaxLocalPoint];
                double *Buffer_Recv_Y = NULL;

                double *Buffer_Send_Z = new double[MaxLocalPoint];
                double *Buffer_Recv_Z = NULL;

                unsigned long *Buffer_Send_Point = new unsigned long[MaxLocalPoint];
                unsigned long *Buffer_Recv_Point = NULL;

                unsigned short *Buffer_Send_MarkerIndex_CfgFile = new unsigned short[MaxLocalPoint];
                unsigned short *Buffer_Recv_MarkerIndex_CfgFile = NULL;

                /*--- Prepare the receive buffers in the master node only. ---*/

                if (rank == TBOX::MASTER_NODE) {

                    Buffer_Recv_X = new double[nProcessor*MaxLocalPoint];
                    Buffer_Recv_Y = new double[nProcessor*MaxLocalPoint];
                    Buffer_Recv_Z = new double[nProcessor*MaxLocalPoint];
                    Buffer_Recv_Point = new unsigned long[nProcessor*MaxLocalPoint];
                    Buffer_Recv_MarkerIndex_CfgFile = new unsigned short[nProcessor*MaxLocalPoint];

                }

                /*--- Main communication routine. Loop over each coordinate and perform
                the MPI comm. Temporary 1-D buffers are used to send the coordinates at
                all nodes on each partition to the master node. These are then unpacked
                by the master and sorted by global index in one large n-dim. array. ---*/

                /*--- Loop over this partition to collect the coords of the local points. ---*/

                jPoint = 0;
                for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {

                    iPointLocal = FFDBox[iFFDBox]->Get_PointIndex(iPoint);

                    if (iPointLocal < geometry->GetnPointDomain()) {

                        /*--- Load local coords into the temporary send buffer. ---*/

                        Buffer_Send_X[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[0];
                        Buffer_Send_Y[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[1];
                        Buffer_Send_Z[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[2];

                        /*--- Store the global index for this local node. ---*/

                        Buffer_Send_Point[jPoint] = geometry->node[FFDBox[iFFDBox]->Get_PointIndex(iPoint)]->GetGlobalIndex();

                        /*--- Marker of the boundary in the local domain. ---*/

                        unsigned short MarkerIndex = FFDBox[iFFDBox]->Get_MarkerIndex(iPoint);
                        string TagBound = config->GetMarker_All_TagBound(MarkerIndex);

                        /*--- Find the Marker of the boundary in the config file.---*/

                        unsigned short MarkerIndex_CfgFile = config->GetMarker_CfgFile_TagBound(TagBound);
                        Buffer_Send_MarkerIndex_CfgFile[jPoint] = MarkerIndex_CfgFile;

                        jPoint++;

                    }

                }

                /*--- Gather the coordinate data on the master node using MPI. ---*/

                MPI_Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, TBOX::MASTER_NODE, MPI_COMM_WORLD);
                MPI_Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, TBOX::MASTER_NODE, MPI_COMM_WORLD);
                MPI_Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, TBOX::MASTER_NODE, MPI_COMM_WORLD);
                MPI_Gather(Buffer_Send_Point, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Point, nBuffer_Scalar, MPI_UNSIGNED_LONG, TBOX::MASTER_NODE, MPI_COMM_WORLD);
                MPI_Gather(Buffer_Send_MarkerIndex_CfgFile, nBuffer_Scalar, MPI_UNSIGNED_SHORT, Buffer_Recv_MarkerIndex_CfgFile, nBuffer_Scalar, MPI_UNSIGNED_SHORT, TBOX::MASTER_NODE, MPI_COMM_WORLD);

                /*--- The master node unpacks and sorts this variable by global index ---*/

                if (rank == TBOX::MASTER_NODE) {

                    jPoint = 0;

                    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

                            /*--- Get global index, then loop over each variable and store ---*/

                            GlobalCoordX[iFFDBox].push_back(Buffer_Recv_X[jPoint]);
                            GlobalCoordY[iFFDBox].push_back(Buffer_Recv_Y[jPoint]);
                            GlobalCoordZ[iFFDBox].push_back(Buffer_Recv_Z[jPoint]);
                            GlobalPoint[iFFDBox].push_back(Buffer_Recv_Point[jPoint]);

                            string TagBound_CfgFile = config->GetMarker_CfgFile_TagBound(Buffer_Recv_MarkerIndex_CfgFile[jPoint]);
                            GlobalTag[iFFDBox].push_back(TagBound_CfgFile);
                            jPoint++;

                        }

                        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

                        jPoint = (iProcessor + 1)*nBuffer_Scalar;

                    }
                }

                /*--- Immediately release the temporary data buffers. ---*/

                delete[] Buffer_Send_X;
                delete[] Buffer_Send_Y;
                delete[] Buffer_Send_Z;
                delete[] Buffer_Send_Point;
                delete[] Buffer_Send_MarkerIndex_CfgFile;

                if (rank == TBOX::MASTER_NODE) {
                    delete[] Buffer_Recv_X;
                    delete[] Buffer_Recv_Y;
                    delete[] Buffer_Recv_Z;
                    delete[] Buffer_Recv_Point;
                    delete[] Buffer_Recv_MarkerIndex_CfgFile;
                }

            }

            if (rank == TBOX::MASTER_NODE) {
                delete[] Buffer_Recv_nPoint;
            }

#endif

        }

        void GRID_SurfaceMovement::WriteFFDInfo(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {


            unsigned short iOrder, jOrder, kOrder, iFFDBox, iCornerPoints, iParentFFDBox, iChildFFDBox;
            unsigned long iSurfacePoints;
            char cstr[TBOX::MAX_STRING_SIZE], mesh_file[TBOX::MAX_STRING_SIZE];
            string str;
            ofstream output_file;
            double *coord;
            string text_line;

            int rank = TBOX::MASTER_NODE;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            unsigned short nDim = geometry->GetnDim();

            /*--- Merge the FFD info ---*/

            MergeFFDInfo(geometry, config);

            /*--- Attach to the mesh file the FFD information ---*/

            if (rank == TBOX::MASTER_NODE) {

                /*--- Read the name of the output file ---*/

                str = config->GetMesh_Out_FileName();
                strcpy(mesh_file, str.c_str());
                strcpy(cstr, mesh_file);

                output_file.precision(15);
                output_file.open(cstr, ios::out | ios::app);

                if (nFFDBox != 0) {
                    output_file << "FFD_NBOX= " << nFFDBox << endl;
                    output_file << "FFD_NLEVEL= " << nLevel << endl;
                }

                for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {

                    output_file << "FFD_TAG= " << FFDBox[iFFDBox]->GetTag() << endl;
                    output_file << "FFD_LEVEL= " << FFDBox[iFFDBox]->GetLevel() << endl;

                    output_file << "FFD_DEGREE_I= " << FFDBox[iFFDBox]->GetlOrder() - 1 << endl;
                    output_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetmOrder() - 1 << endl;
                    if (nDim == 3) output_file << "FFD_DEGREE_K= " << FFDBox[iFFDBox]->GetnOrder() - 1 << endl;

                    output_file << "FFD_PARENTS= " << FFDBox[iFFDBox]->GetnParentFFDBox() << endl;
                    for (iParentFFDBox = 0; iParentFFDBox < FFDBox[iFFDBox]->GetnParentFFDBox(); iParentFFDBox++)
                        output_file << FFDBox[iFFDBox]->GetParentFFDBoxTag(iParentFFDBox) << endl;
                    output_file << "FFD_CHILDREN= " << FFDBox[iFFDBox]->GetnChildFFDBox() << endl;
                    for (iChildFFDBox = 0; iChildFFDBox < FFDBox[iFFDBox]->GetnChildFFDBox(); iChildFFDBox++)
                        output_file << FFDBox[iFFDBox]->GetChildFFDBoxTag(iChildFFDBox) << endl;

                    if (nDim == 2) {
                        output_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() / int(2) << endl;
                        for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints() / int(2); iCornerPoints++) {
                            coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
                            output_file << coord[0] << "\t" << coord[1] << endl;
                        }
                    }
                    else {
                        output_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() << endl;
                        for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints(); iCornerPoints++) {
                            coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
                            output_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
                        }
                    }

                    /*--- Writing control points ---*/

                    if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
                        output_file << "FFD_CONTROL_POINTS= 0" << endl;
                    }
                    else {
                        output_file << "FFD_CONTROL_POINTS= " << FFDBox[iFFDBox]->GetnControlPoints() << endl;
                        for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
                            for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
                                for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
                                    coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
                                    output_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
                                }
                    }

                    /*--- Writing surface points ---*/

                    if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
                        output_file << "FFD_SURFACE_POINTS= 0" << endl;
                    }
                    else {
                        output_file << "FFD_SURFACE_POINTS= " << GlobalTag[iFFDBox].size() << endl;

                        for (iSurfacePoints = 0; iSurfacePoints < GlobalTag[iFFDBox].size(); iSurfacePoints++) {
                            output_file << scientific << GlobalTag[iFFDBox][iSurfacePoints] << "\t" << GlobalPoint[iFFDBox][iSurfacePoints]
                                << "\t" << GlobalCoordX[iFFDBox][iSurfacePoints] << "\t" << GlobalCoordY[iFFDBox][iSurfacePoints]
                                << "\t" << GlobalCoordZ[iFFDBox][iSurfacePoints] << endl;
                        }

                    }

                }

                output_file.close();

            }

        }
    }
}