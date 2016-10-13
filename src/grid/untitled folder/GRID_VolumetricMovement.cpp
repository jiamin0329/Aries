

#include "GRID_VolumetricMovement.hpp"
#include "../Math/Math_Matrix.hpp"
#include "../Math/Math_Vector.hpp"
#include "../Math/Math_LinearSolver.hpp"

using namespace std;
using namespace ARIES::TBOX;

namespace ARIES
{
    namespace GRID
    {
        GRID_VolumetricMovement::GRID_VolumetricMovement(GEOM::GEOM_Geometry *geometry) : GRID_Gridmovement() {

            nDim = geometry->GetnDim();

        }

        GRID_VolumetricMovement::~GRID_VolumetricMovement(void) {

        }


        void GRID_VolumetricMovement::UpdateGridCoord(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            unsigned short iDim;
            unsigned long iPoint, total_index;
            double new_coord;

            /*--- Update the grid coordinates using the solution of the linear system
            after grid deformation (LinSysSol contains the x, y, z displacements). ---*/

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                for (iDim = 0; iDim < nDim; iDim++) {
                    total_index = iPoint*nDim + iDim;
                    new_coord = geometry->node[iPoint]->GetCoord(iDim) + LinSysSol[total_index];
                    if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
                    geometry->node[iPoint]->SetCoord(iDim, new_coord);
                }

        }

        void GRID_VolumetricMovement::UpdateDualGrid(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            /*--- After moving all nodes, update the dual mesh. Recompute the edges and
            dual mesh control volumes in the domain and on the boundaries. ---*/

            geometry->SetCG();
            geometry->SetControlVolume(config, UPDATE);
            geometry->SetBoundControlVolume(config, UPDATE);

        }

        void GRID_VolumetricMovement::UpdateMultiGrid(GEOM::GEOM_Geometry **geometry, TBOX::TBOX_Config *config) {

            unsigned short iMGfine, iMGlevel, nMGlevel = config->GetnMGLevels();

            /*--- Update the multigrid structure after moving the finest grid,
            including computing the grid velocities on the coarser levels. ---*/

            for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
                iMGfine = iMGlevel - 1;
                geometry[iMGlevel]->SetControlVolume(config, geometry[iMGfine], UPDATE);
                geometry[iMGlevel]->SetBoundControlVolume(config, geometry[iMGfine], UPDATE);
                geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
                if (config->GetGrid_Movement())
                    geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine], config);
            }

        }

        void GRID_VolumetricMovement::SetVolume_Deformation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, bool UpdateGeo) {

            unsigned long IterLinSol = 0, Smoothing_Iter, iNonlinear_Iter, MaxIter = 0, RestartIter = 50, Tot_Iter = 0;
            double MinVolume, NumError, Tol_Factor, Residual = 0.0, Residual_Init = 0.0;
            bool Screen_Output;

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Retrieve number or iterations, tol, output, etc. from config ---*/

            Smoothing_Iter = config->GetGridDef_Linear_Iter();
            Screen_Output = config->GetDeform_Output();
            Tol_Factor = config->GetDeform_Tol_Factor();

            /*--- Disable the screen output if we're running SU2_CFD ---*/

            if (config->GetKind_SU2() == SU2_CFD) Screen_Output = false;

            /*--- Initialize the number of spatial dimensions, length of the state
            vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

            nDim = geometry->GetnDim();
            nVar = geometry->GetnDim();
            nPoint = geometry->GetnPoint();
            nPointDomain = geometry->GetnPointDomain();

            /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/

            config->SetKind_Linear_Solver_Prec(LU_SGS);
            LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
            LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
            StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

            /*--- Loop over the total number of grid deformation iterations. The surface
            deformation can be divided into increments to help with stability. In
            particular, the linear elasticity equations hold only for small deformations. ---*/

            for (iNonlinear_Iter = 0; iNonlinear_Iter < config->GetGridDef_Nonlinear_Iter(); iNonlinear_Iter++) {

                /*--- Initialize vector and sparse matrix ---*/

                LinSysSol.SetValZero();
                LinSysRes.SetValZero();
                StiffMatrix.SetValZero();

                /*--- Compute the stiffness matrix entries for all nodes/elements in the
                mesh. FEA uses a finite element method discretization of the linear
                elasticity equations (transfers element stiffnesses to point-to-point). ---*/

                MinVolume = SetFEAMethodContributions_Elem(geometry, config);

                /*--- Compute the tolerance of the linear solver using MinLength ---*/

                NumError = MinVolume * Tol_Factor;

                /*--- Set the boundary displacements (as prescribed by the design variable
                perturbations controlling the surface shape) as a Dirichlet BC. ---*/

                SetBoundaryDisplacements(geometry, config);

                /*--- Fix the location of any points in the domain, if requested. ---*/

                if (config->GetHold_GridFixed())
                    SetDomainDisplacements(geometry, config);

                /*--- Communicate any prescribed boundary displacements via MPI,
                so that all nodes have the same solution and r.h.s. entries
                across all partitions. ---*/

                StiffMatrix.SendReceive_Solution(LinSysSol, geometry, config);
                StiffMatrix.SendReceive_Solution(LinSysRes, geometry, config);

                /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/

                MATH::MATH_MatrixVectorProduct* mat_vec = new MATH::MATH_Matrix_MatrixVectorProduct(StiffMatrix, geometry, config);
                MATH::MATH_Preconditioner* precond = new MATH::MATH_LUSGSPreconditioner(StiffMatrix, geometry, config);
                MATH::MATH_LinearSolver *system = new MATH::MATH_LinearSolver();

                switch (config->GetDeform_Linear_Solver()) {

                    /*--- Solve the linear system (GMRES with restart) ---*/

                case RESTARTED_FGMRES:

                    Tot_Iter = 0; MaxIter = RestartIter;

                    system->FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, 1, &Residual_Init, false);

                    if ((rank == MASTER_NODE) && Screen_Output) {
                        cout << "\n# FGMRES (with restart) residual history" << endl;
                        cout << "# Residual tolerance target = " << NumError << endl;
                        cout << "# Initial residual norm     = " << Residual_Init << endl;
                    }

                    if (rank == MASTER_NODE) { cout << "     " << Tot_Iter << "     " << Residual_Init / Residual_Init << endl; }

                    while (Tot_Iter < Smoothing_Iter) {

                        if (IterLinSol + RestartIter > Smoothing_Iter)
                            MaxIter = Smoothing_Iter - IterLinSol;

                        IterLinSol = system->FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, MaxIter, &Residual, false);
                        Tot_Iter += IterLinSol;

                        if ((rank == MASTER_NODE) && Screen_Output) { cout << "     " << Tot_Iter << "     " << Residual / Residual_Init << endl; }

                        if (Residual < Residual_Init*NumError) { break; }

                    }

                    if ((rank == MASTER_NODE) && Screen_Output) {
                        cout << "# FGMRES (with restart) final (true) residual:" << endl;
                        cout << "# Iteration = " << Tot_Iter << ": |res|/|res0| = " << Residual / Residual_Init << ".\n" << endl;
                    }

                    break;

                    /*--- Solve the linear system (GMRES) ---*/

                case FGMRES:

                    Tot_Iter = system->FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, Smoothing_Iter, &Residual, Screen_Output);

                    break;

                    /*--- Solve the linear system (BCGSTAB) ---*/

                case BCGSTAB:

                    Tot_Iter = system->BCGSTAB_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, Smoothing_Iter, &Residual, Screen_Output);

                    break;

                }

                /*--- Deallocate memory needed by the Krylov linear solver ---*/

                delete system;
                delete mat_vec;
                delete precond;

                /*--- Update the grid coordinates and cell volumes using the solution
                of the linear system (usol contains the x, y, z displacements). ---*/

                UpdateGridCoord(geometry, config);
                if (UpdateGeo)
                    UpdateDualGrid(geometry, config);

                /*--- Check for failed deformation (negative volumes). ---*/

                MinVolume = Check_Grid(geometry);

                if (rank == MASTER_NODE) {
                    cout << "Non-linear iter.: " << iNonlinear_Iter + 1 << "/" << config->GetGridDef_Nonlinear_Iter()
                        << ". Linear iter.: " << Tot_Iter << ". ";
                    if (nDim == 2) cout << "Min. area: " << MinVolume << ". Error: " << Residual << "." << endl;
                    else cout << "Min. volume: " << MinVolume << ". Error: " << Residual << "." << endl;
                }

            }

            /*--- Deallocate vectors for the linear system. ---*/

            LinSysSol.~MATH_Vector();
            LinSysRes.~MATH_Vector();
            StiffMatrix.~MATH_Matrix();

        }

        double GRID_VolumetricMovement::Check_Grid(GEOM::GEOM_Geometry *geometry) {

            unsigned long iElem, ElemCounter = 0, PointCorners[8];
            double Area = 0.0, Volume = 0.0, MaxArea = -1E22, MaxVolume = -1E22, MinArea = 1E22, MinVolume = 1E22, CoordCorners[8][3];
            unsigned short nNodes = 0, iNodes, iDim;
            bool RightVol = true;

            int rank = MASTER_NODE;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Load up each triangle and tetrahedron to check for negative volumes. ---*/

            for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

                if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
                if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;
                if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
                if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
                if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        nNodes = 6;
                if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;

                for (iNodes = 0; iNodes < nNodes; iNodes++) {
                    PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
                    for (iDim = 0; iDim < nDim; iDim++) {
                        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
                    }
                }

                /*--- Triangles ---*/

                if (nDim == 2) {

                    if (nNodes == 3) Area = GetTriangle_Area(CoordCorners);
                    if (nNodes == 4) Area = GetRectangle_Area(CoordCorners);

                    if (Area >= -EPS) RightVol = true;
                    else RightVol = false;;

                    MaxArea = max(MaxArea, Area);
                    MinArea = min(MinArea, Area);

                }

                /*--- Tetrahedra ---*/

                if (nDim == 3) {

                    if (nNodes == 4) Volume = GetTetra_Volume(CoordCorners);
                    if (nNodes == 5) Volume = GetPyram_Volume(CoordCorners);
                    if (nNodes == 6) Volume = GetPrism_Volume(CoordCorners);
                    if (nNodes == 8) Volume = GetHexa_Volume(CoordCorners);

                    if (Volume >= -EPS) RightVol = true;
                    else RightVol = false;;

                    MaxVolume = max(MaxVolume, Volume);
                    MinVolume = min(MinVolume, Volume);

                }

                if (!RightVol) ElemCounter++;

            }

#ifdef HAVE_MPI
            unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
            double MaxVolume_Local = MaxVolume; MaxVolume = 0.0;
            double MinVolume_Local = MinVolume; MinVolume = 0.0;
            MPI_Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

            if ((ElemCounter != 0) && (rank == MASTER_NODE))
                cout << "There are " << ElemCounter << " elements with negative volume.\n" << endl;

            if (nDim == 2) return MinArea;
            else return MinVolume;

        }

        void GRID_VolumetricMovement::ComputeDeforming_Wall_Distance(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            double *coord, dist2, dist;
            unsigned short iDim, iMarker;
            unsigned long iPoint, iVertex, nVertex_SolidWall;

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
            if (rank == MASTER_NODE)
                cout << "Computing distances to the nearest deforming surface." << endl;

            /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
            deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
            meshes after imposing design variable surface deformations (DV_MARKER). ---*/

            unsigned short Kind_SU2 = config->GetKind_SU2();

#ifndef HAVE_MPI

            /*--- Compute the total number of nodes on deforming boundaries ---*/

            nVertex_SolidWall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF)))
                    nVertex_SolidWall += geometry->GetnVertex(iMarker);

            /*--- Allocate an array to hold boundary node coordinates ---*/

            double **Coord_bound;
            Coord_bound = new double*[nVertex_SolidWall];
            for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
                Coord_bound[iVertex] = new double[nDim];

            /*--- Retrieve and store the coordinates of the deforming boundary nodes ---*/

            nVertex_SolidWall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF)))
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Coord_bound[nVertex_SolidWall][iDim] = geometry->node[iPoint]->GetCoord(iDim);
                        nVertex_SolidWall++;
                    }
            }

            /*--- Loop over all interior mesh nodes and compute the distances to each
            of the deforming boundary nodes. Store the minimum distance to the wall for
            each interior mesh node. Store the global minimum distance. ---*/

            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
                coord = geometry->node[iPoint]->GetCoord();
                dist = 1E20;
                for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++) {
                    dist2 = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++)
                        dist2 += (coord[iDim] - Coord_bound[iVertex][iDim])
                        *(coord[iDim] - Coord_bound[iVertex][iDim]);
                    if (dist2 < dist) dist = dist2;
                }
                geometry->node[iPoint]->SetWall_Distance(sqrt(dist));
            }

            /*--- Deallocate the vector of boundary coordinates. ---*/

            for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
                delete[] Coord_bound[iVertex];
            delete[] Coord_bound;


#else

            /*--- Variables and buffers needed for MPI ---*/

            int iProcessor, nProcessor;
            double local_min_dist;
            MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

            unsigned long nLocalVertex_NS = 0, nGlobalVertex_NS = 0, MaxLocalVertex_NS = 0;
            unsigned long *Buffer_Send_nVertex = new unsigned long[1];
            unsigned long *Buffer_Receive_nVertex = new unsigned long[nProcessor];

            /*--- Count the total number of nodes on deforming boundaries within the
            local partition. ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF)))
                    nLocalVertex_NS += geometry->GetnVertex(iMarker);

            /*--- Communicate to all processors the total number of deforming boundary
            nodes, the maximum number of deforming boundary nodes on any single
            partition, and the number of deforming nodes on each partition. ---*/

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

            /*--- Retrieve and store the coordinates of the deforming boundary nodes on
            the local partition and broadcast them to all partitions. ---*/

            nVertex_SolidWall = 0;
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF)))
                    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Buffer_Send_Coord[nVertex_SolidWall*nDim + iDim] = geometry->node[iPoint]->GetCoord(iDim);
                        nVertex_SolidWall++;
                    }

            MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);

            /*--- Loop over all interior mesh nodes on the local partition and compute
            the distances to each of the deforming boundary nodes in the entire mesh.
            Store the minimum distance to the wall for each interior mesh node. ---*/

            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
                coord = geometry->node[iPoint]->GetCoord();
                dist = 1E20;
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
                    for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
                        dist2 = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++)
                            dist2 += (coord[iDim] - Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NS + iVertex)*nDim + iDim])*
                            (coord[iDim] - Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NS + iVertex)*nDim + iDim]);
                        if (dist2 < dist) dist = dist2;
                    }
                geometry->node[iPoint]->SetWall_Distance(sqrt(dist));
            }

            /*--- Deallocate the buffers needed for the MPI communication. ---*/

            delete[] Buffer_Send_Coord;
            delete[] Buffer_Receive_Coord;
            delete[] Buffer_Send_nVertex;
            delete[] Buffer_Receive_nVertex;

#endif

        }

        double GRID_VolumetricMovement::SetFEAMethodContributions_Elem(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            unsigned short iVar, iDim, nNodes = 0, iNodes, StiffMatrix_nElem = 0;
            unsigned long Point_0, Point_1, iElem, iEdge, ElemCounter = 0, PointCorners[8];
            double *Coord_0, *Coord_1, Length, MinLength = 1E10, **StiffMatrix_Elem = NULL, Scale, CoordCorners[8][3];
            double *Edge_Vector = new double[nDim];

            /*--- Allocate maximum size (rectangle and hexahedron) ---*/

            if (nDim == 2) StiffMatrix_nElem = 8;
            else StiffMatrix_nElem = 24;

            StiffMatrix_Elem = new double*[StiffMatrix_nElem];
            for (iVar = 0; iVar < StiffMatrix_nElem; iVar++)
                StiffMatrix_Elem[iVar] = new double[StiffMatrix_nElem];

            /*--- Check the minimum edge length in the entire mesh. ---*/

            for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

                /*--- Points in edge and coordinates ---*/

                Point_0 = geometry->edge[iEdge]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
                Point_1 = geometry->edge[iEdge]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();

                /*--- Compute Edge_Vector ---*/

                Length = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                    Edge_Vector[iDim] = Coord_1[iDim] - Coord_0[iDim];
                    Length += Edge_Vector[iDim] * Edge_Vector[iDim];
                }
                Length = sqrt(Length);
                MinLength = min(Length, MinLength);

            }

            /*--- Compute min volume in the entire mesh. ---*/

            Scale = Check_Grid(geometry);

            /*--- Compute the distance to the nearest deforming surface if needed
            as part of the stiffness calculation. In this case, we can scale based
            on the minimum edge length. ---*/

            if (config->GetDeform_Stiffness_Type() == WALL_DISTANCE) {
                ComputeDeforming_Wall_Distance(geometry, config);
                Scale = MinLength;
            }

            /*--- Compute contributions from each element by forming the stiffness matrix (FEA) ---*/

            for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

                if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
                if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;
                if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
                if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
                if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        nNodes = 6;
                if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;

                for (iNodes = 0; iNodes < nNodes; iNodes++) {
                    PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
                    for (iDim = 0; iDim < nDim; iDim++) {
                        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
                    }
                }

                if (nDim == 2) SetFEA_StiffMatrix2D(geometry, config, StiffMatrix_Elem, PointCorners, CoordCorners, nNodes, Scale);
                if (nDim == 3) SetFEA_StiffMatrix3D(geometry, config, StiffMatrix_Elem, PointCorners, CoordCorners, nNodes, Scale);

                AddFEA_StiffMatrix(geometry, StiffMatrix_Elem, PointCorners, nNodes);

            }

#ifdef HAVE_MPI
            unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
            MPI_Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

            /*--- Deallocate memory and exit ---*/

            for (iVar = 0; iVar < StiffMatrix_nElem; iVar++)
                delete StiffMatrix_Elem[iVar];
            delete[] StiffMatrix_Elem;

            delete[] Edge_Vector;

            /*--- If there are no degenerate cells, use the minimum volume instead ---*/
            if (ElemCounter == 0) MinLength = Scale;

#ifdef HAVE_MPI
            double MinLength_Local = MinLength; MinLength = 0.0;
            MPI_Allreduce(&MinLength_Local, &MinLength, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

            return MinLength;
        }

        double GRID_VolumetricMovement::ShapeFunc_Triangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]) {

            int i, j, k;
            double c0, c1, xsj;
            double xs[3][3], ad[3][3];

            /*--- Shape functions ---*/

            DShapeFunction[0][3] = Xi;
            DShapeFunction[1][3] = Eta;
            DShapeFunction[2][3] = 1 - Xi - Eta;

            /*--- dN/d xi, dN/d eta ---*/

            DShapeFunction[0][0] = 1.0;  DShapeFunction[0][1] = 0.0;
            DShapeFunction[1][0] = 0.0;  DShapeFunction[1][1] = 1.0;
            DShapeFunction[2][0] = -1.0; DShapeFunction[2][1] = -1.0;

            /*--- Jacobian transformation ---*/

            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    xs[i][j] = 0.0;
                    for (k = 0; k < 3; k++) {
                        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
                    }
                }
            }

            /*--- Adjoint to Jacobian ---*/

            ad[0][0] = xs[1][1];
            ad[0][1] = -xs[0][1];
            ad[1][0] = -xs[1][0];
            ad[1][1] = xs[0][0];

            /*--- Determinant of Jacobian ---*/

            xsj = ad[0][0] * ad[1][1] - ad[0][1] * ad[1][0];

            /*--- Jacobian inverse ---*/

            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    xs[i][j] = ad[i][j] / xsj;
                }
            }

            /*--- Derivatives with repect to global coordinates ---*/

            for (k = 0; k < 3; k++) {
                c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1]; // dN/dx
                c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1]; // dN/dy
                DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
                DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
            }

            return xsj;

        }

        double GRID_VolumetricMovement::ShapeFunc_Rectangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]) {

            int i, j, k;
            double c0, c1, xsj;
            double xs[3][3], ad[3][3];

            /*--- Shape functions ---*/

            DShapeFunction[0][3] = 0.25*(1.0 - Xi)*(1.0 - Eta);
            DShapeFunction[1][3] = 0.25*(1.0 + Xi)*(1.0 - Eta);
            DShapeFunction[2][3] = 0.25*(1.0 + Xi)*(1.0 + Eta);
            DShapeFunction[3][3] = 0.25*(1.0 - Xi)*(1.0 + Eta);

            /*--- dN/d xi, dN/d eta ---*/

            DShapeFunction[0][0] = -0.25*(1.0 - Eta); DShapeFunction[0][1] = -0.25*(1.0 - Xi);
            DShapeFunction[1][0] = 0.25*(1.0 - Eta); DShapeFunction[1][1] = -0.25*(1.0 + Xi);
            DShapeFunction[2][0] = 0.25*(1.0 + Eta); DShapeFunction[2][1] = 0.25*(1.0 + Xi);
            DShapeFunction[3][0] = -0.25*(1.0 + Eta); DShapeFunction[3][1] = 0.25*(1.0 - Xi);

            /*--- Jacobian transformation ---*/

            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    xs[i][j] = 0.0;
                    for (k = 0; k < 4; k++) {
                        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
                    }
                }
            }

            /*--- Adjoint to Jacobian ---*/

            ad[0][0] = xs[1][1];
            ad[0][1] = -xs[0][1];
            ad[1][0] = -xs[1][0];
            ad[1][1] = xs[0][0];

            /*--- Determinant of Jacobian ---*/

            xsj = ad[0][0] * ad[1][1] - ad[0][1] * ad[1][0];

            /*--- Jacobian inverse ---*/

            for (i = 0; i < 2; i++) {
                for (j = 0; j < 2; j++) {
                    xs[i][j] = ad[i][j] / xsj;
                }
            }

            /*--- Derivatives with repect to global coordinates ---*/

            for (k = 0; k < 4; k++) {
                c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1]; // dN/dx
                c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1]; // dN/dy
                DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
                DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
            }

            return xsj;

        }

        double GRID_VolumetricMovement::ShapeFunc_Tetra(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {

            int i, j, k;
            double c0, c1, c2, xsj;
            double xs[3][3], ad[3][3];

            /*--- Shape functions ---*/

            DShapeFunction[0][3] = Xi;
            DShapeFunction[1][3] = Zeta;
            DShapeFunction[2][3] = 1.0 - Xi - Eta - Zeta;
            DShapeFunction[3][3] = Eta;

            /*--- dN/d xi, dN/d eta, dN/d zeta ---*/

            DShapeFunction[0][0] = 1.0;  DShapeFunction[0][1] = 0.0;  DShapeFunction[0][2] = 0.0;
            DShapeFunction[1][0] = 0.0;  DShapeFunction[1][1] = 0.0;  DShapeFunction[1][2] = 1.0;
            DShapeFunction[2][0] = -1.0; DShapeFunction[2][1] = -1.0; DShapeFunction[2][2] = -1.0;
            DShapeFunction[3][0] = 0.0;  DShapeFunction[3][1] = 1.0;  DShapeFunction[3][2] = 0.0;

            /*--- Jacobian transformation ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = 0.0;
                    for (k = 0; k < 4; k++) {
                        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
                    }
                }
            }

            /*--- Adjoint to Jacobian ---*/

            ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
            ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
            ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
            ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
            ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
            ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
            ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
            ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
            ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

            /*--- Determinant of Jacobian ---*/

            xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

            /*--- Jacobian inverse ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = ad[i][j] / xsj;
                }
            }

            /*--- Derivatives with repect to global coordinates ---*/

            for (k = 0; k < 4; k++) {
                c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2]; // dN/dx
                c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2]; // dN/dy
                c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2]; // dN/dz
                DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
                DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
                DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
            }

            return xsj;

        }

        double GRID_VolumetricMovement::ShapeFunc_Pyram(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {

            int i, j, k;
            double c0, c1, c2, xsj;
            double xs[3][3], ad[3][3];

            /*--- Shape functions ---*/

            DShapeFunction[0][3] = 0.125*(1.0 - Xi)*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[1][3] = 0.125*(1.0 + Xi)*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[2][3] = 0.125*(1.0 + Xi)*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[3][3] = 0.125*(1.0 - Xi)*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[4][3] = 0.5*(1.0 + Zeta);

            /*--- dN/d xi ---*/

            DShapeFunction[0][0] = -0.125*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[1][0] = 0.125*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[2][0] = 0.125*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[3][0] = -0.125*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[4][0] = 0.0;

            /*--- dN/d eta ---*/

            DShapeFunction[0][1] = -0.125*(1.0 - Xi)*(1.0 - Zeta);
            DShapeFunction[1][1] = -0.125*(1.0 + Xi)*(1.0 - Zeta);
            DShapeFunction[2][1] = 0.125*(1.0 + Xi)*(1.0 - Zeta);
            DShapeFunction[3][1] = 0.125*(1.0 - Xi)*(1.0 - Zeta);
            DShapeFunction[4][1] = 0.0;

            /*--- dN/d zeta ---*/

            DShapeFunction[0][2] = -0.125*(1.0 - Xi)*(1.0 - Eta);
            DShapeFunction[1][2] = -0.125*(1.0 + Xi)*(1.0 - Eta);
            DShapeFunction[2][2] = -0.125*(1.0 + Xi)*(1.0 + Eta);
            DShapeFunction[3][2] = -0.125*(1.0 - Xi)*(1.0 + Eta);
            DShapeFunction[4][2] = 0.5;

            /*--- Jacobian transformation ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = 0.0;
                    for (k = 0; k < 5; k++) {
                        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
                    }
                }
            }

            /*--- Adjoint to Jacobian ---*/

            ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
            ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
            ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
            ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
            ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
            ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
            ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
            ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
            ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

            /*--- Determinant of Jacobian ---*/

            xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

            /*--- Jacobian inverse ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = ad[i][j] / xsj;
                }
            }

            /*--- Derivatives with repect to global coordinates ---*/

            for (k = 0; k < 5; k++) {
                c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2]; // dN/dx
                c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2]; // dN/dy
                c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2]; // dN/dz
                DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
                DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
                DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
            }

            return xsj;

        }

        double GRID_VolumetricMovement::ShapeFunc_Prism(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {

            int i, j, k;
            double c0, c1, c2, xsj;
            double xs[3][3], ad[3][3];

            /*--- Shape functions ---*/

            DShapeFunction[0][3] = 0.5*Eta*(1.0 - Xi);
            DShapeFunction[1][3] = 0.5*Zeta*(1.0 - Xi);
            DShapeFunction[2][3] = 0.5*(1.0 - Eta - Zeta)*(1.0 - Xi);
            DShapeFunction[3][3] = 0.5*Eta*(Xi + 1.0);
            DShapeFunction[4][3] = 0.5*Zeta*(Xi + 1.0);
            DShapeFunction[5][3] = 0.5*(1.0 - Eta - Zeta)*(Xi + 1.0);

            /*--- dN/d Xi, dN/d Eta, dN/d Zeta ---*/

            DShapeFunction[0][0] = -0.5*Eta;            DShapeFunction[0][1] = 0.5*(1.0 - Xi);      DShapeFunction[0][2] = 0.0;
            DShapeFunction[1][0] = -0.5*Zeta;           DShapeFunction[1][1] = 0.0;               DShapeFunction[1][2] = 0.5*(1.0 - Xi);
            DShapeFunction[2][0] = -0.5*(1.0 - Eta - Zeta); DShapeFunction[2][1] = -0.5*(1.0 - Xi);     DShapeFunction[2][2] = -0.5*(1.0 - Xi);
            DShapeFunction[3][0] = 0.5*Eta;             DShapeFunction[3][1] = 0.5*(Xi + 1.0);      DShapeFunction[3][2] = 0.0;
            DShapeFunction[4][0] = 0.5*Zeta;            DShapeFunction[4][1] = 0.0;               DShapeFunction[4][2] = 0.5*(Xi + 1.0);
            DShapeFunction[5][0] = 0.5*(1.0 - Eta - Zeta);  DShapeFunction[5][1] = -0.5*(Xi + 1.0);     DShapeFunction[5][2] = -0.5*(Xi + 1.0);

            /*--- Jacobian transformation ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = 0.0;
                    for (k = 0; k < 6; k++) {
                        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
                    }
                }
            }

            /*--- Adjoint to Jacobian ---*/

            ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
            ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
            ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
            ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
            ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
            ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
            ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
            ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
            ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

            /*--- Determinant of Jacobian ---*/

            xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

            /*--- Jacobian inverse ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = ad[i][j] / xsj;
                }
            }

            /*--- Derivatives with repect to global coordinates ---*/

            for (k = 0; k < 6; k++) {
                c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2]; // dN/dx
                c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2]; // dN/dy
                c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2]; // dN/dz
                DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
                DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
                DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
            }

            return xsj;

        }

        double GRID_VolumetricMovement::ShapeFunc_Hexa(double Xi, double Eta, double Zeta, double CoordCorners[8][3], double DShapeFunction[8][4]) {

            int i, j, k;
            double c0, c1, c2, xsj;
            double xs[3][3], ad[3][3];


            /*--- Shape functions ---*/

            DShapeFunction[0][3] = 0.125*(1.0 - Xi)*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[1][3] = 0.125*(1.0 + Xi)*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[2][3] = 0.125*(1.0 + Xi)*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[3][3] = 0.125*(1.0 - Xi)*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[4][3] = 0.125*(1.0 - Xi)*(1.0 - Eta)*(1.0 + Zeta);
            DShapeFunction[5][3] = 0.125*(1.0 + Xi)*(1.0 - Eta)*(1.0 + Zeta);
            DShapeFunction[6][3] = 0.125*(1.0 + Xi)*(1.0 + Eta)*(1.0 + Zeta);
            DShapeFunction[7][3] = 0.125*(1.0 - Xi)*(1.0 + Eta)*(1.0 + Zeta);

            /*--- dN/d xi ---*/

            DShapeFunction[0][0] = -0.125*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[1][0] = 0.125*(1.0 - Eta)*(1.0 - Zeta);
            DShapeFunction[2][0] = 0.125*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[3][0] = -0.125*(1.0 + Eta)*(1.0 - Zeta);
            DShapeFunction[4][0] = -0.125*(1.0 - Eta)*(1.0 + Zeta);
            DShapeFunction[5][0] = 0.125*(1.0 - Eta)*(1.0 + Zeta);
            DShapeFunction[6][0] = 0.125*(1.0 + Eta)*(1.0 + Zeta);
            DShapeFunction[7][0] = -0.125*(1.0 + Eta)*(1.0 + Zeta);

            /*--- dN/d eta ---*/

            DShapeFunction[0][1] = -0.125*(1.0 - Xi)*(1.0 - Zeta);
            DShapeFunction[1][1] = -0.125*(1.0 + Xi)*(1.0 - Zeta);
            DShapeFunction[2][1] = 0.125*(1.0 + Xi)*(1.0 - Zeta);
            DShapeFunction[3][1] = 0.125*(1.0 - Xi)*(1.0 - Zeta);
            DShapeFunction[4][1] = -0.125*(1.0 - Xi)*(1.0 + Zeta);
            DShapeFunction[5][1] = -0.125*(1.0 + Xi)*(1.0 + Zeta);
            DShapeFunction[6][1] = 0.125*(1.0 + Xi)*(1.0 + Zeta);
            DShapeFunction[7][1] = 0.125*(1.0 - Xi)*(1.0 + Zeta);

            /*--- dN/d zeta ---*/

            DShapeFunction[0][2] = -0.125*(1.0 - Xi)*(1.0 - Eta);
            DShapeFunction[1][2] = -0.125*(1.0 + Xi)*(1.0 - Eta);
            DShapeFunction[2][2] = -0.125*(1.0 + Xi)*(1.0 + Eta);
            DShapeFunction[3][2] = -0.125*(1.0 - Xi)*(1.0 + Eta);
            DShapeFunction[4][2] = 0.125*(1.0 - Xi)*(1.0 - Eta);
            DShapeFunction[5][2] = 0.125*(1.0 + Xi)*(1.0 - Eta);
            DShapeFunction[6][2] = 0.125*(1.0 + Xi)*(1.0 + Eta);
            DShapeFunction[7][2] = 0.125*(1.0 - Xi)*(1.0 + Eta);

            /*--- Jacobian transformation ---*/

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = 0.0;
                    for (k = 0; k < 8; k++) {
                        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
                    }
                }
            }

            /*--- Adjoint to Jacobian ---*/

            ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
            ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
            ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
            ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
            ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
            ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
            ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
            ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
            ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

            /*--- Determinant of Jacobian ---*/

            xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

            /*--- Jacobian inverse ---*/
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    xs[i][j] = ad[i][j] / xsj;
                }
            }

            /*--- Derivatives with repect to global coordinates ---*/

            for (k = 0; k < 8; k++) {
                c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2]; // dN/dx
                c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2]; // dN/dy
                c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2]; // dN/dz
                DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
                DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
                DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
            }

            return xsj;

        }

        double GRID_VolumetricMovement::GetTriangle_Area(double CoordCorners[8][3]) {

            unsigned short iDim;
            double a[3] = { 0.0, 0.0, 0.0 }, b[3] = { 0.0, 0.0, 0.0 };
            double *Coord_0, *Coord_1, *Coord_2, Area;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[1];
            Coord_2 = CoordCorners[2];

            for (iDim = 0; iDim < nDim; iDim++) {
                a[iDim] = Coord_0[iDim] - Coord_2[iDim];
                b[iDim] = Coord_1[iDim] - Coord_2[iDim];
            }

            Area = 0.5*fabs(a[0] * b[1] - a[1] * b[0]);

            return Area;

        }

        double GRID_VolumetricMovement::GetRectangle_Area(double CoordCorners[8][3]) {

            unsigned short iDim;
            double a[3] = { 0.0, 0.0, 0.0 }, b[3] = { 0.0, 0.0, 0.0 };
            double *Coord_0, *Coord_1, *Coord_2, Area;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[1];
            Coord_2 = CoordCorners[2];

            for (iDim = 0; iDim < nDim; iDim++) {
                a[iDim] = Coord_0[iDim] - Coord_2[iDim];
                b[iDim] = Coord_1[iDim] - Coord_2[iDim];
            }

            Area = 0.5*fabs(a[0] * b[1] - a[1] * b[0]);

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[2];
            Coord_2 = CoordCorners[3];

            for (iDim = 0; iDim < nDim; iDim++) {
                a[iDim] = Coord_0[iDim] - Coord_2[iDim];
                b[iDim] = Coord_1[iDim] - Coord_2[iDim];
            }

            Area += 0.5*fabs(a[0] * b[1] - a[1] * b[0]);

            return Area;

        }

        double GRID_VolumetricMovement::GetTetra_Volume(double CoordCorners[8][3]) {

            unsigned short iDim;
            double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
            double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[1];
            Coord_2 = CoordCorners[2];
            Coord_3 = CoordCorners[3];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            return Volume;

        }

        double GRID_VolumetricMovement::GetPyram_Volume(double CoordCorners[8][3]) {

            unsigned short iDim;
            double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
            double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[1];
            Coord_2 = CoordCorners[2];
            Coord_3 = CoordCorners[4];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[2];
            Coord_2 = CoordCorners[3];
            Coord_3 = CoordCorners[4];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            return Volume;

        }

        double GRID_VolumetricMovement::GetPrism_Volume(double CoordCorners[8][3]) {

            unsigned short iDim;
            double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
            double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[2];
            Coord_2 = CoordCorners[1];
            Coord_3 = CoordCorners[5];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[5];
            Coord_2 = CoordCorners[1];
            Coord_3 = CoordCorners[4];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[5];
            Coord_2 = CoordCorners[4];
            Coord_3 = CoordCorners[3];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            return Volume;

        }

        double GRID_VolumetricMovement::GetHexa_Volume(double CoordCorners[8][3]) {

            unsigned short iDim;
            double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
            double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[1];
            Coord_2 = CoordCorners[2];
            Coord_3 = CoordCorners[5];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[2];
            Coord_2 = CoordCorners[7];
            Coord_3 = CoordCorners[5];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[2];
            Coord_2 = CoordCorners[3];
            Coord_3 = CoordCorners[7];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[0];
            Coord_1 = CoordCorners[5];
            Coord_2 = CoordCorners[7];
            Coord_3 = CoordCorners[4];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            Coord_0 = CoordCorners[2];
            Coord_1 = CoordCorners[7];
            Coord_2 = CoordCorners[5];
            Coord_3 = CoordCorners[6];

            for (iDim = 0; iDim < nDim; iDim++) {
                r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
                r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
                r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
            }

            CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
            CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
            CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

            Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

            return Volume;

        }

        void GRID_VolumetricMovement::SetFEA_StiffMatrix2D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, double **StiffMatrix_Elem, unsigned long PointCorners[8], double CoordCorners[8][3], unsigned short nNodes, double scale) {

            double B_Matrix[3][8], D_Matrix[3][3], Aux_Matrix[8][3];
            double Xi = 0.0, Eta = 0.0, Det = 0.0, E, Lambda = 0.0, Nu, Mu = 0.0, Avg_Wall_Dist;
            unsigned short iNode, jNode, iVar, jVar, kVar, iGauss, nGauss = 0;
            double DShapeFunction[8][4] = { { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 } };
            double Location[4][3], Weight[4];
            unsigned short nVar = geometry->GetnDim();

            for (iVar = 0; iVar < nNodes*nVar; iVar++) {
                for (jVar = 0; jVar < nNodes*nVar; jVar++) {
                    StiffMatrix_Elem[iVar][jVar] = 0.0;
                }
            }

            /*--- Each element uses their own stiffness which is inversely
            proportional to the area/volume of the cell. Using Mu = E & Lambda = -E
            is a modification to help allow rigid rotation of elements (see
            "Robust Mesh Deformation using the Linear Elasticity Equations" by
            R. P. Dwight. ---*/

            /*--- Integration formulae from "Shape functions and points of
            integration of the R¨¦sum¨¦" by Josselin DELMAS (2013) ---*/

            /*--- Triangle. Nodes of numerical integration at 1 point (order 1). ---*/

            if (nNodes == 3) {
                nGauss = 1;
                Location[0][0] = 0.333333333333333;  Location[0][1] = 0.333333333333333;  Weight[0] = 0.5;
            }

            /*--- Rectangle. Nodes of numerical integration at 4 points (order 2). ---*/

            if (nNodes == 4) {
                nGauss = 4;
                Location[0][0] = -0.577350269189626;  Location[0][1] = -0.577350269189626;  Weight[0] = 1.0;
                Location[1][0] = 0.577350269189626;   Location[1][1] = -0.577350269189626;  Weight[1] = 1.0;
                Location[2][0] = 0.577350269189626;   Location[2][1] = 0.577350269189626;   Weight[2] = 1.0;
                Location[3][0] = -0.577350269189626;  Location[3][1] = 0.577350269189626;   Weight[3] = 1.0;
            }

            for (iGauss = 0; iGauss < nGauss; iGauss++) {

                Xi = Location[iGauss][0]; Eta = Location[iGauss][1];

                if (nNodes == 3) Det = ShapeFunc_Triangle(Xi, Eta, CoordCorners, DShapeFunction);
                if (nNodes == 4) Det = ShapeFunc_Rectangle(Xi, Eta, CoordCorners, DShapeFunction);

                /*--- Compute the B Matrix ---*/

                for (iVar = 0; iVar < 3; iVar++)
                    for (jVar = 0; jVar < nNodes*nVar; jVar++)
                        B_Matrix[iVar][jVar] = 0.0;

                for (iNode = 0; iNode < nNodes; iNode++) {
                    B_Matrix[0][0 + iNode*nVar] = DShapeFunction[iNode][0];
                    B_Matrix[1][1 + iNode*nVar] = DShapeFunction[iNode][1];

                    B_Matrix[2][0 + iNode*nVar] = DShapeFunction[iNode][1];
                    B_Matrix[2][1 + iNode*nVar] = DShapeFunction[iNode][0];
                }

                /*--- Impose a type of stiffness for each element ---*/

                switch (config->GetDeform_Stiffness_Type()) {

                case INVERSE_VOLUME:
                    E = scale / (Weight[iGauss] * Det);
                    Mu = E;
                    Lambda = -E;
                    break;

                case WALL_DISTANCE:
                    Avg_Wall_Dist = 0.0;
                    for (jNode = 0; jNode < nNodes; jNode++) {
                        Avg_Wall_Dist += geometry->node[PointCorners[jNode]]->GetWall_Distance() / ((double)nNodes);
                    }
                    E = scale / (Weight[iGauss] * Avg_Wall_Dist);
                    Mu = E;
                    Lambda = -E;
                    break;

                case CONSTANT_STIFFNESS:
                    E = config->GetDeform_ElasticityMod();
                    Nu = config->GetDeform_PoissonRatio();
                    //E = 2E11; Nu = 0.30;
                    Mu = E / (2.0*(1.0 + Nu));
                    Lambda = Nu*E / ((1.0 + Nu)*(1.0 - 2.0*Nu));
                    break;
                }

                /*--- Compute the D Matrix (for plane strain and 3-D)---*/

                D_Matrix[0][0] = Lambda + 2.0*Mu;		D_Matrix[0][1] = Lambda;            D_Matrix[0][2] = 0.0;
                D_Matrix[1][0] = Lambda;            D_Matrix[1][1] = Lambda + 2.0*Mu;   D_Matrix[1][2] = 0.0;
                D_Matrix[2][0] = 0.0;               D_Matrix[2][1] = 0.0;               D_Matrix[2][2] = Mu;


                /*--- Compute the BT.D Matrix ---*/

                for (iVar = 0; iVar < nNodes*nVar; iVar++) {
                    for (jVar = 0; jVar < 3; jVar++) {
                        Aux_Matrix[iVar][jVar] = 0.0;
                        for (kVar = 0; kVar < 3; kVar++)
                            Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar] * D_Matrix[kVar][jVar];
                    }
                }

                /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
                matrix using Gauss integration ---*/

                for (iVar = 0; iVar < nNodes*nVar; iVar++) {
                    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
                        for (kVar = 0; kVar < 3; kVar++) {
                            StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar] * B_Matrix[kVar][jVar] * Det;
                        }
                    }
                }

            }

        }

        void GRID_VolumetricMovement::SetFEA_StiffMatrix3D(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, double **StiffMatrix_Elem, unsigned long PointCorners[8], double CoordCorners[8][3], unsigned short nNodes, double scale) {

            double B_Matrix[6][24], D_Matrix[6][6] = { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } }, Aux_Matrix[24][6];
            double Xi = 0.0, Eta = 0.0, Zeta = 0.0, Det = 0.0, Mu = 0.0, E = 0.0, Lambda = 0.0, Nu = 0.0, Avg_Wall_Dist;
            unsigned short iNode, jNode, iVar, jVar, kVar, iGauss, nGauss = 0;
            double DShapeFunction[8][4] = { { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0 } };
            double Location[8][3], Weight[8];
            unsigned short nVar = geometry->GetnDim();

            for (iVar = 0; iVar < nNodes*nVar; iVar++) {
                for (jVar = 0; jVar < nNodes*nVar; jVar++) {
                    StiffMatrix_Elem[iVar][jVar] = 0.0;
                }
            }

            /*--- Each element uses their own stiffness which is inversely
            proportional to the area/volume of the cell. Using Mu = E & Lambda = -E
            is a modification to help allow rigid rotation of elements (see
            "Robust Mesh Deformation using the Linear Elasticity Equations" by
            R. P. Dwight. ---*/

            /*--- Integration formulae from "Shape functions and points of
            integration of the R¨¦sum¨¦" by Josselin Delmas (2013) ---*/

            /*--- Tetrahedrons. Nodes of numerical integration at 1 point (order 1). ---*/

            if (nNodes == 4) {
                nGauss = 1;
                Location[0][0] = 0.25;  Location[0][1] = 0.25;  Location[0][2] = 0.25;  Weight[0] = 0.166666666666666;
            }

            /*--- Pyramids. Nodes numerical integration at 5 points. ---*/

            if (nNodes == 5) {
                nGauss = 5;
                Location[0][0] = 0.5;   Location[0][1] = 0.0;   Location[0][2] = 0.1531754163448146;  Weight[0] = 0.133333333333333;
                Location[1][0] = 0.0;   Location[1][1] = 0.5;   Location[1][2] = 0.1531754163448146;  Weight[1] = 0.133333333333333;
                Location[2][0] = -0.5;  Location[2][1] = 0.0;   Location[2][2] = 0.1531754163448146;  Weight[2] = 0.133333333333333;
                Location[3][0] = 0.0;   Location[3][1] = -0.5;  Location[3][2] = 0.1531754163448146;  Weight[3] = 0.133333333333333;
                Location[4][0] = 0.0;   Location[4][1] = 0.0;   Location[4][2] = 0.6372983346207416;  Weight[4] = 0.133333333333333;
            }

            /*--- Prism. Nodes of numerical integration at 6 points (order 3 in Xi, order 2 in Eta and Mu ). ---*/

            if (nNodes == 6) {
                nGauss = 6;
                Location[0][0] = 0.5;                 Location[0][1] = 0.5;                 Location[0][2] = -0.577350269189626;  Weight[0] = 0.166666666666666;
                Location[1][0] = -0.577350269189626;  Location[1][1] = 0.0;                 Location[1][2] = 0.5;                 Weight[1] = 0.166666666666666;
                Location[2][0] = 0.5;                 Location[2][1] = -0.577350269189626;  Location[2][2] = 0.0;                 Weight[2] = 0.166666666666666;
                Location[3][0] = 0.5;                 Location[3][1] = 0.5;                 Location[3][2] = 0.577350269189626;   Weight[3] = 0.166666666666666;
                Location[4][0] = 0.577350269189626;   Location[4][1] = 0.0;                 Location[4][2] = 0.5;                 Weight[4] = 0.166666666666666;
                Location[5][0] = 0.5;                 Location[5][1] = 0.577350269189626;   Location[5][2] = 0.0;                 Weight[5] = 0.166666666666666;
            }

            /*--- Hexahedrons. Nodes of numerical integration at 6 points (order 3). ---*/

            if (nNodes == 8) {
                nGauss = 8;
                Location[0][0] = -0.577350269189626;  Location[0][1] = -0.577350269189626;  Location[0][2] = -0.577350269189626;  Weight[0] = 1.0;
                Location[1][0] = -0.577350269189626;  Location[1][1] = -0.577350269189626;  Location[1][2] = 0.577350269189626;   Weight[1] = 1.0;
                Location[2][0] = -0.577350269189626;  Location[2][1] = 0.577350269189626;   Location[2][2] = -0.577350269189626;  Weight[2] = 1.0;
                Location[3][0] = -0.577350269189626;  Location[3][1] = 0.577350269189626;   Location[3][2] = 0.577350269189626;   Weight[3] = 1.0;
                Location[4][0] = 0.577350269189626;   Location[4][1] = -0.577350269189626;  Location[4][2] = -0.577350269189626;  Weight[4] = 1.0;
                Location[5][0] = 0.577350269189626;   Location[5][1] = -0.577350269189626;  Location[5][2] = 0.577350269189626;   Weight[5] = 1.0;
                Location[6][0] = 0.577350269189626;   Location[6][1] = 0.577350269189626;   Location[6][2] = -0.577350269189626;  Weight[6] = 1.0;
                Location[7][0] = 0.577350269189626;   Location[7][1] = 0.577350269189626;   Location[7][2] = 0.577350269189626;   Weight[7] = 1.0;
            }

            for (iGauss = 0; iGauss < nGauss; iGauss++) {

                Xi = Location[iGauss][0]; Eta = Location[iGauss][1];  Zeta = Location[iGauss][2];

                if (nNodes == 4) Det = ShapeFunc_Tetra(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
                if (nNodes == 5) Det = ShapeFunc_Pyram(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
                if (nNodes == 6) Det = ShapeFunc_Prism(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
                if (nNodes == 8) Det = ShapeFunc_Hexa(Xi, Eta, Zeta, CoordCorners, DShapeFunction);

                /*--- Compute the B Matrix ---*/

                for (iVar = 0; iVar < 6; iVar++)
                    for (jVar = 0; jVar < nNodes*nVar; jVar++)
                        B_Matrix[iVar][jVar] = 0.0;

                for (iNode = 0; iNode < nNodes; iNode++) {
                    B_Matrix[0][0 + iNode*nVar] = DShapeFunction[iNode][0];
                    B_Matrix[1][1 + iNode*nVar] = DShapeFunction[iNode][1];
                    B_Matrix[2][2 + iNode*nVar] = DShapeFunction[iNode][2];

                    B_Matrix[3][0 + iNode*nVar] = DShapeFunction[iNode][1];
                    B_Matrix[3][1 + iNode*nVar] = DShapeFunction[iNode][0];

                    B_Matrix[4][1 + iNode*nVar] = DShapeFunction[iNode][2];
                    B_Matrix[4][2 + iNode*nVar] = DShapeFunction[iNode][1];

                    B_Matrix[5][0 + iNode*nVar] = DShapeFunction[iNode][2];
                    B_Matrix[5][2 + iNode*nVar] = DShapeFunction[iNode][0];
                }

                /*--- Impose a type of stiffness for each element ---*/

                switch (config->GetDeform_Stiffness_Type()) {

                case INVERSE_VOLUME:
                    E = scale / (Weight[iGauss] * Det);
                    Mu = E;
                    Lambda = -E;
                    break;

                case WALL_DISTANCE:
                    Avg_Wall_Dist = 0.0;
                    for (jNode = 0; jNode < nNodes; jNode++) {
                        Avg_Wall_Dist += geometry->node[PointCorners[jNode]]->GetWall_Distance() / ((double)nNodes);
                    }
                    E = scale / (Weight[iGauss] * Avg_Wall_Dist);
                    Mu = E;
                    Lambda = -E;
                    break;

                case CONSTANT_STIFFNESS:
                    E = config->GetDeform_ElasticityMod();
                    Nu = config->GetDeform_PoissonRatio();
                    //E = 2E11; Nu = 0.30;
                    Mu = E / (2.0*(1.0 + Nu));
                    Lambda = Nu*E / ((1.0 + Nu)*(1.0 - 2.0*Nu));
                    break;
                }

                /*--- Compute the D Matrix (for plane strain and 3-D)---*/

                D_Matrix[0][0] = Lambda + 2.0*Mu;	D_Matrix[0][1] = Lambda;					D_Matrix[0][2] = Lambda;
                D_Matrix[1][0] = Lambda;					D_Matrix[1][1] = Lambda + 2.0*Mu;	D_Matrix[1][2] = Lambda;
                D_Matrix[2][0] = Lambda;					D_Matrix[2][1] = Lambda;					D_Matrix[2][2] = Lambda + 2.0*Mu;
                D_Matrix[3][3] = Mu;
                D_Matrix[4][4] = Mu;
                D_Matrix[5][5] = Mu;


                /*--- Compute the BT.D Matrix ---*/

                for (iVar = 0; iVar < nNodes*nVar; iVar++) {
                    for (jVar = 0; jVar < 6; jVar++) {
                        Aux_Matrix[iVar][jVar] = 0.0;
                        for (kVar = 0; kVar < 6; kVar++)
                            Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar] * D_Matrix[kVar][jVar];
                    }
                }

                /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
                matrix using Gauss integration ---*/

                for (iVar = 0; iVar < nNodes*nVar; iVar++) {
                    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
                        for (kVar = 0; kVar < 6; kVar++) {
                            StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar] * B_Matrix[kVar][jVar] * Det;
                        }
                    }
                }

            }

        }

        void GRID_VolumetricMovement::AddFEA_StiffMatrix(GEOM::GEOM_Geometry *geometry, double **StiffMatrix_Elem, unsigned long PointCorners[8], unsigned short nNodes) {

            unsigned short iVar, jVar, iDim, jDim;

            unsigned short nVar = geometry->GetnDim();

            double **StiffMatrix_Node;
            StiffMatrix_Node = new double*[nVar];
            for (iVar = 0; iVar < nVar; iVar++)
                StiffMatrix_Node[iVar] = new double[nVar];

            for (iVar = 0; iVar < nVar; iVar++)
                for (jVar = 0; jVar < nVar; jVar++)
                    StiffMatrix_Node[iVar][jVar] = 0.0;

            /*--- Transform the stiffness matrix for the hexahedral element into the
            contributions for the individual nodes relative to each other. ---*/

            for (iVar = 0; iVar < nNodes; iVar++) {
                for (jVar = 0; jVar < nNodes; jVar++) {

                    for (iDim = 0; iDim < nVar; iDim++) {
                        for (jDim = 0; jDim < nVar; jDim++) {
                            StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar*nVar) + iDim][(jVar*nVar) + jDim];
                        }
                    }

                    StiffMatrix.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);

                }
            }

            /*--- Deallocate memory and exit ---*/

            for (iVar = 0; iVar < nVar; iVar++)
                delete StiffMatrix_Node[iVar];
            delete[] StiffMatrix_Node;

        }

        void GRID_VolumetricMovement::SetBoundaryDisplacements(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            unsigned short iDim, nDim = geometry->GetnDim(), iMarker, axis = 0;
            unsigned long iPoint, total_index, iVertex;
            double *VarCoord, MeanCoord[3] = { 0.0, 0.0, 0.0 }, VarIncrement = 1.0;

            /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
            deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
            meshes after imposing design variable surface deformations (DV_MARKER). ---*/

            unsigned short Kind_SU2 = config->GetKind_SU2();

            /*--- If requested (no by default) impose the surface deflections in
            increments and solve the grid deformation equations iteratively with
            successive small deformations. ---*/

            VarIncrement = 1.0 / ((double)config->GetGridDef_Nonlinear_Iter());

            /*--- As initialization, set to zero displacements of all the surfaces except the symmetry
            plane and the receive boundaries. ---*/
            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if ((config->GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE)
                    && (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE)) {
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++) {
                            total_index = iPoint*nDim + iDim;
                            LinSysRes[total_index] = 0.0;
                            LinSysSol[total_index] = 0.0;
                            StiffMatrix.DeleteValsRowi(total_index);
                        }
                    }
                }
            }

            /*--- Set to zero displacements of the normal component for the symmetry plane condition ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if ((config->GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) && (nDim == 3)) {

                    for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = 0.0;
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        VarCoord = geometry->node[iPoint]->GetCoord();
                        for (iDim = 0; iDim < nDim; iDim++)
                            MeanCoord[iDim] += VarCoord[iDim] * VarCoord[iDim];
                    }
                    for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = sqrt(MeanCoord[iDim]);

                    if ((MeanCoord[0] <= MeanCoord[1]) && (MeanCoord[0] <= MeanCoord[2])) axis = 0;
                    if ((MeanCoord[1] <= MeanCoord[0]) && (MeanCoord[1] <= MeanCoord[2])) axis = 1;
                    if ((MeanCoord[2] <= MeanCoord[0]) && (MeanCoord[2] <= MeanCoord[1])) axis = 2;

                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        total_index = iPoint*nDim + axis;
                        LinSysRes[total_index] = 0.0;
                        LinSysSol[total_index] = 0.0;
                        StiffMatrix.DeleteValsRowi(total_index);
                    }
                }
            }

            /*--- Set the known displacements, note that some points of the moving surfaces
            could be on on the symmetry plane, we should specify DeleteValsRowi again (just in case) ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
                    ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF))) {
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
                        for (iDim = 0; iDim < nDim; iDim++) {
                            total_index = iPoint*nDim + iDim;
                            LinSysRes[total_index] = VarCoord[iDim] * VarIncrement;
                            LinSysSol[total_index] = VarCoord[iDim] * VarIncrement;
                            StiffMatrix.DeleteValsRowi(total_index);
                        }
                    }
                }
            }

            /*--- Don't move the nearfield plane ---*/

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
                if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++) {
                            total_index = iPoint*nDim + iDim;
                            LinSysRes[total_index] = 0.0;
                            LinSysSol[total_index] = 0.0;
                            StiffMatrix.DeleteValsRowi(total_index);
                        }
                    }
                }
            }

        }

        void GRID_VolumetricMovement::SetDomainDisplacements(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) {

            unsigned short iDim, nDim = geometry->GetnDim();
            unsigned long iPoint, total_index;
            double *Coord, *MinCoordValues, *MaxCoordValues, *Hold_GridFixed_Coord;

            MinCoordValues = new double[nDim];
            MaxCoordValues = new double[nDim];

            for (iDim = 0; iDim < nDim; iDim++) {
                MinCoordValues[iDim] = 0.0;
                MaxCoordValues[iDim] = 0.0;
            }

            Hold_GridFixed_Coord = config->GetHold_GridFixed_Coord();

            MinCoordValues[0] = Hold_GridFixed_Coord[0];
            MinCoordValues[1] = Hold_GridFixed_Coord[1];
            MinCoordValues[2] = Hold_GridFixed_Coord[2];
            MaxCoordValues[0] = Hold_GridFixed_Coord[3];
            MaxCoordValues[1] = Hold_GridFixed_Coord[4];
            MaxCoordValues[2] = Hold_GridFixed_Coord[5];

            /*--- Set to zero displacements of all the points that are not going to be moved
            except the surfaces ---*/

            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
                Coord = geometry->node[iPoint]->GetCoord();
                for (iDim = 0; iDim < nDim; iDim++) {
                    if ((Coord[iDim] < MinCoordValues[iDim]) || (Coord[iDim] > MaxCoordValues[iDim])) {
                        total_index = iPoint*nDim + iDim;
                        LinSysRes[total_index] = 0.0;
                        LinSysSol[total_index] = 0.0;
                        StiffMatrix.DeleteValsRowi(total_index);
                    }
                }
            }

            delete[] MinCoordValues;
            delete[] MaxCoordValues;

        }

        void GRID_VolumetricMovement::Rigid_Rotation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config,
            unsigned short iZone, unsigned long iter) {

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables ---*/
            unsigned short iDim, nDim;
            unsigned long iPoint;
            double r[3] = { 0.0, 0.0, 0.0 }, rotCoord[3] = { 0.0, 0.0, 0.0 }, *Coord, Center[3] = { 0.0, 0.0, 0.0 }, Omega[3] = { 0.0, 0.0, 0.0 }, Lref, dt, Center_Moment[3] = { 0.0, 0.0, 0.0 };
            double *GridVel, newGridVel[3] = { 0.0, 0.0, 0.0 };
            double rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
            double dtheta, dphi, dpsi, cosTheta, sinTheta;
            double cosPhi, sinPhi, cosPsi, sinPsi;
            bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
            bool adjoint = config->GetAdjoint();

            /*--- Problem dimension and physical time step ---*/
            nDim = geometry->GetnDim();
            dt = config->GetDelta_UnstTimeND();
            Lref = config->GetLength_Ref();

            /*--- For time-spectral, motion is the same in each zone (at each instance).
            *    This is used for calls to the config container ---*/
            if (time_spectral)
                iZone = ZONE_0;

            /*--- For the unsteady adjoint, use reverse time ---*/
            if (adjoint) {
                /*--- Set the first adjoint mesh position to the final direct one ---*/
                if (iter == 0) dt = ((double)config->GetnExtIter() - 1)*dt;
                /*--- Reverse the rotation direction for the adjoint ---*/
                else dt = -1.0*dt;
            }
            else {
                /*--- No rotation at all for the first direct solution ---*/
                if (iter == 0) dt = 0;
            }

            /*--- Center of rotation & angular velocity vector from config ---*/

            Center[0] = config->GetMotion_Origin_X(iZone);
            Center[1] = config->GetMotion_Origin_Y(iZone);
            Center[2] = config->GetMotion_Origin_Z(iZone);
            Omega[0] = (config->GetRotation_Rate_X(iZone) / config->GetOmega_Ref());
            Omega[1] = (config->GetRotation_Rate_Y(iZone) / config->GetOmega_Ref());
            Omega[2] = (config->GetRotation_Rate_Z(iZone) / config->GetOmega_Ref());

            /*-- Set dt for time-spectral cases ---*/
            if (time_spectral) {
                /*--- period of oscillation & compute time interval using nTimeInstances ---*/
                double period = config->GetTimeSpectral_Period();
                dt = period * (double)iter / (double)(config->GetnTimeInstances());
            }

            /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

            dtheta = Omega[0] * dt;
            dphi = Omega[1] * dt;
            dpsi = Omega[2] * dt;

            if (rank == MASTER_NODE && iter == 0) {
                cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
                cout << ", " << Omega[2] << ") rad/s." << endl;
            }

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

            /*--- Loop over and rotate each node in the volume mesh ---*/
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

                /*--- Coordinates of the current point ---*/
                Coord = geometry->node[iPoint]->GetCoord();
                GridVel = geometry->node[iPoint]->GetGridVel();

                /*--- Calculate non-dim. position from rotation center ---*/
                for (iDim = 0; iDim < nDim; iDim++)
                    r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;
                if (nDim == 2) r[nDim] = 0.0;

                /*--- Compute transformed point coordinates ---*/
                rotCoord[0] = rotMatrix[0][0] * r[0]
                    + rotMatrix[0][1] * r[1]
                    + rotMatrix[0][2] * r[2];

                rotCoord[1] = rotMatrix[1][0] * r[0]
                    + rotMatrix[1][1] * r[1]
                    + rotMatrix[1][2] * r[2];

                rotCoord[2] = rotMatrix[2][0] * r[0]
                    + rotMatrix[2][1] * r[1]
                    + rotMatrix[2][2] * r[2];

                /*--- Cross Product of angular velocity and distance from center.
                Note that we have assumed the grid velocities have been set to
                an initial value in the plunging routine. ---*/

                newGridVel[0] = GridVel[0] + Omega[1] * rotCoord[2] - Omega[2] * rotCoord[1];
                newGridVel[1] = GridVel[1] + Omega[2] * rotCoord[0] - Omega[0] * rotCoord[2];
                newGridVel[2] = GridVel[2] + Omega[0] * rotCoord[1] - Omega[1] * rotCoord[0];

                /*--- Store new node location & grid velocity. Add center.
                Do not store the grid velocity if this is an adjoint calculation.---*/

                for (iDim = 0; iDim < nDim; iDim++) {
                    geometry->node[iPoint]->SetCoord(iDim, rotCoord[iDim] + Center[iDim]);
                    if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);

                }
            }

            /*--- Set the moment computation center to the new location after
            incrementing the position with the rotation. ---*/

            for (unsigned short jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {

                Center_Moment[0] = config->GetRefOriginMoment_X(jMarker);
                Center_Moment[1] = config->GetRefOriginMoment_Y(jMarker);
                Center_Moment[2] = config->GetRefOriginMoment_Z(jMarker);

                /*--- Calculate non-dim. position from rotation center ---*/

                for (iDim = 0; iDim < nDim; iDim++)
                    r[iDim] = (Center_Moment[iDim] - Center[iDim]) / Lref;
                if (nDim == 2) r[nDim] = 0.0;

                /*--- Compute transformed point coordinates ---*/

                rotCoord[0] = rotMatrix[0][0] * r[0]
                    + rotMatrix[0][1] * r[1]
                    + rotMatrix[0][2] * r[2];

                rotCoord[1] = rotMatrix[1][0] * r[0]
                    + rotMatrix[1][1] * r[1]
                    + rotMatrix[1][2] * r[2];

                rotCoord[2] = rotMatrix[2][0] * r[0]
                    + rotMatrix[2][1] * r[1]
                    + rotMatrix[2][2] * r[2];

                config->SetRefOriginMoment_X(jMarker, Center[0] + rotCoord[0]);
                config->SetRefOriginMoment_Y(jMarker, Center[1] + rotCoord[1]);
                config->SetRefOriginMoment_Z(jMarker, Center[2] + rotCoord[2]);
            }

            /*--- After moving all nodes, update geometry class ---*/

            UpdateDualGrid(geometry, config);

        }

        void GRID_VolumetricMovement::Rigid_Pitching(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter) {

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables ---*/
            double r[3] = { 0.0, 0.0, 0.0 }, rotCoord[3] = { 0.0, 0.0, 0.0 }, *Coord, Center[3] = { 0.0, 0.0, 0.0 },
                Omega[3] = { 0.0, 0.0, 0.0 }, Ampl[3] = { 0.0, 0.0, 0.0 }, Phase[3] = { 0.0, 0.0, 0.0 };
            double Lref, deltaT, alphaDot[3], *GridVel, newGridVel[3] = { 0.0, 0.0, 0.0 };
            double rotMatrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
            double dtheta, dphi, dpsi, cosTheta, sinTheta;
            double cosPhi, sinPhi, cosPsi, sinPsi;
            double time_new, time_old;
            double DEG2RAD = PI_NUMBER / 180.0;
            unsigned short iDim;
            unsigned short nDim = geometry->GetnDim();
            unsigned long iPoint;
            bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
            bool adjoint = config->GetAdjoint();

            /*--- Retrieve values from the config file ---*/
            deltaT = config->GetDelta_UnstTimeND();
            Lref = config->GetLength_Ref();

            /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
            if (time_spectral) {
                iZone = ZONE_0;
            }

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

            if (time_spectral) {
                /*--- period of oscillation & compute time interval using nTimeInstances ---*/
                double period = config->GetTimeSpectral_Period();
                deltaT = period / (double)(config->GetnTimeInstances());
            }

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
                if (time_spectral) {
                    /*--- For time-spectral, begin movement from the zero position ---*/
                    time_old = 0.0;
                }
                else {
                    time_old = time_new;
                    if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
                }
            }

            /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

            dtheta = -Ampl[0] * (sin(Omega[0] * time_new + Phase[0]) - sin(Omega[0] * time_old + Phase[0]));
            dphi = -Ampl[1] * (sin(Omega[1] * time_new + Phase[1]) - sin(Omega[1] * time_old + Phase[1]));
            dpsi = -Ampl[2] * (sin(Omega[2] * time_new + Phase[2]) - sin(Omega[2] * time_old + Phase[2]));

            /*--- Angular velocity at the new time ---*/

            alphaDot[0] = -Omega[0] * Ampl[0] * cos(Omega[0] * time_new);
            alphaDot[1] = -Omega[1] * Ampl[1] * cos(Omega[1] * time_new);
            alphaDot[2] = -Omega[2] * Ampl[2] * cos(Omega[2] * time_new);

            if (rank == MASTER_NODE && iter == 0) {
                cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
                cout << ", " << Omega[2] << ") rad/s." << endl;
                cout << " Pitching amplitude: (" << Ampl[0] / DEG2RAD << ", ";
                cout << Ampl[1] / DEG2RAD << ", " << Ampl[2] / DEG2RAD;
                cout << ") degrees." << endl;
                cout << " Pitching phase lag: (" << Phase[0] / DEG2RAD << ", ";
                cout << Phase[1] / DEG2RAD << ", " << Phase[2] / DEG2RAD;
                cout << ") degrees." << endl;
            }

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

            /*--- Loop over and rotate each node in the volume mesh ---*/
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

                /*--- Coordinates of the current point ---*/
                Coord = geometry->node[iPoint]->GetCoord();
                GridVel = geometry->node[iPoint]->GetGridVel();

                /*--- Calculate non-dim. position from rotation center ---*/
                for (iDim = 0; iDim < nDim; iDim++)
                    r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;
                if (nDim == 2) r[nDim] = 0.0;

                /*--- Compute transformed point coordinates ---*/
                rotCoord[0] = rotMatrix[0][0] * r[0]
                    + rotMatrix[0][1] * r[1]
                    + rotMatrix[0][2] * r[2];

                rotCoord[1] = rotMatrix[1][0] * r[0]
                    + rotMatrix[1][1] * r[1]
                    + rotMatrix[1][2] * r[2];

                rotCoord[2] = rotMatrix[2][0] * r[0]
                    + rotMatrix[2][1] * r[1]
                    + rotMatrix[2][2] * r[2];

                /*--- Cross Product of angular velocity and distance from center.
                Note that we have assumed the grid velocities have been set to
                an initial value in the plunging routine. ---*/

                newGridVel[0] = GridVel[0] + alphaDot[1] * rotCoord[2] - alphaDot[2] * rotCoord[1];
                newGridVel[1] = GridVel[1] + alphaDot[2] * rotCoord[0] - alphaDot[0] * rotCoord[2];
                newGridVel[2] = GridVel[2] + alphaDot[0] * rotCoord[1] - alphaDot[1] * rotCoord[0];

                /*--- Store new node location & grid velocity. Add center location.
                Do not store the grid velocity if this is an adjoint calculation.---*/

                for (iDim = 0; iDim < nDim; iDim++) {
                    geometry->node[iPoint]->SetCoord(iDim, rotCoord[iDim] + Center[iDim]);
                    if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
                }
            }

            /*--- For pitching we don't update the motion origin and moment reference origin. ---*/

            /*--- After moving all nodes, update geometry class ---*/

            UpdateDualGrid(geometry, config);

        }

        void GRID_VolumetricMovement::Rigid_Plunging(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter) {

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables ---*/
            double deltaX[3], newCoord[3], Center[3], *Coord, Omega[3], Ampl[3], Lref;
            double *GridVel, newGridVel[3], xDot[3];
            double deltaT, time_new, time_old;
            unsigned short iDim, nDim = geometry->GetnDim();
            unsigned long iPoint;
            bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
            bool adjoint = config->GetAdjoint();

            /*--- Retrieve values from the config file ---*/
            deltaT = config->GetDelta_UnstTimeND();
            Lref = config->GetLength_Ref();

            /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
            if (time_spectral) {
                iZone = ZONE_0;
            }

            /*--- Plunging frequency and amplitude from config. ---*/
            Center[0] = config->GetMotion_Origin_X(iZone);
            Center[1] = config->GetMotion_Origin_Y(iZone);
            Center[2] = config->GetMotion_Origin_Z(iZone);
            Omega[0] = (config->GetPlunging_Omega_X(iZone) / config->GetOmega_Ref());
            Omega[1] = (config->GetPlunging_Omega_Y(iZone) / config->GetOmega_Ref());
            Omega[2] = (config->GetPlunging_Omega_Z(iZone) / config->GetOmega_Ref());
            Ampl[0] = config->GetPlunging_Ampl_X(iZone) / Lref;
            Ampl[1] = config->GetPlunging_Ampl_Y(iZone) / Lref;
            Ampl[2] = config->GetPlunging_Ampl_Z(iZone) / Lref;

            if (time_spectral) {
                /*--- period of oscillation & time interval using nTimeInstances ---*/
                double period = config->GetTimeSpectral_Period();
                deltaT = period / (double)(config->GetnTimeInstances());
            }

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
                if (time_spectral) {
                    /*--- For time-spectral, begin movement from the zero position ---*/
                    time_old = 0.0;
                }
                else {
                    time_old = time_new;
                    if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
                }
            }

            /*--- Compute delta change in the position in the x, y, & z directions. ---*/
            deltaX[0] = -Ampl[0] * (sin(Omega[0] * time_new) - sin(Omega[0] * time_old));
            deltaX[1] = -Ampl[1] * (sin(Omega[1] * time_new) - sin(Omega[1] * time_old));
            deltaX[2] = -Ampl[2] * (sin(Omega[2] * time_new) - sin(Omega[2] * time_old));

            /*--- Compute grid velocity due to plunge in the x, y, & z directions. ---*/
            xDot[0] = -Ampl[0] * Omega[0] * (cos(Omega[0] * time_new));
            xDot[1] = -Ampl[1] * Omega[1] * (cos(Omega[1] * time_new));
            xDot[2] = -Ampl[2] * Omega[2] * (cos(Omega[2] * time_new));

            if (rank == MASTER_NODE && iter == 0) {
                cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
                cout << ", " << Omega[2] << ") rad/s." << endl;
                cout << " Plunging amplitude: (" << Ampl[0] << ", ";
                cout << Ampl[1] << ", " << Ampl[2] << ") m." << endl;
            }

            /*--- Loop over and move each node in the volume mesh ---*/
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

                /*--- Coordinates of the current point ---*/
                Coord = geometry->node[iPoint]->GetCoord();
                GridVel = geometry->node[iPoint]->GetGridVel();

                /*--- Increment the node position using the delta values. ---*/
                for (iDim = 0; iDim < nDim; iDim++)
                    newCoord[iDim] = Coord[iDim] + deltaX[iDim];

                /*--- Cross Product of angular velocity and distance from center.
                Note that we have assumed the grid velocities have been set to
                an initial value in the plunging routine. ---*/

                newGridVel[0] = GridVel[0] + xDot[0];
                newGridVel[1] = GridVel[1] + xDot[1];
                newGridVel[2] = GridVel[2] + xDot[2];

                /*--- Store new node location & grid velocity. Do not store the grid
                velocity if this is an adjoint calculation. ---*/

                for (iDim = 0; iDim < nDim; iDim++) {
                    geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
                    if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
                }
            }

            /*--- Set the mesh motion center to the new location after
            incrementing the position with the rigid translation. This
            new location will be used for subsequent pitching/rotation.---*/

            config->SetMotion_Origin_X(iZone, Center[0] + deltaX[0]);
            config->SetMotion_Origin_Y(iZone, Center[1] + deltaX[1]);
            config->SetMotion_Origin_Z(iZone, Center[2] + deltaX[2]);

            /*--- As the body origin may have moved, print it to the console ---*/

            //  if (rank == MASTER_NODE) {
            //    cout << " Body origin: (" << Center[0]+deltaX[0];
            //    cout << ", " << Center[1]+deltaX[1] << ", " << Center[2]+deltaX[2];
            //    cout << ")." << endl;
            //  }

            /*--- Set the moment computation center to the new location after
            incrementing the position with the plunging. ---*/

            for (unsigned short jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
                Center[0] = config->GetRefOriginMoment_X(jMarker) + deltaX[0];
                Center[1] = config->GetRefOriginMoment_Y(jMarker) + deltaX[1];
                Center[2] = config->GetRefOriginMoment_Z(jMarker) + deltaX[2];
                config->SetRefOriginMoment_X(jMarker, Center[0]);
                config->SetRefOriginMoment_Y(jMarker, Center[1]);
                config->SetRefOriginMoment_Z(jMarker, Center[2]);
            }

            /*--- After moving all nodes, update geometry class ---*/

            UpdateDualGrid(geometry, config);

        }

        void GRID_VolumetricMovement::Rigid_Translation(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iZone, unsigned long iter) {

            int rank = MASTER_NODE;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Local variables ---*/
            double deltaX[3], newCoord[3], Center[3], *Coord;
            double xDot[3];
            double deltaT, time_new, time_old;
            unsigned short iDim, nDim = geometry->GetnDim();
            unsigned long iPoint;
            bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
            bool adjoint = config->GetAdjoint();

            /*--- Retrieve values from the config file ---*/
            deltaT = config->GetDelta_UnstTimeND();

            /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
            if (time_spectral) {
                iZone = ZONE_0;
            }

            /*--- Get motion center and translation rates from config ---*/
            Center[0] = config->GetMotion_Origin_X(iZone);
            Center[1] = config->GetMotion_Origin_Y(iZone);
            Center[2] = config->GetMotion_Origin_Z(iZone);
            xDot[0] = config->GetTranslation_Rate_X(iZone);
            xDot[1] = config->GetTranslation_Rate_Y(iZone);
            xDot[2] = config->GetTranslation_Rate_Z(iZone);

            if (time_spectral) {
                /*--- period of oscillation & time interval using nTimeInstances ---*/
                double period = config->GetTimeSpectral_Period();
                deltaT = period / (double)(config->GetnTimeInstances());
            }

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
                if (time_spectral) {
                    /*--- For time-spectral, begin movement from the zero position ---*/
                    time_old = 0.0;
                }
                else {
                    time_old = time_new;
                    if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
                }
            }

            /*--- Compute delta change in the position in the x, y, & z directions. ---*/
            deltaX[0] = xDot[0] * (time_new - time_old);
            deltaX[1] = xDot[1] * (time_new - time_old);
            deltaX[2] = xDot[2] * (time_new - time_old);

            if (rank == MASTER_NODE) {
                cout << " New physical time: " << time_new << " seconds." << endl;
                if (iter == 0) {
                    cout << " Translational velocity: (" << xDot[0] << ", " << xDot[1];
                    cout << ", " << xDot[2] << ") m/s." << endl;
                }
            }

            /*--- Loop over and move each node in the volume mesh ---*/
            for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

                /*--- Coordinates of the current point ---*/
                Coord = geometry->node[iPoint]->GetCoord();

                /*--- Increment the node position using the delta values. ---*/
                for (iDim = 0; iDim < nDim; iDim++)
                    newCoord[iDim] = Coord[iDim] + deltaX[iDim];

                /*--- Store new node location & grid velocity. Do not store the grid
                velocity if this is an adjoint calculation. ---*/

                for (iDim = 0; iDim < nDim; iDim++) {
                    geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
                    if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, xDot[iDim]);
                }
            }

            /*--- Set the mesh motion center to the new location after
            incrementing the position with the rigid translation. This
            new location will be used for subsequent pitching/rotation.---*/

            config->SetMotion_Origin_X(iZone, Center[0] + deltaX[0]);
            config->SetMotion_Origin_Y(iZone, Center[1] + deltaX[1]);
            config->SetMotion_Origin_Z(iZone, Center[2] + deltaX[2]);

            /*--- Set the moment computation center to the new location after
            incrementing the position with the translation. ---*/

            for (unsigned short jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
                Center[0] = config->GetRefOriginMoment_X(jMarker) + deltaX[0];
                Center[1] = config->GetRefOriginMoment_Y(jMarker) + deltaX[1];
                Center[2] = config->GetRefOriginMoment_Z(jMarker) + deltaX[2];
                config->SetRefOriginMoment_X(jMarker, Center[0]);
                config->SetRefOriginMoment_Y(jMarker, Center[1]);
                config->SetRefOriginMoment_Z(jMarker, Center[2]);
            }

            /*--- After moving all nodes, update geometry class ---*/

            UpdateDualGrid(geometry, config);

        }
    }
}