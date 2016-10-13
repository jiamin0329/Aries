
#include "GRID_FreeFormDefBox.hpp"

namespace ARIES
{
    namespace GRID
    {

        GRID_FreeFormDefBox::GRID_FreeFormDefBox(void) : GRID_Gridmovement() { }

        GRID_FreeFormDefBox::GRID_FreeFormDefBox(unsigned short val_lDegree, unsigned short val_mDegree, unsigned short val_nDegree) : GRID_Gridmovement() {

            unsigned short iCornerPoints, iOrder, jOrder, kOrder, iDim;

            /*--- FFD is always 3D (even in 2D problems) ---*/

            nDim = 3;
            nCornerPoints = 8;

            /*--- Allocate Corners points ---*/

            Coord_Corner_Points = new double*[nCornerPoints];
            for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
                Coord_Corner_Points[iCornerPoints] = new double[nDim];

            ParamCoord = new double[nDim]; ParamCoord_ = new double[nDim];
            cart_coord = new double[nDim]; cart_coord_ = new double[nDim];
            Gradient = new double[nDim];

            lDegree = val_lDegree; lOrder = lDegree + 1;
            mDegree = val_mDegree; mOrder = mDegree + 1;
            nDegree = val_nDegree; nOrder = nDegree + 1;
            nControlPoints = lOrder*mOrder*nOrder;

            lDegree_Copy = val_lDegree; lOrder_Copy = lDegree + 1;
            mDegree_Copy = val_mDegree; mOrder_Copy = mDegree + 1;
            nDegree_Copy = val_nDegree; nOrder_Copy = nDegree + 1;
            nControlPoints_Copy = lOrder_Copy*mOrder_Copy*nOrder_Copy;

            Coord_Control_Points = new double***[lOrder];
            ParCoord_Control_Points = new double***[lOrder];
            Coord_Control_Points_Copy = new double***[lOrder];
            for (iOrder = 0; iOrder < lOrder; iOrder++) {
                Coord_Control_Points[iOrder] = new double**[mOrder];
                ParCoord_Control_Points[iOrder] = new double**[mOrder];
                Coord_Control_Points_Copy[iOrder] = new double**[mOrder];
                for (jOrder = 0; jOrder < mOrder; jOrder++) {
                    Coord_Control_Points[iOrder][jOrder] = new double*[nOrder];
                    ParCoord_Control_Points[iOrder][jOrder] = new double*[nOrder];
                    Coord_Control_Points_Copy[iOrder][jOrder] = new double*[nOrder];
                    for (kOrder = 0; kOrder < nOrder; kOrder++) {
                        Coord_Control_Points[iOrder][jOrder][kOrder] = new double[nDim];
                        ParCoord_Control_Points[iOrder][jOrder][kOrder] = new double[nDim];
                        Coord_Control_Points_Copy[iOrder][jOrder][kOrder] = new double[nDim];
                        for (iDim = 0; iDim < nDim; iDim++) {
                            Coord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
                            ParCoord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
                            Coord_Control_Points_Copy[iOrder][jOrder][kOrder][iDim] = 0.0;
                        }
                    }
                }
            }

        }

        GRID_FreeFormDefBox::~GRID_FreeFormDefBox(void) {
            unsigned short iOrder, jOrder, kOrder, iCornerPoints;

            for (iOrder = 0; iOrder < lOrder; iOrder++)
                for (jOrder = 0; jOrder < mOrder; jOrder++)
                    for (kOrder = 0; kOrder < nOrder; kOrder++) {
                        delete[] Coord_Control_Points[iOrder][jOrder][kOrder];
                        delete[] ParCoord_Control_Points[iOrder][jOrder][kOrder];
                        delete[] Coord_Control_Points_Copy[iOrder][jOrder][kOrder];
                    }
            delete[] Coord_Control_Points;
            delete[] ParCoord_Control_Points;
            delete[] Coord_Control_Points_Copy;

            delete[] ParamCoord;
            delete[] cart_coord;
            delete[] Gradient;

            for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
                delete[] Coord_Corner_Points[iCornerPoints];
            delete[] Coord_Corner_Points;
        }

        void  GRID_FreeFormDefBox::SetUnitCornerPoints(void) {

            unsigned short iDim;
            double *coord = new double[nDim];

            for (iDim = 0; iDim < nDim; iDim++) coord[iDim] = 0.0;

            coord[0] = 0.0; coord[1] = 0.0; coord[2] = 0.0; this->SetCoordCornerPoints(coord, 0);
            coord[0] = 1.0; coord[1] = 0.0; coord[2] = 0.0; this->SetCoordCornerPoints(coord, 1);
            coord[0] = 1.0; coord[1] = 1.0; coord[2] = 0.0; this->SetCoordCornerPoints(coord, 2);
            coord[0] = 0.0; coord[1] = 1.0; coord[2] = 0.0; this->SetCoordCornerPoints(coord, 3);
            coord[0] = 0.0; coord[1] = 0.0; coord[2] = 1.0; this->SetCoordCornerPoints(coord, 4);
            coord[0] = 1.0; coord[1] = 0.0; coord[2] = 1.0; this->SetCoordCornerPoints(coord, 5);
            coord[0] = 1.0; coord[1] = 1.0; coord[2] = 1.0; this->SetCoordCornerPoints(coord, 6);
            coord[0] = 0.0; coord[1] = 1.0; coord[2] = 1.0; this->SetCoordCornerPoints(coord, 7);

            delete[] coord;

        }

        void GRID_FreeFormDefBox::SetControlPoints_Parallelepiped(void) {
            unsigned short iDim, iDegree, jDegree, kDegree;

            /*--- Set base control points according to the notation of Vtk for hexahedrons ---*/
            for (iDim = 0; iDim < nDim; iDim++) {
                Coord_Control_Points[0][0][0][iDim] = Coord_Corner_Points[0][iDim];
                Coord_Control_Points[lOrder - 1][0][0][iDim] = Coord_Corner_Points[1][iDim];
                Coord_Control_Points[lOrder - 1][mOrder - 1][0][iDim] = Coord_Corner_Points[2][iDim];
                Coord_Control_Points[0][mOrder - 1][0][iDim] = Coord_Corner_Points[3][iDim];
                Coord_Control_Points[0][0][nOrder - 1][iDim] = Coord_Corner_Points[4][iDim];
                Coord_Control_Points[lOrder - 1][0][nOrder - 1][iDim] = Coord_Corner_Points[5][iDim];
                Coord_Control_Points[lOrder - 1][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[6][iDim];
                Coord_Control_Points[0][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[7][iDim];
            }

            /*--- Fill the rest of the cubic matrix of control points with uniform spacing (parallelepiped) ---*/
            for (iDegree = 0; iDegree <= lDegree; iDegree++)
                for (jDegree = 0; jDegree <= mDegree; jDegree++)
                    for (kDegree = 0; kDegree <= nDegree; kDegree++) {
                        Coord_Control_Points[iDegree][jDegree][kDegree][0] = Coord_Corner_Points[0][0]
                            + double(iDegree) / double(lDegree)*(Coord_Corner_Points[1][0] - Coord_Corner_Points[0][0]);
                        Coord_Control_Points[iDegree][jDegree][kDegree][1] = Coord_Corner_Points[0][1]
                            + double(jDegree) / double(mDegree)*(Coord_Corner_Points[3][1] - Coord_Corner_Points[0][1]);
                        Coord_Control_Points[iDegree][jDegree][kDegree][2] = Coord_Corner_Points[0][2]
                            + double(kDegree) / double(nDegree)*(Coord_Corner_Points[4][2] - Coord_Corner_Points[0][2]);
                    }
        }

        void GRID_FreeFormDefBox::SetSupportCP(GRID_FreeFormDefBox *FFDBox) {
            unsigned short iDim, iOrder, jOrder, kOrder;
            unsigned short lOrder = FFDBox->GetlOrder();
            unsigned short mOrder = FFDBox->GetmOrder();
            unsigned short nOrder = FFDBox->GetnOrder();

            Coord_SupportCP = new double***[lOrder];
            for (iOrder = 0; iOrder < lOrder; iOrder++) {
                Coord_SupportCP[iOrder] = new double**[mOrder];
                for (jOrder = 0; jOrder < mOrder; jOrder++) {
                    Coord_SupportCP[iOrder][jOrder] = new double*[nOrder];
                    for (kOrder = 0; kOrder < nOrder; kOrder++)
                        Coord_SupportCP[iOrder][jOrder][kOrder] = new double[nDim];
                }
            }

            /*--- Set base support control points according to the notation of Vtk for hexahedrons ---*/
            for (iDim = 0; iDim < nDim; iDim++) {
                Coord_SupportCP[0][0][0][iDim] = Coord_Corner_Points[0][iDim];
                Coord_SupportCP[lOrder - 1][0][0][iDim] = Coord_Corner_Points[1][iDim];
                Coord_SupportCP[lOrder - 1][mOrder - 1][0][iDim] = Coord_Corner_Points[2][iDim];
                Coord_SupportCP[0][mOrder - 1][0][iDim] = Coord_Corner_Points[3][iDim];
                Coord_SupportCP[0][0][nOrder - 1][iDim] = Coord_Corner_Points[4][iDim];
                Coord_SupportCP[lOrder - 1][0][nOrder - 1][iDim] = Coord_Corner_Points[5][iDim];
                Coord_SupportCP[lOrder - 1][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[6][iDim];
                Coord_SupportCP[0][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[7][iDim];
            }

            /*--- Fill the rest of the cubic matrix of support control points with uniform spacing  ---*/
            for (iOrder = 0; iOrder < lOrder; iOrder++)
                for (jOrder = 0; jOrder < mOrder; jOrder++)
                    for (kOrder = 0; kOrder < nOrder; kOrder++) {
                        Coord_SupportCP[iOrder][jOrder][kOrder][0] = Coord_Corner_Points[0][0]
                            + double(iOrder) / double(lOrder - 1)*(Coord_Corner_Points[1][0] - Coord_Corner_Points[0][0]);
                        Coord_SupportCP[iOrder][jOrder][kOrder][1] = Coord_Corner_Points[0][1]
                            + double(jOrder) / double(mOrder - 1)*(Coord_Corner_Points[3][1] - Coord_Corner_Points[0][1]);
                        Coord_SupportCP[iOrder][jOrder][kOrder][2] = Coord_Corner_Points[0][2]
                            + double(kOrder) / double(nOrder - 1)*(Coord_Corner_Points[4][2] - Coord_Corner_Points[0][2]);
                    }
        }

        void GRID_FreeFormDefBox::SetSupportCPChange(GRID_FreeFormDefBox *FFDBox) {
            unsigned short iDim, iOrder, jOrder, kOrder;
            double *CartCoordNew, *ParamCoord;
            unsigned short lOrder = FFDBox->GetlOrder();
            unsigned short mOrder = FFDBox->GetmOrder();
            unsigned short nOrder = FFDBox->GetnOrder();

            double ****ParamCoord_SupportCP = new double***[lOrder];
            for (iOrder = 0; iOrder < lOrder; iOrder++) {
                ParamCoord_SupportCP[iOrder] = new double**[mOrder];
                for (jOrder = 0; jOrder < mOrder; jOrder++) {
                    ParamCoord_SupportCP[iOrder][jOrder] = new double*[nOrder];
                    for (kOrder = 0; kOrder < nOrder; kOrder++)
                        ParamCoord_SupportCP[iOrder][jOrder][kOrder] = new double[nDim];
                }
            }

            for (iOrder = 0; iOrder < lOrder; iOrder++)
                for (jOrder = 0; jOrder < mOrder; jOrder++)
                    for (kOrder = 0; kOrder < nOrder; kOrder++)
                        for (iDim = 0; iDim < nDim; iDim++)
                            ParamCoord_SupportCP[iOrder][jOrder][kOrder][iDim] =
                            Coord_SupportCP[iOrder][jOrder][kOrder][iDim];

            for (iDim = 0; iDim < nDim; iDim++) {
                Coord_Control_Points[0][0][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 0);
                Coord_Control_Points[1][0][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 1);
                Coord_Control_Points[1][1][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 2);
                Coord_Control_Points[0][1][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 3);
                Coord_Control_Points[0][0][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 4);
                Coord_Control_Points[1][0][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 5);
                Coord_Control_Points[1][1][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 6);
                Coord_Control_Points[0][1][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 7);
            }

            for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++) {
                for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++) {
                    for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
                        ParamCoord = ParamCoord_SupportCP[iOrder][jOrder][kOrder];
                        CartCoordNew = EvalCartesianCoord(ParamCoord);
                        FFDBox->SetCoordControlPoints(CartCoordNew, iOrder, jOrder, kOrder);
                        FFDBox->SetCoordControlPoints_Copy(CartCoordNew, iOrder, jOrder, kOrder);
                    }
                }
            }

        }

        void GRID_FreeFormDefBox::SetTecplot(GEOM::GEOM_Geometry *geometry, unsigned short iFFDBox, bool original) {

            std::ofstream FFDBox_file;
            char FFDBox_filename[TBOX::MAX_STRING_SIZE];
            bool new_file;
            unsigned short iDim, iDegree, jDegree, kDegree;

            nDim = geometry->GetnDim();

            sprintf(FFDBox_filename, "ffd_boxes.dat");

            if ((original) && (iFFDBox == 0)) new_file = true;
            else new_file = false;

            if (new_file) {
                FFDBox_file.open(FFDBox_filename, std::ios::out);
                FFDBox_file << "TITLE = \"Visualization of the FFD boxes generated by SU2_DEF.\"" << std::endl;
                if (nDim == 2) FFDBox_file << "VARIABLES = \"x\", \"y\"" << std::endl;
                else FFDBox_file << "VARIABLES = \"x\", \"y\", \"z\"" << std::endl;
            }
            else FFDBox_file.open(FFDBox_filename, std::ios::out | std::ios::app);

            FFDBox_file << "ZONE T= \"" << Tag;
            if (original) FFDBox_file << " (Original FFD)\"";
            else FFDBox_file << " (Deformed FFD)\"";
            if (nDim == 2) FFDBox_file << ", I=" << lDegree + 1 << ", J=" << mDegree + 1 << ", DATAPACKING=POINT" << std::endl;
            else FFDBox_file << ", I=" << lDegree + 1 << ", J=" << mDegree + 1 << ", K=" << nDegree + 1 << ", DATAPACKING=POINT" << std::endl;

            FFDBox_file.precision(15);

            if (nDim == 2) {
                for (jDegree = 0; jDegree <= mDegree; jDegree++) {
                    for (iDegree = 0; iDegree <= lDegree; iDegree++) {
                        for (iDim = 0; iDim < nDim; iDim++)
                            FFDBox_file << std::scientific << Coord_Control_Points[iDegree][jDegree][0][iDim] << "\t";
                        FFDBox_file << "\n";
                    }
                }
            }
            else {
                for (kDegree = 0; kDegree <= nDegree; kDegree++) {
                    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
                        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
                            for (iDim = 0; iDim < nDim; iDim++)
                                FFDBox_file << std::scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
                            FFDBox_file << "\n";
                        }
                    }
                }
            }

            FFDBox_file.close();
        }


        double *GRID_FreeFormDefBox::GetParametricCoord_Analytical(double *cart_coord) {
            unsigned short iDim;
            double *e1, *e2, *e3, *e12, *e23, *e13, *p;

            /*--- Auxiliary Basis Vectors of the deformed FFDBox ---*/
            e1 = new double[3]; e2 = new double[3]; e3 = new double[3];
            for (iDim = 0; iDim < nDim; iDim++) {
                e1[iDim] = Coord_Corner_Points[1][iDim] - Coord_Corner_Points[0][iDim];
                e2[iDim] = Coord_Corner_Points[3][iDim] - Coord_Corner_Points[0][iDim];
                e3[iDim] = Coord_Corner_Points[4][iDim] - Coord_Corner_Points[0][iDim];
            }

            /*--- Respective Cross-Products ---*/
            e12 = new double[3]; e23 = new double[3]; e13 = new double[3];
            CrossProduct(e1, e2, e12);
            CrossProduct(e1, e3, e13);
            CrossProduct(e2, e3, e23);

            /*--- p is Tranlated vector from the origin ---*/
            p = new double[3];
            for (iDim = 0; iDim < nDim; iDim++)
                p[iDim] = cart_coord[iDim] - Coord_Corner_Points[0][iDim];

            ParamCoord[0] = DotProduct(e23, p) / DotProduct(e23, e1);
            ParamCoord[1] = DotProduct(e13, p) / DotProduct(e13, e2);
            ParamCoord[2] = DotProduct(e12, p) / DotProduct(e12, e3);

            delete[] e1;
            delete[] e2;
            delete[] e3;
            delete[] e12;
            delete[] e23;
            delete[] e13;
            delete[] p;

            return ParamCoord;
        }

        double *GRID_FreeFormDefBox::EvalCartesianCoord(double *ParamCoord) {
            unsigned short iDim, iDegree, jDegree, kDegree;

            for (iDim = 0; iDim < nDim; iDim++)
                cart_coord[iDim] = 0.0;

            for (iDegree = 0; iDegree <= lDegree; iDegree++)
                for (jDegree = 0; jDegree <= mDegree; jDegree++)
                    for (kDegree = 0; kDegree <= nDegree; kDegree++)
                        for (iDim = 0; iDim < nDim; iDim++) {
                            cart_coord[iDim] += Coord_Control_Points[iDegree][jDegree][kDegree][iDim]
                                * GetBernstein(lDegree, iDegree, ParamCoord[0])
                                * GetBernstein(mDegree, jDegree, ParamCoord[1])
                                * GetBernstein(nDegree, kDegree, ParamCoord[2]);
                        }

            return cart_coord;
        }

        double GRID_FreeFormDefBox::GetBernstein(short val_n, short val_i, double val_t) {

            double value = 0.0;

            if (val_i > val_n) { value = 0; return value; }
            if (val_i == 0) {
                if (val_t == 0) value = 1;
                else if (val_t == 1) value = 0;
                else value = Binomial(val_n, val_i)*(pow(val_t, val_i)) * pow(1.0 - val_t, val_n - val_i);
            }
            else if (val_i == val_n) {
                if (val_t == 0) value = 0;
                else if (val_t == 1) value = 1;
                else value = pow(val_t, val_n);
            }
            else value = Binomial(val_n, val_i)*(pow(val_t, val_i)) * pow(1.0 - val_t, val_n - val_i);

            return value;
        }

        double GRID_FreeFormDefBox::GetBernsteinDerivative(short val_n, short val_i,
            double val_t, short val_order) {
            double value = 0.0;

            /*--- Verify this subroutine, it provides negative val_n,
            which is a wrong value for GetBernstein ---*/

            if (val_order == 0) {
                value = GetBernstein(val_n, val_i, val_t); return value;
            }

            if (val_i == 0) {
                value = val_n*(-GetBernsteinDerivative(val_n - 1, val_i, val_t, val_order - 1)); return value;
            }
            else {
                if (val_n == 0) {
                    value = val_t; return value;
                }
                else {
                    value = val_n*(GetBernsteinDerivative(val_n - 1, val_i - 1, val_t, val_order - 1) - GetBernsteinDerivative(val_n - 1, val_i, val_t, val_order - 1));
                    return value;
                }
            }

            return value;
        }

        double *GRID_FreeFormDefBox::GetFFDGradient(double *val_coord, double *xyz) {

            unsigned short iDim, jDim, lmn[3];

            /*--- Set the Degree of the Berstein polynomials ---*/

            lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;

            for (iDim = 0; iDim < nDim; iDim++) Gradient[iDim] = 0.0;

            for (iDim = 0; iDim < nDim; iDim++)
                for (jDim = 0; jDim < nDim; jDim++)
                    Gradient[jDim] += GetDerivative2(val_coord, iDim, xyz, lmn) *
                    GetDerivative3(val_coord, iDim, jDim, lmn);

            return Gradient;

        }

        void GRID_FreeFormDefBox::GetFFDHessian(double *uvw, double *xyz, double **val_Hessian) {

            unsigned short iDim, jDim, lmn[3];

            /*--- Set the Degree of the Berstein polynomials ---*/

            lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;

            for (iDim = 0; iDim < nDim; iDim++)
                for (jDim = 0; jDim < nDim; jDim++)
                    val_Hessian[iDim][jDim] = 0.0;

            /*--- Note that being all the functions linear combinations of polynomials, they are C^\infty,
            and the Hessian will be symmetric; no need to compute the under-diagonal part, for example ---*/

            for (iDim = 0; iDim < nDim; iDim++) {
                val_Hessian[0][0] += 2.0 * GetDerivative3(uvw, iDim, 0, lmn) * GetDerivative3(uvw, iDim, 0, lmn) +
                    GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 0, 0, lmn);

                val_Hessian[1][1] += 2.0 * GetDerivative3(uvw, iDim, 1, lmn) * GetDerivative3(uvw, iDim, 1, lmn) +
                    GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 1, 1, lmn);

                val_Hessian[2][2] += 2.0 * GetDerivative3(uvw, iDim, 2, lmn) * GetDerivative3(uvw, iDim, 2, lmn) +
                    GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 2, 2, lmn);

                val_Hessian[0][1] += 2.0 * GetDerivative3(uvw, iDim, 0, lmn) * GetDerivative3(uvw, iDim, 1, lmn) +
                    GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 0, 1, lmn);

                val_Hessian[0][2] += 2.0 * GetDerivative3(uvw, iDim, 0, lmn) * GetDerivative3(uvw, iDim, 2, lmn) +
                    GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 0, 2, lmn);

                val_Hessian[1][2] += 2.0 * GetDerivative3(uvw, iDim, 1, lmn) * GetDerivative3(uvw, iDim, 2, lmn) +
                    GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 1, 2, lmn);
            }

            val_Hessian[1][0] = val_Hessian[0][1];
            val_Hessian[2][0] = val_Hessian[0][2];
            val_Hessian[2][1] = val_Hessian[1][2];

        }

        double *GRID_FreeFormDefBox::GetParametricCoord_Iterative(unsigned long iPoint, double *xyz, double *ParamCoordGuess, TBOX::TBOX_Config *config) {

            double *IndepTerm, SOR_Factor = 1.0, MinNormError, NormError, Determinant, AdjHessian[3][3], Temp[3] = { 0.0, 0.0, 0.0 };
            unsigned short iDim, jDim, RandonCounter;
            unsigned long iter;

            double tol = config->GetFFD_Tol();
            unsigned short it_max = config->GetnFFD_Iter();
            unsigned short Random_Trials = 500;

            /*--- Allocate the Hessian ---*/

            Hessian = new double*[nDim];
            IndepTerm = new double[nDim];
            for (iDim = 0; iDim < nDim; iDim++)
            {
                Hessian[iDim] = new double[nDim];
                ParamCoord[iDim] = ParamCoordGuess[iDim];
                IndepTerm[iDim] = 0.0;
            }

            RandonCounter = 0; MinNormError = 1E6;

            /*--- External iteration ---*/

            for (iter = 0; iter < it_max*Random_Trials; iter++) {

                /*--- The independent term of the solution of our system is -Gradient(sol_old) ---*/

                Gradient = GetFFDGradient(ParamCoord, xyz);

                for (iDim = 0; iDim < nDim; iDim++) IndepTerm[iDim] = -Gradient[iDim];

                /*--- Hessian = The Matrix of our system, getHessian(sol_old,xyz,...) ---*/

                GetFFDHessian(ParamCoord, xyz, Hessian);

                /*--- Adjoint to Hessian ---*/

                AdjHessian[0][0] = Hessian[1][1] * Hessian[2][2] - Hessian[1][2] * Hessian[2][1];
                AdjHessian[0][1] = Hessian[0][2] * Hessian[2][1] - Hessian[0][1] * Hessian[2][2];
                AdjHessian[0][2] = Hessian[0][1] * Hessian[1][2] - Hessian[0][2] * Hessian[1][1];
                AdjHessian[1][0] = Hessian[1][2] * Hessian[2][0] - Hessian[1][0] * Hessian[2][2];
                AdjHessian[1][1] = Hessian[0][0] * Hessian[2][2] - Hessian[0][2] * Hessian[2][0];
                AdjHessian[1][2] = Hessian[0][2] * Hessian[1][0] - Hessian[0][0] * Hessian[1][2];
                AdjHessian[2][0] = Hessian[1][0] * Hessian[2][1] - Hessian[1][1] * Hessian[2][0];
                AdjHessian[2][1] = Hessian[0][1] * Hessian[2][0] - Hessian[0][0] * Hessian[2][1];
                AdjHessian[2][2] = Hessian[0][0] * Hessian[1][1] - Hessian[0][1] * Hessian[1][0];

                /*--- Determinant of Hessian ---*/

                Determinant = Hessian[0][0] * AdjHessian[0][0] + Hessian[0][1] * AdjHessian[1][0] + Hessian[0][2] * AdjHessian[2][0];

                /*--- Hessian inverse ---*/

                if (Determinant != 0) {
                    for (iDim = 0; iDim < nDim; iDim++) {
                        Temp[iDim] = 0.0;
                        for (jDim = 0; jDim < nDim; jDim++) {
                            Temp[iDim] += AdjHessian[iDim][jDim] * IndepTerm[jDim] / Determinant;
                        }
                    }
                    for (iDim = 0; iDim < nDim; iDim++) {
                        IndepTerm[iDim] = Temp[iDim];
                    }
                }

                /*--- Update with Successive over-relaxation ---*/

                for (iDim = 0; iDim < nDim; iDim++) {
                    ParamCoord[iDim] = (1.0 - SOR_Factor)*ParamCoord[iDim] + SOR_Factor*(ParamCoord[iDim] + IndepTerm[iDim]);
                }

                /*--- If the gradient is small, we have converged ---*/

                if ((fabs(IndepTerm[0]) < tol) && (fabs(IndepTerm[1]) < tol) && (fabs(IndepTerm[2]) < tol))	break;

                /*--- Compute the norm of the error ---*/

                NormError = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                    NormError += IndepTerm[iDim] * IndepTerm[iDim];
                NormError = sqrt(NormError);

                MinNormError = std::min(NormError, MinNormError);

                /*--- If we have no convergence with Random_Trials iterations probably we are in a local minima. ---*/

                if (((iter % it_max) == 0) && (iter != 0)) {

                    RandonCounter++;
                    if (RandonCounter == Random_Trials) {
                        std::cout << std::endl << "Unknown point: " << iPoint << " (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "). Min Error: " << MinNormError << ". Iter: " << iter << "." << std::endl;
                    }
                    else {
                        SOR_Factor = 0.1;
                        for (iDim = 0; iDim < nDim; iDim++)
                            ParamCoord[iDim] = double(rand()) / double(RAND_MAX);
                    }

                }

            }

            for (iDim = 0; iDim < nDim; iDim++)
                delete[] Hessian[iDim];
            delete[] Hessian;
            delete[] IndepTerm;

            /*--- The code has hit the max number of iterations ---*/

            if (iter == it_max*Random_Trials) {
                std::cout << "Unknown point: (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "). Increase the value of FFD_ITERATIONS." << std::endl;
            }

            /*--- Real Solution is now ParamCoord; Return it ---*/

            return ParamCoord;

        }

        unsigned long GRID_FreeFormDefBox::Binomial(unsigned short n, unsigned short m) {

            unsigned short i, j;
            unsigned long binomial[1000];

            binomial[0] = 1;
            for (i = 1; i <= n; ++i) {
                binomial[i] = 1;
                for (j = i - 1U; j > 0; --j) {
                    binomial[j] += binomial[j - 1U];
                }
            }

            return binomial[m];

        }

        bool GRID_FreeFormDefBox::GetPointFFD(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned long iPoint) {
            double Coord[3] = { 0.0, 0.0, 0.0 };
            unsigned short iVar, jVar, iDim;
            bool Inside = false;

            unsigned short Index[5][7] = {
                { 0, 1, 2, 5, 0, 1, 2 },
                { 0, 2, 7, 5, 0, 2, 7 },
                { 0, 2, 3, 7, 0, 2, 3 },
                { 0, 5, 7, 4, 0, 5, 7 },
                { 2, 7, 5, 6, 2, 7, 5 } };
            unsigned short nDim = geometry->GetnDim();

            for (iDim = 0; iDim < nDim; iDim++)
                Coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);

            /*--- 1st tetrahedron {V0, V1, V2, V5}
            2nd tetrahedron {V0, V2, V7, V5}
            3th tetrahedron {V0, V2, V3, V7}
            4th tetrahedron {V0, V5, V7, V4}
            5th tetrahedron {V2, V7, V5, V6} ---*/

            for (iVar = 0; iVar < 5; iVar++) {
                Inside = true;
                for (jVar = 0; jVar < 4; jVar++) {
                    double Distance_Point = geometry->Point2Plane_Distance(Coord,
                        Coord_Corner_Points[Index[iVar][jVar + 1]],
                        Coord_Corner_Points[Index[iVar][jVar + 2]],
                        Coord_Corner_Points[Index[iVar][jVar + 3]]);

                    double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]],
                        Coord_Corner_Points[Index[iVar][jVar + 1]],
                        Coord_Corner_Points[Index[iVar][jVar + 2]],
                        Coord_Corner_Points[Index[iVar][jVar + 3]]);
                    if (Distance_Point*Distance_Vertex < 0.0) Inside = false;
                }
                if (Inside) break;
            }

            return Inside;

        }

        void GRID_FreeFormDefBox::SetDeformationZone(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iFFDBox) {
            double *Coord;
            unsigned short iMarker, iVar, jVar;
            unsigned long iVertex, iPoint;
            bool Inside = false;

            unsigned short Index[5][7] = {
                { 0, 1, 2, 5, 0, 1, 2 },
                { 0, 2, 7, 5, 0, 2, 7 },
                { 0, 2, 3, 7, 0, 2, 3 },
                { 0, 5, 7, 4, 0, 5, 7 },
                { 2, 7, 5, 6, 2, 7, 5 } };

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
                if (config->GetMarker_All_DV(iMarker) == TBOX::YES)
                    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
                        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
                        geometry->node[iPoint]->SetMove(false);

                        Coord = geometry->node[iPoint]->GetCoord();

                        /*--- 1st tetrahedron {V0, V1, V2, V5}
                        2nd tetrahedron {V0, V2, V7, V5}
                        3th tetrahedron {V0, V2, V3, V7}
                        4th tetrahedron {V0, V5, V7, V4}
                        5th tetrahedron {V2, V7, V5, V6} ---*/

                        for (iVar = 0; iVar < 5; iVar++) {
                            Inside = true;
                            for (jVar = 0; jVar < 4; jVar++) {
                                double Distance_Point = geometry->Point2Plane_Distance(Coord,
                                    Coord_Corner_Points[Index[iVar][jVar + 1]],
                                    Coord_Corner_Points[Index[iVar][jVar + 2]],
                                    Coord_Corner_Points[Index[iVar][jVar + 3]]);
                                double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]],
                                    Coord_Corner_Points[Index[iVar][jVar + 1]],
                                    Coord_Corner_Points[Index[iVar][jVar + 2]],
                                    Coord_Corner_Points[Index[iVar][jVar + 3]]);
                                if (Distance_Point*Distance_Vertex < 0.0) Inside = false;
                            }
                            if (Inside) break;
                        }

                        if (Inside) {
                            geometry->node[iPoint]->SetMove(true);
                        }

                    }
        }

        double GRID_FreeFormDefBox::GetDerivative1(double *uvw, unsigned short val_diff, unsigned short *ijk, unsigned short *lmn) {

            unsigned short iDim;
            double value = 0.0;

            value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 1);
            for (iDim = 0; iDim < nDim; iDim++)
                if (iDim != val_diff)
                    value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);

            return value;

        }

        double GRID_FreeFormDefBox::GetDerivative2(double *uvw, unsigned short dim, double *xyz, unsigned short *lmn) {

            unsigned short iDegree, jDegree, kDegree;
            double value = 0.0;

            for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
                for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
                    for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
                        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim]
                            * GetBernstein(lmn[0], iDegree, uvw[0])
                            * GetBernstein(lmn[1], jDegree, uvw[1])
                            * GetBernstein(lmn[2], kDegree, uvw[2]);
                    }

            return 2.0*(value - xyz[dim]);
        }

        double GRID_FreeFormDefBox::GetDerivative3(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short *lmn) {

            unsigned short iDegree, jDegree, kDegree, iDim;
            double value = 0;

            unsigned short *ijk = new unsigned short[nDim];

            for (iDim = 0; iDim < nDim; iDim++) ijk[iDim] = 0;

            for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
                for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
                    for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
                        ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
                        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] *
                            GetDerivative1(uvw, diff_this, ijk, lmn);
                    }

            delete[] ijk;

            return value;
        }

        double GRID_FreeFormDefBox::GetDerivative4(double *uvw, unsigned short val_diff, unsigned short val_diff2,
            unsigned short *ijk, unsigned short *lmn) {
            unsigned short iDim;
            double value = 0.0;

            if (val_diff == val_diff2) {
                value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 2);
                for (iDim = 0; iDim < nDim; iDim++)
                    if (iDim != val_diff)
                        value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
            }
            else {
                value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 1) *
                    GetBernsteinDerivative(lmn[val_diff2], ijk[val_diff2], uvw[val_diff2], 1);
                for (iDim = 0; iDim < nDim; iDim++)
                    if ((iDim != val_diff) && (iDim != val_diff2))
                        value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
            }

            return value;
        }

        double GRID_FreeFormDefBox::GetDerivative5(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also,
            unsigned short *lmn) {

            unsigned short iDegree, jDegree, kDegree, iDim;
            double value = 0.0;

            unsigned short *ijk = new unsigned short[nDim];

            for (iDim = 0; iDim < nDim; iDim++) ijk[iDim] = 0;

            for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
                for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
                    for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
                        ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
                        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] *
                            GetDerivative4(uvw, diff_this, diff_this_also, ijk, lmn);
                    }

            delete[] ijk;

            return value;
        }
    }
}