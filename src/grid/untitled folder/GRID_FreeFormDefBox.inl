#ifndef ARIES_GRID_FREEFORMDEFBOX_INLINE
#define ARIES_GRID_FREEFORMDEFBOX_INLINE

namespace ARIES
{
    namespace GRID
    {

        inline void GRID_FreeFormDefBox::Set_Fix_IPlane(unsigned short val_plane) { Fix_IPlane.push_back(val_plane); }

        inline void GRID_FreeFormDefBox::Set_Fix_JPlane(unsigned short val_plane) { Fix_JPlane.push_back(val_plane); }

        inline void GRID_FreeFormDefBox::Set_Fix_KPlane(unsigned short val_plane) { Fix_KPlane.push_back(val_plane); }

        inline unsigned short GRID_FreeFormDefBox::Get_Fix_IPlane(unsigned short val_index) { return Fix_IPlane[val_index]; }

        inline unsigned short GRID_FreeFormDefBox::Get_Fix_JPlane(unsigned short val_index) { return Fix_JPlane[val_index]; }

        inline unsigned short GRID_FreeFormDefBox::Get_Fix_KPlane(unsigned short val_index) { return Fix_KPlane[val_index]; }

        inline unsigned short GRID_FreeFormDefBox::Get_nFix_IPlane(void) { return Fix_IPlane.size(); }

        inline unsigned short GRID_FreeFormDefBox::Get_nFix_JPlane(void) { return Fix_JPlane.size(); }

        inline unsigned short GRID_FreeFormDefBox::Get_nFix_KPlane(void) { return Fix_KPlane.size(); }

        inline void GRID_FreeFormDefBox::Set_MarkerIndex(unsigned short val_iMarker) { MarkerIndex.push_back(val_iMarker); }

        inline void GRID_FreeFormDefBox::SetParentFFDBox(std::string val_iParentFFDBox) { ParentFFDBox.push_back(val_iParentFFDBox); }

        inline void GRID_FreeFormDefBox::SetChildFFDBox(std::string val_iChildFFDBox) { ChildFFDBox.push_back(val_iChildFFDBox); }

        inline void GRID_FreeFormDefBox::Set_VertexIndex(unsigned long val_iVertex) { VertexIndex.push_back(val_iVertex); }

        inline void GRID_FreeFormDefBox::Set_PointIndex(unsigned long val_iPoint) { PointIndex.push_back(val_iPoint); }

        inline void GRID_FreeFormDefBox::Set_CartesianCoord(double *val_coord) {
            CartesianCoord[0].push_back(val_coord[0]);
            CartesianCoord[1].push_back(val_coord[1]);
            CartesianCoord[2].push_back(val_coord[2]);
        }

        inline void GRID_FreeFormDefBox::Set_CartesianCoord(double *val_coord, unsigned long val_iSurfacePoints) {
            CartesianCoord[0][val_iSurfacePoints] = val_coord[0];
            CartesianCoord[1][val_iSurfacePoints] = val_coord[1];
            CartesianCoord[2][val_iSurfacePoints] = val_coord[2];
        }

        inline void GRID_FreeFormDefBox::Set_ParametricCoord(double *val_coord) {
            ParametricCoord[0].push_back(val_coord[0]);
            ParametricCoord[1].push_back(val_coord[1]);
            ParametricCoord[2].push_back(val_coord[2]);
        }

        inline void GRID_FreeFormDefBox::Set_ParametricCoord(double *val_coord, unsigned long val_iSurfacePoints) {
            ParametricCoord[0][val_iSurfacePoints] = val_coord[0];
            ParametricCoord[1][val_iSurfacePoints] = val_coord[1];
            ParametricCoord[2][val_iSurfacePoints] = val_coord[2];
        }

        inline unsigned short GRID_FreeFormDefBox::Get_MarkerIndex(unsigned long val_iSurfacePoints) { return MarkerIndex[val_iSurfacePoints]; }

        inline unsigned short GRID_FreeFormDefBox::GetnParentFFDBox(void) { return ParentFFDBox.size(); }

        inline unsigned short GRID_FreeFormDefBox::GetnChildFFDBox(void) { return ChildFFDBox.size(); }

        inline std::string GRID_FreeFormDefBox::GetParentFFDBoxTag(unsigned short val_ParentFFDBox) { return ParentFFDBox[val_ParentFFDBox]; }

        inline std::string GRID_FreeFormDefBox::GetChildFFDBoxTag(unsigned short val_ChildFFDBox) { return ChildFFDBox[val_ChildFFDBox]; }

        inline unsigned long GRID_FreeFormDefBox::Get_VertexIndex(unsigned long val_iSurfacePoints) { return VertexIndex[val_iSurfacePoints]; }

        inline unsigned long GRID_FreeFormDefBox::Get_PointIndex(unsigned long val_iSurfacePoints) { return PointIndex[val_iSurfacePoints]; }

        inline double *GRID_FreeFormDefBox::Get_CartesianCoord(unsigned long val_iSurfacePoints) {
            cart_coord_[0] = CartesianCoord[0][val_iSurfacePoints];
            cart_coord_[1] = CartesianCoord[1][val_iSurfacePoints];
            cart_coord_[2] = CartesianCoord[2][val_iSurfacePoints];
            return cart_coord_;
        }

        inline double *GRID_FreeFormDefBox::Get_ParametricCoord(unsigned long val_iSurfacePoints) {
            ParamCoord_[0] = ParametricCoord[0][val_iSurfacePoints];
            ParamCoord_[1] = ParametricCoord[1][val_iSurfacePoints];
            ParamCoord_[2] = ParametricCoord[2][val_iSurfacePoints];
            return ParamCoord_;
        }

        inline unsigned long GRID_FreeFormDefBox::GetnSurfacePoint(void) { return PointIndex.size(); }

        inline void GRID_FreeFormDefBox::SetnCornerPoints(unsigned short val_ncornerpoints) { nCornerPoints = val_ncornerpoints; }

        inline unsigned short GRID_FreeFormDefBox::GetnCornerPoints(void) { return nCornerPoints; }

        inline unsigned short GRID_FreeFormDefBox::GetnControlPoints(void) { return nControlPoints; }

        inline void GRID_FreeFormDefBox::SetnControlPoints(void) { nControlPoints = lOrder*mOrder*nOrder; }

        inline unsigned long GRID_FreeFormDefBox::GetnSurfacePoints(void) { return 0; }

        inline double *GRID_FreeFormDefBox::GetCoordCornerPoints(unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints]; }

        inline double *GRID_FreeFormDefBox::GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return Coord_Control_Points[val_iindex][val_jindex][val_kindex]; }

        inline double *GRID_FreeFormDefBox::GetParCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex) { return ParCoord_Control_Points[val_iindex][val_jindex][val_kindex]; }

        inline double  GRID_FreeFormDefBox::GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints) { return Coord_Corner_Points[val_icornerpoints][val_dim]; }

        inline unsigned short GRID_FreeFormDefBox::GetlOrder(void) { return lOrder; }

        inline unsigned short GRID_FreeFormDefBox::GetmOrder(void) { return mOrder; }

        inline unsigned short GRID_FreeFormDefBox::GetnOrder(void) { return nOrder; }

        inline void GRID_FreeFormDefBox::SetlOrder(unsigned short val_lOrder) { lOrder = val_lOrder; lDegree = lOrder - 1; }

        inline void GRID_FreeFormDefBox::SetmOrder(unsigned short val_mOrder) { mOrder = val_mOrder; mDegree = mOrder - 1; }

        inline void GRID_FreeFormDefBox::SetnOrder(unsigned short val_nOrder) { nOrder = val_nOrder; nDegree = nOrder - 1; }

        inline void  GRID_FreeFormDefBox::SetCoordCornerPoints(double *val_coord, unsigned short val_icornerpoints) {
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
                Coord_Corner_Points[val_icornerpoints][iDim] = val_coord[iDim];
        }

        inline void GRID_FreeFormDefBox::SetCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
            for (unsigned short iDim = 0; iDim < nDim; iDim++) {
                Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
            }
        }

        inline void GRID_FreeFormDefBox::SetCoordControlPoints_Copy(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
            for (unsigned short iDim = 0; iDim < nDim; iDim++) {
                Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
            }
        }

        inline void GRID_FreeFormDefBox::SetParCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree) {
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
                ParCoord_Control_Points[iDegree][jDegree][kDegree][iDim] = val_coord[iDim];
        }

        inline void GRID_FreeFormDefBox::SetCoordCornerPoints(double val_xcoord, double val_ycoord, double val_zcoord, unsigned short val_icornerpoints) {
            Coord_Corner_Points[val_icornerpoints][0] = val_xcoord;
            Coord_Corner_Points[val_icornerpoints][1] = val_ycoord;
            Coord_Corner_Points[val_icornerpoints][2] = val_zcoord;
        }

        inline void GRID_FreeFormDefBox::SetControlPoints(unsigned short *val_index, double *movement) {
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
                Coord_Control_Points[val_index[0]][val_index[1]][val_index[2]][iDim] += movement[iDim];
        }

        inline void GRID_FreeFormDefBox::SetOriginalControlPoints() {
            for (unsigned short iDegree = 0; iDegree <= lDegree_Copy; iDegree++)
                for (unsigned short jDegree = 0; jDegree <= mDegree_Copy; jDegree++)
                    for (unsigned short kDegree = 0; kDegree <= nDegree_Copy; kDegree++)
                        for (unsigned short iDim = 0; iDim < nDim; iDim++)
                            Coord_Control_Points[iDegree][jDegree][kDegree][iDim] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][iDim];

            lDegree = lDegree_Copy; mDegree = mDegree_Copy; nDegree = nDegree_Copy;
            lOrder = lOrder_Copy; mOrder = mOrder_Copy; nOrder = nOrder_Copy;
            nControlPoints = nControlPoints_Copy;
        }

        inline void GRID_FreeFormDefBox::CrossProduct(double *v1, double *v2, double *v3) {
            v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
            v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
            v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
        }

        inline double GRID_FreeFormDefBox::DotProduct(double *v1, double *v2) { double scalar = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]; return scalar; }

        inline double GRID_FreeFormDefBox::GetNorm(double *a) { double  norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); return norm; }

        inline void GRID_FreeFormDefBox::SetTag(std::string val_tag) { Tag = val_tag; }

        inline std::string GRID_FreeFormDefBox::GetTag() { return Tag; }

        inline void GRID_FreeFormDefBox::SetLevel(unsigned short val_level) { Level = val_level; }

        inline unsigned short GRID_FreeFormDefBox::GetLevel() { return Level; }

        inline double GRID_FreeFormDefBox::Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22) {
            return A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20);
        }
    }
}

#endif