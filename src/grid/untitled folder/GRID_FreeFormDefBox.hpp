#ifndef ARIES_GRID_FREEFORMDEFBOX_HPP
#define ARIES_GRID_FREEFORMDEFBOX_HPP


#include "GRID_Gridmovement.hpp"

namespace ARIES
{
    namespace GRID
    {
        /*!
        * \class GRID_FreeFormDefBox
        * \brief Class for defining the free form FFDBox structure.
        * \author F. Palacios & A. Galdran.
        * \version 3.2.9 "eagle"
        */
        class GRID_FreeFormDefBox : public GRID_Gridmovement {
        public:
            unsigned short nDim;                  /*!< \brief Number of dimensions of the problem. */
            unsigned short nCornerPoints,         /*!< \brief Number of corner points of the FFDBox. */
                nControlPoints, nControlPoints_Copy;  /*!< \brief Number of control points of the FFDBox. */
            double **Coord_Corner_Points,		/*!< \brief Coordinates of the corner points. */
                ****Coord_Control_Points,				/*!< \brief Coordinates of the control points. */
                ****ParCoord_Control_Points,		/*!< \brief Coordinates of the control points. */
                ****Coord_Control_Points_Copy,	/*!< \brief Coordinates of the control points (copy). */
                ****Coord_SupportCP;						/*!< \brief Coordinates of the support control points. */
            unsigned short lOrder, lOrder_Copy,	/*!< \brief Order of the FFDBox in the i direction. */
                mOrder, mOrder_Copy, 								/*!< \brief Order of the FFDBox in the j direction. */
                nOrder, nOrder_Copy;									/*!< \brief Order of the FFDBox in the k direction. */
            unsigned short lDegree, lDegree_Copy, /*!< \brief Degree of the FFDBox in the i direction. (lOrder - 1)*/
                mDegree, mDegree_Copy,								/*!< \brief Degree of the FFDBox in the j direction. (mOrder - 1)*/
                nDegree, nDegree_Copy;								/*!< \brief Degree of the FFDBox in the k direction. (nOrder - 1)*/
            double *ParamCoord, *ParamCoord_,	/*!< \brief Parametric coordinates of a point. */
                *cart_coord, *cart_coord_;			/*!< \brief Cartesian coordinates of a point. */
            double ObjFunc;			/*!< \brief Objective function of the point inversion process. */
            double *Gradient;			/*!< \brief Gradient of the point inversion process. */
            double **Hessian;    /*!< \brief Hessian of the point inversion process. */
            double MaxCoord[3];		/*!< \brief Maximum coordinates of the FFDBox. */
            double MinCoord[3];		/*!< \brief Minimum coordinates of the FFDBox. */
            std::string Tag;						/*!< \brief Tag to identify the FFDBox. */
            unsigned short Level;								/*!< \brief Nested level of the FFD box. */
            std::vector<double> CartesianCoord[3];		/*!< \brief Vector with all the cartesian coordinates in the FFD FFDBox. */
            std::vector<double> ParametricCoord[3];	/*!< \brief Vector with all the parametrics coordinates in the FFD FFDBox. */
            std::vector<unsigned short> MarkerIndex;	/*!< \brief Vector with all markers in the FFD FFDBox. */
            std::vector<unsigned long> VertexIndex;	/*!< \brief Vector with all vertex index in the FFD FFDBox. */
            std::vector<unsigned long> PointIndex;		/*!< \brief Vector with all points index in the FFD FFDBox. */
            unsigned long nSurfacePoint;				/*!< \brief Number of surfaces in the FFD FFDBox. */
            std::vector<std::string> ParentFFDBox;					/*!< \brief Vector with all the parent FFD FFDBox. */
            std::vector<std::string> ChildFFDBox;					/*!< \brief Vector with all the child FFD FFDBox. */
            std::vector<unsigned short> Fix_IPlane;  /*!< \brief Fix FFD I plane. */
            std::vector<unsigned short> Fix_JPlane;  /*!< \brief Fix FFD J plane. */
            std::vector<unsigned short> Fix_KPlane;  /*!< \brief Fix FFD K plane. */

        public:

            /*!
            * \brief Constructor of the class.
            */
            GRID_FreeFormDefBox(void);

            /*!
            * \overload
            * \param[in] val_lDegree - Degree of the FFDBox in the i direction.
            * \param[in] val_mDegree - Degree of the FFDBox in the j direction.
            * \param[in] val_nDegree - Degree of the FFDBox in the k direction.
            */
            GRID_FreeFormDefBox(unsigned short val_lDegree, unsigned short val_mDegree, unsigned short val_nDegree);

            /*!
            * \brief Destructor of the class.
            */
            ~GRID_FreeFormDefBox(void);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            void Set_Fix_IPlane(unsigned short val_plane);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            void Set_Fix_JPlane(unsigned short val_plane);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            void Set_Fix_KPlane(unsigned short val_plane);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            unsigned short Get_Fix_IPlane(unsigned short val_index);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            unsigned short Get_Fix_JPlane(unsigned short val_index);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            unsigned short Get_Fix_KPlane(unsigned short val_index);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            unsigned short Get_nFix_IPlane(void);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            unsigned short Get_nFix_JPlane(void);

            /*!
            * \brief Define the I planes to to fix in a FFD box.
            * \param[in] val_plane - Index of the plane to fix.
            */
            unsigned short Get_nFix_KPlane(void);

            /*!
            * \brief Add to the vector of markers a new marker.
            * \param[in] val_iMarker - New marker inside the FFD box.
            */
            void Set_MarkerIndex(unsigned short val_iMarker);

            /*!
            * \brief Add to the vector of vertices a new vertex.
            * \param[in] val_iVertex - New vertex inside the FFD box.
            */
            void Set_VertexIndex(unsigned long val_iVertex);

            /*!
            * \brief Add to the vector of points a new point.
            * \param[in] val_iPoint - New point inside the FFD box.
            */
            void Set_PointIndex(unsigned long val_iPoint);

            /*!
            * \brief Add to the vector of cartesian coordinates a new coordinate.
            * \param[in] val_coord - New coordinate inside the FFD box.
            */
            void Set_CartesianCoord(double *val_coord);

            /*!
            * \brief Add to the vector of parametric coordinates a new coordinate.
            * \param[in] val_coord - New coordinate inside the FFD box.
            */
            void Set_ParametricCoord(double *val_coord);

            /*!
            * \brief Add to the vector of parent FFDBoxes a new FFD FFDBox.
            * \param[in] val_iParentFFDBox - New parent FFDBox in the vector.
            */
            void SetParentFFDBox(std::string val_iParentFFDBox);

            /*!
            * \brief Add to the vector of child FFDBoxes a new FFD FFDBox.
            * \param[in] val_iChildFFDBox - New child FFDBox in the vector.
            */
            void SetChildFFDBox(std::string val_iChildFFDBox);

            /*!
            * \brief _______________.
            * \param[in] val_coord - _______________.
            * \param[in] val_iSurfacePoints - _______________.
            */
            void Set_CartesianCoord(double *val_coord, unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] val_coord - _______________.
            * \param[in] val_iSurfacePoints - _______________.
            */
            void Set_ParametricCoord(double *val_coord, unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] Get_MarkerIndex - _______________.
            * \return _______________.
            */
            unsigned short Get_MarkerIndex(unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] Get_VertexIndex - _______________.
            * \return _______________.
            */
            unsigned long Get_VertexIndex(unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] Get_PointIndex - _______________.
            * \return _______________.
            */
            unsigned long Get_PointIndex(unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] Get_CartesianCoord - _______________.
            * \return _______________.
            */
            double *Get_CartesianCoord(unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] Get_ParametricCoord - _______________.
            * \return _______________.
            */
            double *Get_ParametricCoord(unsigned long val_iSurfacePoints);

            /*!
            * \brief _______________.
            * \param[in] GetnSurfacePoint - _______________.
            * \return _______________.
            */
            unsigned long GetnSurfacePoint(void);

            /*!
            * \brief _______________.
            * \param[in] GetnParentFFDBox - _______________.
            * \return _______________.
            */
            unsigned short GetnParentFFDBox(void);

            /*!
            * \brief _______________.
            * \param[in] GetnChildFFDBox - _______________.
            * \return _______________.
            */
            unsigned short GetnChildFFDBox(void);

            /*!
            * \brief _______________.
            * \param[in] val_ParentFFDBox - _______________.
            * \return _______________.
            */
            std::string GetParentFFDBoxTag(unsigned short val_ParentFFDBox);

            /*!
            * \brief _______________.
            * \param[in] val_ChildFFDBox - _______________.
            * \return _______________.
            */
            std::string GetChildFFDBoxTag(unsigned short val_ChildFFDBox);

            /*!
            * \brief Change the the position of the corners of the unitary FFDBox,
            *        and find the position of the control points for the FFDBox
            * \param[in] FFDBox - Original FFDBox where we want to compute the control points.
            */
            void SetSupportCPChange(GRID_FreeFormDefBox *FFDBox);

            /*!
            * \brief Set the number of corner points.
            * \param[in] val_ncornerpoints - Number of corner points.
            */
            void SetnCornerPoints(unsigned short val_ncornerpoints);

            /*!
            * \brief Get the number of corner points.
            * \return Number of corner points.
            */
            unsigned short GetnCornerPoints(void);

            /*!
            * \brief Get the number of control points.
            * \return Number of control points.
            */
            unsigned short GetnControlPoints(void);

            /*!
            * \brief Get the number of control points.
            * \return Number of control points.
            */
            void SetnControlPoints(void);

            /*!
            * \brief Get the number of numerical points on the surface.
            * \return Number of numerical points on the surface.
            */
            unsigned long GetnSurfacePoints(void);

            /*!
            * \brief Set the corner point for the unitary FFDBox.
            */
            void SetUnitCornerPoints(void);

            /*!
            * \brief Set the coordinates of the corner points.
            * \param[in] val_coord - Coordinates of the corner point with index <i>val_icornerpoints</i>.
            * \param[in] val_icornerpoints - Index of the corner point.
            */
            void SetCoordCornerPoints(double *val_coord, unsigned short val_icornerpoints);

            /*!
            * \overload
            * \param[in] val_xcoord - X coordinate of the corner point with index <i>val_icornerpoints</i>.
            * \param[in] val_ycoord - Y coordinate of the corner point with index <i>val_icornerpoints</i>.
            * \param[in] val_zcoord - Z coordinate of the corner point with index <i>val_icornerpoints</i>.
            * \param[in] val_icornerpoints - Index of the corner point.
            */
            void SetCoordCornerPoints(double val_xcoord, double val_ycoord, double val_zcoord, unsigned short val_icornerpoints);

            /*!
            * \brief Set the coordinates of the control points.
            * \param[in] val_coord - Coordinates of the control point.
            * \param[in] iDegree - Index of the FFDBox, i direction.
            * \param[in] jDegree - Index of the FFDBox, j direction.
            * \param[in] kDegree - Index of the FFDBox, k direction.
            */
            void SetCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree);

            /*!
            * \brief Set the coordinates of the control points.
            * \param[in] val_coord - Coordinates of the control point.
            * \param[in] iDegree - Index of the FFDBox, i direction.
            * \param[in] jDegree - Index of the FFDBox, j direction.
            * \param[in] kDegree - Index of the FFDBox, k direction.
            */
            void SetCoordControlPoints_Copy(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree);

            /*!
            * \brief Set the coordinates of the control points.
            * \param[in] val_coord - Coordinates of the control point.
            * \param[in] iDegree - Index of the FFDBox, i direction.
            * \param[in] jDegree - Index of the FFDBox, j direction.
            * \param[in] kDegree - Index of the FFDBox, k direction.
            */
            void SetParCoordControlPoints(double *val_coord, unsigned short iDegree, unsigned short jDegree, unsigned short kDegree);

            /*!
            * \brief Get the coordinates of the corner points.
            * \param[in] val_dim - Index of the coordinate (x, y, z).
            * \param[in] val_icornerpoints - Index of the corner point.
            * \return Coordinate <i>val_dim</i> of the corner point <i>val_icornerpoints</i>.
            */
            double GetCoordCornerPoints(unsigned short val_dim, unsigned short val_icornerpoints);

            /*!
            * \brief Get the coordinates of the corner points.
            * \param[in] val_icornerpoints - Index of the corner point.
            * \return Pointer to the coordinate vector of the corner point <i>val_icornerpoints</i>.
            */
            double *GetCoordCornerPoints(unsigned short val_icornerpoints);

            /*!
            * \brief Get the coordinates of the control point.
            * \param[in] val_iindex - Value of the local i index of the control point.
            * \param[in] val_jindex - Value of the local j index of the control point.
            * \param[in] val_kindex - Value of the local k index of the control point.
            * \return Pointer to the coordinate vector of the control point with local index (i, j, k).
            */
            double *GetCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex);

            /*!
            * \brief Get the parametric coordinates of the control point.
            * \param[in] val_iindex - Value of the local i index of the control point.
            * \param[in] val_jindex - Value of the local j index of the control point.
            * \param[in] val_kindex - Value of the local k index of the control point.
            * \return Pointer to the coordinate vector of the control point with local index (i, j, k).
            */
            double *GetParCoordControlPoints(unsigned short val_iindex, unsigned short val_jindex, unsigned short val_kindex);

            /*!
            * \brief Set the control points in a parallelepiped (hexahedron).
            */
            void SetControlPoints_Parallelepiped(void);

            /*!
            * \brief Set the control points of the final chuck in a unitary hexahedron free form.
            * \param[in] FFDBox - Original FFDBox where we want to compute the control points.
            */
            void SetSupportCP(GRID_FreeFormDefBox *FFDBox);

            /*!
            * \brief Set the new value of the coordinates of the control points.
            * \param[in] val_index - Local index (i, j, k) of the control point.
            * \param[in] movement - Movement of the control point.
            */
            void SetControlPoints(unsigned short *val_index, double *movement);

            /*!
            * \brief Set the original value of the control points.
            */
            void SetOriginalControlPoints(void);

            /*!
            * \brief Set the tecplot file of the FFD chuck structure.
            * \param[in] iFFDBox - Index of the FFD box.
            * \param[in] original - Original box (before deformation).
            */
            void SetTecplot(GEOM::GEOM_Geometry *geometry, unsigned short iFFDBox, bool original);

            /*!
            * \brief Set the cartesian coords of a point in R^3 and convert them to the parametric coords of
            *        our parametrization of a paralellepiped.
            * \param[in] cart_coord - Cartesian coordinates of a point.
            * \return Pointer to the parametric coordinates of a point.
            */
            double *GetParametricCoord_Analytical(double *cart_coord);

            /*!
            * \brief Iterative strategy for computing the parametric coordinates.
            * \param[in] xyz - Cartesians coordinates of the target point.
            * \param[in] guess - Initial guess for doing the parametric coordinates search.
            * \param[in] tol - Level of convergence of the iterative method.
            * \param[in] it_max - Maximal number of iterations.
            * \return Parametric coordinates of the point.
            */
            double *GetParametricCoord_Iterative(unsigned long iPoint, double *xyz, double *guess, TBOX::TBOX_Config *config);

            /*!
            * \brief Compute the cross product.
            * \param[in] v1 - First input vector.
            * \param[in] v2 - Second input vector.
            * \param[out] v3 - Output vector wuth the cross product.
            */
            void CrossProduct(double *v1, double *v2, double *v3);

            /*!
            * \brief Compute the doc product.
            * \param[in] v1 - First input vector.
            * \param[in] v2 - Sencond input vector.
            * \return Dot product between <i>v1</i>, and <i>v2</i>.
            */
            double DotProduct(double *v1, double *v2);

            /*!
            * \brief Here we take the parametric coords of a point in the box and we convert them to the
            *        physical cartesian coords by plugging the ParamCoords on the Bezier parameterization of our box.
            * \param[in] ParamCoord - Parametric coordinates of a point.
            * \return Pointer to the cartesian coordinates of a point.
            */
            double *EvalCartesianCoord(double *ParamCoord);

            /*!
            * \brief Set the Bernstein polynomial, defined as B_i^n(t) = Binomial(n, i)*t^i*(1-t)^(n-i).
            * \param[in] val_n - Degree of the Bernstein polynomial.
            * \param[in] val_i - Order of the Bernstein polynomial.
            * \param[in] val_t - Value of the parameter where the polynomial is evaluated.
            * \return Value of the Bernstein polynomial.
            */
            double GetBernstein(short val_n, short val_i, double val_t);

            /*!
            * \brief Get the binomial coefficient n over i, defined as n!/(m!(n-m)!)
            * \note If the denominator is 0, the value is 1.
            * \param[in] n - Upper coefficient.
            * \param[in] m - Lower coefficient.
            * \return Value of the binomial coefficient n over m.
            */
            unsigned long Binomial(unsigned short n, unsigned short m);

            /*!
            * \brief Get the order in the l direction of the FFD FFDBox.
            * \return Order in the l direction of the FFD FFDBox.
            */
            unsigned short GetlOrder(void);

            /*!
            * \brief Get the order in the m direction of the FFD FFDBox.
            * \return Order in the m direction of the FFD FFDBox.
            */
            unsigned short GetmOrder(void);

            /*!
            * \brief Get the order in the n direction of the FFD FFDBox.
            * \return Order in the n direction of the FFD FFDBox.
            */
            unsigned short GetnOrder(void);

            /*!
            * \brief Get the order in the l direction of the FFD FFDBox.
            * \return Order in the l direction of the FFD FFDBox.
            */
            void SetlOrder(unsigned short val_lOrder);

            /*!
            * \brief Get the order in the m direction of the FFD FFDBox.
            * \return Order in the m direction of the FFD FFDBox.
            */
            void SetmOrder(unsigned short val_mOrder);

            /*!
            * \brief Get the order in the n direction of the FFD FFDBox.
            * \return Order in the n direction of the FFD FFDBox.
            */
            void SetnOrder(unsigned short val_nOrder);

            /*!
            * \brief Set, at each vertex, the index of the free form FFDBox that contains the vertex.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iFFDBox - Index of the FFDBox.
            */
            bool GetPointFFD(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned long iPoint);

            /*!
            * \brief Set the zone of the computational domain that is going to be deformed.
            * \param[in] geometry - Geometrical definition of the problem.
            * \param[in] config - Definition of the particular problem.
            * \param[in] iFFDBox - Index of the FFDBox.
            */
            void SetDeformationZone(GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config, unsigned short iFFDBox);

            /*!
            * \brief The "order" derivative of the i-th Bernstein polynomial of degree n, evaluated at t,
            *        is calculated as  (B_i^n(t))^{order}(t) = n*(GetBernstein(n-1, i-1, t)-GetBernstein(n-1, i, t)),
            *        having in account that if i=0, GetBernstein(n-1,-1, t) = 0.
            * \param[in] val_n - Degree of the Bernstein polynomial.
            * \param[in] val_i - Order of the Bernstein polynomial.
            * \param[in] val_t - Value of the parameter where the polynomial is evaluated.
            * \param[in] val_order - Order of the derivative.
            * \return Value of the Derivative of the Bernstein polynomial.
            */
            double GetBernsteinDerivative(short val_n, short val_i, double val_t, short val_order);

            /*!
            * \brief The routine computes the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2  evaluated at (u, v, w).
            * \param[in] val_coord - Parametric coordiates of the target point.
            * \param[in] xyz - Cartesians coordinates of the point.
            * \param[in] analytical - Compute the analytical gradient.
            * \return Value of the analytical gradient.
            */
            double *GetFFDGradient(double *val_coord, double *xyz);

            /*!
            * \brief The routine that computes the Hessian of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 evaluated at (u, v, w)
            *        Input: (u, v, w), (x, y, z)
            *        Output: Hessian F (u, v, w).
            * \param[in] uvw - Current value of the parametrics coordinates.
            * \param[in] xyz - Cartesians coordinates of the target point to compose the functional.
            * \param[in] val_Hessian - Value of the hessian.
            */
            void GetFFDHessian(double *uvw, double *xyz, double **val_Hessian);

            /*!
            * \brief An auxiliary routine to help us compute the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
            *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)^2+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)^2
            *        Input: val_t, val_diff (to identify the index of the Bernstein polynomail we differentiate), the i, j, k , l, m, n
            *        E.G.: val_diff=2 => we differentiate w.r.t. w  (val_diff=0,1, or 2) Output: d [B_i^l*B_j^m *B_k^n] / d val_diff
            *        (val_u, val_v, val_w).
            * \param[in] uvw - __________.
            * \param[in] val_diff - __________.
            * \param[in] ijk - __________.
            * \param[in] lmn - Degree of the FFD box.
            * \return __________.
            */
            double GetDerivative1(double *uvw, unsigned short val_diff, unsigned short *ijk, unsigned short *lmn);

            /*!
            * \brief An auxiliary routine to help us compute the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
            *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)^2+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)^2
            *        Input: (u, v, w), dim , xyz=(x, y, z), l, m, n E.G.: dim=2 => we use the third coordinate of the control points,
            *        and the z-coordinate of xyz  (0<=dim<=2) Output: 2* ( (Sum_{i, j, k}^l, m, n P_{ijk}[dim] B_i^l[u] B_j^m[v] B_k^n[w]) -
            *        xyz[dim]).
            * \param[in] uvw - __________.
            * \param[in] dim - __________.
            * \param[in] xyz - __________.
            * \param[in] lmn - Degree of the FFD box.
            * \return __________.
            */
            double GetDerivative2(double *uvw, unsigned short dim, double *xyz, unsigned short *lmn);

            /*!
            * \brief An auxiliary routine to help us compute the gradient of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
            *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
            * \param[in] uvw - Parametric coordiates of the point.
            * \param[in] dim - Value of the coordinate to be differentiate.
            * \param[in] diff_this - Diferentiation with respect this coordinate.
            * \param[in] lmn - Degree of the FFD box.
            * \return Sum_{i, j, k}^{l, m, n} [one of them with -1,
            *        depending on diff_this=0,1 or 2] P_{ijk}[dim] * (B_i^l[u] B_j^m[v] B_k^n[w])--one of them diffrentiated;
            *        which? diff_thiss will tell us ; E.G.: dim=2, diff_this=1 => we use the third coordinate of the control
            *        points, and derivate de v-Bersntein polynomial (use m-1 when summing!!).
            */
            double GetDerivative3(double *uvw, unsigned short dim, unsigned short diff_this,
                unsigned short *lmn);

            /*!
            * \brief An auxiliary routine to help us compute the Hessian of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
            *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
            *        Input: val_t, val_diff, val_diff2 (to identify the index of the Bernstein polynomials we differentiate), the i, j, k , l, m, n
            *        E.G.: val_diff=1, val_diff2=2  =>  we differentiate w.r.t. v and w  (val_diff=0,1, or 2)
            *        E.G.: val_diff=0, val_diff2=0 => we differentiate w.r.t. u two times
            *        Output: [d [B_i^l*B_j^m *B_k^n]/d val_diff *d [B_i^l*B_j^m *B_k^n]/d val_diff2] (val_u, val_v, val_w) .
            * \param[in] uvw - __________.
            * \param[in] val_diff - __________.
            * \param[in] val_diff2 - __________.
            * \param[in] ijk - __________.
            * \param[in] lmn - Degree of the FFD box.
            * \return __________.
            */
            double GetDerivative4(double *uvw, unsigned short val_diff, unsigned short val_diff2,
                unsigned short *ijk, unsigned short *lmn);

            /*!
            * \brief An auxiliary routine to help us compute the Hessian of F(u, v, w) = ||X(u, v, w)-(x, y, z)||^2 =
            *        (Sum_ijk^lmn P1_ijk Bi Bj Bk -x)^2+(Sum_ijk^lmn P2_ijk Bi Bj Bk -y)+(Sum_ijk^lmn P3_ijk Bi Bj Bk -z)
            *        Input: (u, v, w), dim , diff_this, diff_this_also, xyz=(x, y, z), l, m, n
            *        Output:
            *        Sum_{i, j, k}^{l, m, n} [two of them with -1, depending on diff_this, diff_this_also=0,1 or 2]
            *        P_{ijk}[dim] * (B_i^l[u] B_j^m[v] B_k^n[w])--one of them diffrentiated; which? diff_thiss will tell us ;
            *        E.G.: dim=2, diff_this=1 => we use the third coordinate of the control points, and derivate de v-Bersntein
            *        polynomial (use m-1 when summing!!).
            * \param[in] uvw - __________.
            * \param[in] dim - __________.
            * \param[in] diff_this - __________.
            * \param[in] diff_this_also - __________.
            * \param[in] lmn - Degree of the FFD box.
            * \return __________.
            */
            double GetDerivative5(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also,
                unsigned short *lmn);

            /*!
            * \brief Euclidean norm of a vector.
            * \param[in] a - _______.
            * \return __________.
            */
            double GetNorm(double *a);

            /*!
            * \brief Set the tag that identify a FFDBox.
            * \param[in] val_tag - value of the tag.
            */
            void SetTag(std::string val_tag);

            /*!
            * \brief Get the tag that identify a FFDBox.
            * \return Value of the tag that identigy the FFDBox.
            */
            std::string GetTag(void);

            /*!
            * \brief Set the nested level of the FFDBox.
            * \param[in] val_level - value of the level.
            */
            void SetLevel(unsigned short val_level);

            /*!
            * \brief Get the nested level of the FFDBox.
            * \return Value of the nested level of the the FFDBox.
            */
            unsigned short GetLevel(void);

            /*!
            * \brief Compute the determinant of a 3 by 3 matrix.
            * \param[in] val_matrix 3 by 3 matrix.
            * \result Determinant of the matrix
            */
            double Determinant_3x3(double A00, double A01, double A02, double A10, double A11,
                double A12, double A20, double A21, double A22);

        };
    }
}

#include "GRID_FreeFormDefBox.inl"

#endif