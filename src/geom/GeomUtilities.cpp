/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Geometry utilities
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "GeomUtilities.hpp"

#include "const_def.h"
#include "armath.h"

#include <cmath>


namespace ARIES
{
    double GeomUtilities::Point2Plane_Distance(double *coord, double *iCoord, double *jCoord, double *kCoord)
    {
        double crossProduct[3], iVector[3], jVector[3], distance, modulus;
        unsigned short iDim;

        for (iDim = 0; iDim < 3; iDim++)
        {
            iVector[iDim] = jCoord[iDim] - iCoord[iDim];
            jVector[iDim] = kCoord[iDim] - iCoord[iDim];
        }

        crossProduct[0] = iVector[1] * jVector[2] - iVector[2] * jVector[1];
        crossProduct[1] = iVector[2] * jVector[0] - iVector[0] * jVector[2];
        crossProduct[2] = iVector[0] * jVector[1] - iVector[1] * jVector[0];

        modulus = sqrt(crossProduct[0] * crossProduct[0] + crossProduct[1] * crossProduct[1] + crossProduct[2] * crossProduct[2]);

        distance = 0.0;
        for (iDim = 0; iDim < 3; iDim++)
            distance += crossProduct[iDim] * (coord[iDim] - iCoord[iDim]);
        distance /= modulus;

        return distance;
    }

    bool GeomUtilities::RayIntersectsTriangle(double orig[3], double dir[3], double vert0[3], double vert1[3], double vert2[3], double *intersect)
    {
        double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
        double det, inv_det, t, u, v;

        // Find vectors for two edges sharing vert0
        SUB(edge1, vert1, vert0);
        SUB(edge2, vert2, vert0);

        // Begin calculating determinant - also used to calculate U parameter
        CROSS(pvec, dir, edge2);

        // If determinant is near zero, ray lies in plane of triangle
        det = DOT(edge1, pvec);

        if (det > -EPSILON && det < EPSILON)
            return true;
        inv_det = 1.0 / det;

        // Calculate distance from vert0 to ray origin
        SUB(tvec, orig, vert0);

        // Calculate U parameter and test bounds
        u = inv_det * DOT(tvec, pvec);

        if (u < 0.0 || u > 1.0)
            return false;

        // prepare to test V parameter
        CROSS(qvec, tvec, edge1);

        // Calculate V parameter and test bounds
        v = inv_det * DOT(dir, qvec);

        if (v < 0.0 || u + v > 1.0)
            return false;

        // Calculate t, ray intersects triangle
        t = inv_det * DOT(edge2, qvec);

        // Compute the intersection point in cartesian coordinates
        intersect[0] = orig[0] + (t * dir[0]);
        intersect[1] = orig[1] + (t * dir[1]);
        intersect[2] = orig[2] + (t * dir[2]);

        return true;
    }

    bool GeomUtilities::SegmentIntersectsTriangle(double point0[3], double point1[3], double vert0[3], double vert1[3], double vert2[3])
    {
        double dir[3], intersect[3], u[3], v[3], edge1[3], edge2[3], planeNormal[3], denominator, numerator, aux;

        SUB(dir, point1, point0);

        if (RayIntersectsTriangle(point0, dir, vert0, vert1, vert2, intersect))
        {
            // Check that the intersection is in the segment
            SUB(u, point0, intersect);
            SUB(v, point1, intersect);

            SUB(edge1, vert1, vert0);
            SUB(edge2, vert2, vert0);
            CROSS(planeNormal, edge1, edge2);

            denominator = DOT(planeNormal, u);
            numerator   = DOT(planeNormal, v);

            aux = numerator * denominator;

            // Intersection outside the segment
            if (aux > 0.0)
                return false;
        }
        else
        {
            // No intersection with the ray
            return false;
        }

        // Intersection inside the segment
        return true;
    }
    
    bool GeomUtilities::SegmentIntersectsPlane(double *Segment_P0, double *Segment_P1, double Variable_P0, double Variable_P1,
                                               double *Plane_P0, double *Plane_Normal, double *Intersection, double &Variable_Interp)
    {
        double u[3], v[3], Denominator, Numerator, Aux, ModU;
        unsigned short iDim;

        for (iDim = 0; iDim < 3; iDim++)
        {
            u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
            v[iDim] =   Plane_P0[iDim] - Segment_P0[iDim];
        }

        ModU = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

        Numerator   = Plane_Normal[0] * v[0] + Plane_Normal[1] * v[1] + Plane_Normal[2] * v[2];
        Denominator = Plane_Normal[0] * u[0] + Plane_Normal[1] * u[1] + Plane_Normal[2] * u[2];

        if (fabs(Denominator) <= 0.0)
            return false; // No intersection.

        Aux = Numerator / Denominator;

        if (Aux < 0.0 || Aux > 1.0)
            return false; // No intersection.

        for (iDim = 0; iDim < 3; iDim++)
            Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];

        /*--- Check that the intersection is in the segment ---*/
        for (iDim = 0; iDim < 3; iDim++)
        {
            u[iDim] = Segment_P0[iDim] - Intersection[iDim];
            v[iDim] = Segment_P1[iDim] - Intersection[iDim];
        }

        Variable_Interp = Variable_P0 + (Variable_P1 - Variable_P0)*sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) / ModU;

        Denominator = Plane_Normal[0] * u[0] + Plane_Normal[1] * u[1] + Plane_Normal[2] * u[2];
        Numerator   = Plane_Normal[0] * v[0] + Plane_Normal[1] * v[1] + Plane_Normal[2] * v[2];

        Aux = Numerator * Denominator;

        if (Aux > 0.0)
            return false; // Intersection outside the segment.

        return true;
    }

    void GeomUtilities::SetSpline(std::vector<double> &x, std::vector<double> &y, unsigned long n, double yp1, double ypn, std::vector<double> &y2)
    {
        unsigned long i, k;
        double p, qn, sig, un, *u;

        u = new double[n];

        if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
            y2[0] = u[0] = 0.0;			  // -ural"
        else
        {									// or else to have a specified first derivative.
            y2[0] = -0.5;
            u[0] = (3.0 / (x[1] - x[0]))*((y[1] - y[0]) / (x[1] - x[0]) - yp1);
        }

        for (i = 2; i <= n - 1; i++)
        {									//  This is the decomposition loop of the tridiagonal al-
            sig = (x[i - 1] - x[i - 2]) / (x[i] - x[i - 2]);		//	gorithm. y2 and u are used for tem-
            p = sig*y2[i - 2] + 2.0;										//	porary storage of the decomposed
            y2[i - 1] = (sig - 1.0) / p;										//	factors.

            double a1 = (y[i] - y[i - 1]) / (x[i] - x[i - 1]); if (x[i] == x[i - 1]) a1 = 1.0;
            double a2 = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]); if (x[i - 1] == x[i - 2]) a2 = 1.0;
            u[i - 1] = a1 - a2;
            u[i - 1] = (6.0*u[i - 1] / (x[i] - x[i - 2]) - sig*u[i - 2]) / p;
        }

        if (ypn > 0.99e30)						// The upper boundary condition is set either to be
            qn = un = 0.0;									// "natural"
        else
        {												// or else to have a specified first derivative.
            qn = 0.5;
            un = (3.0 / (x[n - 1] - x[n - 2]))*(ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
        }
        y2[n - 1] = (un - qn*u[n - 2]) / (qn*y2[n - 2] + 1.0);
        for (k = n - 1; k >= 1; k--)					// This is the backsubstitution loop of the tridiagonal
            y2[k - 1] = y2[k - 1] * y2[k] + u[k - 1];	  // algorithm.

        delete[] u;
    }

    double GeomUtilities::GetSpline(std::vector<double>&xa, std::vector<double>&ya, std::vector<double>&y2a, unsigned long n, double x)
    {
        unsigned long klo, khi, k;
        double h, b, a, y;

        if (x < xa[0]) x = xa[0];       // Clip max and min values
        if (x > xa[n - 1]) x = xa[n - 1];

        klo = 1;										// We will find the right place in the table by means of
        khi = n;										// bisection. This is optimal if sequential calls to this
        while (khi - klo > 1)
        {			                          // routine are at random values of x. If sequential calls
            k = (khi + klo) >> 1;				// are in order, and closely spaced, one would do better
            if (xa[k - 1] > x)
                khi = k;		// to store previous values of klo and khi and test if
            else
                klo = k;							// they remain appropriate on the next call.
        }								// klo and khi now bracket the input value of x
        h = xa[khi - 1] - xa[klo - 1];
        if (h == 0.0)
            h = EPS; // cout << "Bad xa input to routine splint" << endl;	// The xa’s must be distinct.
        a = (xa[khi - 1] - x) / h;
        b = (x - xa[klo - 1]) / h;				// Cubic spline polynomial is now evaluated.
        y = a*ya[klo - 1] + b*ya[khi - 1] + ((a*a*a - a)*y2a[klo - 1] + (b*b*b - b)*y2a[khi - 1])*(h*h) / 6.0;

        return y;
    }
    
    void GeomUtilities::ComputeAirfoilSection(double *Plane_P0, double *Plane_Normal,
                                              double MinXCoord, double MaxXCoord, double *FlowVariable,
                                              std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil,
                                              std::vector<double> &Zcoord_Airfoil, std::vector<double> &Variable_Airfoil,
                                              bool original_surface, TBOX::TBOX_Config *config)
    {
        unsigned short iMarker, iNode, jNode, iDim;
        bool intersect;
        long MinDist_Point, MinDistAngle_Point;
        unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, jVertex;
        double Segment_P0[3] = { 0.0, 0.0, 0.0 }, Segment_P1[3] = { 0.0, 0.0, 0.0 }, Variable_P0 = 0.0, Variable_P1 = 0.0, Intersection[3] = { 0.0, 0.0, 0.0 }, Trailing_Coord, MinDist_Value, MinDistAngle_Value, Dist_Value,
                                                                                                                                                                        Airfoil_Tangent[3] = { 0.0, 0.0, 0.0 }, Segment[3] = { 0.0, 0.0, 0.0 }, Length, Angle_Value, MaxAngle = 30, *VarCoord = NULL, CosValue, Variable_Interp;
        std::vector<double> Xcoord, Ycoord, Zcoord, Variable;
        std::vector<unsigned long> Duplicate;
        std::vector<unsigned long>::iterator it;
        int rank = TBOX::MASTER_NODE;
        double **Coord_Variation = NULL;

#ifdef HAVE_MPI
        unsigned long nLocalVertex, nGlobalVertex, MaxLocalVertex, *Buffer_Send_nVertex, *Buffer_Receive_nVertex, nBuffer;
        int nProcessor, iProcessor;
        double *Buffer_Send_Coord, *Buffer_Receive_Coord;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

        Xcoord_Airfoil.clear();
        Ycoord_Airfoil.clear();
        Zcoord_Airfoil.clear();
        Variable_Airfoil.clear();

        /*--- Set the right plane in 2D (note the change in Y-Z plane) ---*/
        if (nDim == 2)
        {
            Plane_P0[0] = 0.0;      Plane_P0[1] = 0.0;      Plane_P0[2] = 0.0;
            Plane_Normal[0] = 0.0;  Plane_Normal[1] = 1.0;  Plane_Normal[2] = 0.0;
        }

        /*--- the grid variation is stored using a vertices information,
          we should go from vertex to points ---*/
        if (original_surface == false)
        {
            Coord_Variation = new double *[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                Coord_Variation[iPoint] = new double[nDim];

            for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            {
                if (config->GetMarker_All_GeoEval(iMarker) == TBOX::YES)
                {
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        VarCoord = vertex[iMarker][iVertex]->GetVarCoord();
                        iPoint = vertex[iMarker][iVertex]->GetNode();
                        for (iDim = 0; iDim < nDim; iDim++)
                            Coord_Variation[iPoint][iDim] = VarCoord[iDim];
                    }
                }
            }
        }

        for (iMarker = 0; iMarker < nMarker; iMarker++)
        {
            if (config->GetMarker_All_GeoEval(iMarker) == TBOX::YES)
            {
                for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
                {
                    for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++)
                    {
                        iPoint = bound[iMarker][iElem]->GetNode(iNode);
                        for (jNode = 0; jNode < bound[iMarker][iElem]->GetnNodes(); jNode++)
                        {
                            jPoint = bound[iMarker][iElem]->GetNode(jNode);

                            if ((jPoint > iPoint) && ((node[iPoint]->GetCoord(0) > MinXCoord) && (node[iPoint]->GetCoord(0) < MaxXCoord)))
                            {
                                Segment_P0[0] = 0.0;  Segment_P0[1] = 0.0;  Segment_P0[2] = 0.0;  Variable_P0 = 0.0;
                                Segment_P1[0] = 0.0;  Segment_P1[1] = 0.0;  Segment_P1[2] = 0.0;  Variable_P1 = 0.0;

                                for (iDim = 0; iDim < nDim; iDim++)
                                {
                                    if (original_surface == true)
                                    {
                                        Segment_P0[iDim] = node[iPoint]->GetCoord(iDim);
                                        Segment_P1[iDim] = node[jPoint]->GetCoord(iDim);
                                    }
                                    else
                                    {
                                        Segment_P0[iDim] = node[iPoint]->GetCoord(iDim) + Coord_Variation[iPoint][iDim];
                                        Segment_P1[iDim] = node[jPoint]->GetCoord(iDim) + Coord_Variation[jPoint][iDim];
                                    }
                                }

                                if (FlowVariable != NULL)
                                {
                                    Variable_P0 = FlowVariable[iPoint];
                                    Variable_P1 = FlowVariable[jPoint];
                                }

                                /*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/
                                if (nDim == 2)
                                {
                                    Xcoord.push_back(Segment_P0[0]);    Xcoord.push_back(Segment_P1[0]);
                                    Ycoord.push_back(Segment_P0[2]);    Ycoord.push_back(Segment_P1[2]);
                                    Zcoord.push_back(Segment_P0[1]);    Zcoord.push_back(Segment_P1[1]);
                                    Variable.push_back(Variable_P0);    Variable.push_back(Variable_P1);
                                }
                                /*--- In 3D compute the intersection ---*/
                                else if (nDim == 3)
                                {
                                    intersect = SegmentIntersectsPlane(Segment_P0, Segment_P1, Variable_P0, Variable_P1, Plane_P0, Plane_Normal, Intersection, Variable_Interp);
                                    if (intersect == true)
                                    {
                                        Xcoord.push_back(Intersection[0]);
                                        Ycoord.push_back(Intersection[1]);
                                        Zcoord.push_back(Intersection[2]);
                                        Variable.push_back(Variable_Interp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (original_surface == false)
        {
            for (iPoint = 0; iPoint < nPoint; iPoint++)
                delete[] Coord_Variation[iPoint];
            delete[] Coord_Variation;
        }

#ifdef HAVE_MPI
        /*--- Copy the coordinates of all the points in the plane to the master node ---*/
        nLocalVertex = 0, nGlobalVertex = 0, MaxLocalVertex = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);

        Buffer_Send_nVertex = new unsigned long[1];
        Buffer_Receive_nVertex = new unsigned long[nProcessor];

        nLocalVertex = Xcoord.size();

        Buffer_Send_nVertex[0] = nLocalVertex;

        MPI_Allreduce(&nLocalVertex, &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

        Buffer_Send_Coord = new double[MaxLocalVertex * 4];
        Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex * 4];
        nBuffer = MaxLocalVertex * 4;

        for (iVertex = 0; iVertex < nLocalVertex; iVertex++)
        {
            Buffer_Send_Coord[iVertex * 4 + 0] = Xcoord[iVertex];
            Buffer_Send_Coord[iVertex * 4 + 1] = Ycoord[iVertex];
            Buffer_Send_Coord[iVertex * 4 + 2] = Zcoord[iVertex];
            Buffer_Send_Coord[iVertex * 4 + 3] = Variable[iVertex];
        }

        MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);

        /*--- Clean the vectors before adding the new vertices only to the master node ---*/

        Xcoord.clear();
        Ycoord.clear();
        Zcoord.clear();
        Variable.clear();

        /*--- Copy the boundary to the master node vectors ---*/
        if (rank == MASTER_NODE)
        {
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
            {
                for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++)
                {
                    Xcoord.push_back(Buffer_Receive_Coord[iProcessor*MaxLocalVertex * 4 + iVertex * 4 + 0]);
                    Ycoord.push_back(Buffer_Receive_Coord[iProcessor*MaxLocalVertex * 4 + iVertex * 4 + 1]);
                    Zcoord.push_back(Buffer_Receive_Coord[iProcessor*MaxLocalVertex * 4 + iVertex * 4 + 2]);
                    Variable.push_back(Buffer_Receive_Coord[iProcessor*MaxLocalVertex * 4 + iVertex * 4 + 3]);
                }
            }
        }

        delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
        delete[] Buffer_Send_nVertex; delete[] Buffer_Receive_nVertex;

#endif

        if ((rank == TBOX::MASTER_NODE) && (Xcoord.size() != 0))
        {
            /*--- Create a list with the duplicated points ---*/
            for (iVertex = 0; iVertex < Xcoord.size() - 1; iVertex++)
            {
                for (jVertex = iVertex + 1; jVertex < Xcoord.size(); jVertex++)
                {
                    Segment[0] = Xcoord[jVertex] - Xcoord[iVertex];
                    Segment[1] = Ycoord[jVertex] - Ycoord[iVertex];
                    Segment[2] = Zcoord[jVertex] - Zcoord[iVertex];
                    Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
                    if (Dist_Value < 1E-6)
                    {
                        Duplicate.push_back(jVertex);
                    }
                }
            }

            sort(Duplicate.begin(), Duplicate.end());
            it = unique(Duplicate.begin(), Duplicate.end());
            Duplicate.resize(it - Duplicate.begin());

            /*--- Remove duplicated points (starting from the back) ---*/
            for (iVertex = Duplicate.size(); iVertex > 0; iVertex--)
            {
                Xcoord.erase(Xcoord.begin() + Duplicate[iVertex - 1]);
                Ycoord.erase(Ycoord.begin() + Duplicate[iVertex - 1]);
                Zcoord.erase(Zcoord.begin() + Duplicate[iVertex - 1]);
                Variable.erase(Variable.begin() + Duplicate[iVertex - 1]);
            }

            if (Xcoord.size() != 1)
            {
                /*--- Find the trailing edge ---*/
                Trailing_Point = 0; Trailing_Coord = Xcoord[0];
                for (iVertex = 1; iVertex < Xcoord.size(); iVertex++)
                {
                    if (Xcoord[iVertex] > Trailing_Coord)
                    {
                        Trailing_Point = iVertex; Trailing_Coord = Xcoord[iVertex];
                    }
                }

                /*--- Add the trailing edge to the list, and remove from the original list ---*/
                Xcoord_Airfoil.push_back(Xcoord[Trailing_Point]);
                Ycoord_Airfoil.push_back(Ycoord[Trailing_Point]);
                Zcoord_Airfoil.push_back(Zcoord[Trailing_Point]);
                Variable_Airfoil.push_back(Variable[Trailing_Point]);

                Xcoord.erase(Xcoord.begin() + Trailing_Point);
                Ycoord.erase(Ycoord.begin() + Trailing_Point);
                Zcoord.erase(Zcoord.begin() + Trailing_Point);
                Variable.erase(Variable.begin() + Trailing_Point);

                /*--- Find the next point using the right hand side rule ---*/
                MinDist_Value = 1E6; MinDist_Point = 0;
                for (iVertex = 0; iVertex < Xcoord.size(); iVertex++)
                {
                    Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[0];
                    Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[0];
                    Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[0];
                    Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
                    Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;

                    if ((Dist_Value < MinDist_Value) && (Segment[2] > 0.0))
                    {
                        MinDist_Point = iVertex; MinDist_Value = Dist_Value;
                    }
                }

                Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);
                Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);
                Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);
                Variable_Airfoil.push_back(Variable[MinDist_Point]);

                Xcoord.erase(Xcoord.begin() + MinDist_Point);
                Ycoord.erase(Ycoord.begin() + MinDist_Point);
                Zcoord.erase(Zcoord.begin() + MinDist_Point);
                Variable.erase(Variable.begin() + MinDist_Point);

                /*--- Algorithm for the rest of the points ---*/
                do
                {
                    /*--- Last added point in the list ---*/
                    Airfoil_Point = Xcoord_Airfoil.size() - 1;

                    /*--- Compute the slope of the curve ---*/
                    Airfoil_Tangent[0] = Xcoord_Airfoil[Airfoil_Point] - Xcoord_Airfoil[Airfoil_Point - 1];
                    Airfoil_Tangent[1] = Ycoord_Airfoil[Airfoil_Point] - Ycoord_Airfoil[Airfoil_Point - 1];
                    Airfoil_Tangent[2] = Zcoord_Airfoil[Airfoil_Point] - Zcoord_Airfoil[Airfoil_Point - 1];
                    Length = sqrt(pow(Airfoil_Tangent[0], 2.0) + pow(Airfoil_Tangent[1], 2.0) + pow(Airfoil_Tangent[2], 2.0));
                    Airfoil_Tangent[0] /= Length; Airfoil_Tangent[1] /= Length; Airfoil_Tangent[2] /= Length;

                    /*--- Find the closest point with the right slope ---*/
                    MinDist_Value = 1E6; MinDistAngle_Value = 180;
                    MinDist_Point = -1; MinDistAngle_Point = -1;
                    for (iVertex = 0; iVertex < Xcoord.size(); iVertex++)
                    {
                        Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[Airfoil_Point];
                        Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[Airfoil_Point];
                        Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[Airfoil_Point];

                        /*--- Compute the distance to each point ---*/
                        Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));

                        /*--- Compute the angle of the point ---*/
                        Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;

                        /*--- Clip the value of the cosine, this is important due to the round errors ---*/
                        CosValue = Airfoil_Tangent[0] * Segment[0] + Airfoil_Tangent[1] * Segment[1] + Airfoil_Tangent[2] * Segment[2];
                        if (CosValue >= 1.0) CosValue = 1.0;
                        if (CosValue <= -1.0) CosValue = -1.0;

                        Angle_Value = acos(CosValue) * 180 / TBOX::PI_NUMBER;

                        if (Dist_Value < MinDist_Value)
                        {
                            MinDist_Point = iVertex; MinDist_Value = Dist_Value;
                        }

                        if ((Dist_Value < MinDistAngle_Value) && (Angle_Value < MaxAngle))
                        {
                            MinDistAngle_Point = iVertex;
                            MinDistAngle_Value = Dist_Value;
                        }
                    }

                    if (MinDistAngle_Point != -1) MinDist_Point = MinDistAngle_Point;

                    /*--- Add and remove the min distance to the list ---*/
                    Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);
                    Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);
                    Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);
                    Variable_Airfoil.push_back(Variable[MinDist_Point]);

                    Xcoord.erase(Xcoord.begin() + MinDist_Point);
                    Ycoord.erase(Ycoord.begin() + MinDist_Point);
                    Zcoord.erase(Zcoord.begin() + MinDist_Point);
                    Variable.erase(Variable.begin() + MinDist_Point);
                } while (Xcoord.size() != 0);

                /*--- Clean the vector before using them again for storing the upper or the lower side ---*/
                Xcoord.clear();
                Ycoord.clear();
                Zcoord.clear();
                Variable.clear();
            }
        }
    }
    
    double GeomUtilities::ComputeMaxThickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, AxisOriType val_axisOriType, unsigned short val_nDim,
                                              vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) 
    {
        unsigned long iVertex, jVertex, n, Trailing_Point, Leading_Point;
        double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, zp1, zpn, MaxThickness_Value = 0, Thickness, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, MaxDistance, Distance, AoA;
        vector<double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;

        /*--- Find the leading and trailing edges and compute the angle of attack ---*/
        MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

            if (MaxDistance < Distance) 
            { 
                MaxDistance = Distance; 
                Leading_Point = iVertex; 
            }
        }

        AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point])) * 180 / PI;

        /*--- Translate to the origin ---*/
        Xcoord_Trailing = Xcoord_Airfoil[0];
        Ycoord_Trailing = Ycoord_Airfoil[0];
        Zcoord_Trailing = Zcoord_Airfoil[0];

        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
            Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
            Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
        }

        /*--- Rotate the airfoil ---*/
        ValCos = cos(AoA*PI / 180.0);
        ValSin = sin(AoA*PI / 180.0);
        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            XValue = Xcoord_Airfoil_[iVertex];
            ZValue = Zcoord_Airfoil_[iVertex];

            Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
            Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
        }

        /*--- Identify upper and lower side, and store the value of the normal --*/
        for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) 
        {
            Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex - 1];
            Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex - 1];
            Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex - 1];
            Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
            Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;

            BiNormal[0] = Plane_Normal[0];
            BiNormal[1] = Plane_Normal[1];
            BiNormal[2] = Plane_Normal[2];
            Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
            BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;

            Normal[0] = Tangent[1] * BiNormal[2] - Tangent[2] * BiNormal[1];
            Normal[1] = Tangent[2] * BiNormal[0] - Tangent[0] * BiNormal[2];
            Normal[2] = Tangent[0] * BiNormal[1] - Tangent[1] * BiNormal[0];

            Xcoord_Normal.push_back(Normal[0]);
            Ycoord_Normal.push_back(Normal[1]);
            Zcoord_Normal.push_back(Normal[2]);

            unsigned short index = 2;
            if ((val_axisOriType == Z_AXIS) && (val_nDim == 3))
                index = 0;

            if (Normal[index] >= 0.0) 
            {
                Xcoord.push_back(Xcoord_Airfoil_[iVertex]);
                Ycoord.push_back(Ycoord_Airfoil_[iVertex]);
                Zcoord.push_back(Zcoord_Airfoil_[iVertex]);
            }
        }

        /*--- Order the arrays using the X component ---*/
        for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) 
        {
            for (jVertex = 0; jVertex < Xcoord.size() - 1 - iVertex; jVertex++) 
            {
                if (Xcoord[jVertex] > Xcoord[jVertex + 1]) 
                {
                    auxXCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex + 1]; Xcoord[jVertex + 1] = auxXCoord;
                    auxYCoord = Ycoord[jVertex]; Ycoord[jVertex] = Ycoord[jVertex + 1]; Ycoord[jVertex + 1] = auxYCoord;
                    auxZCoord = Zcoord[jVertex]; Zcoord[jVertex] = Zcoord[jVertex + 1]; Zcoord[jVertex + 1] = auxZCoord;
                }
            }
        }

        n = Xcoord.size();
        zp1 = (Zcoord[1] - Zcoord[0]) / (Xcoord[1] - Xcoord[0]);
        zpn = (Zcoord[n - 1] - Zcoord[n - 2]) / (Xcoord[n - 1] - Xcoord[n - 2]);
        Z2coord.resize(n + 1);
        SetSpline(Xcoord, Zcoord, n, zp1, zpn, Z2coord);

        /*--- Compute the thickness (we add a fabs because we can not guarantee the
          right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/
        MaxThickness_Value = 0.0;
        for (iVertex = 0; iVertex < Xcoord_Airfoil_.size(); iVertex++)
        {
            if (Zcoord_Normal[iVertex] < 0.0) 
            {
                Thickness = fabs(Zcoord_Airfoil_[iVertex] - GetSpline(Xcoord, Zcoord, Z2coord, n, Xcoord_Airfoil_[iVertex]));
                if (Thickness > MaxThickness_Value) 
                { 
                    MaxThickness_Value = Thickness; 
                }
            }
        }

        return MaxThickness_Value;
    }

    double GeomUtilities::ComputeThickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, double Location, AxisOriType val_axisOriType, unsigned short val_nDim,
                                           vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) 
    {
        unsigned long iVertex, jVertex, n_Upper, n_Lower, Trailing_Point, Leading_Point;
        double Thickness_Location, Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Thickness_Value = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, zp1, zpn, Chord, MaxDistance, Distance, AoA;
        vector<double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Z2coord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;

        /*--- Find the leading and trailing edges and compute the angle of attack ---*/
        MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

            if (MaxDistance < Distance) 
            { 
                MaxDistance = Distance; 
                Leading_Point = iVertex; 
            }
        }

        AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point])) * 180 / PI;
        Chord = MaxDistance;

        /*--- Translate to the origin ---*/
        Xcoord_Trailing = Xcoord_Airfoil[0];
        Ycoord_Trailing = Ycoord_Airfoil[0];
        Zcoord_Trailing = Zcoord_Airfoil[0];

        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
            Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
            Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
        }

        /*--- Rotate the airfoil ---*/
        ValCos = cos(AoA*PI / 180.0);
        ValSin = sin(AoA*PI / 180.0);

        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            XValue = Xcoord_Airfoil_[iVertex];
            ZValue = Zcoord_Airfoil_[iVertex];

            Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
            Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
        }

        /*--- Identify upper and lower side, and store the value of the normal --*/
        for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) 
        {
            Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex - 1];
            Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex - 1];
            Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex - 1];
            Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
            Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;

            BiNormal[0] = Plane_Normal[0];
            BiNormal[1] = Plane_Normal[1];
            BiNormal[2] = Plane_Normal[2];
            Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
            BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;

            Normal[0] = Tangent[1] * BiNormal[2] - Tangent[2] * BiNormal[1];
            Normal[1] = Tangent[2] * BiNormal[0] - Tangent[0] * BiNormal[2];
            Normal[2] = Tangent[0] * BiNormal[1] - Tangent[1] * BiNormal[0];

            Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);

            unsigned short index = 2;
            if ((val_axisOriType == Z_AXIS) && (val_nDim == 3)) index = 0;

            if (Normal[index] >= 0.0) 
            {
                Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
                Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
                Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
            }
            else 
            {
                Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
                Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
                Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
            }
        }

        /*--- Order the arrays using the X component ---*/
        for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) 
        {
            for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) 
            {
                if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex + 1]) 
                {
                    auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex + 1]; Xcoord_Upper[jVertex + 1] = auxXCoord;
                    auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex + 1]; Ycoord_Upper[jVertex + 1] = auxYCoord;
                    auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex + 1]; Zcoord_Upper[jVertex + 1] = auxZCoord;
                }
            }
        }

        /*--- Order the arrays using the X component ---*/
        for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) 
        {
            for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) 
            {
                if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex + 1]) 
                {
                    auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex + 1]; Xcoord_Lower[jVertex + 1] = auxXCoord;
                    auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex + 1]; Ycoord_Lower[jVertex + 1] = auxYCoord;
                    auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex + 1]; Zcoord_Lower[jVertex + 1] = auxZCoord;
                }
            }
        }

        n_Upper = Xcoord_Upper.size();
        zp1 = (Zcoord_Upper[1] - Zcoord_Upper[0]) / (Xcoord_Upper[1] - Xcoord_Upper[0]);
        zpn = (Zcoord_Upper[n_Upper - 1] - Zcoord_Upper[n_Upper - 2]) / (Xcoord_Upper[n_Upper - 1] - Xcoord_Upper[n_Upper - 2]);
        Z2coord_Upper.resize(n_Upper + 1);
        SetSpline(Xcoord_Upper, Zcoord_Upper, n_Upper, zp1, zpn, Z2coord_Upper);

        n_Lower = Xcoord_Lower.size();
        zp1 = (Zcoord_Lower[1] - Zcoord_Lower[0]) / (Xcoord_Lower[1] - Xcoord_Lower[0]);
        zpn = (Zcoord_Lower[n_Lower - 1] - Zcoord_Lower[n_Lower - 2]) / (Xcoord_Lower[n_Lower - 1] - Xcoord_Lower[n_Lower - 2]);
        Z2coord_Lower.resize(n_Lower + 1);
        SetSpline(Xcoord_Lower, Zcoord_Lower, n_Lower, zp1, zpn, Z2coord_Lower);

        /*--- Compute the thickness (we add a fabs because we can not guarantee the
          right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/
        Thickness_Location = -Chord*(1.0 - Location);
        Thickness_Value = fabs(GetSpline(Xcoord_Upper, Zcoord_Upper, Z2coord_Upper, n_Upper, Thickness_Location) - GetSpline(Xcoord_Lower, Zcoord_Lower, Z2coord_Lower, n_Lower, Thickness_Location));

        return Thickness_Value;
    }
    
    double GeomUtilities::ComputeAoA(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) 
    {
        unsigned long iVertex, Trailing_Point, Leading_Point;
        double MaxDistance, Distance, AoA = 0.0;

        /*--- Find the leading and trailing edges and compute the angle of attack ---*/
        MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

            if (MaxDistance < Distance) 
            { 
                MaxDistance = Distance; 
                Leading_Point = iVertex; 
            }
        }

        AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point])) * 180 / PI;
        return AoA;
    }

    double GeomUtilities::ComputeChord(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) 
    {
        unsigned long iVertex, Trailing_Point;
        double MaxDistance, Distance, Chord = 0.0;

        /*--- Find the leading and trailing edges and compute the angle of attack ---*/
        MaxDistance = 0.0; Trailing_Point = 0;
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

            if (MaxDistance < Distance) 
            { 
                MaxDistance = Distance; 
            }
        }

        Chord = MaxDistance;
        return Chord;
    }

    double GeomUtilities::ComputeArea(double *Plane_P0, double *Plane_Normal, unsigned short iSection, AxisOriType val_axisOriType, unsigned short val_nDim,
                                      std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil, std::vector<double> &Zcoord_Airfoil, bool original_surface) 
    {
        unsigned long iVertex, jVertex;
        double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Area_Value = 0.0, Area_Value_Upper = 0.0, Area_Value_Lower = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue;
        std::vector<double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
        unsigned long Trailing_Point, Leading_Point;
        double MaxDistance, Distance, AoA;

        /*--- Find the leading and trailing edges and compute the angle of attack ---*/
        MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                            pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

            if (MaxDistance < Distance) 
            { 
                MaxDistance = Distance;
                Leading_Point = iVertex; 
            }
        }
        
        AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point])) * 180 / TBOX::PI_NUMBER;

        /*--- Translate to the origin ---*/
        Xcoord_Trailing = Xcoord_Airfoil[0];
        Ycoord_Trailing = Ycoord_Airfoil[0];
        Zcoord_Trailing = Zcoord_Airfoil[0];

        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
            Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
            Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
        }

        /*--- Rotate the airfoil ---*/
        ValCos = cos(AoA*TBOX::PI_NUMBER / 180.0);
        ValSin = sin(AoA*TBOX::PI_NUMBER / 180.0);
        
        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) 
        {
            XValue = Xcoord_Airfoil_[iVertex];
            ZValue = Zcoord_Airfoil_[iVertex];

            Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
            Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
        }

        /*--- Identify upper and lower side, and store the value of the normal --*/
        for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) 
        {
            Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex - 1];
            Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex - 1];
            Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex - 1];
            Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
            Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;

            BiNormal[0] = Plane_Normal[0];
            BiNormal[1] = Plane_Normal[1];
            BiNormal[2] = Plane_Normal[2];
            Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
            BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;

            Normal[0] = Tangent[1] * BiNormal[2] - Tangent[2] * BiNormal[1];
            Normal[1] = Tangent[2] * BiNormal[0] - Tangent[0] * BiNormal[2];
            Normal[2] = Tangent[0] * BiNormal[1] - Tangent[1] * BiNormal[0];

            Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);

            unsigned short index = 2;
            if ((val_axisOriType == Z_AXIS) && (val_nDim == 3)) index = 0;

            if (Normal[index] >= 0.0) 
            {
                Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
                Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
                Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
            }
            else 
            {
                Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
                Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
                Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
            }

        }

        /*--- Order the arrays using the X component ---*/
        for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) 
        {
            for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) 
            {
                if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex + 1]) 
                {
                    auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex + 1]; Xcoord_Upper[jVertex + 1] = auxXCoord;
                    auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex + 1]; Ycoord_Upper[jVertex + 1] = auxYCoord;
                    auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex + 1]; Zcoord_Upper[jVertex + 1] = auxZCoord;
                }
            }
        }

        /*--- Order the arrays using the X component ---*/
        for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) 
        {
            for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) 
            {
                if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex + 1]) 
                {
                    auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex + 1]; Xcoord_Lower[jVertex + 1] = auxXCoord;
                    auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex + 1]; Ycoord_Lower[jVertex + 1] = auxYCoord;
                    auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex + 1]; Zcoord_Lower[jVertex + 1] = auxZCoord;
                }
            }
        }

        /*--- Compute total area ---*/
        Area_Value_Upper = 0.0;
        Area_Value_Lower = 0.0;

        for (iVertex = 0; iVertex < Xcoord_Upper.size() - 1; iVertex++)
            Area_Value_Upper += (Xcoord_Upper[iVertex + 1] - Xcoord_Upper[iVertex]) * 0.5*(Zcoord_Upper[iVertex + 1] + Zcoord_Upper[iVertex]);
        for (iVertex = 0; iVertex < Xcoord_Lower.size() - 1; iVertex++)
            Area_Value_Lower += (Xcoord_Lower[iVertex + 1] - Xcoord_Lower[iVertex]) * 0.5*(Zcoord_Lower[iVertex + 1] + Zcoord_Lower[iVertex]);

        Area_Value = fabs(Area_Value_Upper - Area_Value_Lower);
        return Area_Value;
    }

    double GeomUtilities::Compute_Volume(TBOX::TBOX_Config *config, bool original_surface)
    {
        int rank = TBOX::MASTER_NODE;

        /*--- MPI initialization ---*/
#ifdef ARIES_HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

        unsigned short iPlane, nPlane = 0;
        double Volume = 0.0, MinPlane, MaxPlane, MinXCoord, MaxXCoord, dPlane, *Area;
        std::vector<double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;

        /*--- Make a large number of section cuts for approximating volume ---*/
        nPlane = config->GetnVolSections();

        /*--- Allocate memory for the section cutting ---*/
        Area = new double[nPlane];

        double **Plane_P0 = new double*[nPlane];
        for (iPlane = 0; iPlane < nPlane; iPlane++)
            Plane_P0[iPlane] = new double[nDim];

        double **Plane_Normal = new double*[nPlane];
        for (iPlane = 0; iPlane < nPlane; iPlane++)
            Plane_Normal[iPlane] = new double[nDim];
        
        MinPlane = config->GetSection_Location(0); MaxPlane = config->GetSection_Location(1);
        MinXCoord = -1E6; MaxXCoord = 1E6;
        dPlane = fabs((MaxPlane - MinPlane) / double(nPlane - 1));
        for (iPlane = 0; iPlane < nPlane; iPlane++)
        {
            Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
            Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
            Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;
            Plane_Normal[iPlane][config->GetAxis_Orientation()] = 1.0;
            Plane_P0[iPlane][config->GetAxis_Orientation()] = MinPlane + iPlane*dPlane;
        }

        /*--- Allocate some std::vectors for storing airfoil coordinates ---*/
        Xcoord_Airfoil   = new std::vector<double>[nPlane];
        Ycoord_Airfoil   = new std::vector<double>[nPlane];
        Zcoord_Airfoil   = new std::vector<double>[nPlane];
        Variable_Airfoil = new std::vector<double>[nPlane];

        /*--- Create the section slices through the geometry ---*/

        for (iPlane = 0; iPlane < nPlane; iPlane++)
        {
            ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], MinXCoord, MaxXCoord, NULL, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], Variable_Airfoil[iPlane], original_surface, config);
        }

        /*--- Compute the area at each section ---*/
        if (rank == TBOX::MASTER_NODE)
        {
            for (iPlane = 0; iPlane < nPlane; iPlane++)
            {
                Area[iPlane] = 0.0;
                if (Xcoord_Airfoil[iPlane].size() != 0)
                {
                    Area[iPlane] = Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane, config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane], original_surface);
                }
            }

            /*--- Compute the volume using a composite Simpson's rule ---*/
            Volume = 0.0;
            for (iPlane = 0; iPlane < nPlane - 2; iPlane += 2)
            {
                if (Xcoord_Airfoil[iPlane].size() != 0)
                {
                    Volume += (1.0 / 3.0)*dPlane*(Area[iPlane] + 4.0*Area[iPlane + 1] + Area[iPlane + 2]);
                }
            }
        }

        /*--- Free memory for the section cuts ---*/
        delete[] Xcoord_Airfoil;
        delete[] Ycoord_Airfoil;
        delete[] Zcoord_Airfoil;
        delete[] Variable_Airfoil;

        for (iPlane = 0; iPlane < nPlane; iPlane++)
            delete[] Plane_P0[iPlane];
        delete[] Plane_P0;

        for (iPlane = 0; iPlane < nPlane; iPlane++)
            delete Plane_Normal[iPlane];
        delete[] Plane_Normal;

        delete[] Area;

        /*--- Return the volume and exit ---*/
        return Volume;
    }

    
    void GEOM_Geometry::ComputeSurf_Curvature(TBOX::TBOX_Config *config)
    {
        unsigned short iMarker, iNeigh_Point, iDim, iNode, iNeighbor_Nodes, Neighbor_Node;
        unsigned long Neighbor_Point, iVertex, iPoint, jPoint, iElem_Bound, iEdge, nLocalVertex, MaxLocalVertex, *Buffer_Send_nVertex, *Buffer_Receive_nVertex, TotalnPointDomain;
        int iProcessor, nProcessor;
        std::vector<unsigned long> Point_NeighborList, Elem_NeighborList, Point_Triangle, Point_Edge, Point_Critical;
        std::vector<unsigned long>::iterator it;
        double U[3] = { 0.0, 0.0, 0.0 }, V[3] = { 0.0, 0.0, 0.0 }, W[3] = { 0.0, 0.0, 0.0 }, Length_U, Length_V, Length_W, CosValue, Angle_Value, *K, *Angle_Defect, *Area_Vertex, *Angle_Alpha, *Angle_Beta, **NormalMeanK, MeanK, GaussK, MaxPrinK, cot_alpha, cot_beta, delta, X1, X2, X3, Y1, Y2, Y3, radius, *Buffer_Send_Coord, *Buffer_Receive_Coord, *Coord, Dist, MinDist, MaxK, MinK, SigmaK;
        bool *Check_Edge;
        int rank;

#ifndef HAVE_MPI
        rank = TBOX::MASTER_NODE;
#else
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

        /*--- Allocate surface curvature ---*/
        K = new double[nPoint];

        if (nDim == 2)
        {
            /*--- Loop over all the markers ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
                {
                    /*--- Loop through all marker vertices again, this time also
                      finding the neighbors of each node.---*/
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();

                        if (node[iPoint]->GetDomain())
                        {
                            /*--- Loop through neighbors. In 2-D, there should be 2 nodes on either
                              side of this vertex that lie on the same surface. ---*/
                            Point_Edge.clear();

                            for (iNeigh_Point = 0; iNeigh_Point < node[iPoint]->GetnPoint(); iNeigh_Point++)
                            {
                                Neighbor_Point = node[iPoint]->GetPoint(iNeigh_Point);

                                /*--- Check if this neighbor lies on the surface. If so,
                                  add to the list of neighbors. ---*/
                                if (node[Neighbor_Point]->GetPhysicalBoundary())
                                {
                                    Point_Edge.push_back(Neighbor_Point);
                                }

                            }

                            if (Point_Edge.size() == 2)
                            {
                                /*--- Compute the curvature using three points ---*/
                                X1 = node[iPoint]->GetCoord(0);
                                X2 = node[Point_Edge[0]]->GetCoord(0);
                                X3 = node[Point_Edge[1]]->GetCoord(0);
                                Y1 = node[iPoint]->GetCoord(1);
                                Y2 = node[Point_Edge[0]]->GetCoord(1);
                                Y3 = node[Point_Edge[1]]->GetCoord(1);

                                radius = sqrt(((X2 - X1)*(X2 - X1) + (Y2 - Y1)*(Y2 - Y1))*
                                              ((X2 - X3)*(X2 - X3) + (Y2 - Y3)*(Y2 - Y3))*
                                              ((X3 - X1)*(X3 - X1) + (Y3 - Y1)*(Y3 - Y1))) /
                                        (2.0*fabs(X1*Y2 + X2*Y3 + X3*Y1 - X1*Y3 - X2*Y1 - X3*Y2));

                                K[iPoint] = 1.0 / radius;
                                node[iPoint]->SetCurvature(K[iPoint]);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            Angle_Defect = new double[nPoint];
            Area_Vertex = new double[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                Angle_Defect[iPoint] = 2 * TBOX::PI_NUMBER;
                Area_Vertex[iPoint] = 0.0;
            }

            Angle_Alpha = new double[nEdge];
            Angle_Beta = new double[nEdge];
            Check_Edge = new bool[nEdge];
            for (iEdge = 0; iEdge < nEdge; iEdge++)
            {
                Angle_Alpha[iEdge] = 0.0;
                Angle_Beta[iEdge] = 0.0;
                Check_Edge[iEdge] = true;
            }

            NormalMeanK = new double *[nPoint];
            for (iPoint = 0; iPoint < nPoint; iPoint++)
            {
                NormalMeanK[iPoint] = new double[nDim];
                for (iDim = 0; iDim < nDim; iDim++)
                {
                    NormalMeanK[iPoint][iDim] = 0.0;
                }
            }

            /*--- Loop over all the markers ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
                {
                    /*--- Loop over all the boundary elements ---*/
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        /*--- Only triangles ---*/
                        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TBOX::TRIANGLE)
                        {
                            /*--- Loop over all the nodes of the boundary element ---*/
                            for (iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++)
                            {
                                iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);

                                Point_Triangle.clear();

                                for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++)
                                {
                                    Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                                    Neighbor_Point = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                                    Point_Triangle.push_back(Neighbor_Point);
                                }

                                iEdge = FindEdge(Point_Triangle[0], Point_Triangle[1]);

                                for (iDim = 0; iDim < nDim; iDim++) {
                                    U[iDim] = node[Point_Triangle[0]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                                    V[iDim] = node[Point_Triangle[1]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                                }

                                W[0] = 0.5*(U[1] * V[2] - U[2] * V[1]); W[1] = -0.5*(U[0] * V[2] - U[2] * V[0]); W[2] = 0.5*(U[0] * V[1] - U[1] * V[0]);

                                Length_U = 0.0, Length_V = 0.0, Length_W = 0.0, CosValue = 0.0;
                                for (iDim = 0; iDim < nDim; iDim++) { Length_U += U[iDim] * U[iDim]; Length_V += V[iDim] * V[iDim]; Length_W += W[iDim] * W[iDim]; }
                                Length_U = sqrt(Length_U); Length_V = sqrt(Length_V); Length_W = sqrt(Length_W);
                                for (iDim = 0; iDim < nDim; iDim++) { U[iDim] /= Length_U; V[iDim] /= Length_V; CosValue += U[iDim] * V[iDim]; }
                                if (CosValue >= 1.0) CosValue = 1.0;
                                if (CosValue <= -1.0) CosValue = -1.0;

                                Angle_Value = acos(CosValue);
                                Area_Vertex[iPoint] += Length_W;
                                Angle_Defect[iPoint] -= Angle_Value;
                                if (Angle_Alpha[iEdge] == 0.0) Angle_Alpha[iEdge] = Angle_Value;
                                else Angle_Beta[iEdge] = Angle_Value;
                            }
                        }
                    }
                }
            }

            /*--- Compute mean curvature ---*/
            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
                {
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
                    {
                        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TBOX::TRIANGLE)
                        {
                            for (iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++)
                            {
                                iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);

                                for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++)
                                {
                                    Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                                    jPoint = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);

                                    iEdge = FindEdge(iPoint, jPoint);

                                    if (Check_Edge[iEdge])
                                    {

                                        Check_Edge[iEdge] = false;

                                        if (tan(Angle_Alpha[iEdge]) != 0.0) cot_alpha = 1.0 / tan(Angle_Alpha[iEdge]); else cot_alpha = 0.0;
                                        if (tan(Angle_Beta[iEdge]) != 0.0) cot_beta = 1.0 / tan(Angle_Beta[iEdge]); else cot_beta = 0.0;

                                        /*--- iPoint, and jPoint ---*/
                                        for (iDim = 0; iDim < nDim; iDim++)
                                        {
                                            if (Area_Vertex[iPoint] != 0.0) NormalMeanK[iPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[iPoint]->GetCoord(iDim) - node[jPoint]->GetCoord(iDim)) / Area_Vertex[iPoint];
                                            if (Area_Vertex[jPoint] != 0.0) NormalMeanK[jPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[jPoint]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim)) / Area_Vertex[jPoint];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /*--- Compute Gauss, mean, max and min principal curvature,
              and set the list of critical points ---*/

            for (iMarker = 0; iMarker < nMarker; iMarker++)
            {
                if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
                {
                    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                    {
                        iPoint = vertex[iMarker][iVertex]->GetNode();

                        if (node[iPoint]->GetDomain())
                        {

                            if (Area_Vertex[iPoint] != 0.0) GaussK = 3.0*Angle_Defect[iPoint] / Area_Vertex[iPoint];
                            else GaussK = 0.0;

                            MeanK = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++)
                                MeanK += NormalMeanK[iPoint][iDim] * NormalMeanK[iPoint][iDim];
                            MeanK = sqrt(MeanK);

                            delta = std::max((MeanK*MeanK - GaussK), 0.0);

                            MaxPrinK = MeanK + sqrt(delta);

                            /*--- Store the curvature value ---*/
                            K[iPoint] = MaxPrinK;
                            node[iPoint]->SetCurvature(K[iPoint]);
                        }
                    }
                }
            }

            delete[] Angle_Defect;
            delete[] Area_Vertex;
            delete[] Angle_Alpha;
            delete[] Angle_Beta;
            delete[] Check_Edge;

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                delete NormalMeanK[iPoint];
            delete[] NormalMeanK;
        }

        /*--- Sharp edge detection is based in the statistical
          distribution of the curvature ---*/

        MaxK = K[0]; MinK = K[0]; MeanK = 0.0; TotalnPointDomain = 0;
        for (iMarker = 0; iMarker < nMarker; iMarker++)
        {
            if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
            {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    if (node[iPoint]->GetDomain())
                    {
                        MaxK = std::max(MaxK, fabs(K[iPoint]));
                        MinK = std::min(MinK, fabs(K[iPoint]));
                        MeanK += fabs(K[iPoint]);
                        TotalnPointDomain++;
                    }
                }
            }
        }

#ifdef HAVE_MPI
        double MyMeanK = MeanK; MeanK = 0.0;
        double MyMaxK = MaxK; MaxK = 0.0;
        unsigned long MynPointDomain = TotalnPointDomain; TotalnPointDomain = 0;
        MPI_Allreduce(&MyMeanK, &MeanK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&MyMaxK, &MaxK, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&MynPointDomain, &TotalnPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

        /*--- Compute the mean ---*/
        MeanK /= double(TotalnPointDomain);

        /*--- Compute the standard deviation ---*/
        SigmaK = 0.0;
        for (iMarker = 0; iMarker < nMarker; iMarker++)
        {
            if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
            {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    if (node[iPoint]->GetDomain())
                    {
                        SigmaK += (fabs(K[iPoint]) - MeanK) * (fabs(K[iPoint]) - MeanK);
                    }
                }
            }
        }

#ifdef HAVE_MPI
        double MySigmaK = SigmaK; SigmaK = 0.0;
        MPI_Allreduce(&MySigmaK, &SigmaK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

        SigmaK = sqrt(SigmaK / double(TotalnPointDomain));

        if (rank == TBOX::MASTER_NODE)
            std::cout << "Max K: " << MaxK << ". Mean K: " << MeanK << ". Standard deviation K: " << SigmaK << "." << std::endl;

        Point_Critical.clear();

        for (iMarker = 0; iMarker < nMarker; iMarker++)
        {
            if (config->GetMarker_All_KindBC(iMarker) != TBOX::SEND_RECEIVE)
            {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
                {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    if (node[iPoint]->GetDomain())
                    {
                        if (fabs(K[iPoint]) > MeanK + config->GetRefSharpEdges()*SigmaK)
                        {
                            Point_Critical.push_back(iPoint);
                        }
                    }
                }
            }
        }

        /*--- Variables and buffers needed for MPI ---*/
#ifdef HAVE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#else
        nProcessor = 1;
#endif
        Buffer_Send_nVertex = new unsigned long[1];
        Buffer_Receive_nVertex = new unsigned long[nProcessor];

        /*--- Count the total number of critical edge nodes. ---*/
        nLocalVertex = Point_Critical.size();
        Buffer_Send_nVertex[0] = nLocalVertex;

        /*--- Communicate to all processors the total number of critical edge nodes. ---*/
#ifdef HAVE_MPI
        MaxLocalVertex = 0;
        MPI_Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
        MaxLocalVertex = nLocalVertex;
        Buffer_Receive_nVertex[0] = nLocalVertex;
#endif


        /*--- Create and initialize to zero some buffers to hold the coordinates
          of the boundary nodes that are communicated from each partition (all-to-all). ---*/

        Buffer_Send_Coord = new double[MaxLocalVertex*nDim];
        Buffer_Receive_Coord = new double[nProcessor*MaxLocalVertex*nDim];

#ifdef HAVE_MPI
        unsigned long nBuffer = MaxLocalVertex*nDim;
#endif

        for (iVertex = 0; iVertex < MaxLocalVertex; iVertex++)
        {
            for (iDim = 0; iDim < nDim; iDim++)
            {
                Buffer_Send_Coord[iVertex*nDim + iDim] = 0.0;
            }
        }

        /*--- Retrieve and store the coordinates of the sharp edges boundary nodes on
          the local partition and broadcast them to all partitions. ---*/

        for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++)
        {
            iPoint = Point_Critical[iVertex];
            for (iDim = 0; iDim < nDim; iDim++)
                Buffer_Send_Coord[iVertex*nDim + iDim] = node[iPoint]->GetCoord(iDim);
        }

#ifdef HAVE_MPI
        MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
#else
        for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++)
        {
            for (iDim = 0; iDim < nDim; iDim++)
            {
                Buffer_Receive_Coord[iVertex*nDim + iDim] = Buffer_Send_Coord[iVertex*nDim + iDim];
            }
        }
#endif

        /*--- Loop over all interior mesh nodes on the local partition and compute
          the distances to each of the no-slip boundary nodes in the entire mesh.
          Store the minimum distance to the wall for each interior mesh node. ---*/

        for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
        {
            Coord = node[iPoint]->GetCoord();

            MinDist = 1E20;
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
            {
                for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++)
                {
                    Dist = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++)
                    {
                        Dist += (Coord[iDim] - Buffer_Receive_Coord[(iProcessor*MaxLocalVertex + iVertex)*nDim + iDim])*
                                (Coord[iDim] - Buffer_Receive_Coord[(iProcessor*MaxLocalVertex + iVertex)*nDim + iDim]);
                    }
                    Dist = sqrt(Dist);
                    if (Dist < MinDist) MinDist = Dist;
                }
            }
            node[iPoint]->SetSharpEdge_Distance(MinDist);
        }

        /*--- Deallocate Max curvature ---*/
        delete[] K;

        /*--- Deallocate the buffers needed for the MPI communication. ---*/
        delete[] Buffer_Send_Coord;
        delete[] Buffer_Receive_Coord;
        delete[] Buffer_Send_nVertex;
        delete[] Buffer_Receive_nVertex;
    }
}
