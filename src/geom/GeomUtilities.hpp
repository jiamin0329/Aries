/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for processor data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_GEOMUTILITIES_HPP
#define ARIES_GEOMUTILITIES_HPP

#include "IProcData.hpp"
#include <vector>

using namespace std;

namespace ARIES
{
    class GeomUtilities
    {
        /*!
         * \brief Get the distance between a plane (defined by three point) and a point.
         * \param[in] Coord - Coordinates of the point.
         * \param[in] iCoord - Coordinates of the first point that defines the plane.
         * \param[in] jCoord - Coordinates of the second point that defines the plane.
         * \param[in] kCoord - Coordinates of the third point that defines the plane.
         * \return Signed distance.
         */
        static double Point2Plane_Distance(double *coord, double *iCoord, double *jCoord, double *kCoord);
        
        /*!
         * \brief Compute the intersection between a segment and a plane.
         * \param[in] Segment_P0 - Definition of the particular problem.
         * \param[in] Segment_P1 - Definition of the particular problem.
         * \param[in] Plane_P0 - Definition of the particular problem.
         * \param[in] Plane_Normal - Definition of the particular problem.
         * \param[in] Intersection - Definition of the particular problem.
         * \returns If the intersection has has been successful.
         */
        static bool RayIntersectsTriangle(double orig[3], double dir[3], double vert0[3], double vert1[3], double vert2[3], double *intersect);
        static bool SegmentIntersectsTriangle(double point0[3], double point1[3], double vert0[3], double vert1[3], double vert2[3]);
        static bool SegmentIntersectsPlane(double *Segment_P0, double *Segment_P1, double Variable_P0, double Variable_P1,
                                           double *Plane_P0, double *Plane_Normal, double *Intersection, double &Variable_Interp);

        static double GetSpline(vector<double> &xa, vector<double> &ya, vector<double>&y2a, unsigned long n, double x);
        static void SetSpline(vector<double> &x,  vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2);
        
        /*
         *  compute parameters of an airfoil
         */
        static void ComputeAirfoilSection(double *Plane_P0, double *Plane_Normal,
                                          double MinXCoord, double MaxXCoord, double *FlowVariable,
                                          std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil,
                                          std::vector<double> &Zcoord_Airfoil, std::vector<double> &Variable_Airfoil,
                                          bool original_surface, TBOX::TBOX_Config *config);
            
        static double ComputeMaxThickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, AxisOriType val_axisOriType, unsigned short val_nDim,
                                          vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface);
        static double ComputeThickness   (double *Plane_P0, double *Plane_Normal, unsigned short iSection, double Location, AxisOriType val_axisOriType, unsigned short val_nDim,
                                          vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface);
        
        static double ComputeAoA  (double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface);
        static double ComputeChord(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface);

        //static double ComputeArea(double *Plane_P0, double *Plane_Normal, unsigned short iSection, IProcData* procData, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface);
        //static double ComputeVolume(IProcData* procData, bool original_surface);
        
        
    };
    
}

#endif

















