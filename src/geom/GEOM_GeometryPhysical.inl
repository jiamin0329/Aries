
#ifndef ARIES_GEOM_GEOMETRYPHYSICAL_INLINE
#define ARIES_GEOM_GEOMETRYPHYSICAL_INLINE

namespace ARIES
{
    namespace GEOM
    {
        inline long GEOM_GeometryPhysical::GetGlobal_to_Local_Point(long val_ipoint) { return Global_to_Local_Point[val_ipoint]; }

        inline unsigned short GEOM_GeometryPhysical::GetGlobal_to_Local_Marker(unsigned short val_imarker) { return Global_to_Local_Marker[val_imarker]; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nPoint(void) { return Global_nPoint; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nPointDomain(void) { return Global_nPointDomain; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElem(void) { return Global_nElem; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemLine(void) { return Global_nelem_edge; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemTria(void) { return Global_nelem_triangle; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemQuad(void) { return Global_nelem_quad; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemTetr(void) { return Global_nelem_tetra; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemHexa(void) { return Global_nelem_hexa; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemPris(void) { return Global_nelem_prism; }

        inline unsigned long GEOM_GeometryPhysical::GetGlobal_nElemPyra(void) { return Global_nelem_pyramid; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemLine(void) { return nelem_edge; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemTria(void) { return nelem_triangle; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemQuad(void) { return nelem_quad; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemTetr(void) { return nelem_tetra; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemHexa(void) { return nelem_hexa; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemPris(void) { return nelem_prism; }

        inline unsigned long GEOM_GeometryPhysical::GetnElemPyra(void) { return nelem_pyramid; }

        inline std::vector<double> GEOM_GeometryPhysical::GetGeometryPlanes() { return XCoordList; }

        inline std::vector<std::vector<double> > GEOM_GeometryPhysical::GetXCoord() { return Xcoord_plane; }

        inline std::vector<std::vector<double> > GEOM_GeometryPhysical::GetYCoord() { return Ycoord_plane; }

        inline std::vector<std::vector<double> > GEOM_GeometryPhysical::GetZCoord() { return Zcoord_plane; }

        inline std::vector<std::vector<unsigned long> > GEOM_GeometryPhysical::GetPlanarPoints() { return Plane_points; }

        inline void GEOM_GeometryPhysical::SetPoint_Connectivity(GEOM_Geometry *geometry) { GEOM_Geometry::SetPoint_Connectivity(geometry); }
    }
}

#endif