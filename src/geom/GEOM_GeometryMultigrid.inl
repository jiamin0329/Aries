





#ifndef ARIES_GEOM_GEOMETRYMULTIGRID_INLINE
#define ARIES_GEOM_GEOMETRYMULTIGRID_INLINE

namespace ARIES
{
    namespace GEOM
    {
        inline void GEOM_GeometryMultigrid::SetPoint_Connectivity(void) { GEOM_Geometry::SetPoint_Connectivity(); }

        inline std::vector<double> GEOM_GeometryMultigrid::GetGeometryPlanes() { return XCoordList; }

        inline std::vector<std::vector<double> > GEOM_GeometryMultigrid::GetXCoord() { return Xcoord_plane; }

        inline std::vector<std::vector<double> > GEOM_GeometryMultigrid::GetYCoord() { return Ycoord_plane; }

        inline std::vector<std::vector<double> > GEOM_GeometryMultigrid::GetZCoord() { return Zcoord_plane; }

        inline std::vector<std::vector<unsigned long> > GEOM_GeometryMultigrid::GetPlanarPoints() { return Plane_points; }
    }
}

#endif