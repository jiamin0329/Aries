/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for primal grid definition, base class of fvm grid classes
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    02-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_GRID_HPP
#define ARIES_GRID_HPP

#include <vector>

using namespace std;

namespace ARIES
{
    class Grid
    {
    public:
        Grid();
        Grid(unsigned short val_nNodes, unsigned short val_nFaces, unsigned short val_VTKType);

        virtual ~Grid();
        
        unsigned short GetNumDim() { return d_numDim; };
        void SetNumDin(unsigned short val_numDim) { d_numDim = val_numDim; };
        
        long GetNeighborElement(unsigned short val_face) { return d_neighborElement[val_face]; };
        void SetNeighborElement(unsigned long val_elem, unsigned short val_face) { d_neighborElement[val_face] = val_elem; };
        void GetAllNeighborElement(); // used for debug

        double GetCG(unsigned short val_iDim) { return d_coordCG[val_iDim]; };
        void SetCG(double **val_coord);
        
        double GetVolume() { return d_volume; };
        void SetVolume(double val_volume) { d_volume = val_volume; };
       
        double GetFaceCG(unsigned short val_face, unsigned short val_dim) { return d_coordFaceElemCG[val_face][val_dim]; };

        bool GetDivide(){ return d_isDivide; };
        void SetDivide(bool val_isDivide) { d_isDivide = val_isDivide; };

        virtual unsigned long GetDomainElement() { return d_domainElement; };
        virtual void SetDomainElement(unsigned long val_domainelement) { d_domainElement = val_domainelement; };

        virtual void ChangeOrientation() = 0;
        virtual unsigned short GetVTKType() = 0;
        
        virtual unsigned short GetRotationType() = 0;
        virtual void SetRotationType(unsigned short val_rotationType) = 0;

        virtual unsigned short GetNumNode() = 0;
        virtual unsigned short GetNumFace() = 0;
        virtual unsigned short GetNumNodeFace(unsigned short val_face) = 0;
        virtual unsigned short GetMaxNodeFace() = 0;
        virtual unsigned short GetNumNeighborNode(unsigned short val_node) = 0;
        virtual unsigned short GetNumNeighborElement() = 0;

        virtual unsigned long GetNode(unsigned short val_node) = 0;
        virtual void SetNode(unsigned short val_node, unsigned long val_point) = 0;

        virtual unsigned short GetFace(unsigned short val_face, unsigned short val_index) = 0;
        virtual unsigned short GetNeighborNode(unsigned short val_node, unsigned short val_index) = 0;

    protected:
        unsigned short d_numDim;		                      /*!< \brief Dimension of the element (2D or 3D) useful for triangles, rectangles and edges. */
        vector<unsigned long> d_node;                         /*!< \brief Vector to store the global nodes of an element. */
        vector<long> d_neighborElement;                       /*!< \brief Vector to store the elements surronding an element. */
        vector<double> d_coordCG;                             /*!< \brief Coordinates of the center-of-gravity of the element. */
        vector<vector<double> > d_coordFaceElemCG;	          /*!< \brief Coordinates of the center-of-gravity of the face of the elements. */
        double d_volume;                                      /*!< \brief Volume of the grid. */
        unsigned long d_domainElement;	                      /*!< \brief Only for boundaries, in this variable the 3D elements which correspond with a boundary element is stored. */
        bool d_isDivide;                                      /*!< \brief Marker used to know if we are going to divide this element in the adaptation proccess. */
    }; // end Grid class
} // end ARIES namespace

#endif



