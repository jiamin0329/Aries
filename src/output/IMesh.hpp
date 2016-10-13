/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for mesh
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_IMESH_HPP
#define ARIES_IMESH_HPP

#include "IDualGrid.hpp"
#include "IElement.hpp"


namespace ARIES
{
    class IMesh
    {
    public:
        explicit IMesh()
        {
        };

        virtual unsigned short GetDim() = 0;



        virtual unsigned long GetNumPoints() = 0;
        virtual unsigned long GetNumElements() = 0;
        
        /*!
         * \brief Set the number of grid points in the domain.
         * \param[in] val_npoint - Number of grid points in the domain.
         */
        virtual unsigned long GenNumPointDomain() = 0;
            
        virtual unsigned short GetNumVerticesOnMarker(unsigned short iMarker) = 0;


        virtual IDualGrid* GetVertex(unsigned short iMarker, unsigned short iVertex) = 0; //dgvertex
          
        virtual IDualGrid* GetNode(unsigned long iPoint) = 0; // dgpoint


        virtual IElement* GetElement(unsigned long iElem) = 0;
        

        virtual unsigned long GetNumElemTria() = 0;
        virtual unsigned long GetNumElemQuad() = 0;
        virtual unsigned long GetNumElemTetr() = 0;
        virtual unsigned long GetNumElemHexa() = 0;
        virtual unsigned long GetNumElemPris() = 0;
        virtual unsigned long GetNumElemPyra() = 0;

        /*!
         * \brief Get the number of boundary elements.
         * \param[in] val_marker - Marker of the boundary.
         */
        virtual unsigned long GetNumElemBoundary(unsigned short iMarker) = 0;


        virtual IElement* GetElementBoundary (unsigned short iMarker, unsigned long iElem) = 0;


            
    protected:
        virtual ~IMesh()
        {
        };
        
    }; // end IMesh class
} // end ARIES namespace

#endif



