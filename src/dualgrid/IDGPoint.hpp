/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for dual grid
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_IDUALGRID_HPP
#define ARIES_IDUALGRID_HPP

namespace ARIES
{

    class IDualGrid
    {
    public:
        explicit IDualGrid()
        {
        };

            
        virtual double GetCoord(int dim) = 0;
        virtual unsigned long GetGlobalIndex() = 0;
        virtual unsigned long GetNode() = 0;

        virtual unsigned short GetRotationType () = 0;

        /*!
         * \brief For parallel computation, its indicates if a point must be computed or not.
         * \return <code>TRUE</code> if the node belong to the physical domain; otherwise <code>FALSE</code>.
         */
        virtual bool GetDomain() = 0; //point

    protected:
        virtual ~IDualGrid()
        {
        };
        
    }; // end IDualGrid class
} // end ARIES namespace

#endif








