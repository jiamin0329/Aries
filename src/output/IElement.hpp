/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for element
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_IELEMENT_HPP
#define ARIES_IELEMENT_HPP


namespace ARIES
{
   
    class IElement
    {
    public:
        explicit IElement()
        {
        };



        virtual unsigned short GetVtkType() = 0;
        
        /*!
         * \brief A pure virtual member.
         * \param[in] val_node - Local index of a node.
         * \return Global index of the node.
         */
        virtual unsigned long GetNode(unsigned short iNode) = 0;

            
    protected:
        virtual ~IElement()
        {
        };
        
    }; // end IElement class
} // end ARIES namespace

#endif
