/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for solution
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_ISOLUTION_HPP
#define ARIES_ISOLUTION_HPP

#include "IVariable.hpp"
#include "IVariableAdj.hpp"

namespace ARIES
{
    class ISolution
    {
    public:
        virtual ~ISolution( )
        {		
        };

        virtual IVariable* GetVariable(unsigned long iPoint) = 0;
        virtual IVariableAdj* GetVariableAdj(unsigned long iPoint) = 0;
            
        virtual double GetCoeffPressure(unsigned short iMarker, unsigned long iVertex) = 0;
        virtual double GetCoeffSkinFriction(unsigned short iMarker, unsigned long iVertex) = 0;
        virtual double GetCoeffSensitivity(unsigned short iMarker, unsigned long iVertex) = 0;
            
        virtual double GetHeatFlux(unsigned short iMarker, unsigned long iVertex) = 0;





            
    protected:
        explicit ISolution()
        {
        };
        
    }; // end ISolution class
} // end ARIES namespace

#endif

