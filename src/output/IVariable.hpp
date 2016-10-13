/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for variable
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    24-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_IVARIABLE_HPP
#define ARIES_IVARIABLE_HPP


namespace ARIES
{
    class IVariable
    {
    public:
        explicit IVariable()
        {
        };

        virtual double GetPressure() = 0;
        virtual double GetVelocity2() = 0;
        virtual double GetSoundSpeed() = 0;

    protected:
        virtual ~IVariable()
        {
        };
        
    }; // end IVariable class
} // end ARIES namespace

#endif

