/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for adjoint variable
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    24-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_IVARIABLEADJ_HPP
#define ARIES_IVARIABLEADJ_HPP


namespace ARIES
{
    class IVariableAdj
    {
    public:
        explicit IVariableAdj()
        {
        };

        virtual double* GetSolution() = 0;

    protected:
        virtual ~IVariableAdj()
        {
        };
        
    }; // end IVariable class
} // end ARIES namespace

#endif

