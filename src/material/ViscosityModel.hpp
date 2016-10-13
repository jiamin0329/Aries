/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for fluid viscosity model
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_VISCOSITYMODEL_HPP
#define ARIES_VISCOSITYMODEL_HPP

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#include "stdio.h"
#include "math.h"

namespace ARIES
{
    class ViscosityModel 
    {
    public:
        ViscosityModel();
        virtual ~ViscosityModel();

        double GetViscosity() { return d_mu; };
        double Getdmudrho_T() { return d_dmudrho_T; };
        double GetdmudT_rho() { return d_dmudT_rho; };

        virtual	void SetViscosity(double T, double rho) = 0;
        virtual	void SetDerViscosity(double T, double rho) = 0;
        
    protected:
        double d_mu;			/*!< \brief Dynamic viscosity. */
        double d_dmudrho_T; 	/*!< \brief DmuDrho_T. */
        double d_dmudT_rho; 	/*!< \brief DmuDT_rho. */
    };
    
    /*!
     * \brief this class defines a constant viscosity
     */
    class ViscosityConstant: public ViscosityModel 
    {
    public:
        ViscosityConstant();
        ViscosityConstant(double val_muConst);
        virtual ~ViscosityConstant();

        void SetViscosity(double T, double rho) {};
        void SetDerViscosity(double T, double rho) {};
    };

    /*!
     * \brief this class defines a constant viscosity
     */
    class ViscositySutherland: public ViscosityModel 
    {
    public:
        ViscositySutherland();
        ViscositySutherland(double mu_ref, double t_ref, double s);
        virtual ~ViscositySutherland();

        void SetViscosity(double T, double rho);
        void SetDerViscosity(double T, double rho);

    protected:
        double d_muRef;
        double d_Tref;	
        double d_S;
    };
}// end namespace ARIES

#endif






