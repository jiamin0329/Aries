/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for conductivity propertis
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_CONDUCTIVITYMODEL_HPP
#define ARIES_CONDUCTIVITYMODEL_HPP

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#include "stdio.h"
#include "math.h"

namespace ARIES
{
    /*!
     * \class CThermalConductivityModel
     * \brief Main class for defining the Transport-Physical Model
     * \version 1.0
     */
    class ConductivityModel
    {
    public:
        ConductivityModel();
        virtual ~ConductivityModel();

        double GetConductivity() { return d_kt; };
        double Getdktdrho_T() { return d_dktdrho_T; };
        double GetdktdT_rho() { return d_dktdT_rho; };

        virtual	void SetConductivity(double T, double rho, double mu, double cp) = 0;
        virtual	void SetDerConductivity(double T, double rho, double dmudrho_T, double dmudT_rho, double cp) = 0;
        
    protected:
        double  d_kt;			/*!< \brief Thermal conductivity. */
        double  d_dktdrho_T; 	/*!< \brief DktDrho_T. */
        double  d_dktdT_rho; 	/*!< \brief DktDT_rho. */
    };


    /*!
     * \class MODEL_ConductivityConstantPrandtl
     * \brief this class defines a constant thermal conductivity using a constant Prandtl's number
     * \author S.Vitale, M.Pini
     * \version 1.0
     */
    class ConductivityConstant : public ConductivityModel 
    {
    public:
        ConductivityConstant();
        ConductivityConstant(double kt_const);
        virtual ~ConductivityConstant();

        void SetConductivity(double T, double rho, double mu, double cp) {};
        void SetDerConductivity(double T, double rho, double dmudrho_T, double dmudT_rho, double cp) {};
    };


    /*!
     * \class MODEL_ConductivityConstantPrandtl
     * \brief this class defines a non-constant thermal conductivity using a constant Prandtl's number
     * \author S.Vitale, M.Pini
     * \version 1.0
     */
    class ConductivityConstantPrandtl : public ConductivityModel 
    {
    public:
        ConductivityConstantPrandtl();
        ConductivityConstantPrandtl(double pr_const);
        virtual ~ConductivityConstantPrandtl();

        void SetConductivity(double T, double rho, double mu, double cp);
        void SetDerConductivity(double T, double rho, double dmudrho_T, double dmudT_rho, double cp);
        
    protected:
        double d_prConst;		/*!< \brief Prandtl's number. */
    };
}

#endif






