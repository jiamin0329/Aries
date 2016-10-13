/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for transport properties
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "ViscosityModel.hpp"

namespace ARIES
{
    ViscosityModel::ViscosityModel() 
    {
        d_mu = 0.0;
        d_dmudrho_T = 0.0;
        d_dmudT_rho = 0.0;
    }

    ViscosityModel::~ViscosityModel()
    { 
    }

    ViscosityConstant::ViscosityConstant() : ViscosityModel() 
    { 
    }

    ViscosityConstant::ViscosityConstant(double val_muConst) : ViscosityModel() 
    {
        d_mu = val_muConst;
        d_dmudrho_T = 0.0;
        d_dmudT_rho = 0.0;
    }

    ViscosityConstant::~ViscosityConstant() 
    { 
    }

    ViscositySutherland::ViscositySutherland() : ViscosityModel() 
    {
        d_muRef = 0.0;
        d_Tref = 0.0;
        d_S = 0.0;
    }

    ViscositySutherland::ViscositySutherland(double val_muRef, double val_Tref, double val_S) : ViscosityModel() 
    {
        d_muRef = val_muRef;
        d_Tref = val_Tref;
        d_S = val_S;
    }

    ViscositySutherland::~ViscositySutherland() 
    { 
    }

    void ViscositySutherland::SetViscosity(double val_T, double val_rho) 
    {
        d_mu = d_muRef*pow((val_T /d_Tref), (3.0 / 2.0))*((d_Tref + d_S) / (val_T + d_S));
    }

    void ViscositySutherland::SetDerViscosity(double val_T, double val_rho)
    {
        d_dmudrho_T = 0.0;
        d_dmudT_rho = d_muRef*((3.0 / 2.0)*pow((val_T / d_Tref), (1.0 / 2.0))*((d_Tref + d_S) / (val_T + d_S))
                            - pow((val_T / d_Tref), (3.0 / 2.0))*(d_Tref + d_S) / (val_T + d_S) / (val_T + d_S));
    }
}
