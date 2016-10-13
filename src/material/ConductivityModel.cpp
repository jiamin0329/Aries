/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for conductivity properties
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "ConductivityModel.hpp"

namespace ARIES
{
    ConductivityModel::ConductivityModel() 
    {
        d_kt = 0.0;
        d_dktdrho_T = 0.0;
        d_dktdT_rho = 0.0;
    }

    ConductivityModel::~ConductivityModel() 
    { 
    }

    ConductivityConstant::ConductivityConstant() : ConductivityModel() 
    { 
    }

    ConductivityConstant::ConductivityConstant(double val_ktConst) : ConductivityModel() 
    {
        d_kt = val_ktConst;
        d_dktdrho_T = 0.0;
        d_dktdT_rho = 0.0;
    }

    ConductivityConstant::~ConductivityConstant() 
    { 
    }

    ConductivityConstantPrandtl::ConductivityConstantPrandtl() : ConductivityModel()
    { 
    }

    ConductivityConstantPrandtl::ConductivityConstantPrandtl(double val_prConst) : ConductivityModel() 
    {
        d_prConst = val_prConst;
    }

    void ConductivityConstantPrandtl::SetConductivity(double val_T, double val_rho, double val_mu, double val_cp) 
    {
        d_kt = val_mu*val_cp / d_prConst;
    }

    void ConductivityConstantPrandtl::SetDerConductivity(double val_T, double val_rho, double val_dmudrho_T, double val_dmudT_rho, double val_cp) 
    {
        d_dktdrho_T = val_dmudrho_T*val_cp / d_prConst;
        d_dktdT_rho = val_dmudT_rho*val_cp / d_prConst;
    }

    ConductivityConstantPrandtl::~ConductivityConstantPrandtl() 
    { 
    }
}
