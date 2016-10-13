/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for fluid models
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    27-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_FLUIDMODEL_HPP
#define ARIES_FLUIDMODEL_HPP

#include "IProcData.hpp"

#include <cmath>

#include "ViscosityModel.hpp"
#include "ConductivityModel.hpp"

namespace ARIES
{
    /*!
     * \brief Main class for defining the Thermo-Physical Model
     */
    class FluidModel 
    {
    public:
        FluidModel();
        virtual ~FluidModel();
 
        double GetPressure() { return d_pressure; };
        double GetTemperature() { return d_temperature; };
        double GetEntropy() { return d_entropy; };
        double GetStaticEnergy() { return d_staticEnergy; };
        double GetDensity() { return d_density; };
        double GetSoundSpeed() { return sqrt(d_soundSpeed2); };
        double GetSoundSpeed2() { return d_soundSpeed2; };
        double GetCp() { return d_cp; };
        
        double GetLaminarViscosity();
        void SetLaminarViscosity(IProcData *procData);
        
        double GetThermalConductivity();
        void SetThermalConductivity(IProcData *procData);
  
        double GetdPdrho_e() { return d_dPdrho_e; };
        double GetdPde_rho() { return d_dPde_rho; };
        double GetdTdrho_e() { return d_dTdrho_e; };
        double GetdTde_rho() { return d_dTde_rho; };
        double Getdmudrho_T() { return d_dmudrho_T; };
        double GetdmudT_rho() { return d_dmudT_rho; };
        double Getdktdrho_T() { return d_dktdrho_T; };
        double GetdktdT_rho() { return d_dktdT_rho; };
        
        virtual void SetTDState_rhoe(double val_rho, double val_e);
        virtual void SetTDState_PT(double val_P, double val_T);
        virtual void SetTDState_Prho(double val_P, double val_rho);
        virtual void SetEnergy_Prho(double val_P, double val_rho);
        virtual void SetTDState_hs(double val_h, double val_s);
        virtual void SetTDState_rhoT(double val_rho, double val_T);

    protected:
        double d_staticEnergy;			/*!< \brief Internal Energy. */
        double d_entropy;  				/*!< \brief Entropy. */
        double d_density;  				/*!< \brief Density. */
        double d_pressure; 				/*!< \brief Pressure. */
        double d_soundSpeed2; 			/*!< \brief SpeedSound. */
        double d_temperature;			/*!< \brief Temperature. */
        double d_dPdrho_e; 				/*!< \brief DpDd_e. */
        double d_dPde_rho; 				/*!< \brief DpDe_d. */
        double d_dTdrho_e; 				/*!< \brief DTDd_e. */
        double d_dTde_rho; 				/*!< \brief DTDe_d. */
        double d_cp;                    /*!< \brief Specific Heat Capacity at constant pressure. */
        double d_mu;					/*!< \brief Specific Heat Capacity at constant pressure. */
        double d_dmudrho_T; 			/*!< \brief Specific Heat Capacity at constant pressure. */
        double d_dmudT_rho;				/*!< \brief Specific Heat Capacity at constant pressure. */
        double d_kt;					/*!< \brief Specific Heat Capacity at constant pressure. */
        double d_dktdrho_T; 			/*!< \brief Specific Heat Capacity at constant pressure. */
        double d_dktdT_rho;				/*!< \brief Specific Heat Capacity at constant pressure. */

        ViscosityModel *d_laminarViscosity;	          /*!< \brief Laminar Viscosity Model */
        ConductivityModel *d_thermalConductivity;	  /*!< \brief Thermal Conductivity Model */
    };

    /*!
     * \brief Child class for defining ideal gas model.
     */
    class FluidIdealGas : public FluidModel 
    {
    public:
        FluidIdealGas();
        FluidIdealGas(double val_gamma, double val_R);
        virtual ~FluidIdealGas(void);

        void SetTDState_rhoe(double val_rho, double val_e);
        void SetTDState_PT(double val_P, double val_T);
        void SetTDState_Prho(double val_P, double val_rho);
        void SetEnergy_Prho(double val_P, double val_rho);
        void SetTDState_hs(double val_h, double val_s);
        void SetTDState_rhoT(double val_rho, double val_T);
        
    protected:
        double d_gamma;				/*!< \brief Heat Capacity Ratio. */
        double d_gammaM1; 			/*!< \brief Heat Capacity Ratio Minus One. */
        double d_gasConstant;		/*!< \brief Gas Constant. */
    };

    /*!
     * \brief Child class for defining the Van der Waals model.
     */
    class FluidVanDerWaalsGas: public FluidIdealGas 
    {
    public:
        FluidVanDerWaalsGas();
        FluidVanDerWaalsGas(double val_gamma, double val_R, double val_Pstar, double val_Tstar);
        virtual ~FluidVanDerWaalsGas(void);

        void SetTDState_rhoe(double val_rho, double val_e);
        void SetTDState_PT(double val_P, double val_T);
        void SetTDState_Prho(double val_P, double val_rho);
        void SetEnergy_Prho(double val_P, double val_rho);
        void SetTDState_hs(double val_h, double val_s);
        void SetTDState_rhoT(double val_rho, double val_T);
        
    protected:
        double d_a, d_b, d_zed;  /*!< \brief Parameters for the Dimensionless Equation. */
    };

    /*!
     * \brief Child class for defining the Peng-Robinson model.
     */
    class FluidPengRobinson : public FluidIdealGas 
    {
    public:
        FluidPengRobinson();
        FluidPengRobinson(double val_gamma, double val_R, double val_Pstar, double val_Tstar, double val_w);
        virtual ~FluidPengRobinson(void);

        void SetTDState_rhoe(double val_rho, double val_e);
        void SetTDState_PT(double val_P, double val_T);
        void SetTDState_Prho(double val_P, double val_rho);
        void SetEnergy_Prho(double val_P, double val_rho);
        void SetTDState_hs(double val_h, double val_s);
        void SetTDState_rhoT(double val_rho, double val_T);
        
    protected:
        double  alpha2(double val_T);
        double  T_v_h(double val_v, double val_h);
        
    protected:
        double d_a; 					/*!< \brief model parameter. */
        double d_b; 					/*!< \brief model parameter. */
        double d_k; 					/*!< \brief model parameter (computed with acentric factor). */
        double d_zed; 					/*!< \brief compressibility factor. */
        double d_TstarCrit;				/*!< \brief Critical temperature. */
    };
}

#endif














