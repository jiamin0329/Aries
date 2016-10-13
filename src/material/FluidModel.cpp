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

#include "FluidModel.hpp"


namespace ARIES
{
    FluidModel::FluidModel() 
    {
        d_staticEnergy = 0.0;
        d_entropy = 0.0;
        d_density = 0.0;
        d_pressure = 0.0;
        d_soundSpeed2 = 0.0;
        d_temperature = 0.0;
        d_dPdrho_e = 0.0;
        d_dPde_rho = 0.0;
        d_dTdrho_e = 0.0;
        d_dTde_rho = 0.0;
        d_cp = 0.0;

        d_laminarViscosity = NULL;
        d_thermalConductivity = NULL;
    }

    FluidModel::~FluidModel() 
    {
    }

    double FluidModel::GetLaminarViscosity() 
    {
        d_laminarViscosity->SetViscosity(d_temperature, d_density);
        d_mu = d_laminarViscosity->GetViscosity();
        d_laminarViscosity->SetDerViscosity(d_temperature, d_density);
        d_dmudrho_T = d_laminarViscosity->Getdmudrho_T();
        d_dmudT_rho = d_laminarViscosity->GetdmudT_rho();
        return d_mu;
    }

    void FluidModel::SetLaminarViscosity(IProcData *procData) 
    {
        switch (procData->GetVisModelType()) 
        {
            case CONSTANT_VISCOSITY:
                d_laminarViscosity = new ViscosityConstant(procData->GetMuConstantND());
                break;
            case SUTHERLAND:
                d_laminarViscosity = new ViscositySutherland(procData->GetMuRefND(),
                                                             procData->GetMuTrefND(),
                                                             procData->GetMuSND());
                break;
        }
    }

    double FluidModel::GetThermalConductivity() 
    {
        d_thermalConductivity->SetConductivity(d_temperature, d_density, d_mu, d_cp);
        d_kt = d_thermalConductivity->GetConductivity();
        d_thermalConductivity->SetDerConductivity(d_temperature, d_density, d_dmudrho_T, d_dmudT_rho, d_cp);
        d_dktdrho_T = d_thermalConductivity->Getdktdrho_T();
        d_dktdT_rho = d_thermalConductivity->GetdktdT_rho();
        return d_kt;
    }
    
    void FluidModel::SetThermalConductivity(IProcData *procData)
    {
        switch (procData->GetConModelType()) 
        {
            case CONSTANT_CONDUCTIVITY:
                d_thermalConductivity = new ConductivityConstant(procData->GetKtConstantND());
                break;
            case CONSTANT_PRANDTL:
                d_thermalConductivity = new ConductivityConstantPrandtl(procData->GetPrandtlLam());
                break;
        }
    }

    /*
     *  Ideal Gas Model
     */
    FluidIdealGas::FluidIdealGas(): FluidModel()
    {
        d_gamma = 0.0;
        d_gammaM1 = 0.0;
        d_gasConstant = 0.0;
        d_cp = 0.0;
    }

    FluidIdealGas::FluidIdealGas(double val_gamma, double val_R): FluidModel()
    {
        d_gamma = val_gamma;
        d_gammaM1 = d_gamma - 1.0;
        d_gasConstant = val_R;
        d_cp = d_gamma/ d_gammaM1*d_gasConstant;
    }

    FluidIdealGas::~FluidIdealGas(void)
    {
    }

    void FluidIdealGas::SetTDState_rhoe(double val_rho, double val_e)
    {
        d_density = val_rho;
        d_staticEnergy = val_e;
        d_pressure = d_gammaM1*d_density*d_staticEnergy;
        d_temperature = d_gammaM1*d_staticEnergy/d_gasConstant;
        d_soundSpeed2 = d_gamma*d_pressure/d_density;
        d_entropy = (1.0/d_gammaM1*log(d_temperature) + log(1.0/d_density))*d_gasConstant;
        d_dPdrho_e = d_gammaM1*d_staticEnergy;
        d_dPde_rho = d_gammaM1*d_density;
        d_dTdrho_e = 0.0;
        d_dTde_rho = d_gammaM1/d_gasConstant;
    }

    void FluidIdealGas::SetTDState_PT(double val_P, double val_T) 
    {
        double e = val_T*d_gasConstant/d_gammaM1;
        double rho = val_P/(val_T*d_gasConstant);
        SetTDState_rhoe(rho, e);
    }

    void FluidIdealGas::SetTDState_Prho(double val_P, double val_rho) 
    {
        double e = val_P/(d_gammaM1*val_rho);
        SetTDState_rhoe(val_rho, e);
    }

    void FluidIdealGas::SetEnergy_Prho(double val_P, double val_rho) 
    {
        d_staticEnergy = val_P / (val_rho*d_gammaM1);
    }

    void FluidIdealGas::SetTDState_hs(double val_h, double val_s) 
    {
        double T = val_h*d_gammaM1/d_gasConstant/d_gamma;
        double e = val_h/d_gamma;
        double v = exp(-1.0/d_gammaM1*log(T) + val_s/d_gasConstant);
        SetTDState_rhoe(1.0/v, e);
    }

    void FluidIdealGas::SetTDState_rhoT(double val_rho, double val_T) 
    {
        double e = val_T*d_gasConstant/d_gammaM1;
        SetTDState_rhoe(val_rho, e);
    }

    /*
     *  van der waal gas
     */
    FluidVanDerWaalsGas::FluidVanDerWaalsGas(): FluidIdealGas()
    {
        d_a = 0.0;
        d_b = 0.0;
    }

    FluidVanDerWaalsGas::FluidVanDerWaalsGas(double val_gamma, double val_R, double val_Pstar, double val_Tstar): FluidIdealGas(val_gamma, val_R)
    {
        d_a = 27.0 / 64.0*d_gasConstant*d_gasConstant*val_Tstar*val_Tstar / val_Pstar;
        d_b = 1.0 / 8.0*d_gasConstant*val_Tstar / val_Pstar;
        d_zed = 1.0;
    }

    FluidVanDerWaalsGas::~FluidVanDerWaalsGas()
    {
    }

    void FluidVanDerWaalsGas::SetTDState_rhoe(double val_rho, double val_e)
    {
        d_density = val_rho;
        d_staticEnergy = val_e;

        d_pressure = d_gammaM1*d_density/(1.0 - d_density*d_b)*(d_staticEnergy + d_density*d_a) - d_a*d_density*d_density;
        d_temperature = (d_pressure + d_density*d_density*d_a)*((1.0 - d_density*d_b) / (d_density*d_gasConstant));
        d_entropy = d_gasConstant *(log(d_temperature) / d_gammaM1 + log(1.0 / d_density - d_b));


        d_dPde_rho = d_density*d_gammaM1 / (1.0 - d_density*d_b);
        d_dPdrho_e = d_gammaM1 / (1.0 - d_density*d_b)*((d_staticEnergy + 2.0 * d_density*d_a) + d_density*d_b*(d_staticEnergy + d_density*d_a) / (1.0 - d_density*d_b)) - 2.0 * d_density*d_a;
        d_dTdrho_e = d_gammaM1 / d_gasConstant*d_a;
        d_dTde_rho = d_gammaM1 / d_gasConstant;

        d_soundSpeed2 = d_dPdrho_e + d_pressure / (d_density*d_density)*d_dPde_rho;

        d_zed = d_pressure / (d_gasConstant*d_temperature*d_density);
    }

    void FluidVanDerWaalsGas::SetTDState_PT(double val_P, double val_T)
    {
        double toll = 1.0e-5;
        unsigned short nmax = 20, count = 0;
        double A, B, Z, DZ = 1.0, F, F1;
        A = d_a*val_P / (val_T*d_gasConstant) / (val_T*d_gasConstant);
        B = d_b*val_P / (val_T*d_gasConstant);

        if (d_zed > 0.1)
            Z = min(d_zed, 0.99);
        else
            Z = 0.99;

        do
        {
            F = Z*Z*Z - Z*Z*(B + 1.0) + Z*A - A*B;
            F1 = 3 * Z*Z - 2 * Z*(B + 1.0) + A;
            DZ = F / F1;
            Z -= 0.7*DZ;
            count++;
        } while (abs(DZ)>toll && count < nmax);

        if (count == nmax)
        {
            cout << "Warning Newton-Raphson exceed number of max iteration in PT" << endl;
            cout << "Compressibility factor  " << Z << " would be substituted with " << d_zed << endl;
        }

        // check if the solution is physical otherwise uses previous point solution
        if (Z <= 1.01 && Z >= 0.05 && count < nmax)
            d_zed = Z;

        d_density = val_P/(d_zed*d_gasConstant*val_T);
        double e = val_T*d_gasConstant / d_gammaM1 - d_a*d_density;
        SetTDState_rhoe(d_density, e);
    }

    void FluidVanDerWaalsGas::SetTDState_Prho(double val_P, double val_rho)
    {
        SetEnergy_Prho(val_P, val_rho);
        SetTDState_rhoe(val_rho, d_staticEnergy);
    }

    void FluidVanDerWaalsGas::SetTDState_hs(double val_h, double val_s)
    {
        double v, T, rho, f, fmid, rtb;
        double x1, x2, xmid, dx, fx1, fx2;
        double toll = 1e-5, FACTOR = 0.2;
        unsigned short count = 0, NTRY = 10, ITMAX = 100;
        double cons_s, cons_h;

        T = 1.0*val_h*d_gammaM1 / d_gasConstant / d_gamma;
        v = exp(-1.0 / d_gammaM1*log(T) + val_s / d_gasConstant);
        if (d_zed<0.9999)
        {
            x1 = d_zed*v;
            x2 = v;
        }
        else
        {
            x1 = 0.5*v;
            x2 = v;
        }
        
        fx1 = log(x1 - d_b) - val_s / d_gasConstant + log((val_h + 2.0 * d_a / x1) / d_gasConstant / (1.0 / d_gammaM1 + x1 / (x1 - d_b))) / d_gammaM1;
        fx2 = log(x2 - d_b) - val_s / d_gasConstant + log((val_h + 2.0 * d_a / x2) / d_gasConstant / (1.0 / d_gammaM1 + x2 / (x2 - d_b))) / d_gammaM1;

        // zbrac algorithm NR
        for (int j = 1; j <= NTRY; j++)
        {
            if (fx1*fx2 > 0.0)
            {
                if (fabs(fx1) < fabs(fx2))
                {
                    x1 += FACTOR*(x1 - x2);
                    fx1 = log(x1 - d_b) - val_s / d_gasConstant + log((val_h + 2.0 * d_a / x1) / d_gasConstant / (1.0 / d_gammaM1 + x1 / (x1 - d_b))) / d_gammaM1;
                }
                else
                {
                    x2 += FACTOR*(x2 - x1);
                    fx2 = log(x2 - d_b) - val_s / d_gasConstant + log((val_h + 2.0 * d_a / x2) / d_gasConstant / (1.0 / d_gammaM1 + x2 / (x2 - d_b))) / d_gammaM1;
                }
            }
        }
        
        // rtbis algorithm NR
        f = fx1;
        fmid = fx2;
        if (f*fmid >= 0.0)
        {
            cout << "Root must be bracketed for bisection in rtbis" << endl;
            SetTDState_rhoT(d_density, d_temperature);
        }
        rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
        do
        {
            xmid = rtb + (dx *= 0.5);
            fmid = log(xmid - d_b) - val_s / d_gasConstant + log((val_h + 2 * d_a / xmid) / d_gasConstant / (1 / d_gammaM1 + xmid / (xmid - d_b))) / d_gammaM1;
            if (fmid <= 0.0)
                rtb = xmid;
            count++;
        } while (abs(fmid) > toll && count<ITMAX);

        v = xmid;
        if (count == ITMAX)
        {
            cout << "Too many bisections in rtbis" << endl;
        }

        rho = 1 / v;
        T = (val_h + 2.0 * d_a / v) / d_gasConstant / (1.0 / d_gammaM1 + v / (v - d_b));
        SetTDState_rhoT(rho, T);

        cons_h = abs(((d_staticEnergy + d_pressure / d_density) - val_h) / val_h);
        cons_s = abs((d_entropy - val_s) / val_s);

        if (cons_h >1e-3 || cons_s >1e-3)
        {
            cout << "TD consistency not verified in hs call" << endl;
        }
    }

    void FluidVanDerWaalsGas::SetEnergy_Prho(double val_P, double val_rho)
    {
        double T = (val_P + val_rho*val_rho*d_a)*(1.0 - val_rho*d_b)/(val_rho*d_gasConstant);
        d_staticEnergy = T*d_gasConstant/d_gammaM1-val_rho*d_a;
    }

    void FluidVanDerWaalsGas::SetTDState_rhoT(double val_rho, double val_T)
    {
        double e = val_T*d_gasConstant / d_gammaM1 - d_a*val_rho;
        SetTDState_rhoe(val_rho, e);
    }
   
    FluidPengRobinson::FluidPengRobinson(): FluidIdealGas() 
    {
        d_a = 0.0;
        d_b = 0.0;
        d_k = 0.0;
        d_TstarCrit = 0.0;
    }

    FluidPengRobinson::FluidPengRobinson(double val_gamma, double val_R, double val_Pstar, double val_Tstar, double val_w) : FluidIdealGas(val_gamma, val_R)
    {
        d_a = 0.45724*d_gasConstant*d_gasConstant*val_Tstar*val_Tstar / val_Pstar;
        d_b = 0.0778*d_gasConstant*val_Tstar/val_Pstar;
        d_TstarCrit = val_Tstar;
        d_zed = 1.0;

        if (val_w <= 0.49)
            d_k = 0.37464 + 1.54226*val_w - 0.26992*val_w*val_w;
        else
            d_k = 0.379642 + 1.48503*val_w - 0.164423*val_w*val_w + 0.016666*val_w*val_w*val_w;
    }

    FluidPengRobinson::~FluidPengRobinson()
    {
    }

    double FluidPengRobinson::alpha2(double val_T)
    {
        return (1.0 + d_k*(1 - sqrt(val_T / d_TstarCrit)))*(1 + d_k*(1 - sqrt(val_T / d_TstarCrit)));
    }

    double FluidPengRobinson::T_v_h(double val_v, double val_h)
    {
        double fv, A, B, C, T, d;
        double sqrt2 = sqrt(2.0);

        d = (val_v*val_v + 2 * d_b*val_v - d_b*d_b);
        fv = atanh(d_b*sqrt2 / (val_v + d_b));

        A = d_gasConstant*(1 / d_gammaM1 + val_v / (val_v - d_b)) - d_a*val_v*d_k*d_k / (d_TstarCrit * d);
        B = d_a*d_k*(d_k + 1) / sqrt(d_TstarCrit) *(fv / (d_b*sqrt2) + 2 * val_v / d);
        C = val_h + d_a*(1 + d_k)*(1 + d_k)*(fv / (d_b*sqrt2) + val_v / d);

        T = (-B + sqrt(B*B + 4 * A*C)) / (2 * A); // Only positive root considered

        return T*T;
    }

    void FluidPengRobinson::SetTDState_rhoe(double val_rho, double val_e)
    {
        double DpDd_T, DpDT_d, DeDd_T, Cv;
        double A, B, C, sqrt2, fv, a2T, rho2;

        d_density = val_rho;
        d_staticEnergy = val_e;

        rho2 = val_rho*val_rho;
        sqrt2 = sqrt(2);

        fv = atanh(val_rho * d_b * sqrt2 / (1 + val_rho*d_b));

        A = d_gasConstant / d_gammaM1;
        B = d_a*d_k*(d_k + 1)*fv / (d_b*sqrt2*sqrt(d_TstarCrit));
        C = d_a*(d_k + 1)*(d_k + 1)*fv / (d_b*sqrt2) + val_e;

        d_temperature = (-B + sqrt(B*B + 4 * A*C)) / (2 * A); /// Only positive root considered
        d_temperature *= d_temperature;

        a2T = alpha2(d_temperature);

        A = (1 / rho2 + 2 * d_b / val_rho - d_b*d_b);
        B = 1 / val_rho - d_b;

        d_pressure = d_temperature*d_gasConstant / B - d_a*a2T / A;

        d_entropy = d_gasConstant / d_gammaM1*log(d_temperature) + d_gasConstant*log(B) - d_a*sqrt(a2T)*d_k*fv / (d_b*sqrt2*sqrt(d_temperature*d_TstarCrit));

        DpDd_T = (d_temperature*d_gasConstant / (B*B) - 2 * d_a*a2T*(1 / val_rho + d_b) / (A*A)) / (rho2);

        DpDT_d = d_gasConstant / B + d_a*d_k / A * sqrt(a2T / (d_temperature*d_TstarCrit));

        Cv = d_gasConstant / d_gammaM1 + (d_a*d_k*(d_k + 1)*fv) / (2 * d_b*sqrt(2 * d_temperature*d_TstarCrit));

        d_dPde_rho = DpDT_d / Cv;

        DeDd_T = -d_a*(1 + d_k) * sqrt(a2T) / A / (rho2);

        d_dPdrho_e = DpDd_T - d_dPde_rho*DeDd_T;

        d_soundSpeed2 = d_dPdrho_e + d_pressure / (rho2)*d_dPde_rho;

        d_dTde_rho = 1 / Cv;

        d_zed = d_pressure / (d_gasConstant*d_temperature*d_density);
    }
    
    void FluidPengRobinson::SetTDState_PT(double val_P, double val_T)
    {
        double toll = 1e-6;
        double A, B, Z, DZ = 1.0, F, F1;
        double rho, fv, e;
        double sqrt2 = sqrt(2);
        unsigned short nmax = 20, count = 0;

        A = d_a*alpha2(val_T)*val_P / (val_T*d_gasConstant) / (val_T*d_gasConstant);
        B = d_b*val_P / (val_T*d_gasConstant);

        if (d_zed > 0.1)
            Z = min(d_zed, 0.99);
        else
            Z = 0.99;
        
        do
        {
            F = Z*Z*Z + Z*Z*(B - 1.0) + Z*(A - 2 * B - 3 * B*B) + (B*B*B + B*B - A*B);
            F1 = 3 * Z*Z + 2 * Z*(B - 1.0) + (A - 2 * B - 3 * B*B);
            DZ = F / F1;
            Z -= DZ;
        } while (abs(DZ)>toll && count < nmax);

        if (count == nmax)
        {
            cout << "Warning Newton-Raphson exceed number of max iteration in PT" << endl;
            cout << "Compressibility factor  " << Z << " would be substituted with " << d_zed << endl;
        }
        // check if the solution is physical otherwise uses previous point  solution
        if (Z <= 1.0001 && Z >= 0.05 && count < nmax)
            d_zed = Z;

        rho = val_P / (d_zed*d_gasConstant*val_T);
        fv = atanh(rho * d_b * sqrt2 / (1 + rho*d_b));

        e = val_T*d_gasConstant / d_gammaM1 - d_a*(d_k + 1)*sqrt(alpha2(val_T))*fv / (d_b*sqrt2);

        SetTDState_rhoe(rho, e);
    }

    void FluidPengRobinson::SetTDState_Prho(double val_P, double val_rho)
    {
        SetEnergy_Prho(val_P, val_rho);
        SetTDState_rhoe(val_rho, d_staticEnergy);
    }

    void FluidPengRobinson::SetTDState_hs(double val_h, double val_s)
    {
        double T, fv, sqrt2 = sqrt(2.0), A;
        double f, v;
        double x1, x2, xmid, dx, fx1, fx2, fmid, rtb;
        double toll = 1e-9, FACTOR = 0.2;
        double cons_s, cons_h;
        unsigned short countrtb = 0, NTRY = 10, ITMAX = 100;

        A = d_gasConstant / d_gammaM1;
        T = val_h*d_gammaM1 / d_gasConstant / d_gamma;
        v = exp(-1 / d_gammaM1*log(T) + val_s / d_gasConstant);

        if (d_zed<0.9999)
        {
            x1 = d_zed*v;
            x2 = v;
        }
        else
        {
            x1 = 0.2*v;
            x2 = v;
        }

        T = T_v_h(x1, val_h);
        fv = atanh(d_b*sqrt2 / (x1 + d_b));
        fx1 = A*log(T) + d_gasConstant*log(x1 - d_b) - d_a*sqrt(alpha2(T)) *d_k*fv / (d_b*sqrt2*sqrt(T*d_TstarCrit)) - val_s;
        T = T_v_h(x2, val_h);
        fv = atanh(d_b*sqrt2 / (x2 + d_b));
        fx2 = A*log(T) + d_gasConstant*log(x2 - d_b) - d_a*sqrt(alpha2(T)) *d_k*fv / (d_b*sqrt2*sqrt(T*d_TstarCrit)) - val_s;

        // zbrac algorithm NR
        for (int j = 1; j <= NTRY; j++)
        {
            if (fx1*fx2 > 0.0)
            {
                if (fabs(fx1) < fabs(fx2))
                {
                    x1 += FACTOR*(x1 - x2);
                    T = T_v_h(x1, val_h);
                    fv = atanh(d_b*sqrt2 / (x1 + d_b));
                    fx1 = A*log(T) + d_gasConstant*log(x1 - d_b) - d_a*sqrt(alpha2(T)) *d_k*fv / (d_b*sqrt2*sqrt(T*d_TstarCrit)) - val_s;
                }
                else
                {
                    x2 += FACTOR*(x2 - x1);
                    T = T_v_h(x2, val_h);
                    fv = atanh(d_b*sqrt2 / (x2 + d_b));
                    fx2 = A*log(T) + d_gasConstant*log(x2 - d_b) - d_a*sqrt(alpha2(T)) *d_k*fv / (d_b*sqrt2*sqrt(T*d_TstarCrit)) - val_s;
                }
            }
        }

        // rtbis algorithm NR
        f = fx1;
        fmid = fx2;
        if (f*fmid >= 0.0) {
            cout << "Root must be bracketed for bisection in rtbis" << endl;
            SetTDState_rhoT(d_density, d_temperature);
        }
        rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
        do{
            xmid = rtb + (dx *= 0.5);
            T = T_v_h(xmid, val_h);
            fv = atanh(d_b* sqrt2 / (xmid + d_b));
            fmid = A*log(T) + d_gasConstant*log(xmid - d_b) - d_a*sqrt(alpha2(T)) *d_k*fv / (d_b*sqrt2*sqrt(T*d_TstarCrit)) - val_s;

            if (fmid <= 0.0)
                rtb = xmid;
            countrtb++;
        } while (abs(fmid) > toll && countrtb<ITMAX);

        v = xmid;
        if (countrtb == ITMAX)
        {
            cout << "Too many bisections in rtbis" << endl;
        }
            
        if (v != v)
        {
            cout << "not physical solution found, h and s input " << val_h << " " << val_s << endl;
            SetTDState_rhoT(d_density, d_temperature);
        }

        T = T_v_h(v, val_h);
        SetTDState_rhoT(1/v, T);

        // consistency check
        cons_h = abs(((d_staticEnergy + d_pressure / d_density) - val_h) / val_h);
        cons_s = abs((d_entropy - val_s) / val_s);

        if (cons_h >1e-4 || cons_s >1e-4)
        {
            cout << "TD consistency not verified in hs call" << endl;

        }
    }

    void FluidPengRobinson::SetEnergy_Prho(double val_P, double val_rho)
    {
        double ad;
        double A, B, C, T, vb1, vb2;
        vb1 = (1 / val_rho - d_b);
        vb2 = (1 / val_rho / val_rho + 2 * d_b / val_rho - d_b*d_b);

        A = d_gasConstant / vb1 - d_a*d_k*d_k / d_TstarCrit / vb2;

        B = 2 * d_a*d_k*(d_k + 1) / sqrt(d_TstarCrit) / vb2;

        C = -val_P - d_a*(1 + d_k)*(1 + d_k) / vb2;

        T = (-B + sqrt(B*B - 4 * A*C)) / (2 * A);
        T *= T;

        ad = d_a*(d_k + 1)*sqrt(alpha2(T)) / (d_b*sqrt(2)) * atanh(val_rho * d_b * sqrt(2) / (1 + val_rho*d_b));

        d_staticEnergy = T * d_gasConstant / d_gammaM1 - ad;

    }

    void FluidPengRobinson::SetTDState_rhoT(double val_rho, double val_T)
    {
        double fv, e;

        fv = atanh(val_rho * d_b * sqrt(2) / (1 + val_rho*d_b));
        e = val_T*d_gasConstant / d_gammaM1 - d_a*(d_k + 1)*sqrt(alpha2(val_T)) / (d_b*sqrt(2)) * fv;
        SetTDState_rhoe(val_rho, e);
    }
}
