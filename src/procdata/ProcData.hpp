/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Class for storing processor data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    01-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_PROCDATA_HPP
#define ARIES_PROCDATA_HPP

#include "IProcData.hpp"

namespace ARIES
{

    class ProcData: virtual public IProcData
    {
    public:
        ProcData();
        ~ProcData();

        unsigned short GetNumDim();
        unsigned short GetNumZone();

        // identifying the types of files to be written
        bool IsWrtVolSol (unsigned short iZone){};
        bool IsWrtSurfSol(unsigned short iZone){};
        bool IsWrtCsvSol (unsigned short iZone){};
        bool IsWrtHalo (){};

        SolverType GetSolverType(unsigned short iZone){};
        OutputType GetOutputType(unsigned short iZone){};
        bool IsFsiSimulation(){ return d_isFsi; };
        
        string GetSurfFlowCoeffFileName(){};
        string GetSurfAdjCoeffFileName(){};
        string GetSurfLinCoeffFileName(){};
        
        UnsteadyType GetUnsteadyType(){};
        bool GetWriteUnsteady(){};
        
        unsigned short GetNumMarkers(){};
        bool IsMarkerPlotted(unsigned short iMarker){};
        BCType GetMarkerBCType(unsigned short iMarker){};
        unsigned short GetMarkerSendRecv(unsigned short iMarker){};



        bool GetSmoothNumGrid(){};

        bool GetGridMovement(){};

        
        
        MeasurementType GetSystemMeasurements(){};


        VisModelType GetVisModelType(){};
        ConModelType GetConModelType(){};

        double GetMuConstantND(){};
        double GetMuRefND(){};
        double GetMuTrefND(){};
        double GetMuSND(){};

        double GetKtConstantND(){};
        double GetPrandtlLam(){};




        AxisOriType GetAxisOriType(){};








        


    private:

        bool d_isFsi;

        
        string d_meshFileName;
        MeshFileType d_meshFileType;

        unsigned short d_numTimeInstance; // for time spectral method

        UnsteadyType d_unsteadyType;



        

    };
}

#endif

