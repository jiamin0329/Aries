/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for processor data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_IPROCDATA_HPP
#define ARIES_IPROCDATA_HPP

#include <string>

using namespace std;

namespace ARIES
{
    /*!
     * \brief different solver types for the CFD component
     */
    typedef enum
    {
        NO_SOLVER = 0,			                    /*!< \brief Definition of no solver. */
        EULER = 1,				                    /*!< \brief Definition of the Euler's solver. */
        NAVIER_STOKES = 2,			                /*!< \brief Definition of the Navier-Stokes' solver. */
        RANS = 3,				                    /*!< \brief Definition of the Reynolds-averaged Navier-Stokes' (RANS) solver. */
        POISSON_EQUATION = 4,       	            /*!< \brief Definition of the poisson potential solver. */
        WAVE_EQUATION = 10,	                        /*!< \brief Definition of the wave solver. */
        HEAT_EQUATION = 29,						    /*!< \brief Definition of the heat solver. */
        LINEAR_ELASTICITY = 11,	                    /*!< \brief Definition of the FEA solver. */
        FLUID_STRUCTURE_EULER = 12,	                /*!< \brief Definition of the FEA solver. */
        FLUID_STRUCTURE_NAVIER_STOKES = 13,	        /*!< \brief Definition of the FEA solver. */
        FLUID_STRUCTURE_RANS = 14,	                /*!< \brief Definition of the FEA solver. */
        ADJ_EULER = 18,			                    /*!< \brief Definition of the continuous adjoint Euler's solver. */
        ADJ_NAVIER_STOKES = 19,		                /*!< \brief Definition of the continuous adjoint Navier-Stokes' solver. */
        ADJ_RANS = 20,				                /*!< \brief Definition of the continuous adjoint Reynolds-averaged Navier-Stokes' (RANS) solver. */
        LIN_EULER = 21,			                    /*!< \brief Definition of the linear Euler's solver. */
        LIN_NAVIER_STOKES = 22,		                /*!< \brief Definition of the linear Navier-Stokes' solver. */
        TEMPLATE_SOLVER = 30,                       /*!< \brief Definition of template solver. */
        TNE2_EULER = 31,
        TNE2_NAVIER_STOKES = 32,
        ADJ_TNE2_EULER = 33,
        ADJ_TNE2_NAVIER_STOKES = 34
    } SolverType;

    /*!
     * \brief type of solution output file formats
     */
    typedef enum
    {
        TECPLOT = 1,  		     /*!< \brief Tecplot format for the solution output. */
        TECPLOT_BINARY = 2,      /*!< \brief Tecplot binary format for the solution output. */
        FIELDVIEW = 3,  		 /*!< \brief FieldView format for the solution output. */
        FIELDVIEW_BINARY = 4,    /*!< \brief FieldView binary format for the solution output. */
        CSV = 5,			     /*!< \brief Comma-separated values format for the solution output. */
        CGNS_SOL = 6,  	     	 /*!< \brief CGNS format for the solution output. */
        PARAVIEW = 7  		     /*!< \brief Paraview format for the solution output. */
    } OutputType;

    /*!
     * \brief types of schemes for unsteady computations
     */
    typedef enum
    {
        STEADY = 0,             /*!< \brief A steady computation. */
        TIME_STEPPING = 1,		/*!< \brief Use a time stepping strategy for unsteady computations. */
        DT_STEPPING_1ST = 2,	/*!< \brief Use a dual time stepping strategy for unsteady computations (1st order). */
        DT_STEPPING_2ND = 3,	/*!< \brief Use a dual time stepping strategy for unsteady computations (2nd order). */
        ROTATIONAL_FRAME = 4,   /*!< \brief Use a rotational source term. */
        TIME_SPECTRAL = 5       /*!< \brief Use a time spectral source term. */
    } UnsteadyType;


    /*!
     * \brief types of boundary conditions
     */
    typedef enum
    {
        EULER_WALL = 1,		            /*!< \brief Boundary Euler wall definition. */
        FAR_FIELD = 2,		            /*!< \brief Boundary far-field definition. */
        SYMMETRY_PLANE = 3,   	        /*!< \brief Boundary symmetry plane definition. */
        INLET_FLOW = 4,		            /*!< \brief Boundary inlet flow definition. */
        OUTLET_FLOW = 5,		        /*!< \brief Boundary outlet flow definition. */
        PERIODIC_BOUNDARY = 6,	        /*!< \brief Periodic boundary definition. */
        NEARFIELD_BOUNDARY = 7,	        /*!< \brief Near-Field boundary definition. */
        ELECTRODE_BOUNDARY = 8,	        /*!< \brief Electrode boundary definition. */
        DIELEC_BOUNDARY = 9,	        /*!< \brief Dipoisson boundary definition. */
        CUSTOM_BOUNDARY = 10,           /*!< \brief custom boundary definition. */
        INTERFACE_BOUNDARY = 11,        /*!< \brief Domain interface boundary definition. */
        DIRICHLET = 12,		            /*!< \brief Boundary Euler wall definition. */
        NEUMANN = 13,		            /*!< \brief Boundary Neumann definition. */
        DISPLACEMENT_BOUNDARY = 14,		/*!< \brief Boundary displacement definition. */
        LOAD_BOUNDARY = 15,		        /*!< \brief Boundary Load definition. */
        FLOWLOAD_BOUNDARY = 16,		    /*!< \brief Boundary Load definition. */
        ELEC_DIELEC_BOUNDARY = 17,	    /*!< \brief Dipoisson boundary definition for the poissonal potential. */
        ELEC_NEUMANN = 18,		        /*!< \brief Boundary Neumann definition. */
        SUPERSONIC_INLET = 19,		    /*!< \brief Boundary supersonic inlet definition. */
        SUPERSONIC_OUTLET = 20,		    /*!< \brief Boundary supersonic inlet definition. */
        ENGINE_INFLOW = 21,		        /*!< \brief Boundary nacelle inflow. */
        ENGINE_EXHAUST = 22,		    /*!< \brief Boundary nacelle exhaust. */
        ENGINE_BLEED = 23,		        /*!< \brief Boundary engine bleed. */
        RIEMANN_BOUNDARY = 24,          /*!< \brief Riemann Boundary definition. */
        ISOTHERMAL = 25,                /*!< \brief No slip isothermal wall boundary condition. */
        HEAT_FLUX = 26,                 /*!< \brief No slip constant heat flux wall boundary condition. */
        PRESSURE_BOUNDARY = 27,   	    /*!< \brief Pressure boundary condition. */
        HEAT_FLUX_NONCATALYTIC = 28,    /*!< \brief No-slip, constant heat flux, noncatalytic bc. */
        HEAT_FLUX_CATALYTIC = 29,       /*!< \brief No-slip, constant heat flux, catalytic bc. */
        ISOTHERMAL_NONCATALYTIC = 30,   /*!< \brief No-slip, constant temperature, noncatalytic bc. */
        ISOTHERMAL_CATALYTIC = 31,      /*!< \brief No-slip, constant temperature, catalytic bc. */
        ACTDISK_INLET = 32,	            /*!< \brief Actuator disk inlet boundary definition. */
        ACTDISK_OUTLET = 33,	        /*!< \brief Actuator disk outlet boundary definition. */
        SEND_RECEIVE = 99		        /*!< \brief Boundary send-receive definition. */
    } BCType;

    /*!
     * \brief different system of measurements
     */
    typedef enum
    {
        SI = 0,			        /*!< \brief Definition of  */
        US = 1,				    /*!< \brief Definition of  */
    } MeasurementType;

    /*!
     * \brief types of viscosity model
     */
    typedef enum
    {
        CONSTANT_VISCOSITY = 0,
        SUTHERLAND = 1
    } VisModelType;

    /*!
     * \brief types of thermal conductivity model
     */
    typedef enum
    {
        CONSTANT_CONDUCTIVITY = 0,
        CONSTANT_PRANDTL = 1
    } ConModelType;

    typedef enum
    {
        AxisOri_x = 0,     /*!< \brief X axis orientation. */
        AxisOri_y = 1, 	/*!< \brief Y axis orientation. */
        AxisOri_z = 2      /*!< \brief Z axis orientation. */
    } AxisOriType;

    typedef enum
    {
        MeshFile_none = 0,
        MeshFile_CGNS = 1,
        MeshFile_SU2  = 2
    } MeshFileType;
 
    class IProcData
    {

    public:
        virtual ~IProcData( )
        {		
        };

        virtual unsigned short GetNumDim() = 0;
        virtual unsigned short GetNumZone() = 0;

        // identifying the types of files to be written
        virtual	bool IsWrtVolSol (unsigned short iZone) = 0;
        virtual bool IsWrtSurfSol(unsigned short iZone) = 0;
        virtual bool IsWrtCsvSol (unsigned short iZone) = 0;
        virtual bool IsWrtHalo () = 0;

        virtual SolverType GetSolverType(unsigned short iZone) = 0;
        virtual OutputType GetOutputType(unsigned short iZone) = 0;
        virtual bool IsFsiSimulation() = 0;
        
        virtual string GetSurfFlowCoeffFileName() = 0;
        virtual string GetSurfAdjCoeffFileName() = 0;
        virtual string GetSurfLinCoeffFileName() = 0;
            
        virtual UnsteadyType GetUnsteadyType() = 0;
        virtual bool GetWriteUnsteady() = 0;
            
        virtual unsigned short GetNumMarkers() = 0;
        virtual bool IsMarkerPlotted(unsigned short iMarker) = 0;
        virtual BCType GetMarkerBCType(unsigned short iMarker) = 0;
        virtual unsigned short GetMarkerSendRecv(unsigned short iMarker) = 0;



        virtual bool GetSmoothNumGrid() = 0;

        virtual bool GetGridMovement() = 0;

            
            
        virtual MeasurementType GetSystemMeasurements() = 0;


        virtual VisModelType GetVisModelType() = 0;
        virtual ConModelType GetConModelType() = 0;

        virtual double GetMuConstantND() = 0;
        virtual double GetMuRefND() = 0;
        virtual double GetMuTrefND() = 0;
        virtual double GetMuSND() = 0;

        virtual double GetKtConstantND() = 0;
        virtual double GetPrandtlLam() = 0;




        virtual AxisOriType GetAxisOriType() = 0;
            
    protected:
        explicit IProcData()
        {
        };
    }; //end IProcData class
} // end ARIES namespace

#endif

