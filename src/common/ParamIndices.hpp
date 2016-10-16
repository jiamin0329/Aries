/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Aries param indices
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    11-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_PARAMINDICES_HPP
#define ARIES_PARAMINDICES_HPP

#define PARAM_paramindices_start             0

//@{
//! @name switches for time manager
#define PARAM_timemgr_start                  1000
#define PARAM_timemgr_print_exclusive        PARAM_timemgr_start+1
#define PARAM_timemgr_print_total            PARAM_timemgr_start+2
#define PARAM_timemgr_print_processor        PARAM_timemgr_start+3
#define PARAM_timemgr_print_max              PARAM_timemgr_start+4
#define PARAM_timemgr_print_summed           PARAM_timemgr_start+5
#define PARAM_timemgr_print_user             PARAM_timemgr_start+6
#define PARAM_timemgr_print_sys              PARAM_timemgr_start+7
#define PARAM_timemgr_print_wall             PARAM_timemgr_start+8
#define PARAM_timemgr_print_percentage       PARAM_timemgr_start+9
#define PARAM_timemgr_print_concurrent       PARAM_timemgr_start+10
#define PARAM_timemgr_print_tiemroverhead    PARAM_timemgr_start+11
#define PARAM_timemgr_print_threshold        PARAM_timemgr_start+12
#define PARAM_timemgr_end                    PARAM_timemgr_start+13
//@}


/*!
 * \brief different software components of ARIES
 */
typedef enum 
{
    ARIES_CFD = 1,	/*!< \brief Running the ARIES_CFD software. */
    ARIES_DEF = 2,	/*!< \brief Running the ARIES_DEF software. */
    ARIES_DOT = 3,	/*!< \brief Running the ARIES_DOT software. */
    ARIES_MSH = 4,	/*!< \brief Running the ARIES_MSH software. */
    ARIES_GEO = 5,	/*!< \brief Running the ARIES_GEO software. */
    ARIES_SOL = 6 	/*!< \brief Running the ARIES_SOL software. */
} PARAM_aries_component_t;

#define PARAM_aries_component 1


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
} PARAM_cfd_solver_type_t;

#define PARAM_cfd_solver_type 2






#define PARAM_paramindices_end                99999

#endif












