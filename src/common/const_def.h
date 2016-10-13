/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Constants used in ARIES are defined
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */


#define MAX_STRING_SIZE  200      

#define MASTER_NODE                       0
#define SINGLE_NODE                       1

#define EPS                         1.0E-16	

#define PI                          3.1415926

#define VERTEX                            1   	 /*!< \brief VTK nomenclature for defining a vertex element. */
#define LINE                              3		 /*!< \brief VTK nomenclature for defining a line element. */
#define TRIANGLE                          5 	 /*!< \brief VTK nomenclature for defining a triangle element. */
#define RECTANGLE                         9		 /*!< \brief VTK nomenclature for defining a rectangle element. */
#define TETRAHEDRON                      10      /*!< \brief VTK nomenclature for defining a tetrahedron element. */
#define HEXAHEDRON                       12      /*!< \brief VTK nomenclature for defining a hexahedron element. */
#define PRISM                            13      /*!< \brief VTK nomenclature for defining a prism element. */
#define PYRAMID                          14  	 /*!< \brief VTK nomenclature for defining a pyramid element. */



#define N_ELEM_TYPES                      7

#define N_POINTS_LINE                     2
#define N_POINTS_TRIANGLE                 3
#define N_POINTS_QUADRILATERAL            4
#define N_POINTS_TETRAHEDRON              4
#define N_POINTS_HEXAHEDRON               8
#define N_POINTS_PYRAMID                  5
#define N_POINTS_PRISM                    6


#define BUFSIZE  3000000		                   /*!< \brief MPI buffer. */
#define MAX_PARAMETERS  10		                   /*!< \brief Maximum number of parameters for a design variable definition. */
#define MAX_NUMBER_PERIODIC  10                   /*!< \brief Maximum number of periodic boundary conditions. */
#define MAX_STRING_SIZE 200                      /*!< \brief Maximum number of domains. */
#define MAX_NUMBER_FFD  10	                       /*!< \brief Maximum number of FFDBoxes for the FFD. */
#define MAX_SOLS  6		                       /*!< \brief Maximum number of solutions at the same time (dimension of solution container array). */
#define MAX_TERMS  6		                       /*!< \brief Maximum number of terms in the numerical equations (dimension of solver container array). */
#define MAX_ZONES  3                             /*!< \brief Maximum number of zones. */
#define NO_RK_ITER  0		                       /*!< \brief No Runge-Kutta iteration. */

#define X_AXIS  0      /*!< \brief X axis orientation. */
#define Y_AXIS  1 	   /*!< \brief Y axis orientation. */
#define Z_AXIS  2      /*!< \brief Z axis orientation. */


