/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    dcomplex class for old-style complex and new complex<double>
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */


#ifndef ARIES_COMPLEX_HPP
#define ARIES_COMPLEX_HPP

#include <complex>

/*!
 * @page Complex Type
 *
 * @brief dcomplex is a typedef to overcome C++ compiler issues with
 * the std::complex type.
 *
 * The std::complex type should be a template however some older C++ compilers
 * implement complex as a double complex.  dcomplex is used to hide this
 * platform issue behind a typedef.
 *
 * @internal NOTE: This should be removed when no longer required.
 *
 */

#ifndef LACKS_TEMPLATE_COMPLEX
typedef std::complex<double> dcomplex;
#else
typedef std::complex dcomplex;
#endif

#endif
