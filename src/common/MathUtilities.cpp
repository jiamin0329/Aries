/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Utilities class to access common POSIX constants and math ops
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#ifndef ARIES_MATHUTILITIES_CPP
#define ARIES_MATHUTILITIES_CPP

#include <limits>
#include <vector>

namespace ARIES
{
    /*
     * Routines to initialize vectors and arrays to signaling NaNs.
     */
    template<class TYPE>
    void MathUtilities<TYPE>::SetVectorToSignalingNaN(std::vector<TYPE>& vector)
    {
        for (int i = 0; i < static_cast<int>(vector.size()); ++i)
        {
            vector[i] = GetSignalingNaN();
        }
    }

    template<class TYPE>
    void MathUtilities<TYPE>::SetArrayToSignalingNaN(TYPE* array, int n)
    {
        for (int i = 0; i < n; ++i)
        {
            array[i] = GetSignalingNaN();
        }
    }
    
    /*
     * Routines to initialize vectors and arrays to max value for type.
     */
    template<class TYPE>
    void MathUtilities<TYPE>::SetVectorToMax(std::vector<TYPE>& vector)
    {
        for (int i = 0; i < static_cast<int>(vector.size()); ++i)
        {
            vector[i] = GetMax();
        }
    }

    template<class TYPE>
    void MathUtilities<TYPE>::SetArrayToMax(TYPE* array, int n)
    {
        for (int i = 0; i < n; ++i)
        {
            array[i] = GetMax();
        }
    }
    
    /*
     * Routines to initialize vectors and arrays to min value for type.
     */
    template<class TYPE>
    void MathUtilities<TYPE>::SetVectorToMin(std::vector<TYPE>& vector)
    {
        for (int i = 0; i < static_cast<int>(vector.size()); ++i)
        {
            vector[i] = GetMin();
        }
    }

    template<class TYPE>
    void MathUtilities<TYPE>::SetArrayToMin(TYPE* array, int n)
    {
        for (int i = 0; i < n; ++i)
        {
            array[i] = GetMin();
        }
    }

    /*
     * Routines to initialize vectors and arrays to epsilon value for type.
     */
    template<class TYPE>
    void MathUtilities<TYPE>::SetVectorToEpsilon(std::vector<TYPE>& vector)
    {
        for (int i = 0; i < static_cast<int>(vector.size()); ++i)
        {
            vector[i] = GetEpsilon();
        }
    }

    template<class TYPE>
    void MathUtilities<TYPE>::SetArrayToEpsilon(TYPE* array, int n)
    {
        for (int i = 0; i < n; ++i)
        {
            array[i] = GetEpsilon();
        }
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::GetZero()
    {
        return d_zero;
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::GetOne()
    {
        return d_one;
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::GetSignalingNaN()
    {
        return std::numeric_limits<TYPE>::signaling_NaN();
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::GetMax()
    {
        return std::numeric_limits<TYPE>::max();
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::GetMin()
    {
        return std::numeric_limits<TYPE>::min();
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::GetEpsilon()
    {
        return std::numeric_limits<TYPE>::epsilon();
    }

    template<class TYPE>
    bool MathUtilities<TYPE>::IsNaN(const TYPE& value)
    {
        NULL_USE(value);
        return false;
    }

    template<class TYPE>
    bool MathUtilities<TYPE>::EqualEps(const TYPE& a, const TYPE& b)
    {
        return a == b;
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::Min(TYPE a, TYPE b)
    {
        return a < b ? a : b;
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::Max(TYPE a, TYPE b)
    {
        return a > b ? a : b;
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::Abs(TYPE value)
    {
        return value;
    }

    template<class TYPE>
    TYPE MathUtilities<TYPE>::Round(TYPE x)
    {
        return x;
    }
} // end namespace ARIES

#endif
