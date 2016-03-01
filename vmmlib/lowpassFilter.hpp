/*
 * Copyright (c) 2006-2016, Visualization and Multimedia Lab,
 *                          University of Zurich <http://vmml.ifi.uzh.ch>,
 *                          Eyescale Software GmbH,
 *                          Blue Brain Project, EPFL
 *
 * This file is part of VMMLib <https://github.com/VMML/vmmlib/>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution.  Neither the name of the Visualization and Multimedia
 * Lab, University of Zurich nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific prior
 * written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* @author Jafet Villafranca
 *
 * Implements low pass filtering on a templated data type, which needs to
 * implement multiplication by a scalar float for smoothing
 */

#ifndef __VMML__LOWPASS_FILTER__HPP__
#define __VMML__LOWPASS_FILTER__HPP__

#include <deque>

namespace vmml
{
template< size_t M, typename T > class LowpassFilter
{
public:
    /** Construct a new lowpass filter with the given smoothing. */
    LowpassFilter( const float F ) : _smoothFactor( F ) {}
    ~LowpassFilter() {}

    /** @return The current filtered output value */
    const T& get() const { return _value; }

    /** Access the filtered value. */
    const T* operator->() const { return &_value; }

    /** Access the filtered value. */
    const T& operator*() const { return _value; }

    /** Add a value to the data set and return the filtered output */
    const T& add( const T& value );

private:
    std::deque< T > _data;
    float _smoothFactor;
    T _value;
};


template< size_t M, typename T >
const T& LowpassFilter< M, T >::add( const T& value )
{
    _data.push_front( value );

    while( _data.size() > M )
        _data.pop_back();

    // update
    typename std::deque< T >::const_iterator i = _data.begin();
    _value = *i;
    double weight = _smoothFactor;

    for( ++i ; i != _data.end(); ++i )
    {
        _value = _value * (1 - weight) + (*i) * weight;
        weight *= _smoothFactor;
    }

    return _value;
}

} // namespace vmml

#endif
