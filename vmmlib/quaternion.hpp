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

#ifndef __VMML__QUATERNION__HPP__
#define __VMML__QUATERNION__HPP__

#include <vmmlib/enable_if.hpp>
#include <vmmlib/types.hpp>
#include <vmmlib/vector.hpp> // used inline

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>


// - declaration - //

#define QUATERNION_TRACE_EPSILON 1e-5

namespace vmml
{

template < typename T > class quaternion
{
public:
    typedef vector< 3, T > vec3;

    /** Construct an identity quaternion */
    quaternion() : array() { array[3] = 1.; }
    quaternion( T x, T y, T z, T w );

    /** Construct a rotation quaternion */
    quaternion( T angle, vec3 axis );

    // uses the top-left 3x3 part of the supplied matrix as rotation matrix
    template< size_t M >
    quaternion( const Matrix< M, M, T >& rotation_matrix_,
        typename enable_if< M >= 3 >::type* = 0 );

    /** @return true if the two quaternion are similar. */
    bool equals( const quaternion& other,
                 T tolerance = std::numeric_limits< T >::epsilon( )) const;

    T x() const { return array[0]; }
    T y() const { return array[1]; }
    T z() const { return array[2]; }
    T w() const { return array[3]; }

    void zero();
    void identity();

    template< size_t D > void set( const Matrix< D, D, T >& rotation_matrix_ );

    bool operator==( const T& a ) const;
    bool operator!=( const T& a ) const;

    bool operator==( const quaternion& a ) const;
    bool operator!=( const quaternion& a ) const;

    bool is_akin( const quaternion& a,
                  const T& delta = std::numeric_limits< T >::epsilon() );

    void conjugate();
    quaternion get_conjugate() const;

    T abs() const;
    T squared_abs() const;

    T normalize();
    quaternion get_normalized() const;

    quaternion negate() const;
    quaternion operator-() const;

    quaternion& operator=(const quaternion& other);

    //
    // quaternion/quaternion operations
    //
    quaternion operator+( const quaternion< T >& a ) const;
    quaternion operator-( const quaternion< T >& a ) const;
    // caution: a * q != q * a in general
    quaternion operator*( const quaternion< T >& a ) const;
    void operator+=( const quaternion< T >& a );
    void operator-=( const quaternion< T >& a );
    // caution: a *= q != q *= a in general
    void operator*=( const quaternion< T >& a );

    //
    // quaternion/scalar operations
    //
    quaternion operator*( T a ) const;
    quaternion operator/( T a ) const;

    void operator*=( T a );
    void operator/=( T a );

    // vec3 = this x b
    vector< 3, T > cross( const quaternion< T >& b ) const;

    T dot( const quaternion< T >& a ) const;
    static T dot( const quaternion< T >& a, const quaternion< T >& b );

    // returns multiplicative inverse
    quaternion inverse();

    void normal( const quaternion& aa, const quaternion& bb, const quaternion& cc,  const quaternion& dd );
    quaternion normal( const quaternion& aa, const quaternion& bb, const quaternion& cc );

    static quaternion slerp( T a, const quaternion& p,
        const quaternion& q, const T epsilon = 1e-13 );

    Matrix< 3, 3, T > get_rotation_matrix() const;

    template< size_t D > void get_rotation_matrix( Matrix< D, D, T >& result ) const;

    static const quaternion IDENTITY;
    static const quaternion QUATERI;
    static const quaternion QUATERJ;
    static const quaternion QUATERK;

    friend std::ostream& operator<< ( std::ostream& os, const quaternion& q )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );

        os << "( "
           << std::setw(10) << q.x() << " "
           << std::setw(10) << q.y() << " "
           << std::setw(10) << q.z() << " "
           << std::setw(10) << q.w() << " "
           << ")";

        os.precision( prec );
        os.setf( flags );
        return os;
    }

private:
    T array[4];
};
}

#include <vmmlib/matrix.hpp>

namespace vmml
{
// - implementation - //

template < typename T >
const quaternion< T > quaternion< T >::IDENTITY( 0, 0, 0, 1 );

template < typename T >
const quaternion< T > quaternion< T >::QUATERI( 1, 0, 0, 0 );

template < typename T >
const quaternion< T > quaternion< T >::QUATERJ( 0, 1, 0, 0 );

template < typename T >
const quaternion< T > quaternion< T >::QUATERK( 0, 0, 1, 0 );

template < typename T >
quaternion< T >::quaternion( T x_, T y_, T z_, T w_ )
{
    array[0] = x_;
    array[1] = y_;
    array[2] = z_;
    array[3] = w_;
}


template< typename T > template< size_t M >
quaternion< T >::quaternion( const Matrix< M, M, T >& rotation_matrix_,
                             typename enable_if< M >= 3 >::type* )
{
    this->template set< M >( rotation_matrix_ );
}

template< typename T >
quaternion< T >::quaternion( const T angle, vec3 axis )
{
    axis.normalize();
    const T sinAngle2 = std::sin( angle * 0.5 );
    array[0] = axis.x() * sinAngle2;
    array[1] = axis.y() * sinAngle2;
    array[2] = axis.z() * sinAngle2;
    array[3] = std::cos( angle * 0.5 );
}

template< typename T >
bool quaternion< T >::equals( const quaternion& other, const T tolerance ) const
{
    return std::abs( array[0] - other.array[0] ) <= tolerance &&
           std::abs( array[1] - other.array[1] ) <= tolerance &&
          std::abs( array[2] - other.array[2] ) <= tolerance &&
          std::abs( array[3] - other.array[3] ) <= tolerance;
}


// top-left 3x3 is interpreted as rot matrix.
template < typename T > template< size_t D >
void quaternion< T >::set( const Matrix< D, D, T >& M )
{
    T trace = M( 0, 0 ) + M( 1, 1 ) + M( 2,2 ) + 1.0;

    // very small traces may introduce a big numerical error
    if( trace > QUATERNION_TRACE_EPSILON )
    {
        T s = 0.5 / std::sqrt( trace );
        array[0] = M( 2, 1 ) - M( 1, 2 );
        array[0] *= s;

        array[1] = M( 0, 2 ) - M( 2, 0 );
        array[1] *= s;

        array[2] = M( 1, 0 ) - M( 0, 1 );
        array[2] *= s;

        array[3] = 0.25 / s;
    }
    else
    {
        vector< 3, T > diag( M( 0, 0 ), M( 1, 1 ), M( 2, 2 ) );
        size_t largest = diag.find_max_index();

        // 0, 0 is largest
        if ( largest == 0 )
        {
            T s = 0.5 / std::sqrt( 1.0 + M( 0, 0 ) - M( 1, 1 ) - M( 2, 2 ) );
            array[0] = 0.25 / s;

            array[1] = M( 0,1 ) + M( 1,0 );
            array[1] *= s;

            array[2] = M( 0,2 ) + M( 2,0 );
            array[2] *= s;

            array[3] = M( 1,2 ) - M( 2,1 );
            array[3] *= s;
        }
        else if ( largest == 1 )
        {
            T s = 0.5 / std::sqrt( 1.0 + M( 1,1 ) - M( 0,0 ) - M( 2,2 ) );
            array[0] = M( 0,1 ) + M( 1,0 );
            array[0] *= s;

            array[1] = 0.25 / s;

            array[2] = M( 1,2 ) + M( 2,1 );
            array[2] *= s;

            array[3] = M( 0,2 ) - M( 2,0 );
            array[3] *= s;
        }
        // 2, 2 is largest
        else if ( largest == 2 )
        {
            T s = 0.5 / std::sqrt( 1.0 + M( 2,2 ) - M( 0,0 ) - M( 1,1 ) );
            array[0] = M( 0,2 ) + M( 2,0 );
            array[0] *= s;

            array[1] = M( 1,2 ) + M( 2,1 );
            array[1] *= s;

            array[2] = 0.25 / s;

            array[3] = M( 0,1 ) - M( 1,0 );
            array[3] *= s;
        }
        else
        {
            throw std::runtime_error( "no clue why, but should not get here" );
        }
    }
}



template < typename T >
void quaternion< T >::zero()
{
    ::memset( array, 0, sizeof( array ));
}



template < typename T >
void quaternion< T >::identity()
{
    (*this) = IDENTITY;
}

template < typename T >
bool quaternion< T >::operator==( const quaternion& rhs ) const
{
    return (
        array[3] == rhs.array[3] &&
        array[0] == rhs.array[0] &&
        array[1] == rhs.array[1] &&
        array[2] == rhs.array[2]
        );
}



template < typename T >
bool
quaternion< T >::operator!=( const quaternion& a ) const
{
    return ! this->operator==( a );
}

template < typename T >
void quaternion< T >::conjugate()
{
    array[0] = -array[0];
    array[1] = -array[1];
    array[2] = -array[2];
}



template < typename T >
quaternion< T > quaternion< T >::get_conjugate() const
{
    return quaternion< T > ( -array[0], -array[1], -array[2], array[3] );
}



template < typename T >
T
quaternion< T >::abs() const
{
    return std::sqrt( squared_abs() );
}



template < typename T >
T quaternion< T >::squared_abs() const
{
    return array[0] * array[0] + array[1] * array[1] + array[2] * array[2] + array[3] * array[3];
}



template < typename T >
quaternion< T > quaternion< T >::inverse()
{
    quaternion< T > q( *this );
    q.conjugate();

    T tmp = squared_abs();
    tmp = static_cast< T >( 1.0 ) / tmp;
    return q * tmp;
}



template < typename T >
T quaternion< T >::normalize()
{
    T len = abs();
    if( len == 0.0 )
        return 0.0;
    len = 1.0f / len;
    this->operator*=( len );
    return len;
}



template < typename T >
quaternion< T >
quaternion< T >::get_normalized() const
{
    quaternion< T > q( *this );
    q.normalize();
    return q;
}



//
// quaternion/quaternion operations
//

template < typename T >
quaternion< T >
quaternion< T >::operator+( const quaternion< T >& rhs ) const
{
    return quaternion( array[0] + rhs.array[0], array[1] + rhs.array[1], array[2] + rhs.array[2], array[3] + rhs.array[3] );
}



template < typename T >
quaternion< T >
quaternion< T >::operator-( const quaternion< T >& rhs ) const
{
    return quaternion( array[0] - rhs.array[0], array[1] - rhs.array[1], array[2] - rhs.array[2], array[3] - rhs.array[3] );
}



// returns Grasssmann product
template < typename T >
quaternion< T > quaternion< T >::operator*( const quaternion< T >& rhs ) const
{
    quaternion< T > ret( *this );
    ret *= rhs;
    return ret;
}



// Grassmann product
template < typename T >
void quaternion< T >::operator*=( const quaternion< T >& q )
{
    #if 0
    quaternion< T > orig( *this );
    array[0] = orig.array[3] * a.array[0] + orig.array[0] * a.array[3] + orig.array[1] * a.array[2] - orig.array[2] * a.array[1];
    array[1] = orig.array[3] * a.array[1] + orig.array[1] * a.array[3] + orig.array[2] * a.array[0] - orig.array[0] * a.array[2];
    array[2] = orig.array[3] * a.array[2] + orig.array[2] * a.array[3] + orig.array[0] * a.array[1] - orig.array[1] * a.array[0];
    array[3] = orig.array[3] * a.array[3] - orig.array[0] * a.array[0] - orig.array[1] * a.array[1] - orig.array[2] * a.array[2];
    #else

    // optimized version, 7 less mul, but 15 more add/subs
    // after Henrik Engstrom, from a gamedev.net article.

    const T& a_ = array[ 3 ];
    const T& b_ = array[ 0 ];
    const T& c = array[ 1 ];
    const T& d = array[ 2 ];
    const T& _x = q.array[ 3 ];
    const T& _y = q.array[ 0 ];
    const T& _z = q.array[ 1 ];
    const T& _w = q.array[ 2 ];

    const T tmp_00 = (d - c) * (_z - _w);
    const T tmp_01 = (a_ + b_) * (_x + _y);
    const T tmp_02 = (a_ - b_) * (_z + _w);
    const T tmp_03 = (c + d) * (_x - _y);
    const T tmp_04 = (d - b_) * (_y - _z);
    const T tmp_05 = (d + b_) * (_y + _z);
    const T tmp_06 = (a_ + c) * (_x - _w);
    const T tmp_07 = (a_ - c) * (_x + _w);
    const T tmp_08 = tmp_05 + tmp_06 + tmp_07;
    const T tmp_09 = 0.5 * (tmp_04 + tmp_08);

    array[ 3 ] = tmp_00 + tmp_09 - tmp_05;
    array[ 0 ] = tmp_01 + tmp_09 - tmp_08;
    array[ 1 ] = tmp_02 + tmp_09 - tmp_07;
    array[ 2 ] = tmp_03 + tmp_09 - tmp_06;

    #endif
}





template < typename T >
quaternion< T >
quaternion< T >::operator-() const
{
    return quaternion( -array[0], -array[1], -array[2], -array[3] );
}



template < typename T >
void quaternion< T >::operator+=( const quaternion< T >& q )
{
    array[ 0 ] += q.array[ 0 ];
    array[ 1 ] += q.array[ 1 ];
    array[ 2 ] += q.array[ 2 ];
    array[ 3 ] += q.array[ 3 ];
}



template < typename T >
void quaternion< T >::operator-=( const quaternion< T >& q )
{
    array[ 0 ] -= q.array[ 0 ];
    array[ 1 ] -= q.array[ 1 ];
    array[ 2 ] -= q.array[ 2 ];
    array[ 3 ] -= q.array[ 3 ];
}



//
// quaternion/scalar operations
//

template < typename T >
quaternion< T >
quaternion< T >::operator*( const T a_ ) const
{
    return quaternion( array[0] * a_, array[1] * a_, array[2] * a_, array[3] * a_ );
}



template < typename T >
quaternion< T >
quaternion< T >::operator/( T a_ ) const
{
    if ( a_ == 0.0 )
        throw std::runtime_error( "Division by zero." );

    a_ = 1.0 / a_;
    return quaternion( array[0] * a_, array[1] * a_, array[2] * a_, array[3] * a_ );
}



template < typename T >
void quaternion< T >::operator*=( T q )
{
    array[ 0 ] *= q;
    array[ 1 ] *= q;
    array[ 2 ] *= q;
    array[ 3 ] *= q;
}



template < typename T >
void quaternion< T >::operator/=( T q )
{
    if ( q == 0.0 )
        throw std::runtime_error( "Division by zero" );

    q = 1.0f / q;
    this->operator*=( q );
}


template < typename T >
vector< 3, T > quaternion< T >::cross( const quaternion< T >& bb ) const
{
    vector< 3, T > result;

    result.array[ 0 ] = array[1] * bb.array[2] - array[2] * bb.array[1];
    result.array[ 1 ] = array[2] * bb.array[0] - array[0] * bb.array[2];
    result.array[ 2 ] = array[0] * bb.array[1] - array[1] * bb.array[0];

    return result;
}



template < typename T >
T quaternion< T >::dot( const quaternion< T >& q ) const
{
    return array[3] * q.array[3] + array[0] * q.array[0] + array[1] * q.array[1] + array[2] * q.array[2];
}



template < typename T >
T quaternion< T >::
dot( const quaternion< T >& p, const quaternion< T >& q )
{
    return p.array[3] * q.array[3] + p.array[0] * q.array[0] + p.array[1] * q.array[1] + p.array[2] * q.array[2];
}



template < typename T >
void quaternion< T >::normal( const quaternion< T >& aa,
                              const quaternion< T >& bb,
                              const quaternion< T >& cc,
                              const quaternion< T >& dd )
{
    //right hand system, CCW triangle
    const quaternion< T > quat_t = bb - aa;
    const quaternion< T > quat_u = cc - aa;
    const quaternion< T > quat_v = dd - aa;
    cross( quat_t );
    cross( quat_u );
    cross( quat_v );
    normalize();
}



template < typename T >
quaternion< T > quaternion< T >::normal( const quaternion< T >& aa,
                                         const quaternion< T >& bb,
                                         const quaternion< T >& cc )
{
    quaternion< T > tmp;
    tmp.normal( *this, aa, bb, cc );
    return tmp;
}

template < typename T >
Matrix< 3, 3, T > quaternion< T >::get_rotation_matrix() const
{
    Matrix< 3, 3, T > result;
    get_rotation_matrix< 3 >( result );
    return result;
}



template < typename T > template< size_t D >
void quaternion< T >::get_rotation_matrix( Matrix< D, D, T >& M ) const
{
    T w2 = array[3] * array[3];
    T x2 = array[0] * array[0];
    T y2 = array[1] * array[1];
    T z2 = array[2] * array[2];
    T wx = array[3] * array[0];
    T wy = array[3] * array[1];
    T wz = array[3] * array[2];
    T xy = array[0] * array[1];
    T xz = array[0] * array[2];
    T yz = array[1] * array[2];

    M( 0, 0 ) = w2 + x2 - y2 - z2;
    M( 0, 1 ) = 2. * (xy - wz);
    M( 0, 2 ) = 2. * (xz + wy);
    M( 1, 0 ) = 2. * (xy + wz);
    M( 1, 1 ) = w2 - x2 + y2 - z2;
    M( 1, 2 ) = 2. * (yz - wx);
    M( 2, 0 ) = 2. * (xz - wy);
    M( 2, 1 ) = 2. * (yz + wx);
    M( 2, 2 ) = w2 - x2 - y2 + z2;

}

template< typename T >
quaternion< T > quaternion< T >::
slerp( T a, const quaternion< T >& p, const quaternion< T >& q, const T epsilon )
{
    quaternion< T > px = p.get_normalized();
    quaternion< T > qx = q.get_normalized();

    T cosine = px.dot( qx );

    // check if inverted rotation is needed
    if ( cosine < 0.0 )
    {
        cosine = -cosine;
        qx = -qx;
    }

    const T abs_cos = static_cast< T >( fabs(cosine) );
    const T one_x   = static_cast< T >( 1. - epsilon );
    if( abs_cos < one_x )
    {
        // standard slerp
        T sine = std::sqrt( 1. - ( cosine * cosine ) );
        T angle = atan2( sine, cosine );
        T coeff1 = std::sin( ( 1.0 - a ) * angle) / sine;
        T coeff2 = std::sin( a * angle ) / sine;

        qx *= coeff2;
        px *= coeff1;

        px += qx;
    }
    else
    {
        // linear interpolation for very small angles
        px *= 1. - a;
        qx *= a;
        px += qx;
        px.normalize();
    }

    return px;
}


template < typename T >
quaternion< T >& quaternion< T >::operator=(const quaternion& other)
{
    memcpy( array, other.array, 4 * sizeof( T ) );
    return *this;
}

}
#endif
