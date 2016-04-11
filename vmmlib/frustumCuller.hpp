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

#ifndef __VMML__FRUSTUMCULLER__HPP__
#define __VMML__FRUSTUMCULLER__HPP__

#include <vmmlib/aabb.hpp> // inline parameter
#include <vmmlib/matrix.hpp> // inline parameter
#include <vmmlib/vector.hpp> // member
#include <vmmlib/visibility.hpp> // return value

// - declaration -

namespace vmml
{

/** Helper class to perform OpenGL view frustum culling. */
template< typename T > class FrustumCuller
{
public:
    typedef vector< 2, T > vec2;
    typedef vector< 3, T > vec3;
    typedef vector< 4, T > vec4;

    /** Construct a new frustum culler */
    FrustumCuller() {}

    /** Construct a frustum culler using a 4x4 projection*model*view matrix. */
    explicit FrustumCuller( const Matrix< 4, 4, T >& projModelView );

    /**
     * Construct a frustum culler using the eight frustum corner points.
     * Corner naming is n(ear)|f(ar), l(eft)|r(ight), t(op)|b(ottom)
     */
    FrustumCuller( const vec3& nlt, const vec3& nrt,
                   const vec3& nlb, const vec3& nrb,
                   const vec3& flt, const vec3& frt,
                   const vec3& flb, const vec3& frb );

    /** Destruct this frustum culler. */
    ~FrustumCuller(){}

    /** @return the visibility of the given sphere */
    Visibility test( const vec4& sphere ) const;

    /** @return the visibility of the axis-aligned bounding box */
    Visibility test( const AABB< T >& aabb ) const;

    /** @return the plane equation of the current near plane. */
    const vec4& getNearPlane() const { return _nearPlane; }

    friend std::ostream& operator << (std::ostream& os, const FrustumCuller& f)
    {
        return os << "Frustum cull planes: " << std::endl
                  << "    left   " << f._leftPlane << std::endl
                  << "    right  " << f._rightPlane << std::endl
                  << "    top    " << f._topPlane << std::endl
                  << "    bottom " << f._bottomPlane << std::endl
                  << "    near   " << f._nearPlane << std::endl
                  << "    far    " << f._farPlane << std::endl;
    }

private:
    inline void _normalizePlane( vec4& plane ) const;
    inline Visibility _test( const vec4& plane, const vec3& middle,
                             const vec3& size_2 ) const;
    vec4    _leftPlane;
    vec4    _rightPlane;
    vec4    _bottomPlane;
    vec4    _topPlane;
    vec4    _nearPlane;
    vec4    _farPlane;

}; // class FrustumCuller
} // namespace vmml

// - implementation - //


namespace vmml
{

/**
 * Setup the culler by extracting the frustum planes from the projection
 * matrix. The projection matrix should contain the viewing transformation.
 */
template < typename T >
FrustumCuller< T >::FrustumCuller( const Matrix< 4, 4, T >& projModelView )
{
    // See http://www2.ravensoft.com/users/ggribb/plane%20extraction.pdf pp.5
    const vec4& row0 = projModelView.getRow( 0 );
    const vec4& row1 = projModelView.getRow( 1 );
    const vec4& row2 = projModelView.getRow( 2 );
    const vec4& row3 = projModelView.getRow( 3 );

    _leftPlane   = row3 + row0;
    _rightPlane  = row3 - row0;
    _bottomPlane = row3 + row1;
    _topPlane    = row3 - row1;
    _nearPlane   = row3 + row2;
    _farPlane    = row3 - row2;

    _normalizePlane( _leftPlane );
    _normalizePlane( _rightPlane );
    _normalizePlane( _bottomPlane );
    _normalizePlane( _topPlane );
    _normalizePlane( _nearPlane );
    _normalizePlane( _farPlane );
}

template < typename T >
FrustumCuller< T >::FrustumCuller( const vec3& a, const vec3& b,
                                   const vec3& c, const vec3& d,
                                   const vec3& e, const vec3& f,
                                   const vec3& g, const vec3& h )
{
    //   e_____f
    //  /____ /|
    // | a b | |
    // | c d |/h
    //  -----
    // CCW winding
    _leftPlane   = compute_plane( c, a, e );
    _rightPlane  = compute_plane( f, b, d );
    _bottomPlane = compute_plane( h, d, c );
    _topPlane    = compute_plane( a, b, f );
    _nearPlane   = compute_plane( b, a, c );
    _farPlane    = compute_plane( g, e, f );
}

template < typename T >
inline void
FrustumCuller< T >::_normalizePlane( vector< 4, T >& plane ) const
{
    const vec3& v3 = plane.template get_sub_vector< 3, 0 >();
    const T len_i = 1.0 / v3.length();
    plane.x() *= len_i;
    plane.y() *= len_i;
    plane.z() *= len_i;
    plane.w() *= len_i;
}


template < typename T > Visibility
FrustumCuller< T >::test( const vector< 4, T >& sphere ) const
{
    Visibility visibility = VISIBILITY_FULL;

    // see http://www.flipcode.com/articles/article_frustumculling.shtml
    // distance = plane.normal . sphere.center + plane.distance
    // Test all planes:
    // - if sphere behind plane: not visible
    // - if sphere intersects one plane: partially visible
    // - else: fully visible

    T distance = _leftPlane.x() * sphere.x() + _leftPlane.y() * sphere.y() +
                 _leftPlane.z() * sphere.z() + _leftPlane.w();
    if( distance <= -sphere.w() )
        return VISIBILITY_NONE;
    if( distance < sphere.w() )
        visibility = VISIBILITY_PARTIAL;

    distance = _rightPlane.x() * sphere.x() + _rightPlane.y() * sphere.y() +
               _rightPlane.z() * sphere.z() + _rightPlane.w();
    if( distance <= -sphere.w() )
        return VISIBILITY_NONE;
    if( distance < sphere.w() )
        visibility = VISIBILITY_PARTIAL;

    distance = _bottomPlane.x() * sphere.x() + _bottomPlane.y() * sphere.y() +
               _bottomPlane.z() * sphere.z() + _bottomPlane.w();
    if( distance <= -sphere.w() )
        return VISIBILITY_NONE;
    if( distance < sphere.w() )
        visibility = VISIBILITY_PARTIAL;

    distance = _topPlane.x() * sphere.x() + _topPlane.y() * sphere.y() +
               _topPlane.z() * sphere.z() + _topPlane.w();
    if( distance <= -sphere.w() )
        return VISIBILITY_NONE;
    if( distance < sphere.w() )
        visibility = VISIBILITY_PARTIAL;

    distance = _nearPlane.x() * sphere.x() + _nearPlane.y() * sphere.y() +
               _nearPlane.z() * sphere.z() + _nearPlane.w();

    if( distance <= -sphere.w() )
        return VISIBILITY_NONE;
    if( distance < sphere.w() )
        visibility = VISIBILITY_PARTIAL;

    distance = _farPlane.x() * sphere.x() + _farPlane.y() * sphere.y() +
               _farPlane.z() * sphere.z() + _farPlane.w();
    if( distance <= -sphere.w() )
        return VISIBILITY_NONE;
    if( distance < sphere.w() )
        visibility = VISIBILITY_PARTIAL;

    return visibility;
}

template < typename T >
Visibility FrustumCuller< T >::_test( const vec4& plane, const vec3& middle,
                                      const vec3& extent ) const
{
    // http://www.cescg.org/CESCG-2002/DSykoraJJelinek/index.html
    const T d = plane.dot( middle );
    const T n = extent.x() * fabs( plane.x( )) +
                extent.y() * fabs( plane.y( )) +
                extent.z() * fabs( plane.z( ));

    if( d - n >= 0 )
        return VISIBILITY_FULL;
    if( d + n > 0 )
        return VISIBILITY_PARTIAL;
    return VISIBILITY_NONE;
}

template < typename T >
Visibility FrustumCuller< T >::test( const AABB< T >& aabb ) const
{
    Visibility result = VISIBILITY_FULL;
    const vec3& middle = aabb.getCenter();
    const vec3& extent = aabb.getSize() * 0.5f;
    switch( _test( _leftPlane, middle, extent ))
    {
        case VISIBILITY_FULL: break;
        case VISIBILITY_PARTIAL: result = VISIBILITY_PARTIAL; break;
        case VISIBILITY_NONE: return VISIBILITY_NONE;
    }

    switch( _test( _rightPlane, middle, extent ))
    {
        case VISIBILITY_FULL: break;
        case VISIBILITY_PARTIAL: result = VISIBILITY_PARTIAL; break;
        case VISIBILITY_NONE: return VISIBILITY_NONE;
    }

    switch( _test( _bottomPlane, middle, extent ))
    {
        case VISIBILITY_FULL: break;
        case VISIBILITY_PARTIAL: result = VISIBILITY_PARTIAL; break;
        case VISIBILITY_NONE: return VISIBILITY_NONE;
    }

    switch( _test( _topPlane, middle, extent ))
    {
        case VISIBILITY_FULL: break;
        case VISIBILITY_PARTIAL: result = VISIBILITY_PARTIAL; break;
        case VISIBILITY_NONE: return VISIBILITY_NONE;
    }

    switch( _test( _nearPlane, middle, extent ))
    {
        case VISIBILITY_FULL: break;
        case VISIBILITY_PARTIAL: result = VISIBILITY_PARTIAL; break;
        case VISIBILITY_NONE: return VISIBILITY_NONE;
    }

    switch( _test( _farPlane, middle, extent ))
    {
        case VISIBILITY_FULL: break;
        case VISIBILITY_PARTIAL: result = VISIBILITY_PARTIAL; break;
        case VISIBILITY_NONE: return VISIBILITY_NONE;
    }

    return result;
}

} // namespace vmml

#endif // include protection
