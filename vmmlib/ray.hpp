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
 * Implements ray-object intersection, storing the ray parameters as attributes
 * and implementing intersection tests against different types of objects
 * (e.g. sphere)
 */

#ifndef __VMML__RAY__HPP__
#define __VMML__RAY__HPP__

#include <vmmlib/vector.hpp>

namespace vmml
{
template< typename T > class Ray
{
public:
    typedef vector< 3, T > vec3;
    typedef vector< 4, T > vec4;

    /**
     * Constructor
     *
     * @param origin Ray origin
     * @param direction Ray direction
     */
    Ray( const vec3& origin, const vec3& direction )
        : _origin ( origin )
        , _direction ( vmml::normalize( direction ))
    {}
    ~Ray() {}

    /**
     * Ray-Sphere Intersection.
     * Optimized solution from "Real-time Rendering 3rd Edition"
     *
     * @param center Sphere center
     * @param radius Sphere radius
     * @param distance Intersection distance
     *
     * @return The distance from the ray origin to the intersection, or a
     *         negative value if there is no intersection.
     */
    T test( const vec4& sphere ) const;

private:
    const vec3 _origin;
    const vec3 _direction;

};


template< typename T >
T Ray< T >::test( const vec4& sphere ) const
{
    const vec3 center = vec3(sphere.x(), sphere.y(), sphere.z());
    const T radius = sphere.w();

    const vec3 centerVec = center - _origin;
    const T vecProjection = centerVec.dot(_direction);

    const T sqDistance = centerVec.squared_length();
    const T sqRadius = radius * radius;

    /** Sphere behind the ray origin and ray origin outside the sphere */
    if( vecProjection < 0 && sqDistance > sqRadius )
        return -1.f;

    /** Squared distance from sphere center to the projection */
    const T sqCenterToProj = sqDistance - vecProjection * vecProjection;

    if( sqCenterToProj > sqRadius )
        return -1.f;

    /** Distance from the sphere center to the surface along the ray direction*/
    const T distSurface = sqrt( sqRadius - sqCenterToProj );

    if(sqDistance > sqRadius)
        return vecProjection - distSurface;

    return vecProjection + distSurface;
}

} // namespace vmml

#endif
