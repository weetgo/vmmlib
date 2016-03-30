/*
 * Copyright (c) 2016 Stefan.Eilemann@epfl.ch
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
#ifndef __VMML__VMMLIB_TYPES__HPP__
#define __VMML__VMMLIB_TYPES__HPP__

#include <sys/types.h>
#ifdef _MSC_VER
#  include <basetsd.h>
typedef UINT8 uint8_t;
#else
#  include <stdint.h>
#endif

namespace vmml
{
template< size_t M, size_t N, typename T > class Matrix;
template< size_t M, typename T > class vector;
template< typename T > class AABB;
template< typename T > class Frustum;
template< typename T > class FrustumCuller;
template< typename T > class Quaternion;
template< typename T > class Ray;

typedef Matrix< 3, 3, double > Matrix3d; //!< A 3x3 double matrix
typedef Matrix< 3, 3, float >  Matrix3f; //!< A 3x3 float matrix
typedef Matrix< 4, 4, double > Matrix4d; //!< A 4x4 double matrix
typedef Matrix< 4, 4, float >  Matrix4f; //!< A 4x4 float matrix

typedef vector< 2, float > Vector2f; //!< A 2-component float vector
typedef vector< 2, int > Vector2i; //!< A 2-component int vector
typedef vector< 2, unsigned > Vector2ui; //!< A 2-component unsigned int vector
typedef vector< 3, double > Vector3d; //!< A 3-component double vector
typedef vector< 3, float > Vector3f; //!< A 3-component float vector
typedef vector< 3, int > Vector3i; //!< A 3-component int vector
typedef vector< 3, uint8_t > Vector3ub; //!< A 3-component byte vector
typedef vector< 3, unsigned > Vector3ui; //!< A 3-component unsigned int vector
typedef vector< 4, double > Vector4d; //!< A 4-component double vector
typedef vector< 4, float > Vector4f; //!< A 4-component float vector
typedef vector< 4, int > Vector4i; //!< A 4-component int vector
typedef vector< 4, uint8_t > Vector4ub; //!< A 4-component byte vector
typedef vector< 4, unsigned > Vector4ui; //!< A 4-component unsigned int vector

typedef Quaternion< double > Quaterniond; //!< A double quaternion
typedef Quaternion< float >  Quaternionf; //!< A float quaternion

typedef Frustum< double > Frustumd; //!< A double frustum
typedef Frustum< float >  Frustumf; //!< A float frustum

typedef FrustumCuller< double > FrustumCullerd; //!< A double frustum culler
typedef FrustumCuller< float >  FrustumCullerf; //!< A float frustum culler

typedef AABB< double > AABBd; //!< A double bounding box
typedef AABB< float >  AABBf; //!< A float bounding box

typedef Ray< double > Rayd; //!< A double ray
typedef Ray< float > Rayf; //!< A float ray

}

#endif
