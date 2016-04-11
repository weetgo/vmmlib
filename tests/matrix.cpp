
/* Copyright (c) 2016, Stefan.Eilemann@epfl.ch
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - Neither the name of Eyescale Software GmbH nor the names of its
 *   contributors may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <vmmlib/matrix.hpp>
#include <vmmlib/types.hpp>

#define BOOST_TEST_MODULE matrix
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(base)
{
    vmml::Matrix4d matrix;
    for( size_t i = 0; i < 16; ++i )
        if( (i%5) == 0 )
            BOOST_CHECK_EQUAL( matrix.array[ i ], 1. );
        else
            BOOST_CHECK_EQUAL( matrix.array[ i ], 0. );

    double data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
    const vmml::Matrix4d matrix2( data, data + 16 );
    matrix = matrix2;
    for( size_t i = 0; i < 16; ++i )
    {
        BOOST_CHECK_EQUAL( matrix2.array[ i ], i+1. );
        BOOST_CHECK_EQUAL( matrix.array[ i ], matrix2.array[ i ] );
    }
}

BOOST_AUTO_TEST_CASE( construction )
{
    vmml::Matrix4f m1( vmml::Quaternionf(), vmml::Vector3f( 1.f, 2.f, 3.f ));
    const float data[] = { 1, 0, 0, 0,
                           0, 1, 0, 0,
                           0, 0, 1, 0,
                           1.f, 2.f, 3.f, 1 };
    BOOST_CHECK_EQUAL( m1, vmml::Matrix4f( data, data + 16 ));
}

BOOST_AUTO_TEST_CASE( lookat )
{
    const vmml::Vector3f eye( 1, 2, 3 );
    const vmml::Vector3f lookAt( eye + vmml::Vector3f( 0, 1, 0 ));
    const vmml::Vector3f up( 0, 0, 1 );

    const vmml::Matrix4f m1( eye, lookAt, up );

    vmml::Vector3f newEye, newLookAt, newUp;
    m1.getLookAt( newEye, newLookAt, newUp );
    BOOST_CHECK_EQUAL( eye, newEye );
    BOOST_CHECK_EQUAL( lookAt, newLookAt );
    BOOST_CHECK_EQUAL( up, newUp );
}

// Verify code by instantiating some templates:
template class vmml::Matrix< 1, 1, float >;
template class vmml::Matrix< 2, 2, double >;
template class vmml::Matrix< 3, 3, short >;
template class vmml::Matrix< 3, 4, int >;
template class vmml::Matrix< 4, 4, long >;
