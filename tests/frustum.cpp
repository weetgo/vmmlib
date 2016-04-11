
/* Copyright (c) 2014-2016, Stefan.Eilemann@epfl.ch
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

#include <vmmlib/frustum.hpp>
#include <vmmlib/frustumCuller.hpp>
#include <vmmlib/types.hpp>

#define BOOST_TEST_MODULE frustum
#include <boost/test/unit_test.hpp>

static void _testCull( const vmml::FrustumCullerf& fc )
{
    const vmml::vector< 4, float > sphereIn( 0.f, 0.f, -10.f, 1.f );
    const vmml::vector< 4, float > sphereOut( 0.f, 0.f, 0.f, .5f );
    const vmml::vector< 4, float > sphereBorder( 0.f, 0.f, -1.f, 1.f );

    BOOST_CHECK_EQUAL( fc.test( sphereIn ), vmml::VISIBILITY_FULL );
    BOOST_CHECK_EQUAL( fc.test( sphereOut ), vmml::VISIBILITY_NONE );
    BOOST_CHECK_EQUAL( fc.test( sphereBorder ),
                       vmml::VISIBILITY_PARTIAL );

    vmml::AABBf aabbIn( vmml::Vector3f( -1.f, -1.f, -4.f ),
                        vmml::Vector3f( 1.f, 1.f, -2.f ));
    vmml::AABBf aabbOut( vmml::Vector3f( -1.f, -1.f, -.5f ),
                         vmml::Vector3f( 1.f, 1.f, 0.f ));
    vmml::AABBf aabbBorder( vmml::Vector3f( -1.f, -1.f, -1.5f ),
                            vmml::Vector3f( 1.f, 1.f, -.5f ));

    BOOST_CHECK_EQUAL( fc.test( aabbIn ), vmml::VISIBILITY_FULL );
    BOOST_CHECK_EQUAL( fc.test( aabbOut ), vmml::VISIBILITY_NONE );
    BOOST_CHECK_EQUAL( fc.test( aabbBorder ), vmml::VISIBILITY_PARTIAL );
}

BOOST_AUTO_TEST_CASE( convert )
{
    const vmml::Frustumf f1;
    const vmml::Frustumf f2( f1.computePerspectiveMatrix( ));
    BOOST_CHECK_MESSAGE( f1.equals( f2, 0.0001f ), f2 );
}

BOOST_AUTO_TEST_CASE( base )
{
    const vmml::Frustumf frustum( -1.f, 1., -1.f, 1., 1.f, 100.f );
    const vmml::Matrix< 4, 4, float > mvp = frustum.computePerspectiveMatrix();

    const vmml::FrustumCullerf fc( mvp );
    _testCull( fc );

    //   e_____f
    //  /     /|
    // | a b | |
    // | c d |/h
    //  -----
    const vmml::Vector3f a( -1.f,  1.f, -1.f );
    const vmml::Vector3f b(  1.f,  1.f, -1.f );
    const vmml::Vector3f c( -1.f, -1.f, -1.f );
    const vmml::Vector3f d(  1.f, -1.f, -1.f );
    const vmml::Vector3f e( -100.f,  100.f, -100.f );
    const vmml::Vector3f f(  100.f,  100.f, -100.f );
    const vmml::Vector3f g( -100.f, -100.f, -100.f );
    const vmml::Vector3f h(  100.f, -100.f, -100.f );

    const vmml::FrustumCullerf fc2( a, b, c, d, e, f, g, h );
    _testCull( fc2 );
}
