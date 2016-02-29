
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

#include <vmmlib/quaternion.hpp>

#define BOOST_TEST_MODULE quaternion
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( construction )
{
    const vmml::Quaternionf q1;
    BOOST_CHECK_EQUAL( q1, vmml::Quaternionf( 0, 0, 0, 1 ));

    const vmml::Quaternionf q2( q1 );
    BOOST_CHECK_EQUAL( q2, vmml::Quaternionf( 0, 0, 0, 1 ));
    BOOST_CHECK_EQUAL( q2, q1 );

    const vmml::Quaternionf q3( 0.3, vmml::Vector3f( .5, .5, 1 ));
    BOOST_CHECK_MESSAGE( q3.equals( vmml::Quaternionf( 0.0610079f, 0.0610079f,
                                                       0.122016f, 0.988771f ),
                                    0.00001 ), q3 );
}

BOOST_AUTO_TEST_CASE( rotate )
{
    const vmml::Quaternionf q1( 0.3, vmml::Vector3f( .5, .5, 1 ));
    BOOST_CHECK_MESSAGE( (q1 * q1).equals(
                vmml::Quaternionf( 0.6, vmml::Vector3f( .5, .5, 1 )), 0.00001 ),
                         vmml::Quaternionf( 0.6, vmml::Vector3f( .5, .5, 1 )));
}
