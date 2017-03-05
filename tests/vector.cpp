
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

#include <vmmlib/vector.hpp>
#include <vmmlib/types.hpp>

#define BOOST_TEST_MODULE vector
#include <boost/test/unit_test.hpp>

#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327
#endif

using namespace vmml;

BOOST_AUTO_TEST_CASE(base)
{
    vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    v.iter_set( data, data+4 );

    // tests copyFrom1DimCArray function
    size_t tmp = 1;
    for( size_t index = 0; index < 4; ++index, ++tmp )
    {
        BOOST_CHECK(v.at( index ) == tmp);
    }

    tmp = 4;
    float dataf[] = { 4, 3, 2, 1 };
    v.iter_set( dataf, dataf + 4 );
    for( size_t index = 0; index < 4; ++index, --tmp )
    {
        BOOST_CHECK(v.at( index ) == tmp);
    }
}

BOOST_AUTO_TEST_CASE(plus)
{
    double data[] = { 1, 2, 3, 4 };
    double datad[] = { 4, 3, 2, 1 };

    // tests operator+ function
    vector< 4, double > v( data );
    vector< 4, double > v_other( datad );
    vector< 4, double > v_result = v + v_other;

    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == 5);

    v_result = v;
    v_result += v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == 5);

    v = vector< 4, double >( data );
    v_result = v + 2.;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == index + 3);

    v_result = v;
    v_result += 2;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == index + 3);
}

BOOST_AUTO_TEST_CASE(minus)
{
    double data[] = { 1, 2, 3, 4 };
    double datad[] = { 1, 2, 3, 4 };

    // tests operator- function
    vector< 4, double > v( data );
    vector< 4, double > v_other( datad );
    vector< 4, double > v_result;

    v_result = v - v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == 0);

    v_result = v;
    v_result -= v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == 0);

    v_result = v - 1.0;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == index);

    v_result = v;
    v_result -= 1.0;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == index);
}

BOOST_AUTO_TEST_CASE(times)
{
    double data[] = { 1, 2, 3, 4 };
    double datad[] = { 24, 12, 8, 6 };

    // tests operator* function
    vector< 4, double > v( data );
    vector< 4, double > v_other( datad );
    vector< 4, double > v_result;

    v_result = v * v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == 24);

    v_result = v;
    v_result *= v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == 24);

    v_result = v * 2.0;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == v.at( index ) * 2.0);

    v_result = v;
    v_result *= 2.0;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(v_result.at( index ) == v.at( index ) * 2.0);
}

BOOST_AUTO_TEST_CASE(divide)
{
    double data[] = { 1, 2, 3, 4 };
    double datad[] = { 2, 4, 6, 8 };

    // tests operator/ function
    vector< 4, double > v( data );
    vector< 4, double > v_other( datad );
    vector< 4, double > v_result;

    v_result = v / v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(( v_result.at( index ) - 0.5 ) < 1e-12);

    v_result = v;
    v_result /= v_other;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(( v_result.at( index ) - 0.5 ) < 1e-12);

    v_result = v / 1.5;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(( v_result.at( index ) - ( v.at( index ) / 1.5 ) ) < 1e-12);

    v_result = v;
    v_result /= 1.5;
    for( size_t index = 0; index < 4; ++index )
        BOOST_CHECK(( v_result.at( index ) - ( v.at( index ) / 1.5 ) ) < 1e-12);
}

BOOST_AUTO_TEST_CASE(vec_norm)
{
    vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    // tests norm / normSquared (length/lengthSquared) computation
    vector< 4, double > vec( data );

    double normSquared = vec.squared_length();
    BOOST_CHECK(normSquared == 1 * 1 + 2 * 2 + 3 * 3 + 4 * 4);

    double norm = vec.length();
    BOOST_CHECK(sqrt( normSquared ) == norm);

    // tests normalize
    vec = vector< 4, double >( data );
    vec.normalize();
    BOOST_CHECK_CLOSE( vec.length(), 1.0, 0.0000001 );

    // constructor tests
    double vData[] = { 1, 2, 3, 4 };
    vector< 4, double > v4( 1, 2, 3, 4 );

    vector< 2, double > v2C( vData );
    vector< 2, double > v2( 1, 2 );

    BOOST_CHECK(v2 == v2C );

    vector< 3, double > v3C( vData );
    vector< 3, double > v3( 1, 2, 3 );

    BOOST_CHECK(v3 == v3C );

    vector< 4, double > v4C( vData );

    BOOST_CHECK(v4 == v4C);

    double vData2[] = { 23, 23, 23, 23 };
    v4C = vector< 4, double >( vData2 );

    vector< 4, double > v4_( 23 );
    BOOST_CHECK(v4_ == v4C);

    v3 = vector< 3, double >( vData );
    v4C = vector< 4, double >( vData );

    vector< 4, double > v4from3_1( v3, vData[ 3 ] );
    BOOST_CHECK(v4from3_1 == v4C);

    double hvData[] = { 1., 2., 3., 0.25 };
    double xvData[] = { 4.0, 8.0, 12.0 };

    vector< 4, double > homogenous;
    homogenous.iter_set( hvData, hvData + 4 );
    vector< 3, double > nonh;
    nonh.iter_set( xvData, xvData + 3 );

    vector< 4, double > htest( nonh );

    // to-homogenous-coordinates ctor
    BOOST_CHECK((htest == vector< 4, double >( 4, 8., 12., 1. ) ));

    vector< 3, double > nhtest( homogenous );

    // from homogenous-coordiates ctor
    BOOST_CHECK(nhtest == nonh );

    // set tests
    double vCData[] = { 2, 3, 4, 5 };
    vec.set( 2, 3, 4, 5 );
    vector< 4, double > vecCorrect( vCData );
    BOOST_CHECK(vec == vecCorrect);

    vec.set( 2 );

    double vCData2[] = { 2, 2, 2, 2 };
    vecCorrect = vector< 4, double >( vCData2 );
    BOOST_CHECK( vec == vecCorrect );

    vector< 3, double > v1( 2, 3, 4 );

    // component accessors
    vector< 4, double > vd( 1, 2, 3, 4 );
    BOOST_CHECK( vd.x() == 1 && vd.y() == 2 && vd.z() == 3 && vd.w() == 4 );
}

BOOST_AUTO_TEST_CASE(dotprod)
{
    // dot product
    vector< 3, float > v0( 1, 2, 3 );
    vector< 3, float > v1( -6, 5, -4 );
    BOOST_CHECK( v0.dot( v1 ) == -8 );
}

BOOST_AUTO_TEST_CASE(crossprod)
{
    vector< 3, float > v0( 1, 2, 3 );
    const vector< 3, float > v1( -6, 5, -4 );
    const vector< 3, float > vcorrect( -23, -14, 17 );

    BOOST_CHECK_EQUAL( cross( v0, v1 ), vcorrect );
    BOOST_CHECK_EQUAL( v0.cross( v1 ), vcorrect );
}

BOOST_AUTO_TEST_CASE(normal)
{
    const vmml::vector< 3, float > v1( 0.f, 0.f, 0.f );
    const vmml::vector< 3, float > v2( 0.f, 1.f, 0.f );
    const vmml::vector< 3, float > v3( 1.f, 0.f, 0.f );
    const vmml::vector< 3, float > n( 0.f, 0.f, -1.f );
    BOOST_CHECK_EQUAL( vmml::compute_normal( v1, v2, v3 ), n );
}

BOOST_AUTO_TEST_CASE(minMax)
{
    vector< 4, float > vf( -1.0f, 3.0f, -99.0f, -0.9f );
    vector< 4, size_t > vui( 0, 5, 2, 4 );

    size_t index = vf.find_min_index();
    float f = vf.find_min();

    BOOST_CHECK( index == 2 && f == -99.0f );

    index = vf.find_max_index();
    f = vf.find_max();
    BOOST_CHECK( index == 1 && f == 3.0f );

    index = vui.find_min_index();
    size_t ui = vui.find_min();
    BOOST_CHECK( index == 0 && ui == 0 );

    index = vui.find_max_index();
    ui = vui.find_max();
    BOOST_CHECK( index == 1 && ui == 5 );
}

BOOST_AUTO_TEST_CASE(product)
{
    const vector< 3, float > v0( 1, 2, 3 );
    BOOST_CHECK_EQUAL( v0.product(), 6 );

    const vector< 3, int > v1( -6, 5, -4 );
    BOOST_CHECK_EQUAL( v1.product(), 120 );
}

BOOST_AUTO_TEST_CASE(tbd1)
{
    vector< 4, float > v1( -1.0f, 3.0f, -99.0f, -0.9f );
    float f = 4.0f;
    vector< 4, float > v_scaled = f * v1;

    BOOST_CHECK(v_scaled == (vector< 4, float >( -4.0f, 12.0f, -396.0f, -3.6f ) ));


    // ???
    vector< 3, float > vf( 3.0, 2.0, 1.0 );
    vector< 3, double > vd( vf );
    vector< 3, double >::const_iterator it = vd.begin(), it_end = vd.end();
    vector< 3, float >::const_iterator fit = vf.begin();
    for( ; it != it_end; ++it, ++fit )
    {
        BOOST_CHECK(*it == *fit);
    }
    vd = vf;
    for( ; it != it_end; ++it, ++fit )
    {
        BOOST_CHECK(*it == *fit);
    }
}

BOOST_AUTO_TEST_CASE(subVector)
{
    vector< 4, float > v4( 3.0, 2.0, 1.0, 1.0 );
    vector< 3, float > v3 = v4.get_sub_vector< 3, 0 >();
    BOOST_CHECK(v3.x() == v4.x() && v3.y() == v4.y());
    v3.normalize();

    BOOST_CHECK_NE( v3.x(), v4.x( ));
    BOOST_CHECK_NE( v3.y(), v4.y( ));

    v4.set_sub_vector< 3, 1 >( v3 );
    BOOST_CHECK_EQUAL( v3.x(), v4.y( ));
    BOOST_CHECK_EQUAL( v3.y(), v4.z( ));
}

BOOST_AUTO_TEST_CASE(tbd2)
{
    //elementwise sqrt
    vector< 4, double > vsq(9.0, 4.0, 1.0, 2.0);
    vector< 4, double > vsq_check( 3.0, 2.0, 1.0, std::sqrt( 2.0 ));
    vsq.sqrt_elementwise();
    BOOST_CHECK_EQUAL( vsq, vsq_check );

    //elementwise sqrt
    vector< 4, float > vr( 9.0, 4.0, 1.0, 2.0 );
    vector< 4, float > vr_check( 0.1111111119389534, 0.25, 1, 0.5 );
    vr.reciprocal();
    BOOST_CHECK(vr == vr_check);
}

BOOST_AUTO_TEST_CASE(l2norm)
{
    vector< 4, float > vr( 9.0, 4.0, 1.0, 2.0 );
    double v_norm_check = 10.09950493836208;
    double v_norm = vr.norm();

    BOOST_CHECK((v_norm - v_norm_check) < 0.0001);
}

BOOST_AUTO_TEST_CASE(rotateVec)
{
    vmml::Vector3f vector = vmml::Vector3f::forward();
    vector.rotate( float( M_PI ), vmml::Vector3f::up( ));
    BOOST_CHECK_MESSAGE( vector.equals( vmml::Vector3f::backward( )), vector );
    BOOST_CHECK_MESSAGE( vmml::rotate( vector, float( M_PI ),
                                       vmml::Vector3f::left( )).equals(
                                           vmml::Vector3f::forward( )),
                         vmml::rotate( vector, float( M_PI ),
                                       vmml::Vector3f::left( )));
}
