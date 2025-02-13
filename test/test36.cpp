#include "test.h"
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::constants;

template<typename X, typename Y, typename Z>
static auto f(const X&_rX, const Y&_rY, const Z&_rZ)
{	return log(sqr(_rX) + sqr(_rY) + sqr(_rZ));
}
static constexpr auto X0 = 1.1;
static constexpr auto Y0 = 1.2;
static constexpr auto Z0 = 1.3;

BOOST_AUTO_TEST_CASE(taylor_36)
{	using namespace taylor;
	using namespace boost::mp11;
	const auto sMap = read("data36.txt");
	constexpr std::size_t MAX = 3;
	auto const sX = ctaylor<makeIndependent<0>, MAX>(X0, true);
	auto const sY = ctaylor<makeIndependent<1>, MAX>(Y0, true);
	auto const sZ = ctaylor<makeIndependent<2>, MAX>(Z0, true);
	const auto s5 = f(sX, sY, sZ);
	__EQUAL_TAYLOR__();
}
BOOST_AUTO_TEST_CASE(jacobian_36)
{	using namespace jacobian;
	using namespace boost::mp11;
	const auto sMap = read("data36.txt", true);
	//double const K = 100.0;  // Strike price.
	//constexpr std::size_t MAX = 2;
	auto const sX = cjacobian<mp_list<mp_size_t<0> > >(X0, true);
	auto const sY = cjacobian<mp_list<mp_size_t<1> > >(Y0, true);
	auto const sZ = cjacobian<mp_list<mp_size_t<2> > >(Z0, true);
	const auto s5 = f(sX, sY, sZ);
	__EQUAL_JACOBIAN__();
}
